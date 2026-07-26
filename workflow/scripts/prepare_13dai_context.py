#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


def sample_name(column: str) -> str:
    name = Path(column).name
    name = re.sub(r"\.(?:FT|MP)\.mapq20\.primary\.bam$", "", name)
    return re.sub(r"\.bam$", "", name)


def read_featurecounts_counts(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", comment="#")
    if "Geneid" not in frame.columns:
        fallback = pd.read_csv(path, sep="\t", index_col=0)
        fallback.columns = [sample_name(column) for column in fallback.columns]
        return fallback.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
    fixed = {"Geneid", "Chr", "Start", "End", "Strand", "Length"}
    sample_columns = [column for column in frame.columns if column not in fixed]
    renamed = {column: sample_name(column) for column in sample_columns}
    frame = frame.rename(columns=renamed).set_index("Geneid")
    return frame[[renamed[column] for column in sample_columns]].apply(
        pd.to_numeric, errors="coerce"
    ).fillna(0).astype(int)


def read_featurecounts_summary(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    first = frame.columns[0]
    frame = frame.rename(columns={first: "Status"}).set_index("Status")
    frame = frame.rename(columns={column: sample_name(column) for column in frame.columns})
    return frame.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)


def boolean_series(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame.columns:
        return pd.Series(False, index=frame.index)
    return frame[column].astype(str).str.lower().isin({"1", "true", "yes"})


def robust_upper_threshold(values: pd.Series, multiplier: float = 3.0) -> float:
    numeric = pd.to_numeric(values, errors="coerce").dropna().astype(float)
    if numeric.empty:
        return 0.0
    median = float(numeric.median())
    mad = float((numeric - median).abs().median())
    robust_sd = 1.4826 * mad
    return median + multiplier * robust_sd if robust_sd > 0 else max(median * 3.0, median + 1.0)


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare the 13-day fungal context analysis.")
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--featurecounts-summary", type=Path, required=True)
    parser.add_argument("--focus-genes", type=Path, required=True)
    parser.add_argument("--prior-background", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--minimum-assigned-reads", type=int, default=5000)
    parser.add_argument("--minimum-detected-genes", type=int, default=300)
    parser.add_argument("--minimum-target-reads", type=int, default=5000)
    parser.add_argument("--minimum-specificity", type=float, default=0.85)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(args.manifest, sep="\t", dtype=str, keep_default_na=False)
    mapping = pd.read_csv(args.mapping, sep="\t")
    counts = read_featurecounts_counts(args.counts)
    summary = read_featurecounts_summary(args.featurecounts_summary)
    focus_table = pd.read_csv(args.focus_genes, sep="\t", dtype=str, keep_default_na=False)
    focus = focus_table["gene_id"].astype(str).tolist()

    sample_ids = manifest["sample_id"].astype(str).tolist()
    missing_counts = sorted(set(sample_ids) - set(counts.columns))
    if missing_counts:
        raise RuntimeError(f"Count matrix is missing samples: {missing_counts}")
    counts = counts.reindex(columns=sample_ids).fillna(0).astype(int)
    if "Assigned" not in summary.index:
        raise RuntimeError("FeatureCounts summary lacks an Assigned row")
    assigned = pd.to_numeric(summary.loc["Assigned"], errors="coerce").reindex(sample_ids).fillna(0)

    samples = manifest.merge(mapping, on="sample_id", how="left", validate="one_to_one")
    samples["gene_assigned_reads"] = samples["sample_id"].map(assigned).fillna(0).astype(int)
    samples["target_genes_ge5"] = [int((counts[sample] >= 5).sum()) for sample in sample_ids]

    for column in ("target_unique", "target_unique_rpm", "target_specificity", "post_fastp_reads"):
        if column not in samples.columns:
            samples[column] = np.nan
        samples[column] = pd.to_numeric(samples[column], errors="coerce")

    control_mask = samples["pathogen"].str.lower().eq("pdb")
    inoculated_mask = samples["pathogen"].str.lower().eq("fus")
    control_ids = samples.loc[control_mask, "sample_id"].tolist()
    inoculated_ids = samples.loc[inoculated_mask, "sample_id"].tolist()
    if not control_ids or not inoculated_ids:
        raise RuntimeError("Both inoculated and control samples are required")

    background_threshold = robust_upper_threshold(samples.loc[control_mask, "target_unique_rpm"])
    samples["background_burden_outlier"] = control_mask & samples["target_unique_rpm"].gt(background_threshold)
    samples["gene_qc_pass"] = ~inoculated_mask | (
        samples["gene_assigned_reads"].ge(args.minimum_assigned_reads)
        & samples["target_genes_ge5"].ge(args.minimum_detected_genes)
        & samples["target_unique"].ge(args.minimum_target_reads)
        & samples["target_specificity"].fillna(0).ge(args.minimum_specificity)
    )

    inoculated = samples.loc[inoculated_mask].copy()
    supported_counts = inoculated.loc[inoculated["gene_qc_pass"]].groupby(["genotype", "water"]).size()
    required_cells = [("wt", "dry"), ("wt", "wet"), ("bmr12", "dry"), ("bmr12", "wet")]
    qc_balanced = all(int(supported_counts.get(cell, 0)) >= 2 for cell in required_cells)
    primary_ids = (
        inoculated.loc[inoculated["gene_qc_pass"], "sample_id"].tolist()
        if qc_balanced else inoculated_ids
    )
    samples["primary_fus_analysis"] = samples["sample_id"].isin(primary_ids)

    control_counts = counts[control_ids]
    background = pd.DataFrame(index=counts.index)
    background["control_samples_ge5"] = (control_counts >= 5).sum(axis=1)
    background["control_total_counts"] = control_counts.sum(axis=1)
    background["control_max_count"] = control_counts.max(axis=1)
    background["exclude_context"] = (
        background["control_samples_ge5"].ge(2) | background["control_max_count"].ge(50)
    )
    background["exclude_prior"] = False
    if args.prior_background and args.prior_background.exists():
        prior = pd.read_csv(args.prior_background, sep="\t", dtype=str, keep_default_na=False)
        if "gene_id" in prior.columns:
            flag_column = "exclude" if "exclude" in prior.columns else "blacklist" if "blacklist" in prior.columns else None
            if flag_column:
                prior_flags = prior.set_index("gene_id")[flag_column].str.lower().isin({"1", "true", "yes"})
                background["exclude_prior"] = prior_flags.reindex(background.index).fillna(False).astype(bool)
    background["exclude_union"] = background["exclude_context"] | background["exclude_prior"]
    background.index.name = "gene_id"

    inoculated_counts = counts[inoculated_ids]
    detectable = (inoculated_counts >= 5).sum(axis=1).ge(2)
    usable = detectable & ~background["exclude_union"]
    missing_focus = sorted(set(focus) - set(counts.index))
    if missing_focus:
        raise RuntimeError(f"Count matrix is missing focus genes: {missing_focus}")

    samples.to_csv(args.output_dir / "sample_summary.tsv", sep="\t", index=False, lineterminator="\n", na_rep="NA")
    background.reset_index().to_csv(
        args.output_dir / "background_genes.tsv", sep="\t", index=False, lineterminator="\n"
    )
    counts.loc[:, primary_ids].to_csv(
        args.output_dir / "counts_primary.tsv", sep="\t", index=True, index_label="gene_id"
    )
    counts.loc[:, inoculated_ids].to_csv(
        args.output_dir / "counts_all.tsv", sep="\t", index=True, index_label="gene_id"
    )
    samples.loc[samples["sample_id"].isin(primary_ids)].to_csv(
        args.output_dir / "samples_primary.tsv", sep="\t", index=False, lineterminator="\n", na_rep="NA"
    )
    inoculated.to_csv(
        args.output_dir / "samples_all.tsv", sep="\t", index=False, lineterminator="\n", na_rep="NA"
    )
    pd.DataFrame({"gene_id": counts.index, "usable_primary": usable.reindex(counts.index).fillna(False).to_numpy()}).to_csv(
        args.output_dir / "gene_filter.tsv", sep="\t", index=False, lineterminator="\n"
    )
    pd.DataFrame({"module": "aromatic_phenolic_metabolism", "gene_id": focus}).to_csv(
        args.output_dir / "focus_genes.tsv", sep="\t", index=False, lineterminator="\n"
    )

    payload = {
        "libraries": int(len(samples)),
        "inoculated_libraries": int(len(inoculated_ids)),
        "control_libraries": int(len(control_ids)),
        "primary_inoculated_libraries": int(len(primary_ids)),
        "qc_balanced": bool(qc_balanced),
        "background_threshold_rpm": float(background_threshold),
        "focus_genes": int(len(focus)),
        "focus_genes_background_flagged": int(background.reindex(focus)["exclude_union"].fillna(False).sum()),
    }
    (args.output_dir / "context_inputs.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
