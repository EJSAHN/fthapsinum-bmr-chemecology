#!/usr/bin/env python3
from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def main() -> None:
    parser = argparse.ArgumentParser(description="Assess fungal gene-count quality and mock-associated loci.")
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--sample-metrics", type=Path, required=True)
    parser.add_argument("--gene-metadata", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-count", type=int, default=5)
    parser.add_argument("--min-samples", type=int, default=3)
    parser.add_argument("--inoculated-label", default="fus")
    parser.add_argument("--control-label", default="pdb")
    args = parser.parse_args()

    manifest = pd.read_csv(args.manifest, sep="\t", dtype=str, keep_default_na=False)
    counts = pd.read_csv(args.counts, sep="\t", index_col=0).apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
    metrics = pd.read_csv(args.sample_metrics, sep="\t")
    gene_meta = pd.read_csv(args.gene_metadata, sep="\t", dtype=str, keep_default_na=False)
    samples = manifest["sample_id"].tolist()
    missing = sorted(set(samples) - set(counts.columns))
    if missing:
        raise RuntimeError(f"Count matrix is missing samples: {missing}")
    counts = counts.loc[:, samples]

    inoculated = manifest.loc[manifest["pathogen"].str.lower().eq(args.inoculated_label.lower()), "sample_id"].tolist()
    controls = manifest.loc[manifest["pathogen"].str.lower().eq(args.control_label.lower()), "sample_id"].tolist()
    assigned = metrics.set_index("sample_id")["assigned_reads"].reindex(samples).astype(float)
    cpm = counts.div(assigned.replace(0, np.nan), axis=1) * 1e6

    control_counts = counts[controls]
    inoculated_counts = counts[inoculated]
    control_cpm = cpm[controls]
    inoculated_cpm = cpm[inoculated]
    recurrent_cutoff = max(3, math.ceil(0.15 * len(controls)))

    background = pd.DataFrame(index=counts.index)
    background["control_n_ge1"] = (control_counts >= 1).sum(axis=1)
    background["control_n_ge3"] = (control_counts >= 3).sum(axis=1)
    background["control_n_ge5"] = (control_counts >= 5).sum(axis=1)
    background["control_total_counts"] = control_counts.sum(axis=1)
    background["control_prevalence_ge3"] = background["control_n_ge3"] / len(controls)
    background["control_median_cpm"] = control_cpm.median(axis=1)
    background["control_max_cpm"] = control_cpm.max(axis=1)
    background["inoculated_n_ge5"] = (inoculated_counts >= 5).sum(axis=1)
    background["inoculated_prevalence_ge5"] = background["inoculated_n_ge5"] / len(inoculated)
    background["inoculated_median_cpm"] = inoculated_cpm.median(axis=1)
    background["recurrent_control"] = background["control_n_ge3"] >= recurrent_cutoff
    background["high_control_burden"] = (background["control_total_counts"] >= 25) & (background["control_n_ge1"] >= 4)
    background["exclude"] = background["recurrent_control"] | background["high_control_burden"]
    background["reason"] = np.select(
        [background["recurrent_control"] & background["high_control_burden"], background["recurrent_control"], background["high_control_burden"]],
        ["recurrent_and_high_control", "recurrent_control", "high_control_total"],
        default="",
    )
    background = background.reset_index().rename(columns={background.index.name or "index": "gene_id"})
    background = background.merge(gene_meta, on="gene_id", how="left")

    top_rows = []
    for sample in samples:
        vector = counts[sample].sort_values(ascending=False)
        total = float(vector.sum())
        top_rows.append({"sample_id": sample, "top10_fraction": float(vector.head(10).sum() / total) if total else np.nan})
    sample_qc = metrics.merge(pd.DataFrame(top_rows), on="sample_id", how="left")
    sample_qc = sample_qc.merge(
        manifest[["sample_id", "genotype", "pathogen", "water", "dai", "replicate", "sequencing_run"]],
        on="sample_id", how="left",
    )
    sample_qc["pass"] = (
        pd.to_numeric(sample_qc["assigned_reads"], errors="coerce").ge(5000)
        & pd.to_numeric(sample_qc["detected_ge5"], errors="coerce").ge(200)
        & pd.to_numeric(sample_qc["top10_fraction"], errors="coerce").le(0.70)
    )

    excluded = set(background.loc[background["exclude"], "gene_id"])
    usable = [
        gene for gene in counts.index
        if gene not in excluded and (inoculated_counts.loc[gene] >= args.min_count).sum() >= args.min_samples
    ]

    correlations = []
    fus_meta = manifest[manifest["pathogen"].str.lower().eq(args.inoculated_label.lower())]
    for (genotype, water), subset in fus_meta.groupby(["genotype", "water"]):
        cell_samples = subset["sample_id"].tolist()
        for first, second in itertools.combinations(cell_samples, 2):
            if len(usable) >= 20:
                rho, p = spearmanr(np.log1p(cpm.loc[usable, first]), np.log1p(cpm.loc[usable, second]))
            else:
                rho, p = np.nan, np.nan
            correlations.append({
                "genotype": genotype, "water": water, "sample_a": first, "sample_b": second,
                "genes": len(usable), "rho": rho, "p": p,
            })

    args.output_dir.mkdir(parents=True, exist_ok=True)
    background.to_csv(args.output_dir / "background_genes.tsv", sep="\t", index=False, lineterminator="\n", na_rep="NA")
    sample_qc.to_csv(args.output_dir / "sample_gene_qc.tsv", sep="\t", index=False, lineterminator="\n", na_rep="NA")
    pd.DataFrame({"gene_id": usable}).to_csv(args.output_dir / "usable_genes.tsv", sep="\t", index=False, lineterminator="\n")
    pd.DataFrame(correlations).to_csv(args.output_dir / "within_cell_correlations.tsv", sep="\t", index=False, lineterminator="\n", na_rep="NA")

    fus_qc = sample_qc[sample_qc["pathogen"].str.lower().eq(args.inoculated_label.lower())]
    summary = {
        "usable_genes": len(usable),
        "background_genes": int(background["exclude"].sum()),
        "primary_inoculated_samples": int(fus_qc["pass"].sum()),
        "all_inoculated_samples": len(fus_qc),
        "excluded_samples": fus_qc.loc[~fus_qc["pass"], "sample_id"].tolist(),
        "median_assigned_reads": float(pd.to_numeric(fus_qc["assigned_reads"]).median()),
        "median_detected_ge5": float(pd.to_numeric(fus_qc["detected_ge5"]).median()),
    }
    (args.output_dir / "gene_count_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
