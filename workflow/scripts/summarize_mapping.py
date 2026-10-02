#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import mannwhitneyu, spearmanr
from rna_metrics import calculate_metrics, fastp_read_counts


def one_row(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8") as handle:
        return next(csv.DictReader(handle, delimiter="\t"))


def alignment_rate(path: Path) -> float:
    text = path.read_text(encoding="utf-8", errors="replace")
    match = re.search(r"Overall alignment rate:\s*([0-9.]+)%", text, re.I)
    return float(match.group(1)) if match else math.nan


def coverage(path: Path) -> dict[str, float]:
    total = 0
    covered = [0, 0, 0]
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 7:
                continue
            total += int(fields[2]) - int(fields[1])
            for index, value in enumerate(map(int, fields[-3:])):
                covered[index] += value
    return {
        "target_breadth_1x": covered[0] / total if total else math.nan,
        "target_breadth_3x": covered[1] / total if total else math.nan,
        "target_breadth_5x": covered[2] / total if total else math.nan,
    }


def robust_threshold(values: pd.Series) -> float:
    values = pd.to_numeric(values, errors="coerce").dropna()
    median = float(values.median())
    mad = float((values - median).abs().median())
    return median + 3 * 1.4826 * mad if mad > 0 else float(values.quantile(0.95))


def cliffs_delta(first: np.ndarray, second: np.ndarray) -> float:
    if not len(first) or not len(second):
        return math.nan
    greater = sum(x > y for x in first for y in second)
    lower = sum(x < y for x in first for y in second)
    return (greater - lower) / (len(first) * len(second))


def bh(values: pd.Series) -> np.ndarray:
    p = values.to_numpy(float)
    order = np.argsort(p)
    adjusted = p[order] * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty(len(p), dtype=float)
    result[order] = np.minimum(adjusted, 1)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize competitive mapping and pathogen RNA burden.")
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--summary-table", type=Path, help="Prepared per-sample input with explicit pre/post read counts")
    parser.add_argument("--rpm-basis", choices=["input", "post_filter"], default="input")
    parser.add_argument("--specificity-groups", nargs="+", default=None)
    parser.add_argument("--species-dir", type=Path)
    parser.add_argument("--summary-dir", type=Path)
    parser.add_argument("--fastp-dir", type=Path)
    parser.add_argument("--coverage-dir", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--target-group", default="FT")
    parser.add_argument("--inoculated-label", default="fus")
    args = parser.parse_args()

    if args.summary_table:
        if any([args.manifest,args.species_dir,args.summary_dir,args.fastp_dir,args.coverage_dir]):
            parser.error("Use summary-table OR raw summary directories, not both")
        table = pd.read_csv(args.summary_table, sep="\t", keep_default_na=True)
    else:
        if not all([args.manifest,args.species_dir,args.summary_dir,args.fastp_dir,args.coverage_dir]):
            parser.error("All five raw-summary inputs are required without summary-table")
        manifest = pd.read_csv(args.manifest, sep="\t", dtype=str, keep_default_na=False)
        rows = []
        for record in manifest.to_dict("records"):
            sample = record["sample_id"]
            species = dict(record)
            species.update(one_row(args.species_dir / f"{sample}.tsv"))
            species.update(coverage(args.coverage_dir / f"{sample}.thresholds.bed.gz"))
            species["overall_alignment_pct"] = alignment_rate(args.summary_dir / f"{sample}.txt")
            payload = json.loads((args.fastp_dir / f"{sample}.json").read_text())
            species["pre_fastp_reads"], species["post_fastp_reads"] = fastp_read_counts(payload)
            rows.append(species)
        table = pd.DataFrame(rows)
    if table["sample_id"].duplicated().any():
        raise ValueError("Duplicate sample identifiers")

    numeric_columns = {"dai", "sequencing_run", "expected_reads", "expected_bytes"}
    numeric_columns.update(column for column in table.columns if column.endswith(("_mapped", "_q20", "_unique")))
    numeric_columns.update({"all_primary", "unmapped_primary", "mapped_primary", "mapq_pass", "unique_pass",
                            "secondary_alignments", "supplementary_alignments", "duplicate_primary", "min_mapq",
                            "overall_alignment_pct", "pre_fastp_reads", "post_fastp_reads", "target_breadth_1x",
                            "target_breadth_3x", "target_breadth_5x"})
    for column in numeric_columns & set(table.columns):
        table[column] = pd.to_numeric(table[column], errors="coerce")

    group = args.target_group.lower()
    fungi = ("ft", "fv", "fpro", "ffuj", "mp")
    decoys = args.specificity_groups or (["fv", "fpro", "ffuj"] if group == "ft" else [g for g in fungi if g != group])
    metrics = pd.DataFrame([calculate_metrics(r, group, decoys, fungi, args.rpm_basis)
                            for r in table.to_dict("records")])
    for column in metrics:
        table[column] = metrics[column].to_numpy()
    table["host_unique_pct"] = 100 * pd.to_numeric(table.get("sb_unique", 0), errors="coerce") / pd.to_numeric(table["all_primary"], errors="coerce").replace(0, np.nan)

    controls = table[table["pathogen"].str.lower().eq("pdb")]
    threshold = robust_threshold(controls["target_unique_rpm"])
    table["background_outlier"] = table["pathogen"].str.lower().eq("pdb") & table["target_unique_rpm"].gt(threshold)
    table["low_pathogen_signal"] = ~table["pathogen"].str.lower().eq("pdb") & table["target_unique_rpm"].le(threshold)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output_dir / "sample_mapping.tsv", sep="\t", index=False, lineterminator="\n", na_rep="NA")
    grouping = (
        table.groupby(["pathogen", "genotype", "water", "dai"], dropna=False)
        .agg(
            n=("sample_id", "size"),
            median_alignment=("overall_alignment_pct", "median"),
            median_host_unique_pct=("host_unique_pct", "median"),
            median_target_rpm=("target_unique_rpm", "median"),
            median_target_specificity=("target_specificity", "median"),
            median_target_breadth=("target_breadth_1x", "median"),
        ).reset_index()
    )
    grouping.to_csv(args.output_dir / "group_mapping.tsv", sep="\t", index=False, lineterminator="\n", na_rep="NA")

    inoculated = table[table["pathogen"].str.lower().eq(args.inoculated_label.lower())].copy()
    if not inoculated.empty and set(inoculated["genotype"]) >= {"wt", "bmr6", "bmr12"}:
        inoculated["log_target_rpm"] = np.log1p(inoculated["target_unique_rpm"])
        fit = smf.ols(
            'log_target_rpm ~ C(genotype, Treatment(reference="wt")) * C(water, Treatment(reference="dry"))',
            data=inoculated,
        ).fit(cov_type="HC3")
        ci = fit.conf_int()
        coefficients = pd.DataFrame({
            "term": fit.params.index,
            "estimate": fit.params.values,
            "se_hc3": fit.bse.values,
            "statistic": fit.tvalues.values,
            "p_hc3": fit.pvalues.values,
            "ci_low": ci[0].values,
            "ci_high": ci[1].values,
        })
        coefficients.to_csv(args.output_dir / "burden_model.tsv", sep="\t", index=False, lineterminator="\n")

        contrasts: list[dict[str, object]] = []
        for water in ("dry", "wet"):
            subset = inoculated[inoculated["water"].eq(water)]
            wt = subset.loc[subset["genotype"].eq("wt"), "target_unique_rpm"].to_numpy(float)
            for genotype in ("bmr6", "bmr12"):
                values = subset.loc[subset["genotype"].eq(genotype), "target_unique_rpm"].to_numpy(float)
                stat, p = mannwhitneyu(values, wt, alternative="two-sided")
                contrasts.append({
                    "contrast": f"{genotype}_vs_wt_{water}",
                    "median_first": np.median(values),
                    "median_second": np.median(wt),
                    "cliffs_delta": cliffs_delta(values, wt),
                    "u": stat,
                    "p": p,
                })
        for genotype in ("wt", "bmr6", "bmr12"):
            subset = inoculated[inoculated["genotype"].eq(genotype)]
            wet = subset.loc[subset["water"].eq("wet"), "target_unique_rpm"].to_numpy(float)
            dry = subset.loc[subset["water"].eq("dry"), "target_unique_rpm"].to_numpy(float)
            stat, p = mannwhitneyu(wet, dry, alternative="two-sided")
            contrasts.append({
                "contrast": f"wet_vs_dry_{genotype}",
                "median_first": np.median(wet),
                "median_second": np.median(dry),
                "cliffs_delta": cliffs_delta(wet, dry),
                "u": stat,
                "p": p,
            })
        contrast_table = pd.DataFrame(contrasts)
        contrast_table["fdr"] = bh(contrast_table["p"])
        contrast_table.to_csv(args.output_dir / "burden_contrasts.tsv", sep="\t", index=False, lineterminator="\n")

        trait_rows: list[dict[str, object]] = []
        traits = [
            "lesion_mm", "av_moist", "ja", "sa", "aba", "opda",
            "sol_syringic_acid", "sol_sinapic_acid", "wb_syringic_acid", "wb_ferulic_acid",
        ]
        for trait in traits:
            if trait not in inoculated:
                continue
            first = pd.to_numeric(inoculated["target_unique_rpm"], errors="coerce")
            second = pd.to_numeric(inoculated[trait], errors="coerce")
            keep = first.notna() & second.notna()
            if keep.sum() >= 4:
                rho, p = spearmanr(first[keep], second[keep])
                trait_rows.append({"trait": trait, "n": int(keep.sum()), "rho": rho, "p": p})
        pd.DataFrame(trait_rows).to_csv(args.output_dir / "burden_traits.tsv", sep="\t", index=False, lineterminator="\n", na_rep="NA")

    summary = {
        "libraries": len(table),
        "rpm_basis": args.rpm_basis,
        "specificity_comparators": list(decoys),
        "background_threshold_rpm": threshold,
        "median_inoculated_rpm": float(inoculated["target_unique_rpm"].median()) if not inoculated.empty else None,
        "median_control_rpm": float(controls["target_unique_rpm"].median()) if not controls.empty else None,
    }
    (args.output_dir / "mapping_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
