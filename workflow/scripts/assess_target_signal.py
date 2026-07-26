#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


def sample_name(column: str, suffix: str) -> str:
    name = Path(column).name
    stripped = re.sub(r"\.(?:FT|MP)\.mapq20\.primary\.bam$", "", name)
    if stripped != name:
        return stripped
    if suffix and name.endswith(suffix):
        return name[:-len(suffix)]
    return name


def main() -> None:
    parser = argparse.ArgumentParser(description="Assess whether pathogen-side RNA supports gene-level analysis.")
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--featurecounts-summary", type=Path, required=True)
    parser.add_argument("--inoculated-label", required=True)
    parser.add_argument("--control-label", default="pdb")
    parser.add_argument("--target-rpm-column", default="target_unique_rpm")
    parser.add_argument("--target-specificity-column", default="target_specificity")
    parser.add_argument("--bam-suffix", default=".bam")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--minimum-median-assigned", type=int, default=5000)
    parser.add_argument("--minimum-median-ratio", type=float, default=2.0)
    parser.add_argument("--minimum-specificity", type=float, default=0.5)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(args.manifest, sep="\t", dtype=str, keep_default_na=False)
    mapping = pd.read_csv(args.mapping, sep="\t")
    counts = pd.read_csv(args.counts, sep="\t", index_col=0).apply(pd.to_numeric, errors="coerce").fillna(0)
    summary = pd.read_csv(args.featurecounts_summary, sep="\t", index_col=0)

    assigned = pd.to_numeric(summary.loc["Assigned"], errors="coerce")
    assigned.index = [sample_name(column, args.bam_suffix) for column in assigned.index]
    genes_ge5 = (counts >= 5).sum(axis=0)
    genes_ge5.index = [sample_name(column, args.bam_suffix) for column in genes_ge5.index]

    samples = manifest.merge(mapping, on="sample_id", how="left", validate="one_to_one")
    samples["assigned_reads"] = samples["sample_id"].map(assigned).fillna(0).astype(int)
    samples["genes_ge5"] = samples["sample_id"].map(genes_ge5).fillna(0).astype(int)
    samples.to_csv(args.output_dir / "sample_signal.tsv", sep="\t", index=False, lineterminator="\n")

    inoculated = samples[samples["pathogen"].str.lower() == args.inoculated_label.lower()].copy()
    controls = samples[samples["pathogen"].str.lower() == args.control_label.lower()].copy()
    if inoculated.empty or controls.empty:
        raise RuntimeError("Both inoculated and control samples are required")

    rpm_inoculated = pd.to_numeric(inoculated[args.target_rpm_column], errors="coerce")
    rpm_control = pd.to_numeric(controls[args.target_rpm_column], errors="coerce")
    specificity = pd.to_numeric(inoculated[args.target_specificity_column], errors="coerce")
    assigned_inoculated = inoculated["assigned_reads"].astype(float)
    assigned_control = controls["assigned_reads"].astype(float)
    ratio = (float(rpm_inoculated.median()) + 1) / (float(rpm_control.median()) + 1)

    cell_support = (
        inoculated.assign(count_positive=inoculated["assigned_reads"] >= args.minimum_median_assigned)
        .groupby(["genotype", "water"], observed=True)["count_positive"]
        .sum()
        .rename("libraries_meeting_assignment_threshold")
        .reset_index()
    )
    cell_support.to_csv(args.output_dir / "cell_support.tsv", sep="\t", index=False, lineterminator="\n")

    sufficient = bool(
        float(assigned_inoculated.median()) >= args.minimum_median_assigned
        and ratio >= args.minimum_median_ratio
        and float(specificity.median()) >= args.minimum_specificity
        and (cell_support["libraries_meeting_assignment_threshold"] >= 2).all()
    )
    summary_payload = {
        "inoculated_samples": int(len(inoculated)),
        "control_samples": int(len(controls)),
        "median_target_rpm_inoculated": float(rpm_inoculated.median()),
        "median_target_rpm_control": float(rpm_control.median()),
        "inoculated_control_rpm_ratio": ratio,
        "median_target_specificity": float(specificity.median()),
        "median_assigned_reads_inoculated": float(assigned_inoculated.median()),
        "median_assigned_reads_control": float(assigned_control.median()),
        "median_genes_with_at_least_five_reads": float(inoculated["genes_ge5"].median()),
        "gene_level_analysis_supported": sufficient,
    }
    (args.output_dir / "target_signal.json").write_text(json.dumps(summary_payload, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
