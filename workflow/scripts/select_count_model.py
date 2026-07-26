#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import re
import shutil
from pathlib import Path

import numpy as np
import pandas as pd


def sample_name(column: str) -> str:
    name = Path(column).name
    name = re.sub(r"\.(?:FT|MP)\.mapq20\.primary\.bam$", "", name)
    return re.sub(r"\.bam$", "", name)


def read_counts(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", comment="#")
    if "Geneid" not in df.columns:
        raise ValueError(f"FeatureCounts output lacks Geneid: {path}")
    fixed = ["Geneid", "Chr", "Start", "End", "Strand", "Length"]
    sample_cols = [c for c in df.columns if c not in fixed]
    rename = {c: sample_name(c) for c in sample_cols}
    df = df.rename(columns=rename).set_index("Geneid")
    return df[[rename[c] for c in sample_cols]].apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)


def read_summary(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    first = df.columns[0]
    df = df.rename(columns={first: "Status"}).set_index("Status")
    df = df.rename(columns={c: sample_name(c) for c in df.columns})
    return df.apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)


def med(series: pd.Series) -> float:
    s = pd.to_numeric(series, errors="coerce").dropna()
    return float(s.median()) if len(s) else math.nan


def main() -> None:
    ap = argparse.ArgumentParser(description="Select the QuantSeq gene-count configuration.")
    ap.add_argument("--counts-dir", type=Path, required=True)
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--extension", type=int, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--inoculated-label", default="fus")
    ap.add_argument("--control-set-label", default="PDB_PRIMARY")
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(args.manifest, sep="\t", dtype=str).fillna("")
    manifest = manifest.set_index("sample_id", drop=False)
    combos: dict[str, dict[str, pd.DataFrame]] = {}
    metrics_rows: list[dict[str, object]] = []
    sample_rows: list[dict[str, object]] = []

    annotations = ["native", f"ext{args.extension}"]
    for annotation in annotations:
        for strand in (0, 1, 2):
            label = f"{annotation}_s{strand}"
            count_path = args.counts_dir / f"counts_{label}.txt"
            summary_path = args.counts_dir / f"counts_{label}.txt.summary"
            if not count_path.exists() or not summary_path.exists():
                raise SystemExit(f"Missing featureCounts output for {label}")
            counts = read_counts(count_path)
            summary = read_summary(summary_path)
            missing = sorted(set(manifest.index) - set(counts.columns))
            if missing:
                raise SystemExit(f"{label} is missing samples: {missing}")
            counts = counts.loc[:, manifest.index]
            summary = summary.loc[:, manifest.index]
            total = summary.sum(axis=0).replace(0, np.nan)
            assigned = summary.loc["Assigned"] if "Assigned" in summary.index else pd.Series(0, index=summary.columns)
            assigned_rate = assigned / total
            detected1 = (counts >= 1).sum(axis=0)
            detected3 = (counts >= 3).sum(axis=0)
            detected5 = (counts >= 5).sum(axis=0)
            detected10 = (counts >= 10).sum(axis=0)
            combos[label] = {"counts": counts, "summary": summary}

            sample_metrics = pd.DataFrame({
                "sample_id": manifest.index,
                "assigned_reads": assigned.reindex(manifest.index).values,
                "featurecounts_total_reads": total.reindex(manifest.index).values,
                "assigned_rate": assigned_rate.reindex(manifest.index).values,
                "detected_ge1": detected1.reindex(manifest.index).values,
                "detected_ge3": detected3.reindex(manifest.index).values,
                "detected_ge5": detected5.reindex(manifest.index).values,
                "detected_ge10": detected10.reindex(manifest.index).values,
            })
            sample_metrics["combo"] = label
            sample_metrics = sample_metrics.merge(manifest.reset_index(drop=True)[["sample_id", "pathogen", "analysis_set"]], on="sample_id", how="left")
            sample_rows.extend(sample_metrics.to_dict("records"))

            inoculated = sample_metrics["pathogen"].str.lower().eq(args.inoculated_label.lower())
            controls = sample_metrics["analysis_set"].eq(args.control_set_label)
            inoculated_assigned = med(sample_metrics.loc[inoculated, "assigned_reads"])
            control_assigned = med(sample_metrics.loc[controls, "assigned_reads"])
            metrics_rows.append({
                "combo": label,
                "annotation": annotation,
                "strand": strand,
                "median_inoculated_assigned_reads": inoculated_assigned,
                "median_control_assigned_reads": control_assigned,
                "median_inoculated_assigned_rate": med(sample_metrics.loc[inoculated, "assigned_rate"]),
                "median_control_assigned_rate": med(sample_metrics.loc[controls, "assigned_rate"]),
                "median_inoculated_detected_ge1": med(sample_metrics.loc[inoculated, "detected_ge1"]),
                "median_inoculated_detected_ge5": med(sample_metrics.loc[inoculated, "detected_ge5"]),
                "median_inoculated_detected_ge10": med(sample_metrics.loc[inoculated, "detected_ge10"]),
                "inoculated_to_control_assigned_ratio": (inoculated_assigned + 1.0) / (control_assigned + 1.0),
            })

    combo_metrics = pd.DataFrame(metrics_rows)
    all_sample_metrics = pd.DataFrame(sample_rows)
    combo_metrics.to_csv(args.output_dir / "featurecounts_combo_metrics.tsv", sep="\t", index=False, na_rep="NA")
    all_sample_metrics.to_csv(args.output_dir / "featurecounts_all_combo_sample_metrics.tsv", sep="\t", index=False, na_rep="NA")

    def choose_strand(annotation: str) -> tuple[int, str]:
        sub = combo_metrics[combo_metrics["annotation"].eq(annotation)].set_index("strand")
        s1 = float(sub.loc[1, "median_inoculated_assigned_rate"])
        s2 = float(sub.loc[2, "median_inoculated_assigned_rate"])
        s0 = float(sub.loc[0, "median_inoculated_assigned_rate"])
        best_stranded = 1 if s1 >= s2 else 2
        best_val = max(s1, s2)
        if s1 >= 0.90 * best_val:
            return 1, "strand_1_within_10pct_of_best_stranded"
        if best_val < 0.20 and s0 > 1.8 * max(best_val, 1e-9):
            return 0, "unstranded_rescue_required"
        return best_stranded, "best_stranded_assignment"

    native_strand, native_note = choose_strand("native")
    ext_name = f"ext{args.extension}"
    ext_strand, ext_note = choose_strand(ext_name)
    native_row = combo_metrics[combo_metrics["combo"].eq(f"native_s{native_strand}")].iloc[0]
    ext_row = combo_metrics[combo_metrics["combo"].eq(f"{ext_name}_s{ext_strand}")].iloc[0]

    native_rate = float(native_row["median_inoculated_assigned_rate"])
    ext_rate = float(ext_row["median_inoculated_assigned_rate"])
    native_det = float(native_row["median_inoculated_detected_ge5"])
    ext_det = float(ext_row["median_inoculated_detected_ge5"])
    native_ratio = float(native_row["inoculated_to_control_assigned_ratio"])
    ext_ratio = float(ext_row["inoculated_to_control_assigned_ratio"])
    extension_gain = (ext_rate >= native_rate * 1.05) or (ext_det >= native_det * 1.10)
    contrast_preserved = ext_ratio >= native_ratio * 0.75
    if extension_gain and contrast_preserved:
        selected_annotation, selected_strand = ext_name, ext_strand
        annotation_note = "QuantSeq extension increased pathogen assignment or detection without collapsing inoculated-to-control contrast"
    else:
        selected_annotation, selected_strand = "native", native_strand
        annotation_note = "Native models retained because extension gain was small or background contrast declined"
    selected_label = f"{selected_annotation}_s{selected_strand}"

    selected_counts = combos[selected_label]["counts"]
    selected_summary = combos[selected_label]["summary"]
    selected_counts.to_csv(args.output_dir / "selected_counts.tsv", sep="\t", index=True, index_label="gene_id")
    selected_summary.to_csv(args.output_dir / "selected_featurecounts_summary.tsv", sep="\t", index=True, index_label="Status")
    selected_sample = all_sample_metrics[all_sample_metrics["combo"].eq(selected_label)].copy()
    selected_sample.to_csv(args.output_dir / "selected_sample_metrics.tsv", sep="\t", index=False, na_rep="NA")

    selection = {
        "selected_combo": selected_label,
        "selected_annotation": selected_annotation,
        "selected_featurecounts_strand": selected_strand,
        "extension_bp": args.extension,
        "native_strand_note": native_note,
        "extended_strand_note": ext_note,
        "annotation_selection_note": annotation_note,
        "native_candidate": native_row.to_dict(),
        "extended_candidate": ext_row.to_dict(),
        "orientation_warning": selected_strand != 1,
    }
    (args.output_dir / "count_selection.json").write_text(json.dumps(selection, indent=2), encoding="utf-8")
    print(args.output_dir / "count_selection.json")
    print(f"Selected featureCounts combination: {selected_label}")


if __name__ == "__main__":
    main()
