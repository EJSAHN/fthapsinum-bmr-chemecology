#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def read_de(path: Path, prefix: str) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t")
    rename = {column: f"{prefix}_{column}" for column in table.columns if column != "gene_id"}
    return table.rename(columns=rename)


def main() -> None:
    parser = argparse.ArgumentParser(description="Join fungal expression models and functional annotations.")
    parser.add_argument("--model-dir", type=Path, required=True)
    parser.add_argument("--annotation", type=Path, required=True)
    parser.add_argument("--gene-sets", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--min-abs-log2fc", type=float, default=1.0)
    parser.add_argument("--focus-module", default="aromatic_phenolic_metabolism")
    args = parser.parse_args()

    annotation = pd.read_csv(args.annotation, sep="\t", dtype=str, keep_default_na=False)
    sets = pd.read_csv(args.gene_sets, sep="\t", dtype=str, keep_default_na=False)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    combined_dir = args.output_dir / "contrasts"
    combined_dir.mkdir(parents=True, exist_ok=True)

    contrasts = pd.read_csv(args.model_dir / "contrast_summary.tsv", sep="\t")["contrast"].drop_duplicates().tolist()
    sensitivity_rows = []
    for contrast in contrasts:
        tables = [
            read_de(args.model_dir / "primary_state" / "de" / f"{contrast}.tsv", "state"),
            read_de(args.model_dir / "primary_tissue" / "de" / f"{contrast}.tsv", "tissue"),
            read_de(args.model_dir / "all_state" / "de" / f"{contrast}.tsv", "all"),
            read_de(args.model_dir / "primary_state_inclusive" / "de" / f"{contrast}.tsv", "inclusive"),
        ]
        merged = tables[0]
        for table in tables[1:]:
            merged = merged.merge(table, on="gene_id", how="outer")
        merged = merged.merge(annotation, on="gene_id", how="left")
        merged["state_supported"] = (
            pd.to_numeric(merged["state_FDR"], errors="coerce").lt(args.fdr)
            & pd.to_numeric(merged["state_logFC"], errors="coerce").abs().ge(args.min_abs_log2fc)
        )
        merged["tissue_supported"] = (
            pd.to_numeric(merged["tissue_FDR"], errors="coerce").lt(args.fdr)
            & pd.to_numeric(merged["tissue_logFC"], errors="coerce").abs().ge(args.min_abs_log2fc)
        )
        merged.to_csv(combined_dir / f"{contrast}.tsv", sep="\t", index=False, lineterminator="\n", na_rep="NA")
        sensitivity_rows.append({
            "contrast": contrast,
            "state_fdr05_lfc": int(merged["state_supported"].sum()),
            "tissue_fdr05_lfc": int(merged["tissue_supported"].sum()),
            "state_all_logfc_r": pd.to_numeric(merged["state_logFC"], errors="coerce").corr(pd.to_numeric(merged["all_logFC"], errors="coerce")),
            "state_inclusive_logfc_r": pd.to_numeric(merged["state_logFC"], errors="coerce").corr(pd.to_numeric(merged["inclusive_logFC"], errors="coerce")),
        })

    pd.DataFrame(sensitivity_rows).to_csv(args.output_dir / "contrast_sensitivity.tsv", sep="\t", index=False, lineterminator="\n", na_rep="NA")

    camera_tables = []
    for model in ("primary_state", "primary_tissue", "all_state", "primary_state_inclusive"):
        path = args.model_dir / model / "camera.tsv"
        table = pd.read_csv(path, sep="\t")
        table["model"] = model
        camera_tables.append(table)
    camera = pd.concat(camera_tables, ignore_index=True)
    camera.to_csv(args.output_dir / "camera.tsv", sep="\t", index=False, lineterminator="\n", na_rep="NA")
    camera[pd.to_numeric(camera["FDR"], errors="coerce").lt(args.fdr)].to_csv(
        args.output_dir / "camera_fdr05.tsv", sep="\t", index=False, lineterminator="\n", na_rep="NA"
    )

    focus_ids = sets.loc[sets["module"].eq(args.focus_module), "gene_id"].drop_duplicates().tolist()
    focus = annotation[annotation["gene_id"].isin(focus_ids)].copy()
    focus.to_csv(args.output_dir / "focus_genes.tsv", sep="\t", index=False, lineterminator="\n", na_rep="NA")


if __name__ == "__main__":
    main()
