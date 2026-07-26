#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


SAMPLE_COLUMNS = {
    "sample_id", "canonical_sample_id", "genotype", "pathogen", "water", "dai",
    "replicate", "sequencing_run", "run_accession", "fastq_url", "fastq_md5",
    "expected_reads", "expected_bytes",
}
REFERENCE_COLUMNS = {
    "label", "accession", "species_group", "include_in_competitive",
    "require_annotation", "require_protein",
}
SET_COLUMNS = {"set_name", "sample_id"}


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate sample, analysis-set, and reference manifests.")
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--analysis-sets", type=Path, required=True)
    parser.add_argument("--references", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    samples = pd.read_csv(args.samples, sep="\t", dtype=str, keep_default_na=False)
    sets = pd.read_csv(args.analysis_sets, sep="\t", dtype=str, keep_default_na=False)
    references = pd.read_csv(args.references, sep="\t", dtype=str, keep_default_na=False)

    missing_samples = sorted(SAMPLE_COLUMNS - set(samples.columns))
    missing_sets = sorted(SET_COLUMNS - set(sets.columns))
    missing_references = sorted(REFERENCE_COLUMNS - set(references.columns))
    if missing_samples or missing_sets or missing_references:
        raise SystemExit(json.dumps({
            "missing_sample_columns": missing_samples,
            "missing_analysis_set_columns": missing_sets,
            "missing_reference_columns": missing_references,
        }, indent=2))

    if samples["sample_id"].duplicated().any():
        raise SystemExit("Sample IDs are not unique")
    if samples["run_accession"].duplicated().any():
        raise SystemExit("Run accessions are not unique")
    if references["label"].duplicated().any():
        raise SystemExit("Reference labels are not unique")
    if sets.duplicated(["set_name", "sample_id"]).any():
        raise SystemExit("Analysis-set assignments are not unique")

    missing_set_samples = sorted(set(sets["sample_id"]) - set(samples["sample_id"]))
    if missing_set_samples:
        raise SystemExit(f"Analysis sets contain unknown samples: {missing_set_samples}")

    for column in ("dai", "sequencing_run", "expected_reads", "expected_bytes"):
        values = pd.to_numeric(samples[column], errors="coerce")
        if values.isna().any():
            bad = samples.loc[values.isna(), "sample_id"].tolist()
            raise SystemExit(f"Column {column} is missing or nonnumeric for: {bad}")
    for column in ("include_in_competitive", "require_annotation", "require_protein"):
        values = pd.to_numeric(references[column], errors="coerce")
        if values.isna().any() or not values.isin([0, 1]).all():
            raise SystemExit(f"Reference column {column} must contain only 0 or 1")

    if not samples["fastq_url"].str.startswith("https://").all():
        raise SystemExit("FASTQ URLs must use HTTPS")
    if not samples["fastq_md5"].str.fullmatch(r"[0-9a-fA-F]{32}").all():
        raise SystemExit("FASTQ MD5 values are malformed")

    set_sizes = sets.groupby("set_name").size().sort_index().astype(int).to_dict()
    payload = {
        "samples": int(len(samples)),
        "references": int(len(references)),
        "analysis_sets": set_sizes,
        "competitive_references": references.loc[
            pd.to_numeric(references["include_in_competitive"]) == 1, "label"
        ].tolist(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
