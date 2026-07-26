#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_assignment(value: str) -> tuple[str, str]:
    if ":" in value:
        set_name, label = value.split(":", 1)
    else:
        set_name, label = value, value
    if not set_name or not label:
        raise argparse.ArgumentTypeError("Set assignments must be SET or SET:LABEL")
    return set_name, label


def main() -> None:
    parser = argparse.ArgumentParser(description="Create an ordered sample manifest from named analysis sets.")
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--analysis-sets", type=Path, required=True)
    parser.add_argument("--set", dest="assignments", action="append", type=parse_assignment, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    samples = pd.read_csv(args.samples, sep="\t", dtype=str, keep_default_na=False)
    sets = pd.read_csv(args.analysis_sets, sep="\t", dtype=str, keep_default_na=False)
    sample_index = samples.set_index("sample_id", drop=False)
    rows: list[pd.DataFrame] = []
    seen: set[str] = set()
    for set_name, label in args.assignments:
        members = sets.loc[sets["set_name"].eq(set_name), "sample_id"].tolist()
        if not members:
            raise RuntimeError(f"Analysis set is empty or absent: {set_name}")
        overlap = sorted(seen.intersection(members))
        if overlap:
            raise RuntimeError(f"Samples occur in more than one selected set: {overlap}")
        seen.update(members)
        missing = sorted(set(members) - set(sample_index.index))
        if missing:
            raise RuntimeError(f"Unknown samples in {set_name}: {missing}")
        table = sample_index.loc[members].copy()
        table["analysis_set"] = label
        rows.append(table)
    output = pd.concat(rows, ignore_index=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(args.output, sep="\t", index=False, lineterminator="\n")


if __name__ == "__main__":
    main()
