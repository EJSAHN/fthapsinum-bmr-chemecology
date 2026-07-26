#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description="Count genes with a previously selected annotation and strand setting.")
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--bam-dir", type=Path, required=True)
    parser.add_argument("--gene-model-dir", type=Path, required=True)
    parser.add_argument("--model-prefix", required=True)
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    selection = json.loads(args.selection.read_text(encoding="utf-8"))
    annotation = str(selection["selected_annotation"])
    strand = int(selection["selected_featurecounts_strand"])
    saf = args.gene_model_dir / f"{args.model_prefix}_{annotation}.saf"
    if not saf.exists():
        raise RuntimeError(f"Missing gene model: {saf}")

    manifest = pd.read_csv(args.manifest, sep="\t", dtype=str, keep_default_na=False)
    bams = [args.bam_dir / f"{sample}.bam" for sample in manifest["sample_id"]]
    missing = [str(path) for path in bams if not path.exists() or path.stat().st_size == 0]
    if missing:
        raise RuntimeError(f"Missing BAM files: {missing}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    command = [
        "featureCounts", "-T", str(args.threads), "-F", "SAF", "-a", str(saf),
        "-o", str(args.output), "-s", str(strand), "--primary",
        *map(str, bams),
    ]
    subprocess.run(command, check=True)


if __name__ == "__main__":
    main()
