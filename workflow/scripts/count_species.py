#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser(description="Count primary competitive alignments by species group.")
    ap.add_argument("--bam", type=Path, required=True)
    ap.add_argument("--contig-map", type=Path, required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--min-mapq", type=int, default=20)
    args = ap.parse_args()

    import pysam

    contig_group: dict[str, str] = {}
    groups: set[str] = set()
    with args.contig_map.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            contig_group[row["contig"]] = row["species_group"]
            groups.add(row["species_group"])

    counts: Counter[str] = Counter()
    mapped_any: Counter[str] = Counter()
    qpass: Counter[str] = Counter()
    unique: Counter[str] = Counter()
    with pysam.AlignmentFile(str(args.bam), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_secondary:
                counts["secondary_alignments"] += 1
                continue
            if read.is_supplementary:
                counts["supplementary_alignments"] += 1
                continue
            counts["all_primary"] += 1
            if read.is_duplicate:
                counts["duplicate_primary"] += 1
            if read.is_unmapped:
                counts["unmapped_primary"] += 1
                continue
            counts["mapped_primary"] += 1
            ref = bam.get_reference_name(read.reference_id)
            group = contig_group.get(ref, "OTHER")
            mapped_any[group] += 1
            if read.mapping_quality >= args.min_mapq:
                counts["mapq_pass"] += 1
                qpass[group] += 1
                nh = read.get_tag("NH") if read.has_tag("NH") else 1
                if nh == 1:
                    counts["unique_pass"] += 1
                    unique[group] += 1

    row: dict[str, int | str] = {
        "sample_id": args.sample,
        "all_primary": counts["all_primary"],
        "unmapped_primary": counts["unmapped_primary"],
        "mapped_primary": counts["mapped_primary"],
        "mapq_pass": counts["mapq_pass"],
        "unique_pass": counts["unique_pass"],
        "secondary_alignments": counts["secondary_alignments"],
        "supplementary_alignments": counts["supplementary_alignments"],
        "duplicate_primary": counts["duplicate_primary"],
        "min_mapq": args.min_mapq,
    }
    for group in sorted(groups | {"OTHER"}):
        key = group.lower()
        row[f"{key}_mapped"] = mapped_any[group]
        row[f"{key}_q20"] = qpass[group]
        row[f"{key}_unique"] = unique[group]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)
    print(args.output)


if __name__ == "__main__":
    main()
