#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import pandas as pd


def clean_id(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def read_fasta(path: Path):
    name = None
    description = ""
    chunks: list[str] = []
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, description, "".join(chunks)
                header = line[1:].strip()
                fields = header.split(maxsplit=1)
                name = fields[0]
                description = fields[1] if len(fields) > 1 else ""
                chunks = []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        yield name, description, "".join(chunks)


def write_record(handle, name: str, description: str, sequence: str) -> None:
    handle.write(f">{name}{(' ' + description) if description else ''}\n")
    for start in range(0, len(sequence), 80):
        handle.write(sequence[start:start + 80] + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a prefixed host-pathogen competitive reference.")
    parser.add_argument("--references", type=Path, required=True)
    parser.add_argument("--reference-root", type=Path, required=True)
    parser.add_argument("--output-fasta", type=Path, required=True)
    parser.add_argument("--contig-map", type=Path, required=True)
    parser.add_argument("--ft-bed", type=Path, required=True)
    parser.add_argument("--mp-bed", type=Path, required=True)
    args = parser.parse_args()

    references = pd.read_csv(args.references, sep="\t", dtype=str, keep_default_na=False)
    selected = references[pd.to_numeric(references["include_in_competitive"]).eq(1)].copy()
    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    args.contig_map.parent.mkdir(parents=True, exist_ok=True)

    seen: set[str] = set()
    map_rows: list[dict[str, object]] = []
    bed_rows: dict[str, list[tuple[str, int]]] = {"FT": [], "MP": []}
    with args.output_fasta.open("wt", encoding="utf-8", newline="\n") as fasta:
        for row in selected.itertuples(index=False):
            source = args.reference_root / row.label / "genome.fna"
            if not source.exists() or source.stat().st_size == 0:
                raise RuntimeError(f"Missing reference genome: {source}")
            for original, description, sequence in read_fasta(source):
                contig = f"{row.label}__{clean_id(original)}"
                if contig in seen:
                    raise RuntimeError(f"Duplicate prefixed contig: {contig}")
                seen.add(contig)
                write_record(fasta, contig, description, sequence)
                map_rows.append({
                    "contig": contig,
                    "species_group": row.species_group,
                    "reference_label": row.label,
                    "source_contig": original,
                    "length": len(sequence),
                })
                if row.species_group in bed_rows:
                    bed_rows[row.species_group].append((contig, len(sequence)))

    with args.contig_map.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["contig", "species_group", "reference_label", "source_contig", "length"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(map_rows)

    for group, path in (("FT", args.ft_bed), ("MP", args.mp_bed)):
        with path.open("wt", encoding="utf-8", newline="\n") as handle:
            for contig, length in bed_rows[group]:
                handle.write(f"{contig}\t0\t{length}\n")
        if not bed_rows[group]:
            raise RuntimeError(f"No {group} contigs were included in the competitive reference")


if __name__ == "__main__":
    main()
