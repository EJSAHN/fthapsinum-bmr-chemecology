#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import json
import shutil
import subprocess
import tempfile
import zipfile
from pathlib import Path


def fasta_bases(path: Path) -> int:
    opener = gzip.open if path.suffix == ".gz" else open
    total = 0
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(">"):
                total += len(line.strip())
    return total


def choose_genome(files: list[Path]) -> Path:
    candidates = [
        path for path in files
        if path.name.endswith("_genomic.fna")
        and not path.name.startswith(("cds_from_genomic", "rna_from_genomic"))
    ]
    if not candidates:
        candidates = [
            path for path in files
            if path.suffix in {".fna", ".fa", ".fasta"}
            and "cds" not in path.name.lower()
            and "rna" not in path.name.lower()
        ]
    if not candidates:
        raise RuntimeError("No genomic FASTA was found in the NCBI data package")
    return max(candidates, key=fasta_bases)


def choose_optional(files: list[Path], names: tuple[str, ...], suffixes: tuple[str, ...]) -> Path | None:
    preferred = [path for path in files if path.name in names]
    if preferred:
        return max(preferred, key=lambda path: path.stat().st_size)
    candidates = [path for path in files if path.name.endswith(suffixes)]
    return max(candidates, key=lambda path: path.stat().st_size) if candidates else None


def copy_text_or_gzip(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if source.suffix == ".gz":
        with gzip.open(source, "rb") as src, destination.open("wb") as dst:
            shutil.copyfileobj(src, dst)
    else:
        shutil.copy2(source, destination)


def main() -> None:
    parser = argparse.ArgumentParser(description="Download one NCBI genome data package.")
    parser.add_argument("--accession", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--require-annotation", type=int, choices=(0, 1), default=0)
    parser.add_argument("--require-protein", type=int, choices=(0, 1), default=0)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    genome_out = args.output_dir / "genome.fna"
    gff_out = args.output_dir / "annotation.gff3"
    protein_out = args.output_dir / "protein.faa"
    metadata_out = args.output_dir / "metadata.json"

    with tempfile.TemporaryDirectory(prefix=f"{args.label}.") as temp_name:
        temp = Path(temp_name)
        archive = temp / "ncbi_dataset.zip"
        command = [
            "datasets", "download", "genome", "accession", args.accession,
            "--include", "genome,gff3,protein", "--filename", str(archive),
        ]
        subprocess.run(command, check=True)
        with zipfile.ZipFile(archive) as bundle:
            bundle.extractall(temp / "package")
        files = [path for path in (temp / "package").rglob("*") if path.is_file()]

        genome = choose_genome(files)
        gff = choose_optional(files, ("genomic.gff", "genomic.gff3"), (".gff", ".gff3"))
        protein = choose_optional(files, ("protein.faa",), (".faa",))
        genome_bp = fasta_bases(genome)
        if genome_bp < 1_000_000:
            raise RuntimeError(f"Selected genome is unexpectedly small: {genome_bp} bp")
        if args.require_annotation and gff is None:
            raise RuntimeError(f"No annotation was available for {args.accession}")
        if args.require_protein and protein is None:
            raise RuntimeError(f"No protein FASTA was available for {args.accession}")

        copy_text_or_gzip(genome, genome_out)
        if gff is None:
            gff_out.write_text("", encoding="utf-8")
        else:
            copy_text_or_gzip(gff, gff_out)
        if protein is None:
            protein_out.write_text("", encoding="utf-8")
        else:
            copy_text_or_gzip(protein, protein_out)

        metadata = {
            "label": args.label,
            "accession": args.accession,
            "genome_source": str(genome.relative_to(temp / "package")),
            "genome_bases": genome_bp,
            "annotation_source": str(gff.relative_to(temp / "package")) if gff else None,
            "protein_source": str(protein.relative_to(temp / "package")) if protein else None,
        }
        metadata_out.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
