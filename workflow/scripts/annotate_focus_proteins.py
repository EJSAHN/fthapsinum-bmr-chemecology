#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import os
import re
import shlex
import subprocess
from pathlib import Path
from typing import Iterable

import pandas as pd


def require(path: Path) -> Path:
    if not path.exists() or path.stat().st_size == 0:
        raise SystemExit(f"Required file missing or empty: {path}")
    return path


def run(cmd: list[str], stdout: Path | None = None) -> None:
    print("RUN:", " ".join(shlex.quote(x) for x in cmd), flush=True)
    if stdout is None:
        subprocess.run(cmd, check=True)
    else:
        stdout.parent.mkdir(parents=True, exist_ok=True)
        with stdout.open("w", encoding="utf-8", newline="") as handle:
            subprocess.run(cmd, check=True, stdout=handle)


def read_fasta(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    seqs: dict[str, str] = {}
    desc: dict[str, str] = {}
    current = ""
    chunks: list[str] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                if current:
                    seqs[current] = "".join(chunks)
                header = line[1:].strip()
                parts = header.split(maxsplit=1)
                current = parts[0]
                desc[current] = parts[1] if len(parts) > 1 else ""
                chunks = []
            elif current:
                chunks.append(line.strip())
    if current:
        seqs[current] = "".join(chunks)
    return seqs, desc


def write_fasta(records: Iterable[tuple[str, str, str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as out:
        for identifier, description, sequence in records:
            out.write(f">{identifier} {description}\n")
            for i in range(0, len(sequence), 60):
                out.write(sequence[i : i + 60] + "\n")


def parse_blast(path: Path, label: str) -> pd.DataFrame:
    cols = [
        "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend",
        "sstart", "send", "evalue", "bitscore", "qlen", "slen",
    ]
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame(columns=cols + ["target_label", "qcov", "scov"])
    df = pd.read_csv(path, sep="\t", names=cols)
    df["qcov"] = 100.0 * df["length"] / df["qlen"].clip(lower=1)
    df["scov"] = 100.0 * df["length"] / df["slen"].clip(lower=1)
    df["target_label"] = label
    return df


def parse_miniprot(path: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if not path.exists():
        return pd.DataFrame()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("##PAF"):
                continue
            text = line.rstrip("\n")
            fields = text.split("\t")
            if fields[0] == "##PAF":
                f = fields[1:]
            else:
                f = text.removeprefix("##PAF").lstrip().split("\t")
            if len(f) < 12:
                continue
            try:
                qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend, nmatch, blen, mapq = f[:12]
                qlen_i, qs_i, qe_i = int(qlen), int(qstart), int(qend)
                nm_i, bl_i = int(nmatch), int(blen)
                tags = {x.split(":", 2)[0]: x.split(":", 2)[-1] for x in f[12:] if x.count(":") >= 2}
                rows.append(
                    {
                        "gene_id": qname,
                        "query_length_aa": qlen_i,
                        "query_start": qs_i,
                        "query_end": qe_i,
                        "query_coverage_pct": 100.0 * (qe_i - qs_i) / max(1, qlen_i),
                        "target_contig": tname,
                        "target_start": int(tstart),
                        "target_end": int(tend),
                        "strand": strand,
                        "aligned_matches": nm_i,
                        "alignment_block": bl_i,
                        "alignment_identity_pct": 100.0 * nm_i / max(1, bl_i),
                        "mapq": int(mapq),
                        "cigar": tags.get("cg", ""),
                        "cs": tags.get("cs", ""),
                    }
                )
            except (ValueError, IndexError):
                continue
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    df = df.sort_values(
        ["gene_id", "query_coverage_pct", "alignment_identity_pct", "mapq"],
        ascending=[True, False, False, False],
    ).drop_duplicates("gene_id", keep="first")
    return df


def parse_pfam(path: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split(maxsplit=22)
            if len(f) < 23:
                continue
            try:
                rows.append(
                    {
                        "pfam_name": f[0],
                        "pfam_accession": f[1],
                        "pfam_length": int(f[2]),
                        "gene_id": f[3],
                        "protein_length": int(f[5]),
                        "full_evalue": float(f[6]),
                        "full_score": float(f[7]),
                        "domain_i_evalue": float(f[12]),
                        "domain_score": float(f[13]),
                        "hmm_from": int(f[15]),
                        "hmm_to": int(f[16]),
                        "ali_from": int(f[17]),
                        "ali_to": int(f[18]),
                        "description": f[22],
                    }
                )
            except ValueError:
                continue
    return pd.DataFrame(rows)



def main() -> None:
    parser = argparse.ArgumentParser(description="Validate projected focus proteins.")
    parser.add_argument("--focus-genes", type=Path, required=True)
    parser.add_argument("--source-proteins", type=Path, required=True)
    parser.add_argument("--fl4-genome", type=Path, required=True)
    parser.add_argument("--comparison", action="append", default=[], help="LABEL=PROTEIN_FASTA")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--pfam-hmm", type=Path)
    parser.add_argument("--threads", type=int, default=8)
    args = parser.parse_args()

    focus = pd.read_csv(require(args.focus_genes), sep="\t", dtype=str).fillna("")
    source_faa = require(args.source_proteins)
    fl4_genome = require(args.fl4_genome)
    out = args.output_dir
    work = args.work_dir
    out.mkdir(parents=True, exist_ok=True)
    work.mkdir(parents=True, exist_ok=True)

    seqs, descriptions = read_fasta(source_faa)
    missing = sorted(set(focus["gene_id"]) - set(seqs))
    if missing:
        raise SystemExit(f"Focus proteins missing from source FASTA: {missing}")
    focus_faa = out / "focus_proteins.faa"
    write_fasta(
        ((row.gene_id, f"product={getattr(row, 'product', '')}", seqs[row.gene_id]) for row in focus.itertuples(index=False)),
        focus_faa,
    )

    projection_file = out / "focus_projection.gff3"
    run(["miniprot", "-t", str(args.threads), "-I", "--gff", "--aln", str(fl4_genome), str(focus_faa)], projection_file)
    projection = parse_miniprot(projection_file)
    projection.to_csv(out / "projection_metrics.tsv", sep="\t", index=False, lineterminator="\n")

    pfam = pd.DataFrame()
    if args.pfam_hmm:
        require(args.pfam_hmm)
        domtbl = out / "pfam.domtblout"
        run([
            "hmmscan", "--cpu", str(args.threads), "--cut_ga", "--noali", "--domtblout", str(domtbl),
            str(args.pfam_hmm), str(focus_faa),
        ], out / "pfam.txt")
        pfam = parse_pfam(domtbl)
    pfam.to_csv(out / "pfam_domains.tsv", sep="\t", index=False, lineterminator="\n")

    comparison_paths = {}
    for item in args.comparison:
        if "=" not in item:
            raise SystemExit(f"Invalid --comparison value: {item}")
        label, path = item.split("=", 1)
        comparison_paths[label] = require(Path(path))

    hit_tables = []
    db_dir = work / "diamond"
    db_dir.mkdir(parents=True, exist_ok=True)
    for label, proteome in comparison_paths.items():
        db = db_dir / label
        run(["diamond", "makedb", "--in", str(proteome), "--db", str(db), "--threads", str(args.threads)])
        hits_path = out / f"hits_{label}.tsv"
        run([
            "diamond", "blastp", "--query", str(focus_faa), "--db", str(db), "--out", str(hits_path),
            "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen",
            "--max-target-seqs", "10", "--evalue", "1e-8", "--sensitive", "--threads", str(args.threads),
        ])
        hit_tables.append(parse_blast(hits_path, label))
    hits = pd.concat(hit_tables, ignore_index=True) if hit_tables else pd.DataFrame()
    hits.to_csv(out / "cross_species_hits.tsv", sep="\t", index=False, lineterminator="\n")
    best = (
        hits.sort_values(["qseqid", "target_label", "bitscore"], ascending=[True, True, False])
        .drop_duplicates(["qseqid", "target_label"], keep="first")
        if not hits.empty else pd.DataFrame()
    )
    best.to_csv(out / "cross_species_best_hits.tsv", sep="\t", index=False, lineterminator="\n")

    projection_map = projection.set_index("gene_id").to_dict("index") if not projection.empty else {}
    pfam_map = {}
    if not pfam.empty:
        for gene, subset in pfam.sort_values("domain_i_evalue").groupby("gene_id"):
            pfam_map[str(gene)] = "; ".join(f"{row.pfam_accession}:{row.pfam_name}" for row in subset.itertuples(index=False))
    best_groups = {gene: subset for gene, subset in best.groupby("qseqid")} if not best.empty else {}

    rows = []
    for row in focus.itertuples(index=False):
        projected = projection_map.get(row.gene_id, {})
        homologs = best_groups.get(row.gene_id, pd.DataFrame())
        strong = int(((homologs["pident"] >= 50) & (homologs["qcov"] >= 70)).sum()) if not homologs.empty else 0
        qcov = float(projected.get("query_coverage_pct", 0) or 0)
        identity = float(projected.get("alignment_identity_pct", 0) or 0)
        domains = pfam_map.get(row.gene_id, "")
        if qcov >= 80 and identity >= 50 and strong >= max(1, len(comparison_paths) - 1) and (domains or not args.pfam_hmm):
            confidence = "high"
        elif qcov >= 80 and identity >= 50 and strong >= 1:
            confidence = "moderate"
        else:
            confidence = "low"
        rows.append({
            "gene_id": row.gene_id,
            "product": getattr(row, "product", ""),
            "source_length_aa": len(seqs[row.gene_id]),
            "projection_query_coverage_pct": qcov,
            "projection_identity_pct": identity,
            "projection_contig": projected.get("target_contig", getattr(row, "contig", "")),
            "projection_start": projected.get("target_start", getattr(row, "start", "")),
            "projection_end": projected.get("target_end", getattr(row, "end", "")),
            "pfam_domains": domains,
            "strong_cross_species_hits": strong,
            "annotation_confidence": confidence,
        })
    pd.DataFrame(rows).to_csv(out / "focus_annotation.tsv", sep="\t", index=False, lineterminator="\n")


if __name__ == "__main__":
    main()
