#!/usr/bin/env python3
from __future__ import annotations

import argparse
import bisect
import csv
import json
import re
import shutil
import tempfile
import urllib.parse
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path


def safe(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def parse_attrs(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    text = text.strip().rstrip(";")
    if not text or text == ".":
        return out
    for piece in text.split(";"):
        piece = piece.strip()
        if not piece:
            continue
        if "=" in piece:
            key, value = piece.split("=", 1)
        elif " " in piece:
            key, value = piece.split(" ", 1)
            value = value.strip().strip('"')
        else:
            continue
        out[key.strip()] = urllib.parse.unquote(value.strip().strip('"'))
    return out


@dataclass
class Feature:
    seqid: str
    kind: str
    start: int
    end: int
    strand: str
    attrs: dict[str, str]


@dataclass
class Gene:
    gene_id: str
    original_id: str
    contig: str
    start: int
    end: int
    strand: str
    attrs: dict[str, str] = field(default_factory=dict)
    source_features: set[str] = field(default_factory=set)

    def absorb(self, start: int, end: int, strand: str, attrs: dict[str, str], kind: str) -> None:
        self.start = min(self.start, start)
        self.end = max(self.end, end)
        if self.strand not in {"+", "-"} and strand in {"+", "-"}:
            self.strand = strand
        for key in ("Name", "gene", "locus_tag", "product", "protein_id", "description", "Target", "Dbxref"):
            if key in attrs and attrs[key] and key not in self.attrs:
                self.attrs[key] = attrs[key]
        self.source_features.add(kind)


def read_fai(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            lengths[f[0]] = int(f[1])
    return lengths


def parse_gff(path: Path, prefix: str, already_prefixed: bool, fai: dict[str, int]) -> tuple[list[Gene], dict[str, int]]:
    features: list[Feature] = []
    kind_counts: defaultdict[str, int] = defaultdict(int)
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                if line.startswith("##FASTA"):
                    break
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            seqid, _, kind, start, end, _, strand, _, attr_text = fields
            if strand not in {"+", "-", ".", "?"}:
                continue
            try:
                s, e = int(start), int(end)
            except ValueError:
                continue
            contig = seqid if already_prefixed else f"{prefix}{safe(seqid)}"
            if not contig.startswith(prefix):
                continue
            if contig not in fai:
                continue
            attrs = parse_attrs(attr_text)
            features.append(Feature(contig, kind.lower(), s, e, strand, attrs))
            kind_counts[kind.lower()] += 1

    transcript_to_gene: dict[str, str] = {}
    gene_id_by_raw: dict[str, str] = {}
    genes: dict[str, Gene] = {}

    def normalized_gene_id(raw: str) -> str:
        raw = raw.strip().split(",")[0]
        if raw in gene_id_by_raw:
            return gene_id_by_raw[raw]
        gid = raw if already_prefixed and raw.startswith(prefix.rstrip("_")) else f"{prefix}{safe(raw)}"
        base = gid
        i = 2
        while gid in genes and genes[gid].original_id != raw:
            gid = f"{base}_{i}"
            i += 1
        gene_id_by_raw[raw] = gid
        return gid

    gene_types = {"gene", "pseudogene"}
    transcript_types = {"mrna", "transcript", "ncrna", "trna", "rrna", "lnc_rna", "primary_transcript"}
    child_types = {"exon", "cds", "three_prime_utr", "five_prime_utr", "utr"}

    for f in features:
        if f.kind in gene_types:
            raw = f.attrs.get("ID") or f.attrs.get("gene_id") or f.attrs.get("locus_tag") or f.attrs.get("Name")
            if not raw:
                raw = f"gene_{f.seqid}_{f.start}_{f.end}_{f.strand}"
            gid = normalized_gene_id(raw)
            g = genes.get(gid)
            if g is None:
                g = Gene(gid, raw, f.seqid, f.start, f.end, f.strand, dict(f.attrs), {f.kind})
                genes[gid] = g
            else:
                g.absorb(f.start, f.end, f.strand, f.attrs, f.kind)

    for f in features:
        if f.kind in transcript_types:
            tid = f.attrs.get("ID") or f.attrs.get("transcript_id")
            parent = f.attrs.get("Parent") or f.attrs.get("gene_id") or f.attrs.get("locus_tag")
            if not tid:
                tid = f"tx_{f.seqid}_{f.start}_{f.end}_{f.strand}"
            if parent:
                raw_gene = parent.split(",")[0]
            else:
                raw_gene = tid
            gid = normalized_gene_id(raw_gene)
            transcript_to_gene[tid] = gid
            g = genes.get(gid)
            if g is None:
                g = Gene(gid, raw_gene, f.seqid, f.start, f.end, f.strand, dict(f.attrs), {f.kind})
                genes[gid] = g
            else:
                g.absorb(f.start, f.end, f.strand, f.attrs, f.kind)

    orphan_counter = 0
    for f in features:
        if f.kind not in child_types:
            continue
        parents = (f.attrs.get("Parent") or f.attrs.get("transcript_id") or f.attrs.get("gene_id") or "").split(",")
        parents = [p for p in parents if p]
        if not parents:
            orphan_counter += 1
            parents = [f"orphan_{f.seqid}_{f.start}_{f.end}_{orphan_counter}"]
        for parent in parents:
            gid = transcript_to_gene.get(parent)
            if gid is None:
                gid = normalized_gene_id(parent)
            g = genes.get(gid)
            if g is None:
                g = Gene(gid, parent, f.seqid, f.start, f.end, f.strand, dict(f.attrs), {f.kind})
                genes[gid] = g
            else:
                if g.contig != f.seqid:
                    continue
                g.absorb(f.start, f.end, f.strand, f.attrs, f.kind)

    result = [
        g for g in genes.values()
        if g.contig in fai and g.start >= 1 and g.end >= g.start and g.strand in {"+", "-"}
    ]
    result.sort(key=lambda g: (g.contig, g.start, g.end, g.gene_id))
    return result, dict(kind_counts)


def extend_genes(genes: list[Gene], fai: dict[str, int], extension: int) -> dict[str, tuple[int, int]]:
    by_contig: defaultdict[str, list[Gene]] = defaultdict(list)
    for g in genes:
        by_contig[g.contig].append(g)
    ext: dict[str, tuple[int, int]] = {}
    for contig, items in by_contig.items():
        items.sort(key=lambda g: (g.start, g.end))
        starts = sorted(g.start for g in items)
        ends = sorted(g.end for g in items)
        length = fai[contig]
        for g in items:
            s, e = g.start, g.end
            if extension > 0 and g.strand == "+":
                idx = bisect.bisect_right(starts, e)
                limit = starts[idx] - 1 if idx < len(starts) else length
                e = min(length, e + extension, limit)
            elif extension > 0 and g.strand == "-":
                idx = bisect.bisect_left(ends, s) - 1
                limit = ends[idx] + 1 if idx >= 0 else 1
                s = max(1, s - extension, limit)
            ext[g.gene_id] = (s, e)
    return ext


def write_outputs(genes: list[Gene], kind_counts: dict[str, int], fai: dict[str, int], extension: int, gff: Path, outdir: Path, output_prefix: str) -> None:
    if len(genes) < 1000:
        raise SystemExit(f"Only {len(genes)} usable gene models were parsed from {gff}; at least 1000 are required")
    ext = extend_genes(genes, fai, extension)
    outdir.parent.mkdir(parents=True, exist_ok=True)
    temp = Path(tempfile.mkdtemp(prefix="gene_models_", dir=str(outdir.parent)))
    try:
        native = temp / f"{output_prefix}_native.saf"
        extended = temp / f"{output_prefix}_ext{extension}.saf"
        metadata = temp / f"{output_prefix}_metadata.tsv"
        metrics = temp / "gene_models_metrics.json"
        for path, use_ext in ((native, False), (extended, True)):
            with path.open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["GeneID", "Chr", "Start", "End", "Strand"])
                for g in genes:
                    s, e = ext[g.gene_id] if use_ext else (g.start, g.end)
                    writer.writerow([g.gene_id, g.contig, s, e, g.strand])
        fields = [
            "gene_id", "original_id", "contig", "start", "end", "strand", "native_length",
            "extended_start", "extended_end", "extended_length", "extension_bp", "source_features",
            "Name", "gene", "locus_tag", "product", "protein_id", "description", "Target", "Dbxref",
        ]
        with metadata.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for g in genes:
                es, ee = ext[g.gene_id]
                row = {
                    "gene_id": g.gene_id,
                    "original_id": g.original_id,
                    "contig": g.contig,
                    "start": g.start,
                    "end": g.end,
                    "strand": g.strand,
                    "native_length": g.end - g.start + 1,
                    "extended_start": es,
                    "extended_end": ee,
                    "extended_length": ee - es + 1,
                    "extension_bp": extension,
                    "source_features": ",".join(sorted(g.source_features)),
                }
                for key in fields[12:]:
                    row[key] = g.attrs.get(key, "")
                writer.writerow(row)
        metrics.write_text(json.dumps({
            "source_gff": str(gff),
            "gene_count": len(genes),
            "contig_count": len({g.contig for g in genes}),
            "extension_bp": extension,
            "feature_type_counts": kind_counts,
            "native_total_bp": sum(g.end - g.start + 1 for g in genes),
            "extended_total_bp": sum(ext[g.gene_id][1] - ext[g.gene_id][0] + 1 for g in genes),
        }, indent=2), encoding="utf-8")
        if outdir.exists():
            shutil.rmtree(outdir)
        temp.rename(outdir)
    except Exception:
        shutil.rmtree(temp, ignore_errors=True)
        raise
    print(outdir / "gene_models_metrics.json")


def main() -> None:
    ap = argparse.ArgumentParser(description="Build native and QuantSeq-aware SAF gene models.")
    ap.add_argument("--gff", type=Path, required=True)
    ap.add_argument("--fai", type=Path, required=True)
    ap.add_argument("--prefix", default="FTH_FL4__")
    ap.add_argument("--already-prefixed", type=int, choices=(0, 1), default=0)
    ap.add_argument("--extension", type=int, default=1000)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--output-prefix", default="target_genes")
    args = ap.parse_args()
    if not args.gff.exists() or args.gff.stat().st_size == 0:
        raise SystemExit(f"GFF is missing or empty: {args.gff}")
    fai = read_fai(args.fai)
    genes, kind_counts = parse_gff(args.gff, args.prefix, bool(args.already_prefixed), fai)
    write_outputs(genes, kind_counts, fai, args.extension, args.gff, args.output_dir, args.output_prefix)


if __name__ == "__main__":
    main()
