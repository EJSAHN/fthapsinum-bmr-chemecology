#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
import shlex
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


def revcomp(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def canonical_kmer(seq: str) -> str:
    seq = seq.upper()
    rc = revcomp(seq)
    return seq if seq <= rc else rc


def bh_adjust(values: Sequence[float]) -> np.ndarray:
    p = np.asarray(values, dtype=float)
    out = np.full(p.shape, np.nan, dtype=float)
    good = np.isfinite(p)
    if not good.any():
        return out
    vals = p[good]
    order = np.argsort(vals)
    ranked = vals[order]
    n = len(ranked)
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    inv = np.empty(n, dtype=int)
    inv[order] = np.arange(n)
    out[good] = adj[inv]
    return out


def require(path: Path) -> Path:
    if not path.exists() or path.stat().st_size == 0:
        raise RuntimeError(f"Required file missing or empty: {path}")
    return path


def read_fasta(path: Path) -> Dict[str, str]:
    seqs: Dict[str, List[str]] = {}
    current: str | None = None
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                if current in seqs:
                    raise RuntimeError(f"Duplicate FASTA record: {current} in {path}")
                seqs[current] = []
            else:
                if current is None:
                    raise RuntimeError(f"Sequence before FASTA header in {path}")
                seqs[current].append(line.upper())
    return {key: "".join(parts) for key, parts in seqs.items()}


def write_fasta(records: Iterable[Tuple[str, str]], path: Path, width: int = 80) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wt", encoding="utf-8", newline="\n") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            for i in range(0, len(seq), width):
                handle.write(seq[i:i + width] + "\n")


def run_command(command: List[str], stdout_path: Path) -> None:
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    print("RUN:", " ".join(shlex.quote(x) for x in command), flush=True)
    with stdout_path.open("wt", encoding="utf-8", newline="\n") as out:
        subprocess.run(command, check=True, stdout=out)


def parse_miniprot_paf(path: Path, species: str) -> pd.DataFrame:
    rows: List[dict] = []
    with path.open("rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw.startswith("##PAF"):
                continue
            fields = raw.rstrip("\n").split("\t")
            f = fields[1:] if fields and fields[0] == "##PAF" else fields
            if len(f) < 12:
                continue
            try:
                qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend, nmatch, blen, mapq = f[:12]
                tags = {}
                for token in f[12:]:
                    if token.count(":") >= 2:
                        key, _kind, value = token.split(":", 2)
                        tags[key] = value
                qlen_i = int(qlen)
                qstart_i = int(qstart)
                qend_i = int(qend)
                nmatch_i = int(nmatch)
                blen_i = int(blen)
                rows.append({
                    "species": species,
                    "gene_id": qname,
                    "query_length_aa": qlen_i,
                    "query_start": qstart_i,
                    "query_end": qend_i,
                    "query_coverage_pct": 100.0 * (qend_i - qstart_i) / max(1, qlen_i),
                    "target_contig": tname,
                    "target_length": int(tlen),
                    "target_start_0based": int(tstart),
                    "target_end_0based_exclusive": int(tend),
                    "strand": strand,
                    "aligned_matches": nmatch_i,
                    "alignment_block": blen_i,
                    "alignment_identity_pct": 100.0 * nmatch_i / max(1, blen_i),
                    "mapq": int(mapq),
                    "alignment_score": float(tags.get("AS", "nan")),
                    "matching_score": float(tags.get("ms", "nan")),
                    "cigar": tags.get("cg", ""),
                    "cs": tags.get("cs", ""),
                })
            except (ValueError, IndexError):
                continue
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    df["alignment_score_sort"] = pd.to_numeric(df["alignment_score"], errors="coerce").fillna(-np.inf)
    df = df.sort_values(
        ["gene_id", "query_coverage_pct", "alignment_score_sort", "alignment_identity_pct", "alignment_block"],
        ascending=[True, False, False, False, False],
    )
    output_rows: List[dict] = []
    for gene_id, sub in df.groupby("gene_id", sort=False):
        sub = sub.reset_index(drop=True)
        top = sub.iloc[0].to_dict()
        top_score = float(top.get("alignment_score", np.nan))
        second_score = float(sub.iloc[1]["alignment_score"]) if len(sub) > 1 else np.nan
        if np.isfinite(top_score) and np.isfinite(second_score) and second_score > 0:
            score_ratio = top_score / second_score
        else:
            score_ratio = np.inf if len(sub) == 1 else np.nan
        qcov = float(top["query_coverage_pct"])
        ident = float(top["alignment_identity_pct"])
        if qcov >= 80 and ident >= 50:
            confidence = "HIGH" if (len(sub) == 1 or (np.isfinite(score_ratio) and score_ratio >= 1.02)) else "MODERATE_PARALOG_AMBIGUITY"
        elif qcov >= 70 and ident >= 40:
            confidence = "MODERATE"
        else:
            confidence = "LOW"
        top["alignment_count"] = int(len(sub))
        top["second_alignment_score"] = second_score
        top["top_second_score_ratio"] = score_ratio
        top["mapping_confidence"] = confidence
        top.pop("alignment_score_sort", None)
        output_rows.append(top)
    return pd.DataFrame(output_rows)


def extract_oriented_promoter(row: pd.Series, genome: Dict[str, str], promoter_bp: int) -> Tuple[str, int, int, str]:
    contig = str(row["target_contig"])
    if contig not in genome:
        return "", 0, 0, "contig_not_found"
    seq = genome[contig]
    start = int(row["target_start_0based"])
    end = int(row["target_end_0based_exclusive"])
    strand = str(row["strand"])
    if strand == "+":
        left = max(0, start - promoter_bp)
        right = start
        promoter = seq[left:right]
    elif strand == "-":
        left = end
        right = min(len(seq), end + promoter_bp)
        promoter = revcomp(seq[left:right])
    else:
        return "", 0, 0, "invalid_strand"
    if not promoter:
        return "", left, right, "empty"
    n_fraction = promoter.count("N") / len(promoter)
    status = "ok" if len(promoter) >= 200 and n_fraction <= 0.2 else "low_quality"
    return promoter, left, right, status


def scan_motif(sequence: str, motif: str) -> dict:
    motif = motif.upper()
    rc = revcomp(motif)
    hits: List[Tuple[int, str]] = []
    for pattern, orientation in ((motif, "+"), (rc, "-")):
        for match in re.finditer(f"(?={re.escape(pattern)})", sequence):
            hits.append((match.start(), orientation))
    hits.sort(key=lambda x: (x[0], x[1]))
    positions = [start + 1 for start, _ in hits]
    distances = [len(sequence) - (start + len(motif)) for start, _ in hits]
    return {
        "motif_present": bool(hits),
        "motif_count": len(hits),
        "motif_positions_1based": ",".join(map(str, positions)),
        "motif_orientations": ",".join(ori for _, ori in hits),
        "motif_distances_to_gene_bp": ",".join(map(str, distances)),
        "nearest_gene_distance_bp": min(distances) if distances else np.nan,
        "nearest_hit_start_0based": hits[int(np.argmin(distances))][0] if hits else np.nan,
        "nearest_hit_orientation": hits[int(np.argmin(distances))][1] if hits else "",
    }


def best_window_identity(query: str, target: str) -> float:
    query = query.upper()
    target = target.upper()
    k = len(query)
    if not query or len(target) < k:
        return float("nan")
    q_rc = revcomp(query)
    best = 0
    for i in range(0, len(target) - k + 1):
        window = target[i:i + k]
        score1 = sum(a == b for a, b in zip(query, window))
        score2 = sum(a == b for a, b in zip(q_rc, window))
        best = max(best, score1, score2)
    return best / k


def flanking_window(sequence: str, hit_start: int, motif_len: int, flank_each_side: int = 12) -> str:
    center_left = max(0, hit_start - flank_each_side)
    center_right = min(len(sequence), hit_start + motif_len + flank_each_side)
    fragment = sequence[center_left:center_right]
    desired = motif_len + 2 * flank_each_side
    if len(fragment) < desired:
        fragment = ("N" * max(0, flank_each_side - hit_start)) + fragment
        fragment = fragment + ("N" * max(0, desired - len(fragment)))
    return fragment[:desired]


def draw_matched_set(pools: Dict[str, List[str]], valid_genes: set[str], rng: np.random.Generator, max_tries: int = 200) -> List[str]:
    keys = list(pools)
    for _ in range(max_tries):
        selected: List[str] = []
        used: set[str] = set()
        order = keys.copy()
        rng.shuffle(order)
        ok = True
        for focus_gene in order:
            candidates = [x for x in pools[focus_gene] if x in valid_genes and x not in used]
            if not candidates:
                ok = False
                break
            pick = str(rng.choice(candidates))
            selected.append(pick)
            used.add(pick)
        if ok and len(selected) == len(keys):
            return selected
    union = sorted({x for key in keys for x in pools[key] if x in valid_genes})
    if len(union) < len(keys):
        return []
    return list(rng.choice(union, size=len(keys), replace=False))


def empirical_species_p(
    focus_genes: List[str],
    motif_map: Dict[str, bool],
    pools: Dict[str, List[str]],
    iterations: int,
    rng: np.random.Generator,
) -> Tuple[float, int, int]:
    valid_focus = [g for g in focus_genes if g in motif_map]
    local_pools = {g: pools[g] for g in valid_focus if g in pools}
    valid_focus = [g for g in valid_focus if g in local_pools]
    if len(valid_focus) < 5:
        return float("nan"), sum(bool(motif_map.get(g, False)) for g in valid_focus), len(valid_focus)
    observed = sum(bool(motif_map[g]) for g in valid_focus)
    valid_background = set(motif_map) - set(focus_genes)
    extreme = 0
    successful = 0
    for _ in range(iterations):
        chosen = draw_matched_set(local_pools, valid_background, rng)
        if len(chosen) != len(valid_focus):
            continue
        successful += 1
        score = sum(bool(motif_map[g]) for g in chosen)
        extreme += int(score >= observed)
    if successful == 0:
        return float("nan"), observed, len(valid_focus)
    return (extreme + 1) / (successful + 1), observed, len(valid_focus)



def main() -> None:
    parser = argparse.ArgumentParser(description="Test a fixed promoter motif in ortholog promoters.")
    parser.add_argument("--motif-table", type=Path, required=True)
    parser.add_argument("--focus-promoters", type=Path, required=True)
    parser.add_argument("--matched-pools", type=Path, required=True)
    parser.add_argument("--focus-promoter-fasta", type=Path, required=True)
    parser.add_argument("--background-promoter-fasta", type=Path, required=True)
    parser.add_argument("--source-proteins", type=Path, required=True)
    parser.add_argument("--gene-annotation", type=Path, required=True)
    parser.add_argument("--genome", action="append", default=[], help="LABEL=GENOME_FASTA")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=12)
    parser.add_argument("--promoter-bp", type=int, default=1000)
    parser.add_argument("--matched-iterations", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=20260716)
    parser.add_argument("--motif", default="AUTO")
    parser.add_argument("--miniprot", default="miniprot")
    args = parser.parse_args()

    out = args.output_dir
    tables = out / "tables"
    sequences = out / "sequences"
    mappings = out / "mappings"
    for directory in (out, args.work_dir, tables, sequences, mappings):
        directory.mkdir(parents=True, exist_ok=True)

    motif_table = pd.read_csv(require(args.motif_table), sep="\t")
    if args.motif.upper() == "AUTO":
        robust = motif_table[motif_table["robust"].astype(str).str.lower().isin({"true", "1", "yes"})]
        if len(robust) != 1:
            raise RuntimeError(f"Expected one robust motif, found {len(robust)}")
        motif = canonical_kmer(str(robust.iloc[0]["motif"]))
    else:
        motif = canonical_kmer(args.motif)

    focus_info = pd.read_csv(require(args.focus_promoters), sep="\t", dtype=str).fillna("")
    pools_table = pd.read_csv(require(args.matched_pools), sep="\t", dtype=str).fillna("")
    focus_fasta = read_fasta(require(args.focus_promoter_fasta))
    background_fasta = read_fasta(require(args.background_promoter_fasta))
    proteins = read_fasta(require(args.source_proteins))
    annotation = pd.read_csv(require(args.gene_annotation), sep="\t", dtype=str).fillna("")

    focus_genes = focus_info["gene_id"].astype(str).tolist()
    pools = {
        str(gene): subset.sort_values("rank")["matched_gene"].astype(str).tolist()
        for gene, subset in pools_table.groupby("focus_gene")
    }
    background_genes = sorted(set(pools_table["matched_gene"].astype(str)))
    query_genes = focus_genes + [gene for gene in background_genes if gene not in focus_genes]
    missing = [gene for gene in query_genes if gene not in proteins]
    if missing:
        raise RuntimeError(f"Missing source proteins: {missing[:5]}")
    query_faa = sequences / "query_proteins.faa"
    write_fasta(((gene, proteins[gene]) for gene in query_genes), query_faa)

    genomes = {"FTH_FL4": None}
    for item in args.genome:
        if "=" not in item:
            raise RuntimeError(f"Invalid --genome value: {item}")
        label, path = item.split("=", 1)
        genomes[label] = require(Path(path))

    promoter_rows = []
    promoter_sequences: Dict[Tuple[str, str], str] = {}
    focus_map = focus_info.set_index("gene_id").to_dict("index")
    background_info_path = args.background_promoter_fasta.parent.parent / "tables" / "matched_background_promoters.tsv"
    background_map = pd.read_csv(background_info_path, sep="\t", dtype=str).fillna("").set_index("gene_id").to_dict("index") if background_info_path.exists() else {}
    for gene in query_genes:
        sequence = focus_fasta.get(gene, background_fasta.get(gene, ""))
        if not sequence:
            continue
        info = focus_map.get(gene, background_map.get(gene, {}))
        call = scan_motif(sequence, motif)
        promoter_sequences[("FTH_FL4", gene)] = sequence
        promoter_rows.append({
            "species": "FTH_FL4", "gene_id": gene, "is_focus": gene in focus_genes,
            "mapping_confidence": "reference", "query_coverage_pct": 100.0,
            "alignment_identity_pct": 100.0, "target_contig": info.get("fasta_contig", info.get("contig", "")),
            "strand": info.get("strand", ""), "promoter_length": len(sequence),
            "promoter_status": "ok", **call,
        })

    mapping_tables = []
    for species, genome_path in genomes.items():
        if species == "FTH_FL4":
            continue
        mapping_file = mappings / f"{species}.miniprot.gff3"
        run_command([args.miniprot, "-t", str(args.threads), "-I", "--gff", "--aln", str(genome_path), str(query_faa)], mapping_file)
        mapping = parse_miniprot_paf(mapping_file, species)
        if mapping.empty:
            raise RuntimeError(f"No mappings for {species}")
        mapping_tables.append(mapping)
        genome = read_fasta(genome_path)
        for _, row in mapping.iterrows():
            gene = str(row["gene_id"])
            if gene not in query_genes:
                continue
            promoter, left, right, status = extract_oriented_promoter(row, genome, args.promoter_bp)
            call = scan_motif(promoter, motif) if promoter else {
                "motif_present": False, "motif_count": 0, "motif_positions_1based": "",
                "motif_orientations": "", "motif_distances_to_gene_bp": "",
                "nearest_gene_distance_bp": np.nan, "nearest_hit_start_0based": np.nan,
                "nearest_hit_orientation": "",
            }
            if promoter:
                promoter_sequences[(species, gene)] = promoter
            promoter_rows.append({
                "species": species, "gene_id": gene, "is_focus": gene in focus_genes,
                "mapping_confidence": row["mapping_confidence"],
                "query_coverage_pct": row["query_coverage_pct"],
                "alignment_identity_pct": row["alignment_identity_pct"],
                "target_contig": row["target_contig"], "strand": row["strand"],
                "promoter_start": left, "promoter_end": right, "promoter_length": len(promoter),
                "promoter_status": status, **call,
            })
    calls = pd.DataFrame(promoter_rows)
    valid_confidence = {"reference", "HIGH", "MODERATE", "MODERATE_PARALOG_AMBIGUITY"}
    calls["analysis_valid"] = calls["promoter_status"].eq("ok") & calls["mapping_confidence"].isin(valid_confidence)
    calls.to_csv(tables / "ortholog_promoter_motif_calls.tsv", sep="\t", index=False, lineterminator="\n", na_rep="NA")
    if mapping_tables:
        pd.concat(mapping_tables, ignore_index=True).to_csv(tables / "ortholog_mapping.tsv", sep="\t", index=False, lineterminator="\n")
    write_fasta(
        ((f"{species}|{gene}", sequence) for (species, gene), sequence in promoter_sequences.items() if gene in focus_genes),
        sequences / "focus_ortholog_promoters.fasta",
    )

    rng = np.random.default_rng(args.seed)
    species_rows = []
    motif_maps: Dict[str, Dict[str, bool]] = {}
    for species in genomes:
        subset = calls[(calls["species"].eq(species)) & calls["analysis_valid"]].copy()
        motif_map = dict(zip(subset["gene_id"].astype(str), subset["motif_present"].astype(bool)))
        motif_maps[species] = motif_map
        focus_valid = [gene for gene in focus_genes if gene in motif_map]
        background_valid = set(gene for gene in background_genes if gene in motif_map)
        observed = sum(motif_map[gene] for gene in focus_valid)
        p, successful = empirical_species_p(
            motif_map, focus_valid, pools, background_valid, args.matched_iterations, rng
        )
        species_rows.append({
            "species": species,
            "focus_promoters_valid": len(focus_valid),
            "focus_motif_positive": observed,
            "focus_fraction": observed / len(focus_valid) if focus_valid else np.nan,
            "matched_empirical_p": p,
            "successful_iterations": successful,
        })
    species_table = pd.DataFrame(species_rows)

    fl4 = calls[(calls["species"].eq("FTH_FL4")) & calls["is_focus"]].set_index("gene_id")
    conservation_rows = []
    for gene in focus_genes:
        row = {"gene_id": gene, "fl4_motif_present": bool(fl4.loc[gene, "motif_present"])}
        fl4_distance = float(fl4.loc[gene, "nearest_gene_distance_bp"]) if pd.notna(fl4.loc[gene, "nearest_gene_distance_bp"]) else np.nan
        row["fl4_distance_bp"] = fl4_distance
        for species in genomes:
            subset = calls[(calls["species"].eq(species)) & calls["gene_id"].eq(gene)]
            valid = bool(len(subset) and subset.iloc[0]["analysis_valid"])
            present = bool(valid and subset.iloc[0]["motif_present"])
            distance = float(subset.iloc[0]["nearest_gene_distance_bp"]) if valid and pd.notna(subset.iloc[0]["nearest_gene_distance_bp"]) else np.nan
            row[f"{species}_valid"] = valid
            row[f"{species}_motif_present"] = present
            row[f"{species}_distance_bp"] = distance
            row[f"{species}_position_delta_bp"] = abs(distance - fl4_distance) if present and np.isfinite(fl4_distance) and np.isfinite(distance) else np.nan
        conservation_rows.append(row)
    conservation = pd.DataFrame(conservation_rows).merge(
        annotation[[column for column in ("gene_id", "product", "description_text") if column in annotation]],
        on="gene_id", how="left",
    )
    conservation.to_csv(tables / "focus_motif_conservation.tsv", sep="\t", index=False, lineterminator="\n", na_rep="NA")

    fl4_positive = set(conservation.loc[conservation["fl4_motif_present"], "gene_id"])
    retention_rows = []
    for species in genomes:
        if species == "FTH_FL4":
            continue
        valid = conservation[conservation["gene_id"].isin(fl4_positive) & conservation[f"{species}_valid"]]
        retention_rows.append({
            "species": species,
            "fl4_positive_orthologs_valid": len(valid),
            "motif_retained": int(valid[f"{species}_motif_present"].sum()),
            "retention_fraction": float(valid[f"{species}_motif_present"].mean()) if len(valid) else np.nan,
        })
    retention = pd.DataFrame(retention_rows)
    species_table = species_table.merge(retention, on="species", how="left")
    species_table.to_csv(tables / "species_motif_enrichment.tsv", sep="\t", index=False, lineterminator="\n", na_rep="NA")

    summary = {
        "motif": motif,
        "reverse_complement": revcomp(motif),
        "species": species_table.to_dict("records"),
    }
    (tables / "motif_conservation.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
