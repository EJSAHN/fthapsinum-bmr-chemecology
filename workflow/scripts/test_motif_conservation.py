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
from motif_inference import choose_follow_up


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
    patterns = [(motif, "+")] if motif == rc else [(motif, "+"), (rc, "-")]
    for pattern, orientation in patterns:
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



def parse_label_paths(values: Sequence[str]) -> dict[str, Path]:
    paths = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"Expected LABEL=PATH, received {value}")
        label, text = value.split("=", 1)
        if not label or label in paths:
            raise ValueError(f"Empty or duplicate reference label: {label}")
        paths[label] = require(Path(text))
    return paths


def write_no_selection(tables: Path, decision: dict) -> None:
    # Empty tables replace stale canonical summaries, including in reused output folders.
    schemas = {
        "ortholog_promoter_motif_calls.tsv": ["species", "gene_id", "is_focus", "analysis_valid", "motif_present"],
        "ortholog_mapping_summary.tsv": ["species", "gene_id", "target_contig"],
        "species_motif_enrichment.tsv": ["species", "focus_motif_positive", "focus_promoters_valid", "matched_empirical_p"],
        "focus_gene_motif_conservation.tsv": ["gene_id", "fl4_motif_present"],
        "focus_motif_conservation.tsv": ["gene_id", "fl4_motif_present"],
        "focus_motif_position_conservation.tsv": ["gene_id", "species", "valid", "motif_present"],
        "pooled_cross_species_motif_tests.tsv": ["test", "species", "matched_empirical_p"],
        "sequence_record_layout.tsv": ["species", "projected_genes", "unique_sequence_records"],
        "fixed_motif_definition.tsv": ["motif", "selection_mode", "full_family_q"],
    }
    for name, columns in schemas.items():
        pd.DataFrame(columns=columns).to_csv(tables / name, sep="\t", index=False)
    (tables / "motif_conservation.json").write_text(json.dumps(decision, indent=2) + "\n", encoding="utf-8")


def sequence_record_layout(calls: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for species, sub in calls[calls["is_focus"] & calls["analysis_valid"]].groupby("species", sort=True):
        sub = sub.drop_duplicates("gene_id")
        gaps = []
        records = sub.to_dict("records")
        for i, left in enumerate(records):
            for right in records[i + 1:]:
                if left["target_contig"] != right["target_contig"]:
                    continue
                a, b = sorted((left, right), key=lambda r: int(r["target_start_0based"]))
                gaps.append(max(0, int(b["target_start_0based"]) - int(a["target_end_0based_exclusive"])))
        rows.append({
            "species": species, "projected_genes": len(sub),
            "unique_sequence_records": sub["target_contig"].nunique(),
            "same_sequence_pairs": len(gaps),
            "minimum_within_sequence_intergenic_gap_bp": min(gaps) if gaps else np.nan,
            "pairs_gap_le_20kb": sum(g <= 20000 for g in gaps),
            "pairs_gap_le_50kb": sum(g <= 50000 for g in gaps),
            "pairs_gap_le_100kb": sum(g <= 100000 for g in gaps),
            "between_record_distances": "unknown",
        })
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description="Compare an explicitly selected promoter sequence across projected loci.")
    ap.add_argument("--motif-table", type=Path, required=True)
    ap.add_argument("--focus-promoters", type=Path, required=True)
    ap.add_argument("--matched-pools", type=Path, required=True)
    ap.add_argument("--focus-promoter-fasta", type=Path, required=True)
    ap.add_argument("--background-promoter-fasta", type=Path, required=True)
    ap.add_argument("--source-proteins", type=Path, required=True)
    ap.add_argument("--gene-annotation", type=Path, required=True)
    ap.add_argument("--genome", action="append", default=[], help="LABEL=GENOME_FASTA")
    ap.add_argument("--mapping", action="append", default=[], help="LABEL=EXISTING_MINIPROT_GFF; optional upstream projection input")
    ap.add_argument("--reference-label", default="FTH_FL4")
    ap.add_argument("--within-species-label", default="FTH_NRRL22049")
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--work-dir", type=Path, required=True)
    ap.add_argument("--threads", type=int, default=12)
    ap.add_argument("--promoter-bp", type=int, default=1000)
    ap.add_argument("--matched-iterations", type=int, default=20000)
    ap.add_argument("--position-tolerance-bp", type=int, default=150)
    ap.add_argument("--seed", type=int, default=20260716)
    ap.add_argument("--selection-mode", choices=["auto", "exploratory"], default="auto")
    ap.add_argument("--selection-fdr", type=float, default=0.10)
    ap.add_argument("--motif", default="AUTO")
    ap.add_argument("--miniprot", default="miniprot")
    args = ap.parse_args()
    if args.matched_iterations < 1 or args.promoter_bp < 1:
        raise ValueError("Iteration count and promoter window must be positive")
    out, work = args.output_dir, args.work_dir
    tables, sequences, mappings = out / "tables", out / "sequences", out / "mappings"
    for directory in (out, work, tables, sequences, mappings):
        directory.mkdir(parents=True, exist_ok=True)
    motif_table = pd.read_csv(require(args.motif_table), sep="\t")
    selection = choose_follow_up(motif_table, args.selection_mode, args.motif, args.selection_fdr)
    (tables / "candidate_selection.json").write_text(json.dumps(selection, indent=2) + "\n", encoding="utf-8")
    if selection["motif"] is None:
        write_no_selection(tables, selection)
        print("MOTIF_FOLLOW_UP=" + selection["status"], flush=True)
        return
    motif, motif_rc = selection["motif"], selection["reverse_complement"]
    focus_promoters = pd.read_csv(require(args.focus_promoters), sep="\t")
    matched_pools_df = pd.read_csv(require(args.matched_pools), sep="\t")
    focus_fasta = read_fasta(require(args.focus_promoter_fasta))
    background_fasta = read_fasta(require(args.background_promoter_fasta))
    source_proteins = read_fasta(require(args.source_proteins))
    gene_annotation = pd.read_csv(require(args.gene_annotation), sep="\t")
    projection_inputs = parse_label_paths(args.mapping)
    species_genomes = {args.reference_label: None, **parse_label_paths(args.genome)}
    if args.reference_label in parse_label_paths(args.genome):
        raise ValueError("The focal reference uses the supplied promoters; do not repeat it in --genome")
    if set(projection_inputs) - set(species_genomes):
        raise ValueError("Every supplied mapping must have a matching --genome")
    matched_pools_df["rank"] = pd.to_numeric(matched_pools_df["rank"], errors="raise")
    focus_genes = focus_promoters["gene_id"].astype(str).tolist()
    if len(set(focus_genes)) != len(focus_genes) or not focus_genes:
        raise ValueError("Focus identifiers must be nonempty and unique")
    pools: Dict[str, List[str]] = {
        str(g): sub.sort_values("rank")["matched_gene"].astype(str).tolist()
        for g, sub in matched_pools_df.groupby("focus_gene")
    }
    background_genes = sorted(set(matched_pools_df["matched_gene"].astype(str)))
    query_genes = focus_genes + [g for g in background_genes if g not in focus_genes]
    missing_proteins = [g for g in query_genes if g not in source_proteins]
    if missing_proteins:
        raise RuntimeError(f"Source proteins missing for {len(missing_proteins)} queried genes; first={missing_proteins[:5]}")

    query_faa = sequences / "focus_and_matched_source_proteins.faa"
    write_fasta(((g, source_proteins[g]) for g in query_genes), query_faa)

    promoter_rows: List[dict] = []
    mapping_rows: List[pd.DataFrame] = []
    promoter_sequences: Dict[Tuple[str, str], str] = {}

    if set(focus_genes) != set(focus_fasta) or set(background_genes) != set(background_fasta):
        raise ValueError("FASTA records do not match focus and matched-background tables")
    if set(focus_genes) & set(background_genes):
        raise ValueError("Foreground and background overlap")
    # Use the focal-reference promoter sequences from the upstream analysis.
    fl4_info = focus_promoters.set_index("gene_id").to_dict("index")
    background_info = pd.read_csv(require(args.background_promoter_fasta.parent.parent / "tables" / "matched_background_promoters.tsv"), sep="\t").set_index("gene_id").to_dict("index")
    for gene_id in query_genes:
        seq = focus_fasta.get(gene_id, background_fasta.get(gene_id, ""))
        info = fl4_info.get(gene_id, background_info.get(gene_id, {}))
        if not seq:
            continue
        scan = scan_motif(seq, motif)
        promoter_sequences[(args.reference_label, gene_id)] = seq
        promoter_rows.append({
            "species": args.reference_label,
            "gene_id": gene_id,
            "is_focus": gene_id in focus_genes,
            "mapping_confidence": "reference",
            "query_coverage_pct": 100.0,
            "alignment_identity_pct": 100.0,
            "alignment_count": 1,
            "top_second_score_ratio": np.inf,
            "target_contig": info.get("fasta_contig", info.get("contig", "")),
            "target_start_0based": int(info.get("gene_start", 1)) - 1 if info else np.nan,
            "target_end_0based_exclusive": int(info.get("gene_end", 0)) if info else np.nan,
            "strand": info.get("strand", ""),
            "promoter_length": len(seq),
            "promoter_gc": (seq.count("G") + seq.count("C")) / max(1, len(seq)),
            "promoter_n_fraction": seq.count("N") / max(1, len(seq)),
            "promoter_status": "ok",
            **scan,
        })

    for species, genome_path in species_genomes.items():
        if species == args.reference_label:
            continue
        map_path = mappings / f"{species}.miniprot.gff3"
        if species in projection_inputs:
            # Explicit supplied upstream alignment; never download or overwrite it.
            map_path = projection_inputs[species]
        else:
            run_command([args.miniprot, "-t", str(args.threads), "-I", "--gff", "--aln",
                         str(genome_path), str(query_faa)], map_path)
        mapping = parse_miniprot_paf(map_path, species)
        if mapping.empty:
            raise RuntimeError(f"No miniprot mappings parsed for {species}")
        mapping_rows.append(mapping)
        genome = read_fasta(genome_path)
        for _, row in mapping.iterrows():
            gene_id = str(row["gene_id"])
            if gene_id not in query_genes:
                continue
            if int(row["query_length_aa"]) != len(source_proteins[gene_id]):
                raise ValueError(f"Projection query length mismatch: {species} {gene_id}")
            contig = str(row["target_contig"])
            if contig not in genome or int(row["target_length"]) != len(genome[contig]):
                raise ValueError(f"Projection reference length mismatch: {species} {contig}")
            promoter, left, right, status = extract_oriented_promoter(row, genome, args.promoter_bp)
            scan = scan_motif(promoter, motif) if promoter else {
                "motif_present": False, "motif_count": 0, "motif_positions_1based": "",
                "motif_orientations": "", "motif_distances_to_gene_bp": "",
                "nearest_gene_distance_bp": np.nan, "nearest_hit_start_0based": np.nan,
                "nearest_hit_orientation": "",
            }
            if promoter:
                promoter_sequences[(species, gene_id)] = promoter
            promoter_rows.append({
                "species": species,
                "gene_id": gene_id,
                "is_focus": gene_id in focus_genes,
                "mapping_confidence": row["mapping_confidence"],
                "query_coverage_pct": row["query_coverage_pct"],
                "alignment_identity_pct": row["alignment_identity_pct"],
                "alignment_count": row["alignment_count"],
                "top_second_score_ratio": row["top_second_score_ratio"],
                "target_contig": row["target_contig"],
                "target_start_0based": row["target_start_0based"],
                "target_end_0based_exclusive": row["target_end_0based_exclusive"],
                "strand": row["strand"],
                "promoter_genomic_left_0based": left,
                "promoter_genomic_right_0based_exclusive": right,
                "promoter_length": len(promoter),
                "promoter_gc": (promoter.count("G") + promoter.count("C")) / max(1, len(promoter)) if promoter else np.nan,
                "promoter_n_fraction": promoter.count("N") / max(1, len(promoter)) if promoter else np.nan,
                "promoter_status": status,
                **scan,
            })

    promoter_df = pd.DataFrame(promoter_rows)
    valid_conf = {"reference", "HIGH", "MODERATE_PARALOG_AMBIGUITY", "MODERATE"}
    promoter_df["analysis_valid"] = promoter_df["promoter_status"].eq("ok") & promoter_df["mapping_confidence"].isin(valid_conf)
    promoter_df.to_csv(tables / "ortholog_promoter_motif_calls.tsv", sep="\t", index=False)
    if mapping_rows:
        pd.concat(mapping_rows, ignore_index=True).to_csv(tables / "ortholog_mapping_summary.tsv", sep="\t", index=False)

    # Write focus ortholog promoter FASTA.
    focus_records = []
    for (species, gene_id), seq in promoter_sequences.items():
        if gene_id in focus_genes:
            focus_records.append((f"{species}|{gene_id}", seq))
    write_fasta(focus_records, sequences / "focus_ortholog_promoters.fasta")

    fl4_focus = promoter_df[(promoter_df["species"] == args.reference_label) & promoter_df["is_focus"]]
    observed_focus = int(fl4_focus["motif_present"].sum())
    if observed_focus != selection["foreground_present"] or len(fl4_focus) != selection["foreground_total"]:
        raise ValueError("Follow-up focal counts disagree with the full-family discovery table")

    rng = np.random.default_rng(args.seed)
    species_rows: List[dict] = []
    motif_maps_by_species: Dict[str, Dict[str, bool]] = {}
    for species in species_genomes:
        sub = promoter_df[(promoter_df["species"] == species) & promoter_df["analysis_valid"]].copy()
        motif_map = dict(zip(sub["gene_id"].astype(str), sub["motif_present"].astype(bool)))
        motif_maps_by_species[species] = motif_map
        focus_sub = sub[sub["gene_id"].isin(focus_genes)]
        bg_sub = sub[sub["gene_id"].isin(background_genes)]
        focus_pos = int(focus_sub["motif_present"].sum())
        focus_total = int(len(focus_sub))
        bg_pos = int(bg_sub["motif_present"].sum())
        bg_total = int(len(bg_sub))
        if focus_total > 0 and bg_total > 0:
            odds, fisher_p = fisher_exact([[focus_pos, focus_total - focus_pos], [bg_pos, bg_total - bg_pos]], alternative="greater")
        else:
            odds, fisher_p = np.nan, np.nan
        empirical_p, observed, empirical_n = empirical_species_p(focus_genes, motif_map, pools, args.matched_iterations, rng)
        species_rows.append({
            "species": species,
            "focus_motif_positive": focus_pos,
            "focus_promoters_valid": focus_total,
            "background_motif_positive": bg_pos,
            "background_promoters_valid": bg_total,
            "focus_fraction": focus_pos / focus_total if focus_total else np.nan,
            "background_fraction": bg_pos / bg_total if bg_total else np.nan,
            "odds_ratio": odds,
            "fisher_p": fisher_p,
            "matched_empirical_p": empirical_p,
            "empirical_focus_positive": observed,
            "empirical_focus_n": empirical_n,
        })
    species_df = pd.DataFrame(species_rows)
    species_df["fisher_fdr"] = bh_adjust(species_df["fisher_p"])
    species_df["matched_empirical_fdr"] = bh_adjust(species_df["matched_empirical_p"])

    # Gene-wise conservation relative to the fixed FL-4 motif call.
    fl4_lookup = fl4_focus.set_index("gene_id")
    conservation_rows: List[dict] = []
    position_rows: List[dict] = []
    for gene_id in focus_genes:
        fl4_row = fl4_lookup.loc[gene_id]
        fl4_present = bool(fl4_row["motif_present"])
        fl4_distance = float(fl4_row["nearest_gene_distance_bp"]) if pd.notna(fl4_row["nearest_gene_distance_bp"]) else np.nan
        fl4_seq = promoter_sequences.get((args.reference_label, gene_id), "")
        ref_fragment = ""
        if fl4_present and pd.notna(fl4_row["nearest_hit_start_0based"]):
            ref_fragment = flanking_window(fl4_seq, int(fl4_row["nearest_hit_start_0based"]), len(motif), 12)
        row_out = {"gene_id": gene_id, "fl4_motif_present": fl4_present, "fl4_nearest_gene_distance_bp": fl4_distance}
        valid_species = 0
        positive_species = 0
        external_positive = 0
        external_valid = 0
        for species in species_genomes:
            sub = promoter_df[(promoter_df["species"] == species) & (promoter_df["gene_id"] == gene_id)]
            if sub.empty:
                valid = False
                present = False
                distance = np.nan
                confidence = "MISSING"
                fragment_identity = np.nan
            else:
                r = sub.iloc[0]
                valid = bool(r["analysis_valid"])
                present = bool(r["motif_present"]) if valid else False
                distance = float(r["nearest_gene_distance_bp"]) if valid and pd.notna(r["nearest_gene_distance_bp"]) else np.nan
                confidence = str(r["mapping_confidence"])
                target_seq = promoter_sequences.get((species, gene_id), "")
                fragment_identity = best_window_identity(ref_fragment, target_seq) if ref_fragment and target_seq else np.nan
            if valid:
                valid_species += 1
                positive_species += int(present)
                if species not in {args.reference_label, args.within_species_label}:
                    external_valid += 1
                    external_positive += int(present)
            pos_delta = abs(distance - fl4_distance) if fl4_present and present and np.isfinite(distance) and np.isfinite(fl4_distance) else np.nan
            position_conserved = bool(np.isfinite(pos_delta) and pos_delta <= args.position_tolerance_bp)
            row_out[f"{species}_valid"] = valid
            row_out[f"{species}_motif_present"] = present
            row_out[f"{species}_nearest_gene_distance_bp"] = distance
            row_out[f"{species}_mapping_confidence"] = confidence
            row_out[f"{species}_flank31_best_identity"] = fragment_identity
            position_rows.append({
                "gene_id": gene_id,
                "species": species,
                "valid": valid,
                "motif_present": present,
                "nearest_gene_distance_bp": distance,
                "fl4_distance_bp": fl4_distance,
                "position_delta_bp": pos_delta,
                "position_conserved_within_tolerance": position_conserved,
                "flank31_best_identity": fragment_identity,
            })
        row_out["valid_species_count"] = valid_species
        row_out["motif_positive_species_count"] = positive_species
        row_out["external_valid_species_count"] = external_valid
        row_out["external_motif_positive_species_count"] = external_positive
        conservation_rows.append(row_out)
    conservation_df = pd.DataFrame(conservation_rows)
    conservation_df = conservation_df.merge(
        gene_annotation[[c for c in ["gene_id", "source_accession", "source_product", "annotation_text"] if c in gene_annotation.columns]],
        on="gene_id", how="left",
    )
    conservation_df.to_csv(tables / "focus_gene_motif_conservation.tsv", sep="\t", index=False)
    pd.DataFrame(position_rows).to_csv(tables / "focus_motif_position_conservation.tsv", sep="\t", index=False)

    fl4_positive_genes = set(conservation_df.loc[conservation_df["fl4_motif_present"], "gene_id"])
    retention_rows = []
    for species in species_genomes:
        if species == args.reference_label:
            continue
        valid_col = f"{species}_valid"
        motif_col = f"{species}_motif_present"
        sub = conservation_df[conservation_df["gene_id"].isin(fl4_positive_genes) & conservation_df[valid_col].astype(bool)]
        retention_rows.append({
            "species": species,
            "fl4_positive_orthologs_valid": int(len(sub)),
            "motif_retained": int(sub[motif_col].astype(bool).sum()),
            "retention_fraction": float(sub[motif_col].astype(bool).mean()) if len(sub) else np.nan,
        })
    retention_df = pd.DataFrame(retention_rows)
    species_df = species_df.merge(retention_df, on="species", how="left")

    # Fixed-candidate pooled summaries preserve reference correspondence, not independent species replication.
    def pooled_empirical(species_list: List[str]) -> Tuple[float, float, int, int]:
        focus_valid_pairs = []
        for g in focus_genes:
            for s in species_list:
                if g in motif_maps_by_species.get(s, {}):
                    focus_valid_pairs.append((g, s))
        if len(focus_valid_pairs) < 15:
            return np.nan, np.nan, 0, len(focus_valid_pairs)
        observed = sum(int(motif_maps_by_species[s][g]) for g, s in focus_valid_pairs) / len(focus_valid_pairs)
        extreme = 0
        successful = 0
        for _ in range(args.matched_iterations):
            chosen_by_focus: Dict[str, str] = {}
            used: set[str] = set()
            ok = True
            for g in focus_genes:
                candidates = [x for x in pools[g] if x not in used and any(x in motif_maps_by_species.get(s, {}) for s in species_list)]
                if not candidates:
                    ok = False
                    break
                pick = str(rng.choice(candidates))
                chosen_by_focus[g] = pick
                used.add(pick)
            if not ok:
                continue
            vals = []
            for g, s in focus_valid_pairs:
                bg = chosen_by_focus[g]
                if bg in motif_maps_by_species.get(s, {}):
                    vals.append(int(motif_maps_by_species[s][bg]))
            if not vals:
                continue
            successful += 1
            random_fraction = sum(vals) / len(vals)
            extreme += int(random_fraction >= observed - 1e-12)
        p = (extreme + 1) / (successful + 1) if successful else np.nan
        return observed, p, successful, len(focus_valid_pairs)

    ext_species = [label for label in species_genomes if label not in {args.reference_label, args.within_species_label}]
    external_fraction, external_p, external_successful, external_pairs = pooled_empirical(ext_species)
    nonfl4_fraction, nonfl4_p, nonfl4_successful, nonfl4_pairs = pooled_empirical(([args.within_species_label] if args.within_species_label in species_genomes else []) + ext_species)
    pooled_df = pd.DataFrame([
        {"test": "external_species", "species": ",".join(ext_species), "focus_motif_fraction": external_fraction, "matched_empirical_p": external_p, "successful_iterations": external_successful, "valid_focus_species_pairs": external_pairs},
        {"test": "all_non_reference", "species": ",".join(([args.within_species_label] if args.within_species_label in species_genomes else []) + ext_species), "focus_motif_fraction": nonfl4_fraction, "matched_empirical_p": nonfl4_p, "successful_iterations": nonfl4_successful, "valid_focus_species_pairs": nonfl4_pairs},
    ])
    pooled_df.to_csv(tables / "pooled_cross_species_motif_tests.tsv", sep="\t", index=False)
    species_df["inference_scope"] = "fixed_candidate_exploratory_not_scan_adjusted"
    species_df["discovery_full_family_q"] = selection["full_family_q"]
    species_df.to_csv(tables / "species_motif_enrichment.tsv", sep="\t", index=False)

    conservation_df.to_csv(tables / "focus_motif_conservation.tsv", sep="\t", index=False, na_rep="NA")
    sequence_record_layout(promoter_df).to_csv(tables / "sequence_record_layout.tsv", sep="\t", index=False, na_rep="NA")
    summary = {
        **selection,
        "execution_status": "COMPLETED",
        "matched_iterations": args.matched_iterations,
        "fisher_alternative": "greater",
        "matched_tail": "upper",
        "external_pooled_matched_p": external_p if np.isfinite(external_p) else None,
        "external_pooled_fraction": external_fraction if np.isfinite(external_fraction) else None,
        "genomic_scope": "between_sequence_record_distances_unknown",
        "comparison_scope": "related_promoters_not_independent_infection_or_binding_validation",
        "projection_mode": "explicit_upstream_mappings" if projection_inputs else "new_miniprot_projections",
    }
    (tables / "motif_conservation.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    pd.DataFrame([selection]).to_csv(tables / "fixed_motif_definition.tsv", sep="\t", index=False, na_rep="NA")
    print("MOTIF_FOLLOW_UP=" + selection["status"], flush=True)


if __name__ == "__main__":
    main()
