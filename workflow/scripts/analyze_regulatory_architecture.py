from __future__ import annotations

import argparse
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, spearmanr
import statsmodels.api as sm


def bh_adjust(pvals: Sequence[float]) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
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


def revcomp(seq: str) -> str:
    table = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    return seq.translate(table)[::-1]


def canonical_kmer(kmer: str) -> str:
    rc = revcomp(kmer)
    return kmer if kmer <= rc else rc


def read_fasta(path: Path) -> Dict[str, str]:
    seqs: Dict[str, List[str]] = {}
    name = None
    with path.open('rt', encoding='utf-8', errors='replace') as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                name = line[1:].split()[0]
                if name in seqs:
                    raise RuntimeError(f'Duplicate FASTA record: {name}')
                seqs[name] = []
            else:
                if name is None:
                    raise RuntimeError('FASTA sequence before first header')
                seqs[name].append(line.upper())
    return {k: ''.join(v) for k, v in seqs.items()}


def find_contig_key(contig: str, seqs: Dict[str, str]) -> str | None:
    candidates = [contig]
    for prefix in ('FTH_FL4__', 'FT_FL4__', 'FTH__'):
        if contig.startswith(prefix):
            candidates.append(contig[len(prefix):])
    candidates.extend([x.split('__', 1)[-1] for x in candidates if '__' in x])
    for c in candidates:
        if c in seqs:
            return c
    # last-resort suffix match, only if unique
    hits = [k for k in seqs if contig.endswith(k) or k.endswith(contig)]
    return hits[0] if len(hits) == 1 else None


def extract_promoter(row: pd.Series, seqs: Dict[str, str], promoter_bp: int) -> Tuple[str, str, int, int, str]:
    contig = str(row['contig'])
    key = find_contig_key(contig, seqs)
    if key is None:
        return '', '', 0, 0, 'contig_not_found'
    seq = seqs[key]
    start = int(row['start'])
    end = int(row['end'])
    strand = str(row['strand'])
    if strand == '+':
        left = max(1, start - promoter_bp)
        right = start - 1
        prom = seq[left - 1:right]
    elif strand == '-':
        left = end + 1
        right = min(len(seq), end + promoter_bp)
        prom = revcomp(seq[left - 1:right])
    else:
        return '', key, 0, 0, 'invalid_strand'
    if not prom:
        return '', key, left, right, 'empty'
    n_frac = prom.count('N') / len(prom)
    status = 'ok' if len(prom) >= 200 and n_frac <= 0.2 else 'low_quality'
    return prom, key, left, right, status


def write_fasta(records: Iterable[Tuple[str, str]], path: Path, width: int = 80) -> None:
    with path.open('wt', encoding='utf-8', newline='\n') as fh:
        for name, seq in records:
            fh.write(f'>{name}\n')
            for i in range(0, len(seq), width):
                fh.write(seq[i:i+width] + '\n')


def design_matrix(meta: pd.DataFrame) -> np.ndarray:
    groups = pd.get_dummies(meta['group'].astype(str), prefix='group', drop_first=True, dtype=float)
    burden = np.log1p(pd.to_numeric(meta['ft_unique_rpm_input'], errors='coerce').fillna(0).to_numpy())
    X = np.column_stack([np.ones(len(meta)), groups.to_numpy(), burden])
    return X.astype(float)


def residualize_matrix(Y: np.ndarray, X: np.ndarray) -> np.ndarray:
    beta, *_ = np.linalg.lstsq(X, Y, rcond=None)
    return Y - X @ beta


def upper_mean_corr(mat: np.ndarray) -> Tuple[float, float, float, np.ndarray]:
    # mat: samples x genes
    if mat.shape[1] < 2:
        return float('nan'), float('nan'), float('nan'), np.eye(mat.shape[1])
    z = mat.copy().astype(float)
    z -= np.nanmean(z, axis=0, keepdims=True)
    sd = np.nanstd(z, axis=0, ddof=1, keepdims=True)
    sd[~np.isfinite(sd) | (sd == 0)] = np.nan
    z = z / sd
    corr = np.corrcoef(z, rowvar=False)
    iu = np.triu_indices(corr.shape[0], 1)
    vals = corr[iu]
    vals = vals[np.isfinite(vals)]
    if len(vals) == 0:
        return float('nan'), float('nan'), float('nan'), corr
    return float(np.mean(vals)), float(np.median(vals)), float(np.mean(vals > 0)), corr


def draw_matched_set(pools: Dict[str, List[str]], rng: np.random.Generator, max_tries: int = 200) -> List[str]:
    keys = list(pools)
    for _ in range(max_tries):
        chosen: List[str] = []
        used = set()
        ok = True
        order = keys.copy()
        rng.shuffle(order)
        for k in order:
            candidates = pools[k].copy()
            rng.shuffle(candidates)
            pick = next((x for x in candidates if x not in used), None)
            if pick is None:
                ok = False
                break
            chosen.append(pick)
            used.add(pick)
        if ok and len(chosen) == len(keys):
            return chosen
    # fallback: sample from union, preserving set size
    union = sorted(set(x for v in pools.values() for x in v))
    if len(union) < len(keys):
        raise RuntimeError('Matched background union is smaller than the focus set')
    return rng.choice(union, size=len(keys), replace=False).tolist()


def cluster_metrics(ids: Sequence[str], ann_idx: pd.DataFrame) -> Dict[str, float]:
    sub = ann_idx.loc[list(ids)]
    contigs = sub['contig'].astype(str).tolist()
    starts = pd.to_numeric(sub['start']).astype(int).to_numpy()
    ends = pd.to_numeric(sub['end']).astype(int).to_numpy()
    centers = (starts + ends) / 2
    same_pairs = 0
    within20 = 0
    within50 = 0
    within100 = 0
    distances = []
    for i in range(len(ids)):
        for j in range(i + 1, len(ids)):
            if contigs[i] == contigs[j]:
                same_pairs += 1
                d = abs(centers[i] - centers[j])
                distances.append(d)
                within20 += int(d <= 20_000)
                within50 += int(d <= 50_000)
                within100 += int(d <= 100_000)
    counts = Counter(contigs)
    return {
        'unique_contigs': float(len(counts)),
        'max_genes_one_contig': float(max(counts.values()) if counts else 0),
        'same_contig_pairs': float(same_pairs),
        'pairs_within_20kb': float(within20),
        'pairs_within_50kb': float(within50),
        'pairs_within_100kb': float(within100),
        'median_same_contig_distance': float(np.median(distances)) if distances else float('nan'),
    }


def permutation_p_beta(y_res: np.ndarray, x_res: np.ndarray, groups: Sequence[str], n_perm: int, rng: np.random.Generator) -> Tuple[float, float]:
    denom = float(np.dot(x_res, x_res))
    if denom <= 0:
        return float('nan'), float('nan')
    obs = float(np.dot(x_res, y_res) / denom)
    idx_by_group = [np.where(np.asarray(groups) == g)[0] for g in pd.unique(groups)]
    extreme = 0
    for _ in range(n_perm):
        yp = y_res.copy()
        for idx in idx_by_group:
            yp[idx] = rng.permutation(yp[idx])
        beta = float(np.dot(x_res, yp) / denom)
        if abs(beta) >= abs(obs) - 1e-12:
            extreme += 1
    return obs, (extreme + 1) / (n_perm + 1)


def loo_stability(y: np.ndarray, x: np.ndarray, X0: np.ndarray) -> float:
    signs = []
    full_x = residualize_matrix(x[:, None], X0).ravel()
    full_y = residualize_matrix(y[:, None], X0).ravel()
    denom = np.dot(full_x, full_x)
    if denom <= 0:
        return float('nan')
    full_sign = np.sign(np.dot(full_x, full_y) / denom)
    for i in range(len(y)):
        keep = np.ones(len(y), dtype=bool)
        keep[i] = False
        xr = residualize_matrix(x[keep, None], X0[keep, :]).ravel()
        yr = residualize_matrix(y[keep, None], X0[keep, :]).ravel()
        d = np.dot(xr, xr)
        if d <= 0:
            continue
        signs.append(np.sign(np.dot(xr, yr) / d) == full_sign)
    return float(np.mean(signs)) if signs else float('nan')


def true_tf_mask(text: pd.Series) -> pd.Series:
    include = re.compile(
        r'transcription factor|transcriptional regulator|mads[- ]box|myb[- ]like dna-binding|homeobox|\bste12\b|'
        r'\bpro-?1\b|cutinase transcription factor|zn\(2\)[- ]?cys\(6\)|zn2cys6|zinc cluster|\bbzip\b|'
        r'\bgata\b|forkhead|velvet', re.I
    )
    exclude = re.compile(r'proteasome|protease regulatory|regulatory subunit|mTOR|cytoskeleton|transposition|IWS1', re.I)
    return text.fillna('').map(lambda s: bool(include.search(str(s))) and not bool(exclude.search(str(s))))



def main() -> None:
    parser = argparse.ArgumentParser(description="Analyze genomic dispersion, residual coexpression, promoter k-mers, and candidate regulators.")
    parser.add_argument("--annotation", type=Path, required=True)
    parser.add_argument("--focus-genes", type=Path, required=True)
    parser.add_argument("--module-scores", type=Path, required=True)
    parser.add_argument("--logcpm", type=Path, required=True)
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--samples", type=Path, required=True)
    parser.add_argument("--interaction-de", type=Path, required=True)
    parser.add_argument("--genome", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260716)
    parser.add_argument("--matched-iterations", type=int, default=20000)
    parser.add_argument("--regulator-permutations", type=int, default=50000)
    parser.add_argument("--promoter-bp", type=int, default=1000)
    args = parser.parse_args()

    out = args.output_dir
    tables = out / "tables"
    sequences = out / "sequences"
    for directory in (out, tables, sequences):
        directory.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(args.seed)
    ann = pd.read_csv(args.annotation, sep="\t", dtype=str).fillna("")
    focus = pd.read_csv(args.focus_genes, sep="\t", dtype=str).fillna("")
    samples = pd.read_csv(args.samples, sep="\t")
    logcpm = pd.read_csv(args.logcpm, sep="\t", index_col=0)
    counts = pd.read_csv(args.counts, sep="\t", index_col=0)
    interaction = pd.read_csv(args.interaction_de, sep="\t")
    scores = pd.read_csv(args.module_scores, sep="\t")
    scores = scores[scores["score_type"].eq("zmean")].copy()

    if "group" not in samples:
        samples["group"] = samples["genotype"].astype(str) + "_" + samples["water"].astype(str)
    sample_ids = [sample for sample in samples["sample_id"].astype(str) if sample in logcpm.columns]
    samples = samples[samples["sample_id"].isin(sample_ids)].copy()
    logcpm = logcpm.loc[:, sample_ids]
    counts = counts.reindex(columns=sample_ids).fillna(0)
    scores = samples[["sample_id", "genotype", "water", "group", "ft_unique_rpm_input"]].merge(
        scores[["sample_id", "score"]], on="sample_id", how="inner", validate="one_to_one"
    ).set_index("sample_id").loc[sample_ids].reset_index()

    focus_ids = [gene for gene in focus["gene_id"].astype(str) if gene in logcpm.index]
    if len(focus_ids) < 8:
        raise RuntimeError(f"Only {len(focus_ids)} focus genes are present in logCPM")
    focus_set = set(focus_ids)

    seqs = read_fasta(args.genome)
    ann["gene_id"] = ann["gene_id"].astype(str)
    if "annotation_text" not in ann:
        parts = [ann[column].astype(str) for column in ("product", "description_text", "description", "Name") if column in ann]
        ann["annotation_text"] = pd.concat(parts, axis=1).agg(" ".join, axis=1) if parts else ""
    ann = ann.drop_duplicates("gene_id").set_index("gene_id", drop=False)
    usable_ids = [
        gene for gene in ann.index
        if gene in logcpm.index and str(ann.loc[gene, "contig"]) and str(ann.loc[gene, "start"]) and str(ann.loc[gene, "end"])
    ]

    promoter_rows = []
    promoters: Dict[str, str] = {}
    for gene in usable_ids:
        row = ann.loc[gene]
        promoter, fasta_contig, left, right, status = extract_promoter(row, seqs, args.promoter_bp)
        if promoter:
            promoters[gene] = promoter
        promoter_rows.append({
            "gene_id": gene,
            "annotation_text": row.get("annotation_text", ""),
            "contig": row["contig"],
            "fasta_contig": fasta_contig,
            "gene_start": int(row["start"]),
            "gene_end": int(row["end"]),
            "strand": row["strand"],
            "promoter_start": left,
            "promoter_end": right,
            "promoter_length": len(promoter),
            "promoter_gc": ((promoter.count("G") + promoter.count("C")) / len(promoter)) if promoter else np.nan,
            "promoter_n_fraction": (promoter.count("N") / len(promoter)) if promoter else np.nan,
            "promoter_status": status,
            "is_focus": gene in focus_set,
        })
    promoter_qc = pd.DataFrame(promoter_rows)
    promoter_qc.to_csv(tables / "promoter_qc.tsv", sep="\t", index=False, lineterminator="\n")
    valid = promoter_qc[promoter_qc["promoter_status"].eq("ok")].copy()
    focus_valid = [gene for gene in focus_ids if gene in set(valid["gene_id"])]
    if len(focus_valid) < 8:
        raise RuntimeError(f"Only {len(focus_valid)} focus promoters passed sequence QC")

    mean_logcpm = logcpm.mean(axis=1)
    detection = (counts.reindex(logcpm.index).fillna(0) >= 5).mean(axis=1)
    features = valid.set_index("gene_id").copy()
    features["mean_logcpm"] = mean_logcpm.reindex(features.index)
    features["detection_fraction"] = detection.reindex(features.index)
    features = features.replace([np.inf, -np.inf], np.nan).dropna(
        subset=["mean_logcpm", "detection_fraction", "promoter_gc", "promoter_length"]
    )
    background_ids = [gene for gene in features.index if gene not in focus_set]
    feature_columns = ["mean_logcpm", "detection_fraction", "promoter_gc", "promoter_length"]
    mean = features.loc[background_ids, feature_columns].mean()
    sd = features.loc[background_ids, feature_columns].std(ddof=1).replace(0, 1)
    scaled = (features[feature_columns] - mean) / sd
    pools: Dict[str, List[str]] = {}
    pool_rows = []
    for gene in focus_valid:
        distances = ((scaled.loc[background_ids] - scaled.loc[gene]) ** 2).sum(axis=1).pow(0.5).sort_values()
        pool = distances.head(min(60, len(distances))).index.tolist()
        pools[gene] = pool
        for rank, matched in enumerate(pool, 1):
            pool_rows.append({"focus_gene": gene, "matched_gene": matched, "rank": rank, "distance": float(distances.loc[matched])})
    pool_table = pd.DataFrame(pool_rows)
    pool_table.to_csv(tables / "matched_background_pools.tsv", sep="\t", index=False, lineterminator="\n")
    background_union = sorted(set(pool_table["matched_gene"]))
    features.loc[focus_valid].reset_index().to_csv(tables / "focus_promoters.tsv", sep="\t", index=False, lineterminator="\n")
    features.loc[background_union].reset_index().to_csv(tables / "matched_background_promoters.tsv", sep="\t", index=False, lineterminator="\n")
    write_fasta([(gene, promoters[gene]) for gene in focus_valid], sequences / "focus_promoters.fasta")
    write_fasta([(gene, promoters[gene]) for gene in background_union], sequences / "matched_background_promoters.fasta")

    matched_sets = [draw_matched_set(pools, rng) for _ in range(args.matched_iterations)]
    annotation_index = ann.loc[[gene for gene in ann.index if gene in features.index]]
    observed_cluster = cluster_metrics(focus_valid, annotation_index)
    null_cluster = defaultdict(list)
    for matched in matched_sets:
        metrics = cluster_metrics(matched, annotation_index)
        for key, value in metrics.items():
            null_cluster[key].append(value)
    cluster_rows = []
    for metric, observed in observed_cluster.items():
        null = np.asarray(null_cluster[metric], dtype=float)
        if metric in ("unique_contigs", "median_same_contig_distance"):
            finite = null[np.isfinite(null)]
            p = (1 + np.sum(finite <= observed)) / (1 + len(finite)) if np.isfinite(observed) and len(finite) else np.nan
        else:
            p = (1 + np.sum(null >= observed)) / (1 + len(null))
        cluster_rows.append({
            "metric": metric, "observed": observed, "null_mean": float(np.nanmean(null)),
            "null_median": float(np.nanmedian(null)), "empirical_p": p,
        })
    pd.DataFrame(cluster_rows).to_csv(tables / "genomic_clustering.tsv", sep="\t", index=False, lineterminator="\n", na_rep="NA")

    design = design_matrix(scores)
    all_expression = logcpm.loc[features.index, sample_ids].T.to_numpy(float)
    residuals = residualize_matrix(all_expression, design)
    gene_column = {gene: index for index, gene in enumerate(features.index)}
    focus_columns = [gene_column[gene] for gene in focus_valid]
    focus_residuals = residuals[:, focus_columns]
    mean_corr, median_corr, proportion_positive, corr = upper_mean_corr(focus_residuals)
    standardized = focus_residuals - focus_residuals.mean(axis=0, keepdims=True)
    sds = focus_residuals.std(axis=0, ddof=1, keepdims=True)
    sds[sds == 0] = 1
    standardized /= sds
    singular = np.linalg.svd(standardized, full_matrices=False, compute_uv=False)
    pc1_variance = float(singular[0] ** 2 / np.sum(singular ** 2))
    null_corr = []
    null_pc1 = []
    for matched in matched_sets:
        columns = [gene_column[gene] for gene in matched]
        value, _, _, _ = upper_mean_corr(residuals[:, columns])
        null_corr.append(value)
        matrix = residuals[:, columns]
        matrix = matrix - matrix.mean(axis=0, keepdims=True)
        sd_matrix = matrix.std(axis=0, ddof=1, keepdims=True)
        sd_matrix[sd_matrix == 0] = 1
        matrix /= sd_matrix
        sv = np.linalg.svd(matrix, full_matrices=False, compute_uv=False)
        null_pc1.append(float(sv[0] ** 2 / np.sum(sv ** 2)))
    corr_p = (1 + np.sum(np.asarray(null_corr) >= mean_corr)) / (1 + len(null_corr))
    pc1_p = (1 + np.sum(np.asarray(null_pc1) >= pc1_variance)) / (1 + len(null_pc1))
    pd.DataFrame([{
        "focus_genes": len(focus_valid),
        "mean_residual_pairwise_correlation": mean_corr,
        "median_residual_pairwise_correlation": median_corr,
        "proportion_positive_pairs": proportion_positive,
        "pc1_variance_explained": pc1_variance,
        "matched_p_mean_correlation": corr_p,
        "matched_p_pc1": pc1_p,
        "matched_iterations": len(null_corr),
    }]).to_csv(tables / "coexpression_coherence.tsv", sep="\t", index=False, lineterminator="\n")
    pd.DataFrame(corr, index=focus_valid, columns=focus_valid).to_csv(tables / "residual_correlation_matrix.tsv", sep="\t", lineterminator="\n")

    motif_presence: Dict[str, set] = {}
    for gene in sorted(set(focus_valid + background_union)):
        sequence = promoters[gene]
        present = set()
        for length in (6, 7):
            for start in range(0, len(sequence) - length + 1):
                kmer = sequence[start:start + length]
                if set(kmer) <= set("ACGT"):
                    present.add(canonical_kmer(kmer))
        motif_presence[gene] = present
    focus_counter = Counter(motif for gene in focus_valid for motif in motif_presence[gene])
    background_counter = Counter(motif for gene in background_union for motif in motif_presence[gene])
    motif_rows = []
    for motif, a in focus_counter.items():
        if a < 4:
            continue
        b = len(focus_valid) - a
        c = background_counter.get(motif, 0)
        d = len(background_union) - c
        odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        motif_rows.append({
            "motif": motif, "k": len(motif), "focus_present": a, "focus_total": len(focus_valid),
            "background_present": c, "background_total": len(background_union), "odds_ratio": odds, "fisher_p": p,
        })
    motif_table = pd.DataFrame(motif_rows)
    if not motif_table.empty:
        motif_table["fisher_fdr"] = bh_adjust(motif_table["fisher_p"])
        motif_table = motif_table.sort_values(["fisher_p", "motif"]).reset_index(drop=True)
        top = motif_table.head(50)["motif"].tolist()
        empirical = {}
        for motif in top:
            observed = int(motif_table.loc[motif_table["motif"].eq(motif), "focus_present"].iloc[0])
            exceed = sum(sum(motif in motif_presence.get(gene, set()) for gene in matched) >= observed for matched in matched_sets)
            empirical[motif] = (exceed + 1) / (len(matched_sets) + 1)
        motif_table["matched_empirical_p"] = motif_table["motif"].map(empirical)
        motif_table["robust"] = (
            motif_table["fisher_fdr"].le(0.10)
            & motif_table["matched_empirical_p"].fillna(1).le(0.01)
            & motif_table["focus_present"].ge(5)
        )
    motif_table.to_csv(tables / "promoter_kmer_enrichment.tsv", sep="\t", index=False, lineterminator="\n", na_rep="NA")

    tf_mask = true_tf_mask(ann["annotation_text"])
    candidate_ids = [gene for gene in ann.index[tf_mask] if gene in logcpm.index and gene not in focus_set]
    response = scores["score"].to_numpy(float)
    groups = scores["group"].astype(str).to_numpy()
    base_design = design_matrix(scores)
    response_residual = residualize_matrix(response[:, None], base_design).ravel()
    interaction_index = interaction.set_index("gene_id") if "gene_id" in interaction else pd.DataFrame()
    regulator_rows = []
    for gene in candidate_ids:
        predictor = logcpm.loc[gene, sample_ids].to_numpy(float)
        fit = sm.OLS(response, np.column_stack([base_design, predictor])).fit(cov_type="HC3")
        predictor_residual = residualize_matrix(predictor[:, None], base_design).ravel()
        beta_perm, p_perm = permutation_p_beta(
            response_residual, predictor_residual, groups, args.regulator_permutations, rng
        )
        leave_one_out = loo_stability(response, predictor, base_design)
        interaction_row = interaction_index.loc[gene] if gene in interaction_index.index else None
        regulator_rows.append({
            "gene_id": gene,
            "annotation_text": ann.loc[gene, "annotation_text"],
            "partial_beta": float(fit.params[-1]),
            "partial_se_hc3": float(fit.bse[-1]),
            "partial_p_hc3": float(fit.pvalues[-1]),
            "permutation_beta": beta_perm,
            "permutation_p": p_perm,
            "leave_one_out_sign_stability": leave_one_out,
            "interaction_logFC": float(interaction_row["state_logFC"]) if interaction_row is not None and pd.notna(interaction_row.get("state_logFC")) else np.nan,
            "interaction_FDR": float(interaction_row["state_FDR"]) if interaction_row is not None and pd.notna(interaction_row.get("state_FDR")) else np.nan,
        })
    regulator_table = pd.DataFrame(regulator_rows)
    if not regulator_table.empty:
        regulator_table["partial_fdr"] = bh_adjust(regulator_table["partial_p_hc3"])
        regulator_table["permutation_fdr"] = bh_adjust(regulator_table["permutation_p"])
        regulator_table["robust"] = (
            regulator_table["permutation_fdr"].le(0.10)
            & regulator_table["partial_p_hc3"].le(0.05)
            & regulator_table["leave_one_out_sign_stability"].ge(0.90)
        )
    regulator_table.to_csv(tables / "candidate_regulators.tsv", sep="\t", index=False, lineterminator="\n", na_rep="NA")

    summary = {
        "focus_genes": len(focus_valid),
        "unique_contigs": int(observed_cluster["unique_contigs"]),
        "same_contig_pairs": int(observed_cluster["same_contig_pairs"]),
        "mean_residual_correlation": mean_corr,
        "matched_p_mean_correlation": corr_p,
        "pc1_variance_explained": pc1_variance,
        "matched_p_pc1": pc1_p,
        "robust_motifs": int(motif_table.get("robust", pd.Series(dtype=bool)).fillna(False).sum()),
        "robust_regulators": int(regulator_table.get("robust", pd.Series(dtype=bool)).fillna(False).sum()),
    }
    (tables / "regulatory_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
