#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from itertools import combinations, product
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import scipy.stats as st
import statsmodels.api as sm


GROUP_ORDER = ["wt_dry", "bmr6_dry", "bmr12_dry", "wt_wet", "bmr6_wet", "bmr12_wet"]

CONTRASTS: dict[str, dict[str, float]] = {
    "bmr6_vs_wt_dry": {"bmr6_dry": 1.0, "wt_dry": -1.0},
    "bmr12_vs_wt_dry": {"bmr12_dry": 1.0, "wt_dry": -1.0},
    "bmr6_vs_wt_wet": {"bmr6_wet": 1.0, "wt_wet": -1.0},
    "bmr12_vs_wt_wet": {"bmr12_wet": 1.0, "wt_wet": -1.0},
    "wet_vs_dry_wt": {"wt_wet": 1.0, "wt_dry": -1.0},
    "wet_vs_dry_bmr6": {"bmr6_wet": 1.0, "bmr6_dry": -1.0},
    "wet_vs_dry_bmr12": {"bmr12_wet": 1.0, "bmr12_dry": -1.0},
    "bmr6_water_interaction": {"bmr6_wet": 1.0, "bmr6_dry": -1.0, "wt_wet": -1.0, "wt_dry": 1.0},
    "bmr12_water_interaction": {"bmr12_wet": 1.0, "bmr12_dry": -1.0, "wt_wet": -1.0, "wt_dry": 1.0},
    "bmr12_vs_bmr6_dry": {"bmr12_dry": 1.0, "bmr6_dry": -1.0},
    "bmr12_vs_bmr6_wet": {"bmr12_wet": 1.0, "bmr6_wet": -1.0},
}

KEY_AROMATIC_CONTRASTS = [
    "bmr6_vs_wt_wet",
    "bmr12_vs_bmr6_wet",
    "bmr6_water_interaction",
    "wet_vs_dry_wt",
]


def require(path: Path) -> Path:
    if not path.exists() or (path.is_file() and path.stat().st_size == 0):
        raise FileNotFoundError(f"Required input is missing or empty: {path}")
    return path


def bh_adjust(values: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(values), dtype=float)
    out = np.full(p.shape, np.nan, dtype=float)
    finite = np.isfinite(p)
    if not finite.any():
        return out
    vals = p[finite]
    order = np.argsort(vals)
    ranked = vals[order]
    n = len(ranked)
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    temp = np.empty_like(adjusted)
    temp[order] = adjusted
    out[finite] = temp
    return out


def zscore_rows(df: pd.DataFrame) -> pd.DataFrame:
    arr = df.to_numpy(dtype=float, copy=True)
    means = np.nanmean(arr, axis=1, keepdims=True)
    sds = np.nanstd(arr, axis=1, ddof=0, keepdims=True)
    sds[~np.isfinite(sds) | (sds == 0)] = np.nan
    z = (arr - means) / sds
    return pd.DataFrame(z, index=df.index, columns=df.columns)


def pca1_score(z: pd.DataFrame, orient: pd.Series) -> pd.Series:
    x = z.T.to_numpy(dtype=float, copy=True)
    x = np.nan_to_num(x, nan=0.0)
    x -= x.mean(axis=0, keepdims=True)
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    score = u[:, 0] * s[0] if len(s) else np.zeros(x.shape[0])
    if np.std(score) > 0 and np.std(orient.to_numpy()) > 0:
        corr = np.corrcoef(score, orient.to_numpy(dtype=float))[0, 1]
        if np.isfinite(corr) and corr < 0:
            score = -score
    return pd.Series(score, index=z.columns, name="pc1")


def group_effect(values: pd.Series, sample_table: pd.DataFrame, coeffs: dict[str, float]) -> float:
    merged = sample_table.set_index("sample_id").join(values.rename("score"), how="inner")
    effect = 0.0
    for group, coefficient in coeffs.items():
        sub = merged.loc[merged["group"] == group, "score"]
        if sub.empty:
            return float("nan")
        effect += coefficient * float(sub.mean())
    return effect


def exact_pair_permutation(
    values: pd.Series,
    sample_table: pd.DataFrame,
    group_a: str,
    group_b: str,
) -> tuple[float, float, int]:
    merged = sample_table.set_index("sample_id").join(values.rename("score"), how="inner")
    x = merged.loc[merged["group"] == group_a, "score"].to_numpy(dtype=float)
    y = merged.loc[merged["group"] == group_b, "score"].to_numpy(dtype=float)
    all_values = np.concatenate([x, y])
    n_a = len(x)
    observed = float(np.mean(x) - np.mean(y))
    stats = []
    for idx in combinations(range(len(all_values)), n_a):
        mask = np.zeros(len(all_values), dtype=bool)
        mask[list(idx)] = True
        stats.append(float(all_values[mask].mean() - all_values[~mask].mean()))
    arr = np.asarray(stats)
    p = float(np.mean(np.abs(arr) >= abs(observed) - 1e-12))
    return observed, p, len(arr)


def exact_interaction_permutation(
    values: pd.Series,
    sample_table: pd.DataFrame,
    genotype: str,
) -> tuple[float, float, int]:
    merged = sample_table.set_index("sample_id").join(values.rename("score"), how="inner").reset_index()
    sub = merged[merged["genotype"].isin(["wt", genotype])].copy()

    def wet_minus_dry(frame: pd.DataFrame, geno: str) -> float:
        wet = frame[(frame["genotype"] == geno) & (frame["water"] == "wet")]["score"]
        dry = frame[(frame["genotype"] == geno) & (frame["water"] == "dry")]["score"]
        return float(wet.mean() - dry.mean())

    observed = wet_minus_dry(sub, genotype) - wet_minus_dry(sub, "wt")
    data: dict[str, np.ndarray] = {}
    assignments: dict[str, list[tuple[int, ...]]] = {}
    for geno in ["wt", genotype]:
        frame = sub[sub["genotype"] == geno]
        data[geno] = frame["score"].to_numpy(dtype=float)
        n_wet = int((frame["water"] == "wet").sum())
        assignments[geno] = list(combinations(range(len(frame)), n_wet))

    stats = []
    for wt_idx, g_idx in product(assignments["wt"], assignments[genotype]):
        def diff(arr: np.ndarray, wet_idx: tuple[int, ...]) -> float:
            mask = np.zeros(len(arr), dtype=bool)
            mask[list(wet_idx)] = True
            return float(arr[mask].mean() - arr[~mask].mean())

        stats.append(diff(data[genotype], g_idx) - diff(data["wt"], wt_idx))
    arr = np.asarray(stats)
    p = float(np.mean(np.abs(arr) >= abs(observed) - 1e-12))
    return observed, p, len(arr)


def exact_permutation(
    values: pd.Series,
    sample_table: pd.DataFrame,
    contrast_name: str,
) -> tuple[float, float, int]:
    coeffs = CONTRASTS[contrast_name]
    positive = [g for g, v in coeffs.items() if v == 1]
    negative = [g for g, v in coeffs.items() if v == -1]
    if len(coeffs) == 2 and len(positive) == 1 and len(negative) == 1:
        return exact_pair_permutation(values, sample_table, positive[0], negative[0])
    if contrast_name == "bmr6_water_interaction":
        return exact_interaction_permutation(values, sample_table, "bmr6")
    if contrast_name == "bmr12_water_interaction":
        return exact_interaction_permutation(values, sample_table, "bmr12")
    raise ValueError(f"Unsupported exact contrast: {contrast_name}")


def fit_group_model(
    values: pd.Series,
    sample_table: pd.DataFrame,
    coeffs: dict[str, float],
    burden_adjusted: bool,
) -> dict[str, float]:
    merged = sample_table.set_index("sample_id").join(values.rename("score"), how="inner").reset_index()
    group_matrix = pd.get_dummies(merged["group"], dtype=float).reindex(columns=GROUP_ORDER, fill_value=0.0)
    x_parts = [group_matrix.to_numpy(dtype=float, copy=True)]
    names = list(GROUP_ORDER)
    if burden_adjusted:
        log_ft = np.log10(pd.to_numeric(merged["ft_unique_rpm_input"], errors="coerce").fillna(0).to_numpy() + 1.0)
        x_parts.append(log_ft[:, None])
        names.append("log10_ft_rpm_plus1")
    x = np.column_stack(x_parts)
    y = merged["score"].to_numpy(dtype=float, copy=True)
    fit = sm.OLS(y, x).fit(cov_type="HC3")
    c = np.zeros(x.shape[1], dtype=float)
    for group, coef in coeffs.items():
        c[names.index(group)] = coef
    test = fit.t_test(c)
    ci = np.asarray(test.conf_int(alpha=0.05)).reshape(-1)
    return {
        "estimate": float(np.asarray(test.effect).reshape(-1)[0]),
        "se_hc3": float(np.asarray(test.sd).reshape(-1)[0]),
        "statistic": float(np.asarray(test.tvalue).reshape(-1)[0]),
        "p_hc3": float(np.asarray(test.pvalue).reshape(-1)[0]),
        "ci_low": float(ci[0]),
        "ci_high": float(ci[1]),
        "r_squared": float(fit.rsquared),
        "n": int(fit.nobs),
    }


def make_quantile_bins(feature: pd.Series, q: int) -> pd.Series:
    ranks = feature.rank(method="first")
    try:
        return pd.qcut(ranks, q=q, labels=False, duplicates="drop").astype(int)
    except ValueError:
        return pd.Series(np.zeros(len(feature), dtype=int), index=feature.index)


def choose_matched_controls(
    module_genes: list[str],
    universe: list[str],
    mean_logcpm: pd.Series,
    detection: pd.Series,
    controls_per_gene: int,
) -> tuple[list[str], pd.DataFrame]:
    feature = pd.DataFrame({
        "mean_logcpm": mean_logcpm,
        "detection": detection,
    }).loc[universe]
    scaled = (feature - feature.mean()) / feature.std(ddof=0).replace(0, 1)
    candidates = [g for g in universe if g not in module_genes]
    selected: list[str] = []
    rows: list[dict[str, object]] = []
    used: set[str] = set()
    for gene in module_genes:
        dist = ((scaled.loc[candidates] - scaled.loc[gene]) ** 2).sum(axis=1).sort_values()
        taken = 0
        for candidate, value in dist.items():
            if candidate in used:
                continue
            used.add(candidate)
            selected.append(candidate)
            rows.append({
                "module_gene": gene,
                "control_gene": candidate,
                "squared_feature_distance": float(value),
                "module_mean_logcpm": float(mean_logcpm.loc[gene]),
                "control_mean_logcpm": float(mean_logcpm.loc[candidate]),
                "module_detection": float(detection.loc[gene]),
                "control_detection": float(detection.loc[candidate]),
            })
            taken += 1
            if taken >= controls_per_gene:
                break
    return selected, pd.DataFrame(rows)


def contrast_direction_counts(
    de_table: pd.DataFrame,
    members: list[str],
    logfc_col: str,
    expected_direction: str,
) -> dict[str, float | int]:
    sub = de_table[de_table["gene_id"].isin(members)].copy()
    vals = pd.to_numeric(sub[logfc_col], errors="coerce").dropna()
    positive = int((vals > 0).sum())
    negative = int((vals < 0).sum())
    n = int(len(vals))
    matching = positive if expected_direction.lower() == "up" else negative
    p_binom = float(st.binomtest(matching, n, p=0.5, alternative="greater").pvalue) if n else float("nan")
    return {
        "n_genes": n,
        "positive": positive,
        "negative": negative,
        "direction_matching": matching,
        "direction_matching_fraction": matching / n if n else float("nan"),
        "median_logfc": float(np.median(vals)) if n else float("nan"),
        "mean_logfc": float(np.mean(vals)) if n else float("nan"),
        "descriptive_binomial_p": p_binom,
    }


def matched_set_empirical(
    z_matrix: pd.DataFrame,
    module_genes: list[str],
    sample_table: pd.DataFrame,
    contrast_name: str,
    mean_logcpm: pd.Series,
    detection: pd.Series,
    iterations: int,
    bins: int,
    seed: int,
) -> dict[str, float | int]:
    universe = list(z_matrix.index)
    module_set = set(module_genes)
    mean_bin = make_quantile_bins(mean_logcpm.loc[universe], bins)
    detect_bin = make_quantile_bins(detection.loc[universe], bins)
    labels = mean_bin.astype(str) + "_" + detect_bin.astype(str)
    module_counts = labels.loc[module_genes].value_counts().to_dict()
    pools: dict[str, np.ndarray] = {}
    index_map = {g: i for i, g in enumerate(universe)}
    for label, count in module_counts.items():
        genes = [g for g in labels.index[(labels == label)] if g not in module_set]
        if len(genes) < count:
            genes = [g for g in universe if g not in module_set]
        pools[label] = np.asarray([index_map[g] for g in genes], dtype=int)

    z_arr = z_matrix.to_numpy(dtype=float, copy=True)
    sample_ids = list(z_matrix.columns)
    sample_pos = {sid: i for i, sid in enumerate(sample_ids)}
    group_positions = {
        group: np.asarray([sample_pos[sid] for sid in sample_table.loc[sample_table["group"] == group, "sample_id"] if sid in sample_pos], dtype=int)
        for group in GROUP_ORDER
    }
    coeffs = CONTRASTS[contrast_name]

    observed_score = z_matrix.loc[module_genes].mean(axis=0)
    observed = group_effect(observed_score, sample_table, coeffs)

    rng = np.random.default_rng(seed)
    null = np.empty(iterations, dtype=float)
    for i in range(iterations):
        selected: list[int] = []
        for label, count in module_counts.items():
            pool = pools[label]
            replace = len(pool) < count
            selected.extend(rng.choice(pool, size=count, replace=replace).tolist())
        score = np.nanmean(z_arr[np.asarray(selected, dtype=int), :], axis=0)
        effect = 0.0
        for group, coefficient in coeffs.items():
            effect += coefficient * float(np.mean(score[group_positions[group]]))
        null[i] = effect
    empirical_p = float((1 + np.sum(np.abs(null) >= abs(observed))) / (iterations + 1))
    return {
        "observed_effect": float(observed),
        "empirical_p_two_sided": empirical_p,
        "iterations": int(iterations),
        "null_mean": float(np.mean(null)),
        "null_sd": float(np.std(null, ddof=1)),
        "null_q025": float(np.quantile(null, 0.025)),
        "null_q975": float(np.quantile(null, 0.975)),
        "null": null,
    }


def markdown_table(df: pd.DataFrame, columns: list[str], max_rows: int = 30) -> str:
    if df.empty:
        return "No rows met the reporting criterion."
    view = df.loc[:, [c for c in columns if c in df.columns]].head(max_rows).copy()
    for col in view.columns:
        if pd.api.types.is_float_dtype(view[col]):
            view[col] = view[col].map(lambda x: "NA" if not np.isfinite(x) else f"{x:.4g}")
    header = "| " + " | ".join(view.columns) + " |"
    sep = "| " + " | ".join(["---"] * len(view.columns)) + " |"
    rows = ["| " + " | ".join(str(v).replace("|", "/") for v in row) + " |" for row in view.itertuples(index=False, name=None)]
    return "\n".join([header, sep, *rows])



def main() -> None:
    parser = argparse.ArgumentParser(description="Audit the aromatic/phenolic-processing gene set.")
    parser.add_argument("--samples-primary", type=Path, required=True)
    parser.add_argument("--samples-all", type=Path, required=True)
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--logcpm-primary", type=Path, required=True)
    parser.add_argument("--logcpm-all", type=Path, required=True)
    parser.add_argument("--focus-genes", type=Path, required=True)
    parser.add_argument("--de-contrast", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260715)
    parser.add_argument("--matched-iterations", type=int, default=10000)
    parser.add_argument("--match-bins", type=int, default=5)
    parser.add_argument("--controls-per-gene", type=int, default=3)
    parser.add_argument("--focus-module", default="aromatic_phenolic_metabolism")
    args = parser.parse_args()

    samples = pd.read_csv(args.samples_primary, sep="\t")
    samples_all = pd.read_csv(args.samples_all, sep="\t")
    counts = pd.read_csv(args.counts, sep="\t", index_col=0)
    logcpm = pd.read_csv(args.logcpm_primary, sep="\t", index_col=0)
    logcpm_all = pd.read_csv(args.logcpm_all, sep="\t", index_col=0)
    focus_table = pd.read_csv(args.focus_genes, sep="\t", dtype=str).fillna("")
    focus = [gene for gene in focus_table["gene_id"].astype(str) if gene in logcpm.index]
    if len(focus) < 5:
        raise RuntimeError(f"Too few focus genes: {len(focus)}")

    sample_ids = [sample for sample in samples["sample_id"] if sample in logcpm.columns]
    sample_ids_all = [sample for sample in samples_all["sample_id"] if sample in logcpm_all.columns]
    samples = samples[samples["sample_id"].isin(sample_ids)].copy()
    samples_all = samples_all[samples_all["sample_id"].isin(sample_ids_all)].copy()
    if "group" not in samples:
        samples["group"] = samples["genotype"].astype(str) + "_" + samples["water"].astype(str)
    if "group" not in samples_all:
        samples_all["group"] = samples_all["genotype"].astype(str) + "_" + samples_all["water"].astype(str)
    logcpm = logcpm.loc[:, sample_ids]
    logcpm_all = logcpm_all.loc[:, sample_ids_all]
    counts_primary = counts.reindex(index=logcpm.index, columns=sample_ids).fillna(0)
    mean_logcpm = logcpm.mean(axis=1)
    detection = (counts_primary > 0).mean(axis=1)

    expression = logcpm.loc[focus]
    z = zscore_rows(expression)
    zmean = z.mean(axis=0)
    rankmean = expression.rank(axis=1, method="average", pct=True).mean(axis=0) - 0.5
    pc1 = pca1_score(z, zmean)
    controls, control_map = choose_matched_controls(
        focus, list(logcpm.index), mean_logcpm, detection, args.controls_per_gene
    )
    matched_logratio = expression.mean(axis=0) - logcpm.loc[controls].mean(axis=0)
    scores = {"zmean": zmean, "rankmean": rankmean, "pc1": pc1, "matched_logratio": matched_logratio}

    args.output_dir.mkdir(parents=True, exist_ok=True)
    control_map.to_csv(args.output_dir / "matched_controls.tsv", sep="\t", index=False, lineterminator="\n")
    score_rows = []
    test_rows = []
    for score_name, values in scores.items():
        joined = samples.set_index("sample_id").join(values.rename("score"), how="inner").reset_index()
        joined["module"] = args.focus_module
        joined["score_type"] = score_name
        score_rows.append(joined)
        for contrast, coefficients in CONTRASTS.items():
            effect, p_exact, permutations = exact_permutation(values, samples, contrast)
            unadjusted = fit_group_model(values, samples, coefficients, burden_adjusted=False)
            adjusted = fit_group_model(values, samples, coefficients, burden_adjusted=True)
            row = {
                "score_type": score_name, "contrast": contrast, "effect": effect,
                "p_exact": p_exact, "exact_permutations": permutations,
            }
            row.update({f"unadjusted_{key}": value for key, value in unadjusted.items()})
            row.update({f"burden_adjusted_{key}": value for key, value in adjusted.items()})
            test_rows.append(row)
    pd.concat(score_rows, ignore_index=True).to_csv(args.output_dir / "module_scores.tsv", sep="\t", index=False, lineterminator="\n")
    tests = pd.DataFrame(test_rows)
    tests["p_exact_fdr"] = bh_adjust(tests["p_exact"])
    tests.to_csv(args.output_dir / "module_tests.tsv", sep="\t", index=False, lineterminator="\n")

    z_all = zscore_rows(logcpm_all.loc[focus]).mean(axis=0)
    all_rows = []
    for contrast in KEY_AROMATIC_CONTRASTS:
        effect, p_exact, permutations = exact_permutation(z_all, samples_all, contrast)
        all_rows.append({"contrast": contrast, "effect": effect, "p_exact": p_exact, "exact_permutations": permutations})
    all_sensitivity = pd.DataFrame(all_rows)
    all_sensitivity.to_csv(args.output_dir / "all_sample_sensitivity.tsv", sep="\t", index=False, lineterminator="\n")

    leave_sample = []
    involved = {
        "bmr6_vs_wt_wet": ["bmr6_wet", "wt_wet"],
        "bmr6_water_interaction": ["bmr6_wet", "bmr6_dry", "wt_wet", "wt_dry"],
    }
    for contrast, groups in involved.items():
        for omitted in samples.loc[samples["group"].isin(groups), "sample_id"]:
            keep = zmean.index[zmean.index != omitted]
            effect, p_exact, permutations = exact_permutation(
                zmean.loc[keep], samples[samples["sample_id"].isin(keep)], contrast
            )
            leave_sample.append({
                "contrast": contrast, "omitted_sample": omitted, "effect": effect,
                "p_exact": p_exact, "exact_permutations": permutations,
                "expected_sign": effect > 0,
            })
    leave_sample_table = pd.DataFrame(leave_sample)
    leave_sample_table.to_csv(args.output_dir / "leave_one_sample.tsv", sep="\t", index=False, lineterminator="\n")

    leave_gene = []
    product_map = focus_table.set_index("gene_id").get("product", pd.Series(dtype=str)).to_dict()
    for omitted in focus:
        retained = [gene for gene in focus if gene != omitted]
        score = zscore_rows(logcpm.loc[retained]).mean(axis=0)
        for contrast in ("bmr6_vs_wt_wet", "bmr12_vs_bmr6_wet", "bmr6_water_interaction"):
            effect, p_exact, permutations = exact_permutation(score, samples, contrast)
            expected = effect < 0 if contrast == "bmr12_vs_bmr6_wet" else effect > 0
            leave_gene.append({
                "contrast": contrast, "omitted_gene": omitted,
                "omitted_product": product_map.get(omitted, ""), "effect": effect,
                "p_exact": p_exact, "exact_permutations": permutations,
                "expected_sign": expected,
            })
    leave_gene_table = pd.DataFrame(leave_gene)
    leave_gene_table.to_csv(args.output_dir / "leave_one_gene.tsv", sep="\t", index=False, lineterminator="\n")

    empirical = matched_set_empirical(
        zscore_rows(logcpm), focus, samples, "bmr6_vs_wt_wet",
        mean_logcpm, detection, args.matched_iterations, args.match_bins, args.seed,
    )
    empirical.pop("null", None)
    pd.DataFrame([empirical]).to_csv(args.output_dir / "matched_set_test.tsv", sep="\t", index=False, lineterminator="\n")

    de = pd.read_csv(args.de_contrast, sep="\t")
    logfc = pd.to_numeric(de.set_index("gene_id")["state_logFC"], errors="coerce").reindex(focus)
    direction_fraction = float((logfc > 0).mean())
    direction = pd.DataFrame({"gene_id": focus, "logFC": logfc.values, "positive": logfc.values > 0})
    direction.to_csv(args.output_dir / "gene_direction.tsv", sep="\t", index=False, lineterminator="\n")

    def test_value(score_type: str, contrast: str, column: str) -> float:
        row = tests[(tests["score_type"] == score_type) & (tests["contrast"] == contrast)]
        return float(row.iloc[0][column])

    criteria = {
        "zmean_pair": test_value("zmean", "bmr6_vs_wt_wet", "p_exact") <= 0.05,
        "rank_pair": test_value("rankmean", "bmr6_vs_wt_wet", "p_exact") <= 0.05,
        "logratio_pair": test_value("matched_logratio", "bmr6_vs_wt_wet", "p_exact") <= 0.05,
        "interaction": test_value("zmean", "bmr6_water_interaction", "p_exact") <= 0.05,
        "burden_adjusted_interaction": test_value("zmean", "bmr6_water_interaction", "burden_adjusted_p_hc3") <= 0.05,
        "all_pair": float(all_sensitivity.loc[all_sensitivity["contrast"] == "bmr6_vs_wt_wet", "p_exact"].iloc[0]) <= 0.05,
        "all_interaction": float(all_sensitivity.loc[all_sensitivity["contrast"] == "bmr6_water_interaction", "p_exact"].iloc[0]) <= 0.05,
        "gene_direction": direction_fraction >= 0.8,
        "leave_sample": bool(leave_sample_table.loc[leave_sample_table["contrast"] == "bmr6_vs_wt_wet", "expected_sign"].all()),
        "leave_gene": bool(leave_gene_table.loc[leave_gene_table["contrast"] == "bmr6_vs_wt_wet", "expected_sign"].all()),
        "matched_set": float(empirical["empirical_p_two_sided"]) <= 0.05,
    }
    summary = {
        "focus_genes": len(focus),
        "criteria_passed": int(sum(criteria.values())),
        "criteria_total": len(criteria),
        "criteria": criteria,
        "bmr6_vs_wt_wet_p_exact": test_value("zmean", "bmr6_vs_wt_wet", "p_exact"),
        "bmr6_water_interaction_p_exact": test_value("zmean", "bmr6_water_interaction", "p_exact"),
        "burden_adjusted_interaction_p": test_value("zmean", "bmr6_water_interaction", "burden_adjusted_p_hc3"),
        "matched_set_p": float(empirical["empirical_p_two_sided"]),
    }
    (args.output_dir / "module_audit.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
