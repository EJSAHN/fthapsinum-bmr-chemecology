#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from itertools import combinations, product
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf


def exact_two_group(values: np.ndarray, labels: np.ndarray, first: str, second: str) -> tuple[float, float, int]:
    keep = np.isin(labels, [first, second])
    y = values[keep]
    groups = labels[keep]
    first_n = int((groups == first).sum())
    observed = float(y[groups == first].mean() - y[groups == second].mean())
    statistics = []
    for indexes in combinations(range(len(y)), first_n):
        mask = np.zeros(len(y), dtype=bool)
        mask[list(indexes)] = True
        statistics.append(float(y[mask].mean() - y[~mask].mean()))
    array = np.asarray(statistics)
    p_value = float(np.mean(np.abs(array) >= abs(observed) - 1e-12))
    return observed, p_value, len(array)


def exact_interaction(frame: pd.DataFrame) -> tuple[float, float, int]:
    def water_difference(data: pd.DataFrame) -> float:
        return float(data.loc[data["water"].eq("wet"), "score"].mean() - data.loc[data["water"].eq("dry"), "score"].mean())

    wt = frame[frame["genotype"].eq("wt")].copy()
    bmr12 = frame[frame["genotype"].eq("bmr12")].copy()
    observed = water_difference(bmr12) - water_difference(wt)
    arrays = {"wt": wt["score"].to_numpy(float), "bmr12": bmr12["score"].to_numpy(float)}
    wet_counts = {"wt": int(wt["water"].eq("wet").sum()), "bmr12": int(bmr12["water"].eq("wet").sum())}
    statistics = []
    for wt_indexes, bmr_indexes in product(
        combinations(range(len(arrays["wt"])), wet_counts["wt"]),
        combinations(range(len(arrays["bmr12"])), wet_counts["bmr12"]),
    ):
        def difference(values: np.ndarray, indexes: tuple[int, ...]) -> float:
            mask = np.zeros(len(values), dtype=bool)
            mask[list(indexes)] = True
            return float(values[mask].mean() - values[~mask].mean())
        statistics.append(difference(arrays["bmr12"], bmr_indexes) - difference(arrays["wt"], wt_indexes))
    array = np.asarray(statistics)
    return observed, float(np.mean(np.abs(array) >= abs(observed) - 1e-12)), len(array)


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize the fixed fungal module in the 13-day context.")
    parser.add_argument("--sample-summary", type=Path, required=True)
    parser.add_argument("--state-logcpm", type=Path, required=True)
    parser.add_argument("--focus-genes", type=Path, required=True)
    parser.add_argument("--state-camera", type=Path, required=True)
    parser.add_argument("--tissue-camera", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    samples = pd.read_csv(args.sample_summary, sep="\t")
    logcpm = pd.read_csv(args.state_logcpm, sep="\t").set_index("gene_id")
    focus = pd.read_csv(args.focus_genes, sep="\t", dtype=str)["gene_id"].astype(str).tolist()
    primary_flag = samples["primary_fus_analysis"].astype(str).str.lower().isin({"1", "true", "yes"})
    primary = samples.loc[primary_flag].copy()
    focus_ids = [gene for gene in focus if gene in logcpm.index]
    if len(focus_ids) < 5:
        raise RuntimeError(f"Only {len(focus_ids)} focus genes are represented in the state model")

    expression = logcpm.loc[focus_ids, primary["sample_id"]].T.astype(float)
    standardized = (expression - expression.mean(axis=0)) / expression.std(axis=0, ddof=1).replace(0, np.nan)
    primary["score"] = standardized.mean(axis=1, skipna=True).reindex(primary["sample_id"]).to_numpy()
    primary["group"] = primary["genotype"].astype(str) + "_" + primary["water"].astype(str)
    primary.to_csv(args.output_dir / "module_scores.tsv", sep="\t", index=False, lineterminator="\n", na_rep="NA")

    rows = []
    labels = primary["group"].to_numpy(str)
    values = primary["score"].to_numpy(float)
    pairs = [
        ("bmr12_vs_wt_dry", "bmr12_dry", "wt_dry"),
        ("bmr12_vs_wt_wet", "bmr12_wet", "wt_wet"),
        ("wet_vs_dry_wt", "wt_wet", "wt_dry"),
        ("wet_vs_dry_bmr12", "bmr12_wet", "bmr12_dry"),
    ]
    for name, first, second in pairs:
        effect, p_value, permutations = exact_two_group(values, labels, first, second)
        rows.append({"contrast": name, "effect": effect, "p_exact": p_value, "permutations": permutations})
    effect, p_value, permutations = exact_interaction(primary)
    rows.append({"contrast": "bmr12_water_interaction", "effect": effect, "p_exact": p_value, "permutations": permutations})
    tests = pd.DataFrame(rows)
    tests.to_csv(args.output_dir / "module_tests.tsv", sep="\t", index=False, lineterminator="\n")

    primary["log_target_rpm"] = np.log1p(pd.to_numeric(primary["target_unique_rpm"], errors="coerce").fillna(0))
    fit = smf.ols("score ~ C(genotype) * C(water) + log_target_rpm", data=primary).fit(cov_type="HC3")
    confidence = fit.conf_int()
    coefficients = pd.DataFrame({
        "term": fit.params.index,
        "estimate": fit.params.values,
        "p_hc3": fit.pvalues.values,
        "ci_low": confidence[0].values,
        "ci_high": confidence[1].values,
    })
    coefficients.to_csv(args.output_dir / "burden_adjusted_model.tsv", sep="\t", index=False, lineterminator="\n")

    camera = pd.concat([
        pd.read_csv(args.state_camera, sep="\t").assign(scale="state"),
        pd.read_csv(args.tissue_camera, sep="\t").assign(scale="tissue"),
    ], ignore_index=True)
    camera.to_csv(args.output_dir / "camera.tsv", sep="\t", index=False, lineterminator="\n")

    lookup = tests.set_index("contrast")
    median_target_reads = float(pd.to_numeric(primary["target_unique"], errors="coerce").median())
    minimum_group = int(primary.groupby(["genotype", "water"]).size().min())
    if median_target_reads < 5000 or minimum_group < 2:
        verdict = "INSUFFICIENT_CONTEXT_SIGNAL"
    elif lookup.loc["wet_vs_dry_bmr12", "effect"] > 0 and lookup.loc["wet_vs_dry_bmr12", "p_exact"] <= 0.05:
        verdict = "LATE_BMR12_AROMATIC_PROGRAM"
    elif (
        lookup.loc["wet_vs_dry_wt", "effect"] > 0 and lookup.loc["wet_vs_dry_wt", "p_exact"] <= 0.05
    ):
        verdict = "LATE_AROMATIC_PROGRAM_GENERALIZES"
    elif (
        lookup.loc["wet_vs_dry_wt", "effect"] <= 0.25
        and lookup.loc["wet_vs_dry_bmr12", "effect"] <= 0.25
        and lookup.loc["wet_vs_dry_wt", "p_exact"] > 0.05
        and lookup.loc["wet_vs_dry_bmr12", "p_exact"] > 0.05
    ):
        verdict = "LATE_CONTEXT_DOES_NOT_RECAPITULATE_BMR6_PATTERN"
    else:
        verdict = "LATE_CONTEXT_AMBIGUOUS"
    payload = {
        "verdict": verdict,
        "direct_bmr6_test_available": False,
        "primary_inoculated_libraries": int(len(primary)),
        "focus_genes_in_state_model": int(len(focus_ids)),
        "median_target_reads": median_target_reads,
        "minimum_group_size": minimum_group,
    }
    (args.output_dir / "context_summary.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
