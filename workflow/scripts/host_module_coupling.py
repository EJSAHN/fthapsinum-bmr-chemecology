#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm


METADATA_COLUMNS = {
    "sample_id", "canonical_sample_id", "genotype", "pathogen", "water",
    "dai", "replicate", "sequencing_run", "group", "cell",
}


def normalize_module_name(value: object) -> str:
    text = str(value).strip().lower()
    text = re.sub(r"^me[._ -]*", "", text)
    return re.sub(r"[^a-z0-9]+", "", text)


def canonical_sample_id(value: object) -> str | None:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return None
    text = str(value).strip().strip('"').strip("'").lower()
    replacements = {
        "f. thapsinum": "fus", "f_thapsinum": "fus",
        "fusarium_thapsinum": "fus", "fusarium": "fus",
        "macrophomina_phaseolina": "macro", "macrophomina": "macro",
        "well-watered": "wet", "well_watered": "wet",
        "water-limited": "dry", "water_limited": "dry",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    text = re.sub(r"[ .:/\\-]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    full = re.fullmatch(r"(bmr12|bmr6|wt)_(fus|macro|pdb)_(dry|wet)_(0|3|13)_([a-h])", text)
    if full:
        return "_".join(full.groups())
    abbreviated = re.fullmatch(r"(bmr12|bmr6|wt)_(fus|macro|pdb)_(dry|wet)_([a-h])", text)
    if abbreviated:
        genotype, pathogen, water, replicate = abbreviated.groups()
        return f"{genotype}_{pathogen}_{water}_3_{replicate}"
    return None


def split_sample(sample_id: str) -> dict[str, str | int]:
    match = re.fullmatch(r"(bmr12|bmr6|wt)_(fus|macro|pdb)_(dry|wet)_(0|3|13)_([a-h])", sample_id)
    if not match:
        raise ValueError(f"Unrecognized sample ID: {sample_id}")
    genotype, pathogen, water, dai, replicate = match.groups()
    return {
        "genotype": genotype,
        "pathogen": pathogen,
        "water": water,
        "dai": int(dai),
        "replicate": replicate,
    }


def read_delimited(path: Path) -> pd.DataFrame:
    errors: list[str] = []
    for options in ({"sep": ","}, {"sep": None, "engine": "python"}, {"sep": "\t"}):
        try:
            table = pd.read_csv(path, **options)
            if table.shape[0] > 1 and table.shape[1] > 1:
                return table
        except Exception as error:
            errors.append(str(error))
    raise RuntimeError(f"Could not parse {path}: {' | '.join(errors)}")


def load_eigengenes(path: Path) -> tuple[pd.DataFrame, dict[str, object]]:
    raw = read_delimited(path)
    diagnostics: dict[str, object] = {"rows": int(raw.shape[0]), "columns": int(raw.shape[1])}
    matches = {str(column): sum(canonical_sample_id(value) is not None for value in raw[column]) for column in raw.columns}
    sample_column = max(matches, key=matches.get)
    row_matches = matches[sample_column]

    if row_matches >= max(12, int(0.2 * len(raw))):
        work = raw.copy()
        work["sample_id"] = work[sample_column].map(canonical_sample_id)
        diagnostics["orientation"] = "samples_in_rows"
        diagnostics["sample_column"] = sample_column
    else:
        sample_columns = [column for column in raw.columns if canonical_sample_id(column) is not None]
        if len(sample_columns) < 12:
            raise RuntimeError("Sample IDs were not detected in the eigengene table")
        label_columns = [column for column in raw.columns if column not in sample_columns]
        if not label_columns:
            raise RuntimeError("No module-label column was detected in the transposed eigengene table")
        work = raw.set_index(label_columns[0])[sample_columns].T
        work.index = [canonical_sample_id(value) for value in work.index]
        work.index.name = "sample_id"
        work = work.reset_index()
        diagnostics["orientation"] = "samples_in_columns"
        diagnostics["module_label_column"] = str(label_columns[0])

    work = work[work["sample_id"].notna()].copy()
    modules: dict[str, str] = {}
    for column in work.columns:
        if column == "sample_id" or str(column) == str(sample_column):
            continue
        numeric = pd.to_numeric(work[column], errors="coerce")
        if numeric.notna().mean() < 0.7:
            continue
        name = normalize_module_name(column)
        if not name or name in {"sample", "sampleid", "name", "rowname", "rownames", "x"}:
            continue
        if name in modules.values():
            raise RuntimeError(f"Duplicate normalized module name: {name}")
        modules[str(column)] = name

    if len(modules) < 5:
        raise RuntimeError(f"Only {len(modules)} numeric modules were detected")

    output = pd.DataFrame({"sample_id": work["sample_id"].astype(str)})
    for original, normalized in modules.items():
        output[normalized] = pd.to_numeric(work[original], errors="coerce")
    output["_complete"] = output.drop(columns="sample_id").notna().sum(axis=1)
    duplicates = output.loc[output.duplicated("sample_id", keep=False), "sample_id"].unique().tolist()
    output = (
        output.sort_values(["sample_id", "_complete"], ascending=[True, False])
        .drop_duplicates("sample_id")
        .drop(columns="_complete")
        .sort_values("sample_id")
        .reset_index(drop=True)
    )
    diagnostics.update({
        "canonical_samples": int(len(output)),
        "modules": int(len(modules)),
        "duplicate_samples_resolved": duplicates,
        "module_names": list(output.columns[1:]),
    })
    return output, diagnostics


def bh_adjust(values: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(values), dtype=float)
    output = np.full_like(p, np.nan)
    valid = np.isfinite(p)
    if not valid.any():
        return output
    observed = p[valid]
    order = np.argsort(observed)
    ranked = observed[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    inverse = np.empty_like(order)
    inverse[order] = np.arange(len(order))
    output[np.where(valid)[0]] = np.clip(adjusted[inverse], 0, 1)
    return output


def spearman_safe(first: pd.Series, second: pd.Series) -> tuple[float, float, int]:
    data = pd.DataFrame({"first": pd.to_numeric(first, errors="coerce"), "second": pd.to_numeric(second, errors="coerce")}).dropna()
    if len(data) < 5 or data["first"].nunique() < 2 or data["second"].nunique() < 2:
        return np.nan, np.nan, len(data)
    rho, pvalue = stats.spearmanr(data["first"], data["second"])
    return float(rho), float(pvalue), int(len(data))


def covariate_matrix(frame: pd.DataFrame) -> np.ndarray:
    groups = pd.get_dummies(frame["group"].astype(str), drop_first=True, dtype=float)
    matrix = pd.concat([
        pd.Series(1.0, index=frame.index, name="intercept"),
        groups,
        np.log1p(pd.to_numeric(frame["target_unique_rpm"], errors="coerce")).rename("log_target_rpm"),
    ], axis=1)
    return matrix.to_numpy(dtype=float)


def partial_association(frame: pd.DataFrame, module: str) -> dict[str, float]:
    data = frame[["fungal_score", module, "group", "target_unique_rpm"]].dropna().copy()
    if len(data) < 12:
        return {"partial_beta": np.nan, "partial_se_hc3": np.nan, "partial_p_hc3": np.nan, "n": len(data)}
    covariates = covariate_matrix(data)
    host = pd.to_numeric(data[module], errors="coerce").to_numpy(float)
    fungal = pd.to_numeric(data["fungal_score"], errors="coerce").to_numpy(float)
    fit = sm.OLS(fungal, np.column_stack([covariates, host])).fit(cov_type="HC3")
    return {
        "partial_beta": float(fit.params[-1]),
        "partial_se_hc3": float(fit.bse[-1]),
        "partial_p_hc3": float(fit.pvalues[-1]),
        "n": int(len(data)),
    }


def stratified_permutation_indices(groups: pd.Series, iterations: int, rng: np.random.Generator) -> np.ndarray:
    values = groups.astype(str).to_numpy()
    output = np.tile(np.arange(len(values), dtype=int), (iterations, 1))
    for group in pd.unique(values):
        positions = np.where(values == group)[0]
        order = np.argsort(rng.random((iterations, len(positions))), axis=1)
        output[:, positions] = positions[order]
    return output


def permutation_partial_p(frame: pd.DataFrame, module: str, permutations: np.ndarray) -> tuple[float, float]:
    data = frame[["fungal_score", module, "group", "target_unique_rpm"]].dropna().copy()
    if len(data) != len(frame):
        return np.nan, np.nan
    covariates = covariate_matrix(data)
    host = pd.to_numeric(data[module], errors="coerce").to_numpy(float)
    fungal = pd.to_numeric(data["fungal_score"], errors="coerce").to_numpy(float)
    host_residual = sm.OLS(host, covariates).fit().resid
    fungal_residual = sm.OLS(fungal, covariates).fit().resid
    denominator = float(np.dot(host_residual, host_residual))
    if denominator <= 0:
        return np.nan, np.nan
    observed = float(np.dot(host_residual, fungal_residual) / denominator)
    null = fungal_residual[permutations] @ host_residual / denominator
    pvalue = (1 + np.sum(np.abs(null) >= abs(observed))) / (len(null) + 1)
    return observed, float(pvalue)


def leave_one_out_stability(frame: pd.DataFrame, module: str, beta: float) -> float:
    if not np.isfinite(beta) or beta == 0:
        return np.nan
    stable: list[bool] = []
    for index in frame.index:
        estimate = partial_association(frame.drop(index), module)["partial_beta"]
        if np.isfinite(estimate) and estimate != 0:
            stable.append(np.sign(estimate) == np.sign(beta))
    return float(np.mean(stable)) if stable else np.nan


def zmean_score(logcpm_path: Path, focus_genes: list[str]) -> pd.DataFrame:
    table = pd.read_csv(logcpm_path, sep="\t", compression="infer")
    available = [gene for gene in focus_genes if gene in set(table["gene_id"].astype(str))]
    if len(available) < 5:
        raise RuntimeError(f"Only {len(available)} focus genes were found in {logcpm_path}")
    matrix = table.set_index("gene_id").loc[available].T.apply(pd.to_numeric, errors="coerce")
    zscores = (matrix - matrix.mean(axis=0)) / matrix.std(axis=0, ddof=0).replace(0, np.nan)
    return zscores.mean(axis=1).rename("fungal_score").reset_index().rename(columns={"index": "sample_id"})


def fit_host_factorial(host: pd.DataFrame, modules: list[str], iterations: int, rng: np.random.Generator) -> pd.DataFrame:
    cells = sorted(host["cell"].unique())
    design = pd.get_dummies(host["cell"], dtype=float).reindex(columns=cells, fill_value=0.0)
    matrix = design.to_numpy(float)
    contrasts = {
        "bmr6_fus_water_interaction": {"bmr6_fus_wet": 1, "bmr6_fus_dry": -1, "wt_fus_wet": -1, "wt_fus_dry": 1},
        "bmr6_pdb_water_interaction": {"bmr6_pdb_wet": 1, "bmr6_pdb_dry": -1, "wt_pdb_wet": -1, "wt_pdb_dry": 1},
        "bmr6_three_way": {
            "bmr6_fus_wet": 1, "bmr6_fus_dry": -1, "wt_fus_wet": -1, "wt_fus_dry": 1,
            "bmr6_pdb_wet": -1, "bmr6_pdb_dry": 1, "wt_pdb_wet": 1, "wt_pdb_dry": -1,
        },
        "bmr12_fus_water_interaction": {"bmr12_fus_wet": 1, "bmr12_fus_dry": -1, "wt_fus_wet": -1, "wt_fus_dry": 1},
        "bmr12_pdb_water_interaction": {"bmr12_pdb_wet": 1, "bmr12_pdb_dry": -1, "wt_pdb_wet": -1, "wt_pdb_dry": 1},
        "bmr12_three_way": {
            "bmr12_fus_wet": 1, "bmr12_fus_dry": -1, "wt_fus_wet": -1, "wt_fus_dry": 1,
            "bmr12_pdb_wet": -1, "bmr12_pdb_dry": 1, "wt_pdb_wet": 1, "wt_pdb_dry": -1,
        },
    }

    difference_weights: dict[tuple[str, str], np.ndarray] = {}
    for genotype in ("wt", "bmr6", "bmr12"):
        for water in ("dry", "wet"):
            positions = np.where((host["genotype"].to_numpy() == genotype) & (host["water"].to_numpy() == water))[0]
            if len(positions) != 8:
                raise RuntimeError(f"Expected eight samples in {genotype}/{water}; found {len(positions)}")
            selected = np.argsort(rng.random((iterations, len(positions))), axis=1)[:, :4]
            assignment = np.zeros((iterations, len(positions)), dtype=bool)
            assignment[np.arange(iterations)[:, None], selected] = True
            weights = np.zeros((iterations, len(host)), dtype=np.float32)
            weights[:, positions] = np.where(assignment, 0.25, -0.25)
            difference_weights[(genotype, water)] = weights
    permutation_weights = {
        "bmr6_three_way": difference_weights[("bmr6", "wet")] - difference_weights[("bmr6", "dry")] - difference_weights[("wt", "wet")] + difference_weights[("wt", "dry")],
        "bmr12_three_way": difference_weights[("bmr12", "wet")] - difference_weights[("bmr12", "dry")] - difference_weights[("wt", "wet")] + difference_weights[("wt", "dry")],
    }

    rows: list[dict[str, object]] = []
    for module in modules:
        response = pd.to_numeric(host[module], errors="coerce").to_numpy(float)
        if not np.isfinite(response).all():
            continue
        fit = sm.OLS(response, matrix).fit(cov_type="HC3")
        for contrast_name, weights in contrasts.items():
            vector = np.array([weights.get(cell, 0.0) for cell in cells], dtype=float)
            test = fit.t_test(vector)
            effect = float(np.asarray(test.effect).ravel()[0])
            pvalue = float(np.asarray(test.pvalue).ravel()[0])
            permutation_p = np.nan
            if contrast_name in permutation_weights:
                null = permutation_weights[contrast_name] @ response
                permutation_p = (1 + np.sum(np.abs(null) >= abs(effect))) / (len(null) + 1)
            rows.append({
                "module": module,
                "contrast": contrast_name,
                "effect": effect,
                "se_hc3": float(np.asarray(test.sd).ravel()[0]),
                "p_hc3": pvalue,
                "permutation_p": float(permutation_p) if np.isfinite(permutation_p) else np.nan,
            })
    output = pd.DataFrame(rows)
    if not output.empty:
        output["p_hc3_fdr"] = output.groupby("contrast")["p_hc3"].transform(bh_adjust)
        output["permutation_fdr"] = output.groupby("contrast")["permutation_p"].transform(bh_adjust)
    return output


def module_trait_correlations(eigengenes: pd.DataFrame, traits: pd.DataFrame, modules: list[str]) -> pd.DataFrame:
    metadata = eigengenes["sample_id"].map(split_sample).apply(pd.Series)
    host = pd.concat([eigengenes, metadata], axis=1)
    host = host[(host["dai"] == 3) & (host["pathogen"] == "fus")].copy()
    traits = traits.copy()
    if "canonical_sample_id" in traits.columns:
        traits["sample_id"] = traits["canonical_sample_id"].replace("", np.nan).fillna(traits["sample_id"])
    numeric_traits = []
    for column in traits.columns:
        if column in METADATA_COLUMNS or column.endswith("accession") or column.endswith("url") or column.endswith("md5"):
            continue
        converted = pd.to_numeric(traits[column], errors="coerce")
        if converted.notna().sum() >= 8:
            traits[column] = converted
            numeric_traits.append(column)
    merged = host.merge(traits[["sample_id"] + numeric_traits], on="sample_id", how="left")
    rows: list[dict[str, object]] = []
    for module in modules:
        for trait in numeric_traits:
            rho, pvalue, n = spearman_safe(merged[module], merged[trait])
            if n >= 8 and np.isfinite(pvalue):
                rows.append({"module": module, "trait": trait, "rho": rho, "p": pvalue, "n": n})
    output = pd.DataFrame(rows)
    if not output.empty:
        output["fdr"] = bh_adjust(output["p"])
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description="Test associations between published host eigengenes and a fixed fungal module.")
    parser.add_argument("--eigengenes", type=Path, required=True)
    parser.add_argument("--module-scores", type=Path, required=True)
    parser.add_argument("--focus-genes", type=Path, required=True)
    parser.add_argument("--logcpm-all", type=Path, required=True)
    parser.add_argument("--samples-primary", type=Path, required=True)
    parser.add_argument("--samples-all", type=Path, required=True)
    parser.add_argument("--sample-traits", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--focus-module", default="aromatic_phenolic_metabolism")
    parser.add_argument("--permutations", type=int, default=50000)
    parser.add_argument("--seed", type=int, default=20260716)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    eigengenes, diagnostics = load_eigengenes(args.eigengenes)
    modules = list(eigengenes.columns[1:])
    eigengenes.to_csv(args.output_dir / "host_eigengenes.tsv", sep="\t", index=False, lineterminator="\n")
    (args.output_dir / "eigengene_import.json").write_text(json.dumps(diagnostics, indent=2), encoding="utf-8")

    host_metadata = eigengenes["sample_id"].map(split_sample).apply(pd.Series)
    host = pd.concat([eigengenes, host_metadata], axis=1)
    host["cell"] = host["genotype"] + "_" + host["pathogen"] + "_" + host["water"]
    host48 = host[(host["dai"] == 3) & host["pathogen"].isin(["fus", "pdb"]) & host["replicate"].isin(list("abcd"))].copy()
    counts = host48.groupby(["genotype", "pathogen", "water"], observed=True).size()
    if len(host48) != 48 or len(counts) != 12 or not (counts == 4).all():
        raise RuntimeError("The eigengene table did not yield the expected balanced 48-sample design")
    host48 = host48.sort_values(["genotype", "water", "pathogen", "replicate"]).reset_index(drop=True)

    scores = pd.read_csv(args.module_scores, sep="\t")
    primary_scores = scores[(scores["module"] == args.focus_module) & (scores["score_type"] == "zmean")].copy()
    primary_scores = primary_scores.rename(columns={"score": "fungal_score"})
    samples_primary = pd.read_csv(args.samples_primary, sep="\t")
    if "target_unique_rpm" not in samples_primary.columns and "ft_unique_rpm_input" in samples_primary.columns:
        samples_primary["target_unique_rpm"] = samples_primary["ft_unique_rpm_input"]
    primary = primary_scores.copy()
    if "group" not in primary.columns:
        primary = primary.merge(
            samples_primary[["sample_id", "genotype", "water", "group"]],
            on="sample_id", how="left", validate="one_to_one",
        )
    primary = primary.merge(
        samples_primary[["sample_id", "target_unique_rpm"]],
        on="sample_id", how="inner", validate="one_to_one",
    )
    primary = primary.merge(eigengenes, on="sample_id", how="inner", validate="one_to_one")
    if len(primary) != 22:
        raise RuntimeError(f"Expected 22 primary samples; found {len(primary)}")

    focus_genes = pd.read_csv(args.focus_genes, sep="\t")["gene_id"].astype(str).tolist()
    all_scores = zmean_score(args.logcpm_all, focus_genes)
    samples_all = pd.read_csv(args.samples_all, sep="\t")
    if "target_unique_rpm" not in samples_all.columns and "ft_unique_rpm_input" in samples_all.columns:
        samples_all["target_unique_rpm"] = samples_all["ft_unique_rpm_input"]
    all_samples = all_scores.merge(samples_all[["sample_id", "genotype", "water", "group", "target_unique_rpm"]], on="sample_id", how="inner")
    all_samples = all_samples.merge(eigengenes, on="sample_id", how="inner")
    if len(all_samples) != 24:
        raise RuntimeError(f"Expected 24 sensitivity samples; found {len(all_samples)}")

    rng = np.random.default_rng(args.seed)
    permutation_indices = stratified_permutation_indices(primary["group"], args.permutations, rng)
    rows: list[dict[str, object]] = []
    for module in modules:
        rho, raw_p, n = spearman_safe(primary[module], primary["fungal_score"])
        adjusted = partial_association(primary, module)
        permutation_beta, permutation_p = permutation_partial_p(primary, module, permutation_indices)
        all_adjusted = partial_association(all_samples, module)
        rows.append({
            "module": module,
            "n": n,
            "raw_spearman_rho": rho,
            "raw_spearman_p": raw_p,
            **adjusted,
            "permutation_beta": permutation_beta,
            "permutation_p": permutation_p,
            "leave_one_out_sign_stability": leave_one_out_stability(primary, module, adjusted["partial_beta"]),
            "all_samples_partial_beta": all_adjusted["partial_beta"],
            "all_samples_partial_p_hc3": all_adjusted["partial_p_hc3"],
            "all_samples_direction_consistent": bool(
                np.isfinite(adjusted["partial_beta"]) and np.isfinite(all_adjusted["partial_beta"]) and
                np.sign(adjusted["partial_beta"]) == np.sign(all_adjusted["partial_beta"])
            ),
        })
    associations = pd.DataFrame(rows)
    associations["raw_spearman_fdr"] = bh_adjust(associations["raw_spearman_p"])
    associations["partial_fdr_hc3"] = bh_adjust(associations["partial_p_hc3"])
    associations["permutation_fdr"] = bh_adjust(associations["permutation_p"])
    associations["all_samples_partial_fdr_hc3"] = bh_adjust(associations["all_samples_partial_p_hc3"])
    associations = associations.sort_values(["permutation_fdr", "permutation_p", "partial_p_hc3"], na_position="last")
    associations.to_csv(args.output_dir / "host_fungal_associations.tsv", sep="\t", index=False, lineterminator="\n")

    factorial = fit_host_factorial(host48, modules, args.permutations, np.random.default_rng(args.seed + 1))
    factorial.to_csv(args.output_dir / "host_factorial_contrasts.tsv", sep="\t", index=False, lineterminator="\n")
    traits = pd.read_csv(args.sample_traits, sep="\t", dtype=str, keep_default_na=False)
    trait_associations = module_trait_correlations(eigengenes, traits, modules)
    trait_associations.to_csv(args.output_dir / "host_module_trait_associations.tsv", sep="\t", index=False, lineterminator="\n")

    summary = {
        "host_modules": len(modules),
        "primary_samples": len(primary),
        "sensitivity_samples": len(all_samples),
        "balanced_host_samples": len(host48),
        "permutations": args.permutations,
        "modules_with_permutation_fdr_below_0_05": int((associations["permutation_fdr"] < 0.05).sum()),
        "modules_with_three_way_permutation_fdr_below_0_05": int((factorial.loc[factorial["contrast"].str.endswith("three_way"), "permutation_fdr"] < 0.05).sum()),
    }
    (args.output_dir / "host_module_coupling.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
