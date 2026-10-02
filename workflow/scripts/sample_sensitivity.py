#!/usr/bin/env python3
"""Prepare fixed-count sample subsets and summarise their refitted models."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from score_sensitivity import analyse, standardised_score


def read(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def save(table: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, sep="\t", index=False, na_rep="NA", lineterminator="\n")


def unique(table: pd.DataFrame, column: str) -> None:
    if column not in table or table[column].isna().any() or table[column].duplicated().any():
        raise ValueError(f"Missing or duplicate {column}")


def prepare(args: argparse.Namespace) -> None:
    cfg = json.loads(args.settings.read_text(encoding="utf-8"))
    counts, samples, sets, mapping = map(read, (args.counts, args.samples, args.gene_sets, args.mapping))
    unique(counts, "gene_id")
    unique(samples, "sample_id")
    unique(mapping, "sample_id")
    genes = counts.gene_id.astype(str).tolist()
    ids = samples.sample_id.tolist()
    if len(ids) != cfg["expected_primary_n"] or len(genes) != cfg["expected_genes"]:
        raise ValueError("Primary sample count or fixed gene universe differs from the settings")
    excluded = cfg["exclude_samples"]
    if len(set(excluded)) != len(excluded) or not set(excluded).issubset(ids):
        raise ValueError("Invalid explicitly configured exclusions")
    if len(ids) - len(excluded) != cfg["expected_reduced_n"]:
        raise ValueError("Reduced sample count differs from the settings")
    if not set(ids).issubset(counts.columns) or not set(ids).issubset(mapping.sample_id):
        raise ValueError("Counts or mapping data lack a required sample")
    values = counts.loc[:, ids].to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any() or (values != np.floor(values)).any():
        raise ValueError("Counts must be finite non-negative integers")
    group = samples.genotype.astype(str) + "_" + samples.water.astype(str)
    if "group" in samples and not samples.group.eq(group).all():
        raise ValueError("Group labels disagree with genotype and water")
    samples["group"] = group
    ordered = mapping.set_index("sample_id").loc[ids]
    for column in ("all_primary", "target_unique_rpm_input", "target_unique_rpm_post_filter"):
        value = pd.to_numeric(ordered[column], errors="raise").to_numpy(float)
        if not np.isfinite(value).all() or (value < 0).any():
            raise ValueError(f"Invalid {column}")
        if column in samples and not np.allclose(samples[column].to_numpy(float), value, rtol=1e-12, atol=1e-12):
            raise ValueError(f"Existing sample values disagree with mapping: {column}")
        samples[column] = value
    if not np.array_equal(samples.all_primary.to_numpy(float), ordered.post_fastp_reads.to_numpy(float)):
        raise ValueError("Tissue offsets must equal post-filter primary-read counts")
    if not (samples.all_primary > 0).all():
        raise ValueError("Tissue offsets must be positive")
    focus = sets.loc[sets.module.eq(cfg["focus_module"]) & sets.gene_id.isin(genes), "gene_id"].drop_duplicates().tolist()
    if len(focus) != cfg["expected_focus_genes"]:
        raise ValueError("Fixed functional set has unexpected membership")
    cfg["focus_genes"] = focus
    reduced = samples.loc[~samples.sample_id.isin(excluded)]
    for frame in (samples, reduced):
        if set(frame.group) != {"wt_dry", "bmr6_dry", "bmr12_dry", "wt_wet", "bmr6_wet", "bmr12_wet"} or frame.groupby("group").size().min() < 2:
            raise ValueError("Missing or under-replicated genotype-water cell")
    args.output.mkdir(parents=True, exist_ok=True)
    save(counts, args.output / "counts.tsv")
    save(samples, args.output / "samples.tsv")
    save(pd.DataFrame({"gene_id": genes}), args.output / "genes.tsv")
    save(sets, args.output / "sets.tsv")
    qc = ordered.reset_index().copy()
    qc["excluded_in_sensitivity"] = qc.sample_id.isin(excluded)
    qc["assignment_fraction_below_half"] = qc.target_specificity.lt(0.5)
    save(qc, args.output / "species_assignment_qc.tsv")
    (args.output / "settings.json").write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")


def summarise(args: argparse.Namespace) -> None:
    cfg = json.loads((args.prepared / "settings.json").read_text(encoding="utf-8"))
    if not (args.refit / "REFIT_COMPLETE.txt").is_file():
        raise FileNotFoundError("The four-model refit has not completed")
    primary = read(args.prepared / "samples.tsv")
    reduced = primary.loc[~primary.sample_id.isin(cfg["exclude_samples"])].copy()
    tables, scores, camera, fixed = [], [], [], []
    for label, samples in (("baseline22_state", primary), ("exclude21_state", reduced)):
        lc = read(args.refit / label / "logCPM.tsv.gz").set_index("gene_id")
        if list(lc.columns) != samples.sample_id.tolist():
            raise ValueError("Refit sample order differs from the prepared subset")
        score = standardised_score(lc, cfg["focus_genes"])
        tables.append(analyse(score, samples, label))
        scores.append(pd.DataFrame({"analysis": label, "sample_id": score.index, "score": score.values}))
        if label == "baseline22_state":
            fixed.extend((analyse(score, primary, "baseline_fixed_scores"),
                          analyse(score, reduced, "excluded_fixed_scores_diagnostic")))
    for label in ("baseline22_state", "baseline22_tissue", "exclude21_state", "exclude21_tissue"):
        info = json.loads((args.refit / label / "model.json").read_text())
        expected = cfg["expected_primary_n"] if label.startswith("baseline") else cfg["expected_reduced_n"]
        if info["n_samples"] != expected or info["n_genes"] != cfg["expected_genes"]:
            raise ValueError("Refit dimensions disagree with the fixed subsets")
        table = read(args.refit / label / "camera.tsv")
        table = table.loc[table.module.eq(cfg["focus_module"])].copy()
        table.insert(0, "analysis", label)
        camera.append(table)
    args.output.mkdir(parents=True, exist_ok=True)
    save(pd.concat(tables, ignore_index=True), args.output / "score_sensitivity.tsv")
    save(pd.concat(scores, ignore_index=True), args.output / "refitted_gene_set_scores.tsv")
    save(pd.concat(camera, ignore_index=True), args.output / "fixed_set_camera_sensitivity.tsv")
    save(pd.concat(fixed, ignore_index=True), args.output / "fixed_score_diagnostic.tsv")
    save(read(args.prepared / "species_assignment_qc.tsv"), args.output / "species_assignment_qc.tsv")
    status = {"status": "COMPLETED", "primary_n": cfg["expected_primary_n"],
              "reduced_n": cfg["expected_reduced_n"], "excluded_samples": cfg["exclude_samples"],
              "fixed_gene_universe": cfg["expected_genes"], "focus_genes": cfg["focus_genes"],
              "significance_required_for_completion": False}
    (args.output / "sample_sensitivity.json").write_text(json.dumps(status, indent=2) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    prep = sub.add_parser("prepare")
    for name in ("counts", "samples", "gene-sets", "mapping", "settings", "output"):
        prep.add_argument("--" + name, type=Path, required=True)
    summary = sub.add_parser("summarise")
    for name in ("prepared", "refit", "output"):
        summary.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    (prepare if args.command == "prepare" else summarise)(args)


if __name__ == "__main__":
    main()
