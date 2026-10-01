"""Complete-family promoter-presence tests and explicit follow-up selection."""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from collections import Counter
from functools import lru_cache
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

FAMILY_NAME = "all_canonical_acgt_kmers"


def reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def canonical(sequence: str) -> str:
    sequence = sequence.upper()
    if not sequence or set(sequence) - set("ACGT"):
        raise ValueError("A motif must be a nonempty A/C/G/T sequence")
    return min(sequence, reverse_complement(sequence))


def canonical_family(lengths: Sequence[int]) -> list[str]:
    lengths = sorted(set(int(k) for k in lengths))
    if not lengths or any(k < 1 or k > 8 for k in lengths):
        raise ValueError("Motif lengths must be integers from 1 through 8")
    return sorted({canonical("".join(chars)) for k in lengths
                   for chars in itertools.product("ACGT", repeat=k)})


def read_sequences(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    name = None
    with path.open(encoding="utf-8-sig") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                fields = line[1:].split()
                if not fields or fields[0] in records:
                    raise ValueError(f"Empty or duplicate FASTA identifier in {path}")
                name = fields[0]
                records[name] = ""
            else:
                if name is None:
                    raise ValueError(f"Sequence before FASTA header in {path}")
                records[name] += line.upper()
    if not records or any(not value for value in records.values()):
        raise ValueError(f"Empty FASTA input or sequence: {path}")
    return records


def sequence_digest(records: Mapping[str, str]) -> str:
    text = "".join(f">{key}\n{records[key].upper()}\n" for key in sorted(records))
    return hashlib.sha256(text.encode("ascii")).hexdigest()


def promoter_presence(sequence: str, lengths: Sequence[int]) -> set[str]:
    sequence = sequence.upper()
    result: set[str] = set()
    for k in lengths:
        for offset in range(len(sequence) - k + 1):
            word = sequence[offset:offset + k]
            if not set(word) - set("ACGT"):
                result.add(canonical(word))
    return result


@lru_cache(maxsize=None)
def fisher_presence(a: int, n: int, c: int, m: int) -> tuple[float, float]:
    result = fisher_exact([[a, n - a], [c, m - c]], alternative="greater")
    return float(result.statistic), float(result.pvalue)


def analyse_family(
    promoters: Mapping[str, str], foreground: Sequence[str], background: Sequence[str],
    lengths: Sequence[int] = (6, 7), matched_sets: Sequence[Sequence[str]] | None = None,
    empirical_top: int = 50, alpha: float = 0.10,
) -> tuple[pd.DataFrame, dict]:
    """BH includes every canonical word; foreground presence never filters tests.

    Fixed-word matched-set tail probabilities are ancillary descriptive quantities,
    not scan-wide adjusted probabilities and not gates for BH eligibility.
    """
    if not 0 < alpha < 1 or empirical_top < 0:
        raise ValueError("Invalid FDR threshold or empirical-top count")
    foreground, background = list(foreground), list(background)
    if not foreground or not background:
        raise ValueError("Foreground and background must both be nonempty")
    if len(set(foreground)) != len(foreground) or len(set(background)) != len(background):
        raise ValueError("Promoter identifiers must be unique within each group")
    if set(foreground) & set(background):
        raise ValueError("Foreground and background must be disjoint")
    if (set(foreground) | set(background)) - set(promoters):
        raise ValueError("Missing promoter sequences")
    family = canonical_family(lengths)
    presence = {gene: promoter_presence(promoters[gene], lengths)
                for gene in foreground + background}
    fc = Counter(word for gene in foreground for word in presence[gene])
    bc = Counter(word for gene in background for word in presence[gene])
    rows = []
    for word in family:
        a, c = fc.get(word, 0), bc.get(word, 0)
        odds, p = fisher_presence(a, len(foreground), c, len(background))
        rows.append({"motif": word, "k": len(word), "focus_present": a,
                     "focus_total": len(foreground), "background_present": c,
                     "background_total": len(background), "odds_ratio": odds,
                     "fisher_p": p})
    table = pd.DataFrame(rows)
    table["fisher_fdr"] = multipletests(table["fisher_p"], method="fdr_bh")[1]
    table["bh_full_family_q"] = table["fisher_fdr"]
    table["hypothesis_family"] = FAMILY_NAME
    table["hypotheses_in_family"] = len(family)
    table["passes_full_family_fdr"] = table["fisher_fdr"] <= alpha
    table["fdr_threshold"] = alpha
    table = table.sort_values(["fisher_p", "motif"], kind="stable").reset_index(drop=True)
    table["matched_empirical_p"] = np.nan
    if matched_sets is not None:
        if not matched_sets:
            raise ValueError("Matched-set collection is empty")
        for chosen in matched_sets:
            if len(chosen) != len(foreground) or len(set(chosen)) != len(chosen):
                raise ValueError("Matched sets must be unique-gene sets of foreground size")
            if set(chosen) - set(background):
                raise ValueError("Matched sets must contain background promoters only")
        for i, row in table.head(empirical_top).iterrows():
            word, observed = str(row["motif"]), int(row["focus_present"])
            extreme = sum(sum(word in presence[gene] for gene in chosen) >= observed
                          for chosen in matched_sets)
            table.loc[i, "matched_empirical_p"] = (extreme + 1) / (len(matched_sets) + 1)
    table["empirical_scope"] = np.where(table["matched_empirical_p"].notna(),
                                         "fixed_motif_upper_tail_not_scan_adjusted", "not_evaluated")
    summary = {
        "hypothesis_family": FAMILY_NAME, "motif_lengths": sorted(set(lengths)),
        "hypotheses_in_family": len(family),
        "observed_in_pooled_sequences": len(set(fc) | set(bc)),
        "foreground_promoters": len(foreground), "background_promoters": len(background),
        "zero_foreground_hypotheses": int((table["focus_present"] == 0).sum()),
        "fisher_alternative": "greater", "adjustment": "Benjamini-Hochberg",
        "foreground_frequency_filter": False, "fdr_threshold": alpha,
        "discovery_motifs": int(table["passes_full_family_fdr"].sum()),
        "matched_iterations": len(matched_sets) if matched_sets is not None else 0,
        "matched_empirical_scope": "fixed_motif_not_scan_wide",
        "foreground_sequence_sha256": sequence_digest({g: promoters[g] for g in foreground}),
        "background_sequence_sha256": sequence_digest({g: promoters[g] for g in background}),
        "interpretation": "No passing motif is a valid result. BH adjustment does not establish regulatory function.",
    }
    return table, summary


def choose_follow_up(table: pd.DataFrame, mode: str, motif: str, alpha: float = 0.10) -> dict:
    """Do not silently turn a failed discovery into a positive follow-up."""
    required = {"motif", "k", "fisher_p", "bh_full_family_q", "hypothesis_family", "hypotheses_in_family"}
    if not required <= set(table):
        raise ValueError("Follow-up requires a complete-family motif table, not a foreground-filtered table")
    if table.empty or set(table["hypothesis_family"]) != {FAMILY_NAME}:
        raise ValueError("Invalid motif hypothesis-family metadata")
    if not 0 < alpha < 1:
        raise ValueError("Invalid selection FDR threshold")
    lengths = sorted(set(pd.to_numeric(table["k"]).astype(int)))
    family = canonical_family(lengths)
    if table["motif"].duplicated().any() or set(table["motif"]) != set(family):
        raise ValueError("Incomplete or duplicated canonical motif family")
    if not (pd.to_numeric(table["hypotheses_in_family"]) == len(family)).all():
        raise ValueError("Hypothesis count does not match the complete family")
    p = pd.to_numeric(table["fisher_p"], errors="raise").to_numpy(float)
    q = pd.to_numeric(table["bh_full_family_q"], errors="raise").to_numpy(float)
    if not np.isfinite(p).all() or not np.isfinite(q).all() or np.any((p < 0) | (p > 1) | (q < 0) | (q > 1)):
        raise ValueError("Invalid P or adjusted P values")
    if not np.allclose(q, multipletests(p, method="fdr_bh")[1], rtol=1e-10, atol=1e-12):
        raise ValueError("BH values do not reproduce from the supplied full-family P values")
    supported = table.loc[q <= alpha]
    result = {"selection_mode": mode, "fdr_threshold": alpha, "hypotheses_in_family": len(family),
              "discovery_motifs": len(supported), "motif": None,
              "full_family_q": None, "regulatory_function_validated": False,
              "comparative_inference": "selected_sequence_follow_up_not_independent_discovery_validation"}
    if mode == "auto":
        if motif.upper() != "AUTO":
            raise ValueError("An explicit motif requires --selection-mode exploratory")
        if len(supported) != 1:
            result["status"] = "NO_DISCOVERY_MOTIF" if not len(supported) else "MULTIPLE_DISCOVERY_MOTIFS_NO_SELECTION"
            return result
        chosen = supported.iloc[0]
        result["status"] = "DISCOVERY_SUPPORTED_CANDIDATE_SELECTED"
    elif mode == "exploratory":
        if motif.upper() == "AUTO":
            raise ValueError("Exploratory follow-up requires an explicitly named motif")
        word = canonical(motif)
        selected = table.loc[table["motif"].eq(word)]
        if len(selected) != 1:
            raise ValueError("Explicit motif is not in the tested family")
        chosen = selected.iloc[0]
        result["status"] = "EXPLORATORY_CANDIDATE_SELECTED"
    else:
        raise ValueError("Selection mode must be auto or exploratory")
    result.update(motif=str(chosen["motif"]), reverse_complement=reverse_complement(str(chosen["motif"])),
                  full_family_q=float(chosen["bh_full_family_q"]),
                  passes_full_family_fdr=bool(float(chosen["bh_full_family_q"]) <= alpha),
                  foreground_present=int(chosen["focus_present"]),
                  foreground_total=int(chosen["focus_total"]))
    return result


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--foreground", type=Path, required=True)
    ap.add_argument("--background", type=Path, required=True)
    ap.add_argument("--lengths", type=int, nargs="+", default=[6, 7])
    ap.add_argument("--fdr", type=float, default=0.10)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    fg, bg = read_sequences(args.foreground), read_sequences(args.background)
    table, summary = analyse_family({**fg, **bg}, list(fg), list(bg), args.lengths, alpha=args.fdr)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output_dir / "promoter_kmer_enrichment.tsv", sep="\t", index=False, na_rep="NA")
    (args.output_dir / "motif_search_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
