"""Regression tests for complete-family inference and no-hit workflow behaviour."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))
import motif_inference as mi


def load(name):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class MotifTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.fg = {f"f{i}": "ACGTACGT" for i in range(5)}
        cls.bg = {f"b{i}": "ACGTACGT" for i in range(12)}
        cls.table, cls.summary = mi.analyse_family({**cls.fg, **cls.bg}, list(cls.fg), list(cls.bg), (2,))
        cls.h = load("test_motif_conservation")

    def test_canonical_reverse_complement(self):
        self.assertEqual(mi.canonical("ATGCTAA"), mi.canonical("TTAGCAT"))

    def test_full_family_sizes(self):
        self.assertEqual(len(mi.canonical_family((6,))), 2080)
        self.assertEqual(len(mi.canonical_family((7,))), 8192)
        self.assertEqual(len(mi.canonical_family((6, 7))), 10272)

    def test_repeat_count_is_presence(self):
        self.assertEqual(mi.promoter_presence("AAAAAAAAA", (2,)), {"AA"})

    def test_ambiguity_is_skipped(self):
        self.assertEqual(mi.promoter_presence("ANNNT", (2,)), set())

    def test_zero_foreground_in_family(self):
        row = self.table.set_index("motif").loc["AA"]
        self.assertEqual(int(row.focus_present), 0)
        self.assertEqual(float(row.fisher_p), 1)
        self.assertEqual(len(self.table), 10)

    def test_fisher_known_table(self):
        odds, p = mi.fisher_presence(6, 10, 38, 557)
        self.assertAlmostEqual(odds, 20.486842105263158)
        self.assertAlmostEqual(p, 2.5964516432639312e-5, places=15)

    def test_bh_matches_scipy(self):
        np.testing.assert_allclose(self.table.bh_full_family_q,
                                   false_discovery_control(self.table.fisher_p, method="bh"))

    def test_disjoint_groups_required(self):
        with self.assertRaises(ValueError):
            mi.analyse_family({"f": "ACGT"}, ["f"], ["f"], (2,))

    def test_duplicates_rejected(self):
        with self.assertRaises(ValueError):
            mi.analyse_family({**self.fg, **self.bg}, ["f0", "f0"], list(self.bg), (2,))

    def test_missing_sequences_rejected(self):
        with self.assertRaises(ValueError):
            mi.analyse_family(self.fg, list(self.fg), ["missing"], (2,))

    def test_auto_no_hit_is_valid(self):
        decision = mi.choose_follow_up(self.table, "auto", "AUTO")
        self.assertEqual(decision["status"], "NO_DISCOVERY_MOTIF")
        self.assertIsNone(decision["motif"])

    def test_explicit_followup_requires_exploratory_mode(self):
        with self.assertRaises(ValueError):
            mi.choose_follow_up(self.table, "auto", "AA")

    def test_exploratory_candidate_not_promoted(self):
        decision = mi.choose_follow_up(self.table, "exploratory", "AA")
        self.assertEqual(decision["status"], "EXPLORATORY_CANDIDATE_SELECTED")
        self.assertFalse(decision["passes_full_family_fdr"])

    def test_exploratory_auto_rejected(self):
        with self.assertRaises(ValueError):
            mi.choose_follow_up(self.table, "exploratory", "AUTO")

    def test_legacy_table_rejected(self):
        with self.assertRaises(ValueError):
            mi.choose_follow_up(self.table.drop(columns="hypothesis_family"), "auto", "AUTO")

    def test_truncated_family_rejected(self):
        with self.assertRaises(ValueError):
            mi.choose_follow_up(self.table.iloc[:-1], "auto", "AUTO")

    def test_forged_adjustment_rejected(self):
        t = self.table.copy()
        t["bh_full_family_q"] = 0
        with self.assertRaises(ValueError):
            mi.choose_follow_up(t, "auto", "AUTO")

    def test_multiple_auto_hits_are_not_cherry_picked(self):
        fg = {f"f{i}": "AAC" for i in range(10)}
        bg = {f"b{i}": "GGG" for i in range(40)}
        t, _ = mi.analyse_family({**fg, **bg}, list(fg), list(bg), (2,))
        self.assertEqual(mi.choose_follow_up(t, "auto", "AUTO")["status"],
                         "MULTIPLE_DISCOVERY_MOTIFS_NO_SELECTION")

    def test_palindromic_site_counted_once(self):
        self.assertEqual(self.h.scan_motif("GGATATCC", "ATAT")["motif_count"], 1)

    def test_plus_minus_promoter_coordinates(self):
        g = {"c": "ACGT" * 100}
        plus = pd.Series({"target_contig": "c", "target_start_0based": 280,
                          "target_end_0based_exclusive": 300, "strand": "+"})
        seq, left, right, _ = self.h.extract_oriented_promoter(plus, g, 200)
        self.assertEqual((left, right, seq), (80, 280, g["c"][80:280]))
        minus = plus.copy(); minus["strand"] = "-"
        seq, left, right, _ = self.h.extract_oriented_promoter(minus, g, 200)
        self.assertEqual((left, right, seq), (300, 400, self.h.revcomp(g["c"][300:400])))

    def test_unknown_between_record_distance(self):
        calls = pd.DataFrame([
            {"species": "s", "gene_id": "a", "is_focus": True, "analysis_valid": True,
             "target_contig": "c1", "target_start_0based": 1, "target_end_0based_exclusive": 10},
            {"species": "s", "gene_id": "b", "is_focus": True, "analysis_valid": True,
             "target_contig": "c2", "target_start_0based": 1, "target_end_0based_exclusive": 10}])
        result = self.h.sequence_record_layout(calls).iloc[0]
        self.assertTrue(pd.isna(result.minimum_within_sequence_intergenic_gap_bp))
        self.assertEqual(result.between_record_distances, "unknown")

    def test_same_record_gap(self):
        calls = pd.DataFrame([
            {"species": "s", "gene_id": "a", "is_focus": True, "analysis_valid": True,
             "target_contig": "c", "target_start_0based": 100, "target_end_0based_exclusive": 200},
            {"species": "s", "gene_id": "b", "is_focus": True, "analysis_valid": True,
             "target_contig": "c", "target_start_0based": 500, "target_end_0based_exclusive": 600}])
        self.assertEqual(self.h.sequence_record_layout(calls).iloc[0].minimum_within_sequence_intergenic_gap_bp, 300)

    def test_conservation_helper_signature_and_result(self):
        fg = [f"f{i}" for i in range(5)]
        bg = [f"b{i}" for i in range(10)]
        m = {**dict.fromkeys(fg, True), **dict.fromkeys(bg, False)}
        p, observed, n = self.h.empirical_species_p(fg, m, {g: bg for g in fg}, 9, np.random.default_rng(42))
        self.assertEqual((p, observed, n), (.1, 5, 5))

    def test_nohit_cli_does_not_require_genomes_and_clears_stale_table(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp); table = d / "motifs.tsv"
            self.table.to_csv(table, sep="\t", index=False)
            out = d / "out"; (out / "tables").mkdir(parents=True)
            (out / "tables/species_motif_enrichment.tsv").write_text("stale positive result")
            args = [sys.executable, str(SCRIPTS / "test_motif_conservation.py"),
                    "--motif-table", str(table), "--output-dir", str(out), "--work-dir", str(d / "work")]
            for key in ["focus-promoters", "matched-pools", "focus-promoter-fasta", "background-promoter-fasta", "source-proteins", "gene-annotation"]:
                args += ["--" + key, str(d / "not_present")]
            result = subprocess.run(args, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            status = json.loads((out / "tables/motif_conservation.json").read_text())
            self.assertEqual(status["status"], "NO_DISCOVERY_MOTIF")
            self.assertTrue(pd.read_csv(out / "tables/species_motif_enrichment.tsv", sep="\t").empty)


if __name__ == "__main__":
    unittest.main(verbosity=2)
