from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]


def test_sample_manifest():
    samples = pd.read_csv(ROOT / "resources" / "samples.tsv", sep="\t", dtype=str, keep_default_na=False)
    assert len(samples) == 129
    assert samples["sample_id"].is_unique
    assert samples["run_accession"].is_unique
    assert samples["fastq_url"].str.startswith("https://").all()
    assert samples["fastq_md5"].str.fullmatch(r"[0-9a-fA-F]{32}").all()
    inferred = samples[samples["timepoint_basis"].eq("inferred_from_published_notebook")]
    assert set(inferred["canonical_sample_id"]) == {
        "bmr6_pdb_wet_3_a", "bmr6_pdb_wet_3_b", "bmr6_pdb_wet_3_c"
    }


def test_analysis_sets():
    sets = pd.read_csv(ROOT / "resources" / "analysis_sets.tsv", sep="\t", dtype=str)
    observed = sets.groupby("set_name").size().to_dict()
    assert observed == {
        "background_3dai": 21,
        "context_13dai": 29,
        "fusarium_3dai": 24,
        "host_3dai": 48,
        "macrophomina_3dai": 24,
        "macrophomina_mapping_3dai": 45,
        "mapping_3dai": 45,
    }


def test_reference_manifest():
    references = pd.read_csv(ROOT / "resources" / "references.tsv", sep="\t", dtype=str)
    assert len(references) == 7
    assert references["label"].is_unique
    assert set(references.loc[references["include_in_competitive"].eq("1"), "species_group"]) == {
        "SB", "FT", "MP", "FV", "FPRO", "FFUJ"
    }
