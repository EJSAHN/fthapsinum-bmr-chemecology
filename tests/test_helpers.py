import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "workflow" / "scripts"


def load(name):
    path = SCRIPTS / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_reverse_complement_and_canonical_kmer():
    regulatory = load("analyze_regulatory_architecture")
    assert regulatory.revcomp("ATGCTAA") == "TTAGCAT"
    assert regulatory.canonical_kmer("ATGCTAA") == "ATGCTAA"
    assert regulatory.canonical_kmer("TTAGCAT") == "ATGCTAA"


def test_sample_name_normalization():
    coupling = load("host_module_coupling")
    assert coupling.canonical_sample_id("bmr6_pdb_wet_a") == "bmr6_pdb_wet_3_a"
    assert coupling.canonical_sample_id("WT-FUS-Dry-3-A") == "wt_fus_dry_3_a"
    assert coupling.split_sample("bmr12_fus_wet_13_c") == {
        "genotype": "bmr12", "pathogen": "fus", "water": "wet", "dai": 13, "replicate": "c"
    }


def test_featurecounts_sample_names_and_raw_parser(tmp_path):
    selector = load("select_count_model")
    context = load("prepare_13dai_context")
    assert selector.sample_name("sample.FT.mapq20.primary.bam") == "sample"
    assert selector.sample_name("sample.MP.mapq20.primary.bam") == "sample"
    assert context.sample_name("sample.MP.mapq20.primary.bam") == "sample"

    counts_path = tmp_path / "counts.txt"
    counts_path.write_text(
        "# Program: featureCounts\n"
        "Geneid\tChr\tStart\tEnd\tStrand\tLength\t/a/sample.FT.mapq20.primary.bam\t/b/other.MP.mapq20.primary.bam\n"
        "gene1\tcontig1\t1\t10\t+\t10\t4\t7\n",
        encoding="utf-8",
    )
    parsed = context.read_featurecounts_counts(counts_path)
    assert parsed.columns.tolist() == ["sample", "other"]
    assert parsed.loc["gene1"].tolist() == [4, 7]
