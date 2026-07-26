import json
import subprocess
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "workflow" / "scripts"


def test_manifest_validation(tmp_path):
    output = tmp_path / "validation.json"
    subprocess.run([
        sys.executable, str(SCRIPTS / "validate_inputs.py"),
        "--samples", str(ROOT / "resources" / "samples.tsv"),
        "--analysis-sets", str(ROOT / "resources" / "analysis_sets.tsv"),
        "--references", str(ROOT / "resources" / "references.tsv"),
        "--output", str(output),
    ], check=True)
    payload = json.loads(output.read_text())
    assert payload["samples"] == 129
    assert payload["analysis_sets"]["mapping_3dai"] == 45


def test_sample_selection(tmp_path):
    output = tmp_path / "selected.tsv"
    subprocess.run([
        sys.executable, str(SCRIPTS / "select_samples.py"),
        "--samples", str(ROOT / "resources" / "samples.tsv"),
        "--analysis-sets", str(ROOT / "resources" / "analysis_sets.tsv"),
        "--set", "fusarium_3dai:FUS",
        "--set", "background_3dai:PDB_PRIMARY",
        "--output", str(output),
    ], check=True)
    table = pd.read_csv(output, sep="\t")
    assert len(table) == 45
    assert table["analysis_set"].value_counts().to_dict() == {"FUS": 24, "PDB_PRIMARY": 21}
