from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_ROOT = ROOT / "workflow" / "scripts"


def test_scripts_are_portable_and_table_only():
    forbidden = [
        "C:" + "\\",
        "/" + "mnt/",
        "/" + "home/",
        "matplot" + "lib",
        "py" + "plot",
        "save" + "fig",
        "gg" + "plot",
        "plot" + "ly",
        "this " + "revision",
        "Chat" + "GPT",
        "Open" + "AI",
    ]
    for path in SCRIPT_ROOT.iterdir():
        if path.suffix not in {".py", ".R", ".sh"}:
            continue
        text = path.read_text(encoding="utf-8")
        for token in forbidden:
            assert token not in text, f"{token!r} found in {path.relative_to(ROOT)}"


def test_no_binary_files_are_tracked():
    for path in ROOT.rglob("*"):
        if not path.is_file() or ".git" in path.parts:
            continue
        assert path.stat().st_size < 5_000_000, f"Unexpected large tracked file: {path.relative_to(ROOT)}"


def test_text_files_use_lf_line_endings():
    suffixes = {".py", ".R", ".sh", ".md", ".yaml", ".yml", ".tsv", ".json"}
    candidates = [ROOT / "Snakefile", ROOT / ".gitattributes", ROOT / ".gitignore"]
    candidates.extend(path for path in ROOT.rglob("*") if path.is_file() and path.suffix in suffixes)
    for path in candidates:
        if not path.exists() or any(part in {".git", ".pytest_cache", "__pycache__"} for part in path.parts):
            continue
        assert b"\r\n" not in path.read_bytes(), f"CRLF line endings found in {path.relative_to(ROOT)}"
