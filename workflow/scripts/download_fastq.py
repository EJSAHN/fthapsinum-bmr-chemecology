#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import shutil
import subprocess
from pathlib import Path


def md5sum(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def gzip_ok(path: Path) -> bool:
    try:
        with gzip.open(path, "rb") as handle:
            while handle.read(8 * 1024 * 1024):
                pass
        return True
    except (OSError, EOFError):
        return False


def count_reads(path: Path) -> int:
    command = ["seqkit", "stats", "-T", str(path)]
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    lines = [line for line in completed.stdout.splitlines() if line.strip()]
    if len(lines) < 2:
        raise RuntimeError(f"seqkit did not return a data row for {path}")
    header = lines[0].split("\t")
    values = lines[1].split("\t")
    record = dict(zip(header, values))
    return int(record["num_seqs"].replace(",", ""))


def main() -> None:
    parser = argparse.ArgumentParser(description="Download and validate one ENA FASTQ file.")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--url", required=True)
    parser.add_argument("--md5", required=True)
    parser.add_argument("--expected-reads", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metrics", type=Path, required=True)
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.metrics.parent.mkdir(parents=True, exist_ok=True)
    part = args.output.with_suffix(args.output.suffix + ".part")

    valid_existing = (
        args.output.exists()
        and args.output.stat().st_size > 0
        and gzip_ok(args.output)
        and md5sum(args.output).lower() == args.md5.lower()
    )
    if not valid_existing:
        if args.output.exists():
            args.output.unlink()
        command = [
            "curl", "--fail", "--location", "--retry", "8", "--retry-all-errors",
            "--continue-at", "-", "--output", str(part), args.url,
        ]
        subprocess.run(command, check=True)
        if not gzip_ok(part):
            raise RuntimeError(f"Downloaded FASTQ failed gzip validation: {part}")
        observed_md5 = md5sum(part)
        if observed_md5.lower() != args.md5.lower():
            raise RuntimeError(f"MD5 mismatch for {args.sample}: {observed_md5} != {args.md5}")
        shutil.move(part, args.output)

    observed_reads = count_reads(args.output)
    if observed_reads != args.expected_reads:
        raise RuntimeError(
            f"Read-count mismatch for {args.sample}: {observed_reads} != {args.expected_reads}"
        )
    metrics = {
        "sample_id": args.sample,
        "url": args.url,
        "md5": args.md5.lower(),
        "expected_reads": args.expected_reads,
        "observed_reads": observed_reads,
        "bytes": args.output.stat().st_size,
    }
    args.metrics.write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
