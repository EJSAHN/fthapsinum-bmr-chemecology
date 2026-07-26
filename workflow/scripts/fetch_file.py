#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import os
import shutil
import tempfile
from pathlib import Path
from urllib.request import Request, urlopen


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description="Download a public file with atomic replacement.")
    parser.add_argument("--url", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sha256")
    parser.add_argument("--minimum-bytes", type=int, default=1)
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.output.exists() and args.output.stat().st_size >= args.minimum_bytes:
        if not args.sha256 or sha256(args.output).lower() == args.sha256.lower():
            return

    descriptor, temporary_name = tempfile.mkstemp(prefix=args.output.name + ".", suffix=".part", dir=args.output.parent)
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        request = Request(args.url, headers={"User-Agent": "Mozilla/5.0"})
        with urlopen(request, timeout=120) as response, temporary.open("wb") as output:
            shutil.copyfileobj(response, output, length=1024 * 1024)
        if temporary.stat().st_size < args.minimum_bytes:
            raise RuntimeError(f"Downloaded file is smaller than {args.minimum_bytes} bytes")
        if args.sha256 and sha256(temporary).lower() != args.sha256.lower():
            raise RuntimeError("SHA-256 checksum mismatch")
        temporary.replace(args.output)
    finally:
        temporary.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
