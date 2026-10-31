#!/usr/bin/env python3
"""Package compact analysis results and a SHA256 checksum."""

from __future__ import annotations

import argparse
import hashlib
import tarfile
from pathlib import Path


EXCLUDED_SUFFIXES = {"CHGCAR", "AECCAR0", "AECCAR1", "AECCAR2", "WAVECAR", "vasprun.xml", "POTCAR"}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True, help="Results root to archive")
    parser.add_argument("--output", type=Path, required=True, help="Output .tgz path")
    parser.add_argument("--include-snapshots", action="store_true", help="Include selected ACF.dat and POSCAR files")
    args = parser.parse_args()
    root = args.root.expanduser().resolve()
    archive = args.output.expanduser().resolve()
    if not root.is_dir():
        parser.error(f"results root is not a directory: {root}")
    archive.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "w:gz") as tar:
        for path in sorted(root.rglob("*")):
            if not path.is_file() or path == archive:
                continue
            if path.name in EXCLUDED_SUFFIXES:
                continue
            if not args.include_snapshots and path.name in {"ACF.dat", "POSCAR"} and any(part.casefold().startswith("bader") for part in path.parts):
                continue
            relative = path.relative_to(root)
            tar.add(path, arcname=relative.as_posix(), recursive=False)
    checksum = _sha256(archive)
    checksum_path = Path(str(archive) + ".sha256")
    checksum_path.write_text(f"{checksum}  {archive.name}\n", encoding="utf-8")
    print(f"Archive: {archive}\nSHA256: {checksum_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
