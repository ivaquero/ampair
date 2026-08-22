#!/usr/bin/env python3
"""Verify the cached functional-test archive and its snapshot metadata."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from _shared import sha256_file, snapshot_metadata_path


def verify(archive: Path) -> dict[str, object]:
    metadata_path = snapshot_metadata_path(archive)
    if not archive.is_file():
        raise FileNotFoundError(f"Test dataset archive not found: {archive}")
    if not metadata_path.is_file():
        raise FileNotFoundError(f"Test dataset metadata not found: {metadata_path}")

    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    required = {
        "schema_version",
        "snapshot_id",
        "archive_sha256",
        "archive_size",
        "genus",
        "assembly_level",
    }
    missing = sorted(required - metadata.keys())
    if missing:
        raise ValueError(
            f"Test dataset metadata is missing fields: {', '.join(missing)}"
        )
    if metadata["schema_version"] != 1:
        raise ValueError(
            f"Unsupported test dataset metadata schema: {metadata['schema_version']}"
        )

    actual_sha256 = sha256_file(archive)
    if actual_sha256 != metadata["archive_sha256"]:
        raise ValueError(
            f"Test dataset checksum mismatch: expected {metadata['archive_sha256']}, "
            f"got {actual_sha256}"
        )
    actual_size = archive.stat().st_size
    if actual_size != metadata["archive_size"]:
        raise ValueError(
            f"Test dataset size mismatch: expected {metadata['archive_size']}, "
            f"got {actual_size}"
        )
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    args = parser.parse_args()
    metadata = verify(args.archive)
    print(
        "test dataset ok: "
        f"snapshot={metadata['snapshot_id']} "
        f"genus={metadata['genus']} "
        f"assembly_level={metadata['assembly_level']}"
    )


if __name__ == "__main__":
    main()
