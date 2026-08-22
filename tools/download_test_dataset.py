#!/usr/bin/env python3
"""Download the functional-test genomes and package them as a tar archive."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
import tarfile
import tempfile
import time
from pathlib import Path

import yaml
from _shared import extend_pythonpath, sha256_file, snapshot_metadata_path

ROOT = Path(__file__).resolve().parents[1]
DOWNLOAD_SCRIPT = ROOT / "workflow" / "scripts" / "genomes_download.py"
DEFAULT_CONFIG = ROOT / "config" / "config.yaml"
DEFAULT_OUTPUT = ROOT / "data" / "borrelia-genomes.tar.gz"
SNAPSHOT_SCHEMA_VERSION = 1
DOWNLOAD_ATTEMPTS = 3
DOWNLOAD_RETRY_DELAY_SECONDS = 5


def _snapshot_id(manifest: Path, genus: str, assembly_level: str) -> str:
    with manifest.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    stable_rows = [
        {
            key: row.get(key, "")
            for key in ("label", "format", "n_fna", "total_bytes", "data_fingerprint")
        }
        for row in rows
    ]
    payload = json.dumps(
        {"genus": genus, "assembly_level": assembly_level, "rows": stable_rows},
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _write_snapshot_metadata(
    archive: Path, manifest: Path, genus: str, assembly_level: str, config: Path
) -> None:
    metadata = {
        "schema_version": SNAPSHOT_SCHEMA_VERSION,
        "snapshot_id": _snapshot_id(manifest, genus, assembly_level),
        "archive_sha256": sha256_file(archive),
        "archive_size": archive.stat().st_size,
        "manifest_sha256": sha256_file(manifest),
        "config_sha256": sha256_file(config),
        "genus": genus,
        "assembly_level": assembly_level,
    }
    destination = snapshot_metadata_path(archive)
    partial = destination.with_name(destination.name + ".partial")
    partial.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    os.replace(partial, destination)


def download_archive(output: Path, genus: str, config: Path) -> None:
    if not DOWNLOAD_SCRIPT.is_file():
        raise FileNotFoundError(
            f"Genome downloader script not found: {DOWNLOAD_SCRIPT}. "
            "Check out the current repository revision before running CI."
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    partial_output = output.with_name(output.name + ".partial")
    if partial_output.exists():
        partial_output.unlink()

    # Keep downloader paths on the same filesystem as the checkout. On
    # Windows, ncbi-genome-download computes relative paths internally and
    # fails when the system temp directory is on another drive.
    with tempfile.TemporaryDirectory(
        prefix=".ampair-test-data-", dir=ROOT
    ) as temp_dir:
        root = Path(temp_dir)
        genomes = root / "genomes"
        log_path = root / "download.log"
        command = [
            sys.executable,
            str(DOWNLOAD_SCRIPT),
            "--config",
            str(config),
            "--genus",
            genus,
            "--genomic-dir",
            str(genomes / "genomic"),
            "--cds-dir",
            str(genomes / "cds"),
            "--rna-dir",
            str(genomes / "rna"),
            "--manifest",
            str(genomes / "download_manifest.tsv"),
            "--log",
            str(log_path),
        ]
        completed = None
        for attempt in range(1, DOWNLOAD_ATTEMPTS + 1):
            try:
                completed = subprocess.run(
                    command, cwd=ROOT, check=False, env=_project_env()
                )
            except OSError as exc:
                if attempt == DOWNLOAD_ATTEMPTS:
                    raise RuntimeError(
                        f"Genome test-data downloader could not start: {exc}"
                    ) from exc
                print(
                    f"Downloader start failed ({exc}); retrying "
                    f"{attempt}/{DOWNLOAD_ATTEMPTS}...",
                    file=sys.stderr,
                )
                time.sleep(DOWNLOAD_RETRY_DELAY_SECONDS)
                continue

            if completed.returncode == 0:
                break

            if log_path.is_file():
                print(log_path.read_text(encoding="utf-8"), file=sys.stderr)
            if attempt < DOWNLOAD_ATTEMPTS:
                print(
                    f"Genome test-data download failed with exit code "
                    f"{_format_returncode(completed.returncode)}; retrying "
                    f"{attempt}/{DOWNLOAD_ATTEMPTS}...",
                    file=sys.stderr,
                )
                time.sleep(DOWNLOAD_RETRY_DELAY_SECONDS)

        if completed is None or completed.returncode != 0:
            returncode = (
                _format_returncode(completed.returncode)
                if completed is not None
                else "unknown"
            )
            raise RuntimeError(
                f"Genome test-data download failed with exit code {returncode}. "
                "The downloader log was printed above."
            )

        manifest = genomes / "download_manifest.tsv"
        with tarfile.open(partial_output, "w:gz") as archive:
            archive.add(genomes, arcname="genomes")
        os.replace(partial_output, output)
        _write_snapshot_metadata(
            output, manifest, genus, _assembly_level(config), config
        )


def _assembly_level(config: Path) -> str:
    values = yaml.safe_load(config.read_text(encoding="utf-8")) or {}
    return str(values.get("assembly_level", "complete"))


def _format_returncode(returncode: int) -> str:
    """Make Windows' unsigned representation of negative exits diagnosable."""
    if returncode >= 0x80000000:
        signed = returncode - 0x100000000
        return f"{signed} (Windows 0x{returncode:08X})"
    return str(returncode)


def _project_env() -> dict[str, str]:
    return extend_pythonpath(DOWNLOAD_SCRIPT.parent)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genus", default="Borrelia")
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    download_archive(args.output, args.genus, args.config)
    print(f"downloaded test dataset: {args.output}")


if __name__ == "__main__":
    main()
