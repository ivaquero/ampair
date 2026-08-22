#!/usr/bin/env python3
# =============================================================================
# genomes_download.py
#
# Download genomic, CDS, and RNA FASTA files for a bacterial genus with
# ncbi-genome-download.
# =============================================================================

import argparse
import csv
import gzip
import logging
import shutil
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path
from time import perf_counter

from common import configure_logging, fasta_directory_summary, sha256_file
from config_schema import load_config_file

log = logging.getLogger(__name__)
NCBI_ENTRYPOINT = (
    "import sys; from ncbi_genome_download import __main__ as n; sys.exit(n.main())"
)


def run_download(genus, assembly_level, fmt, out_dir):
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    cmd = [
        # Invoke the package through the active Python environment.  This is
        # reliable on Windows, where the generated console-script name and
        # PATH resolution differ from Unix.
        sys.executable,
        "-c",
        NCBI_ENTRYPOINT,
        "bacteria",
        "--genera",
        genus,
        "--assembly-levels",
        assembly_level,
        "--formats",
        fmt,
        "--flat-output",
        "--output-folder",
        out_dir,
        "--retries",
        "5",
        "--verbose",
    ]
    log.info("Running: %s", subprocess.list2cmdline(cmd))
    result = subprocess.run(  # noqa: S603 - fixed downloader with shell=False.
        cmd, capture_output=True, text=True
    )
    if result.stdout:
        log.info(result.stdout.rstrip())
    if result.stderr:
        log.info(result.stderr.rstrip())
    if result.returncode != 0:
        log.error("ncbi-genome-download failed with exit code %d", result.returncode)
        return result.returncode
    return 0


def reset_download_outputs(downloads, manifest):
    """Remove outputs managed by this rule before a fresh download."""
    for _, _, out_dir in downloads:
        path = Path(out_dir)
        if path.exists():
            log.info("Removing stale download directory: %s", path)
            shutil.rmtree(path)
        path.mkdir(parents=True, exist_ok=True)

    if manifest:
        manifest_path = Path(manifest)
        if manifest_path.exists():
            log.info("Removing stale download manifest: %s", manifest_path)
            manifest_path.unlink()


def decompress_gzip_files(*directories):
    for directory in directories:
        for gz_path in Path(directory).rglob("*.gz"):
            out_path = gz_path.with_suffix("")
            log.info("Decompressing %s -> %s", gz_path, out_path)
            with gzip.open(gz_path, "rb") as src, open(out_path, "wb") as dst:
                shutil.copyfileobj(src, dst)
            gz_path.unlink()


def write_manifest(path, genus, assembly_level, rows, config_sha256=""):
    if not path:
        return

    manifest_path = Path(path)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    generated_at = datetime.now(UTC).isoformat(timespec="seconds")
    fieldnames = [
        "generated_at",
        "genus",
        "assembly_level",
        "label",
        "format",
        "output_dir",
        "n_fna",
        "total_bytes",
        "data_fingerprint",
        "config_sha256",
    ]
    with manifest_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "generated_at": generated_at,
                    "genus": genus,
                    "assembly_level": assembly_level,
                    "config_sha256": config_sha256,
                    **row,
                }
            )


def parse_args():
    parser = argparse.ArgumentParser(description="Download genus FASTA files.")
    parser.add_argument("--config", help="Optional AmPrime config.yaml")
    parser.add_argument("--genus")
    parser.add_argument("--assembly-level")
    parser.add_argument("--genomic-dir", required=True)
    parser.add_argument("--cds-dir", required=True)
    parser.add_argument("--rna-dir", required=True)
    parser.add_argument("--manifest")
    parser.add_argument("--log", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    configure_logging(args.log)
    started = perf_counter()

    cfg = load_config_file(args.config) if args.config else {}
    config_sha256 = sha256_file(args.config) if args.config else ""
    genus = args.genus or cfg.get("genus")
    assembly_level = args.assembly_level or cfg.get("assembly_level")
    if not genus:
        raise SystemExit("missing --genus or config setting: genus")
    if not assembly_level:
        raise SystemExit("missing --assembly-level or config setting: assembly_level")

    downloads = [
        ("genomic", "fasta", args.genomic_dir),
        ("cds", "cds-fasta", args.cds_dir),
        ("rna", "rna-fasta", args.rna_dir),
    ]
    reset_download_outputs(downloads, args.manifest)

    for _, fmt, out_dir in downloads:
        code = run_download(genus, assembly_level, fmt, out_dir)
        if code != 0:
            return code

    decompress_gzip_files(args.genomic_dir, args.cds_dir, args.rna_dir)

    manifest_rows = []
    for label, fmt, out_dir in downloads:
        summary = fasta_directory_summary(out_dir)
        manifest_rows.append(
            {"label": label, "format": fmt, "output_dir": out_dir, **summary}
        )
    write_manifest(
        args.manifest,
        genus,
        assembly_level,
        manifest_rows,
        config_sha256=config_sha256,
    )

    n_gen = manifest_rows[0]["n_fna"]
    n_cds = manifest_rows[1]["n_fna"]
    n_rna = manifest_rows[2]["n_fna"]
    total_bytes = sum(int(row["total_bytes"]) for row in manifest_rows)
    elapsed = perf_counter() - started
    log.info(
        "Downloaded: %d genomic, %d CDS, %d RNA files for genus %s "
        "(%.1f MB total) in %.2f s",
        n_gen,
        n_cds,
        n_rna,
        genus,
        total_bytes / 1_000_000,
        elapsed,
    )
    if args.manifest:
        log.info("Wrote download manifest to %s", args.manifest)

    if n_gen == 0:
        log.error("No genomic FASTA downloaded. Check genus name and assembly level.")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
