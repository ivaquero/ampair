#!/usr/bin/env python3
# =============================================================================
# fasta_cluster.py
#
# Dereplicate FASTA records with VSEARCH. On Windows, a missing executable is
# installed through Scoop when Scoop is available. Empty input produces an
# empty centroid FASTA so downstream per-gene reporting can continue.
# =============================================================================

import argparse
import logging
import os
import subprocess
import sys
from time import perf_counter

from common import configure_logging, log_process_output
from dependencies import ensure_tool
from fasta_io import count_fasta_records, parse_fasta

log = logging.getLogger(__name__)

WARN_SEQUENCE_COUNT = 1000
WARN_CENTROID_COUNT = 500


def cluster_with_vsearch(executable, input_path, output_path, identity, threads):
    command = [
        executable,
        "--cluster_fast",
        input_path,
        "--id",
        f"{identity:.6f}",
        "--centroids",
        output_path,
        "--minseqlength",
        "1",
        "--threads",
        str(threads),
    ]
    log.info("Running VSEARCH: %s", " ".join(command))
    completed = subprocess.run(  # noqa: S603 - executable came from PATH lookup.
        command,
        capture_output=True,
        text=True,
        check=False,
    )
    log_process_output(completed, log)
    if completed.returncode != 0:
        raise RuntimeError(
            f"VSEARCH clustering failed with exit code {completed.returncode}"
        )


def parse_args():
    parser = argparse.ArgumentParser(description="Cluster FASTA records.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--identity", required=True, type=float)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--log", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    if not 0 < args.identity <= 1:
        raise ValueError("--identity must be greater than 0 and at most 1")
    if args.threads < 1:
        raise ValueError("--threads must be a positive integer")
    configure_logging(args.log)
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    started = perf_counter()
    executable = ensure_tool("vsearch")
    n_in = 0
    total_bp = 0
    for _, seq in parse_fasta(args.input):
        n_in += 1
        total_bp += len(seq)
    if n_in == 0:
        open(args.output, "w", encoding="utf-8").close()
        log.info("No input sequences; wrote empty centroids FASTA")
        return 0

    log.info("Input size: %d sequence(s), %d bp total", n_in, total_bp)

    if n_in > WARN_SEQUENCE_COUNT:
        log.info("Large cluster input (%d sequences); using VSEARCH", n_in)

    cluster_with_vsearch(
        executable, args.input, args.output, args.identity, args.threads
    )
    n_centroids = count_fasta_records(args.output)

    if n_centroids > WARN_CENTROID_COUNT:
        log.warning(
            "Many centroids retained (%d). Alignment and primer design may be slow.",
            n_centroids,
        )

    n_out = count_fasta_records(args.output)
    elapsed = perf_counter() - started
    log.info(
        "Clustered %d sequences into %d centroids with VSEARCH in %.2f s",
        n_in,
        n_out,
        elapsed,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
