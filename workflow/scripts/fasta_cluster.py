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
import subprocess
import sys
from time import perf_counter

from common import configure_logging, ensure_parent, log_process_output
from config_schema import load_config_file, resolve_clusterer
from dependencies import ensure_tool
from fasta_io import count_fasta_records, parse_fasta

log = logging.getLogger(__name__)

WARN_SEQUENCE_COUNT = 1000
WARN_CENTROID_COUNT = 500
# Generous wall-clock guard for a single VSEARCH invocation. Normal clustering
# finishes well within this; it exists to convert a wedged child process into a
# clear failure rather than an indefinite hang.
CLUSTER_TIMEOUT_SECONDS = 7_200


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
    try:
        completed = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=False,
            timeout=CLUSTER_TIMEOUT_SECONDS,
        )
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(
            f"VSEARCH clustering timed out after {CLUSTER_TIMEOUT_SECONDS} s"
        ) from exc
    log_process_output(completed, log)
    if completed.returncode != 0:
        raise RuntimeError(
            f"VSEARCH clustering failed with exit code {completed.returncode}"
        )


# Dispatch table from a configured backend name to its runner. A runner takes
# ``(executable, input_path, output_path, identity, threads)``. Only VSEARCH is
# implemented today; adding a backend means registering one entry here plus
# declaring its name in config_schema.SUPPORTED_CLUSTERERS, with no changes to
# the surrounding flow.
BACKEND_RUNNERS = {"vsearch": cluster_with_vsearch}


def parse_args():
    parser = argparse.ArgumentParser(description="Cluster FASTA records.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--identity", required=True, type=float)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument(
        "--config", help="AmPair config.yaml; selects the clusterer backend"
    )
    parser.add_argument("--log", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    if not 0 < args.identity <= 1:
        raise ValueError("--identity must be greater than 0 and at most 1")
    if args.threads < 1:
        raise ValueError("--threads must be a positive integer")
    configure_logging(args.log)
    ensure_parent(args.output)

    started = perf_counter()
    cfg = load_config_file(args.config) if args.config else {}
    clusterer = resolve_clusterer(cfg)
    executable = ensure_tool(clusterer)
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
        log.info("Large cluster input (%d sequences); using %s", n_in, clusterer)

    run_backend = BACKEND_RUNNERS[clusterer]
    run_backend(executable, args.input, args.output, args.identity, args.threads)
    n_out = count_fasta_records(args.output)

    if n_out > WARN_CENTROID_COUNT:
        log.warning(
            "Many centroids retained (%d). Alignment and primer design may be slow.",
            n_out,
        )

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
