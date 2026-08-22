#!/usr/bin/env python3
# =============================================================================
# fasta_align.py
#
# Align centroid FASTA records with MUSCLE. On Windows, a missing executable is
# installed through Scoop when Scoop is available. A single input sequence is
# copied unchanged because no multiple sequence alignment is required.
# =============================================================================

import argparse
import csv
import logging
import os
import shutil
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path
from time import perf_counter

from common import configure_logging, log_process_output
from dependencies import ensure_tool
from fasta_io import count_fasta_records, parse_fasta

log = logging.getLogger(__name__)
WARN_ALIGNMENT_SEQUENCE_COUNT = 500
WARN_ALIGNMENT_BP = 2_000_000


def _backend_version(executable):
    try:
        result = subprocess.run(
            [executable, "--version"],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return ""

    text = (result.stdout or result.stderr).strip().splitlines()
    return text[0] if text else ""


def write_alignment_metadata(path, row):
    if not path:
        return

    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "generated_at",
        "requested_backend",
        "backend_used",
        "fallback_used",
        "backend_executable",
        "backend_version",
        "n_input_sequences",
        "n_output_sequences",
        "input_total_bp",
        "elapsed_seconds",
    ]
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerow(
            {"generated_at": datetime.now(UTC).isoformat(timespec="seconds"), **row}
        )


def run_muscle(executable, input_path, output_path, threads):
    """Run MUSCLE and return ``(succeeded, used_legacy_syntax)``."""
    output_tmp = Path(output_path).with_suffix(Path(output_path).suffix + ".tmp")
    commands = [
        [
            executable,
            "-align",
            input_path,
            "-output",
            str(output_tmp),
            "-threads",
            str(threads),
        ],
        [
            executable,
            "-in",
            input_path,
            "-out",
            str(output_tmp),
            "-threads",
            str(threads),
        ],
    ]
    for command_index, command in enumerate(commands):
        if output_tmp.exists():
            output_tmp.unlink()
        log.info("Running MUSCLE: %s", " ".join(str(part) for part in command))
        completed = subprocess.run(command, capture_output=True, text=True, check=False)
        log_process_output(completed, log)
        if (
            completed.returncode == 0
            and output_tmp.exists()
            and output_tmp.stat().st_size
        ):
            os.replace(output_tmp, output_path)
            return True, command_index == 1

    if output_tmp.exists():
        output_tmp.unlink()
    return False, False


def parse_args():
    parser = argparse.ArgumentParser(description="Align FASTA records with MUSCLE.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--config", help="Accepted for workflow compatibility")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--metadata", help="Optional alignment metadata TSV")
    parser.add_argument("--log", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.threads < 1:
        raise ValueError("--threads must be a positive integer")
    configure_logging(args.log)
    os.makedirs(os.path.dirname(args.output), exist_ok=True)

    started = perf_counter()
    executable = ensure_tool("muscle")
    n_in = 0
    total_bp = 0
    min_len = None
    max_len = None
    for _, seq in parse_fasta(args.input):
        n_in += 1
        seq_len = len(seq)
        total_bp += seq_len
        min_len = seq_len if min_len is None else min(min_len, seq_len)
        max_len = seq_len if max_len is None else max(max_len, seq_len)

    if n_in < 2:
        shutil.copyfile(args.input, args.output)
        log.info("Only %d sequence(s); skipped MUSCLE alignment", n_in)
        elapsed = perf_counter() - started
        write_alignment_metadata(
            args.metadata,
            {
                "requested_backend": "muscle",
                "backend_used": "skipped",
                "fallback_used": False,
                "backend_executable": executable,
                "backend_version": _backend_version(executable),
                "n_input_sequences": n_in,
                "n_output_sequences": n_in,
                "input_total_bp": total_bp,
                "elapsed_seconds": round(elapsed, 3),
            },
        )
        return 0

    log.info(
        "Input size: %d sequence(s), %d bp total, length range %d-%d bp",
        n_in,
        total_bp,
        min_len,
        max_len,
    )
    if n_in > WARN_ALIGNMENT_SEQUENCE_COUNT or total_bp > WARN_ALIGNMENT_BP:
        log.info("Large alignment input; using MUSCLE")

    alignment_ok, fallback_used = run_muscle(
        executable, args.input, args.output, args.threads
    )
    if not alignment_ok:
        raise SystemExit("MUSCLE alignment failed; see log for command output")

    n_out = count_fasta_records(args.output)
    elapsed = perf_counter() - started
    write_alignment_metadata(
        args.metadata,
        {
            "requested_backend": "muscle",
            "backend_used": "muscle",
            "fallback_used": fallback_used,
            "backend_executable": executable,
            "backend_version": _backend_version(executable),
            "n_input_sequences": n_in,
            "n_output_sequences": n_out,
            "input_total_bp": total_bp,
            "elapsed_seconds": round(elapsed, 3),
        },
    )
    log.info("Aligned %d sequences with MUSCLE in %.2f s", n_out, elapsed)
    return 0


if __name__ == "__main__":
    sys.exit(main())
