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

from common import configure_logging, ensure_parent, log_process_output
from config_schema import load_config_file, resolve_aligner
from dependencies import ensure_tool
from fasta_io import count_fasta_records, parse_fasta

log = logging.getLogger(__name__)
WARN_ALIGNMENT_SEQUENCE_COUNT = 500
WARN_ALIGNMENT_BP = 2_000_000
# Generous wall-clock guard for a single MUSCLE invocation. Legitimate
# alignments finish within this; it converts a wedged child process into a
# clean failure instead of an indefinite hang.
ALIGN_TIMEOUT_SECONDS = 7_200


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
        try:
            completed = subprocess.run(
                command,
                capture_output=True,
                text=True,
                check=False,
                timeout=ALIGN_TIMEOUT_SECONDS,
            )
        except subprocess.TimeoutExpired:
            log.info(
                "MUSCLE attempt timed out after %d s: %s",
                ALIGN_TIMEOUT_SECONDS,
                " ".join(str(part) for part in command),
            )
            continue
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


# Dispatch table from a configured backend name to its runner. A runner takes
# ``(executable, input_path, output_path, threads)`` and returns
# ``(succeeded, fallback_used)``. Only MUSCLE is implemented today; adding a
# backend means registering one entry here plus declaring its name in
# config_schema.SUPPORTED_ALIGNERS, with no changes to the surrounding flow.
BACKEND_RUNNERS = {"muscle": run_muscle}


def parse_args():
    parser = argparse.ArgumentParser(description="Align FASTA records with MUSCLE.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--config", help="AmPair config.yaml; selects the aligner backend"
    )
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--metadata", help="Optional alignment metadata TSV")
    parser.add_argument("--log", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    if args.threads < 1:
        raise ValueError("--threads must be a positive integer")
    configure_logging(args.log)
    ensure_parent(args.output)

    started = perf_counter()
    cfg = load_config_file(args.config) if args.config else {}
    aligner = resolve_aligner(cfg)
    executable = ensure_tool(aligner)
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
        log.info("Only %d sequence(s); skipped %s alignment", n_in, aligner)
        elapsed = perf_counter() - started
        write_alignment_metadata(
            args.metadata,
            {
                "requested_backend": aligner,
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
        log.info("Large alignment input; using %s", aligner)

    run_backend = BACKEND_RUNNERS[aligner]
    alignment_ok, fallback_used = run_backend(
        executable, args.input, args.output, args.threads
    )
    if not alignment_ok:
        raise SystemExit(f"{aligner} alignment failed; see log for command output")

    n_out = count_fasta_records(args.output)
    elapsed = perf_counter() - started
    write_alignment_metadata(
        args.metadata,
        {
            "requested_backend": aligner,
            "backend_used": aligner,
            "fallback_used": fallback_used,
            "backend_executable": executable,
            "backend_version": _backend_version(executable),
            "n_input_sequences": n_in,
            "n_output_sequences": n_out,
            "input_total_bp": total_bp,
            "elapsed_seconds": round(elapsed, 3),
        },
    )
    log.info("Aligned %d sequences with %s in %.2f s", n_out, aligner, elapsed)
    return 0


if __name__ == "__main__":
    sys.exit(main())
