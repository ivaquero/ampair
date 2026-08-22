#!/usr/bin/env python3
"""Small shared helpers for AmPair workflow command-line tools."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from ampair.provenance import fasta_directory_summary, sha256_file

__all__ = [
    "IUPAC_COMPLEMENT_TABLE",
    "config_param",
    "configure_logging",
    "fasta_directory_summary",
    "log_process_output",
    "required_param",
    "reverse_complement",
    "sha256_file",
]
IUPAC_COMPLEMENT_TABLE = str.maketrans(
    "ACGTRYMKSWHBVDNacgtrymkswhbvdn", "TGCAYRKMSWDVBHNtgcayrkmswdvbhn"
)


def configure_logging(log_path: str) -> None:
    Path(log_path).parent.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        filename=log_path,
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def config_param(cli_value, cfg: dict, key: str):
    return cli_value if cli_value is not None else cfg.get(key)


def required_param(name: str, value):
    if value is None:
        raise SystemExit(
            f"missing --{name.replace('_', '-')} or config setting: {name}"
        )
    return value


def reverse_complement(seq: str) -> str:
    return seq.translate(IUPAC_COMPLEMENT_TABLE)[::-1]


def log_process_output(completed, logger) -> None:
    """Log a completed subprocess's stdout and stderr at INFO level."""
    if completed.stdout:
        logger.info(completed.stdout.rstrip())
    if completed.stderr:
        logger.info(completed.stderr.rstrip())
