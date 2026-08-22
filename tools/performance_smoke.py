#!/usr/bin/env python3
"""Run a small deterministic performance guard for batch PCR scanning."""

from __future__ import annotations

import csv
import os
import random
import subprocess
import sys
import tempfile
from pathlib import Path
from time import perf_counter

from _shared import ROOT, SCRIPTS, load_script_module, smoke_env

reverse_complement = load_script_module("common").reverse_complement

GENOME_COUNT = 4
GENOME_LENGTH = 20_000
PRIMER_COUNT = 4
DEFAULT_BUDGET_SECONDS = 30.0


def write_fixture(genome_dir: Path, primers_path: Path) -> None:
    target_fwd = "ACGTTGCAGTACCTGATCGA"
    target_rev = "TGCACGATCGGATCAGGTCA"
    rows = [
        {"primer_id": "primer_pair_1", "fwd": target_fwd, "rev": target_rev},
        {
            "primer_id": "primer_pair_2",
            "fwd": "GATCTAGCGTACGATCGTAC",
            "rev": "CAGTTCGATGGCATCGATGC",
        },
        {
            "primer_id": "primer_pair_3",
            "fwd": "TTACCGGATCGTACGGTACA",
            "rev": "AGTCGATCCGATGCTAGTCA",
        },
        {
            "primer_id": "primer_pair_4",
            "fwd": "CGATGACCTAGCTAGGCTAA",
            "rev": "TTCGGCATAGCTACGATGGA",
        },
    ]
    with primers_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["primer_id", "fwd", "rev", "combined_score"], delimiter="\t"
        )
        writer.writeheader()
        for rank, row in enumerate(rows, 1):
            writer.writerow({**row, "combined_score": 1.0 / rank})

    rng = random.Random(1729)
    binding = reverse_complement(target_rev)
    for index in range(GENOME_COUNT):
        sequence: list[str] = [str(rng.choice("ACGT")) for _ in range(GENOME_LENGTH)]
        start = 2_000 + index * 100
        sequence[start : start + len(target_fwd)] = list(target_fwd)
        reverse_start = start + 650
        sequence[reverse_start : reverse_start + len(binding)] = list(binding)
        sequence[8_000:9_000] = "N" * 1_000
        (genome_dir / f"genome_{index}.fna").write_text(
            f">genome_{index} [organism=Performance test species]\n"
            + "".join(sequence)
            + "\n",
            encoding="utf-8",
        )


def main() -> int:
    budget = float(os.environ.get("AMPAIR_PERFORMANCE_BUDGET", DEFAULT_BUDGET_SECONDS))
    with tempfile.TemporaryDirectory(prefix="ampair-performance-") as temp_dir:
        root = Path(temp_dir)
        genome_dir = root / "genomes"
        genome_dir.mkdir()
        primers = root / "primers.tsv"
        output = root / "amplicons.tsv"
        species_summary = root / "species_summary.tsv"
        species_table = root / "species.tsv"
        log_path = root / "in_silico_pcr.log"
        write_fixture(genome_dir, primers)

        command = [
            sys.executable,
            str(SCRIPTS / "in_silico_pcr.py"),
            "--primers-tsv",
            str(primers),
            "--genome-dir",
            str(genome_dir),
            "--out-tsv",
            str(output),
            "--gene",
            "performance",
            "--mismatch",
            "0",
            "--amplicon-min-len",
            "500",
            "--amplicon-max-len",
            "800",
            "--top-n",
            str(PRIMER_COUNT),
            "--workers",
            "2",
            "--species-summary",
            str(species_summary),
            "--species-tsv",
            str(species_table),
            "--log",
            str(log_path),
        ]
        env = smoke_env()
        started = perf_counter()
        completed = subprocess.run(
            command, cwd=ROOT, env=env, capture_output=True, text=True, check=False
        )
        elapsed = perf_counter() - started
        if completed.returncode != 0:
            print(completed.stdout, file=sys.stdout, end="")
            print(completed.stderr, file=sys.stderr, end="")
            raise RuntimeError(
                f"Performance fixture failed with exit code {completed.returncode}"
            )
        if not output.is_file() or output.stat().st_size == 0:
            raise RuntimeError("Performance fixture did not produce amplicon output")
        if elapsed > budget:
            raise RuntimeError(
                f"Performance smoke exceeded {budget:.1f} s: {elapsed:.2f} s"
            )
        print(
            f"performance smoke ok: genomes={GENOME_COUNT} "
            f"length={GENOME_LENGTH} candidates={PRIMER_COUNT} "
            f"elapsed={elapsed:.2f}s budget={budget:.1f}s"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
