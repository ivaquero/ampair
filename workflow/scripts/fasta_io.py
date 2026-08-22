#!/usr/bin/env python3
"""Small FASTA I/O helpers shared by AmPrime command-line tools."""

from functools import lru_cache
from pathlib import Path


def parse_fasta(path):
    """Yield (header, sequence) tuples from a FASTA file.

    Header strings keep the leading ">" to preserve the legacy script behavior.
    """
    header = None
    seq_parts = []
    with open(path, encoding="utf-8-sig") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        yield header, "".join(seq_parts)


def write_fasta(records, path):
    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as fh:
        for header, seq in records:
            fh.write(header + "\n")
            fh.writelines(seq[i : i + 80] + "\n" for i in range(0, len(seq), 80))


@lru_cache(maxsize=128)
def _count_cached(path: str, mtime_ns: int, size: int) -> int:
    return sum(1 for _ in parse_fasta(path))


def count_fasta_records(path):
    """Count FASTA records, caching by path + mtime + size.

    Files produced by this process (cluster/align outputs) are counted right
    after being written; the cache avoids re-parsing the same file multiple
    times while still invalidating on any content change.
    """
    stat = Path(path).stat()
    return _count_cached(str(path), stat.st_mtime_ns, stat.st_size)
