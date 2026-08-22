"""Content fingerprints shared by the API and workflow scripts."""

from __future__ import annotations

import hashlib
from pathlib import Path


def sha256_file(path: str | Path) -> str:
    """Return the SHA-256 digest of a file without loading it into memory."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fasta_directory_summary(directory: str | Path) -> dict[str, int | str]:
    """Summarize and fingerprint all FASTA files in a directory.

    The fingerprint includes each relative path and file content, so a rename,
    addition, deletion, or sequence change is visible in the result provenance.
    """
    root = Path(directory)
    paths = sorted(path for path in root.rglob("*.fna") if path.is_file())
    digest = hashlib.sha256()
    total_bytes = 0
    for path in paths:
        relative = path.relative_to(root).as_posix().encode("utf-8")
        content_digest = sha256_file(path).encode("ascii")
        digest.update(relative)
        digest.update(b"\0")
        digest.update(content_digest)
        digest.update(b"\n")
        total_bytes += path.stat().st_size
    return {
        "n_fna": len(paths),
        "total_bytes": total_bytes,
        "data_fingerprint": digest.hexdigest(),
    }
