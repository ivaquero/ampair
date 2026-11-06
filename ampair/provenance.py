"""Content fingerprints shared by the API and workflow scripts."""

from __future__ import annotations

import hashlib
import subprocess
from pathlib import Path


def sha256_file(path: str | Path) -> str:
    """Return the SHA-256 digest of a file without loading it into memory."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _run_git(repo_dir: Path, args: list[str]) -> subprocess.CompletedProcess | None:
    """Run a git command in *repo_dir*, or return None on any failure."""
    try:
        return subprocess.run(
            ["git", *args],
            cwd=repo_dir,
            capture_output=True,
            text=True,
            timeout=10,
        )
    except (OSError, subprocess.SubprocessError):
        return None


def git_revision(repo_dir: str | Path | None = None) -> str:
    """Return a short git revision string for provenance, or "" if unavailable.

    The result looks like ``a1b2c3d`` for a clean checkout, or ``a1b2c3d-dirty``
    when the working tree has uncommitted changes.  This is purely informational:
    if git is missing, the directory is not a repository, or the command fails or
    times out, an empty string is returned so provenance recording can never break
    the pipeline.
    """
    root = Path(repo_dir) if repo_dir is not None else Path(__file__).resolve().parent
    revision = _run_git(root, ["rev-parse", "--short", "HEAD"])
    if revision is None or revision.returncode != 0:
        return ""
    rev = revision.stdout.strip()
    if not rev:
        return ""
    status = _run_git(root, ["status", "--porcelain"])
    if status is not None and status.returncode == 0 and status.stdout.strip():
        return f"{rev}-dirty"
    return rev


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
