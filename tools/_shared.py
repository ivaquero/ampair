"""Shared helpers for the tools/ helper scripts.

These scripts run outside the installed package, so this module exposes a
single source of truth for the small utilities they would otherwise duplicate.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from amprime.provenance import sha256_file  # noqa: E402  (sys.path set above)


def snapshot_metadata_path(archive: Path) -> Path:
    """Return the JSON sidecar path paired with a test-data archive."""
    return archive.with_name(archive.name + ".json")


def extend_pythonpath(scripts_dir: str | Path, env: dict[str, str] | None = None) -> dict[str, str]:
    """Return ``env`` with ``scripts_dir`` prepended to PYTHONPATH."""
    env = dict(os.environ) if env is None else dict(env)
    pythonpath = env.get("PYTHONPATH")
    env["PYTHONPATH"] = (
        str(scripts_dir)
        if not pythonpath
        else f"{scripts_dir}{os.pathsep}{pythonpath}"
    )
    return env


__all__ = ["extend_pythonpath", "sha256_file", "snapshot_metadata_path"]
