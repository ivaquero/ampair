"""Shared helpers for the tools/ helper scripts.

These scripts run outside the installed package, so this module exposes a
single source of truth for the small utilities they would otherwise duplicate.
"""

from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "workflow" / "scripts"

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

from ampair.provenance import sha256_file  # noqa: E402  (sys.path set above)


def snapshot_metadata_path(archive: Path) -> Path:
    """Return the JSON sidecar path paired with a test-data archive."""
    return archive.with_name(archive.name + ".json")


def extend_pythonpath(
    scripts_dir: str | Path, env: dict[str, str] | None = None
) -> dict[str, str]:
    """Return ``env`` with ``scripts_dir`` prepended to PYTHONPATH."""
    env = dict(os.environ) if env is None else dict(env)
    pythonpath = env.get("PYTHONPATH")
    env["PYTHONPATH"] = (
        str(scripts_dir) if not pythonpath else f"{scripts_dir}{os.pathsep}{pythonpath}"
    )
    return env


def smoke_env() -> dict[str, str]:
    """Environment for running workflow scripts as subprocesses."""
    return extend_pythonpath(SCRIPTS)


def load_script_module(module_name: str):
    """Load a workflow/scripts module by name, registering it in ``sys.modules``."""
    module_path = SCRIPTS / f"{module_name}.py"
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {module_name} from {module_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


__all__ = [
    "extend_pythonpath",
    "load_script_module",
    "sha256_file",
    "smoke_env",
    "snapshot_metadata_path",
]
