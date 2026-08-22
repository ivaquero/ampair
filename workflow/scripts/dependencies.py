#!/usr/bin/env python3
"""Detect and install external AmPair command-line dependencies."""

from __future__ import annotations

import argparse
import logging
import os
import re
import shutil
import subprocess
import tarfile
from pathlib import Path

from common import log_process_output

log = logging.getLogger(__name__)
REQUIRED_TOOLS = ("vsearch", "muscle", "seqkit")
SCOOP_BUCKET_NAME = "main-plus"
SCOOP_BUCKET_URL = "https://github.com/Scoopforge/Main-Plus"


def _run_scoop(scoop: str, arguments: list[str]) -> subprocess.CompletedProcess[str]:
    command = [scoop, *arguments]
    return subprocess.run(command, capture_output=True, text=True, check=False)


def _seqkit_runs(executable: str) -> bool:
    completed = subprocess.run(
        [executable, "version"], capture_output=True, text=True, check=False
    )
    return completed.returncode == 0


def _repair_seqkit_archive(scoop: str) -> str | None:
    """Repair Scoop's historical tarball-as-exe SeqKit installation."""
    scoop_path = Path(scoop).resolve()
    scoop_root = Path(os.environ.get("SCOOP", scoop_path.parent.parent))
    installed = scoop_root / "apps" / "seqkit" / "current" / "seqkit.exe"
    if not installed.is_file():
        return None

    with installed.open("rb") as fh:
        if fh.read(2) != b"\x1f\x8b":
            return None

    with tarfile.open(installed, mode="r:gz") as archive:
        member = next(
            (
                item
                for item in archive.getmembers()
                if item.isfile() and Path(item.name).name.lower() == "seqkit.exe"
            ),
            None,
        )
        if member is None:
            return None
        source = archive.extractfile(member)
        if source is None:
            return None
        repaired = installed.with_suffix(".repaired")
        repaired.write_bytes(source.read())

    os.replace(repaired, installed)
    log.info("Repaired Scoop SeqKit archive at %s", installed)
    return str(installed)


def _validate_tool(tool: str, executable: str, scoop: str | None = None) -> str:
    if tool != "seqkit" or _seqkit_runs(executable):
        return executable
    if scoop is not None:
        repaired = _repair_seqkit_archive(scoop)
        if repaired is not None and _seqkit_runs(repaired):
            return repaired
    raise RuntimeError(
        "SeqKit was found on PATH but could not run 'seqkit version'. "
        "Repair or reinstall it with Scoop, then retry."
    )


def ensure_scoop_bucket(scoop: str) -> None:
    """Ensure the Scoop bucket containing AmPair's Windows tools is loaded."""
    listed = _run_scoop(scoop, ["bucket", "list"])
    log_process_output(listed, log)
    if listed.returncode != 0:
        raise RuntimeError(
            f"Scoop bucket listing failed with exit code {listed.returncode}."
        )

    bucket_output = re.sub(
        r"\x1b\[[0-?]*[ -/]*[@-~]", "", f"{listed.stdout}\n{listed.stderr}"
    )
    bucket_loaded = any(
        line.strip().split(maxsplit=1)[0:1] == [SCOOP_BUCKET_NAME]
        for line in bucket_output.splitlines()
        if line.strip()
    )
    if bucket_loaded:
        log.info("Scoop bucket already loaded: %s", SCOOP_BUCKET_NAME)
        return

    log.info("Loading Scoop bucket: %s", SCOOP_BUCKET_NAME)
    added = _run_scoop(scoop, ["bucket", "add", SCOOP_BUCKET_NAME, SCOOP_BUCKET_URL])
    log_process_output(added, log)
    if added.returncode != 0:
        raise RuntimeError(
            f"Scoop failed to add bucket {SCOOP_BUCKET_NAME} with exit code "
            f"{added.returncode}."
        )


def _ensure_windows_tool(tool: str) -> str:
    """Install and return a missing tool through Scoop on Windows."""
    scoop = shutil.which("scoop")
    if scoop is None:
        raise RuntimeError(
            f"{tool} is required on Windows, but Scoop was not found. "
            f"Install Scoop, then run 'scoop install {tool}'."
        )
    assert scoop is not None

    log.info("%s not found; installing it with Scoop", tool)
    ensure_scoop_bucket(scoop)
    completed = _run_scoop(scoop, ["install", tool])
    log_process_output(completed, log)
    if completed.returncode != 0:
        raise RuntimeError(
            f"Scoop failed to install {tool} with exit code "
            f"{completed.returncode}. Run 'scoop install {tool}' manually."
        )

    executable = shutil.which(tool)
    if executable is None:
        raise RuntimeError(
            f"Scoop reported a successful {tool} installation, but {tool} was "
            "not found on PATH. Restart the shell and retry."
        )
    assert executable is not None
    return _validate_tool(tool, executable, scoop)


def ensure_tool(tool: str) -> str:
    """Return a tool path without probing Scoop on Unix platforms."""
    executable = shutil.which(tool)
    if executable is not None:
        scoop = shutil.which("scoop") if os.name == "nt" else None
        return _validate_tool(tool, executable, scoop)

    if os.name == "nt":
        return _ensure_windows_tool(tool)

    raise RuntimeError(
        f"{tool} is required but was not found on PATH. Install it with "
        "Pixi/Conda and retry."
    )


def ensure_required_tools() -> dict[str, str]:
    """Ensure every external executable required by the workflow is present."""
    return {tool: ensure_tool(tool) for tool in REQUIRED_TOOLS}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    tools = ensure_required_tools()
    for name, executable in tools.items():
        log.info("%s: %s", name, executable)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
