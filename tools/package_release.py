#!/usr/bin/env python3
"""Create release archives from the repository source tree."""

import os
import tarfile
import tomllib
import zipfile

from _shared import ROOT

DIST = ROOT / "dist"
INCLUDE_PATHS = [
    "README.md",
    "LICENSE",
    "Snakefile",
    "pyproject.toml",
    "pixi.toml",
    "pixi.lock",
    "amprime",
    "config",
    "docs",
    "workflow",
    "tools",
]


def project_version():
    with open(ROOT / "pixi.toml", "rb") as fh:
        manifest = tomllib.load(fh)
    return os.environ.get("GITHUB_REF_NAME") or f"v{manifest['workspace']['version']}"


def iter_files():
    for rel in INCLUDE_PATHS:
        path = ROOT / rel
        if not path.exists():
            continue
        if path.is_file():
            yield path
            continue
        for child in sorted(path.rglob("*")):
            if child.is_file() and "__pycache__" not in child.parts:
                yield child


def main():
    version = project_version()
    archive_stem = f"amprime-{version}"
    DIST.mkdir(exist_ok=True)

    zip_path = DIST / f"{archive_stem}.zip"
    tar_path = DIST / f"{archive_stem}.tar.gz"

    files = list(iter_files())

    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in files:
            zf.write(path, f"{archive_stem}/{path.relative_to(ROOT).as_posix()}")

    with tarfile.open(tar_path, "w:gz") as tf:
        for path in files:
            tf.add(path, f"{archive_stem}/{path.relative_to(ROOT).as_posix()}")

    print(zip_path.relative_to(ROOT).as_posix())
    print(tar_path.relative_to(ROOT).as_posix())


if __name__ == "__main__":
    main()
