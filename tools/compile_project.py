#!/usr/bin/env python3
"""Compile all project Python files."""

import py_compile

from _shared import ROOT

PYTHON_FILES = sorted(
    [
        *(ROOT / "workflow" / "scripts").glob("*.py"),
        *(ROOT / "tools").glob("*.py"),
        *(ROOT / "amprime").glob("*.py"),
    ]
)


def main():
    for path in PYTHON_FILES:
        py_compile.compile(str(path), doraise=True)
        print(f"compiled {path.relative_to(ROOT).as_posix()}")


if __name__ == "__main__":
    main()
