"""Shared helpers for the AmPair HTML report builders.

The per-gene and cross-gene reports both read tab-separated inputs and render
an HTML template, so those small utilities live here instead of being copied.
"""

from __future__ import annotations

import csv
from pathlib import Path

_HTML_DIR = Path(__file__).resolve().parent


def read_tsv(path) -> list[dict]:
    """Read a tab-separated file into a list of dict rows."""
    with open(path, newline="", encoding="utf-8") as fh:
        return [dict(row) for row in csv.DictReader(fh, delimiter="\t")]


def load_template(name: str) -> str:
    """Return the raw contents of an HTML template in this directory."""
    return (_HTML_DIR / name).read_text(encoding="utf-8")


def render_page(template: str, **placeholders) -> str:
    """Replace ``{KEY}`` placeholders in ``template`` with the given values.

    Values are inserted verbatim; callers are responsible for escaping when
    needed (e.g. ``html.escape`` on user-controlled text).
    """
    html = template
    for key, value in placeholders.items():
        html = html.replace("{" + key + "}", str(value))
    return html


__all__ = ["load_template", "read_tsv", "render_page"]
