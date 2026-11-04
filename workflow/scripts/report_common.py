"""Shared helpers for the AmPair HTML report builders.

The per-gene and cross-gene reports both read tab-separated inputs and render
an HTML template, so those small utilities live here instead of being copied.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

_HTML_DIR = Path(__file__).resolve().parent

# Matches the report templates' relative stylesheet link, e.g.
# ``<link rel="stylesheet" type="text/css" href="report.css">``.
_CSS_LINK = re.compile(r"<link[^>]*href=\"report\.css\"[^>]*>", re.IGNORECASE)


def read_tsv(path) -> list[dict]:
    """Read a tab-separated file into a list of dict rows."""
    with open(path, newline="", encoding="utf-8") as fh:
        return [dict(row) for row in csv.DictReader(fh, delimiter="\t")]


def load_template(name: str) -> str:
    """Return the raw contents of an HTML template in this directory."""
    return (_HTML_DIR / name).read_text(encoding="utf-8")


def load_css(name: str = "report.css") -> str:
    """Return the bundled stylesheet used by the HTML reports."""
    return (_HTML_DIR / name).read_text(encoding="utf-8")


def inline_css(html: str, css_text: str) -> str:
    """Replace a relative ``report.css`` link with an inline ``<style>`` block.

    Reports are written under ``results/{genus}/reports/`` where ``report.css``
    is not copied, so a relative link leaves the delivered HTML unstyled.
    Inlining keeps each report a single self-contained file.
    """
    style = f"<style>\n{css_text}\n</style>"
    # A function replacement avoids re-interpreting backslashes in the CSS body.
    return _CSS_LINK.sub(lambda _match: style, html)


def render_page(template: str, **placeholders) -> str:
    """Replace ``{KEY}`` placeholders in ``template`` with the given values.

    Values are inserted verbatim; callers are responsible for escaping when
    needed (e.g. ``html.escape`` on user-controlled text).
    """
    html = template
    for key, value in placeholders.items():
        html = html.replace("{" + key + "}", str(value))
    return html


__all__ = [
    "inline_css",
    "load_css",
    "load_template",
    "read_tsv",
    "render_page",
]
