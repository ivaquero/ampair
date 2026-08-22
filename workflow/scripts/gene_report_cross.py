#!/usr/bin/env python3
"""Build a compact cross-gene comparison report from species summaries."""

from __future__ import annotations

import argparse
from html import escape
from pathlib import Path

from report_common import load_template, read_tsv, render_page

_HTML_TEMPLATE = load_template("gene_report_cross.html")


def _read_metrics(path: str) -> dict[str, str]:
    return {row["metric"]: row["value"] for row in read_tsv(path) if row.get("metric")}


def _percent(metrics: dict[str, str], key: str) -> str:
    try:
        return f"{float(metrics.get(key, 0)) * 100:.1f}%"
    except ValueError:
        return "0.0%"


def _gene_from_path(path: str) -> str:
    name = Path(path).name
    suffix = "_species_summary.tsv"
    return name.removesuffix(suffix) if name.endswith(suffix) else Path(name).stem


def _build_html(summary_paths: list[str]) -> str:
    rows = []
    for path in summary_paths:
        metrics = _read_metrics(path)
        rows.append(
            {
                "gene": _gene_from_path(path),
                "amplified_genomes": metrics.get("amplified_genomes", "0"),
                "amplification_rate": _percent(metrics, "amplification_rate"),
                "amplified_species": metrics.get("amplified_species", "0"),
                "multi_allele_genomes": metrics.get("multi_allele_genomes", "0"),
                "overlap_species": metrics.get("overlap_species", "0"),
                "overlap_rate": _percent(metrics, "overlap_rate"),
                "unique_amplicon_alleles": metrics.get("unique_amplicon_alleles", "0"),
            }
        )
    rows.sort(key=lambda row: row["gene"])
    body_rows = "\n".join(
        "<tr>" + "".join(f"<td>{escape(row[key])}</td>" for key in row) + "</tr>"
        for row in rows
    )
    headers = [
        ("gene", "Gene"),
        ("amplified_genomes", "Amplified genomes"),
        ("amplification_rate", "Amplification rate"),
        ("amplified_species", "Amplified species"),
        ("multi_allele_genomes", "Multiple-allele genomes"),
        ("overlap_species", "Overlap species"),
        ("overlap_rate", "Overlap rate"),
        ("unique_amplicon_alleles", "Unique alleles"),
    ]
    header_html = "".join(f"<th>{escape(label)}</th>" for _, label in headers)
    return render_page(_HTML_TEMPLATE, HEADER_HTML=header_html, BODY_ROWS=body_rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary-tsv", nargs="+", required=True)
    parser.add_argument("--out-html", required=True)
    args = parser.parse_args()
    Path(args.out_html).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_html).write_text(_build_html(args.summary_tsv), encoding="utf-8")


if __name__ == "__main__":
    main()
