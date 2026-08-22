#!/usr/bin/env python3
# =============================================================================
# gene_report.py - self-contained HTML report for a single gene
#
# Standalone CLI used by Snakemake. Requires markdown (conda-forge).
# Builds the body in Markdown, converts to HTML, and wraps in a styled page.
# =============================================================================

import argparse
import base64
import logging
import os
from datetime import UTC, datetime
from html import escape

from common import sha256_file
from config_schema import load_config_file
from report_common import load_template, read_tsv, render_page

# =============================================================================
# HTML page shell - loaded from separate .html file
# =============================================================================
_PAGE = load_template("gene_report.html")

# =============================================================================
# Markdown engine - tables turned on
# =============================================================================
_MD_EXTENSIONS = ["tables"]


# =============================================================================
# Helpers
# =============================================================================
def _b64_png(path):
    with open(path, "rb") as fh:
        return base64.b64encode(fh.read()).decode()


def _num(row, key, default=0.0):
    value = row.get(key, default)
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _md_cell(value):
    return escape(str(value)).replace("|", r"\|")


def _md_code(value):
    return f"`{_md_cell(value)}`"


def _find_primer_by_id(primers, primer_id):
    for primer in primers:
        if primer.get("primer_id") == primer_id:
            return primer
    return None


def _recommended_primer(primers, amplicons):
    if not primers:
        return None
    if not amplicons:
        return primers[0]

    best_pcr = amplicons[0]
    primer = _find_primer_by_id(primers, best_pcr.get("primer_id"))
    return primer or primers[0]


def _build_body(
    gene,
    genus,
    timestamp,
    primers,
    top_primer,
    pcr,
    diversity_img,
    alignment_meta=None,
    config_sha256="",
    manifest_sha256="",
    data_fingerprints=None,
    species_metrics=None,
):
    """Build the report body as a Markdown string, using inline HTML only
    for styled elements (alerts, badges).  All control flow is plain Python."""

    out = []
    gene_text = _md_cell(gene)
    genus_text = _md_cell(genus)

    # ---- title -----------------------------------------------------------
    out.append(f"# {genus_text} - *{gene_text}* primer design")
    out.append("")
    out.append(f"*Generated {timestamp}*")
    out.append("")

    # ---- provenance ------------------------------------------------------
    out.append("## Provenance")
    out.append("")
    out.append("|  |  |")
    out.append("|---|---|")
    out.append(f"| **Config SHA-256** | {_md_code(config_sha256 or 'unavailable')} |")
    out.append(
        f"| **Download manifest SHA-256** | "
        f"{_md_code(manifest_sha256 or 'unavailable')} |"
    )
    for label, fingerprint in sorted((data_fingerprints or {}).items()):
        out.append(
            f"| **{_md_cell(label)} data fingerprint** | {_md_code(fingerprint)} |"
        )
    out.append("")

    # ---- species-level validation ---------------------------------------
    out.append("## Species-level validation")
    out.append("")
    metrics = species_metrics or {}
    out.append("| Metric | Value |")
    out.append("|---|---|")
    metric_labels = [
        ("total_genomes", "Total genomes"),
        ("amplified_genomes", "Amplified genomes"),
        ("amplification_rate", "Amplification rate"),
        ("total_species", "Total species"),
        ("amplified_species", "Amplified species"),
        ("multi_allele_genomes", "Genomes with multiple amplicon alleles"),
        ("multi_allele_rate", "Multiple-allele genome rate"),
        ("overlap_species", "Species with inter-species allele overlap"),
        ("overlap_rate", "Inter-species overlap rate"),
        ("unique_amplicon_alleles", "Unique amplicon alleles"),
    ]
    for key, label in metric_labels:
        value = metrics.get(key, "0")
        if key.endswith("rate"):
            value = f"{_num({key: value}, key) * 100:.1f}%"
        out.append(f"| **{_md_cell(label)}** | {_md_cell(value)} |")
    out.append("")

    # ---- alignment metadata ---------------------------------------------
    if alignment_meta:
        meta = alignment_meta
        out.append("## Alignment")
        out.append("")
        out.append("|  |  |")
        out.append("|---|---|")
        out.append(
            f"| **Requested backend** | {_md_code(meta.get('requested_backend', ''))} |"
        )
        out.append(f"| **Backend used** | {_md_code(meta.get('backend_used', ''))} |")
        if meta.get("fallback_used"):
            out.append(f"| **Fallback used** | {_md_cell(meta.get('fallback_used'))} |")
        if meta.get("n_input_sequences"):
            out.append(
                f"| **Input sequences** | {_md_cell(meta.get('n_input_sequences'))} |"
            )
        if meta.get("n_output_sequences"):
            out.append(
                f"| **Output sequences** | {_md_cell(meta.get('n_output_sequences'))} |"
            )
        if meta.get("elapsed_seconds"):
            out.append(
                f"| **Elapsed seconds** | {_md_cell(meta.get('elapsed_seconds'))} |"
            )
        out.append("")

    # ---- in silico PCR ---------------------------------------------------
    out.append("## In silico PCR validation")
    out.append("")
    if pcr:
        a = pcr
        rate = f"{_num(a, 'amplification_rate') * 100:.1f}%"
        out.append("|  |  |")
        out.append("|---|---|")
        out.append(f"| **Primer pair** | {_md_code(a.get('primer_id', ''))} |")
        if a.get("validation_rank"):
            out.append(
                f"| **Validation rank** | {_md_cell(a.get('validation_rank'))} |"
            )
        if a.get("input_rank"):
            out.append(
                f"| **Original candidate rank** | {_md_cell(a.get('input_rank'))} |"
            )
        out.append(f"| **Forward primer** | {_md_code(a.get('fwd', ''))} |")
        out.append(f"| **Reverse primer** | {_md_code(a.get('rev', ''))} |")
        out.append(
            f"| **Genomes amplified** | {_md_cell(a.get('n_genomes_amplified', ''))} / "
            f"{_md_cell(a.get('total_genomes', ''))}"
            f' <span class="badge badge-green">{rate}</span> |'
        )
        out.append(
            f"| **Mean amplicon length** | {_num(a, 'mean_amplicon_len'):.0f} bp |"
        )
    else:
        out.append('<div class="alert alert-info">')
        out.append(
            "<strong>No primer pair was validated.</strong> "
            "Either no primers passed the design filters, or in silico PCR "
            "could not run. See the candidate table below."
        )
        out.append("</div>")
    out.append("")

    # ---- recommended pair ------------------------------------------------
    out.append("## Recommended primer pair")
    out.append("")
    if top_primer:
        t = top_primer
        out.append("| Property | Forward | Reverse |")
        out.append("|---|---|---|")
        fwd_seq = _md_code(t.get("fwd", ""))
        rev_seq = _md_code(t.get("rev", ""))
        out.append(f"| Sequence | {fwd_seq} | {rev_seq} |")
        out.append(
            f"| Position | {_md_cell(t.get('fwd_pos', ''))} | "
            f"{_md_cell(t.get('rev_pos', ''))} |"
        )
        out.append(
            f"| GC fraction | {_num(t, 'fwd_GC'):.2f} | {_num(t, 'rev_GC'):.2f} |"
        )
        out.append("")
        out.append(
            f"**Amplicon length:** {_md_cell(t.get('amplicon_len', ''))} bp &ensp;"
            f"**Pair diversity (Shannon):** {_num(t, 'pair_diversity'):.3f} &ensp;"
            f"**GC difference:** {_num(t, 'delta_GC'):.3f} &ensp;"
            f"**Combined score:** {_num(t, 'combined_score'):.3f}"
        )
    else:
        out.append('<div class="alert alert-warn">')
        out.append(
            "<strong>No candidate primers were found.</strong> "
            "Try raising <code>div_cut</code> or relaxing <code>GC_tol</code> "
            "in <code>config.yaml</code>."
        )
        out.append("</div>")
    out.append("")

    # ---- diversity plot --------------------------------------------------
    out.append("## Sequence diversity")
    out.append("")
    if diversity_img:
        out.append(
            "Per-position Shannon entropy across the alignment. "
            "Low-entropy (conserved) regions are good primer targets. "
            "Red rectangles mark the binding sites of the top-scoring primer pairs."
        )
        out.append("")
        out.append(
            f'<img src="data:image/png;base64,{diversity_img}" '
            f'alt="Diversity plot for {escape(str(gene), quote=True)}">'
        )
    else:
        out.append('<div class="alert alert-info">Diversity plot not available.</div>')
    out.append("")

    # ---- all candidates --------------------------------------------------
    out.append("## All candidate primer pairs")
    out.append("")
    if primers:
        keys = list(primers[0].keys())
        out.append("| " + " | ".join(_md_cell(key) for key in keys) + " |")
        out.append("|" + "|".join(["---"] * len(keys)) + "|")
        for row in primers:
            vals = [_md_cell(row.get(k, "")) for k in keys]
            out.append("| " + " | ".join(vals) + " |")
    else:
        out.append('<div class="alert alert-info">No candidates to display.</div>')
    out.append("")

    # ---- footer ----------------------------------------------------------
    out.append(f"*Primer pipeline - {genus_text}*")

    return "\n".join(out)


# =============================================================================
# Main
# =============================================================================
def parse_args():
    parser = argparse.ArgumentParser(description="Build a per-gene HTML report.")
    parser.add_argument("--gene", required=True)
    parser.add_argument("--genus")
    parser.add_argument("--config", help="Optional AmPrime config.yaml")
    parser.add_argument("--primers-tsv", required=True)
    parser.add_argument("--amplicons-tsv", required=True)
    parser.add_argument("--diversity-png", required=True)
    parser.add_argument("--alignment-meta")
    parser.add_argument("--download-manifest")
    parser.add_argument("--species-summary")
    parser.add_argument("--out-html", required=True)
    parser.add_argument("--log", required=True)
    return parser.parse_args()


def main():
    args = parse_args()

    import markdown

    log_path = args.log
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    logging.basicConfig(
        filename=log_path,
        level=logging.INFO,
        format="%(asctime)s  %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    log = logging.getLogger()

    gene = args.gene
    cfg = load_config_file(args.config) if args.config else {}
    genus = args.genus or cfg.get("genus")
    if not genus:
        raise SystemExit("missing --genus or config setting: genus")

    # -- read inputs -------------------------------------------------------
    primers = read_tsv(args.primers_tsv) if os.path.isfile(args.primers_tsv) else []
    amplicons = (
        read_tsv(args.amplicons_tsv) if os.path.isfile(args.amplicons_tsv) else []
    )
    alignment_meta = (
        read_tsv(args.alignment_meta)
        if args.alignment_meta and os.path.isfile(args.alignment_meta)
        else []
    )
    manifest_rows = (
        read_tsv(args.download_manifest)
        if args.download_manifest and os.path.isfile(args.download_manifest)
        else []
    )
    species_metrics = {}
    if args.species_summary and os.path.isfile(args.species_summary):
        species_metrics = {
            row.get("metric", ""): row.get("value", "")
            for row in read_tsv(args.species_summary)
            if row.get("metric")
        }
    data_fingerprints = {
        row.get("label", ""): row.get("data_fingerprint", "")
        for row in manifest_rows
        if row.get("label") and row.get("data_fingerprint")
    }
    config_sha256 = (
        sha256_file(args.config) if args.config and os.path.isfile(args.config) else ""
    )
    manifest_sha256 = (
        sha256_file(args.download_manifest)
        if args.download_manifest and os.path.isfile(args.download_manifest)
        else ""
    )
    log.info(
        "primers: %d rows, amplicons: %d rows, alignment metadata: %d rows",
        len(primers),
        len(amplicons),
        len(alignment_meta),
    )

    has_diversity = os.path.isfile(args.diversity_png)

    # -- build markdown body ------------------------------------------------
    body_md = _build_body(
        gene=gene,
        genus=genus,
        timestamp=datetime.now(UTC).strftime("%Y-%m-%d %H:%M"),
        primers=primers,
        top_primer=_recommended_primer(primers, amplicons),
        pcr=amplicons[0] if amplicons else None,
        diversity_img=_b64_png(args.diversity_png) if has_diversity else None,
        alignment_meta=alignment_meta[0] if alignment_meta else None,
        config_sha256=config_sha256,
        manifest_sha256=manifest_sha256,
        data_fingerprints=data_fingerprints,
        species_metrics=species_metrics,
    )

    # -- markdown to HTML, then wrap in page shell --------------------------
    body_html = markdown.markdown(body_md, extensions=_MD_EXTENSIONS)
    title = f"{genus} - {gene} primer design"
    html = render_page(_PAGE, TITLE=escape(title), CONTENT=body_html)

    # -- write --------------------------------------------------------------
    out_html = args.out_html
    os.makedirs(os.path.dirname(out_html), exist_ok=True)
    with open(out_html, "w", encoding="utf-8") as fh:
        fh.write(html)

    log.info("Wrote report -> %s", out_html)


if __name__ == "__main__":
    main()
