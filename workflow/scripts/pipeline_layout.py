"""Single machine-readable source for the per-genus result layout and rules.

Both the Snakemake workflow (``workflow/Snakefile`` + ``workflow/rules/*.smk``)
and the thin Python API (``ampair/api.py``) read these names, so the output
filename contract and the pipeline rule sequence live in exactly one place.
Renaming a result artifact or reordering a stage is a one-line change here.

Keep this module free of third-party and workflow imports so it is safe to load
from the Snakefile, the rule modules, and the packaged API alike.
"""

from __future__ import annotations

# Top-level results subdirectories under ``results/{genus}/``. The Snakefile
# composes full Paths from these; the API joins them the same way.
GENOMES = "genomes"
EXTRACTED = "extracted"
ALIGNED = "aligned"
PRIMERS = "primers"
REPORTS = "reports"

# Genus-level (not per-gene) artifacts.
GENOMIC = "genomic"
CDS = "cds"
RNA = "rna"
DOWNLOAD_MANIFEST = "download_manifest.tsv"
CROSS_REPORT = "gene_report_cross.html"


# Per-gene filename templates. ``gene`` is either a concrete gene name (used by
# list comprehensions over ``config['genes']``) or the Snakemake wildcard
# placeholder string ``"{gene}"`` (used inside rule input/output declarations);
# both produce the identical file name, so the layout never drifts.
def extracted_fasta(gene: str) -> str:
    return f"{gene}.fasta"


def centroids(gene: str) -> str:
    return f"{gene}.centroids.fasta"


def aln(gene: str) -> str:
    return f"{gene}.aln"


def alignment_meta(gene: str) -> str:
    return f"{gene}.alignment.tsv"


def primers_raw(gene: str) -> str:
    return f"{gene}_primers_raw.tsv"


def diversity(gene: str) -> str:
    return f"{gene}_diversity.png"


def primers(gene: str) -> str:
    return f"{gene}_primers.tsv"


def amplicons(gene: str) -> str:
    return f"{gene}_amplicons.tsv"


def species_summary(gene: str) -> str:
    return f"{gene}_species_summary.tsv"


def species(gene: str) -> str:
    return f"{gene}_species.tsv"


def report_html(gene: str) -> str:
    return f"{gene}_report.html"


# Pipeline stage rule names, in dependency order. ``PER_GENE_RULES`` is the
# sequence the API forces a rerun of for a single-gene report; ``CROSS_RULES``
# is the separate cross-gene comparison invocation. Names must match the
# ``rule`` declarations in workflow/rules/*.smk.
PER_GENE_RULES = [
    "gene_extract",
    "cluster",
    "align",
    "primers_design",
    "primers_check",
    "in_silico_pcr",
    "gene_report",
]
CROSS_RULES = ["gene_report_cross"]
