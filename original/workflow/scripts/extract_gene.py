#!/usr/bin/env python3
# =============================================================================
# extract_gene.py
#
# Extracts sequences for a target gene from all CDS and RNA FASTA files
# in the genome download directories, writes them to a single FASTA file.
#
# Strategy:
#   1. Search CDS FASTA files first (protein-coding genes: rpoB, tuf, etc.)
#   2. Search RNA FASTA files second (rRNA genes: 16S, 23S, etc.)
#   3. Warn and skip individual files where the gene is not found
#   4. If no sequences found across all files, exit with error
#
# Header matching: looks for [gene=<name>] tag in NCBI CDS/RNA FASTA headers.
# Also checks [product=...] as fallback for rRNA genes.
#
# Usage (called via Snakemake script: directive):
#   snakemake.params.cds_dir : CDS directory (from checkpoint, via input[0])
#   snakemake.params.rna_dir : RNA directory (from checkpoint, via input[1])
#   snakemake.output[0]      : output FASTA path
#   snakemake.params.gene    : gene name (e.g. "rpoB", "tuf", "16S")
#   snakemake.params.aliases : list of alternative names for the gene
#   snakemake.log[0]         : log file path
# =============================================================================

import os
import sys
import re
import logging

# ---------------------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------------------
log_path = snakemake.log[0]
os.makedirs(os.path.dirname(log_path), exist_ok=True)

logging.basicConfig(
    filename=log_path,
    level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
log = logging.getLogger()

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
cds_dir   = snakemake.params.cds_dir
rna_dir   = snakemake.params.rna_dir
out_fasta = snakemake.output[0]
gene      = snakemake.params.gene
aliases   = snakemake.params.get("aliases", [])

# Build set of all names to search for (gene name + aliases), lowercased
search_names = set([gene.lower()] + [a.lower() for a in aliases])

log.info("Target gene : %s", gene)
log.info("Aliases     : %s", aliases)
log.info("Search names: %s", sorted(search_names))
log.info("CDS dir     : %s", cds_dir)
log.info("RNA dir     : %s", rna_dir)

# ---------------------------------------------------------------------------
# FASTA parser
# ---------------------------------------------------------------------------
def parse_fasta(filepath):
    """Yield (header, sequence) tuples from a FASTA file."""
    header = None
    seq_parts = []
    with open(filepath, "r") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq_parts)
                header = line
                seq_parts = []
            else:
                seq_parts.append(line)
    if header is not None:
        yield header, "".join(seq_parts)

# ---------------------------------------------------------------------------
# Header matching
# ---------------------------------------------------------------------------
def header_matches(header, names):
    """
    Return True if the FASTA header contains any of the target gene names.
    Checks:
      [gene=<name>]           — standard CDS/RNA FASTA tag
      [product=<...name...>]  — fallback for rRNA (e.g. "16S ribosomal RNA")
    """
    h = header.lower()

    # [gene=rpoB] style
    gene_tag = re.search(r'\[gene=([^\]]+)\]', h)
    if gene_tag and gene_tag.group(1).strip() in names:
        return True

    # [product=...] style — check if any search name appears in the product string
    product_tag = re.search(r'\[product=([^\]]+)\]', h)
    if product_tag:
        product = product_tag.group(1).strip()
        if any(name in product for name in names):
            return True

    return False

# ---------------------------------------------------------------------------
# Extract from a directory of FASTA files
# ---------------------------------------------------------------------------
def extract_from_dir(directory, names, label):
    """
    Scan all .fna files in directory, return list of (header, seq) matches.
    Logs a warning for files where the gene is not found.
    """
    fna_files = sorted([
        os.path.join(directory, f)
        for f in os.listdir(directory)
        if f.endswith(".fna")
    ])

    if not fna_files:
        log.warning("No .fna files found in %s dir: %s", label, directory)
        return []

    log.info("Scanning %d %s files ...", len(fna_files), label)

    extracted = []
    for fna in fna_files:
        genome_id = os.path.basename(fna)
        hits = [
            (hdr, seq)
            for hdr, seq in parse_fasta(fna)
            if header_matches(hdr, names)
        ]
        if hits:
            log.info("  %s : %d sequence(s) found", genome_id, len(hits))
            extracted.extend(hits)
        else:
            log.warning("  %s : gene '%s' not found — skipping", genome_id, gene)

    return extracted

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
results = []

# Search CDS first
if os.path.isdir(cds_dir):
    results.extend(extract_from_dir(cds_dir, search_names, "CDS"))
else:
    log.warning("CDS directory does not exist: %s", cds_dir)

# Search RNA if nothing found yet (or always, to catch rRNA genes)
if os.path.isdir(rna_dir):
    results.extend(extract_from_dir(rna_dir, search_names, "RNA"))
else:
    log.warning("RNA directory does not exist: %s", rna_dir)

log.info("Total sequences extracted: %d", len(results))

if len(results) == 0:
    log.error(
        "No sequences found for gene '%s' in any genome. "
        "Check gene name, aliases, or try a less stringent assembly_level.",
        gene
    )
    sys.exit(1)

# Write output FASTA
os.makedirs(os.path.dirname(out_fasta), exist_ok=True)
with open(out_fasta, "w") as fh:
    for header, seq in results:
        fh.write(header + "\n")
        # wrap sequence at 80 chars
        for i in range(0, len(seq), 80):
            fh.write(seq[i:i+80] + "\n")

log.info("Written to %s", out_fasta)
