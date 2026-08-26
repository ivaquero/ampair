#!/usr/bin/env python3
# =============================================================================
# in_silico_pcr.py
#
# Validates the top-ranked primer pair against all full genome assemblies
# using `seqkit amplicon`. Reports how many genomes are successfully
# amplified, the amplification rate, and the mean amplicon length.
#
# Strategy:
#   1. Read {gene}_primers.tsv, take the top row (highest combined_score)
#   2. Write a seqkit primer file: name<tab>fwd<tab>rev
#   3. Run seqkit amplicon on each *.fna genome in genomes/genomic/
#      with -m <mismatch> allowed mismatches
#   4. Count amplified genomes, compute rate and mean amplicon length
#   5. Write summary to {gene}_amplicons.tsv
#
# Handles empty primer TSV gracefully (writes empty summary, exits 0).
#
# Snakemake interface:
#   snakemake.input.primers   : {gene}_primers.tsv
#   snakemake.input.genome_dir: genomes/genomic/ directory
#   snakemake.output[0]       : {gene}_amplicons.tsv
#   snakemake.params.gene     : gene name
#   snakemake.params.mismatch : allowed mismatches (int)
#   snakemake.params.amplicon_min_len / amplicon_max_len : valid product window
#   snakemake.log[0]          : log file
# =============================================================================

import os
import sys
import csv
import glob
import subprocess
import logging
import tempfile

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
log_path = snakemake.log[0]
os.makedirs(os.path.dirname(log_path), exist_ok=True)
logging.basicConfig(
    filename=log_path, level=logging.INFO,
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)
log = logging.getLogger()

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
primers_tsv = snakemake.input.primers
genome_dir  = snakemake.input.genome_dir
out_tsv     = snakemake.output[0]
gene        = snakemake.params.gene
mismatch    = int(snakemake.params.mismatch)

# Valid-product length window: discard spurious long-range "amplicons" that
# arise from mismatched off-target priming. We allow a margin around the
# configured amplicon size range.
min_len     = int(snakemake.params.amplicon_min_len)
max_len_cfg = int(snakemake.params.amplicon_max_len)
len_margin  = 100  # bp of slack on each side
valid_lo    = max(0, min_len - len_margin)
valid_hi    = max_len_cfg + len_margin

log.info("Gene        : %s", gene)
log.info("Primers TSV : %s", primers_tsv)
log.info("Genome dir  : %s", genome_dir)
log.info("Mismatch    : %d", mismatch)

# ---------------------------------------------------------------------------
# Output header
# ---------------------------------------------------------------------------
OUT_COLS = [
    "primer_id", "fwd", "rev",
    "n_genomes_amplified", "total_genomes",
    "amplification_rate", "mean_amplicon_len"
]

def write_summary(rows):
    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    with open(out_tsv, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(OUT_COLS)
        for r in rows:
            w.writerow(r)

# ---------------------------------------------------------------------------
# 1. Read top primer pair
# ---------------------------------------------------------------------------
with open(primers_tsv) as fh:
    reader = csv.DictReader(fh, delimiter="\t")
    primer_rows = list(reader)

if len(primer_rows) == 0:
    log.warning("Primer TSV is empty - no primers to validate. Writing empty summary.")
    write_summary([])
    sys.exit(0)

top = primer_rows[0]   # already sorted by combined_score in design_primers.R
primer_id = top["primer_id"]
fwd = top["fwd"]
rev = top["rev"]
log.info("Top primer pair: %s  fwd=%s  rev=%s", primer_id, fwd, rev)

# ---------------------------------------------------------------------------
# 2. Locate genome files
# ---------------------------------------------------------------------------
genomes = sorted(glob.glob(os.path.join(genome_dir, "*.fna")))
total_genomes = len(genomes)
log.info("Found %d genome files", total_genomes)

if total_genomes == 0:
    log.error("No genome files found in %s", genome_dir)
    write_summary([])
    sys.exit(1)

# ---------------------------------------------------------------------------
# 3. Write seqkit primer file
# ---------------------------------------------------------------------------
with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as pf:
    pf.write("%s\t%s\t%s\n" % (primer_id, fwd, rev))
    primer_file = pf.name
log.info("Wrote seqkit primer file: %s", primer_file)

# ---------------------------------------------------------------------------
# 4. Run seqkit amplicon per genome, collect amplicon lengths
# ---------------------------------------------------------------------------
amplified_genomes = 0
amplicon_lengths  = []

for g in genomes:
    genome_id = os.path.basename(g)
    cmd = [
        "seqkit", "amplicon",
        "-p", primer_file,
        "-m", str(mismatch),
        "--bed",
        g
    ]
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, check=True
        )
    except subprocess.CalledProcessError as e:
        log.warning("  %s : seqkit failed (%s) - skipping", genome_id, e.stderr.strip())
        continue

    lines = [ln for ln in result.stdout.splitlines() if ln.strip()]

    # Keep only products whose length falls in the valid window; long-range
    # products from off-target priming are discarded as spurious.
    valid_here = []
    for ln in lines:
        fields = ln.split("\t")
        if len(fields) >= 3:
            try:
                amp_len = int(fields[2]) - int(fields[1])
            except ValueError:
                continue
            if valid_lo <= amp_len <= valid_hi:
                valid_here.append(amp_len)

    if valid_here:
        amplified_genomes += 1
        amplicon_lengths.extend(valid_here)
        log.info("  %s : amplified (%d valid product[s], %d raw)",
                 genome_id, len(valid_here), len(lines))
    elif lines:
        log.info("  %s : only spurious products (%d raw, none in %d-%d bp) - not counted",
                 genome_id, len(lines), valid_lo, valid_hi)
    else:
        log.info("  %s : no amplification", genome_id)

os.unlink(primer_file)

# ---------------------------------------------------------------------------
# 5. Compute summary and write
# ---------------------------------------------------------------------------
rate = amplified_genomes / total_genomes if total_genomes else 0.0
mean_len = sum(amplicon_lengths) / len(amplicon_lengths) if amplicon_lengths else 0

log.info("Amplified %d / %d genomes (rate %.3f), mean amplicon length %.1f bp",
         amplified_genomes, total_genomes, rate, mean_len)

write_summary([[
    primer_id, fwd, rev,
    amplified_genomes, total_genomes,
    round(rate, 4), round(mean_len, 1)
]])

log.info("Written to %s", out_tsv)
