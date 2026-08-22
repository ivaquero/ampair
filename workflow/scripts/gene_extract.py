#!/usr/bin/env python3
# =============================================================================
# gene_extract.py
#
# Extract target gene sequences from downloaded CDS and RNA FASTA files.
# Batch mode scans both directories once and writes one FASTA per configured
# gene; the single-gene mode remains available for debugging and reuse.
# =============================================================================

import argparse
import logging
import os
import re
from time import perf_counter

from common import configure_logging
from config_schema import load_config_file
from fasta_io import parse_fasta, write_fasta

log = logging.getLogger(__name__)
WARN_FASTA_FILE_COUNT = 1000


def _fna_files(directory, label):
    fna_files = sorted(
        os.path.join(directory, filename)
        for filename in os.listdir(directory)
        if filename.endswith(".fna")
    )
    if not fna_files:
        log.warning("No .fna files found in %s dir: %s", label, directory)
    return fna_files


def _matching_genes(header, names_by_gene):
    """Return all target genes matched by one FASTA header."""
    h = header.lower()
    gene_tag = re.search(r"\[gene=([^\]]+)\]", h)
    product_tag = re.search(r"\[product=([^\]]+)\]", h)
    gene_name = gene_tag.group(1).strip() if gene_tag else ""
    product = product_tag.group(1).strip() if product_tag else ""
    return [
        gene
        for gene, names in names_by_gene.items()
        if gene_name in names or any(name in product for name in names)
    ]


def scan_fasta_dirs(directories, names_by_gene):
    """Scan CDS/RNA FASTA files once and collect records for every gene."""
    started = perf_counter()
    extracted = {gene: [] for gene in names_by_gene}

    for label, directory in directories:
        fna_files = _fna_files(directory, label) if os.path.isdir(directory) else []
        log.info(
            "Scanning %d %s files for %d gene(s) ...",
            len(fna_files),
            label,
            len(names_by_gene),
        )
        if len(fna_files) > WARN_FASTA_FILE_COUNT:
            log.warning(
                "Large %s input (%d FASTA files). Gene extraction may be I/O-bound.",
                label,
                len(fna_files),
            )

        n_records = 0
        n_bp = 0
        hits_by_gene = dict.fromkeys(names_by_gene, 0)
        for fna in fna_files:
            for header, sequence in parse_fasta(fna):
                n_records += 1
                n_bp += len(sequence)
                for gene in _matching_genes(header, names_by_gene):
                    extracted[gene].append((header, sequence))
                    hits_by_gene[gene] += 1

        log.info(
            "Scanned %d %s FASTA record(s), %d bp total; hits by gene: %s",
            n_records,
            label,
            n_bp,
            ", ".join(f"{gene}={count}" for gene, count in hits_by_gene.items()),
        )

    elapsed = perf_counter() - started
    log.info(
        "Scanned %d gene target(s) across %d FASTA directories in %.2f s",
        len(names_by_gene),
        len(directories),
        elapsed,
    )
    return extracted


def _names_by_gene(cfg):
    aliases_by_gene = cfg.get("gene_aliases", {})
    return {
        gene: {
            gene.lower(),
            *(alias.lower() for alias in aliases_by_gene.get(gene, [])),
        }
        for gene in cfg.get("genes", [])
    }


def parse_args():
    parser = argparse.ArgumentParser(description="Extract target gene sequences.")
    parser.add_argument("--cds-dir", required=True)
    parser.add_argument("--rna-dir", required=True)
    parser.add_argument("--out-fasta")
    parser.add_argument("--gene")
    parser.add_argument(
        "--batch",
        action="store_true",
        help="Extract every gene in the config with one FASTA scan.",
    )
    parser.add_argument("--out-dir")
    parser.add_argument("--config", help="Optional AmPair config.yaml")
    parser.add_argument("--alias", action="append", default=[], dest="aliases")
    parser.add_argument("--log", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    configure_logging(args.log)
    cfg = load_config_file(args.config) if args.config else {}

    if args.batch:
        if not args.out_dir:
            raise SystemExit("--out-dir is required with --batch")
        names_by_gene = _names_by_gene(cfg)
        results_by_gene = scan_fasta_dirs(
            [("CDS", args.cds_dir), ("RNA", args.rna_dir)], names_by_gene
        )
        for gene, results in results_by_gene.items():
            if not results:
                log.warning(
                    "No sequences found for gene '%s'; writing empty FASTA", gene
                )
            out_fasta = os.path.join(args.out_dir, f"{gene}.fasta")
            write_fasta(results, out_fasta)
            log.info(
                "Written %d sequence(s) for %s -> %s", len(results), gene, out_fasta
            )
        return 0

    if not args.gene or not args.out_fasta:
        raise SystemExit("--gene and --out-fasta are required unless --batch is used")

    aliases = list(args.aliases)
    gene_aliases = cfg.get("gene_aliases", {})
    if isinstance(gene_aliases, dict):
        configured_aliases = gene_aliases.get(args.gene, [])
        if isinstance(configured_aliases, list):
            aliases.extend(
                alias for alias in configured_aliases if isinstance(alias, str)
            )
    search_names = {args.gene.lower(), *(alias.lower() for alias in aliases)}
    log.info("Target gene : %s", args.gene)
    log.info("Aliases     : %s", aliases)
    log.info("Search names: %s", sorted(search_names))
    log.info("CDS dir     : %s", args.cds_dir)
    log.info("RNA dir     : %s", args.rna_dir)

    results = scan_fasta_dirs(
        [("CDS", args.cds_dir), ("RNA", args.rna_dir)], {args.gene: search_names}
    )[args.gene]
    log.info("Total sequences extracted: %d", len(results))
    if not results:
        log.warning("No sequences found for gene '%s'; writing empty FASTA", args.gene)
    write_fasta(results, args.out_fasta)
    log.info("Written to %s", args.out_fasta)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
