# =============================================================================
# gene_extract.smk
#
# Batch rule: extracts all target gene sequences from the downloaded CDS and
# RNA FASTA files in one pass, producing one merged FASTA per gene.
#
# Depends on: rule genomes_download. Genome directories and config are passed
# to the CLI as plain paths.
#
# Output : results/{genus}/extracted/{gene}.fasta   [temp, one per config gene]
# =============================================================================

rule gene_extract:
    threads: 1
    resources:
        mem_mb=2048
    input:
        cds = str(GENOMES_CDS),
        rna = str(GENOMES_RNA),
        config = CONFIG_FILE
    output:
        fasta = [temp(str(EXTRACTED / layout.extracted_fasta(gene))) for gene in GENES]
    params:
        config_file = CONFIG_FILE,
        cds_dir = str(GENOMES_CDS),
        rna_dir = str(GENOMES_RNA),
        out_dir = str(EXTRACTED),
        script = str(SCRIPTS_DIR / "gene_extract.py")
    log:
        str(RESULTS / "logs" / "gene_extract.log")
    benchmark:
        str(RESULTS / "benchmarks" / "gene_extract.txt")
    shell:
        """
        python {params.script:q} \
            --cds-dir {params.cds_dir:q} \
            --rna-dir {params.rna_dir:q} \
            --out-dir {params.out_dir:q} \
            --batch \
            --config {params.config_file:q} \
            --log {log:q}
        """
