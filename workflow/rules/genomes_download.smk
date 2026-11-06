# =============================================================================
# genomes_download.smk
#
# Downloads three FASTA types for the target genus from NCBI, each into its
# own subdirectory:
#   genomic   - full genome assemblies, used by in_silico_pcr (KEPT)
#   cds       - CDS FASTA, used to extract protein-coding genes
#   rna       - RNA FASTA, used to extract rRNA genes
#
# Uses ncbi-genome-download with --genera (no taxid resolution needed).
#
# Outputs:
#   results/{genus}/genomes/genomic/   - *_genomic.fna           [KEPT]
#   results/{genus}/genomes/cds/       - *_cds_from_genomic.fna
#   results/{genus}/genomes/rna/       - *_rna_from_genomic.fna
#   results/{genus}/genomes/download_manifest.tsv
# =============================================================================


rule genomes_download:
    threads: 1
    resources:
        mem_mb=1024
    input:
        config = CONFIG_FILE
    output:
        genomic = directory(str(GENOMES_GENOMIC)),
        cds     = directory(str(GENOMES_CDS)),
        rna     = directory(str(GENOMES_RNA)),
        manifest = str(RESULTS / layout.GENOMES / layout.DOWNLOAD_MANIFEST)
    params:
        config_file = CONFIG_FILE,
        script = str(SCRIPTS_DIR / "genomes_download.py")
    log:
        str(RESULTS / "logs" / "genomes_download.log")
    benchmark:
        str(RESULTS / "benchmarks" / "genomes_download.txt")
    shell:
        """
        python {params.script:q} \
            --config {params.config_file:q} \
            --genomic-dir {output.genomic:q} \
            --cds-dir {output.cds:q} \
            --rna-dir {output.rna:q} \
            --manifest {output.manifest:q} \
            --log {log:q}
        """
