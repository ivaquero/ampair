# =============================================================================
# reports.smk
#
# Per-gene HTML reports plus a cross-gene comparison report.
# =============================================================================


rule gene_report:
    threads: 1
    resources:
        mem_mb=1024
    input:
        primers = str(PRIMERS / layout.primers("{gene}")),
        amplicons = str(PRIMERS / layout.amplicons("{gene}")),
        species_summary = str(PRIMERS / layout.species_summary("{gene}")),
        diversity = str(PRIMERS / layout.diversity("{gene}")),
        alignment_meta = str(ALIGNED / layout.alignment_meta("{gene}")),
        config = CONFIG_FILE
    output:
        str(REPORTS / layout.report_html("{gene}"))
    params:
        gene = lambda wc: wc.gene,
        config_file = CONFIG_FILE,
        download_manifest = str(RESULTS / layout.GENOMES / layout.DOWNLOAD_MANIFEST),
        script = str(SCRIPTS_DIR / "gene_report.py")
    log:
        str(RESULTS / "logs" / "gene_report" / "{gene}.log")
    benchmark:
        str(RESULTS / "benchmarks" / "gene_report" / "{gene}.txt")
    shell:
        """
        python {params.script:q} \
            --gene {params.gene:q} \
            --config {params.config_file:q} \
            --primers-tsv {input.primers:q} \
            --amplicons-tsv {input.amplicons:q} \
            --species-summary {input.species_summary:q} \
            --diversity-png {input.diversity:q} \
            --alignment-meta {input.alignment_meta:q} \
            --download-manifest {params.download_manifest:q} \
            --out-html {output:q} \
            --log {log:q}
        """


rule gene_report_cross:
    threads: 1
    resources:
        mem_mb=1024
    input:
        summaries = expand(str(PRIMERS / layout.species_summary("{gene}")), gene=GENES),
        config = CONFIG_FILE
    output:
        str(REPORTS / layout.CROSS_REPORT)
    params:
        script = str(SCRIPTS_DIR / "gene_report_cross.py")
    log:
        str(RESULTS / "logs" / "gene_report_cross.log")
    benchmark:
        str(RESULTS / "benchmarks" / "gene_report_cross.txt")
    shell:
        """
        python {params.script:q} \
            --summary-tsv {input.summaries:q} \
            --out-html {output:q}
        """
