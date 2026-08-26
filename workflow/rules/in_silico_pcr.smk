# =============================================================================
# in_silico_pcr.smk
#
# Per-gene rule: validates the top-ranked QC-passed primer candidates against
# the full genome assemblies using SeqKit. Reports are sorted
# by in silico PCR performance so the HTML report can recommend the best hit.
#
# Inputs : results/{genus}/primers/{gene}_primers.tsv   (from primers_design)
#          results/{genus}/genomes/genomic/             (from genomes_download)
# Output : results/{genus}/primers/{gene}_amplicons.tsv
# =============================================================================

rule in_silico_pcr:
    threads: 4
    resources:
        mem_mb=4096
    input:
        primers    = str(PRIMERS / "{gene}_primers.tsv"),
        genome_dir = str(GENOMES_GENOMIC),
        config = CONFIG_FILE
    output:
        summary = str(PRIMERS / "{gene}_amplicons.tsv"),
        species_summary = str(PRIMERS / "{gene}_species_summary.tsv"),
        species = str(PRIMERS / "{gene}_species.tsv")
    params:
        gene             = lambda wc: wc.gene,
        config_file = CONFIG_FILE,
        script = str(SCRIPTS_DIR / "in_silico_pcr.py")
    log:
        str(RESULTS / "logs" / "in_silico_pcr" / "{gene}.log")
    benchmark:
        str(RESULTS / "benchmarks" / "in_silico_pcr" / "{gene}.txt")
    shell:
        """
        python {params.script:q} \
            --primers-tsv {input.primers:q} \
            --genome-dir {input.genome_dir:q} \
            --out-tsv {output.summary:q} \
            --gene {params.gene:q} \
            --config {params.config_file:q} \
            --mismatch {config[pcr_mismatch]} \
            --degenerate \
            --workers {threads} \
            --species-summary {output.species_summary:q} \
            --species-tsv {output.species:q} \
            --log {log:q}
        """
