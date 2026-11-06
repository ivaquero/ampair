# =============================================================================
# primers_design.smk
#
# Per-gene rule: runs primers_design.py on the aligned FASTA to produce a
# ranked primer-pair TSV and a diversity PNG with primer sites overlaid.
#
# Inputs  : results/{genus}/aligned/{gene}.aln
# Outputs : results/{genus}/primers/{gene}_primers_raw.tsv  [temp, feeds primers_check]
#           results/{genus}/primers/{gene}_diversity.png
# =============================================================================


rule primers_design:
    threads: 1
    resources:
        mem_mb=4096
    input:
        aln = str(ALIGNED / layout.aln("{gene}")),
        config = CONFIG_FILE
    output:
        tsv  = temp(str(PRIMERS / layout.primers_raw("{gene}"))),
        plot = str(PRIMERS / layout.diversity("{gene}"))
    params:
        gene = lambda wc: wc.gene,
        config_file = CONFIG_FILE,
        script = str(SCRIPTS_DIR / "primers_design.py")
    log:
        str(RESULTS / "logs" / "primers_design" / "{gene}.log")
    benchmark:
        str(RESULTS / "benchmarks" / "primers_design" / "{gene}.txt")
    shell:
        """
        python {params.script:q} \
            --aln {input.aln:q} \
            --out-tsv {output.tsv:q} \
            --out-plot {output.plot:q} \
            --gene {params.gene:q} \
            --config {params.config_file:q} \
            --log {log:q}
        """
