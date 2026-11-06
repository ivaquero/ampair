# =============================================================================
# primers_check.smk
#
# Per-gene rule: computes hairpin, dimer, and 3'-end stability metrics for
# every primer pair, filters out pairs that fail configurable thresholds,
# and writes the result as the final primer TSV consumed by downstream rules.
#
# Input  : results/{genus}/primers/{gene}_primers_raw.tsv   [temp, from primers_design]
# Output : results/{genus}/primers/{gene}_primers.tsv       [final]
# =============================================================================


rule primers_check:
    threads: 1
    resources:
        mem_mb=2048
    input:
        primer_tsv = str(PRIMERS / layout.primers_raw("{gene}")),
        config = CONFIG_FILE
    output:
        str(PRIMERS / layout.primers("{gene}"))
    params:
        config_file = CONFIG_FILE,
        script = str(SCRIPTS_DIR / "primers_check.py")
    log:
        str(RESULTS / "logs" / "primers_check" / "{gene}.log")
    benchmark:
        str(RESULTS / "benchmarks" / "primers_check" / "{gene}.txt")
    shell:
        """
        python {params.script:q} \
            --in-tsv {input.primer_tsv:q} \
            --out-tsv {output:q} \
            --config {params.config_file:q} \
            --log {log:q}
        """
