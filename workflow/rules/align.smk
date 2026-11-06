# =============================================================================
# align.smk
#
# Per-gene rule: multiple sequence alignment of the clustered centroids with
# MUSCLE. The executable is installed by Pixi/Conda on Unix and by Scoop on
# Windows when it is missing.
#
# Input  : results/{genus}/extracted/{gene}.centroids.fasta  [temp, from cluster]
# Output : results/{genus}/aligned/{gene}.aln                [temp]
#          results/{genus}/aligned/{gene}.alignment.tsv      [backend metadata]
#
# .aln is temp(): consumed only by primers_design; not needed afterwards.
# =============================================================================


rule align:
    threads: 4
    resources:
        mem_mb=4096
    input:
        centroids = str(EXTRACTED / layout.centroids("{gene}")),
        config = CONFIG_FILE
    output:
        aln = temp(str(ALIGNED / layout.aln("{gene}"))),
        metadata = str(ALIGNED / layout.alignment_meta("{gene}"))
    params:
        config_file = CONFIG_FILE,
        script = str(SCRIPTS_DIR / "fasta_align.py")
    log:
        str(RESULTS / "logs" / "align" / "{gene}.log")
    benchmark:
        str(RESULTS / "benchmarks" / "align" / "{gene}.txt")
    shell:
        """
        python {params.script:q} \
            --input {input.centroids:q} \
            --output {output.aln:q} \
            --metadata {output.metadata:q} \
            --threads {threads} \
            --config {params.config_file:q} \
            --log {log:q}
        """
