# =============================================================================
# cluster.smk
#
# Per-gene rule: dereplicates near-identical sequences with VSEARCH.
# This reduces redundancy before alignment.
#
# Input  : results/{genus}/extracted/{gene}.fasta            [temp, from gene_extract]
# Output : results/{genus}/extracted/{gene}.centroids.fasta  [temp]
#
# centroids.fasta is temp(): it feeds align only. in_silico_pcr validates
# against the full genomes (genomes/genomic/), not the centroids.
# =============================================================================

rule cluster:
    threads: 4
    resources:
        mem_mb=4096
    input:
        fasta = str(EXTRACTED / "{gene}.fasta"),
        config = CONFIG_FILE
    output:
        temp(str(EXTRACTED / "{gene}.centroids.fasta"))
    params:
        identity = config["cluster_identity"],
        script = str(SCRIPTS_DIR / "fasta_cluster.py")
    log:
        str(RESULTS / "logs" / "cluster" / "{gene}.log")
    benchmark:
        str(RESULTS / "benchmarks" / "cluster" / "{gene}.txt")
    shell:
        """
        python {params.script:q} \
            --input {input.fasta:q} \
            --output {output:q} \
            --identity {params.identity} \
            --threads {threads} \
            --log {log:q}
        """
