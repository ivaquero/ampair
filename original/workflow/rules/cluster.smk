# =============================================================================
# cluster.smk
#
# Per-gene rule: dereplicates near-identical sequences using vsearch
# cluster_fast at 97% identity, keeping one centroid per cluster.
# This reduces redundancy before alignment and speeds up MUSCLE.
#
# Input  : results/{genus}/extracted/{gene}.fasta            [temp, from extract_gene]
# Output : results/{genus}/extracted/{gene}.centroids.fasta  [temp]
#
# centroids.fasta is temp(): it feeds align only. in_silico_pcr validates
# against the full genomes (genomes/genomic/), not the centroids.
# =============================================================================

rule cluster:
    input:
        str(EXTRACTED / "{gene}.fasta")
    output:
        temp(str(EXTRACTED / "{gene}.centroids.fasta"))
    params:
        identity = 0.97
    log:
        str(RESULTS / "logs" / "cluster" / "{gene}.log")
    benchmark:
        str(RESULTS / "benchmarks" / "cluster" / "{gene}.txt")
    shell:
        """
        vsearch --cluster_fast {input} \
            --strand both \
            --id {params.identity} \
            --centroids {output} \
            2>> {log}

        n_in=$(grep -c "^>" {input} || echo 0)
        n_out=$(grep -c "^>" {output} || echo 0)
        echo "Clustered $n_in sequences into $n_out centroids" >> {log}
        """
