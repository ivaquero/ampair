# =============================================================================
# in_silico_pcr.smk
#
# Per-gene rule: validates the top-ranked primer pair against the full
# genome assemblies using seqkit amplicon. Reports amplification rate
# and mean amplicon length.
#
# Inputs : results/{genus}/primers/{gene}_primers.tsv   (from design_primers)
#          results/{genus}/genomes/genomic/             (from download_genomes)
# Output : results/{genus}/primers/{gene}_amplicons.tsv
# =============================================================================

rule in_silico_pcr:
    input:
        primers    = str(PRIMERS / "{gene}_primers.tsv"),
        genome_dir = str(GENOMES_GENOMIC)
    output:
        str(PRIMERS / "{gene}_amplicons.tsv")
    params:
        gene             = lambda wc: wc.gene,
        mismatch         = config["pcr_mismatch"],
        amplicon_min_len = config["amplicon_min_len"],
        amplicon_max_len = config["amplicon_max_len"]
    log:
        str(RESULTS / "logs" / "in_silico_pcr" / "{gene}.log")
    benchmark:
        str(RESULTS / "benchmarks" / "in_silico_pcr" / "{gene}.txt")
    script:
        "../scripts/in_silico_pcr.py"
