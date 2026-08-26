# =============================================================================
# design_primers.smk
#
# Per-gene rule: runs design_primers.R on the aligned FASTA to produce a
# ranked primer-pair TSV and a diversity PNG with primer sites overlaid.
#
# Inputs  : results/{genus}/aligned/{gene}.aln
# Outputs : results/{genus}/primers/{gene}_primers.tsv
#           results/{genus}/primers/{gene}_diversity.png
# =============================================================================


rule design_primers:
    input:
        aln = str(RESULTS / "aligned" / "{gene}.aln")
    output:
        tsv  = str(PRIMERS / "{gene}_primers.tsv"),
        plot = str(PRIMERS / "{gene}_diversity.png")
    params:
        primer_len       = config["primer_len"],
        amplicon_min_len = config["amplicon_min_len"],
        amplicon_max_len = config["amplicon_max_len"],
        div_cut          = config["div_cut"],
        GC_tol           = config["GC_tol"]
    log:
        str(RESULTS / "logs" / "design_primers" / "{gene}.log")
    benchmark:
        str(RESULTS / "benchmarks" / "design_primers" / "{gene}.txt")
    script:
        "../scripts/design_primers.R"
