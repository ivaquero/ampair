# =============================================================================
# reports.smk
#
# Report generation rules.
#
# gene_report        : per-gene HTML report (primers + PCR + diversity)
# comparison_report  : cross-gene HTML report (only when >1 gene)
#
# Reports are rendered from R Markdown templates via rmarkdown::render().
# All paths are made absolute and knit_root_dir is pinned to the project
# root, because rmarkdown::render() otherwise changes the working directory
# to the Rmd's location and breaks relative output/param paths.
# =============================================================================


rule gene_report:
    input:
        primers   = str(PRIMERS / "{gene}_primers.tsv"),
        amplicons = str(PRIMERS / "{gene}_amplicons.tsv"),
        diversity = str(PRIMERS / "{gene}_diversity.png"),
        rmd       = "workflow/scripts/gene_report.Rmd"
    output:
        str(REPORTS / "{gene}_report.html")
    params:
        gene  = lambda wc: wc.gene,
        genus = config["genus"]
    log:
        str(RESULTS / "logs" / "gene_report" / "{gene}.log")
    benchmark:
        str(RESULTS / "benchmarks" / "gene_report" / "{gene}.txt")
    shell:
        r"""
        Rscript -e '
            root <- getwd()
            rmarkdown::render(
                input         = file.path(root, "{input.rmd}"),
                output_file   = file.path(root, "{output}"),
                knit_root_dir = root,
                params = list(
                    gene          = "{params.gene}",
                    genus         = "{params.genus}",
                    primers_tsv   = file.path(root, "{input.primers}"),
                    amplicons_tsv = file.path(root, "{input.amplicons}"),
                    diversity_png = file.path(root, "{input.diversity}")
                ),
                quiet = TRUE
            )
        ' > {log} 2>&1
        """
