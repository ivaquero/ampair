# =============================================================================
# extract_gene.smk
#
# Per-gene rule: extracts target gene sequences from all downloaded CDS and
# RNA FASTA files, producing a single merged FASTA per gene.
#
# Depends on: checkpoint download_genomes. The input function forces the
# checkpoint to complete (establishing the DAG dependency). The cds/rna
# directory paths are deterministic (results/{genus}/genomes/{cds,rna}) and
# are passed as plain-string params built from the global path variables, so
# the script never has to introspect checkpoint output or resolved inputs.
#
# Output : results/{genus}/extracted/{gene}.fasta   [temp]
# =============================================================================

def trigger_download(wildcards):
    # Force the checkpoint to complete; return its outputs purely to
    # establish the DAG edge. We do not rely on the return value downstream.
    co = checkpoints.download_genomes.get(**wildcards).output
    return [str(co.cds), str(co.rna)]

rule extract_gene:
    input:
        trigger_download
    output:
        temp(str(EXTRACTED / "{gene}.fasta"))
    params:
        gene    = lambda wc: wc.gene,
        aliases = lambda wc: config.get("gene_aliases", {}).get(wc.gene, []),
        cds_dir = str(GENOMES_CDS),
        rna_dir = str(GENOMES_RNA)
    log:
        str(RESULTS / "logs" / "extract_gene" / "{gene}.log")
    benchmark:
        str(RESULTS / "benchmarks" / "extract_gene" / "{gene}.txt")
    script:
        "../scripts/extract_gene.py"
