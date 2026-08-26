# =============================================================================
# download_genomes.smk
#
# Downloads three FASTA types for the target genus from NCBI, each into its
# own subdirectory:
#   genomic   — full genome assemblies, used by in_silico_pcr (KEPT)
#   cds       — CDS FASTA, used to extract protein-coding genes
#   rna       — RNA FASTA, used to extract rRNA genes
#
# Uses ncbi-genome-download with --genera (no taxid resolution needed).
#
# Outputs (checkpoint):
#   results/{genus}/genomes/genomic/   — *_genomic.fna           [KEPT]
#   results/{genus}/genomes/cds/       — *_cds_from_genomic.fna
#   results/{genus}/genomes/rna/       — *_rna_from_genomic.fna
# =============================================================================


checkpoint download_genomes:
    output:
        genomic = directory(str(GENOMES_GENOMIC)),
        cds     = directory(str(GENOMES_CDS)),
        rna     = directory(str(GENOMES_RNA))
    params:
        genus          = config["genus"],
        assembly_level = config["assembly_level"]
    log:
        str(RESULTS / "logs" / "download_genomes.log")
    benchmark:
        str(RESULTS / "benchmarks" / "download_genomes.txt")
    shell:
        """
        mkdir -p {output.genomic} {output.cds} {output.rna}

        # Genomic FASTA (full assemblies, for in silico PCR)
        ncbi-genome-download bacteria \
            --genera "{params.genus}" \
            --assembly-levels {params.assembly_level} \
            --formats fasta \
            --flat-output \
            --output-folder {output.genomic} \
            2>> {log}

        # CDS FASTA (for protein-coding gene extraction)
        ncbi-genome-download bacteria \
            --genera "{params.genus}" \
            --assembly-levels {params.assembly_level} \
            --formats cds-fasta \
            --flat-output \
            --output-folder {output.cds} \
            2>> {log}

        # RNA FASTA (for rRNA gene extraction)
        ncbi-genome-download bacteria \
            --genera "{params.genus}" \
            --assembly-levels {params.assembly_level} \
            --formats rna-fasta \
            --flat-output \
            --output-folder {output.rna} \
            2>> {log}

        # Decompress everything
        find {output.genomic} {output.cds} {output.rna} \
            -name "*.gz" -exec gunzip {{}} \\; 2>> {log}

        # Sanity check
        n_gen=$(find {output.genomic} -name "*.fna" | wc -l)
        n_cds=$(find {output.cds} -name "*.fna" | wc -l)
        n_rna=$(find {output.rna} -name "*.fna" | wc -l)
        echo "Downloaded: $n_gen genomic, $n_cds CDS, $n_rna RNA files for genus {params.genus}" >> {log}

        if [ "$n_gen" -eq 0 ]; then
            echo "ERROR: No genomic FASTA downloaded. Check genus name and assembly level." >> {log}
            exit 1
        fi
        """
