# AmPrime

[![CI](https://github.com/Xinming9606/AmPrime/actions/workflows/ci.yml/badge.svg)](https://github.com/Xinming9606/AmPrime/actions/workflows/ci.yml)
[![Release](https://github.com/Xinming9606/AmPrime/actions/workflows/release.yml/badge.svg)](https://github.com/Xinming9606/AmPrime/actions/workflows/release.yml)
[![Managed with Pixi](https://img.shields.io/badge/managed%20with-pixi-ffcb47)](https://pixi.sh)
![Python](https://img.shields.io/badge/python-3.12-blue)
![Platforms](https://img.shields.io/badge/platforms-linux--64%20%7C%20osx--arm64%20%7C%20win--64-2ea44f)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

AmPrime designs and validates amplicon-sequencing primer pairs for bacterial housekeeping genes using public NCBI genomes.

You give it:

- a bacterial genus, such as `Borrelia`
- one or more target genes, such as `recG`, `clpA`, or an MLST gene set

It returns one self-contained HTML report per gene with the recommended primer pair, in silico PCR and species-level validation, a sequence diversity plot, and all candidate primer pairs. When multiple genes are configured, it also writes a compact cross-gene comparison report.

## When To Use

Use AmPrime when you want a reproducible first-pass primer design workflow for a bacterial genus, especially when you want primers that should work across many genomes within that genus.

It is not a full specificity checker yet. The current pipeline checks whether the best QC-passed primer candidates amplify genomes inside the target genus, but it does not test off-target amplification outside the genus.

## How It Works

For each gene independently, AmPrime runs this pipeline:

```mermaid
flowchart TD
    CFG[config.yaml<br/>genus + genes]

    subgraph S1[Genome Inputs]
        direction LR
        A[Download<br/>NCBI genomes] --> B[Extract gene<br/>CDS + rRNA]
    end

    subgraph S2[Representative Sequences]
        direction LR
        C[Dereplicate] --> D[Align]
    end

    subgraph S3[Primer Selection]
        direction LR
        E[Design<br/>entropy scan] --> F[QC<br/>hairpin + dimer]
    end

    subgraph S4[Validation And Output]
        direction LR
        G[In silico<br/>PCR] --> H[HTML<br/>report]
    end

    CFG --> A
    B --> C
    D --> E
    F --> G
```

The main idea is simple:

1. Download genomic, CDS, and RNA FASTA files for the target genus.
2. Extract the target gene from CDS/RNA annotations.
3. Dereplicate near-identical sequences so redundant strains do not dominate.
4. Align representative sequences with MUSCLE.
5. Find conserved primer windows that flank a variable amplicon region.
6. Filter primer pairs for simple secondary-structure risks.
7. Validate the best QC-passed candidates against full genomes with SeqKit.
8. Write a browsable HTML report.

Missing genes are handled gracefully. If a gene cannot be found, the pipeline continues and writes a report showing that no candidates were available.

## Quick Start

```bash
# 1. Clone the repository
git clone git@github.com:Xinming9606/AmPrime.git
cd AmPrime

# 2. Install Pixi if needed: https://pixi.sh

# 3. Preview the work to be done
pixi run dry-run

# 4. Edit config/config.yaml
# Set your genus, genes, and optional gene aliases.

# 5. Run the pipeline
pixi run pipeline
```

Reports are written to:

```text
results/<genus>/reports/<gene>_report.html
```

Open the HTML file in a browser to inspect the recommended primer pair and validation summary.

## Configuration

All user-facing settings live in [config/config.yaml](config/config.yaml).

Minimal example:

```yaml
genus: Borrelia

genes:
  - recG
  - clpA
  - uvrA

assembly_level: complete
primer_len: 20
amplicon_min_len: 300
amplicon_max_len: 1000
div_cut: 2.0
GC_tol: 0.1

pcr_mismatch: 3
pcr_top_n: 10
max_primer_pairs: 100000
```

Useful options:

| Setting                                | Meaning                                                                                         |
| -------------------------------------- | ----------------------------------------------------------------------------------------------- |
| `genus`                                | Bacterial genus name recognized by NCBI.                                                        |
| `genes`                                | One or more gene names. Each gene is processed independently.                                   |
| `assembly_level`                       | NCBI assembly level: `complete`, `chromosome`, `scaffold`, or `contig`.                         |
| `primer_len`                           | Primer length in bp.                                                                            |
| `amplicon_min_len`, `amplicon_max_len` | Target amplicon size range.                                                                     |
| `div_cut`                              | Maximum Shannon entropy allowed for conserved primer windows. Raise it if no primers are found. |
| `GC_tol`                               | Maximum GC fraction difference between forward and reverse primers.                             |
| `pcr_mismatch`                         | Mismatches allowed per primer during in silico PCR.                                             |
| `pcr_top_n`                            | Number of QC-passed primer pairs to validate before choosing the report recommendation.         |
| `max_primer_pairs`                     | Maximum number of scored primer pairs retained by primer design.                               |

If a gene is annotated under different names across genomes, add aliases:

```yaml
gene_aliases:
  tuf:
    - tsf
  16S:
    - "16S ribosomal RNA"
```

You can also override the diversity cutoff for specific genes:

```yaml
div_cut_per_gene:
  16S: 3.0
  rpoB: 1.5
```

## Output Files

For each gene, outputs are written under `results/<genus>/`.

| File                                 | Contents                                                                             |
| ------------------------------------ | ------------------------------------------------------------------------------------ |
| `reports/<gene>_report.html`         | Main deliverable: recommendation, PCR/species validation, plot, and candidates.      |
| `reports/gene_report_cross.html`     | Cross-gene comparison of amplification and species-level metrics.                    |
| `genomes/download_manifest.tsv`      | Download manifest with FASTA counts, sizes, config SHA-256, and data fingerprints.   |
| `aligned/<gene>.alignment.tsv`       | MUSCLE version and alignment metadata.                                               |
| `primers/<gene>_primers.tsv`         | Filtered candidate primer pairs ranked by score.                                     |
| `primers/<gene>_amplicons.tsv`       | In silico PCR results for the top validated primer candidates, sorted best first.    |
| `primers/<gene>_species_summary.tsv` | Species-level amplification, allele multiplicity, and inter-species overlap metrics. |
| `primers/<gene>_species.tsv`         | Per-species amplicon allele and overlap summary.                                     |
| `primers/<gene>_diversity.png`       | Per-position Shannon entropy plot with top primer sites marked.                      |
| `logs/...`                           | Per-step logs for debugging.                                                         |
| `benchmarks/...`                     | Snakemake benchmark files.                                                           |

## Project Layout

```text
AmPrime/
|-- Snakefile
|-- pixi.toml
|-- pixi.lock
|-- pyproject.toml
|-- .github/
|   `-- workflows/
|-- config/
|   `-- config.yaml
|-- tools/
|   |-- compile_project.py
|   |-- check_metadata.py
|   |-- download_test_dataset.py
|   |-- verify_test_dataset.py
|   |-- smoke_project.py
|   |-- test_built_package.py
|   `-- package_release.py
|-- amprime/
|   |-- api.py
|   |-- cli.py
|   `-- __main__.py
|-- docs/
|   |-- paper.md
|   `-- test-functional.md
|-- workflow/
|   |-- Snakefile
|   |-- rules/
|   |   |-- genomes_download.smk
|   |   |-- gene_extract.smk
|   |   |-- cluster.smk
|   |   |-- align.smk
|   |   |-- primers_design.smk
|   |   |-- primers_check.smk
|   |   |-- in_silico_pcr.smk
|   |   `-- reports.smk
|   |-- scripts/
|   |   |-- genomes_download.py
|   |   |-- gene_extract.py
|   |   |-- fasta_cluster.py
|   |   |-- fasta_align.py
|   |   |-- primers_design.py
|   |   |-- primers_check.py
|   |   |-- in_silico_pcr.py
|   |   |-- gene_report.py
|   |   |-- gene_report_cross.py
|   |   |-- gene_report_cross.html
|   |   |-- report.css
|   |   |-- dependencies.py
|   |   |-- common.py
|   |   |-- fasta_io.py
|   |   `-- gene_report.html
|   `-- envs/
|       `-- environment.yaml
`-- results/
```

## Design Notes

Snakemake is intentionally kept as a thin scheduler. It manages dependencies, parallel execution, logs, benchmarks, and resumability. The actual work is done by standalone Python command-line tools in `workflow/scripts/`. The shared `dependencies.py` module detects VSEARCH, MUSCLE, and SeqKit and installs missing Windows tools through Scoop.

Download outputs are refreshed as a unit when the download rule runs, so stale FASTA files from a previous genus or assembly level do not mix into a new run. Clustering always uses VSEARCH at 97% identity, and representative sequences are aligned with MUSCLE. Alignment runs write a small metadata TSV next to the alignment, recording the MUSCLE executable and version. The same alignment summary is included in each HTML report so runs remain easy to audit.

Gene extraction is a single batch scan for all configured genes, avoiding a full CDS/RNA directory rescan per gene. In-silico PCR batches genomes into SeqKit `locate` invocations, searches both strands for every candidate site, reconstructs every valid primer-site pairing, and aggregates deterministic report metrics in Python. Workflow rules declare thread and memory requirements; the Pixi convenience tasks cap scheduled memory at 8 GB. The `performance-smoke` task guards against large regressions before very large genera are processed.

This keeps each step easy to test and debug. For example:

```bash
python workflow/scripts/gene_extract.py --help
python workflow/scripts/primers_design.py --help
python workflow/scripts/gene_report.py --help
```

## Development

Pixi is the primary project manager. It creates the conda/bioconda environment from [pixi.toml](pixi.toml) and keeps runs reproducible with `pixi.lock`. The locked Pixi platforms are Linux, Apple Silicon macOS, and Windows. Intel macOS is not supported. Pixi installs VSEARCH, MUSCLE, and SeqKit automatically on Ubuntu and Apple Silicon macOS. On Windows, the workflow checks for all three executables and installs missing tools with Scoop when Scoop is available.

Useful commands:

```bash
pixi run compile
pixi run metadata-check
pixi run lint
pixi run format-check
pixi run ty-check
pixi run ensure-dependencies
pixi run performance-smoke
pixi run smoke
pixi run dry-run
pixi run pipeline
pixi run ci
pixi run functional-test
pixi run functional-test-ci
pixi run source-archive
pixi run conda-build
pixi run conda-install-test
```

- `source-archive` writes source `.zip` and `.tar.gz` archives under `dist/`.
- `conda-build` writes a local conda package under `dist/conda/`.
- `conda-install-test` builds `amprime`, publishes it to an indexed local conda channel under `dist/conda-channel/`, installs it into a fresh Pixi consumer project, checks the `amprime` command, verifies the bundled config/workflow resources, runs a Snakemake dry run, and executes the functional test against the repository fixture.
- `metadata-check` keeps mirrored project metadata honest: package names and versions must match across `pixi.toml` and `pyproject.toml`, conda runtime dependencies must stay in `pixi.toml`, and the legacy `environment.yaml` must mirror the default Pixi environment.

For the end-to-end Borrelia functional test, including CI dataset download, see
[docs/test-functional.md](docs/test-functional.md).

You can also call the workflow through the lightweight Python API:

```python
from amprime import AmPrimeProject

project = AmPrimeProject()
result = project.run_functional_test()
print(result.report_html)
```

After installing the conda package, the same API is exposed as the `amprime` command:

```bash
amprime functional-test
amprime verify --genus Borrelia --gene recG
```

The GitHub Actions workflow downloads and caches the Borrelia test archive before running `pixi run ci`; the archive is paired with a SHA-256 sidecar and an explicit cache snapshot version. Bump that version in the workflow when refreshing the reference dataset. Local `pixi run ci` does not access NCBI. Pushing a tag like `v0.1.0` runs the release workflow, verifies a clean conda-package install, builds source archives plus a conda package under `dist/`, uploads them as workflow artifacts, and attaches them to the GitHub Release.

## Troubleshooting

If the dry run fails, check that you are running from the repository root:

```bash
snakemake -n
```

If no primers are found, try one or more of the following:

- raise `div_cut`
- relax `GC_tol`
- add gene aliases under `gene_aliases`
- use a less restrictive `assembly_level`
- inspect the per-gene logs under `results/<genus>/logs/`

If genome download fails, confirm that the genus name is recognized by NCBI and that your internet connection is available.

On Ubuntu and Apple Silicon macOS, `pixi install` supplies VSEARCH and MUSCLE.
Check them with:

```bash
vsearch --version
muscle --version
seqkit version
```

On Windows, VSEARCH, MUSCLE, and SeqKit are required. Dependency setup first checks for each executable, then checks Scoop, loads the `main-plus` bucket with `scoop bucket add main-plus https://github.com/Scoopforge/Main-Plus` when needed, and installs the missing tool with Scoop. The Windows CI and Release workflows perform the same bucket setup explicitly before running the dependency check. If Scoop is not installed, install it first using the instructions below.

If the test dataset check reports a checksum mismatch, remove the archive and its `.json` sidecar, then regenerate them together with `pixi run download-ci-test-data`.

If a batch run is slow, inspect the per-step logs and benchmarks under `results/<genus>/logs/` and `results/<genus>/benchmarks/`. The Python sequence steps log input sequence counts, centroid counts, scanned genome bases, and elapsed time. Start with a stricter assembly level such as `complete`, a smaller gene set, or a lower `pcr_top_n` when first testing a large genus.

Check the report or `results/<genus>/aligned/<gene>.alignment.tsv` to see the MUSCLE executable and version used for alignment.

## Requirements

Use Pixi for the simplest setup:

```bash
pixi run ci
```

The workflow dependencies are declared in:

```text
pixi.toml
```

The legacy micromamba/conda environment file mirrors the default Pixi dependencies for users who prefer that tooling:

```bash
micromamba env create -f workflow/envs/environment.yaml
micromamba activate amprime
snakemake --cores 4
```

That environment file lives at:

```text
workflow/envs/environment.yaml
```

Both dependency files include the same default runtime: Snakemake, Python 3.12, Biopython, NumPy, Matplotlib, `ncbi-genome-download`, PyYAML, Python `markdown`, VSEARCH, MUSCLE, and SeqKit. Pixi installs the command-line tools on Linux and Apple Silicon macOS; Windows installs them through Scoop when the workflow runs.

### Windows command-line tools

On Windows, AmPrime requires VSEARCH, MUSCLE, and SeqKit. Install Scoop first, then run:

```powershell
scoop bucket add main-plus https://github.com/Scoopforge/Main-Plus
scoop install vsearch muscle seqkit
```

## Limitations

- Off-target specificity outside the target genus is not checked yet.
- Primer sequences use a majority-rule consensus and IUPAC representation;
  Shannon entropy is calculated independently from the observed A/C/G/T bases
  in each alignment column.
- MUSCLE is required for multiple sequence alignment.
- Degenerate-base handling is conservative.
- Very large genera can take a long time to download, align, and scan even with batched SeqKit validation.

## License

MIT. See [LICENSE](LICENSE).
