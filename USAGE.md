# AmPair Usage

This document describes how to install, configure, run, and operate AmPair, and
how to read its outputs. The project overview, motivation, architecture, and
limitations live in [README.md](README.md).

## 1. Quick Start

```bash
# 1. Clone the repository
git clone git@github.com:Xinming9606/AmPair.git
cd AmPair

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

Open the HTML file in a browser to inspect the recommended primer pair and
validation summary.

## 2. Setup

See the **Setup** section in [README.md](README.md#setup). It covers Pixi
install, verifying the workflow and external tools, the micromamba/conda legacy
path, and the default vs `dev` environment split.

## 3. Running the Pipeline

The pipeline is a Snakemake workflow. The root `Snakefile` includes
`workflow/Snakefile`; run everything from the repository root.

```bash
# Dry run (show the planned steps without executing)
pixi run dry-run

# Full run with 4 cores, capped at 8 GB scheduled memory
pixi run pipeline
```

Equivalent direct Snakemake invocations:

```bash
snakemake --snakefile Snakefile -n all \
  --resources mem_mb=8192 \
  --shared-fs-usage input-output persistence software-deployment sources

snakemake --snakefile Snakefile --cores 4 all \
  --resources mem_mb=8192 \
  --shared-fs-usage input-output persistence software-deployment sources
```

For very large genera, start with a stricter `assembly_level` (e.g. `complete`),
a smaller gene set, or a lower `pcr_top_n`, and inspect the per-step logs and
benchmarks under `results/<genus>/logs/` and `results/<genus>/benchmarks/`.

## 4. Configuration

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
| `max_primer_pairs`                     | Maximum number of scored primer pairs retained by primer design.                                |

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

### Scientific-refinement options

The following options borrow ideas from the MSA-primers manuscript approach.
They are **disabled by default** so the default behaviour is unchanged; enable
them to test the improvements.

```yaml
# (1) Adaptive div_cut relaxation. Set to a positive integer N to raise div_cut
#     in steps of div_cut_auto_step until at least N candidate kmers pass, or
#     div_cut_auto_max is reached. Prevents empty output on highly diverse genera.
div_cut_auto_min_candidates: null   # e.g. 50
div_cut_auto_step: 0.05
div_cut_auto_max: 3.0

# (2) Degeneracy penalty in the combined score. 0 keeps the score without a fold
#     term. Set >0 (e.g. 0.5) to de-rank high-fold (highly degenerate) primer pairs:
#     score = 1 / (|pair_div| + 10*delta_gc^2 + w_fold*log2(total_fold) + 0.01)
score_weight_fold: 0.0

# (3) Candidate-pair score heatmap. 0 disables. Set e.g. 20 to write a Top-N pair
#     combined_score heatmap next to the diversity plot.
plot_pair_heatmap_top_n: 0
```

> Note: `in silico PCR` validation uses either degenerate matching
> (`seqkit locate --degenerate`, no mismatches) or strict matching
> (`seqkit amplicon -m`, with mismatches). The two modes can change which
> top-scoring pairs pass validation, so the final recommended pair may differ
> between modes even when the designed pairs are identical. This is controlled
> by `pcr_mismatch`/`--max-mismatch` in the PCR rule.

## 5. Reading the Report

The main deliverable is `results/<genus>/reports/<gene>_report.html`. It shows,
for one gene:

- the recommended primer pair, with forward and reverse sequences
- in silico PCR and species-level validation of the recommended pair
- a per-position Shannon entropy plot with the top primer sites marked
- all scored candidate primer pairs

When multiple genes are configured,
`results/<genus>/reports/gene_report_cross.html` compares amplification and
species-level metrics across genes.

All output files for a genus live under `results/<genus>/`:

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

## 6. Standalone Command-Line Tools

The actual work is done by standalone Python command-line tools under
`workflow/scripts/`. Each accepts `--help`:

```bash
pixi run python workflow/scripts/gene_extract.py --help
pixi run python workflow/scripts/primers_design.py --help
pixi run python workflow/scripts/gene_report.py --help
pixi run python workflow/scripts/in_silico_pcr.py --help
pixi run python workflow/scripts/design_compare.py --help
```

Example: run only the primer-design step on an existing alignment

```bash
pixi run python workflow/scripts/primers_design.py \
  --aln test_data/tuf_align.aln \
  --out-tsv results/Borrelia/primers/recG_primers.tsv \
  --out-plot results/Borrelia/primers/recG_diversity.png \
  --config config/config.yaml \
  --log results/Borrelia/logs/primers_design.log
```

Example: compare two design outputs (current vs original)

```bash
pixi run python workflow/scripts/design_compare.py \
  --current results/Borrelia/primers/recG_primers.tsv \
  --original original_run/recG_primers.tsv \
  --top-n 50
```

`design_compare.py` is dependency-free (standard library only) and exits
non-zero when the Top-1 recommended pairs differ, so it can gate CI/diff checks.

## 7. Visualization Module

The plotting routines live in a standalone module
[`workflow/scripts/primers_plot.py`](workflow/scripts/primers_plot.py) so the
primer-design logic stays free of matplotlib concerns. It is imported and called
by `primers_design.py`; you can also reuse it directly.

matplotlib and NumPy are imported lazily inside each plotting function, so no
separate setup call is needed before calling them:

```python
from primers_plot import (
    plot_diversity,
    plot_placeholder,
    plot_pair_heatmap,
)

# Text-only figure when no candidate primers pass the filters:
plot_placeholder(
    message="No conserved primer windows found under the current div_cut",
    aln_file="test_data/tuf_align.aln",
    out_plot="results/Borrelia/primers/recG_diversity.png",
    log=logging.getLogger(__name__),
)

# Per-position Shannon entropy with rolling mean and Top-N primer sites:
plot_diversity(
    divs=divs,                 # list[float] entropy per alignment column
    roll_means=roll_means,     # list[float] rolling mean (NaN where undefined)
    roll_k=roll_k,             # int rolling window size
    results=results,           # list[dict] ranked candidate primer pairs
    top_n=top_n,               # int number of top sites to shade
    primer_len=primer_len,     # int primer length in bp
    aln_file="test_data/tuf_align.aln",
    out_plot="results/Borrelia/primers/recG_diversity.png",
    log=logging.getLogger(__name__),
)

# Top-N candidate-pair combined_score heatmap (no-op when top_n <= 0):
plot_pair_heatmap(
    results=results,           # list[dict] ranked candidate primer pairs
    top_n=plot_pair_heatmap_top_n,  # int from config (0 disables)
    out_plot="results/Borrelia/primers/recG_pair_heatmap.png",
    log=logging.getLogger(__name__),
)
```

Public API:

| Function                                                                                        | Purpose                                                                                      |
| ----------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------- |
| `plot_placeholder(message, aln_file, out_plot, log)`                                            | Text-only figure used when no candidates pass.                                               |
| `plot_diversity(divs, roll_means, roll_k, results, top_n, primer_len, aln_file, out_plot, log)` | Diversity plot with Top-N sites shaded.                                                      |
| `plot_pair_heatmap(results, top_n, out_plot, log)`                                              | Combined-score heatmap for the Top-N pairs; no-op when `top_n <= 0` or fewer than 2 results. |

## 8. Requirements

The workflow dependencies are declared in [pixi.toml](pixi.toml); the legacy
[workflow/envs/environment.yaml](workflow/envs/environment.yaml) mirrors them for
micromamba/conda users.

Both files include the same default runtime: Snakemake (`snakemake-minimal`),
Python, NumPy, Matplotlib (`matplotlib-base`), `ncbi-genome-download`, PyYAML,
Python `markdown`, VSEARCH, MUSCLE, and SeqKit. The development/checking tools
Ruff and `ty` are pinned in the `dev` feature (not installed by default).

The three command-line tools are installed by platform:

| Tool    | Ubuntu / macOS              | Windows                 |
| ------- | --------------------------- | ----------------------- |
| VSEARCH | Installed by `pixi install` | Installed through Scoop |
| MUSCLE  | Installed by `pixi install` | Installed through Scoop |
| SeqKit  | Installed by `pixi install` | Installed through Scoop |

### Windows command-line tools

On Windows, install Scoop first, then run:

```powershell
scoop bucket add main-plus https://github.com/Scoopforge/Main-Plus
scoop install vsearch muscle seqkit
```

## 9. Troubleshooting

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

If genome download fails, confirm that the genus name is recognized by NCBI and
that your internet connection is available.

If the test dataset check reports a checksum mismatch, remove the archive and its
`.json` sidecar, then regenerate them together with
`pixi run download-ci-test-data`.

If a batch run is slow, inspect the per-step logs and benchmarks under
`results/<genus>/logs/` and `results/<genus>/benchmarks/`. The Python sequence
steps log input sequence counts, centroid counts, scanned genome bases, and
elapsed time. Start with a stricter assembly level such as `complete`, a smaller
gene set, or a lower `pcr_top_n` when first testing a large genus.

Check the report or `results/<genus>/aligned/<gene>.alignment.tsv` to see the
MUSCLE executable and version used for alignment.

## 10. Development Commands

Useful Pixi tasks (run with `pixi run <task>`):

```bash
pixi run compile
pixi run metadata-check
pixi run lint                 # dev environment only (ruff)
pixi run format-check         # dev environment only (ruff)
pixi run ty-check             # dev environment only (ty)
pixi run ensure-dependencies
pixi run performance-smoke
pixi run smoke
pixi run dry-run
pixi run pipeline
pixi run ci
pixi run functional-test
pixi run functional-test-ci
pixi run download-ci-test-data   # CI only: fetch and cache the reference test dataset
pixi run source-archive
pixi run conda-build
pixi run conda-install-test
```

- `source-archive` writes source `.zip` and `.tar.gz` archives under `dist/`.
- `conda-build` writes a local conda package under `dist/conda/`.
- `conda-install-test` builds `ampair`, publishes it to an indexed local conda
  channel under `dist/conda-channel/`, installs it into a fresh Pixi consumer
  project, checks the `ampair` command, verifies the bundled config/workflow
  resources, runs a Snakemake dry run, and executes the functional test against
  the repository fixture.
- `metadata-check` keeps mirrored project metadata honest: package names and
  versions must match across `pixi.toml` and `pyproject.toml`, conda runtime
  dependencies must stay in `pixi.toml`, and the legacy `environment.yaml` must
  mirror the default Pixi environment.
- `download-ci-test-data` is used by CI to fetch and cache the Borrelia reference
  dataset; local runs do not access NCBI.

### Python API

You can also call the workflow through the lightweight Python API:

```python
from ampair import AmPairProject

project = AmPairProject()
result = project.run_functional_test()
print(result.report_html)
```

After installing the conda package, the same API is exposed as the `ampair`
command:

```bash
ampair functional-test
ampair verify --genus Borrelia --gene recG
```

## 11. Project Layout

```text
AmPair/
|-- Snakefile
|-- pixi.toml
|-- pixi.lock
|-- pyproject.toml
|-- config/
|   `-- config.yaml
|-- tools/            # development and release helper scripts
|-- ampair/           # the installed Python package
|-- workflow/
|   |-- Snakefile
|   |-- rules/        # Snakemake rules
|   |-- scripts/      # standalone Python command-line tools
|   `-- envs/         # legacy environment.yaml
`-- results/          # generated outputs (not committed)
```

The GitHub Actions workflow downloads and caches the Borrelia test archive
before running `pixi run ci`; the archive is paired with a SHA-256 sidecar and an
explicit cache snapshot version. Bump that version in the workflow when
refreshing the reference dataset. Local `pixi run ci` does not access NCBI.
Pushing a tag like `v0.1.0` runs the release workflow, verifies a clean
conda-package install, builds source archives plus a conda package under `dist/`,
uploads them as workflow artifacts, and attaches them to the GitHub Release.
