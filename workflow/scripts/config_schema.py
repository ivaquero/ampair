#!/usr/bin/env python3
"""Configuration validation shared by the Snakefile and developer checks."""

import math
import warnings

REQUIRED_SETTINGS = [
    "genus",
    "genes",
    "assembly_level",
    "primer_len",
    "amplicon_min_len",
    "amplicon_max_len",
    "div_cut",
    "GC_tol",
    "pcr_mismatch",
]

ALLOWED_ASSEMBLY_LEVELS = {"complete", "chromosome", "scaffold", "contig"}
# Single source of truth for the primer quality-check knobs. Each entry is
# (config key, default threshold kcal/mol, primers_check metric columns gated by
# that threshold). config_schema validation, config.yaml, and primers_check
# (argparse flags, threshold dict, per-column checks, output columns) all derive
# from this table so adding or renaming a dG filter is a one-line change here.
QC_THRESHOLD_SPECS = [
    ("max_hairpin_dg", 0.0, ("hairpin_fwd_dg", "hairpin_rev_dg")),
    ("max_homodimer_dg", -6.0, ("homodimer_fwd_dg", "homodimer_rev_dg")),
    ("max_heterodimer_dg", -6.0, ("heterodimer_dg",)),
    ("max_3end_dg", -8.0, ("end3_fwd_dg", "end3_rev_dg")),
]
QC_THRESHOLD_SETTINGS = [key for key, _default, _cols in QC_THRESHOLD_SPECS]
# Metric column -> owning threshold key, derived from the specs above and used
# by primers_check to gate each computed dG against its configured limit.
QC_METRIC_TO_THRESHOLD = {
    col: key for key, _default, cols in QC_THRESHOLD_SPECS for col in cols
}
# Plausible magnitude range (kcal/mol) for the nearest-neighbour dG thresholds.
# Real primer secondary-structure dG values stay well inside this window, and
# the documented "-99 to disable filtering" convention is preserved; anything
# outside it is treated as a typo.
QC_DG_BOUNDS = (-100.0, 100.0)

DEFAULT_SETTINGS = {
    "gene_aliases": {},
    "div_cut_per_gene": {},
    "min_allele_freq": 0.05,
    "max_degeneracy": 16,
    "max_primer_pairs": 100_000,
    # VSEARCH dereplication identity threshold for the cluster step. Sourced
    # here so config.yaml is the single source of truth (cluster.smk reads it
    # from the validated config instead of a hardcoded literal).
    "cluster_identity": 0.97,
    # Number of QC-passed primer pairs to validate against full genomes. The
    # default matches the shipped config.yaml reference value so a config that
    # omits this key behaves identically to the shipped one (the first-version
    # pipeline validated the top 50 candidates before recommending).
    "pcr_top_n": 50,
    # Amplicon length tolerance (bp) added to amplicon_min_len/amplicon_max_len
    # when accepting a SeqKit in silico PCR hit. Promoted from the previously
    # hardcoded ``len_margin = 100`` in in_silico_pcr.py.
    "pcr_len_margin": 100,
    # Maximum ranked candidate pairs rendered into the per-gene HTML report's
    # candidate table. The complete ranked list always remains in the primers
    # TSV; this only bounds the in-HTML table so a large max_primer_pairs run
    # does not produce an enormous page.
    "report_top_n": 100,
    # Auto-relax div_cut until at least `div_cut_min_candidates` pass, then stop.
    # Disabled by default (None) so behavior is backward-compatible with a fixed
    # div_cut. Set to a number (e.g. 50) to enable adaptive relaxation of the
    # diversity cutoff.
    "div_cut_auto_min_candidates": None,
    # Increment step and upper bound for the adaptive div_cut search.
    "div_cut_auto_step": 0.05,
    "div_cut_auto_max": 3.0,
    # Degeneracy penalty weight in the combined score. 0 (default) keeps the
    # score without a fold term; >0 penalizes high-fold primer pairs so that
    # conserved (low-degeneracy) designs rank higher.
    "score_weight_fold": 0.0,
    # Draw a Top-N candidate-pair score heatmap alongside the diversity plot.
    "plot_pair_heatmap_top_n": 0,
}
# QC dG threshold defaults are sourced from the single QC_THRESHOLD_SPECS table.
DEFAULT_SETTINGS.update({key: default for key, default, _cols in QC_THRESHOLD_SPECS})


def _is_int(value):
    return isinstance(value, int) and not isinstance(value, bool)


def _is_number(value):
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(value)
    )


def apply_config_defaults(cfg):
    """Return a shallow copy with optional settings filled in."""
    normalized = dict(DEFAULT_SETTINGS)
    normalized.update(cfg or {})
    return normalized


def load_config_file(path):
    """Load a YAML config file without making every CLI depend on PyYAML at import."""
    import yaml

    with open(path, encoding="utf-8") as fh:
        cfg = apply_config_defaults(yaml.safe_load(fh) or {})

    validate_config(cfg)
    return cfg


def validate_config(cfg):
    """Raise ValueError with actionable messages when the config is invalid."""
    errors = []

    errors.extend(
        [
            f"missing required setting: {key}"
            for key in REQUIRED_SETTINGS
            if key not in cfg
        ]
    )

    genus = cfg.get("genus")
    if not isinstance(genus, str) or not genus.strip():
        errors.append("genus must be a non-empty string")
    elif "/" in genus or "\\" in genus or "." in genus:
        errors.append(f"genus names cannot contain '/', '\\', or '.': {genus}")

    genes = cfg.get("genes")
    if not isinstance(genes, list) or not genes:
        errors.append("genes must be a non-empty list")
    else:
        for gene in genes:
            if not isinstance(gene, str) or not gene.strip():
                errors.append("each item in genes must be a non-empty string")
            elif "/" in gene or "\\" in gene or "." in gene:
                errors.append(f"gene names cannot contain '/', '\\', or '.': {gene}")
        # Duplicated genes expand to the same Snakemake output twice, which
        # surfaces later as a cryptic "duplicate output file" error.
        seen = set()
        duplicates = set()
        for gene in genes:
            if isinstance(gene, str) and gene in seen:
                duplicates.add(gene)
            seen.add(gene)
        if duplicates:
            errors.append(
                "genes must not contain duplicates: " + ", ".join(sorted(duplicates))
            )

    if cfg.get("assembly_level") not in ALLOWED_ASSEMBLY_LEVELS:
        errors.append(
            "assembly_level must be one of: complete, chromosome, scaffold, contig"
        )

    errors.extend(
        [
            f"{key} must be a positive integer"
            for key in ["primer_len", "amplicon_min_len", "amplicon_max_len"]
            if not _is_int(cfg.get(key)) or cfg.get(key) <= 0
        ]
    )

    if (
        _is_int(cfg.get("amplicon_min_len"))
        and _is_int(cfg.get("amplicon_max_len"))
        and cfg["amplicon_min_len"] > cfg["amplicon_max_len"]
    ):
        errors.append("amplicon_min_len must be <= amplicon_max_len")

    if not _is_number(cfg.get("div_cut")) or cfg.get("div_cut") < 0:
        errors.append("div_cut must be a non-negative number")

    if not _is_number(cfg.get("GC_tol")) or not 0 <= cfg.get("GC_tol") <= 1:
        errors.append("GC_tol must be a number between 0 and 1")

    min_allele_freq = cfg.get("min_allele_freq", 0.05)
    if not _is_number(min_allele_freq) or not 0 <= min_allele_freq <= 1:
        errors.append("min_allele_freq must be a number between 0 and 1")

    max_degeneracy = cfg.get("max_degeneracy", 16)
    if not _is_int(max_degeneracy) or max_degeneracy < 1:
        errors.append("max_degeneracy must be a positive integer")

    max_primer_pairs = cfg.get("max_primer_pairs", 100_000)
    if not _is_int(max_primer_pairs) or max_primer_pairs < 1:
        errors.append("max_primer_pairs must be a positive integer")

    if not _is_int(cfg.get("pcr_mismatch")) or cfg.get("pcr_mismatch") < 0:
        errors.append("pcr_mismatch must be a non-negative integer")

    pcr_top_n = cfg.get("pcr_top_n", 50)
    if not _is_int(pcr_top_n) or pcr_top_n < 1:
        errors.append("pcr_top_n must be a positive integer")

    cluster_identity = cfg.get("cluster_identity", 0.97)
    if not _is_number(cluster_identity) or not 0 < cluster_identity <= 1:
        errors.append("cluster_identity must be a number greater than 0 and at most 1")

    pcr_len_margin = cfg.get("pcr_len_margin", 100)
    if not _is_int(pcr_len_margin) or pcr_len_margin < 0:
        errors.append("pcr_len_margin must be a non-negative integer")

    report_top_n = cfg.get("report_top_n", 100)
    if not _is_int(report_top_n) or report_top_n < 1:
        errors.append("report_top_n must be a positive integer")

    # --- Scientific refinement settings ------------------------------------
    div_cut_auto_min = cfg.get("div_cut_auto_min_candidates", None)
    if div_cut_auto_min is not None and (
        not _is_int(div_cut_auto_min) or div_cut_auto_min < 1
    ):
        errors.append("div_cut_auto_min_candidates must be a positive integer or null")

    div_cut_auto_step = cfg.get("div_cut_auto_step", 0.05)
    if not _is_number(div_cut_auto_step) or div_cut_auto_step <= 0:
        errors.append("div_cut_auto_step must be a positive number")

    div_cut_auto_max = cfg.get("div_cut_auto_max", 3.0)
    if not _is_number(div_cut_auto_max) or div_cut_auto_max < 0:
        errors.append("div_cut_auto_max must be a non-negative number")

    score_weight_fold = cfg.get("score_weight_fold", 0.0)
    if not _is_number(score_weight_fold) or score_weight_fold < 0:
        errors.append("score_weight_fold must be a non-negative number")

    plot_pair_heatmap_top_n = cfg.get("plot_pair_heatmap_top_n", 0)
    if not _is_int(plot_pair_heatmap_top_n) or plot_pair_heatmap_top_n < 0:
        errors.append("plot_pair_heatmap_top_n must be a non-negative integer")

    aliases = cfg.get("gene_aliases", {})
    if not isinstance(aliases, dict):
        errors.append("gene_aliases must be a mapping")
    else:
        for gene, values in aliases.items():
            if not isinstance(gene, str) or not gene.strip():
                errors.append("gene_aliases keys must be non-empty strings")
            if not isinstance(values, list) or not all(
                isinstance(alias, str) and alias.strip() for alias in values
            ):
                errors.append(f"gene_aliases.{gene} must be a list of strings")

    div_overrides = cfg.get("div_cut_per_gene", {})
    if not isinstance(div_overrides, dict):
        errors.append("div_cut_per_gene must be a mapping")
    else:
        for gene, value in div_overrides.items():
            if not isinstance(gene, str) or not gene.strip():
                errors.append("div_cut_per_gene keys must be non-empty strings")
            if not _is_number(value) or value < 0:
                errors.append(f"div_cut_per_gene.{gene} must be non-negative")

    errors.extend(
        [
            f"{key} must be a finite number within {QC_DG_BOUNDS} kcal/mol or null"
            for key in QC_THRESHOLD_SETTINGS
            if key in cfg
            and cfg[key] is not None
            and (
                not _is_number(cfg[key])
                or not QC_DG_BOUNDS[0] <= cfg[key] <= QC_DG_BOUNDS[1]
            )
        ]
    )

    # Non-fatal: flag keys that are neither required nor known defaults so
    # typos and stray settings surface early, while still allowing forward-
    # compatible extras added in future versions.
    known_keys = set(REQUIRED_SETTINGS) | set(DEFAULT_SETTINGS)
    unknown_keys = sorted(set(cfg) - known_keys)
    if unknown_keys:
        warnings.warn(
            "config contains unrecognized key(s): " + ", ".join(unknown_keys),
            stacklevel=2,
        )

    if errors:
        raise ValueError("Invalid config/config.yaml:\n- " + "\n- ".join(errors))
