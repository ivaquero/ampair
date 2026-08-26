#!/usr/bin/env python3
# =============================================================================
# primers_design.py - pure Python primer design implementation
#
# Reads a multiple sequence alignment (FASTA), computes per-position Shannon
# entropy, builds degenerate (IUPAC) consensus kmer windows, evaluates all
# valid primer pairs, scores them, and writes a ranked TSV + diversity PNG.
#
# Dependencies: numpy, matplotlib.
#
# CLI:
#   python primers_design.py --aln input.aln --out-tsv primers.tsv \
#       --out-plot diversity.png --primer-len 20 ...
#
# TSV columns (same as R version):
#   primer_id, fwd, rev, fwd_pos, rev_pos, amplicon_len,
#   fwd_GC, rev_GC, fwd_ndegen, rev_ndegen, fwd_fold, rev_fold, total_fold,
#   pair_diversity, delta_GC, combined_score
# =============================================================================

from __future__ import annotations

import argparse
import csv
import logging
import os
from bisect import bisect_left, bisect_right
from time import perf_counter
from typing import Any

from common import (
    config_param as _param,
    configure_logging,
    required_param as _required_param,
    reverse_complement as _rev_comp,
)
from config_schema import load_config_file
from fasta_io import parse_fasta
from primers_plot import (
    ensure_matplotlib,
    plot_diversity,
    plot_pair_heatmap,
    plot_placeholder,
)

# numpy/matplotlib are imported once, inside main(), and bound to these module
# globals so the computation helpers below can use np without threading it
# through every call. primers_plot keeps its own (lazy) copies for drawing.
np: Any = None
plt: Any = None

_IUPAC_MAP = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "AG": "R",
    "CT": "Y",
    "CG": "S",
    "AT": "W",
    "GT": "K",
    "AC": "M",
    "CGT": "B",
    "AGT": "D",
    "ACT": "H",
    "ACG": "V",
    "ACGT": "N",
}

_TSV_FIELDNAMES = [
    "primer_id",
    "fwd",
    "rev",
    "fwd_pos",
    "rev_pos",
    "amplicon_len",
    "fwd_GC",
    "rev_GC",
    "fwd_ndegen",
    "rev_ndegen",
    "fwd_fold",
    "rev_fold",
    "total_fold",
    "pair_diversity",
    "delta_GC",
    "combined_score",
]

_VALID_BASES = frozenset("ACGT")
WARN_KMER_COUNT = 50_000
WARN_PAIR_COUNT = 200_000


def _iupac_encode(bases: str) -> str:
    """Sorted unique base letters to IUPAC ambiguity code."""
    key = "".join(sorted(set(bases.upper()) & _VALID_BASES))
    return _IUPAC_MAP.get(key, "N")


def _shannon(col) -> float:
    """Shannon entropy of the observed A/C/G/T bases in one column.

    Gaps and unknown bases are missing observations, not additional bases in
    the distribution.
    """
    observed = col[np.isin(col, ["a", "c", "g", "t"])]
    if len(observed) == 0:
        return 0.0
    _, counts = np.unique(observed, return_counts=True)
    probs = counts / len(observed)
    probs = probs[probs > 0]
    return float(-np.sum(probs * np.log(probs)))


def _rollmean(values, k: int):
    """Rolling mean with NaN padding at edges (matches R's zoo::rollmean)."""
    result = np.full(len(values), np.nan)
    if k <= len(values):
        conv = np.convolve(values, np.ones(k) / k, mode="valid")
        start = k // 2
        result[start : start + len(conv)] = conv
    return result


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------
def _compute_column_stats(dna_matrix: np.ndarray, n_seqs: int, min_allele_freq: float):
    """Return (consensus, divs, pos_code, pos_fold) in a single column pass.

    ``consensus`` is the most frequent base (gaps and Ns included), ``divs``
    is the Shannon entropy over observed A/C/G/T, and ``pos_code``/``pos_fold``
    are the IUPAC degeneracy code and fold for each column.
    """
    aln_len = dna_matrix.shape[1]
    consensus = np.empty(aln_len, dtype="<U1")
    divs = np.empty(aln_len, dtype=float)
    pos_code = np.empty(aln_len, dtype="<U1")
    pos_fold = np.ones(aln_len, dtype=int)

    for i in range(aln_len):
        unique, counts = np.unique(dna_matrix[:, i], return_counts=True)

        # consensus: most frequent base, gaps/Ns included.
        consensus[i] = unique[np.argmax(counts)]

        # entropy: observed A/C/G/T only.
        is_acgt = np.isin(unique, ["a", "c", "g", "t"])
        acgt_counts = counts[is_acgt]
        if acgt_counts.size:
            probs = acgt_counts / acgt_counts.sum()
            probs = probs[probs > 0]
            divs[i] = float(-np.sum(probs * np.log(probs)))
        else:
            divs[i] = 0.0

        # annotations: drop gaps and Ns, apply min_allele_freq threshold.
        mask = ~np.isin(unique, ["-", "n"])
        bases = unique[mask]
        base_counts = counts[mask]
        if bases.size == 0:  # all-gap column
            pos_code[i] = "N"
            pos_fold[i] = 1
            continue

        freqs = base_counts / n_seqs
        keep = bases[freqs >= min_allele_freq]
        keep = np.array([b for b in keep if b.upper() in _VALID_BASES])

        if keep.size == 0:  # nothing passes freq threshold
            keep = np.array([bases[np.argmax(base_counts)]])

        pos_code[i] = _iupac_encode("".join(keep))
        pos_fold[i] = len(keep)

    return consensus, divs, pos_code, pos_fold


def _build_kmers(consensus, pos_code, pos_fold, divs, primer_len):
    """Slide a window of length *primer_len* across the alignment."""
    n_kmers = len(divs) - primer_len + 1
    kmers = []
    for j in range(n_kmers):
        idx = slice(j, j + primer_len)
        gc_count = int(np.sum(np.isin(consensus[idx], ["g", "c"])))
        kmers.append(
            {
                "pos": j,
                "degen": "".join(pos_code[idx]),
                "n_degen": int(np.sum(pos_fold[idx] > 1)),
                "fold": int(np.prod(pos_fold[idx])),
                "divs": float(np.sum(divs[idx])),
                "GC": gc_count / primer_len,
            }
        )
    return kmers


def _pair_sort_key(row):
    return (-row["combined_score"], row["total_fold"], row["fwd_pos"], row["rev_pos"])


def _evaluate_pairs(
    candidates,
    amplicon_min_len,
    amplicon_max_len,
    GC_tol,
    max_results,
    score_weight_fold=0.0,
):
    """Return the best primer pairs without retaining an unbounded result set.

    When ``score_weight_fold`` > 0 a degeneracy penalty is added:
        score = 1 / (|pair_div| + 10*delta_gc^2 + w_fold*log2(total_fold) + 0.01)
    so that high-degeneracy (high-fold) primer pairs are de-ranked. With
    ``score_weight_fold == 0`` the score keeps the original (no fold term) form.
    """
    results = []
    candidates = sorted(candidates, key=lambda item: item["pos"])
    positions = [candidate["pos"] for candidate in candidates]

    for i, ci in enumerate(candidates[:-1]):
        lo = ci["pos"] + amplicon_min_len
        hi = ci["pos"] + amplicon_max_len
        start = bisect_left(positions, lo, i + 1)
        end = bisect_right(positions, hi, i + 1)

        for cj in candidates[start:end]:
            amp_len = cj["pos"] - ci["pos"]
            delta_gc = abs(ci["GC"] - cj["GC"])
            if delta_gc >= GC_tol:
                continue

            pair_div = ci["divs"] + cj["divs"]
            total_fold = ci["fold"] * cj["fold"]
            fold_penalty = (
                score_weight_fold * np.log2(total_fold) if score_weight_fold else 0.0
            )
            score = 1.0 / (abs(pair_div) + 10.0 * delta_gc**2 + fold_penalty + 0.01)

            results.append(
                {
                    "fwd": ci["degen"],
                    "rev": _rev_comp(cj["degen"]),
                    "fwd_pos": ci["pos"],
                    "rev_pos": cj["pos"],
                    "amplicon_len": amp_len,
                    "fwd_GC": round(ci["GC"], 4),
                    "rev_GC": round(cj["GC"], 4),
                    "fwd_ndegen": ci["n_degen"],
                    "rev_ndegen": cj["n_degen"],
                    "fwd_fold": ci["fold"],
                    "rev_fold": cj["fold"],
                    "total_fold": total_fold,
                    "pair_diversity": round(pair_div, 4),
                    "delta_GC": round(delta_gc, 4),
                    "combined_score": round(score, 6),
                }
            )
            if len(results) >= max_results * 2:
                results.sort(key=_pair_sort_key)
                del results[max_results:]

    results.sort(key=_pair_sort_key)
    del results[max_results:]
    for idx, row in enumerate(results, 1):
        row["primer_id"] = f"primer_pair_{idx}"
    return results


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------
def _write_tsv(rows, path):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=_TSV_FIELDNAMES, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def _write_empty_tsv(path):
    _write_tsv([], path)


# Visualization routines live in primers_plot (imported above): plot_placeholder,
# plot_diversity, plot_pair_heatmap. matplotlib/numpy are imported lazily there
# via ensure_matplotlib(), which main calls before any plotting.


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser(description="Design conserved primer pairs.")
    parser.add_argument("--aln", required=True, help="Input aligned FASTA")
    parser.add_argument("--out-tsv", required=True, help="Output primer TSV")
    parser.add_argument("--out-plot", required=True, help="Output diversity PNG")
    parser.add_argument("--config", help="Optional AmPair config.yaml")
    parser.add_argument("--gene", help="Gene name, used for per-gene config overrides")
    parser.add_argument("--primer-len", type=int)
    parser.add_argument("--amplicon-min-len", type=int)
    parser.add_argument("--amplicon-max-len", type=int)
    parser.add_argument("--div-cut", type=float)
    parser.add_argument("--gc-tol", type=float)
    parser.add_argument("--min-allele-freq", type=float)
    parser.add_argument("--max-degeneracy", type=int)
    parser.add_argument("--max-primer-pairs", type=int)
    parser.add_argument(
        "--score-weight-fold",
        type=float,
        help="Degeneracy penalty weight in combined_score (0 = no fold term).",
    )
    parser.add_argument("--log", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    started = perf_counter()

    global np, plt
    import matplotlib.pyplot as plt
    import numpy as np

    ensure_matplotlib()  # initialises plt/np inside primers_plot for plotting

    aln_file = args.aln
    out_tsv = args.out_tsv
    out_plot = args.out_plot
    cfg = load_config_file(args.config) if args.config else {}

    configure_logging(args.log)
    log = logging.getLogger(__name__)

    primer_len = _required_param(
        "primer_len", _param(args.primer_len, cfg, "primer_len")
    )
    amplicon_min_len = _required_param(
        "amplicon_min_len", _param(args.amplicon_min_len, cfg, "amplicon_min_len")
    )
    amplicon_max_len = _required_param(
        "amplicon_max_len", _param(args.amplicon_max_len, cfg, "amplicon_max_len")
    )
    if args.div_cut is not None:
        div_cut = args.div_cut
    elif args.gene:
        div_overrides = cfg.get("div_cut_per_gene", {})
        if isinstance(div_overrides, dict):
            div_cut = div_overrides.get(args.gene, cfg.get("div_cut"))
        else:
            div_cut = cfg.get("div_cut")
    else:
        div_cut = cfg.get("div_cut")
    div_cut = _required_param("div_cut", div_cut)
    GC_tol = _required_param("GC_tol", _param(args.gc_tol, cfg, "GC_tol"))
    min_allele_freq = _required_param(
        "min_allele_freq", _param(args.min_allele_freq, cfg, "min_allele_freq")
    )
    max_degeneracy = _required_param(
        "max_degeneracy", _param(args.max_degeneracy, cfg, "max_degeneracy")
    )
    max_primer_pairs = _required_param(
        "max_primer_pairs", _param(args.max_primer_pairs, cfg, "max_primer_pairs")
    )
    if max_primer_pairs < 1:
        raise SystemExit("max_primer_pairs must be a positive integer")

    # --- Scientific refinement settings -----------------------------------
    div_cut_auto_min = cfg.get("div_cut_auto_min_candidates", None)
    div_cut_auto_step = cfg.get("div_cut_auto_step", 0.05)
    div_cut_auto_max = cfg.get("div_cut_auto_max", 3.0)
    score_weight_fold = _param(args.score_weight_fold, cfg, "score_weight_fold") or 0.0
    plot_pair_heatmap_top_n = cfg.get("plot_pair_heatmap_top_n", 0)

    log.info("Parameters:")
    for k, v in [
        ("aln_file", aln_file),
        ("primer_len", primer_len),
        ("amplicon_min_len", amplicon_min_len),
        ("amplicon_max_len", amplicon_max_len),
        ("div_cut", div_cut),
        ("GC_tol", GC_tol),
        ("min_allele_freq", min_allele_freq),
        ("max_degeneracy", max_degeneracy),
        ("max_primer_pairs", max_primer_pairs),
    ]:
        log.info("  %-18s = %s", k, v)

    # --- 1. Load alignment ------------------------------------------------
    records = list(parse_fasta(aln_file))
    if len(records) < 2:
        msg = f"Need at least 2 sequences to estimate diversity; found {len(records)}."
        log.warning(msg)
        _write_empty_tsv(out_tsv)
        plot_placeholder(msg, aln_file, out_plot, log)
        return

    dna_matrix = np.array([list(sequence.lower()) for _, sequence in records])
    n_seqs, aln_len = dna_matrix.shape

    log.info("Loaded %d sequences, alignment length %d bp", n_seqs, aln_len)
    # --- 2. Per-column stats (entropy, consensus, IUPAC code + fold) ------
    consensus, divs, pos_code, pos_fold = _compute_column_stats(
        dna_matrix, n_seqs, min_allele_freq
    )
    log.info("Mean per-position entropy: %.4f", np.mean(divs))

    # --- 3. Degeneracy summary -------------------------------------------
    log.info(
        "Degenerate columns (fold > 1): %d / %d", int(np.sum(pos_fold > 1)), aln_len
    )

    # --- 4. Rolling mean for diversity plot -------------------------------
    roll_k = min(10, aln_len)
    roll_means = _rollmean(divs, roll_k)

    # --- 5. Build kmer table ----------------------------------------------
    if aln_len < primer_len:
        log.warning(
            "Alignment (%d bp) < primer_len (%d bp). Empty TSV.", aln_len, primer_len
        )
        _write_empty_tsv(out_tsv)
        plot_diversity(
            divs, roll_means, roll_k, [], 0, primer_len, aln_file, out_plot, log
        )
        return

    kmers = _build_kmers(consensus, pos_code, pos_fold, divs, primer_len)
    if len(kmers) > WARN_KMER_COUNT:
        log.warning(
            "Large primer search space (%d kmers). Consider stricter div_cut "
            "or a narrower amplicon range for faster batch runs.",
            len(kmers),
        )

    # --- 6. Filter candidates --------------------------------------------
    # Starting cutoff. With adaptive relaxation (div_cut_auto_min_candidates
    # set), the cutoff is progressively raised until enough candidates pass or
    # div_cut_auto_max is reached. Otherwise the fixed div_cut is used.
    effective_div_cut = div_cut
    if div_cut_auto_min is not None:
        while True:
            candidates = [
                k
                for k in kmers
                if k["divs"] <= effective_div_cut and k["fold"] <= max_degeneracy
            ]
            if (
                len(candidates) >= div_cut_auto_min
                or effective_div_cut >= div_cut_auto_max
            ):
                break
            effective_div_cut += div_cut_auto_step
        if effective_div_cut != div_cut:
            log.info(
                "Auto-relaxed div_cut %.2f -> %.2f (%d candidates pass, target %d)",
                div_cut,
                effective_div_cut,
                len(candidates),
                div_cut_auto_min,
            )
    else:
        candidates = [
            k
            for k in kmers
            if k["divs"] <= effective_div_cut and k["fold"] <= max_degeneracy
        ]

    if not candidates:
        log.info(
            "No candidates pass div_cut=%.2f + max_degeneracy=%d. Empty TSV.",
            effective_div_cut,
            max_degeneracy,
        )
        _write_empty_tsv(out_tsv)
        plot_diversity(
            divs, roll_means, roll_k, [], 0, primer_len, aln_file, out_plot, log
        )
        return

    log.info("%d candidate kmers pass filters", len(candidates))

    # --- 7. Evaluate pairs ------------------------------------------------
    results = _evaluate_pairs(
        candidates,
        amplicon_min_len,
        amplicon_max_len,
        GC_tol,
        max_primer_pairs,
        score_weight_fold=score_weight_fold,
    )
    if len(results) == max_primer_pairs:
        log.warning(
            "Primer pair results were limited to the best %d pairs; "
            "raise max_primer_pairs to retain more candidates.",
            max_primer_pairs,
        )

    if not results:
        log.info("No valid primer pairs found. Empty TSV.")
        _write_empty_tsv(out_tsv)
        plot_diversity(
            divs, roll_means, roll_k, [], 0, primer_len, aln_file, out_plot, log
        )
        return

    # --- 8. Write TSV -----------------------------------------------------
    _write_tsv(results, out_tsv)
    log.info("Wrote %d primer pairs -> %s", len(results), out_tsv)

    # --- 8b. Optional candidate-pair score heatmap ------------------------
    if plot_pair_heatmap_top_n and len(results) >= 2:
        heatmap_path = os.path.splitext(out_plot)[0] + "_heatmap.png"
        plot_pair_heatmap(
            results, min(plot_pair_heatmap_top_n, len(results)), heatmap_path, log
        )

    # --- 9. Diversity plot ------------------------------------------------
    plot_diversity(
        divs,
        roll_means,
        roll_k,
        results,
        min(5, len(results)),
        primer_len,
        aln_file,
        out_plot,
        log,
    )
    log.info("Primer design completed in %.2f s", perf_counter() - started)


if __name__ == "__main__":
    main()
