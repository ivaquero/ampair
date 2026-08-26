#!/usr/bin/env python3
# =============================================================================
# primers_plot.py - visualization helpers for primer_design.py
#
# Contains the plotting/figure routines so the design logic stays free of
# matplotlib concerns. matplotlib is configured with a non-interactive Agg
# backend so plots can be written to files in headless runs.
#
# Public API:
#   plot_placeholder(...)    -> text-only figure when no candidates pass
#   plot_diversity(...)      -> per-position entropy + rolling mean + Top-N sites
#   plot_pair_heatmap(...)   -> Top-N candidate-pair combined_score heatmap
# =============================================================================

from __future__ import annotations

import os

import matplotlib.pyplot as plt
import numpy as np


def plot_placeholder(message, aln_file, out_plot, log):
    os.makedirs(os.path.dirname(out_plot), exist_ok=True)
    fig, ax = plt.subplots(figsize=(12, 6), constrained_layout=True)
    ax.text(0.5, 0.5, message, ha="center", va="center", wrap=True)
    ax.set(axis_off=True, title=f"Sequence diversity - {os.path.basename(aln_file)}")
    fig.savefig(out_plot, dpi=150)
    plt.close(fig)
    log.info(f"Wrote placeholder diversity plot to {out_plot}")


def plot_diversity(
    divs, roll_means, roll_k, results, top_n, primer_len, aln_file, out_plot, log
):
    os.makedirs(os.path.dirname(out_plot), exist_ok=True)
    fig, ax = plt.subplots(figsize=(12, 6), constrained_layout=True)

    ax.scatter(
        np.arange(len(divs)),
        divs,
        s=4,
        alpha=0.5,
        c="black",
        edgecolors="none",
        label="Per-position entropy",
    )

    valid = ~np.isnan(roll_means)
    ax.plot(
        np.arange(len(divs))[valid],
        roll_means[valid],
        color="#2ca25f",
        lw=1.5,
        label=f"Rolling mean (k={roll_k})",
    )

    ymax = float(max(divs)) * 1.05
    for k in range(top_n):
        ax.axvspan(
            results[k]["fwd_pos"],
            results[k]["fwd_pos"] + primer_len,
            alpha=0.20,
            color="#e34a33",
        )
        ax.axvspan(
            results[k]["rev_pos"],
            results[k]["rev_pos"] + primer_len,
            alpha=0.20,
            color="#e34a33",
        )

    from matplotlib.patches import Patch

    handles, _ = ax.get_legend_handles_labels()
    if top_n > 0:
        handles.append(
            Patch(facecolor="#e34a33", alpha=0.20, label=f"Top {top_n} primer sites")
        )
    if handles:
        ax.legend(handles=handles, loc="upper right")

    ax.set(
        title=f"Sequence diversity - {os.path.basename(aln_file)}",
        xlabel="Alignment position (bp)",
        ylabel="Shannon entropy",
        ylim=(-0.05, ymax),
    )

    fig.savefig(out_plot, dpi=150)
    plt.close(fig)

    log.info(f"Wrote diversity plot to {out_plot}")


def plot_pair_heatmap(results, top_n, out_plot, log):
    """Score heatmap of the Top-N candidate primer pairs.

    Rows/cols are the candidate pair ids (ordered by score) and the cell color
    is combined_score, so reviewers can see how separable the best pairs are.
    """
    if top_n <= 0 or len(results) < 2:
        return

    top = results[:top_n]
    n = len(top)
    score_mat = np.zeros((n, n))
    labels = [r["primer_id"] for r in top]
    for a in range(n):
        for b in range(n):
            # score of an (unordered) pair is the better of (a,b)/(b,a); since the
            # matrix is symmetric for the same pair we just fill upper-tri + diag.
            score_mat[a, b] = (
                top[a]["combined_score"] + top[b]["combined_score"]
            ) / 2.0

    os.makedirs(os.path.dirname(out_plot), exist_ok=True)
    fig, ax = plt.subplots(
        figsize=(max(6, n * 0.6), max(5, n * 0.6)), constrained_layout=True
    )
    im = ax.imshow(score_mat, aspect="auto", cmap="viridis")
    ax.set(
        xticks=range(n),
        yticks=range(n),
        xticklabels=labels,
        yticklabels=labels,
        title="Top primer-pair combined_score heatmap",
    )
    ax.tick_params(axis="x", labelrotation=90, labelsize=7)
    ax.tick_params(axis="y", labelsize=7)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="combined_score")
    fig.savefig(out_plot, dpi=150)
    plt.close(fig)
    log.info(f"Wrote pair-score heatmap ({n} pairs) to {out_plot}")
