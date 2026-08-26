#!/usr/bin/env python3
# =============================================================================
# design_compare.py
#
# Offline comparison of primer-design outputs between the current AmPair
# pipeline (workflow/scripts/primers_design.py) and the first-version
# "original" pipeline (original/workflow/scripts/design_primers.R).
#
# It does NOT run either pipeline; it only compares two already-produced
# primer TSVs so you can isolate *why* the two code bases diverge.
#
# The script is dependency-free (standard library only) so it runs anywhere.
#
# Usage:
#   python design_compare.py \
#       --current  results/Borrelia/primers/recG_primers.tsv \
#       --original original_run/recG_primers.tsv \
#       [--top-n 50]
#
# The two TSVs may use different column names:
#   current : primer_id, fwd, rev, fwd_pos, rev_pos, amplicon_len, combined_score
#   original: gene, fwd, rev, fwd.pos, rev.pos, amplicon.len, score
# Either naming is accepted; sequences are compared case-insensitively after
# uppercasing.
# =============================================================================

import argparse
import csv
import sys
from collections import OrderedDict


def _read_pairs(path):
    """Return an ordered list of (fwd, rev, score) tuples from a primer TSV."""
    pairs = []
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            return pairs
        # Detect column naming style.
        fields = {f.lower(): f for f in reader.fieldnames}
        fwd_col = fields.get("fwd") or fields.get("forward")
        rev_col = fields.get("rev") or fields.get("reverse")
        score_col = (
            fields.get("combined_score")
            or fields.get("score")
            or fields.get("score_combined")
        )
        if fwd_col is None or rev_col is None:
            raise SystemExit(
                f"Cannot find forward/reverse columns in {path}. "
                f"Found: {reader.fieldnames}"
            )
        for row in reader:
            fwd = (row.get(fwd_col) or "").upper().strip()
            rev = (row.get(rev_col) or "").upper().strip()
            if not fwd or not rev:
                continue
            score = None
            if score_col and row.get(score_col):
                try:
                    score = float(row[score_col])
                except ValueError:
                    score = None
            pairs.append((fwd, rev, score))
    return pairs


def _key(pair):
    return (pair[0], pair[1])


def _report(current, original, top_n):
    cur = current[:top_n]
    org = original[:top_n]

    cur_keys = OrderedDict((_key(p), p) for p in cur)
    org_keys = OrderedDict((_key(p), p) for p in org)

    shared = set(cur_keys) & set(org_keys)
    only_cur = set(cur_keys) - set(org_keys)
    only_org = set(org_keys) - set(cur_keys)

    cur_top = _key(cur[0]) if cur else None
    org_top = _key(org[0]) if org else None

    print("=" * 72)
    print("AmPair primer-design comparison (current vs original)")
    print("=" * 72)
    print(f"Top-N compared            : {top_n}")
    print(f"Current candidates (total): {len(current)}  (in top-N: {len(cur)})")
    print(f"Original candidates (total): {len(original)}  (in top-N: {len(org)})")
    print("-" * 72)
    print(f"Shared primer pairs       : {len(shared)}")
    print(f"Only in current           : {len(only_cur)}")
    print(f"Only in original          : {len(only_org)}")
    denom = max(1, len(set(cur_keys) | set(org_keys)))
    print(f"Jaccard overlap           : {len(shared) / denom:.3f}")
    print("-" * 72)
    print("Top-1 (recommended) pair:")
    print(f"  current : {cur_top}")
    print(f"  original: {org_top}")
    print(f"  match   : {'YES' if cur_top == org_top else 'NO'}")
    print("=" * 72)

    if only_cur:
        print(f"\nPairs only in CURRENT (top-N), first {min(10, len(only_cur))}:")
        for k in list(only_cur)[:10]:
            print(f"  fwd={k[0]} rev={k[1]}")
    if only_org:
        print(f"\nPairs only in ORIGINAL (top-N), first {min(10, len(only_org))}:")
        for k in list(only_org)[:10]:
            print(f"  fwd={k[0]} rev={k[1]}")

    # Exit non-zero when the recommended pairs differ, so CI / diff checks fail.
    return 0 if cur_top == org_top else 1


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--current",
        required=True,
        help="Primer TSV from the current pipeline (primers_design.py).",
    )
    parser.add_argument(
        "--original",
        required=True,
        help="Primer TSV from the original pipeline (design_primers.R).",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=50,
        help="Number of top-ranked pairs to compare (default: 50).",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    current = _read_pairs(args.current)
    original = _read_pairs(args.original)
    return _report(current, original, args.top_n)


if __name__ == "__main__":
    sys.exit(main())
