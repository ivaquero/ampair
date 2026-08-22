#!/usr/bin/env python3
# =============================================================================
# primers_check.py
#
# Primer quality filter: computes hairpin, homodimer, heterodimer, and
# 3'-end stability for every primer pair, then drops pairs that fail any
# configurable threshold.  Writes a filtered TSV with extra quality columns.
#
# All dG calculations use SantaLucia 1998 nearest-neighbour parameters in
# pure Python; no primer3 C library required. Works on any OS.
#
# CLI:
#   python primers_check.py --in-tsv primers_raw.tsv --out-tsv primers.tsv \
#       --max-hairpin-dg 0 --max-homodimer-dg -6 ...
# =============================================================================

import argparse
import csv
import logging
import os

from common import IUPAC_COMPLEMENT_TABLE, reverse_complement
from config_schema import load_config_file

# ---------------------------------------------------------------------------
# IUPAC ambiguity codes
# ---------------------------------------------------------------------------
IUPAC_BASES = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "M": frozenset("AC"),
    "K": frozenset("GT"),
    "S": frozenset("CG"),
    "W": frozenset("AT"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}

# ---------------------------------------------------------------------------
# SantaLucia 1998 unified nearest-neighbour dG (kcal/mol, 37 C, 1 M NaCl)
# ---------------------------------------------------------------------------
NN_DG = {
    "AA": -1.00,
    "TT": -1.00,
    "AT": -0.88,
    "TA": -0.58,
    "CA": -1.45,
    "TG": -1.45,
    "GT": -1.44,
    "AC": -1.44,
    "CT": -1.28,
    "AG": -1.28,
    "GA": -1.30,
    "TC": -1.30,
    "CG": -2.17,
    "GC": -2.24,
    "GG": -1.84,
    "CC": -1.84,
}

SYMMETRY_DG = 0.43  # self-complementary duplex penalty
INIT_AT = 2.30  # initiation, terminal A/T
INIT_GC = 2.55  # initiation, terminal G/C
MIN_STEM = 2  # minimum stem length for hairpin
MIN_LOOP = 3  # minimum loop size for hairpin


# ---------------------------------------------------------------------------
# Building blocks
# ---------------------------------------------------------------------------
def _complement(base: str) -> str:
    return base.translate(IUPAC_COMPLEMENT_TABLE)


def _base_set(base: str):
    return IUPAC_BASES.get(base.upper(), frozenset())


def _compatible_after_reverse_complement(base_a: str, base_b_rc: str) -> bool:
    """True when two 5' to 3' bases can pair in an antiparallel duplex."""
    return bool(_base_set(base_a) & _base_set(base_b_rc))


def _terminal_init_dg(base: str) -> float:
    bases = _base_set(base)
    values = []
    if bases & frozenset("AT"):
        values.append(INIT_AT)
    if bases & frozenset("GC"):
        values.append(INIT_GC)
    return min(values) if values else 0.0


def _nearest_neighbor_dg(dimer: str) -> float:
    """Most stable canonical nearest-neighbour value compatible with dimer."""
    left, right = dimer[0], dimer[1]
    values = [
        dg
        for a in _base_set(left)
        for b in _base_set(right)
        if (dg := NN_DG.get(a + b)) is not None
    ]
    return min(values) if values else 0.0


def _duplex_dg(seq_5p: str) -> float:
    """dG of a duplex scored from one strand (5' to 3'), SantaLucia 1998.

    seq_5p is the sequence of one strand in 5' to 3' orientation. The
    complementary strand is assumed to pair perfectly.
    """
    n = len(seq_5p)
    if n < 2:
        return 0.0
    dg = _terminal_init_dg(seq_5p[0])
    for i in range(n - 1):
        dg += _nearest_neighbor_dg(seq_5p[i : i + 2])
    if seq_5p == reverse_complement(seq_5p):
        dg += SYMMETRY_DG
    return dg


# ---------------------------------------------------------------------------
# Dimer dG (homo- and hetero-)
# ---------------------------------------------------------------------------
def _align_dg(seq_a: str, seq_b_rc: str) -> float:
    """Best (most negative) dG of antiparallel alignment of seq_a vs seq_b_rc.

    seq_b_rc is the reverse complement of some other sequence in 5' to 3'.
    We slide the two 5' to 3' strands across each other and score only
    contiguous compatible stems.
    """
    best = 0.0
    seq_a = seq_a.upper()
    seq_b_rc = seq_b_rc.upper()
    la, lb = len(seq_a), len(seq_b_rc)

    for offset in range(-(lb - 1), la):
        if offset >= 0:
            a0, b0, n = offset, 0, min(la - offset, lb)
        else:
            a0, b0, n = 0, -offset, min(la, lb + offset)

        if n < MIN_STEM:
            continue

        run = []
        for k in range(n):
            base_a = seq_a[a0 + k]
            base_b_rc = seq_b_rc[b0 + k]
            if _compatible_after_reverse_complement(base_a, base_b_rc):
                run.append(base_a)
                continue

            if len(run) >= MIN_STEM:
                best = min(best, _duplex_dg("".join(run)))
            run = []

        if len(run) >= MIN_STEM:
            best = min(best, _duplex_dg("".join(run)))

    return round(best, 2)


def calc_homodimer_dg(seq: str) -> float:
    """Minimum homodimer dG (self-dimer)."""
    return _align_dg(seq, reverse_complement(seq))


def calc_heterodimer_dg(seq_a: str, seq_b: str) -> float:
    """Minimum heterodimer dG (cross-dimer)."""
    return _align_dg(seq_a, reverse_complement(seq_b))


# ---------------------------------------------------------------------------
# Hairpin dG
# ---------------------------------------------------------------------------
def calc_hairpin_dg(seq: str) -> float:
    """Minimum hairpin dG: find the most stable stem within the sequence.

    Enumerates all possible stem-start (i), stem-end-start (j) and stem
    length (k) combinations with a minimum loop of MIN_LOOP bases between the
    two stem halves. Only genuinely complementary stems are scored.
    """
    seq = seq.upper()
    n = len(seq)
    best = 0.0

    for stem_len in range(MIN_STEM, n // 2 + 1):
        max_i = n - 2 * stem_len - MIN_LOOP
        for i in range(max_i + 1):
            j_min = i + stem_len + MIN_LOOP
            for j in range(j_min, n - stem_len + 1):
                stem5 = seq[i : i + stem_len]
                stem3_rc = reverse_complement(seq[j : j + stem_len])
                if all(
                    _compatible_after_reverse_complement(a, b)
                    for a, b in zip(stem5, stem3_rc, strict=True)
                ):
                    best = min(best, _duplex_dg(stem5))

    return round(best, 2)


# ---------------------------------------------------------------------------
# 3'-end stability
# ---------------------------------------------------------------------------
def calc_3end_dg(seq: str, n_bases: int = 5) -> float:
    """dG of the 3'-most *n_bases* of *seq* (kcal/mol, SantaLucia 1998)."""
    end = seq[-n_bases:].upper()
    if len(end) < 2:
        return 0.0

    dg = INIT_GC if end[0] in "GC" else INIT_AT
    for i in range(len(end) - 1):
        dg += NN_DG.get(end[i : i + 2], 0)

    rc = "".join(_complement(b) for b in reversed(end))
    if end == rc:
        dg += SYMMETRY_DG

    return round(dg, 2)


# ---------------------------------------------------------------------------
# Quality-check a single primer pair
# ---------------------------------------------------------------------------
def check_pair(row: dict, thresholds: dict) -> dict:
    fwd = row["fwd"]
    rev = row["rev"]

    metrics = {
        "hairpin_fwd_dg": calc_hairpin_dg(fwd),
        "hairpin_rev_dg": calc_hairpin_dg(rev),
        "homodimer_fwd_dg": calc_homodimer_dg(fwd),
        "homodimer_rev_dg": calc_homodimer_dg(rev),
        "heterodimer_dg": calc_heterodimer_dg(fwd, rev),
        "end3_fwd_dg": calc_3end_dg(fwd),
        "end3_rev_dg": calc_3end_dg(rev),
    }

    # thresholds: more negative dG = stronger secondary structure = worse
    checks = [
        ("hairpin_fwd_dg", thresholds.get("max_hairpin_dg")),
        ("hairpin_rev_dg", thresholds.get("max_hairpin_dg")),
        ("homodimer_fwd_dg", thresholds.get("max_homodimer_dg")),
        ("homodimer_rev_dg", thresholds.get("max_homodimer_dg")),
        ("heterodimer_dg", thresholds.get("max_heterodimer_dg")),
        ("end3_fwd_dg", thresholds.get("max_3end_dg")),
        ("end3_rev_dg", thresholds.get("max_3end_dg")),
    ]

    fails = []
    for key, limit in checks:
        if limit is not None and metrics[key] < limit:
            fails.append(key)

    result: dict[str, float | bool | str] = dict(metrics)
    result["qc_pass"] = len(fails) == 0
    result["qc_fail_reasons"] = ";".join(fails) if fails else ""

    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def _nullable_float(value):
    if value.lower() in {"none", "null", "na"}:
        return None
    return float(value)


def parse_args():
    parser = argparse.ArgumentParser(description="Filter primer pairs by QC metrics.")
    parser.add_argument("--in-tsv", required=True)
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--config", help="Optional AmPrime config.yaml for thresholds")
    parser.add_argument(
        "--max-hairpin-dg", type=_nullable_float, default=argparse.SUPPRESS
    )
    parser.add_argument(
        "--max-homodimer-dg", type=_nullable_float, default=argparse.SUPPRESS
    )
    parser.add_argument(
        "--max-heterodimer-dg", type=_nullable_float, default=argparse.SUPPRESS
    )
    parser.add_argument(
        "--max-3end-dg", type=_nullable_float, default=argparse.SUPPRESS
    )
    parser.add_argument("--log", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    log_path = args.log
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    logging.basicConfig(
        filename=log_path,
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    log = logging.getLogger()

    in_tsv = args.in_tsv
    out_tsv = args.out_tsv

    thresholds = {
        "max_hairpin_dg": None,
        "max_homodimer_dg": None,
        "max_heterodimer_dg": None,
        "max_3end_dg": None,
    }
    if args.config:
        cfg = load_config_file(args.config)
        thresholds.update({key: cfg.get(key) for key in thresholds})

    cli_thresholds = {
        "max_hairpin_dg": "max_hairpin_dg",
        "max_homodimer_dg": "max_homodimer_dg",
        "max_heterodimer_dg": "max_heterodimer_dg",
        "max_3end_dg": "max_3end_dg",
    }
    for threshold_key, arg_name in cli_thresholds.items():
        if hasattr(args, arg_name):
            thresholds[threshold_key] = getattr(args, arg_name)

    log.info("Input  : %s", in_tsv)
    log.info("Output : %s", out_tsv)
    log.info("Thresholds: %s", {k: v for k, v in thresholds.items() if v is not None})

    # Read input
    with open(in_tsv, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)

    if not rows:
        log.warning("Input TSV is empty; writing empty output.")
        os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
        with open(out_tsv, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            # minimal header when there's no data
            if reader.fieldnames:
                extra = [
                    "hairpin_fwd_dg",
                    "hairpin_rev_dg",
                    "homodimer_fwd_dg",
                    "homodimer_rev_dg",
                    "heterodimer_dg",
                    "end3_fwd_dg",
                    "end3_rev_dg",
                    "qc_pass",
                    "qc_fail_reasons",
                ]
                w.writerow(list(reader.fieldnames) + extra)
        return

    in_fields = list(reader.fieldnames or [])
    qc_fields = [
        "hairpin_fwd_dg",
        "hairpin_rev_dg",
        "homodimer_fwd_dg",
        "homodimer_rev_dg",
        "heterodimer_dg",
        "end3_fwd_dg",
        "end3_rev_dg",
        "qc_pass",
        "qc_fail_reasons",
    ]
    out_fields = in_fields + qc_fields

    checked = []
    for r in rows:
        m = check_pair(r, thresholds)
        checked.append({**r, **m})

    passed = [r for r in checked if r["qc_pass"]]
    n_dropped = len(checked) - len(passed)
    log.info(
        "Checked %d pairs; %d passed, %d dropped.", len(checked), len(passed), n_dropped
    )

    if not passed:
        log.warning(
            "All %d pairs failed quality filters. Writing empty TSV.", len(checked)
        )
        os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
        with open(out_tsv, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(out_fields)
        return

    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    with open(out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=out_fields, delimiter="\t", extrasaction="ignore"
        )
        w.writeheader()
        w.writerows(passed)

    log.info("Written %d filtered pairs to %s", len(passed), out_tsv)


if __name__ == "__main__":
    main()
