"""§2.5 measurement helper — profile per-Series byte footprint on a real
mapping-result pickle.

Usage:

    python measure_series_overhead.py /path/to/result.pickle

Prints a table of (tier, strategy) → bytes (values + index), the total over
all non-None Series, and the pickle size. Point this at the 11 k-guide chrX
result (or the 45 M-read AVITI result, if available) to decide whether the
54-Series → DataFrame dedup proposed in IMPROVEMENTS.md §2.5 is worth the
shape-changing refactor.

Decision threshold suggested by the Phase 6 measurement (~100 KB for a
45-guide sim):
- If total < 10 MB on real-scale data: close §2.5 as not worth the refactor.
- If total > 100 MB: schedule the refactor in Phase 8.
- In between: judgment call based on per-user RSS headroom.
"""
from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path


def _series_bytes(s) -> int:
    """Approximate in-memory byte cost of a pandas Series (values + index)."""
    total = int(getattr(s, "values", b"").nbytes or 0)
    idx = getattr(s, "index", None)
    if idx is None:
        return total
    if hasattr(idx, "levels"):
        for lv in idx.levels:
            total += int(lv.values.nbytes)
    else:
        total += int(idx.values.nbytes)
    return total


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("pickle_path", type=Path, help="Path to a WhitelistReporterCountsResult pickle.")
    args = ap.parse_args()

    path = args.pickle_path
    with path.open("rb") as fh:
        result = pickle.load(fh)

    allw = result.all_match_set_whitelist_reporter_counter_series_results
    tiers = [
        "protospacer_match",
        "protospacer_match_surrogate_match",
        "protospacer_match_barcode_match",
        "protospacer_match_surrogate_match_barcode_match",
        "protospacer_mismatch_surrogate_match",
        "protospacer_mismatch_surrogate_match_barcode_match",
    ]

    print(f"{'tier':50s} {'strategy':45s} {'rows':>8s} {'bytes':>12s}")
    print("-" * 120)
    total = 0
    n_series = 0
    for tn in tiers:
        t = getattr(allw, tn, None)
        if t is None:
            continue
        for attr in sorted(a for a in dir(t) if a.endswith("counterseries") and not a.startswith("_")):
            s = getattr(t, attr, None)
            if s is None:
                continue
            b = _series_bytes(s)
            total += b
            n_series += 1
            print(f"{tn:50s} {attr:45s} {len(s):>8d} {b:>12,d}")

    print("-" * 120)
    print(f"{n_series} non-None Series; total = {total/1024/1024:.2f} MB")

    try:
        pkl_size = path.stat().st_size
        print(f"Source pickle size on disk: {pkl_size/1024/1024:.2f} MB")
    except OSError:
        pass


if __name__ == "__main__":
    sys.exit(main() or 0)
