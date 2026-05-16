"""Ground-truth expected counter Series computed from the simulation truth_df.

We re-derive what each of CRISPR-Correct's output Series *should* be, purely
from the truth dataframe + library — without using the package itself.
Compared against CRISPR-Correct's actual output, any mismatch is a bug.

The logic mirrors the inference pipeline at a high level:

  1. For each read, the *observed tuple* is the subset of
     (protospacer, surrogate, barcode) requested by the parsing mode.
  2. Reads are grouped by the observed tuple; within a group, UMIs form a
     Counter.
  3. For each group, library matches are computed independently per component
     using strict-less-than Hamming (same as the package). The per-component
     match sets are combined by tier:

       - protospacer_match                             = PM_match_set
       - protospacer_match_surrogate_match             = PM ∩ SM        (non-empty)
       - protospacer_match_barcode_match               = PM ∩ BM        (non-empty)
       - protospacer_match_surrogate_match_barcode_match = PM ∩ SM ∩ BM
       - protospacer_mismatch_surrogate_match          = PM / SM  when PM ∩ SM = ∅ and both non-empty  → mismatch tier
       - protospacer_mismatch_surrogate_match_barcode_match = similar, AND BM ⊇ SM (i.e. SM ⊂ BM? the tool requires
                                                              SM and BM both non-empty, and PM ∩ SM ∩ BM empty)

  4. For each (tier, ambiguity-strategy, umi-strategy) combination, increment
     per-guide counters accordingly:

       ambiguous_accepted  → every match gets the group's weight
       ambiguous_ignored   → only groups with unique match (|matches|==1) contribute
       ambiguous_spread    → every match gets weight / |matches|

       weight(umi_collapsed)     = len(umi_counter)
       weight(umi_noncollapsed)  = sum(umi_counter.values())
       weight(no_umi)            = total reads in group (== noncollapsed when UMI is off)
"""
from __future__ import annotations

from collections import Counter, defaultdict
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


def _hamming_against(obs: str, lib_seqs_np: np.ndarray) -> np.ndarray:
    """lib_seqs_np is a (N, L) uint8 array. obs is a str of length L.
    Returns a length-N numpy int array of hamming distances.

    Assumes all sequences in lib_seqs_np have the same length as obs;
    truncates/pads if not (rare: only protospacer lengths may differ by 1-2 bp).
    """
    L = lib_seqs_np.shape[1]
    obs_padded = obs[:L].ljust(L, "N")
    obs_np = np.frombuffer(obs_padded.encode(), dtype=np.uint8)
    return (lib_seqs_np != obs_np).sum(axis=1)


def _match_set(obs: str, lib_seqs_np: np.ndarray, threshold: int, subset_idx: Optional[List[int]] = None) -> List[int]:
    """Strict-less-than: distances < threshold. Return indices of *all* rows at
    the minimum distance (not just rows below threshold). If subset_idx is
    provided, search is restricted to those library indices (mirrors the
    nested-subset logic used in crispr_guide_inference.infer_whitelist_sequence).
    """
    if lib_seqs_np is None:
        return []
    if subset_idx is not None:
        if not subset_idx:
            return []
        sub_np = lib_seqs_np[subset_idx]
    else:
        sub_np = lib_seqs_np
    dists = _hamming_against(obs, sub_np)
    mn = int(dists.min())
    if mn >= threshold:
        return []
    local_matches = np.where(dists == mn)[0].tolist()
    if subset_idx is not None:
        return [subset_idx[i] for i in local_matches]
    return local_matches


def _pack_library(lib_df: pd.DataFrame, col: str) -> np.ndarray:
    strs = lib_df[col].tolist()
    L = max(len(s) for s in strs)
    padded = [s.ljust(L, "N") for s in strs]
    return np.frombuffer("".join(padded).encode(), dtype=np.uint8).reshape(len(strs), L).copy()


def compute_expected_counts(
    truth_df: pd.DataFrame,
    library_df: pd.DataFrame,
    contains_surrogate: bool,
    contains_barcode: bool,
    contains_umi: bool,
    protospacer_hamming_threshold: int = 7,
    surrogate_hamming_threshold: int = 10,
    barcode_hamming_threshold: int = 2,
) -> Dict[Tuple[str, str, str], pd.Series]:
    """Return a dict keyed by (tier, ambiguity_strategy, umi_strategy) → Series.

    Series are MultiIndexed on the library's (protospacer, surrogate, barcode)
    rows, matching the package's output index. Only guides with nonzero counts
    appear in the Series — again matching the package after our perf fix.

    `umi_strategy` is one of 'collapsed', 'noncollapsed', ''. The '' key is
    used for the non-UMI flavor of each tier (the plain "counterseries").

    Mismatch tiers use a different key-shape (pair of library tuples); those
    Series are also MultiIndexed accordingly.
    """
    if contains_umi and "umi" not in truth_df.columns:
        raise ValueError("contains_umi=True but truth_df has no 'umi' column")

    # Pack library sequences
    proto_np = _pack_library(library_df, "protospacer")
    surr_np  = _pack_library(library_df, "surrogate") if contains_surrogate else None
    bc_np    = _pack_library(library_df, "barcode")   if contains_barcode   else None

    # Previously the §1.1 surrogate-truncation bug was active upstream, so we
    # mirrored it here by clamping `surr_np = surr_np[:, :proto_np.shape[1]]`.
    # The bug has now been fixed on perf/roadmap-phase1 — the surrogate uses
    # its full 32 bp length — so the mirror has been removed. If this test
    # ever runs against an older package where the bug is live, reinstate the
    # clamp.

    # Which tiers are produced given the mode? Mismatch tiers get split into
    # two pseudo-tiers (:MATCH / :MISMATCH) because the package emits two
    # separate Series per (ambiguity, umi) combo in those tiers.
    tiers = ["protospacer_match"]
    if contains_surrogate:
        tiers.append("protospacer_match_surrogate_match")
        tiers.append("protospacer_mismatch_surrogate_match:MATCH")
        tiers.append("protospacer_mismatch_surrogate_match:MISMATCH")
        if contains_barcode:
            tiers.append("protospacer_match_surrogate_match_barcode_match")
            tiers.append("protospacer_mismatch_surrogate_match_barcode_match:MATCH")
            tiers.append("protospacer_mismatch_surrogate_match_barcode_match:MISMATCH")
    if contains_barcode:
        tiers.append("protospacer_match_barcode_match")

    umi_strategies = ["collapsed", "noncollapsed"] if contains_umi else [""]
    ambig_strategies = ["accepted", "ignored", "spread"]

    counters: Dict[Tuple[str, str, str], Dict] = {}
    for tier in tiers:
        for amb in ambig_strategies:
            for umi in umi_strategies:
                counters[(tier, amb, umi)] = defaultdict(float)

    # Group reads by observed tuple (according to parsing mode)
    key_cols = ["observed_protospacer"]
    if contains_surrogate:
        key_cols.append("observed_surrogate")
    if contains_barcode:
        key_cols.append("observed_guide_barcode")

    # Precompute library tuple for each library index
    lib_cols = ["protospacer", "surrogate", "barcode"]
    lib_tuples = list(library_df[lib_cols].itertuples(index=False, name=None))

    groups = truth_df.groupby(key_cols, sort=False)

    for key_vals, sub in groups:
        obs_proto = key_vals if not isinstance(key_vals, tuple) else key_vals[0]
        obs_surr = key_vals[1] if (contains_surrogate and isinstance(key_vals, tuple)) else None
        obs_bc   = key_vals[-1] if contains_barcode else None
        if not contains_surrogate and not contains_barcode:
            obs_proto = key_vals  # single-column groupby returns scalar

        # Read-level weights for this group
        total_reads = len(sub)
        if contains_umi:
            umi_counts = Counter(sub["umi"])
            weight_by_umi = {"noncollapsed": sum(umi_counts.values()), "collapsed": len(umi_counts)}
        else:
            weight_by_umi = {"": total_reads}

        # Library match sets — mirror the nested-subset logic in
        # crispr_guide_inference.infer_whitelist_sequence. The order matters:
        # BM first (full library), then PM within BM subset, then SM within
        # the PM-within-BM subset for the full-triplet tier. PM+SM (no BM)
        # uses full-library PM then SM within PM subset.

        pm_full = _match_set(obs_proto, proto_np, protospacer_hamming_threshold)

        # BM full-library (independent of PM)
        bm_full = _match_set(obs_bc, bc_np, barcode_hamming_threshold) if contains_barcode else None

        # PM+BM tier: restrict PM to BM subset
        if contains_barcode:
            pm_in_bm = _match_set(obs_proto, proto_np, protospacer_hamming_threshold, subset_idx=bm_full) if bm_full else []
        else:
            pm_in_bm = []

        # SM surrogate-only (full library)
        sm_full = _match_set(obs_surr, surr_np, surrogate_hamming_threshold) if contains_surrogate else None

        # PM+SM tier: SM restricted to PM subset
        if contains_surrogate:
            sm_in_pm = _match_set(obs_surr, surr_np, surrogate_hamming_threshold, subset_idx=pm_full) if pm_full else []
        else:
            sm_in_pm = []

        # PM+SM+BM tier: SM restricted to (PM∩BM) subset
        if contains_surrogate and contains_barcode:
            sm_in_pm_bm = _match_set(obs_surr, surr_np, surrogate_hamming_threshold, subset_idx=pm_in_bm) if pm_in_bm else []
        else:
            sm_in_pm_bm = []

        per_tier_matches: Dict[str, List] = {}
        per_tier_matches["protospacer_match"] = pm_full
        if contains_surrogate:
            per_tier_matches["protospacer_match_surrogate_match"] = sm_in_pm
        if contains_barcode:
            per_tier_matches["protospacer_match_barcode_match"] = pm_in_bm
        if contains_surrogate and contains_barcode:
            per_tier_matches["protospacer_match_surrogate_match_barcode_match"] = sm_in_pm_bm

        # Mismatch tiers actually emit TWO series per ambiguity×UMI combo:
        #
        #   *_match_counterseries    (3-level key) — populated when PM_full and
        #                            SM_full (or PM_full and SM_in_BM) intersect
        #                            non-trivially. Counts each shared guide.
        #   *_mismatch_counterseries (6-level key = (proto_tup, surr_tup)) —
        #                            populated when the intersection is empty:
        #                            every (p, s) cross pair is counted.
        #
        # We track them separately using tier keys with a suffix "_MATCH" or
        # "_MISMATCH" so the same per-tier-matches iteration loop can drive both.
        if contains_surrogate:
            # plain PM-mismatch-SM tier
            inter = [i for i in pm_full if i in set(sm_full)] if (pm_full and sm_full) else []
            if inter:
                per_tier_matches["protospacer_mismatch_surrogate_match:MATCH"] = inter
                per_tier_matches["protospacer_mismatch_surrogate_match:MISMATCH"] = []
            elif pm_full and sm_full:
                per_tier_matches["protospacer_mismatch_surrogate_match:MATCH"] = []
                per_tier_matches["protospacer_mismatch_surrogate_match:MISMATCH"] = [(p, s) for p in pm_full for s in sm_full]
            else:
                per_tier_matches["protospacer_mismatch_surrogate_match:MATCH"] = []
                per_tier_matches["protospacer_mismatch_surrogate_match:MISMATCH"] = []

            # PM-mismatch-SM-BM tier: SM is taken within BM subset; PM stays full-library.
            if contains_barcode:
                sm_in_bm = _match_set(obs_surr, surr_np, surrogate_hamming_threshold, subset_idx=bm_full) if bm_full else []
                inter_bm = [i for i in pm_full if i in set(sm_in_bm)] if (pm_full and sm_in_bm) else []
                if inter_bm:
                    per_tier_matches["protospacer_mismatch_surrogate_match_barcode_match:MATCH"] = inter_bm
                    per_tier_matches["protospacer_mismatch_surrogate_match_barcode_match:MISMATCH"] = []
                elif pm_full and sm_in_bm:
                    per_tier_matches["protospacer_mismatch_surrogate_match_barcode_match:MATCH"] = []
                    per_tier_matches["protospacer_mismatch_surrogate_match_barcode_match:MISMATCH"] = [(p, s) for p in pm_full for s in sm_in_bm]
                else:
                    per_tier_matches["protospacer_mismatch_surrogate_match_barcode_match:MATCH"] = []
                    per_tier_matches["protospacer_mismatch_surrogate_match_barcode_match:MISMATCH"] = []

        # Now populate counters. The :MISMATCH pseudo-tiers have 6-level pair
        # keys; all others have 3-level single-guide keys.
        for tier, matches in per_tier_matches.items():
            if not matches:
                continue
            n = len(matches)
            is_pair_keyed = tier.endswith(":MISMATCH")
            for amb in ambig_strategies:
                if amb == "ignored" and n != 1:
                    continue
                for umi_strat, weight in weight_by_umi.items():
                    per_match = (weight / n) if amb == "spread" else weight
                    for m in matches:
                        if is_pair_keyed:
                            key = (lib_tuples[m[0]], lib_tuples[m[1]])
                        else:
                            key = lib_tuples[m]
                        counters[(tier, amb, umi_strat)][key] += per_match

    # Build pd.Series. Match tiers: 3-level MultiIndex. :MISMATCH tiers: 6-level.
    result: Dict[Tuple[str, str, str], pd.Series] = {}
    for (tier, amb, umi_strat), cdict in counters.items():
        if not cdict:
            result[(tier, amb, umi_strat)] = pd.Series(dtype=float)
            continue

        is_pair_keyed = tier.endswith(":MISMATCH")
        if is_pair_keyed:
            records = [(*k[0], *k[1], v) for k, v in cdict.items()]
            cols = [c + "_ProtospacerMatch" for c in lib_cols] + [c + "_SurrogateMatch" for c in lib_cols] + ["value"]
            df = pd.DataFrame.from_records(records, columns=cols)
            result[(tier, amb, umi_strat)] = df.set_index(cols[:-1])["value"]
        else:
            records = [(*k, v) for k, v in cdict.items()]
            cols = lib_cols + ["value"]
            df = pd.DataFrame.from_records(records, columns=cols)
            result[(tier, amb, umi_strat)] = df.set_index(cols[:-1])["value"]

    return result
