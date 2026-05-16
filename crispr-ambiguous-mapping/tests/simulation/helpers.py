"""Shared helpers for the test_local / test_0142 notebooks.

Two key functions:
  - run_local(...) drives the multi-sample-support API
    (`get_whitelist_reporter_counts_from_fastq`).
  - run_0142(...) drives the legacy 0.0.142 API
    (`get_whitelist_reporter_counts_from_umitools_output` — the only entry
    point available in that version). The umitools entry point always assumes
    protospacer at R1[0:20] and surrogate at R2[0:32] revcomp, so "surrogate off"
    is not a supported mode there.
  - compare_series(observed, expected) compares pandas Series on aligned
    MultiIndex values with an exact-equality check, returning a structured
    pass/fail record.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd


# Regexes that the simulator's header format satisfies
BARCODE_REGEX = r"_([^_ ]+)[\s+]"
UMI_REGEX = r":([^+:]{6})(.{2})\+"


# ---------------------------------------------------------------------------
# LOCAL (multi-sample-support) — 8 modes
# ---------------------------------------------------------------------------

def run_local(
    cam_module,
    library_df: pd.DataFrame,
    r1_path: str,
    r2_path: str,
    contains_surrogate: bool,
    contains_barcode: bool,
    contains_umi: bool,
    cores: int = 4,
):
    """Call `cam.mapping.get_whitelist_reporter_counts_from_fastq` with the
    parsing knobs turned on/off per the flags. Returns the full result."""

    kwargs: Dict[str, Any] = dict(
        whitelist_guide_reporter_df=library_df,
        fastq_r1_fns=[r1_path],
        fastq_r2_fns=[r2_path] if contains_surrogate else None,

        protospacer_start_position=0, protospacer_length=20,
        is_protospacer_r1=True, is_protospacer_header=False, revcomp_protospacer=False,

        protospacer_hamming_threshold_strict=7,
        # Keep the inference dict so the ground-truth tests can exercise
        # downstream allele paths (default is slim / None).
        retain_inference_results=True,
        cores=cores,
    )
    if contains_surrogate:
        kwargs.update(
            surrogate_start_position=0, surrogate_length=32,
            is_surrogate_r1=False, is_surrogate_header=False, revcomp_surrogate=True,
            surrogate_hamming_threshold_strict=10,
        )
    if contains_barcode:
        kwargs.update(
            guide_barcode_pattern_regex=BARCODE_REGEX,
            is_guide_barcode_r1=False, is_guide_barcode_header=True, revcomp_guide_barcode=True,
            guide_barcode_hamming_threshold_strict=2,
        )
    if contains_umi:
        kwargs.update(
            guide_umi_pattern_regex=UMI_REGEX,
            is_guide_umi_r1=False, is_guide_umi_header=True, revcomp_guide_umi=False,
        )

    return cam_module.mapping.get_whitelist_reporter_counts_from_fastq(**kwargs)


# ---------------------------------------------------------------------------
# 0142 (UMI-tools entry point) — 4 modes (surrogate always on)
# ---------------------------------------------------------------------------

def run_0142(
    cam_module,
    library_df: pd.DataFrame,
    r1_path: str,
    r2_path: str,
    contains_barcode: bool,
    contains_umi: bool,
    cores: int = 4,
):
    """The 0142 umitools entry point has fixed protospacer/surrogate positions;
    barcode + UMI come from the header via the same regexes."""
    kwargs: Dict[str, Any] = dict(
        whitelist_guide_reporter_df=library_df,
        fastq_r1_fn=r1_path,
        fastq_r2_fn=r2_path,
        protospacer_hamming_threshold_strict=7,
        surrogate_hamming_threshold_strict=10,
        barcode_hamming_threshold_strict=2,
        cores=cores,
    )
    if contains_barcode:
        kwargs["barcode_pattern_regex"] = BARCODE_REGEX
    if contains_umi:
        kwargs["umi_pattern_regex"] = UMI_REGEX

    return cam_module.mapping.get_whitelist_reporter_counts_from_umitools_output(**kwargs)


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------

def compare_series(
    observed: Optional[pd.Series],
    expected: Optional[pd.Series],
    *,
    tolerance: float = 1e-9,
) -> Dict[str, Any]:
    """Compare two pandas Series by aligning on index and doing exact value
    equality (within tolerance for float spread values). Returns
    {status: PASS|FAIL|SKIP, detail: str, ...}."""
    # 0142 sometimes emits a raw defaultdict instead of a pd.Series for some
    # tier/ambiguity combos; convert if so.
    if observed is not None and not isinstance(observed, pd.Series):
        try:
            observed = pd.Series(dict(observed))
        except Exception:
            return {"status": "ERROR", "detail": f"observed has unexpected type {type(observed).__name__}",
                    "n_obs": 0, "n_exp": len(expected) if expected is not None else 0}

    obs_empty = observed is None or len(observed) == 0
    exp_empty = expected is None or len(expected) == 0

    if obs_empty and exp_empty:
        return {"status": "PASS", "detail": "both empty", "n_obs": 0, "n_exp": 0}
    if obs_empty != exp_empty:
        return {
            "status": "FAIL",
            "detail": f"emptiness mismatch: obs empty={obs_empty} exp empty={exp_empty}",
            "n_obs": 0 if obs_empty else len(observed),
            "n_exp": 0 if exp_empty else len(expected),
        }

    # Harmonize index names (observed / expected may differ on MultiIndex level
    # names). Work on positionally-stripped indices to avoid name-comparison
    # failures inside pandas' index.union implementation.
    obs_c = observed.copy()
    exp_c = expected.copy()
    obs_c.index = obs_c.index.set_names([None] * obs_c.index.nlevels) if isinstance(obs_c.index, pd.MultiIndex) else obs_c.index
    exp_c.index = exp_c.index.set_names([None] * exp_c.index.nlevels) if isinstance(exp_c.index, pd.MultiIndex) else exp_c.index

    if obs_c.index.nlevels != exp_c.index.nlevels:
        return {
            "status": "FAIL",
            "detail": f"index nlevels mismatch: obs={obs_c.index.nlevels} exp={exp_c.index.nlevels}",
            "n_obs": len(observed), "n_exp": len(expected),
        }

    all_idx = obs_c.index.union(exp_c.index)
    obs_aligned = obs_c.reindex(all_idx, fill_value=0.0).astype(float)
    exp_aligned = exp_c.reindex(all_idx, fill_value=0.0).astype(float)

    diff = (obs_aligned - exp_aligned).abs()
    max_diff = float(diff.max())
    n_diff = int((diff > tolerance).sum())

    if n_diff == 0:
        return {
            "status": "PASS",
            "detail": f"all {len(all_idx)} keys match (max_diff={max_diff})",
            "n_obs": len(observed), "n_exp": len(expected),
        }
    # Top differing rows for diagnostics
    top = diff.sort_values(ascending=False).head(5)
    return {
        "status": "FAIL",
        "detail": f"{n_diff}/{len(all_idx)} keys differ; max_diff={max_diff}; top diffs: {top.to_dict()}",
        "n_obs": len(observed), "n_exp": len(expected),
    }


# ---------------------------------------------------------------------------
# Result-extraction helpers (walk the nested dataclasses)
# ---------------------------------------------------------------------------

_MATCH_TIER_FIELDS = [
    "ambiguous_ignored_umi_noncollapsed_counterseries",
    "ambiguous_ignored_umi_collapsed_counterseries",
    "ambiguous_ignored_counterseries",
    "ambiguous_accepted_umi_noncollapsed_counterseries",
    "ambiguous_accepted_umi_collapsed_counterseries",
    "ambiguous_accepted_counterseries",
    "ambiguous_spread_umi_noncollapsed_counterseries",
    "ambiguous_spread_umi_collapsed_counterseries",
    "ambiguous_spread_counterseries",
]

_MISMATCH_TIER_FIELDS = [
    # the mismatch-tier dataclass uses "_match_" in field names
    "ambiguous_ignored_umi_noncollapsed_match_counterseries",
    "ambiguous_ignored_umi_collapsed_match_counterseries",
    "ambiguous_ignored_match_counterseries",
    "ambiguous_accepted_umi_noncollapsed_match_counterseries",
    "ambiguous_accepted_umi_collapsed_match_counterseries",
    "ambiguous_accepted_match_counterseries",
    "ambiguous_spread_umi_noncollapsed_match_counterseries",
    "ambiguous_spread_umi_collapsed_match_counterseries",
    "ambiguous_spread_match_counterseries",
]

TIER_TO_ATTR = {
    "protospacer_match": "protospacer_match",
    "protospacer_match_surrogate_match": "protospacer_match_surrogate_match",
    "protospacer_match_barcode_match": "protospacer_match_barcode_match",
    "protospacer_match_surrogate_match_barcode_match": "protospacer_match_surrogate_match_barcode_match",
    "protospacer_mismatch_surrogate_match": "protospacer_mismatch_surrogate_match",
    "protospacer_mismatch_surrogate_match_barcode_match": "protospacer_mismatch_surrogate_match_barcode_match",
}


def extract_series(result, tier: str, ambiguity: str, umi_strategy: str):
    """Pull the named Series from result.all_match_set_whitelist_reporter_counter_series_results.

    For mismatch tiers, `tier` must carry a `:MATCH` or `:MISMATCH` suffix
    selecting which of the two Series (3-level `_match_counterseries` vs
    6-level `_mismatch_counterseries`) to fetch. Non-mismatch tiers ignore the
    suffix.
    """
    allw = result.all_match_set_whitelist_reporter_counter_series_results
    base_tier, _, mismatch_suffix = tier.partition(":")
    tier_wrap = getattr(allw, TIER_TO_ATTR[base_tier], None)
    if tier_wrap is None:
        return None

    is_mismatch_tier = base_tier.startswith("protospacer_mismatch")
    if is_mismatch_tier:
        series_kind = "match" if (mismatch_suffix or "MATCH") == "MATCH" else "mismatch"
        if umi_strategy == "":
            field = f"ambiguous_{ambiguity}_{series_kind}_counterseries"
        else:
            field = f"ambiguous_{ambiguity}_umi_{umi_strategy}_{series_kind}_counterseries"
    else:
        if umi_strategy == "":
            field = f"ambiguous_{ambiguity}_counterseries"
        else:
            field = f"ambiguous_{ambiguity}_umi_{umi_strategy}_counterseries"
    return getattr(tier_wrap, field, None)
