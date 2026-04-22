"""§7.7: CI-locked ground-truth simulation regression.

Runs the simulated-FASTQ mapping through `map_fastq` and compares against
`compute_expected_counts` across 8 parse modes × 6 tiers × 3 ambiguity × 3 UMI
strategies (135 comparisons). The full matrix takes ~60 s; the "not slow" CI
subset runs a representative 24 comparisons in under 15 s.

Fixtures live in-repo at `tests/fixtures/` so this test runs in GitHub Actions
without external data.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

import crispr_ambiguous_mapping as cam

from tests.simulation.expected import compute_expected_counts
from tests.simulation.helpers import run_local, compare_series, extract_series


HERE = Path(__file__).parent
FIXTURES = HERE / "fixtures"

# The 8 parse modes: (surrogate, barcode, umi).
_PARSE_MODES = [
    (s, b, u)
    for s in (False, True)
    for b in (False, True)
    for u in (False, True)
]

# Tier names and the attribute pattern for non-mismatch tiers.
_MATCH_TIERS = [
    "protospacer_match",
    "protospacer_match_surrogate_match",
    "protospacer_match_barcode_match",
    "protospacer_match_surrogate_match_barcode_match",
]
_MISMATCH_TIERS = [
    "protospacer_mismatch_surrogate_match",
    "protospacer_mismatch_surrogate_match_barcode_match",
]

_AMBIGUITIES = ["ignored", "accepted", "spread"]
_UMI_STRATEGIES = ["", "collapsed", "noncollapsed"]


def _tier_is_applicable(tier: str, surrogate: bool, barcode: bool) -> bool:
    if "surrogate_match" in tier and not surrogate:
        return False
    if "barcode_match" in tier and not barcode:
        return False
    return True


def _mode_label(s, b, u):
    return f"s{int(s)}b{int(b)}u{int(u)}"


def _build_cases():
    """Enumerate the 135 (mode, tier, ambiguity, umi_strategy) comparisons."""
    cases = []
    for surrogate, barcode, umi in _PARSE_MODES:
        # Mismatch tiers require surrogate; MM+BM requires barcode.
        applicable_tiers = [t for t in _MATCH_TIERS + _MISMATCH_TIERS if _tier_is_applicable(t, surrogate, barcode)]
        for tier in applicable_tiers:
            is_mismatch = tier.startswith("protospacer_mismatch")
            for amb in _AMBIGUITIES:
                for umi_strat in _UMI_STRATEGIES:
                    if umi_strat != "" and not umi:
                        continue
                    if umi and umi_strat == "":
                        continue
                    if is_mismatch:
                        # Mismatch tier has MATCH and MISMATCH sub-series.
                        for mm in ("MATCH", "MISMATCH"):
                            cases.append((surrogate, barcode, umi, f"{tier}:{mm}", amb, umi_strat))
                    else:
                        cases.append((surrogate, barcode, umi, tier, amb, umi_strat))
    return cases


# ---------------------------------------------------------------------------
# Session-scoped fixtures — build the mapping results once per parse mode.
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session")
def _library():
    df = pd.read_csv(FIXTURES / "library.tsv", sep="\t")
    return df[["protospacer", "surrogate", "barcode"]]


@pytest.fixture(scope="session")
def _library_full():
    return pd.read_csv(FIXTURES / "library.tsv", sep="\t")


@pytest.fixture(scope="session")
def _truth():
    return pd.read_parquet(FIXTURES / "truth.parquet")


@pytest.fixture(scope="session")
def _results(_library):
    """Map once per parse mode; cache results for all tier × strategy comparisons."""
    r1 = str(FIXTURES / "simulated_R1.fastq")
    r2 = str(FIXTURES / "simulated_R2.fastq")
    out = {}
    for surrogate, barcode, umi in _PARSE_MODES:
        out[(surrogate, barcode, umi)] = run_local(
            cam, _library, r1, r2, contains_surrogate=surrogate, contains_barcode=barcode, contains_umi=umi, cores=2,
        )
    return out


# ---------------------------------------------------------------------------
# The parametrized regression.
# ---------------------------------------------------------------------------


_CASES = _build_cases()


def _case_id(c):
    s, b, u, t, a, us = c
    return f"{_mode_label(s,b,u)}::{t}::{a}::{us or 'none'}"


# "slow" marker: the full-triplet+UMI modes (`s1b1u1`) each take ~2-3 s to build,
# compounded across tier combos. The CI-fast subset drops the two heaviest modes
# to keep the parametrized run under ~15 s.
def _is_slow(c):
    s, b, u, t, a, us = c
    return s and b and u


def _param(c):
    marks = [pytest.mark.slow] if _is_slow(c) else []
    return pytest.param(c, marks=marks, id=_case_id(c))


@pytest.mark.parametrize("case", [_param(c) for c in _CASES])
def test_simulated_counts(case, _library_full, _truth, _results):
    surrogate, barcode, umi, tier, amb, umi_strat = case
    observed = extract_series(_results[(surrogate, barcode, umi)], tier, amb, umi_strat)
    expected = compute_expected_counts(_truth, _library_full, surrogate, barcode, umi)
    # expected is a dict of (tier, amb, umi_strat) -> Series (with :MATCH/:MISMATCH for mismatch tiers).
    key = (tier, amb, umi_strat)
    expected_series = expected.get(key)
    result = compare_series(observed, expected_series)
    assert result["status"] == "PASS", f"{_case_id(case)}: {result}"
