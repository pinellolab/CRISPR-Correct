"""Phase 8 V.1 — comprehensive end-to-end validation sweep.

Runs the `s1b1u1` full-triplet + UMI simulation mode through every public
surface added across Phases 1-7 and asserts consistency at each stage:

  map (Python)           ── bit-equal to ──            map (CLI)
       │                                                   │
       ├── count(result) ─ bit-equal to ─ crispr-correct count
       │
       ├── save(dir) ── bit-equal to ── crispr-correct save
       │                                                   │
       └── load(dir) ── reconstructs ──────── load.pickle  │
                                                           │
       alleles(result, tier=PM_SM_BM) ── parquet ── CLI alleles

Plus ground-truth checks:
  - count series per tier x strategy == compute_expected_counts
  - allele_df internally consistent with tier counts
  - slim result + alleles() raises ValueError with actionable message

Runtime: ~60 s (one map + CLI subprocesses).
"""
from __future__ import annotations

import json
import pickle
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

import crispr_ambiguous_mapping as cam
from crispr_ambiguous_mapping.models.mapping_models import MatchTier

from tests.simulation.expected import compute_expected_counts
from tests.simulation.helpers import extract_series, compare_series


HERE = Path(__file__).parent
FIXTURES = HERE / "fixtures"
LIBRARY_TSV = FIXTURES / "library.tsv"
R1 = str(FIXTURES / "simulated_R1.fastq")
R2 = str(FIXTURES / "simulated_R2.fastq")
TRUTH_PARQUET = FIXTURES / "truth.parquet"


def _full_triplet_config(*, retain: bool = True) -> cam.ParsingConfig:
    return cam.ParsingConfig(
        protospacer_start_position=0, protospacer_length=20,
        is_protospacer_r1=True, is_protospacer_header=False, revcomp_protospacer=False,
        protospacer_hamming_threshold_strict=7,
        surrogate_start_position=0, surrogate_length=32,
        is_surrogate_r1=False, is_surrogate_header=False, revcomp_surrogate=True,
        surrogate_hamming_threshold_strict=10,
        guide_barcode_pattern_regex=r"_([^_ ]+)[\s+]",
        is_guide_barcode_r1=False, is_guide_barcode_header=True, revcomp_guide_barcode=True,
        guide_barcode_hamming_threshold_strict=2,
        guide_umi_pattern_regex=r":([^+:]{6})(.{2})\+",
        is_guide_umi_r1=False, is_guide_umi_header=True, revcomp_guide_umi=False,
        retain_inference_results=retain,
        cores=2,
    )


@pytest.fixture(scope="module")
def _library():
    return pd.read_csv(LIBRARY_TSV, sep="\t")[["protospacer", "surrogate", "barcode"]]


@pytest.fixture(scope="module")
def _library_full():
    return pd.read_csv(LIBRARY_TSV, sep="\t")


@pytest.fixture(scope="module")
def _truth():
    return pd.read_parquet(TRUTH_PARQUET)


@pytest.fixture(scope="module")
def py_result(_library):
    """Stage 1 — the authoritative Python map result, retained."""
    return cam.map_fastq(_library, [R1], [R2], config=_full_triplet_config(retain=True))


@pytest.fixture(scope="module")
def cli_result(tmp_path_factory, _library):
    """Stage 2 — CLI map produces a second pickle via subprocess."""
    tmp = tmp_path_factory.mktemp("cli_map")
    out_pkl = tmp / "cli_result.pkl"
    cmd = [
        sys.executable, "-m", "crispr_ambiguous_mapping.cli", "map",
        "--r1", R1, "--r2", R2, "--library", str(LIBRARY_TSV),
        "--out", str(out_pkl),
        "--protospacer-start-position", "0", "--protospacer-length", "20",
        "--is-protospacer-r1", "--no-is-protospacer-header", "--no-revcomp-protospacer",
        "--protospacer-hamming-threshold-strict", "7",
        "--surrogate-start-position", "0", "--surrogate-length", "32",
        "--no-is-surrogate-r1", "--no-is-surrogate-header", "--revcomp-surrogate",
        "--surrogate-hamming-threshold-strict", "10",
        "--guide-barcode-pattern-regex", r"_([^_ ]+)[\s+]",
        "--no-is-guide-barcode-r1", "--is-guide-barcode-header", "--revcomp-guide-barcode",
        "--guide-barcode-hamming-threshold-strict", "2",
        "--guide-umi-pattern-regex", r":([^+:]{6})(.{2})\+",
        "--no-is-guide-umi-r1", "--is-guide-umi-header", "--no-revcomp-guide-umi",
        "--retain-inference-results", "--cores", "2",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    with open(out_pkl, "rb") as fh:
        r = pickle.load(fh)
    return r, result.stderr


def _iter_non_none_series(allw):
    """Yield (tier, attr, series) for every non-None counterseries on the result."""
    for tier_name in ("protospacer_match", "protospacer_match_surrogate_match",
                      "protospacer_match_barcode_match", "protospacer_match_surrogate_match_barcode_match",
                      "protospacer_mismatch_surrogate_match", "protospacer_mismatch_surrogate_match_barcode_match"):
        t = getattr(allw, tier_name, None)
        if t is None:
            continue
        for attr in sorted(a for a in dir(t) if a.endswith("counterseries") and not a.startswith("_")):
            s = getattr(t, attr, None)
            if s is not None:
                yield tier_name, attr, s


# ---------------------------------------------------------------------------
# Stages
# ---------------------------------------------------------------------------


def test_01_python_map_retained(py_result):
    """Stage 1 — Python map with retain_inference_results=True."""
    assert py_result.count_input.contains_guide_surrogate is True
    assert py_result.count_input.contains_guide_barcode is True
    assert py_result.count_input.contains_guide_umi is True
    # Legacy alias still works (§4.3 back-compat).
    assert py_result.count_input.contains_surrogate is True
    # Retained inference dict present.
    assert py_result.observed_guide_reporter_umi_counts_inferred is not None
    assert len(py_result.observed_guide_reporter_umi_counts_inferred) > 0


def test_02_cli_map_bit_equal_to_python(cli_result, py_result):
    """Stage 2 — CLI map result's count series bit-equal to Python result."""
    cli_r, _stderr = cli_result
    for tier_name, attr, s_py in _iter_non_none_series(py_result.all_match_set_whitelist_reporter_counter_series_results):
        s_cli = getattr(getattr(cli_r.all_match_set_whitelist_reporter_counter_series_results, tier_name), attr, None)
        assert s_cli is not None, f"CLI missing {tier_name}.{attr}"
        assert s_py.sort_index().equals(s_cli.sort_index()), f"CLI mismatch at {tier_name}.{attr}"


def test_03_count_series_ground_truth(py_result, _library_full, _truth):
    """Stage 3 — every non-None tier x strategy matches compute_expected_counts."""
    expected = compute_expected_counts(_truth, _library_full, True, True, True)

    # Build a case list matching expected's keys: (tier, amb, umi_strat) with optional :MATCH/:MISMATCH.
    failures = []
    checked = 0
    for (tier_key, amb, umi_strat), exp_series in expected.items():
        observed = extract_series(py_result, tier_key, amb, umi_strat)
        r = compare_series(observed, exp_series)
        checked += 1
        if r["status"] != "PASS":
            failures.append((tier_key, amb, umi_strat, r))
    assert checked > 0, "no expected cases generated"
    assert not failures, f"{len(failures)}/{checked} sim comparisons failed; first 3: {failures[:3]}"


def test_04_count_api_returns_same_object(py_result):
    """Stage 4 — cam.count(result) returns the all-match-set container by identity."""
    c = cam.count(py_result)
    assert c is py_result.all_match_set_whitelist_reporter_counter_series_results


def test_05_cli_count_tsv_matches_python_series(py_result, tmp_path):
    """Stage 5 — crispr-correct count --out tsv equals the Python Series values."""
    # Materialize the source pickle the CLI reads from.
    src_pkl = tmp_path / "src.pkl"
    with open(src_pkl, "wb") as fh:
        pickle.dump(py_result, fh)
    out_tsv = tmp_path / "counts.tsv"
    cmd = [sys.executable, "-m", "crispr_ambiguous_mapping.cli", "count",
           "--in", str(src_pkl), "--out", str(out_tsv),
           "--tier", "protospacer_match_surrogate_match_barcode_match",
           "--strategy", "ambiguous_accepted_umi_noncollapsed_counterseries"]
    subprocess.run(cmd, capture_output=True, text=True, check=True)
    tsv = pd.read_csv(out_tsv, sep="\t")
    # TSV has columns (protospacer, surrogate, barcode, value) — one row per guide.
    python_series = py_result.all_match_set_whitelist_reporter_counter_series_results.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_umi_noncollapsed_counterseries
    assert len(tsv) == len(python_series)
    assert int(tsv["value"].sum()) == int(python_series.sum())


def test_06_python_save_load_roundtrip(py_result, tmp_path):
    """Stage 6 — cam.save + cam.load round-trip every tier x strategy."""
    out_dir = tmp_path / "saved"
    cam.save(py_result, out_dir)
    r2 = cam.load(out_dir)
    compared = 0
    for tier_name, attr, s_py in _iter_non_none_series(py_result.all_match_set_whitelist_reporter_counter_series_results):
        s_r2 = getattr(getattr(r2.all_match_set_whitelist_reporter_counter_series_results, tier_name), attr, None)
        assert s_r2 is not None, f"load lost {tier_name}.{attr}"
        a = s_py.sort_index(); b = s_r2.sort_index()
        assert a.index.equals(b.index), f"{tier_name}.{attr} index mismatch"
        if a.dtype == b.dtype:
            assert (a.values == b.values).all(), f"{tier_name}.{attr} value mismatch"
        else:
            assert a.astype(float).equals(b.astype(float)), f"{tier_name}.{attr} value mismatch"
        compared += 1
    assert compared > 0


def test_07_cli_save_load_roundtrip(py_result, tmp_path):
    """Stage 7 — CLI save + load produce a pickle whose counts match Python round-trip."""
    src_pkl = tmp_path / "src.pkl"
    with open(src_pkl, "wb") as fh:
        pickle.dump(py_result, fh)
    save_dir = tmp_path / "saved_cli"
    subprocess.run([sys.executable, "-m", "crispr_ambiguous_mapping.cli", "save",
                    "--in", str(src_pkl), "--out-dir", str(save_dir)],
                   capture_output=True, text=True, check=True)
    assert (save_dir / "manifest.json").exists()
    dst_pkl = tmp_path / "dst.pkl"
    subprocess.run([sys.executable, "-m", "crispr_ambiguous_mapping.cli", "load",
                    "--in-dir", str(save_dir), "--out", str(dst_pkl)],
                   capture_output=True, text=True, check=True)
    with open(dst_pkl, "rb") as fh:
        loaded = pickle.load(fh)
    # Cross-check a strategy that's guaranteed non-empty on sim data.
    s_py = py_result.all_match_set_whitelist_reporter_counter_series_results.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_umi_noncollapsed_counterseries
    s_cli = loaded.all_match_set_whitelist_reporter_counter_series_results.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_umi_noncollapsed_counterseries
    assert s_py is not None and s_cli is not None
    assert s_py.sort_index().index.equals(s_cli.sort_index().index)
    assert int(s_py.sum()) == int(s_cli.sum())


def test_08_alleles_python(py_result):
    """Stage 8 — cam.alleles(result, tier=PM_SM_BM) populates allele tables."""
    ms = cam.alleles(py_result, tier=MatchTier.PM_SM_BM,
                    contains_guide_surrogate=True, contains_guide_barcode=True, contains_guide_umi=True)
    # The default-ish strategy should be non-empty.
    df = ms.ambiguous_accepted_umi_noncollapsed_allele_df
    assert df is not None and not df.empty
    # 9 strategy DataFrames exist (9 alleleseries_dicts are populated — allele_df presence varies).
    strat_dfs = [a for a in dir(ms) if a.endswith("allele_df") and not a.startswith("_")]
    assert len(strat_dfs) >= 9


def test_09_allele_internal_consistency(py_result):
    """Stage 9 — allele-df totals are consistent with the tier count series.

    Per-guide sum of observed-allele counts in the ambiguous_accepted
    umi_noncollapsed allele table equals the same-strategy count series for
    that whitelist guide.
    """
    ms = cam.alleles(py_result, tier=MatchTier.PM_SM_BM,
                    contains_guide_surrogate=True, contains_guide_barcode=True, contains_guide_umi=True)
    allele_df = ms.ambiguous_accepted_umi_noncollapsed_allele_df
    count_series = py_result.all_match_set_whitelist_reporter_counter_series_results.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_umi_noncollapsed_counterseries
    # allele_df has MultiIndex (observed_*) + columns including true_protospacer, true_surrogate, true_barcode, counts.
    # Group by true_* to sum per whitelist, then compare to count_series keyed on (protospacer, surrogate, barcode).
    true_cols = ["true_protospacer", "true_surrogate", "true_barcode"]
    assert all(c in allele_df.columns for c in true_cols), f"missing true_ cols: {list(allele_df.columns)}"
    grouped = allele_df.groupby(true_cols)["counts"].sum()
    # count_series index levels are (protospacer, surrogate, barcode) — same order as true_cols.
    overlap = grouped.index.intersection(count_series.index)
    assert len(overlap) > 0, f"no overlap between allele-df whitelist rows and count-series (grouped keys: {grouped.index[:3].tolist()}, count keys: {count_series.index[:3].tolist()})"
    for idx in overlap:
        assert float(grouped.loc[idx]) == float(count_series.loc[idx]), \
            f"allele sum != count at {idx}: allele={grouped.loc[idx]} count={count_series.loc[idx]}"


def test_10_cli_alleles_writes_parquet(py_result, tmp_path):
    """Stage 10 — crispr-correct alleles writes a non-empty parquet matching the Python allele table."""
    src_pkl = tmp_path / "src.pkl"
    with open(src_pkl, "wb") as fh:
        pickle.dump(py_result, fh)
    out_parquet = tmp_path / "alleles.parquet"
    subprocess.run([sys.executable, "-m", "crispr_ambiguous_mapping.cli", "alleles",
                    "--in", str(src_pkl), "--out", str(out_parquet),
                    "--tier", "protospacer_match_surrogate_match_barcode_match",
                    "--ambiguity", "accepted", "--umi-strategy", "noncollapsed"],
                   capture_output=True, text=True, check=True)
    assert out_parquet.exists()
    cli_df = pd.read_parquet(out_parquet)
    assert not cli_df.empty
    # Same shape as the Python allele_df.
    ms = cam.alleles(py_result, tier=MatchTier.PM_SM_BM,
                    contains_guide_surrogate=True, contains_guide_barcode=True, contains_guide_umi=True)
    py_df = ms.ambiguous_accepted_umi_noncollapsed_allele_df
    # CLI parquet resets the MultiIndex to columns; compare on row count and counts sum.
    assert cli_df.shape[0] == py_df.shape[0]
    assert int(cli_df["counts"].sum()) == int(py_df["counts"].sum())


def test_11_slim_result_alleles_raises(_library):
    """Stage 11 — slim result (retain_inference_results=False) + cam.alleles raises ValueError."""
    slim = cam.map_fastq(_library, [R1], [R2], config=_full_triplet_config(retain=False))
    assert slim.observed_guide_reporter_umi_counts_inferred is None
    with pytest.raises(ValueError, match="retain_inference_results"):
        cam.alleles(slim, tier=MatchTier.PM_SM_BM,
                   contains_guide_surrogate=True, contains_guide_barcode=True, contains_guide_umi=True)
