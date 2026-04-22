"""Smoke tests that run in CI — they cover the API surface but don't require
real FASTQ fixtures. Full ground-truth regression lives in the wrapper folder
(`CRISPR-Correct-Folder/tests/simulation/`) and is run out-of-band."""
from pathlib import Path

import pandas as pd
import pytest


def test_package_imports():
    import crispr_ambiguous_mapping as cam
    assert hasattr(cam, "mapping")
    assert hasattr(cam, "processing")
    assert hasattr(cam, "models")
    assert hasattr(cam, "quality_control")


def test_match_tier_enum_is_str_backed():
    from crispr_ambiguous_mapping.models.mapping_models import MatchTier
    # §4.2 / §7.5: subclassing `str` keeps magic-string callers working.
    assert MatchTier.PM_SM_BM == "protospacer_match_surrogate_match_barcode_match"
    assert MatchTier.PM == "protospacer_match"
    assert MatchTier.PM_MISMATCH_SM == "protospacer_mismatch_surrogate_match"


def test_guide_count_error_has_suggestion():
    from crispr_ambiguous_mapping.models.error_models import (
        ProtospacerHammingThresholdGuideCountError, BarcodeMissingInfoGuideCountError,
    )
    e = ProtospacerHammingThresholdGuideCountError(hamming_min=8, hamming_threshold=7)
    assert "threshold" in e.suggestion.lower()
    assert len(e.suggestion) > 20
    # Missing-info errors should also have suggestions.
    m = BarcodeMissingInfoGuideCountError()
    assert "extraction" in m.suggestion.lower()


def test_whitelist_validation_catches_missing_columns():
    from crispr_ambiguous_mapping.processing.crispr_guide_counting import (
        get_whitelist_reporter_counts_with_umi,
    )
    # `spacer` instead of `protospacer` — should raise ValueError with a clear message.
    from collections import Counter
    bad_df = pd.DataFrame({"spacer": ["ACGT"], "surrogate": ["ACGTACGT"]})
    obs: Counter = Counter({("ACGT", "ACGTACGT"): 1})
    with pytest.raises(ValueError, match="missing required columns"):
        get_whitelist_reporter_counts_with_umi(
            observed_guide_reporter_umi_counts=obs,
            whitelist_guide_reporter_df=bad_df,
            contains_guide_surrogate=True,
        )


def test_slim_result_raises_on_postproc():
    # §2.2: post-proc functions raise ValueError when called on a slim result
    # (retain_inference_results=False default).
    from crispr_ambiguous_mapping.processing.crispr_editing_processing import (
        get_matchset_alleleseries,
    )
    with pytest.raises(ValueError, match="retain_inference_results"):
        get_matchset_alleleseries(
            observed_guide_reporter_umi_counts_inferred=None,
            attribute_name="protospacer_match",
            contains_guide_surrogate=False,
            contains_guide_barcode=False,
            contains_guide_umi=False,
        )


def test_v0_1_0_public_api_importable():
    # §4.5 / §7.1: map_fastq / count / alleles / ParsingConfig exposed at
    # package root via api.py.
    import crispr_ambiguous_mapping as cam
    assert callable(cam.map_fastq)
    assert callable(cam.count)
    assert callable(cam.alleles)
    assert cam.ParsingConfig is not None

    cfg = cam.ParsingConfig(protospacer_length=20, cores=2)
    kw = cfg.to_kwargs()
    assert kw["protospacer_length"] == 20
    assert kw["cores"] == 2
    # None fields are dropped.
    assert "protospacer_pattern_regex" not in kw


def test_crispr_correct_shim_package():
    # §4.11: `import crispr_correct as cc` aliases the canonical name.
    import crispr_correct as cc
    assert callable(cc.map_fastq)
    assert cc.ParsingConfig is not None
    assert cc.MatchTier.PM_SM_BM == "protospacer_match_surrogate_match_barcode_match"


def test_count_input_contains_surrogate_backcompat():
    # §4.3 / K.1: legacy `contains_surrogate` attribute still reads the new
    # `contains_guide_surrogate` field.
    import pandas as pd
    from crispr_ambiguous_mapping.models.mapping_models import CountInput
    ci = CountInput(
        whitelist_guide_reporter_df=pd.DataFrame(),
        contains_guide_surrogate=True,
        contains_guide_barcode=False,
        contains_guide_umi=False,
        contains_sample_barcode=False,
        protospacer_hamming_threshold_strict=7,
        surrogate_hamming_threshold_strict=10,
        guide_barcode_hamming_threshold_strict=2,
    )
    assert ci.contains_guide_surrogate is True
    assert ci.contains_surrogate is True  # deprecated alias


def test_cli_help_runs():
    # §4.5: basic CLI smoke — the map subcommand's --help renders.
    from click.testing import CliRunner
    from crispr_ambiguous_mapping.cli import main
    runner = CliRunner()
    r = runner.invoke(main, ["map", "--help"])
    assert r.exit_code == 0, r.output
    # Flat-arg design check: every ParsingConfig field has a flag.
    assert "--protospacer-start-position" in r.output
    assert "--is-protospacer-r1 / --no-is-protospacer-r1" in r.output
    assert "--retain-inference-results" in r.output


def test_cli_parsingconfig_option_count_matches_fields():
    # Regression: every ParsingConfig field should become a CLI flag.
    from dataclasses import fields as _f
    from click.testing import CliRunner
    from crispr_ambiguous_mapping.cli import main
    from crispr_ambiguous_mapping.api import ParsingConfig
    runner = CliRunner()
    r = runner.invoke(main, ["map", "--help"])
    # At least one --flag per ParsingConfig field (boolean ones produce --flag/--no-flag).
    for f in _f(ParsingConfig):
        flag = "--" + f.name.replace("_", "-")
        assert flag in r.output, f"missing CLI flag for ParsingConfig field `{f.name}`"


def test_sample_barcode_inferred_value_is_shared():
    """§2.6: in sample-barcode (cell-barcode) mode, two cells that saw the
    same observed guide should share the `inferred_value` by reference — not
    get independent CompleteInferenceMatchResult copies. Verifies on the real
    scCRISPR result if available; otherwise skips."""
    import pickle
    import os
    golden = "tests/golden/sccrispr_cell_barcode.pkl"  # built by test_sccrispr_cell_barcode.py on first run
    if not os.path.exists(golden):
        pytest.skip("scCRISPR golden pickle not built yet")
    with open(golden, "rb") as fh:
        result = pickle.load(fh)
    inf = getattr(result, "observed_guide_reporter_umi_counts_inferred", None)
    if inf is None:
        pytest.skip("golden pickle was built with retain_inference_results=False")
    # Nested dict: {cell: {observed_tuple: InferenceResult}}. Find two cells that share an observed.
    cells = list(inf.keys())[:500]
    for c1 in cells:
        for obs in inf[c1]:
            for c2 in cells:
                if c2 == c1:
                    continue
                if obs in inf[c2]:
                    assert id(inf[c1][obs].inferred_value) == id(inf[c2][obs].inferred_value), \
                        "two cells saw the same observed tuple but got independent inferred_value objects"
                    return
    pytest.skip("no cross-cell shared observed tuple found in sampled cells")


def test_save_load_roundtrip_bit_equal(tmp_path):
    """§7.4: map -> save -> load round-trip must preserve count Series values +
    index across every non-None tier × strategy."""
    from collections import Counter
    import pandas as pd
    import crispr_ambiguous_mapping as cam

    HERE = Path(__file__).parent
    library = pd.read_csv(HERE / "fixtures" / "library.tsv", sep="\t")[["protospacer", "surrogate", "barcode"]]
    r1 = str(HERE / "fixtures" / "simulated_R1.fastq")
    r2 = str(HERE / "fixtures" / "simulated_R2.fastq")
    cfg = cam.ParsingConfig(
        protospacer_start_position=0, protospacer_length=20,
        is_protospacer_r1=True, is_protospacer_header=False, revcomp_protospacer=False,
        protospacer_hamming_threshold_strict=7,
        cores=2,
    )
    r = cam.map_fastq(library, [r1], [r2], config=cfg)
    out_dir = tmp_path / "result"
    cam.save(r, out_dir)
    r2_loaded = cam.load(out_dir)

    a1 = r.all_match_set_whitelist_reporter_counter_series_results
    a2 = r2_loaded.all_match_set_whitelist_reporter_counter_series_results
    compared = 0
    for tn in ("protospacer_match", "protospacer_match_surrogate_match", "protospacer_match_barcode_match",
               "protospacer_match_surrogate_match_barcode_match"):
        t1 = getattr(a1, tn, None); t2 = getattr(a2, tn, None)
        if t1 is None or t2 is None:
            continue
        for fname in dir(t1):
            if not fname.endswith("counterseries") or fname.startswith("_"):
                continue
            s1 = getattr(t1, fname, None); s2 = getattr(t2, fname, None)
            if s1 is None and s2 is None:
                continue
            compared += 1
            s1s = s1.sort_index(); s2s = s2.sort_index()
            assert s1s.index.equals(s2s.index), f"{tn}.{fname}: index mismatch"
            # Compare values loosely on dtype — parquet can widen int64→float64
            # when NaN alignment is involved; values themselves must match.
            if s1s.dtype == s2s.dtype:
                assert (s1s.values == s2s.values).all(), f"{tn}.{fname}: value mismatch"
            else:
                assert s1s.astype(float).equals(s2s.astype(float)), f"{tn}.{fname}: value mismatch"
    assert compared > 0


def test_cli_save_load_roundtrip(tmp_path):
    """§7.4: `crispr-correct save` + `crispr-correct load` preserve the result."""
    import pickle
    from click.testing import CliRunner
    from crispr_ambiguous_mapping.cli import main
    import crispr_ambiguous_mapping as cam

    HERE = Path(__file__).parent
    library = pd.read_csv(HERE / "fixtures" / "library.tsv", sep="\t")[["protospacer", "surrogate", "barcode"]]
    r1 = str(HERE / "fixtures" / "simulated_R1.fastq")
    r2 = str(HERE / "fixtures" / "simulated_R2.fastq")
    cfg = cam.ParsingConfig(protospacer_start_position=0, protospacer_length=20,
                            is_protospacer_r1=True, is_protospacer_header=False,
                            revcomp_protospacer=False, protospacer_hamming_threshold_strict=7, cores=2)
    result = cam.map_fastq(library, [r1], [r2], config=cfg)

    src_pkl = tmp_path / "src.pkl"
    with open(src_pkl, "wb") as fh:
        pickle.dump(result, fh)

    save_dir = tmp_path / "saved"
    runner = CliRunner()
    r_save = runner.invoke(main, ["save", "--in", str(src_pkl), "--out-dir", str(save_dir)])
    assert r_save.exit_code == 0, r_save.output
    manifest = save_dir / "manifest.json"
    assert manifest.exists()

    dst_pkl = tmp_path / "dst.pkl"
    r_load = runner.invoke(main, ["load", "--in-dir", str(save_dir), "--out", str(dst_pkl)])
    assert r_load.exit_code == 0, r_load.output
    assert dst_pkl.exists()
    with open(dst_pkl, "rb") as fh:
        loaded = pickle.load(fh)
    # Verify one Series matches after the CLI round-trip.
    s1 = result.all_match_set_whitelist_reporter_counter_series_results.protospacer_match.ambiguous_accepted_counterseries
    s2 = loaded.all_match_set_whitelist_reporter_counter_series_results.protospacer_match.ambiguous_accepted_counterseries
    if s1 is not None:
        assert s2 is not None
        assert s1.sort_index().index.equals(s2.sort_index().index)


def test_revcomp_translate_matches_biopython():
    # §3.10: translate-based revcomp must produce the same result as the
    # previous Bio.Seq.reverse_complement() on IUPAC bases we actually see.
    from crispr_ambiguous_mapping.parsing.reporter_standard_fastq_parsing import _RCMAP
    for seq, expected in [
        ("ACGT", "ACGT"),
        ("AAAA", "TTTT"),
        ("GATTACA", "TGTAATC"),
        ("N", "N"),
        ("acgtN", "Nacgt"),  # revcomp of "acgtN" = translate + reverse
    ]:
        assert seq.translate(_RCMAP)[::-1] == expected
