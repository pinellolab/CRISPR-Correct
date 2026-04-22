"""Smoke tests that run in CI — they cover the API surface but don't require
real FASTQ fixtures. Full ground-truth regression lives in the wrapper folder
(`CRISPR-Correct-Folder/tests/simulation/`) and is run out-of-band."""
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
