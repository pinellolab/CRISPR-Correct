"""Public API surface for CRISPR-Correct v0.1.0.

Exposes the three-stage workflow (map / count / alleles) cleanly, with
dataclass configuration to replace the 50-kwarg entry point. The legacy
`mapping.get_whitelist_reporter_counts_from_fastq(...)` entry point continues
to work and is what `map_fastq` delegates to — this module is currently a
thin re-packaging with IDE-friendly signatures; a future release can split
the stages more aggressively once drivers are migrated.

Typical usage:

```python
import crispr_ambiguous_mapping as cam
from crispr_ambiguous_mapping.api import map_fastq, count, alleles, ParsingConfig
from crispr_ambiguous_mapping.models import MatchTier

cfg = ParsingConfig(
    protospacer_start_position=0, protospacer_length=20,
    is_protospacer_r1=True, is_protospacer_header=False, revcomp_protospacer=False,
    protospacer_hamming_threshold_strict=7,
    surrogate_start_position=0, surrogate_length=32,
    is_surrogate_r1=False, is_surrogate_header=False, revcomp_surrogate=True,
    surrogate_hamming_threshold_strict=10,
    retain_inference_results=True,
    cores=4,
)

result = map_fastq(library, fastq_r1_fns=["R1.fq.gz"], fastq_r2_fns=["R2.fq.gz"], config=cfg)
counts_per_tier = count(result)
allele_df = alleles(result, tier=MatchTier.PM_SM_BM, contains_guide_surrogate=True, contains_guide_barcode=False, contains_guide_umi=False)
```
"""
from __future__ import annotations
from dataclasses import dataclass, asdict, fields
from typing import Any, Dict, List, Optional
import pandas as pd

from .models.mapping_models import (
    WhitelistReporterCountsResult,
    AllMatchSetWhitelistReporterCounterSeriesResults,
    MatchTier,
)


@dataclass
class ParsingConfig:
    """IDE-friendly bundle of the ~50 parsing/threshold kwargs.

    Every field mirrors a kwarg on
    `mapping.get_whitelist_reporter_counts_from_fastq`. Pass a ParsingConfig
    to `map_fastq` and it unpacks into the legacy signature.

    §4.1 / §7.3: this is the minimum-viable dataclass replacement — a single
    flat struct so IDE autocomplete works. A future release can decompose
    into ComponentConfig / ThresholdConfig subtrees once users have migrated.
    """
    # Regex / flank / position kwargs — protospacer
    protospacer_pattern_regex: Optional[str] = None
    protospacer_left_flank: Optional[str] = None
    protospacer_right_flank: Optional[str] = None
    protospacer_start_position: Optional[int] = None
    protospacer_end_position: Optional[int] = None
    protospacer_length: Optional[int] = None
    is_protospacer_r1: Optional[bool] = None
    is_protospacer_header: Optional[bool] = None
    revcomp_protospacer: Optional[bool] = None
    protospacer_hamming_threshold_strict: Optional[int] = None

    # surrogate
    surrogate_pattern_regex: Optional[str] = None
    surrogate_left_flank: Optional[str] = None
    surrogate_right_flank: Optional[str] = None
    surrogate_start_position: Optional[int] = None
    surrogate_end_position: Optional[int] = None
    surrogate_length: Optional[int] = None
    is_surrogate_r1: Optional[bool] = None
    is_surrogate_header: Optional[bool] = None
    revcomp_surrogate: Optional[bool] = None
    surrogate_hamming_threshold_strict: Optional[int] = None

    # guide barcode
    guide_barcode_pattern_regex: Optional[str] = None
    guide_barcode_left_flank: Optional[str] = None
    guide_barcode_right_flank: Optional[str] = None
    guide_barcode_start_position: Optional[int] = None
    guide_barcode_end_position: Optional[int] = None
    guide_barcode_length: Optional[int] = None
    is_guide_barcode_r1: Optional[bool] = None
    is_guide_barcode_header: Optional[bool] = None
    revcomp_guide_barcode: Optional[bool] = None
    guide_barcode_hamming_threshold_strict: Optional[int] = None

    # guide UMI
    guide_umi_pattern_regex: Optional[str] = None
    guide_umi_left_flank: Optional[str] = None
    guide_umi_right_flank: Optional[str] = None
    guide_umi_start_position: Optional[int] = None
    guide_umi_end_position: Optional[int] = None
    guide_umi_length: Optional[int] = None
    is_guide_umi_r1: Optional[bool] = None
    is_guide_umi_header: Optional[bool] = None
    revcomp_guide_umi: Optional[bool] = None

    # sample barcode
    sample_barcode_pattern_regex: Optional[str] = None
    sample_barcode_left_flank: Optional[str] = None
    sample_barcode_right_flank: Optional[str] = None
    sample_barcode_start_position: Optional[int] = None
    sample_barcode_end_position: Optional[int] = None
    sample_barcode_length: Optional[int] = None
    is_sample_barcode_r1: Optional[bool] = None
    is_sample_barcode_header: Optional[bool] = None
    revcomp_sample_barcode: Optional[bool] = None

    # misc
    retain_inference_results: bool = False
    cores: int = 1

    def to_kwargs(self) -> Dict[str, Any]:
        """Dict of non-None fields suitable for splatting into the legacy entry point."""
        return {f.name: getattr(self, f.name) for f in fields(self) if getattr(self, f.name) is not None}


def map_fastq(
    whitelist_guide_reporter_df: pd.DataFrame,
    fastq_r1_fns: List[str],
    fastq_r2_fns: Optional[List[str]] = None,
    *,
    config: Optional[ParsingConfig] = None,
    **kwargs: Any,
) -> WhitelistReporterCountsResult:
    """Map FASTQs to a whitelist library and return a `WhitelistReporterCountsResult`.

    §4.5 / §7.1: public-API wrapper around
    `mapping.get_whitelist_reporter_counts_from_fastq`. Accepts a
    `ParsingConfig` and/or flat kwargs (flat kwargs override config fields).

    Parameters
    ----------
    whitelist_guide_reporter_df, fastq_r1_fns, fastq_r2_fns
        Same semantics as the legacy entry point.
    config
        Optional `ParsingConfig`. Fields with value `None` are dropped; provide
        any kwargs not in the config separately via `**kwargs`.
    **kwargs
        Flat kwargs; override any corresponding `config` field.

    Returns
    -------
    WhitelistReporterCountsResult
    """
    # Import here to avoid a circular import at module load time.
    from .mapping.main_mapping import get_whitelist_reporter_counts_from_fastq
    merged: Dict[str, Any] = {}
    if config is not None:
        merged.update(config.to_kwargs())
    merged.update(kwargs)
    return get_whitelist_reporter_counts_from_fastq(
        whitelist_guide_reporter_df=whitelist_guide_reporter_df,
        fastq_r1_fns=fastq_r1_fns,
        fastq_r2_fns=fastq_r2_fns,
        **merged,
    )


def count(result: WhitelistReporterCountsResult) -> AllMatchSetWhitelistReporterCounterSeriesResults:
    """Return the per-tier count Series container from a mapping result.

    §4.5 / §7.1: thin accessor that surfaces the already-built counts. In the
    current implementation the mapping call builds counts eagerly; a future
    release can split mapping + counting into separate stages so counts are
    lazily evaluated only when asked for.
    """
    return result.all_match_set_whitelist_reporter_counter_series_results


def alleles(
    result: WhitelistReporterCountsResult,
    tier: str,
    *,
    contains_guide_surrogate: bool,
    contains_guide_barcode: bool,
    contains_guide_umi: bool,
):
    """Build allele count Series for a given mapping tier from a retained mapping result.

    §4.5 / §7.1: wraps `processing.get_matchset_alleleseries`. Requires the
    mapping call to have been run with `retain_inference_results=True` (the
    default is slim); raises `ValueError` with a clear remediation message
    otherwise.

    Parameters
    ----------
    result
        A mapping result from `map_fastq` built with
        `retain_inference_results=True`.
    tier
        A `MatchTier` enum member or its string value.
    contains_guide_surrogate, contains_guide_barcode, contains_guide_umi
        Must match the mapping configuration.

    Returns
    -------
    MatchSetWhitelistReporterObservedSequenceCounterSeriesResults
    """
    from .processing.crispr_editing_processing import get_matchset_alleleseries
    return get_matchset_alleleseries(
        observed_guide_reporter_umi_counts_inferred=result.observed_guide_reporter_umi_counts_inferred,
        attribute_name=str(tier),
        contains_guide_surrogate=contains_guide_surrogate,
        contains_guide_barcode=contains_guide_barcode,
        contains_guide_umi=contains_guide_umi,
    )
