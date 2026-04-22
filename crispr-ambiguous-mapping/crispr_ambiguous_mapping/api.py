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
from pathlib import Path
from typing import Any, Dict, List, Optional
import json
import pandas as pd

from .models.mapping_models import (
    WhitelistReporterCountsResult,
    AllMatchSetWhitelistReporterCounterSeriesResults,
    MatchSetWhitelistReporterCounterSeriesResults,
    SurrogateProtospacerMismatchSetWhitelistReporterCounterSeriesResults,
    CountInput,
    MatchTier,
)
from .models.quality_control_models import (
    QualityControlResult,
    MatchSetSingleInferenceQualityControlResult,
    SurrogateProtospacerMismatchSingleInferenceQualityControlResult,
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


# ---------------------------------------------------------------------------
# §7.4 / §7.1: parquet-backed save/load (cross-language durable serialization)
# ---------------------------------------------------------------------------

_MATCH_STRATEGY_ATTRS = [
    ("ambiguous_ignored", ""),
    ("ambiguous_ignored", "umi_collapsed"),
    ("ambiguous_ignored", "umi_noncollapsed"),
    ("ambiguous_accepted", ""),
    ("ambiguous_accepted", "umi_collapsed"),
    ("ambiguous_accepted", "umi_noncollapsed"),
    ("ambiguous_spread", ""),
    ("ambiguous_spread", "umi_collapsed"),
    ("ambiguous_spread", "umi_noncollapsed"),
]

_MATCH_TIERS = ["protospacer_match", "protospacer_match_surrogate_match",
                "protospacer_match_barcode_match", "protospacer_match_surrogate_match_barcode_match"]
_MISMATCH_TIERS = ["protospacer_mismatch_surrogate_match", "protospacer_mismatch_surrogate_match_barcode_match"]


def _strategy_attr(amb: str, umi: str, kind: str = "") -> str:
    """Resolve (ambiguity, umi_strategy, kind) to the tier-object attribute name.

    kind in {'', 'match', 'mismatch'} — empty for the non-mismatch tiers;
    'match' / 'mismatch' for SurrogateProtospacerMismatch tiers.
    """
    pieces = [amb]
    if umi:
        pieces.append(umi)
    if kind:
        pieces.append(kind)
    pieces.append("counterseries")
    return "_".join(pieces)


def _tier_to_dataframe(tier_obj, kind: str = "") -> pd.DataFrame:
    """Convert a tier dataclass (9 Series) into a single DataFrame (strategies = columns)."""
    series_by_col: Dict[str, pd.Series] = {}
    for amb, umi in _MATCH_STRATEGY_ATTRS:
        attr = _strategy_attr(amb, umi, kind)
        s = getattr(tier_obj, attr, None)
        if s is None:
            continue
        col = f"{amb}" + (f"__{umi}" if umi else "")
        series_by_col[col] = s
    if not series_by_col:
        return pd.DataFrame()
    # pandas aligns on a union of indices automatically.
    return pd.DataFrame(series_by_col)


def _column_to_series(col_values: pd.Series) -> pd.Series:
    """Normalize one DataFrame column back to a Series the package expects:
    drop NaN (from the alignment union), downcast to int64 when every value
    is whole (parquet widens to float64 when any NaN is present), strip the
    column name so callers aren't suprised by dtype / name drift."""
    s = col_values.dropna()
    if s.empty:
        return s
    # If all values are whole, restore int64 — matches the pre-save dtype of
    # the non-spread counters. Ambiguous-spread can produce fractional
    # counts (/N_matches) so we keep those as float.
    vals = s.values
    try:
        if (vals == vals.astype("int64")).all():
            s = s.astype("int64")
    except (TypeError, ValueError):
        pass
    s.name = None
    return s


def _dataframe_to_tier(df: pd.DataFrame, tier_class, kind: str = ""):
    """Reverse of `_tier_to_dataframe`: populate a tier dataclass from a DataFrame."""
    obj = tier_class()
    if df.empty:
        return obj
    for col in df.columns:
        amb, _, umi = col.partition("__")
        attr = _strategy_attr(amb, umi, kind)
        # The __setattr__ guard on the tier class nulls all-zero Series.
        setattr(obj, attr, _column_to_series(df[col]))
    return obj


def _qc_tier_to_dict(qc_tier) -> Dict[str, Any]:
    """Serialize one QC tier to a JSON-friendly dict. Enum error keys → .name."""
    out: Dict[str, Any] = {}
    for f in fields(qc_tier):
        v = getattr(qc_tier, f.name)
        if v is None:
            out[f.name] = None
        elif f.name.startswith("guide_count_error_type") and hasattr(v, "items"):
            out[f.name] = {str(k): int(cnt) for k, cnt in v.items()}
        else:
            out[f.name] = v
    return out


def save(result: WhitelistReporterCountsResult, directory: str | Path, *, overwrite: bool = False) -> Path:
    """Write a mapping result to a directory as parquet + JSON.

    §7.4 / §7.1: cross-language durable serialization. The directory ends up with:

        result_dir/
            manifest.json           # schema version + tier list + timestamp
            counts_<tier>.parquet   # match tiers: 1 file; mismatch tiers: _match + _mismatch
            qc.json                 # QualityControlResult summary counts
            count_input.json        # echo of parsing flags (whitelist DF not round-tripped)

    Complement to pickling — does NOT include the per-observation inference
    dict. If you need that, keep a pickle alongside.

    Returns the directory path.
    """
    directory = Path(directory)
    if directory.exists():
        if not overwrite and any(directory.iterdir()):
            raise FileExistsError(f"{directory} is not empty; pass overwrite=True to replace.")
    directory.mkdir(parents=True, exist_ok=True)

    allw = result.all_match_set_whitelist_reporter_counter_series_results
    manifest = {
        "schema_version": "1",
        "match_tiers": [],
        "mismatch_tiers": [],
    }

    # Match tiers.
    for tier_name in _MATCH_TIERS:
        t = getattr(allw, tier_name, None)
        if t is None:
            continue
        df = _tier_to_dataframe(t)
        if df.empty:
            continue
        path = directory / f"counts_{tier_name}.parquet"
        df.to_parquet(path)
        manifest["match_tiers"].append(tier_name)

    # Mismatch tiers — two DataFrames per tier.
    for tier_name in _MISMATCH_TIERS:
        t = getattr(allw, tier_name, None)
        if t is None:
            continue
        wrote_any = False
        for kind in ("match", "mismatch"):
            df = _tier_to_dataframe(t, kind=kind)
            if df.empty:
                continue
            path = directory / f"counts_{tier_name}_{kind}.parquet"
            df.to_parquet(path)
            wrote_any = True
        if wrote_any:
            manifest["mismatch_tiers"].append(tier_name)

    # QC.
    qc_obj = result.quality_control_result
    qc_out = {}
    for f in fields(qc_obj):
        tier_qc = getattr(qc_obj, f.name)
        qc_out[f.name] = _qc_tier_to_dict(tier_qc) if tier_qc is not None else None
    (directory / "qc.json").write_text(json.dumps(qc_out, indent=2, default=str))

    # CountInput (drop the DataFrame — user can reload the library separately).
    ci = result.count_input
    ci_out = {f.name: getattr(ci, f.name) for f in fields(ci) if f.name != "whitelist_guide_reporter_df"}
    (directory / "count_input.json").write_text(json.dumps(ci_out, indent=2, default=str))

    (directory / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return directory


def load(directory: str | Path) -> WhitelistReporterCountsResult:
    """Reconstruct a mapping result from a directory written by `save`.

    §7.4 / §7.1: pairs with `save`. Reconstructs `all_match_set_whitelist_reporter_counter_series_results`,
    `quality_control_result` summary counts, and `count_input`. The per-
    observation inference dict (`observed_guide_reporter_umi_counts_inferred`)
    is left as `None` — it is not round-tripped (use pickle for that).
    """
    directory = Path(directory)
    manifest_path = directory / "manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(f"No manifest.json in {directory}; not a mapping-result directory.")
    manifest = json.loads(manifest_path.read_text())

    allw = AllMatchSetWhitelistReporterCounterSeriesResults()
    for tier_name in manifest.get("match_tiers", []):
        path = directory / f"counts_{tier_name}.parquet"
        df = pd.read_parquet(path)
        setattr(allw, tier_name, _dataframe_to_tier(df, MatchSetWhitelistReporterCounterSeriesResults))
    for tier_name in manifest.get("mismatch_tiers", []):
        obj = SurrogateProtospacerMismatchSetWhitelistReporterCounterSeriesResults()
        for kind in ("match", "mismatch"):
            path = directory / f"counts_{tier_name}_{kind}.parquet"
            if path.exists():
                df = pd.read_parquet(path)
                for col in df.columns:
                    amb, _, umi = col.partition("__")
                    attr = _strategy_attr(amb, umi, kind)
                    setattr(obj, attr, _column_to_series(df[col]))
        setattr(allw, tier_name, obj)

    # QC — reconstruct the tier objects from the dict.
    qc_raw = json.loads((directory / "qc.json").read_text())
    qc_obj = QualityControlResult()
    for f in fields(qc_obj):
        raw = qc_raw.get(f.name)
        if raw is None:
            continue
        cls = SurrogateProtospacerMismatchSingleInferenceQualityControlResult if f.name.startswith("protospacer_mismatch") else MatchSetSingleInferenceQualityControlResult
        inst = cls()
        for sub_f in fields(inst):
            setattr(inst, sub_f.name, raw.get(sub_f.name))
        setattr(qc_obj, f.name, inst)

    # CountInput (whitelist DF intentionally None — caller can reload separately).
    ci_raw = json.loads((directory / "count_input.json").read_text())
    count_input = CountInput(
        whitelist_guide_reporter_df=pd.DataFrame(),  # stub; user reloads if needed
        contains_guide_surrogate=ci_raw.get("contains_guide_surrogate", False),
        contains_guide_barcode=ci_raw.get("contains_guide_barcode", False),
        contains_guide_umi=ci_raw.get("contains_guide_umi", False),
        contains_sample_barcode=ci_raw.get("contains_sample_barcode", False),
        protospacer_hamming_threshold_strict=ci_raw.get("protospacer_hamming_threshold_strict"),
        surrogate_hamming_threshold_strict=ci_raw.get("surrogate_hamming_threshold_strict"),
        guide_barcode_hamming_threshold_strict=ci_raw.get("guide_barcode_hamming_threshold_strict"),
    )

    return WhitelistReporterCountsResult(
        all_match_set_whitelist_reporter_counter_series_results=allw,
        observed_guide_reporter_umi_counts_inferred=None,
        quality_control_result=qc_obj,
        count_input=count_input,
    )
