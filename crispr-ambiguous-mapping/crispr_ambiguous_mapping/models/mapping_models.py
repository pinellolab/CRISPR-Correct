from dataclasses import dataclass
from typing import Union, List, Mapping, Tuple, Optional, Any, DefaultDict, Dict
from typing import Counter as CounterType

from collections import Counter, defaultdict
import pandas as pd

from .error_models import GuideCountError
from .quality_control_models import QualityControlResult
from .types import * 

@dataclass
class SingleInferenceMatchResultValue:
    pass

@dataclass
class MatchSetSingleInferenceMatchResultValue(SingleInferenceMatchResultValue):
    matches: Optional[pd.DataFrame] = None

@dataclass
class SurrogateProtospacerMismatchSingleInferenceMatchResultValue(SingleInferenceMatchResultValue):
    mismatched: bool
    surrogate_matches: pd.DataFrame
    protospacer_matches: pd.DataFrame 
    protospacer_surrogate_matches: pd.DataFrame

@dataclass
class SingleInferenceMatchResult:
    error: Optional[GuideCountError] = None
    value: Optional[SingleInferenceMatchResultValue] = None

@dataclass
class MatchSetSingleInferenceMatchResult(SingleInferenceMatchResult):
    value: Optional[MatchSetSingleInferenceMatchResultValue] = None

@dataclass
class SurrogateProtospacerMismatchSingleInferenceMatchResult(SingleInferenceMatchResult):
    value: Optional[SurrogateProtospacerMismatchSingleInferenceMatchResultValue] = None

@dataclass
class CompleteInferenceMatchResult:
    protospacer_match: Optional[MatchSetSingleInferenceMatchResult] = None
    protospacer_match_surrogate_match_barcode_match: Optional[MatchSetSingleInferenceMatchResult] = None
    protospacer_match_surrogate_match: Optional[MatchSetSingleInferenceMatchResult] = None
    protospacer_match_barcode_match: Optional[MatchSetSingleInferenceMatchResult] =None
    protospacer_mismatch_surrogate_match_barcode_match: Optional[SurrogateProtospacerMismatchSingleInferenceMatchResult] = None
    protospacer_mismatch_surrogate_match: Optional[SurrogateProtospacerMismatchSingleInferenceMatchResult] = None

@dataclass
class InferenceResult:
    observed_value: Union[int, CounterType[Optional[str]]]
    inferred_value: CompleteInferenceMatchResult


@dataclass
class WhitelistReporterCountsResult:
    observed_guide_reporter_umi_counts_inferred: DefaultDict[Tuple[str,Optional[str],Optional[str]], dict]
    quality_control_result: QualityControlResult


@dataclass
class MatchSetWhitelistReporterCounterSeriesResults:
    ambiguous_ignored_umi_noncollapsed_counterseries : Optional[pd.Series] = None
    ambiguous_ignored_umi_collapsed_counterseries : Optional[pd.Series] = None
    ambiguous_ignored_counterseries : Optional[pd.Series] = None

    ambiguous_accepted_umi_noncollapsed_counterseries : Optional[pd.Series] = None
    ambiguous_accepted_umi_collapsed_counterseries : Optional[pd.Series] = None
    ambiguous_accepted_counterseries : Optional[pd.Series] = None

    ambiguous_spread_umi_noncollapsed_counterseries : Optional[pd.Series] = None
    ambiguous_spread_umi_collapsed_counterseries : Optional[pd.Series] = None
    ambiguous_spread_counterseries : Optional[pd.Series] = None

    # This ensures that any empty series are kept at None
    def __setattr__(self, name, value):
        if isinstance(value, pd.Series):
            if value.sum() == 0: # If sum of pandas series is 0, just return None
                super().__setattr__(name, None)
            else:
                super().__setattr__(name, value)
        else:
            # Set the attribute
            super().__setattr__(name, value)

@dataclass
class SurrogateProtospacerMismatchSetWhitelistReporterCounterSeriesResults:
    
    # MATCH counters
    ambiguous_ignored_umi_noncollapsed_match_counterseries : Optional[pd.Series]  = None
    ambiguous_ignored_umi_collapsed_match_counterseries : Optional[pd.Series]  = None
    ambiguous_ignored_match_counterseries : Optional[pd.Series] = None

    ambiguous_accepted_umi_noncollapsed_match_counterseries : Optional[pd.Series] = None
    ambiguous_accepted_umi_collapsed_match_counterseries : Optional[pd.Series] = None
    ambiguous_accepted_match_counterseries : Optional[pd.Series] = None

    ambiguous_spread_umi_noncollapsed_match_counterseries : Optional[pd.Series] = None
    ambiguous_spread_umi_collapsed_match_counterseries : Optional[pd.Series] = None
    ambiguous_spread_match_counterseries : Optional[pd.Series] = None

    # MISMATCH counters (keys are PAIRS of indices, representing the protospacer and surrogate match separately), [str, Optional[str], Optional[str], str, Optional[str], Optional[str]], first three are Protospacer_Match, last three is Surrogate_Match
    ambiguous_ignored_umi_noncollapsed_mismatch_counterseries : Optional[pd.Series] = None
    ambiguous_ignored_umi_collapsed_mismatch_counterseries : Optional[pd.Series] = None
    ambiguous_ignored_mismatch_counterseries : Optional[pd.Series] = None

    ambiguous_accepted_umi_noncollapsed_mismatch_counterseries : Optional[pd.Series] = None
    ambiguous_accepted_umi_collapsed_mismatch_counterseries : Optional[pd.Series] = None
    ambiguous_accepted_mismatch_counterseries : Optional[pd.Series] = None

    ambiguous_spread_umi_noncollapsed_mismatch_counterseries : Optional[pd.Series] = None
    ambiguous_spread_umi_collapsed_mismatch_counterseries : Optional[pd.Series] = None
    ambiguous_spread_mismatch_counterseries : Optional[pd.Series] = None

    # This ensures that any empty series are kept at None
    def __setattr__(self, name, value):
        if isinstance(value, pd.Series):
            if value.sum() == 0: # If sum of pandas series is 0, just return None
                super().__setattr__(name, None)
            else:
                super().__setattr__(name, value)
        else:
            # Set the attribute
            super().__setattr__(name, value)

@dataclass
class AllMatchSetWhitelistReporterCounterSeriesResults:
    protospacer_match: Optional[MatchSetWhitelistReporterCounterSeriesResults] = None
    protospacer_match_surrogate_match_barcode_match: Optional[MatchSetWhitelistReporterCounterSeriesResults] = None
    protospacer_match_surrogate_match: Optional[MatchSetWhitelistReporterCounterSeriesResults] = None
    protospacer_match_barcode_match: Optional[MatchSetWhitelistReporterCounterSeriesResults] = None
    protospacer_mismatch_surrogate_match_barcode_match: Optional[SurrogateProtospacerMismatchSetWhitelistReporterCounterSeriesResults] = None
    protospacer_mismatch_surrogate_match: Optional[SurrogateProtospacerMismatchSetWhitelistReporterCounterSeriesResults] = None

@dataclass
class CountInput:
    whitelist_guide_reporter_df: pd.DataFrame
    contains_surrogate:bool
    contains_barcode:bool
    contains_umi:bool
    protospacer_hamming_threshold_strict: Optional[int]
    surrogate_hamming_threshold_strict: Optional[int]
    barcode_hamming_threshold_strict: Optional[int]

@dataclass
class WhitelistReporterCountsResult:
    all_match_set_whitelist_reporter_counter_series_results: AllMatchSetWhitelistReporterCounterSeriesResults
    observed_guide_reporter_umi_counts_inferred: GeneralMappingInferenceDict
    quality_control_result: QualityControlResult
    count_input: CountInput
