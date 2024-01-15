from dataclasses import dataclass
from typing import Union, List, Mapping, Tuple, Optional, Any, DefaultDict, Dict
from typing import Counter as CounterType
from enum import Enum
from collections import Counter, defaultdict

import pandas as pd

class GuideCountErrorType(Enum):
    NO_MATCH = "No guide found within hamming distance"
    MULTIPLE_MATCH = "Multiple exact matches found for guide (likely a truncated guide read assuming guide series is unique)"
    NO_PROTOSPACER_WITH_SAME_LENGTH = "No whitelisted guides with same length as observed guide - maybe try enabling truncating whitelisted guides"
    NO_BARCODE_WITH_SAME_LENGTH = "No whitelisted barcode with the same length as the observed"
    NO_SURROGATE_WITH_SAME_LENGTH = "No whitelisted surrogate with the same length as the observed"
    NO_MATCH_PROTOSPACER_HAMMING_THRESHOLD = "No whitelisted guide that is below the provided guide hamming distance threshold"
    NO_MATCH_BARCODE_HAMMING_THRESHOLD = "No whitelisted barcode that is below the provided barcode hamming distance threshold"
    NO_MATCH_SURROGATE_HAMMING_THRESHOLD = "The inferred whitelisted surrogate does not match with the observed surrogate below the given hamming distance threshold"
    MULTIPLE_MATCH_EXACT = "Multiple exact guide matches, double check that there are no duplicates in your guide library (especially if truncation is enabled)"
    NO_PROTOSPACER_MATCH_MISSING_INFO = "The protospacer was not provided"
    NO_SURROGATE_MATCH_MISSING_INFO = "The surrogate was not provided"
    NO_BARCODE_MATCH_MISSING_INFO = "The barcode was not provided"
    NO_MATCH_OBSERVED_SURROGATE_SHORTER = "The observed surrogate is shorter than the inferred surrogate"
    NULL_MATCH_RESULT = "Match result was null, likely since match result type is not supported given the inputs (i.e. a barcode_match result when barcode is not specified)"

@dataclass
class GuideCountError:
    guide_count_error_type: GuideCountErrorType
    miscellaneous_info_dict: dict = None

@dataclass
class InsufficientLengthGuideCountError(GuideCountError):
    # NOTE: Added default argument to prevent error of having non-default after defaults in child classes
    sequence_length: Optional[int]  = None
    minimum_length: Optional[int] = None

@dataclass
class ProtospacerInsufficientLengthGuideCountError(InsufficientLengthGuideCountError):
    guide_count_error_type: GuideCountErrorType = GuideCountErrorType.NO_PROTOSPACER_WITH_SAME_LENGTH

@dataclass
class SurrogateInsufficientLengthGuideCountError(InsufficientLengthGuideCountError):
    guide_count_error_type: GuideCountErrorType = GuideCountErrorType.NO_SURROGATE_WITH_SAME_LENGTH

@dataclass
class BarcodeInsufficientLengthGuideCountError(InsufficientLengthGuideCountError):
    guide_count_error_type: GuideCountErrorType = GuideCountErrorType.NO_BARCODE_WITH_SAME_LENGTH

@dataclass
class MissingInfoGuideCountError(GuideCountError):
    sequence_value: Optional[Any] = None

@dataclass
class ProtospacerMissingInfoGuideCountError(MissingInfoGuideCountError):
    guide_count_error_type: GuideCountErrorType = GuideCountErrorType.NO_PROTOSPACER_MATCH_MISSING_INFO

@dataclass
class SurrogateMissingInfoGuideCountError(MissingInfoGuideCountError):
    guide_count_error_type: GuideCountErrorType = GuideCountErrorType.NO_SURROGATE_MATCH_MISSING_INFO

@dataclass
class BarcodeMissingInfoGuideCountError(MissingInfoGuideCountError):
    guide_count_error_type: GuideCountErrorType = GuideCountErrorType.NO_BARCODE_MATCH_MISSING_INFO


@dataclass
class HammingThresholdGuideCountError(GuideCountError):
    hamming_min: Optional[int] = None
    hamming_threshold: Optional[int] = None
    original_df: Optional[pd.DataFrame] = None
    hamming_min_match_df: Optional[pd.DataFrame] = None

@dataclass
class ProtospacerHammingThresholdGuideCountError(HammingThresholdGuideCountError):
    guide_count_error_type: GuideCountErrorType = GuideCountErrorType.NO_MATCH_PROTOSPACER_HAMMING_THRESHOLD
    barcode_subsetted: Optional[bool] = None

@dataclass
class SurrogateHammingThresholdGuideCountError(HammingThresholdGuideCountError):
    guide_count_error_type: GuideCountErrorType = GuideCountErrorType.NO_MATCH_SURROGATE_HAMMING_THRESHOLD
    barcode_subsetted: Optional[bool] = None

@dataclass
class BarcodeHammingThresholdGuideCountError(HammingThresholdGuideCountError):
    guide_count_error_type: GuideCountErrorType = GuideCountErrorType.NO_MATCH_BARCODE_HAMMING_THRESHOLD


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
class SingleInferenceQualityControlResult:
    non_error_dict: Optional[defaultdict] = None
    num_non_error_umi_noncollapsed_counts: Optional[int] = None
    num_non_error_umi_collapsed_counts: Optional[int] = None
    num_non_error_counts: Optional[int] = None
    num_total_umi_noncollapsed_counts: Optional[int] = None
    num_total_umi_collapsed_counts: Optional[int] = None
    num_total_counts: Optional[int] = None
    guide_count_error_type_umi_noncollapsed_count: DefaultDict[GuideCountError, int] = None
    guide_count_error_type_umi_collapsed_count: DefaultDict[GuideCountError, int] = None
    guide_count_error_type_count: DefaultDict[GuideCountError, int] = None

@dataclass
class MatchSetSingleInferenceQualityControlResult(SingleInferenceQualityControlResult):
    # TODO: In the future, can add any QC that is specific to the singe match set result type
    pass

@dataclass
class SurrogateProtospacerMismatchSingleInferenceQualityControlResult(SingleInferenceQualityControlResult):
    # TODO: In the future, can add any QC that is specific to the surrogate_protospacer mismatch result type
    pass

@dataclass
class QualityControlResult:
    protospacer_match: Optional[MatchSetSingleInferenceQualityControlResult] = None
    protospacer_match_surrogate_match_barcode_match: Optional[MatchSetSingleInferenceQualityControlResult] = None
    protospacer_match_surrogate_match: Optional[MatchSetSingleInferenceQualityControlResult] = None
    protospacer_match_barcode_match: Optional[MatchSetSingleInferenceQualityControlResult] = None
    protospacer_mismatch_surrogate_match_barcode_match: Optional[SurrogateProtospacerMismatchSingleInferenceQualityControlResult] = None
    protospacer_mismatch_surrogate_match: Optional[SurrogateProtospacerMismatchSingleInferenceQualityControlResult] = None


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
class MatchSetWhitelistReporterObservedSequenceCounterSeriesResults:
    
    # Storing as a dictionary 
    ambiguous_ignored_umi_noncollapsed_alleleseries_dict : Optional[DefaultDict[Tuple[str, Optional[str], Optional[str]], pd.Series]] = None
    ambiguous_ignored_umi_collapsed_alleleseries_dict : Optional[DefaultDict[Tuple[str, Optional[str], Optional[str]], pd.Series]] = None
    ambiguous_ignored_alleleseries_dict : Optional[DefaultDict[Tuple[str, Optional[str], Optional[str]], pd.Series]] = None

    ambiguous_accepted_umi_noncollapsed_alleleseries_dict : Optional[DefaultDict[Tuple[str, Optional[str], Optional[str]], pd.Series]] = None
    ambiguous_accepted_umi_collapsed_alleleseries_dict : Optional[DefaultDict[Tuple[str, Optional[str], Optional[str]], pd.Series]] = None
    ambiguous_accepted_alleleseries_dict : Optional[DefaultDict[Tuple[str, Optional[str], Optional[str]], pd.Series]] = None

    ambiguous_spread_umi_noncollapsed_alleleseries_dict : Optional[DefaultDict[Tuple[str, Optional[str], Optional[str]], pd.Series]] = None
    ambiguous_spread_umi_collapsed_alleleseries_dict : Optional[DefaultDict[Tuple[str, Optional[str], Optional[str]], pd.Series]] = None
    ambiguous_spread_alleleseries_dict : Optional[DefaultDict[Tuple[str, Optional[str], Optional[str]], pd.Series]] = None
        
    # Storing as a dataframe
    ambiguous_ignored_umi_noncollapsed_allele_df : pd.DataFrame = None
    ambiguous_ignored_umi_collapsed_allele_df : pd.DataFrame = None
    ambiguous_ignored_allele_df : pd.DataFrame = None

    ambiguous_accepted_umi_noncollapsed_allele_df : pd.DataFrame = None
    ambiguous_accepted_umi_collapsed_allele_df : pd.DataFrame = None
    ambiguous_accepted_allele_df : pd.DataFrame = None

    ambiguous_spread_umi_noncollapsed_allele_df : pd.DataFrame = None
    ambiguous_spread_umi_collapsed_allele_df : pd.DataFrame = None
    ambiguous_spread_allele_df : pd.DataFrame = None

    # This ensures that any empty series are kept at None
    def __setattr__(self, name, value):
        if isinstance(value, dict):
            if len(value) == 0: # If no items in dictionary, just return None
                super().__setattr__(name, None)
            else:
                super().__setattr__(name, value)
        elif isinstance(value, pd.DataFrame):
            if value.empty(): # If no items in dataframe, just return None
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
    observed_guide_reporter_umi_counts_inferred: DefaultDict[Tuple[str,Optional[str],Optional[str]], dict]
    quality_control_result: QualityControlResult
    count_input: CountInput

@dataclass
class ObservedSequenceMutations:
    linked_mutations_whitelist_reporter_dict: Dict[Tuple[str, Optional[str], Optional[str]], pd.DataFrame] = None
    all_observed_protospacer_unlinked_mutations_df: pd.DataFrame = None
    all_observed_surrogate_unlinked_mutations_df: Optional[pd.DataFrame] = None
    all_observed_barcode_unlinked_mutations_df: Optional[pd.DataFrame] = None
    
@dataclass
class MatchSetWhitelistReporterObservedSequenceMutationsResults:
    ambiguous_ignored_umi_noncollapsed_mutations : Optional[ObservedSequenceMutations] = None
    ambiguous_ignored_umi_collapsed_mutations : Optional[ObservedSequenceMutations] = None
    ambiguous_ignored_mutations : Optional[ObservedSequenceMutations] = None

    ambiguous_accepted_umi_noncollapsed_mutations : Optional[ObservedSequenceMutations] = None
    ambiguous_accepted_umi_collapsed_mutations : Optional[ObservedSequenceMutations] = None
    ambiguous_accepted_mutations : Optional[ObservedSequenceMutations] = None

    ambiguous_spread_umi_noncollapsed_mutations : Optional[ObservedSequenceMutations] = None
    ambiguous_spread_umi_collapsed_mutations : Optional[ObservedSequenceMutations] = None
    ambiguous_spread_mutations : Optional[ObservedSequenceMutations] = None