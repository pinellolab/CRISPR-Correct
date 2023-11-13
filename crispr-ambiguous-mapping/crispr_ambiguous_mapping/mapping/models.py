from dataclasses import dataclass
from typing import Union, List, Mapping, Tuple, Optional, Any, DefaultDict
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

@dataclass
class UnequalLengthGuideCountError(GuideCountError):
    # NOTE: Added default argument to prevent error of having non-default after defaults in child classes
    sequence_length: Optional[int]  = None
    expected_length: Optional[int] = None

@dataclass
class ProtospacerUnequalLengthGuideCountError(UnequalLengthGuideCountError):
    guide_count_error_type: GuideCountErrorType = GuideCountErrorType.NO_PROTOSPACER_WITH_SAME_LENGTH

@dataclass
class SurrogateUnequalLengthGuideCountError(UnequalLengthGuideCountError):
    guide_count_error_type: GuideCountErrorType = GuideCountErrorType.NO_SURROGATE_WITH_SAME_LENGTH

@dataclass
class BarcodeUnequalLengthGuideCountError(UnequalLengthGuideCountError):
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

@dataclass
class ProtospacerHammingThresholdGuideCountError(HammingThresholdGuideCountError):
    guide_count_error_type: GuideCountErrorType = GuideCountErrorType.NO_MATCH_PROTOSPACER_HAMMING_THRESHOLD

@dataclass
class SurrogateHammingThresholdGuideCountError(HammingThresholdGuideCountError):
    guide_count_error_type: GuideCountErrorType = GuideCountErrorType.NO_MATCH_SURROGATE_HAMMING_THRESHOLD

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


@dataclass
class AllMatchSetWhitelistReporterCounterSeriesResults:
    protospacer_match: Optional[MatchSetWhitelistReporterCounterSeriesResults] = None
    protospacer_match_surrogate_match_barcode_match: Optional[MatchSetWhitelistReporterCounterSeriesResults] = None
    protospacer_match_surrogate_match: Optional[MatchSetWhitelistReporterCounterSeriesResults] = None
    protospacer_match_barcode_match: Optional[MatchSetWhitelistReporterCounterSeriesResults] = None
    #protospacer_mismatch_surrogate_match_barcode_match: Optional[SurrogateProtospacerMismatchSingleInferenceQualityControlResult] = None
    #protospacer_mismatch_surrogate_match: Optional[SurrogateProtospacerMismatchSingleInferenceQualityControlResult] = None
