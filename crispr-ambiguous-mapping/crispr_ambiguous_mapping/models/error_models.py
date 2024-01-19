
from enum import Enum
from typing import Optional, Any
from dataclasses import dataclass
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