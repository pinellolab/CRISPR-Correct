
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
    BARCODE_MISC = "Miscellaneous barcode error"
    PROTOSPACER_MISC = "Miscellaneous protospacer error"
    SURROGATE_MISC = "Miscellaneous surrogate error"

_DEFAULT_SUGGESTIONS = {
    GuideCountErrorType.NO_MATCH_PROTOSPACER_HAMMING_THRESHOLD: (
        "Observed protospacer did not match any whitelist guide within the strict Hamming threshold. "
        "Consider relaxing `protospacer_hamming_threshold_strict` (e.g. 7 -> 10) if your read quality is low, "
        "or verify the protospacer extraction parameters (position/flanks/length, revcomp, r1 vs r2)."
    ),
    GuideCountErrorType.NO_MATCH_SURROGATE_HAMMING_THRESHOLD: (
        "The inferred surrogate did not match within `surrogate_hamming_threshold_strict`. "
        "Base-editing drives surrogate Hamming distances higher than protospacer distances — "
        "threshold ~10 for 32bp surrogates is typical. If you're seeing many of these, raise the threshold."
    ),
    GuideCountErrorType.NO_MATCH_BARCODE_HAMMING_THRESHOLD: (
        "Observed barcode did not match any whitelist barcode within `guide_barcode_hamming_threshold_strict`. "
        "4bp barcodes typically use threshold 2. Check the barcode extraction regex and revcomp setting."
    ),
    GuideCountErrorType.NO_PROTOSPACER_WITH_SAME_LENGTH: (
        "The observed protospacer is shorter than the whitelist minimum. Check extraction length/position."
    ),
    GuideCountErrorType.NO_BARCODE_WITH_SAME_LENGTH: (
        "The observed barcode is shorter than the whitelist minimum. Check extraction length/position."
    ),
    GuideCountErrorType.NO_SURROGATE_WITH_SAME_LENGTH: (
        "The observed surrogate is shorter than the whitelist minimum. Check extraction length/position."
    ),
    GuideCountErrorType.NO_PROTOSPACER_MATCH_MISSING_INFO: (
        "Protospacer extraction returned None — the regex/flank did not match the read. "
        "Inspect a few reads vs. your `protospacer_pattern_regex` / flanks."
    ),
    GuideCountErrorType.NO_SURROGATE_MATCH_MISSING_INFO: (
        "Surrogate extraction returned None. Inspect reads vs. your surrogate pattern."
    ),
    GuideCountErrorType.NO_BARCODE_MATCH_MISSING_INFO: (
        "Barcode extraction returned None. Inspect reads vs. your barcode pattern."
    ),
    GuideCountErrorType.MULTIPLE_MATCH_EXACT: (
        "Multiple whitelist rows share an exact sequence — deduplicate the whitelist or disable truncation."
    ),
}


@dataclass
class GuideCountError:
    guide_count_error_type: GuideCountErrorType
    miscellaneous_info_dict: dict = None

    @property
    def suggestion(self) -> str:
        """User-facing remediation hint for this error type.

        §4.10: errors now expose a `suggestion` property instead of printing a
        generic enum description — pair error-type filtering downstream with a
        human-readable next step.
        """
        return _DEFAULT_SUGGESTIONS.get(self.guide_count_error_type, "")

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
    # PERF/MEM: previously carried `original_df` and `hamming_min_match_df`
    # pandas DataFrames attached to every error record — roughly 150 KB per
    # instance for a 2 k-guide whitelist. For 45 M-read runs with a 30 %
    # error rate this single field accounted for O(10 GB) of peak RSS.
    # Replaced with lightweight counts; reconstruct the DataFrame at report
    # time if ever needed.
    n_whitelist_candidates: Optional[int] = None
    n_hamming_min_match: Optional[int] = None

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