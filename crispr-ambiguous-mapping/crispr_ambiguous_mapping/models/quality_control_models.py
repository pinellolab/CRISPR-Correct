from dataclasses import dataclass
from typing import Optional, DefaultDict
from .error_models import GuideCountError

@dataclass
class SingleInferenceQualityControlResult:
    non_error_dict: Optional[DefaultDict] = None
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