from typeguard import typechecked
from ..models.error_models import GuideCountErrorType, GuideCountError
from ..models.quality_control_models import SingleInferenceQualityControlResult, MatchSetSingleInferenceQualityControlResult, QualityControlResult
from typing import DefaultDict, Tuple, Optional
from collections import defaultdict
from ..processing.crispr_editing_processing import get_non_error_dict

#
# PERFORM THE QC
#
def perform_counts_quality_control(observed_guide_reporter_umi_counts_inferred: DefaultDict[Tuple[str,Optional[str],Optional[str]], dict], contains_umi: bool, contains_surrogate: bool, contains_barcode: bool):
    get_umi_noncollapsed_counts = lambda counts_inferred_dict : sum([sum(counts_inferred_value.observed_value.values()) for counts_inferred_value in counts_inferred_dict.values()]) # NOTE: for UMI, observed_value is a Counter, so need to use .values()
    get_umi_collapsed_counts = lambda counts_inferred_dict : sum([len(counts_inferred_value.observed_value.values()) for counts_inferred_value in counts_inferred_dict.values()]) # NOTE: for UMI, observed_value is a Counter, so need to use .values()
    get_counts = lambda counts_inferred_dict : sum([counts_inferred_value.observed_value for counts_inferred_value in counts_inferred_dict.values()]) # NOTE: for UMI, observed_value is an int, so NO NEED to use .values()

    @typechecked
    def set_num_non_error_counts(single_inference_quality_control_result: SingleInferenceQualityControlResult, counts_inferred_dict: dict) -> SingleInferenceQualityControlResult:
        # TODO: Look into the MatchSetSingleInferenceMatchResultValue.matches, and calculate how many duplicates, no matches, and single matches there are. Are no_matches thrown as error?
        if contains_umi:
            single_inference_quality_control_result.num_non_error_umi_noncollapsed_counts = get_umi_noncollapsed_counts(counts_inferred_dict)
            single_inference_quality_control_result.num_non_error_umi_collapsed_counts = get_umi_collapsed_counts(counts_inferred_dict)

            single_inference_quality_control_result.num_total_umi_noncollapsed_counts = get_umi_noncollapsed_counts(observed_guide_reporter_umi_counts_inferred)
            single_inference_quality_control_result.num_total_umi_collapsed_counts = get_umi_collapsed_counts(observed_guide_reporter_umi_counts_inferred)
        else:
            single_inference_quality_control_result.num_non_error_counts = get_counts(counts_inferred_dict)
            single_inference_quality_control_result.num_total_counts = get_counts(observed_guide_reporter_umi_counts_inferred)
        return single_inference_quality_control_result

    @typechecked
    def set_num_guide_count_error_types(single_inference_quality_control_result: SingleInferenceQualityControlResult, attribute_name: str) -> SingleInferenceQualityControlResult:
        guide_count_error_type_umi_noncollapsed_count: DefaultDict[GuideCountErrorType, int] = defaultdict(int)
        guide_count_error_type_umi_collapsed_count: DefaultDict[GuideCountErrorType, int] = defaultdict(int)
        guide_count_error_type_count: DefaultDict[GuideCountErrorType, int] = defaultdict(int)
        for _, observed_guide_reporter_umi_counts_inferred_value in observed_guide_reporter_umi_counts_inferred.items():
            error_result: Optional[GuideCountError] = getattr(observed_guide_reporter_umi_counts_inferred_value.inferred_value, attribute_name).error if getattr(observed_guide_reporter_umi_counts_inferred_value.inferred_value, attribute_name) is not None else GuideCountError(guide_count_error_type=GuideCountErrorType.NULL_MATCH_RESULT)
            if error_result is not None:
                if contains_umi:
                    guide_count_error_type_umi_noncollapsed_count[error_result.guide_count_error_type] += sum(observed_guide_reporter_umi_counts_inferred_value.observed_value.values())
                    guide_count_error_type_umi_collapsed_count[error_result.guide_count_error_type] += len(observed_guide_reporter_umi_counts_inferred_value.observed_value.values())
                else:
                    guide_count_error_type_count[error_result.guide_count_error_type] += observed_guide_reporter_umi_counts_inferred_value.observed_value

        single_inference_quality_control_result.guide_count_error_type_umi_noncollapsed_count = guide_count_error_type_umi_noncollapsed_count
        single_inference_quality_control_result.guide_count_error_type_umi_collapsed_count = guide_count_error_type_umi_collapsed_count
        single_inference_quality_control_result.guide_count_error_type_count = guide_count_error_type_count

        return single_inference_quality_control_result

    @typechecked
    def set_match_set_single_inference_quality_control_results(attribute_name: str) -> MatchSetSingleInferenceQualityControlResult:
        single_inference_quality_control_result = MatchSetSingleInferenceQualityControlResult()
        single_inference_quality_control_result.non_error_dict = get_non_error_dict(observed_guide_reporter_umi_counts_inferred, attribute_name)
        
        single_inference_quality_control_result = set_num_non_error_counts(single_inference_quality_control_result, single_inference_quality_control_result.non_error_dict)
        single_inference_quality_control_result = set_num_guide_count_error_types(single_inference_quality_control_result, attribute_name)
        return single_inference_quality_control_result


    quality_control_result = QualityControlResult()
    quality_control_result.protospacer_match = set_match_set_single_inference_quality_control_results("protospacer_match")
    if contains_barcode:
        quality_control_result.protospacer_match_barcode_match = set_match_set_single_inference_quality_control_results("protospacer_match_barcode_match")
        if contains_surrogate:
            quality_control_result.protospacer_match_surrogate_match_barcode_match = set_match_set_single_inference_quality_control_results("protospacer_match_surrogate_match_barcode_match")
            protospacer_mismatch_surrogate_match_barcode_match_nonerror_dict = get_non_error_dict("protospacer_mismatch_surrogate_match_barcode_match")
            quality_control_result.protospacer_mismatch_surrogate_match_barcode_match = set_match_set_single_inference_quality_control_results("protospacer_mismatch_surrogate_match_barcode_match")
    if contains_surrogate:
        quality_control_result.protospacer_match_surrogate_match = set_match_set_single_inference_quality_control_results("protospacer_match_surrogate_match")
        quality_control_result.protospacer_mismatch_surrogate_match = set_match_set_single_inference_quality_control_results("protospacer_mismatch_surrogate_match")
    
    return quality_control_result