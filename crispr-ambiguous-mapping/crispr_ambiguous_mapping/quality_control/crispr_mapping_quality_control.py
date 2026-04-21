from typeguard import typechecked
from ..models.error_models import GuideCountErrorType, GuideCountError
from ..models.quality_control_models import SingleInferenceQualityControlResult, MatchSetSingleInferenceQualityControlResult, QualityControlResult
from ..models.mapping_models import GeneralMappingInferenceDict
from typing import DefaultDict, Tuple, Optional, Union
from collections import defaultdict
from ..processing.crispr_editing_processing import get_non_error_dict

# -----------------------------
# QC helpers
# -----------------------------
def iter_all_obs_values(d, contains_sample_barcode: bool):
    """
    Iterate over all observed values safely.
    Works for:
    - Flat dict: {guide: InferenceResult}
    - Flat dict with tuple keys: {(sample, guide): InferenceResult}
    - Nested dict: {sample: {guide: InferenceResult}}
    """
    for val in d.values():
        if isinstance(val, dict):
            # Nested dictionary: sample -> {guide: InferenceResult}
            for inner_val in val.values():
                yield inner_val
        else:
            # Flat dictionary: directly InferenceResult
            yield val

def get_umi_noncollapsed_counts(d, contains_guide_umi: bool, contains_sample_barcode: bool):
    return sum(
        sum(v.observed_value.values()) if contains_guide_umi else v.observed_value
        for v in iter_all_obs_values(d, contains_sample_barcode)
    )


def get_umi_collapsed_counts(d, contains_guide_umi: bool, contains_sample_barcode: bool):
    return sum(
        len(v.observed_value) if contains_guide_umi else v.observed_value
        for v in iter_all_obs_values(d, contains_sample_barcode)
    )


# -----------------------------
# Update num counts in QC results
# -----------------------------

def set_num_non_error_counts(
    qc_result: SingleInferenceQualityControlResult,
    counts_inferred_dict: dict,
    observed_guide_reporter_umi_counts_inferred: Union[GeneralMappingInferenceDict, DefaultDict[str, GeneralMappingInferenceDict]],
    contains_guide_umi: bool,
    contains_sample_barcode: bool
) -> SingleInferenceQualityControlResult:
    """
    Update a QC result object with both non-error counts and total counts.
    """
    if contains_guide_umi:
        # Non-error counts
        qc_result.num_non_error_umi_noncollapsed_counts = get_umi_noncollapsed_counts(
            counts_inferred_dict, contains_guide_umi=True, contains_sample_barcode=contains_sample_barcode
        )
        qc_result.num_non_error_umi_collapsed_counts = get_umi_collapsed_counts(
            counts_inferred_dict, contains_guide_umi=True, contains_sample_barcode=contains_sample_barcode
        )

        # Total counts (all sequences, including errors)
        qc_result.num_total_umi_noncollapsed_counts = get_umi_noncollapsed_counts(
            observed_guide_reporter_umi_counts_inferred, contains_guide_umi=True, contains_sample_barcode=contains_sample_barcode
        )
        qc_result.num_total_umi_collapsed_counts = get_umi_collapsed_counts(
            observed_guide_reporter_umi_counts_inferred, contains_guide_umi=True, contains_sample_barcode=contains_sample_barcode
        )
    else:
        # Non-error counts
        qc_result.num_non_error_counts = get_umi_noncollapsed_counts(
            counts_inferred_dict, contains_guide_umi=False, contains_sample_barcode=contains_sample_barcode
        )
        # Total counts
        qc_result.num_total_counts = get_umi_noncollapsed_counts(
            observed_guide_reporter_umi_counts_inferred, contains_guide_umi=False, contains_sample_barcode=contains_sample_barcode
        )

    return qc_result



def set_num_guide_count_error_types(qc_result, attribute_name: str, observed_guide_reporter_umi_counts_inferred, contains_guide_umi: bool, contains_sample_barcode: bool):
    guide_count_error_type_umi_noncollapsed_count: DefaultDict = defaultdict(int)
    guide_count_error_type_umi_collapsed_count: DefaultDict = defaultdict(int)
    guide_count_error_type_count: DefaultDict = defaultdict(int)

    for obs_value in iter_all_obs_values(observed_guide_reporter_umi_counts_inferred, contains_sample_barcode):
        attr = getattr(obs_value.inferred_value, attribute_name, None)
        error = getattr(attr, "error", None)

        # Only count errors, skip NO_ERROR
        if error is not None and error.guide_count_error_type.name != "NO_ERROR":
            if contains_guide_umi:
                guide_count_error_type_umi_noncollapsed_count[error.guide_count_error_type] += sum(obs_value.observed_value.values())
                guide_count_error_type_umi_collapsed_count[error.guide_count_error_type] += len(obs_value.observed_value)
            else:
                guide_count_error_type_count[error.guide_count_error_type] += obs_value.observed_value

    qc_result.guide_count_error_type_umi_noncollapsed_count = guide_count_error_type_umi_noncollapsed_count
    qc_result.guide_count_error_type_umi_collapsed_count = guide_count_error_type_umi_collapsed_count
    qc_result.guide_count_error_type_count = guide_count_error_type_count
    return qc_result


# -----------------------------
# Main QC
# -----------------------------

def perform_counts_quality_control(
    observed_guide_reporter_umi_counts_inferred:  Union[GeneralMappingInferenceDict, DefaultDict[str, GeneralMappingInferenceDict]],
    contains_guide_umi: bool,
    contains_guide_surrogate: bool,
    contains_guide_barcode: bool,
    contains_sample_barcode: bool,
):
    quality_control_result = QualityControlResult()

    # Helper to process a single attribute
    def set_match_set_single_inference_quality_control_results(attribute_name: str):
        qc = MatchSetSingleInferenceQualityControlResult()
        # PERF/MEM: compute the non-error dict locally — don't retain it as a
        # QC attribute. It's only needed here for counting, and storing it
        # previously duplicated ~78 % of the full result pickle.
        non_error_dict = get_non_error_dict(observed_guide_reporter_umi_counts_inferred, attribute_name)
        qc = set_num_non_error_counts(
            qc,
            non_error_dict,
            observed_guide_reporter_umi_counts_inferred,
            contains_guide_umi,
            contains_sample_barcode
        )
        qc = set_num_guide_count_error_types(
            qc, attribute_name, observed_guide_reporter_umi_counts_inferred, contains_guide_umi, contains_sample_barcode
        )
        return qc

    # Assign QC results
    quality_control_result = QualityControlResult()
    quality_control_result.protospacer_match = set_match_set_single_inference_quality_control_results("protospacer_match")
    if contains_guide_barcode:
        quality_control_result.protospacer_match_barcode_match = set_match_set_single_inference_quality_control_results("protospacer_match_barcode_match")
        if contains_guide_surrogate:
            quality_control_result.protospacer_match_surrogate_match_barcode_match = set_match_set_single_inference_quality_control_results("protospacer_match_surrogate_match_barcode_match")
            quality_control_result.protospacer_mismatch_surrogate_match_barcode_match = set_match_set_single_inference_quality_control_results("protospacer_mismatch_surrogate_match_barcode_match")
    if contains_guide_surrogate:
        quality_control_result.protospacer_match_surrogate_match = set_match_set_single_inference_quality_control_results("protospacer_match_surrogate_match")
        quality_control_result.protospacer_mismatch_surrogate_match = set_match_set_single_inference_quality_control_results("protospacer_mismatch_surrogate_match")
    
    return quality_control_result