from typeguard import typechecked
from typing import Union, List, Mapping, Tuple, Optional, Any, DefaultDict, Dict
from typing import Counter as CounterType
from ..models.mapping_models import GeneralMappingInferenceDict, GeneralMatchCountDict, GeneralMismatchCountDict, WhitelistReporterObservedSequenceMapping

from collections import Counter
from collections import defaultdict
from ..models.mapping_models import MatchSetWhitelistReporterCounterSeriesResults, CompleteInferenceMatchResult, MatchSetSingleInferenceMatchResult, SurrogateProtospacerMismatchSetWhitelistReporterCounterSeriesResults, SurrogateProtospacerMismatchSingleInferenceMatchResult, AllMatchSetWhitelistReporterCounterSeriesResults, InferenceResult
from .crispr_editing_processing import check_match_result_non_error, get_non_error_dict
import pandas as pd

@typechecked
def helper_get_observed_values_given_whitelist_value(whitelist_sequence_list: List[Tuple[str, Optional[str], Optional[str]]], observed_guide_reporter_umi_counts_inferred: GeneralMappingInferenceDict, attribute_name:str, contains_guide_umi: bool, ambiguous_accepted:bool = True) -> WhitelistReporterObservedSequenceMapping:
    if observed_guide_reporter_umi_counts_inferred is None:
        raise ValueError(
            "helper_get_observed_values_given_whitelist_value requires "
            "`observed_guide_reporter_umi_counts_inferred` but the result was built "
            "with `retain_inference_results=False` (the default). Re-run mapping with "
            "`retain_inference_results=True` to enable allele post-processing."
        )
    """
    Provides the set of observed sequences for a given whitelist sequence

    Args:
        whitelist_sequence_list (List[Tuple[str, Optional[str], Optional[str]]]): List of whitelist sequence to retrieve observed sequences
        observed_guide_reporter_umi_counts_inferred (GeneralMappingInferenceDict): Datastructure contain the observed sequence to mapped whitelist sequence
        attribute_name (str): string specifying white mapping match type to consider.
        contains_guide_umi (bool): specifying whether UMI is used
        ambiguous_accepted (bool): specifying whether to consider ambiguous mapping

    Returns:
        WhitelistReporterObservedSequenceMapping: dictionary mapping between the whitelist sequence and the list of observed sequence. Each observed sequence has the corresponding count Union[int, Dict [ str, int ]] (either int if no UMI, or Dict[str, int] for UMI-collapsed and non-collpased count)
    """
    whitelist_sequence_mapping_list: WhitelistReporterObservedSequenceMapping = defaultdict(list)


    # Iterate through the inference results (will attempt to do retrieval for all requested whitelist sequences at once for optimization)
    inferred_value_results: InferenceResult
    observed_sequence: Tuple[str, Optional[str], Optional[str]]
    for observed_sequence, inferred_value_results in get_non_error_dict(observed_guide_reporter_umi_counts_inferred, attribute_name).items():
        #
        #   Get the relevant attributes
        #
        observed_value_counts: Union[int, CounterType[Optional[str]]] = inferred_value_results.observed_value # Read count of the observed sequence
        inferred_value_result: CompleteInferenceMatchResult =  inferred_value_results.inferred_value  # Whitelist inference results
        match_set_single_inference_match_result : Optional[MatchSetSingleInferenceMatchResult] = getattr(inferred_value_result, attribute_name) # Get inference result for specific mapping strategy
        assert match_set_single_inference_match_result is not None, "match_set_single_inference_match_result should not be none since this is from the non error list. Developer error."
        
        matches: pd.DataFrame = match_set_single_inference_match_result.value.matches # Get the list of ambiguous
        
        if not matches.empty:

            # Skip observed sequence if there are multiple matches when not accepting ambiguous mapping
            if (ambiguous_accepted is False) and (matches.shape[0] > 1):
                continue

            # PERF §3.7: itertuples over the underlying numpy values is 5-10x
            # faster than iterrows (no per-row Series construction).
            _all_match_tuples = [tuple(row) for row in matches.itertuples(index=False, name=None)]
            for dict_index in _all_match_tuples:
                if dict_index in whitelist_sequence_list:
                    all_match_sequences = _all_match_tuples
                    total_match_counts = len(all_match_sequences)

                    if contains_guide_umi:
                        assert isinstance(observed_value_counts, Counter), f"For UMI, expecting observed value is a Counter, but type is {type(observed_value_counts)}"

                        # Calculate both UMI noncollapsed and collapsed count
                        observed_sequence_count = {
                            "umi_noncollapsed_count": sum(observed_value_counts.values()),
                            "umi_collapsed_count": len(observed_value_counts.values())
                        }

                        # Add observed sequence and count to mapping
                        whitelist_sequence_mapping_list[dict_index].append( (observed_sequence, observed_sequence_count, total_match_counts, all_match_sequences) )
                    
                    # STANDARD NON-UMI BASED COUNTING
                    else:
                        assert isinstance(observed_value_counts, int), f"For non UMI, expecting observed value is an int, but type is {type(observed_value_counts)}"
                        
                        observed_sequence_count = observed_value_counts
                        
                        # Add observed sequence and count to mapping
                        whitelist_sequence_mapping_list[dict_index].append( (observed_sequence, observed_sequence_count, total_match_counts, all_match_sequences) )

    return whitelist_sequence_mapping_list


#
#   GET THE WHITELIST COUNT PANDAS SERIES
#
# HELPER FUNCTION GETS COUNTS FOR THE THE MATCHES - defined in-function to reduce arguments being passed (NOTE: There is some duplicate code with mismatch counts function - keep in mind if making modifications)
@typechecked
def get_matchset_counterseries(
    observed_guide_reporter_umi_counts_inferred: Union[GeneralMappingInferenceDict, DefaultDict[str, GeneralMappingInferenceDict]],
    whitelist_guide_reporter_df: pd.DataFrame,
    contains_guide_umi: bool,
    contains_sample_barcode: bool,
    attribute_name: str
) -> MatchSetWhitelistReporterCounterSeriesResults:
    
    #
    #   DEFINE THE DEFAULTDICTS FOR COUNTING
    #   Key is the guide, value is the count
    #
    if contains_sample_barcode:
        ambiguous_ignored_umi_noncollapsed_counterdict  = defaultdict(lambda: defaultdict(int))
        ambiguous_ignored_umi_collapsed_counterdict     = defaultdict(lambda: defaultdict(int))
        ambiguous_ignored_counterdict                   = defaultdict(lambda: defaultdict(int))

        ambiguous_accepted_umi_noncollapsed_counterdict = defaultdict(lambda: defaultdict(int))
        ambiguous_accepted_umi_collapsed_counterdict    = defaultdict(lambda: defaultdict(int))
        ambiguous_accepted_counterdict                  = defaultdict(lambda: defaultdict(int))

        ambiguous_spread_umi_noncollapsed_counterdict   = defaultdict(lambda: defaultdict(float))
        ambiguous_spread_umi_collapsed_counterdict      = defaultdict(lambda: defaultdict(float))
        ambiguous_spread_counterdict                    = defaultdict(lambda: defaultdict(float))
    else:
        ambiguous_ignored_umi_noncollapsed_counterdict  = defaultdict(int)
        ambiguous_ignored_umi_collapsed_counterdict     = defaultdict(int)
        ambiguous_ignored_counterdict                   = defaultdict(int)

        ambiguous_accepted_umi_noncollapsed_counterdict = defaultdict(int)
        ambiguous_accepted_umi_collapsed_counterdict    = defaultdict(int)
        ambiguous_accepted_counterdict                  = defaultdict(int)

        ambiguous_spread_umi_noncollapsed_counterdict   = defaultdict(float)
        ambiguous_spread_umi_collapsed_counterdict      = defaultdict(float)
        ambiguous_spread_counterdict                    = defaultdict(float)


    #
    # ITERATE THROUGH THE NON-ERROR INFERRED RESULTS AND FILL THE COUNTS
    #
    def process_sample(observed_dict, cell_barcode=None):
        inferred_value_results: InferenceResult
        for inferred_value_results in get_non_error_dict(observed_dict, attribute_name).values():
            #
            #   Get the relevant attributes
            #
            observed_value_counts: Union[int, CounterType[Optional[str]]] = inferred_value_results.observed_value
            inferred_value_result: CompleteInferenceMatchResult =  inferred_value_results.inferred_value 
            match_set_single_inference_match_result : Optional[MatchSetSingleInferenceMatchResult] = getattr(inferred_value_result, attribute_name)
            assert match_set_single_inference_match_result is not None, "match_set_single_inference_match_result should not be none since this is from the non error list. Developer error."
            
            matches: pd.DataFrame = match_set_single_inference_match_result.value.matches
            if not matches.empty:
                # ITERATE THROUGH MATCHE(S) TO PERFORM COUNTS
                # PERF: itertuples is 5-10x faster than iterrows
                for dict_index in matches.itertuples(index=False, name=None):

                    # Helper for incrementing either flat or nested dicts
                    def add_count(counterdict, value, spread=False):
                        if contains_sample_barcode:
                            if spread:
                                counterdict[cell_barcode][dict_index] += float(value)
                            else:
                                counterdict[cell_barcode][dict_index] += value
                        else:
                            if spread:
                                counterdict[dict_index] += float(value)
                            else:
                                counterdict[dict_index] += value

                    # UMI-BASED COUNTING
                    if contains_guide_umi:
                        assert isinstance(observed_value_counts, Counter), f"For UMI, expecting observed value is a Counter, but type is {type(observed_value_counts)}"
                        total_reads = sum(observed_value_counts.values())
                        collapsed = len(observed_value_counts.values())
                        num_matches = float(matches.shape[0])

                        add_count(ambiguous_accepted_umi_noncollapsed_counterdict, total_reads)
                        add_count(ambiguous_accepted_umi_collapsed_counterdict, collapsed)
                        add_count(ambiguous_spread_umi_noncollapsed_counterdict, total_reads / num_matches, spread=True)
                        add_count(ambiguous_spread_umi_collapsed_counterdict, collapsed / num_matches, spread=True)

                        # If there is no ambiguous matches, then add to ambiguous_ignored counter
                        if matches.shape[0] == 1:
                            add_count(ambiguous_ignored_umi_noncollapsed_counterdict, total_reads)
                            add_count(ambiguous_ignored_umi_collapsed_counterdict, collapsed)

                    # STANDARD NON-UMI BASED COUNTING
                    else:
                        assert isinstance(observed_value_counts, int), f"For non UMI, expecting observed value is an int, but type is {type(observed_value_counts)}"
                        num_matches = float(matches.shape[0])
                        add_count(ambiguous_accepted_counterdict, observed_value_counts)
                        add_count(ambiguous_spread_counterdict, observed_value_counts / num_matches, spread=True)
                        
                        # If there is no ambiguous matches, then add to ambiguous_ignored counter
                        if matches.shape[0] == 1:
                            add_count(ambiguous_ignored_counterdict, observed_value_counts)

    # Process all samples
    if contains_sample_barcode:
        for cell_barcode, observed_dict in observed_guide_reporter_umi_counts_inferred.items():
            process_sample(observed_dict, cell_barcode)
    else:
        process_sample(observed_guide_reporter_umi_counts_inferred)


    # Helper function that converts defaultdict to series
    def create_counterseries(counterdict):
        if contains_sample_barcode:
            # Flatten nested dict
            records = [
                (cell_barcode, *key_tuple, value)
                for cell_barcode, inner in counterdict.items()
                for key_tuple, value in inner.items()
            ]
            columns = ["CellBarcode"] + list(whitelist_guide_reporter_df.columns) + ["value"]
            df = pd.DataFrame.from_records(records, columns=columns)
            return df.set_index(["CellBarcode"] + list(whitelist_guide_reporter_df.columns))["value"]
        else:
            # PERF: build directly from counterdict.items() instead of iterating
            # every whitelist row via .apply(axis=1) (O(N_guides) per series ×
            # 9 series per tier was dominating the run, especially since
            # defaultdict access on a miss inserts a zero entry). Only non-zero
            # keys end up in the Series; mapping_models.__setattr__ already
            # coerces all-zero Series to None, so downstream is unaffected.
            cols = list(whitelist_guide_reporter_df.columns)
            if not counterdict:
                return pd.Series(
                    dtype=float,
                    index=pd.MultiIndex.from_frame(whitelist_guide_reporter_df.iloc[0:0]),
                )
            records = [(*key, value) for key, value in counterdict.items()]
            df = pd.DataFrame.from_records(records, columns=cols + ["value"])
            return df.set_index(cols)["value"]


    #
    #   CONVERT THE COUNT DICTS INTO PANDAS SERIES, since this is a more ideal structure.
    #
    match_set_whitelist_reporter_counter_series_results = MatchSetWhitelistReporterCounterSeriesResults()
    match_set_whitelist_reporter_counter_series_results.ambiguous_ignored_umi_noncollapsed_counterseries = create_counterseries(ambiguous_ignored_umi_noncollapsed_counterdict)
    match_set_whitelist_reporter_counter_series_results.ambiguous_ignored_umi_collapsed_counterseries = create_counterseries(ambiguous_ignored_umi_collapsed_counterdict)
    match_set_whitelist_reporter_counter_series_results.ambiguous_ignored_counterseries = create_counterseries(ambiguous_ignored_counterdict)

    match_set_whitelist_reporter_counter_series_results.ambiguous_accepted_umi_noncollapsed_counterseries = create_counterseries(ambiguous_accepted_umi_noncollapsed_counterdict)
    match_set_whitelist_reporter_counter_series_results.ambiguous_accepted_umi_collapsed_counterseries = create_counterseries(ambiguous_accepted_umi_collapsed_counterdict)
    match_set_whitelist_reporter_counter_series_results.ambiguous_accepted_counterseries = create_counterseries(ambiguous_accepted_counterdict)

    match_set_whitelist_reporter_counter_series_results.ambiguous_spread_umi_noncollapsed_counterseries = create_counterseries(ambiguous_spread_umi_noncollapsed_counterdict)
    match_set_whitelist_reporter_counter_series_results.ambiguous_spread_umi_collapsed_counterseries = create_counterseries(ambiguous_spread_umi_collapsed_counterdict)
    match_set_whitelist_reporter_counter_series_results.ambiguous_spread_counterseries = create_counterseries(ambiguous_spread_counterdict)

    return match_set_whitelist_reporter_counter_series_results



#
# HELPER FUNCTION GETS COUNTS FOR THE THE MISMATCHES - defined in-function to reduce arguments being passed (NOTE: There is some duplicate code with match counts function - keep in mind if making modifications)
#
@typechecked
def get_mismatchset_counterseries(observed_guide_reporter_umi_counts_inferred: Union[GeneralMappingInferenceDict, DefaultDict[str, GeneralMappingInferenceDict]] , whitelist_guide_reporter_df: pd.DataFrame, contains_guide_umi: bool, contains_sample_barcode: bool, attribute_name: str) -> SurrogateProtospacerMismatchSetWhitelistReporterCounterSeriesResults:
    #
    #   DEFINE THE DEFAULTDICTS FOR COUNTING
    #

    if contains_sample_barcode:
        ambiguous_ignored_umi_noncollapsed_match_counterdict  = defaultdict(lambda: defaultdict(int))
        ambiguous_ignored_umi_collapsed_match_counterdict  = defaultdict(lambda: defaultdict(int))
        ambiguous_ignored_match_counterdict  = defaultdict(lambda: defaultdict(int))

        ambiguous_accepted_umi_noncollapsed_match_counterdict = defaultdict(lambda: defaultdict(int))
        ambiguous_accepted_umi_collapsed_match_counterdict  = defaultdict(lambda: defaultdict(int))
        ambiguous_accepted_match_counterdict  = defaultdict(lambda: defaultdict(int))

        ambiguous_spread_umi_noncollapsed_match_counterdict  = defaultdict(lambda: defaultdict(float))
        ambiguous_spread_umi_collapsed_match_counterdict  = defaultdict(lambda: defaultdict(float))
        ambiguous_spread_match_counterdict  = defaultdict(lambda: defaultdict(float))

        # MISMATCH counters (keys are PAIRS of indices, representing the protospacer and surrogate match separately)
        ambiguous_ignored_umi_noncollapsed_mismatch_counterdict  = defaultdict(lambda: defaultdict(int))
        ambiguous_ignored_umi_collapsed_mismatch_counterdict  = defaultdict(lambda: defaultdict(int))
        ambiguous_ignored_mismatch_counterdict = defaultdict(lambda: defaultdict(int))

        ambiguous_accepted_umi_noncollapsed_mismatch_counterdict  =defaultdict(lambda:  defaultdict(int))
        ambiguous_accepted_umi_collapsed_mismatch_counterdict = defaultdict(lambda: defaultdict(int))
        ambiguous_accepted_mismatch_counterdict  = defaultdict(lambda: defaultdict(int))

        ambiguous_spread_umi_noncollapsed_mismatch_counterdict  = defaultdict(lambda: defaultdict(float))
        ambiguous_spread_umi_collapsed_mismatch_counterdict  = defaultdict(lambda: defaultdict(float))
        ambiguous_spread_mismatch_counterdict  = defaultdict(lambda: defaultdict(float))
    else:
        # MATCH counters
        ambiguous_ignored_umi_noncollapsed_match_counterdict  = defaultdict(int)
        ambiguous_ignored_umi_collapsed_match_counterdict  = defaultdict(int)
        ambiguous_ignored_match_counterdict  = defaultdict(int)

        ambiguous_accepted_umi_noncollapsed_match_counterdict  = defaultdict(int)
        ambiguous_accepted_umi_collapsed_match_counterdict  = defaultdict(int)
        ambiguous_accepted_match_counterdict  = defaultdict(int)

        ambiguous_spread_umi_noncollapsed_match_counterdict  = defaultdict(float)
        ambiguous_spread_umi_collapsed_match_counterdict = defaultdict(float)
        ambiguous_spread_match_counterdict  = defaultdict(float)

        # MISMATCH counters (keys are PAIRS of indices, representing the protospacer and surrogate match separately)
        ambiguous_ignored_umi_noncollapsed_mismatch_counterdict  = defaultdict(int)
        ambiguous_ignored_umi_collapsed_mismatch_counterdict  = defaultdict(int)
        ambiguous_ignored_mismatch_counterdict  = defaultdict(int)

        ambiguous_accepted_umi_noncollapsed_mismatch_counterdict   = defaultdict(int)
        ambiguous_accepted_umi_collapsed_mismatch_counterdict = defaultdict(int)
        ambiguous_accepted_mismatch_counterdict  = defaultdict(int)

        ambiguous_spread_umi_noncollapsed_mismatch_counterdict  = defaultdict(float)
        ambiguous_spread_umi_collapsed_mismatch_counterdict = defaultdict(float)
        ambiguous_spread_mismatch_counterdict = defaultdict(float)



    #
    # ITERATE THROUGH THE NON-ERROR INFERRED RESULTS (NOTE: If only one of the protospacer or surrogate matches but not the other, this is treated as an error, and will NOT be counted or considered (even in the single match series). For those counts, should just use the protospacer-only or surrogate-only match results)
    #
    def process_sample(observed_dict, cell_barcode=None):
        inferred_value_results: InferenceResult
        for inferred_value_results in get_non_error_dict(observed_dict, attribute_name).values():
            #
            #   Get the relevant attributes
            #
            observed_value_counts: Union[int, CounterType[Optional[str]]] = inferred_value_results.observed_value
            inferred_value_result: CompleteInferenceMatchResult = inferred_value_results.inferred_value
            surrogate_protospacer_mismatch_single_inference_match_result: Optional[SurrogateProtospacerMismatchSingleInferenceMatchResult] = getattr(inferred_value_result, attribute_name)
            assert surrogate_protospacer_mismatch_single_inference_match_result is not None, "surrogate_protospacer_mismatch_single_inference_match_result should not be none since this is from the non error list. Developer error."
            
            mismatched: bool = surrogate_protospacer_mismatch_single_inference_match_result.value.mismatched
            surrogate_matches: pd.DataFrame = surrogate_protospacer_mismatch_single_inference_match_result.value.surrogate_matches
            protospacer_matches: pd.DataFrame = surrogate_protospacer_mismatch_single_inference_match_result.value.protospacer_matches
            protospacer_surrogate_matches: pd.DataFrame = surrogate_protospacer_mismatch_single_inference_match_result.value.protospacer_surrogate_matches

            # Helper for incrementing nested or flat dicts
            def add_count(counterdict, key, value):
                if contains_sample_barcode:
                    counterdict[cell_barcode][key] += value
                else:
                    counterdict[key] += value

            if mismatched:  # If mismatched, Add count to protospacer/surrogate pair
                assert not surrogate_matches.empty, "Developer error: to be called a mismatch, there must be both separate protospacer and surrogate. No surrogate match (possible no protospacer match)"
                assert not protospacer_matches.empty, "Developer error: to be called a mismatch, there must be both separate protospacer and surrogate. No surrogate match (possible no protospacer match)"
                assert protospacer_surrogate_matches.empty, "Developer error: to be called a mismatch, matches dataframe must be empty."

                # PERF: itertuples is 5-10x faster than iterrows
                for proto_row in protospacer_matches.itertuples(index=False, name=None):
                    for surr_row in surrogate_matches.itertuples(index=False, name=None):
                        dict_index = (proto_row, surr_row)
                        
                        if contains_guide_umi:
                            assert isinstance(observed_value_counts, Counter), f"For UMI, expecting observed value is a Counter, but type is {type(observed_value_counts)}"
                            total_reads = sum(observed_value_counts.values())
                            collapsed = len(observed_value_counts.values())
                            denom = float(protospacer_matches.shape[0] * surrogate_matches.shape[0])

                            add_count(ambiguous_accepted_umi_noncollapsed_mismatch_counterdict, dict_index, total_reads)
                            add_count(ambiguous_accepted_umi_collapsed_mismatch_counterdict, dict_index, collapsed)
                            add_count(ambiguous_spread_umi_noncollapsed_mismatch_counterdict, dict_index, total_reads / denom)
                            add_count(ambiguous_spread_umi_collapsed_mismatch_counterdict, dict_index, collapsed / denom)

                            if (protospacer_matches.shape[0] == 1) and (surrogate_matches.shape[0] == 1):
                                add_count(ambiguous_ignored_umi_noncollapsed_mismatch_counterdict, dict_index, total_reads)
                                add_count(ambiguous_ignored_umi_collapsed_mismatch_counterdict, dict_index, collapsed)
                        else:
                            assert isinstance(observed_value_counts, int), f"For non UMI, expecting observed value is an int, but type is {type(observed_value_counts)}"
                            denom = float(protospacer_matches.shape[0] * surrogate_matches.shape[0])

                            add_count(ambiguous_accepted_mismatch_counterdict, dict_index, observed_value_counts)
                            add_count(ambiguous_spread_mismatch_counterdict, dict_index, observed_value_counts / denom)

                            if (protospacer_matches.shape[0] == 1) and (surrogate_matches.shape[0] == 1):
                                add_count(ambiguous_ignored_mismatch_counterdict, dict_index, observed_value_counts)

            else:  # Matched case
                assert not protospacer_surrogate_matches.empty, f"mismatched==false, but the match dataframe is empty, unexpected paradox. Developer error"
                
                matches = protospacer_surrogate_matches
                # PERF: itertuples is 5-10x faster than iterrows
                for dict_index in matches.itertuples(index=False, name=None):
                    if contains_guide_umi:
                        assert isinstance(observed_value_counts, Counter), f"For UMI, expecting observed value is a Counter, but type is {type(observed_value_counts)}"
                        total_reads = sum(observed_value_counts.values())
                        collapsed = len(observed_value_counts.values())
                        denom = float(matches.shape[0])

                        add_count(ambiguous_accepted_umi_noncollapsed_match_counterdict, dict_index, total_reads)
                        add_count(ambiguous_accepted_umi_collapsed_match_counterdict, dict_index, collapsed)
                        add_count(ambiguous_spread_umi_noncollapsed_match_counterdict, dict_index, total_reads / denom)
                        add_count(ambiguous_spread_umi_collapsed_match_counterdict, dict_index, collapsed / denom)

                        if matches.shape[0] == 1:
                            add_count(ambiguous_ignored_umi_noncollapsed_match_counterdict, dict_index, total_reads)
                            add_count(ambiguous_ignored_umi_collapsed_match_counterdict, dict_index, collapsed)
                    else:
                        assert isinstance(observed_value_counts, int), f"For non UMI, expecting observed value is an int, but type is {type(observed_value_counts)}"
                        denom = float(matches.shape[0])
                        add_count(ambiguous_accepted_match_counterdict, dict_index, observed_value_counts)
                        add_count(ambiguous_spread_match_counterdict, dict_index, observed_value_counts / denom)

                        if matches.shape[0] == 1:
                            add_count(ambiguous_ignored_match_counterdict, dict_index, observed_value_counts)

    if contains_sample_barcode:
        for cell_barcode, observed_dict in observed_guide_reporter_umi_counts_inferred.items():
            process_sample(observed_dict, cell_barcode)
    else:
        process_sample(observed_guide_reporter_umi_counts_inferred)

    #
    # Series creation helpers (now aware of sample barcode)
    #
    def create_match_counterseries(counterdict):
        if contains_sample_barcode:
            records = [
                (cell_barcode, *key_tuple, value)
                for cell_barcode, inner in counterdict.items()
                for key_tuple, value in inner.items()
            ]
            columns = ["CellBarcode"] + list(whitelist_guide_reporter_df.columns) + ["value"]
            df = pd.DataFrame.from_records(records, columns=columns)
            return df.set_index(["CellBarcode"] + list(whitelist_guide_reporter_df.columns))["value"]
        else:
            # PERF: same fix as create_counterseries above (avoid .apply over whitelist).
            cols = list(whitelist_guide_reporter_df.columns)
            if not counterdict:
                return pd.Series(
                    dtype=float,
                    index=pd.MultiIndex.from_frame(whitelist_guide_reporter_df.iloc[0:0]),
                )
            records = [(*key, value) for key, value in counterdict.items()]
            df = pd.DataFrame.from_records(records, columns=cols + ["value"])
            return df.set_index(cols)["value"]

    def create_mismatch_counterseries(counterdict):
        protospacer_match_suffix = "_ProtospacerMatch"
        surrogate_match_suffix = "_SurrogateMatch"
        proto_cols = [c + protospacer_match_suffix for c in whitelist_guide_reporter_df.columns]
        surr_cols  = [c + surrogate_match_suffix  for c in whitelist_guide_reporter_df.columns]

        if contains_sample_barcode:
            records = []
            for cell_barcode, inner in counterdict.items():
                for (protospacer_key, surrogate_key), value in inner.items():
                    records.append((cell_barcode, *protospacer_key, *surrogate_key, value))
            columns = ["CellBarcode"] + proto_cols + surr_cols + ["value"]
            df = pd.DataFrame.from_records(records, columns=columns)
            return df.set_index(["CellBarcode"] + proto_cols + surr_cols)["value"]
        else:
            # PERF: previously this built a cross-product of the whitelist
            # (N_guides^2 rows — 1.4M for a 1186-guide library) and then
            # .apply(axis=1) looked up every (proto, surr) pair in counterdict.
            # Because defaultdict access inserts a zero entry on miss, the dict
            # also grew O(N^2). Build directly from counterdict.items() instead;
            # rows for unseen pairs simply don't appear (they would have been
            # zero anyway, and mapping_models.__setattr__ strips all-zero Series).
            if not counterdict:
                empty_df = pd.DataFrame(columns=proto_cols + surr_cols)
                return pd.Series(
                    dtype=float,
                    index=pd.MultiIndex.from_frame(empty_df),
                )
            records = [(*proto_key, *surr_key, value)
                       for (proto_key, surr_key), value in counterdict.items()]
            df = pd.DataFrame.from_records(records, columns=proto_cols + surr_cols + ["value"])
            return df.set_index(proto_cols + surr_cols)["value"]


    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results = SurrogateProtospacerMismatchSetWhitelistReporterCounterSeriesResults()
    
    # MATCH counters
    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_ignored_umi_noncollapsed_match_counterseries = create_match_counterseries(ambiguous_ignored_umi_noncollapsed_match_counterdict)
    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_ignored_umi_collapsed_match_counterseries = create_match_counterseries(ambiguous_ignored_umi_collapsed_match_counterdict)
    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_ignored_match_counterseries = create_match_counterseries(ambiguous_ignored_match_counterdict)

    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_accepted_umi_noncollapsed_match_counterseries = create_match_counterseries(ambiguous_accepted_umi_noncollapsed_match_counterdict)
    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_accepted_umi_collapsed_match_counterseries = create_match_counterseries(ambiguous_accepted_umi_collapsed_match_counterdict)
    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_accepted_match_counterseries = create_match_counterseries(ambiguous_accepted_match_counterdict)

    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_spread_umi_noncollapsed_match_counterseries = create_match_counterseries(ambiguous_spread_umi_noncollapsed_match_counterdict)
    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_spread_umi_collapsed_match_counterseries = create_match_counterseries(ambiguous_spread_umi_collapsed_match_counterdict)
    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_spread_match_counterseries = create_match_counterseries(ambiguous_spread_match_counterdict)

    # MISMATCH counters
    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_ignored_umi_noncollapsed_mismatch_counterseries = create_mismatch_counterseries(ambiguous_ignored_umi_noncollapsed_mismatch_counterdict)
    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_ignored_umi_collapsed_mismatch_counterseries = create_mismatch_counterseries(ambiguous_ignored_umi_collapsed_mismatch_counterdict)
    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_ignored_mismatch_counterseries = create_mismatch_counterseries(ambiguous_ignored_mismatch_counterdict)

    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_accepted_umi_noncollapsed_mismatch_counterseries = create_mismatch_counterseries(ambiguous_accepted_umi_noncollapsed_mismatch_counterdict)
    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_accepted_umi_collapsed_mismatch_counterseries = create_mismatch_counterseries(ambiguous_accepted_umi_collapsed_mismatch_counterdict)
    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_accepted_mismatch_counterseries = create_mismatch_counterseries(ambiguous_accepted_mismatch_counterdict)

    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_spread_umi_noncollapsed_mismatch_counterseries = create_mismatch_counterseries(ambiguous_spread_umi_noncollapsed_mismatch_counterdict)
    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_spread_umi_collapsed_mismatch_counterseries = create_mismatch_counterseries(ambiguous_spread_umi_collapsed_mismatch_counterdict)
    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_spread_mismatch_counterseries = create_mismatch_counterseries(ambiguous_spread_mismatch_counterdict)

    return surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results


#
# CALLED FUNCTION TO RETRIEVE ALL MATCHSET AND MISMATCHSET COUNTERSERIES
#
@typechecked
def get_counterseries_all_results(observed_guide_reporter_umi_counts_inferred: Union[GeneralMappingInferenceDict, DefaultDict[str, GeneralMappingInferenceDict]], whitelist_guide_reporter_df: pd.DataFrame, contains_guide_barcode: bool, contains_guide_surrogate: bool,  contains_guide_umi: bool, contains_sample_barcode: bool) -> AllMatchSetWhitelistReporterCounterSeriesResults:
    all_match_set_whitelist_reporter_counter_series_results = AllMatchSetWhitelistReporterCounterSeriesResults()
    all_match_set_whitelist_reporter_counter_series_results.protospacer_match = get_matchset_counterseries(observed_guide_reporter_umi_counts_inferred, whitelist_guide_reporter_df, contains_guide_umi, contains_sample_barcode, "protospacer_match")
    if contains_guide_barcode:
        all_match_set_whitelist_reporter_counter_series_results.protospacer_match_barcode_match = get_matchset_counterseries(observed_guide_reporter_umi_counts_inferred, whitelist_guide_reporter_df, contains_guide_umi, contains_sample_barcode, "protospacer_match_barcode_match")
        if contains_guide_surrogate:
            all_match_set_whitelist_reporter_counter_series_results.protospacer_match_surrogate_match_barcode_match = get_matchset_counterseries(observed_guide_reporter_umi_counts_inferred, whitelist_guide_reporter_df, contains_guide_umi, contains_sample_barcode, "protospacer_match_surrogate_match_barcode_match")
            
            all_match_set_whitelist_reporter_counter_series_results.protospacer_mismatch_surrogate_match_barcode_match = get_mismatchset_counterseries(observed_guide_reporter_umi_counts_inferred, whitelist_guide_reporter_df, contains_guide_umi, contains_sample_barcode, "protospacer_mismatch_surrogate_match_barcode_match")
    if contains_guide_surrogate:
        all_match_set_whitelist_reporter_counter_series_results.protospacer_match_surrogate_match = get_matchset_counterseries(observed_guide_reporter_umi_counts_inferred, whitelist_guide_reporter_df, contains_guide_umi, contains_sample_barcode, "protospacer_match_surrogate_match")
        all_match_set_whitelist_reporter_counter_series_results.protospacer_mismatch_surrogate_match = get_mismatchset_counterseries(observed_guide_reporter_umi_counts_inferred, whitelist_guide_reporter_df, contains_guide_umi, contains_sample_barcode, "protospacer_mismatch_surrogate_match")
    
    return all_match_set_whitelist_reporter_counter_series_results
