from typeguard import typechecked
from typing import Union, List, Mapping, Tuple, Optional, Any, DefaultDict
from typing import Counter as CounterType
from collections import Counter
from collections import defaultdict
from ..models.mapping_models import MatchSetWhitelistReporterCounterSeriesResults, CompleteInferenceMatchResult, MatchSetSingleInferenceMatchResult, SurrogateProtospacerMismatchSetWhitelistReporterCounterSeriesResults, SurrogateProtospacerMismatchSingleInferenceMatchResult, AllMatchSetWhitelistReporterCounterSeriesResults
from .crispr_editing_processing import check_match_result_non_error, get_non_error_dict
import pandas as pd



#
#   GET THE WHITELIST COUNT PANDAS SERIES
#
# HELPER FUNCTION GETS COUNTS FOR THE THE MATCHES - defined in-function to reduce arguments being passed (NOTE: There is some duplicate code with mismatch counts function - keep in mind if making modifications)
@typechecked
def get_matchset_counterseries(observed_guide_reporter_umi_counts_inferred: DefaultDict[Tuple[str,Optional[str],Optional[str]], dict], whitelist_guide_reporter_df: pd.DataFrame, contains_umi: bool, attribute_name: str) -> MatchSetWhitelistReporterCounterSeriesResults: 
    #
    #   DEFINE THE DEFAULTDICTS FOR COUNTING
    #
    ambiguous_ignored_umi_noncollapsed_counterdict : DefaultDict[Tuple[str, Optional[str], Optional[str]], int]  = defaultdict(int)
    ambiguous_ignored_umi_collapsed_counterdict : DefaultDict[Tuple[str, Optional[str], Optional[str]], int]  = defaultdict(int)
    ambiguous_ignored_counterdict : DefaultDict[Tuple[str, Optional[str], Optional[str]], int]  = defaultdict(int)

    ambiguous_accepted_umi_noncollapsed_counterdict : DefaultDict[Tuple[str, Optional[str], Optional[str]], int]  = defaultdict(int)
    ambiguous_accepted_umi_collapsed_counterdict : DefaultDict[Tuple[str, Optional[str], Optional[str]], int]  = defaultdict(int)
    ambiguous_accepted_counterdict : DefaultDict[Tuple[str, Optional[str], Optional[str]], int]  = defaultdict(int)

    ambiguous_spread_umi_noncollapsed_counterdict : DefaultDict[Tuple[str, Optional[str], Optional[str]], float]  = defaultdict(float)
    ambiguous_spread_umi_collapsed_counterdict : DefaultDict[Tuple[str, Optional[str], Optional[str]], float]  = defaultdict(float)
    ambiguous_spread_counterdict : DefaultDict[Tuple[str, Optional[str], Optional[str]], float]  = defaultdict(float)

    #
    # ITERATE THROUGH THE NON-ERROR INFERRED RESULTS AND FILL THE COUNTS
    #
    for inferred_value_results in get_non_error_dict(observed_guide_reporter_umi_counts_inferred, attribute_name).values():
        #
        #   Get the relevant attributes
        #
        observed_value_counts: Union[CounterType[Optional[str]], int] = inferred_value_results.observed_value
        inferred_value_result: CompleteInferenceMatchResult =  inferred_value_results.inferred_value 
        match_set_single_inference_match_result : Optional[MatchSetSingleInferenceMatchResult] = getattr(inferred_value_result, attribute_name)
        assert match_set_single_inference_match_result is not None, "match_set_single_inference_match_result should not be none since this is from the non error list. Developer error."
        
        matches: pd.DataFrame = match_set_single_inference_match_result.value.matches
        if not matches.empty:
            # ITERATE THROUGH MATCHE(S) TO PERFORM COUNTS
            for whitelist_reporter_series in matches.iterrows(): 
                # UMI-BASED COUNTING
                dict_index = tuple(whitelist_reporter_series[1])
                if contains_umi:
                    assert isinstance(observed_value_counts, Counter), f"For UMI, expecting observed value is a Counter, but type is {type(observed_value_counts)}"
                    ambiguous_accepted_umi_noncollapsed_counterdict[dict_index] += sum(observed_value_counts.values())
                    ambiguous_accepted_umi_collapsed_counterdict[dict_index] += len(observed_value_counts.values())

                    ambiguous_spread_umi_noncollapsed_counterdict[dict_index] += sum(observed_value_counts.values()) / float(matches.shape[0])
                    ambiguous_spread_umi_collapsed_counterdict[dict_index] += len(observed_value_counts.values()) / float(matches.shape[0])
                    
                    # If there is no ambiguous matches, then add to ambiguous_ignored counter
                    if matches.shape[0] == 1:
                        ambiguous_ignored_umi_noncollapsed_counterdict[dict_index] += sum(observed_value_counts.values())
                        ambiguous_ignored_umi_collapsed_counterdict[dict_index] += len(observed_value_counts.values())
                
                # STANDARD NON-UMI BASED COUNTING
                else:
                    assert isinstance(observed_value_counts, int), f"For non UMI, expecting observed value is an int, but type is {type(observed_value_counts)}"
                    ambiguous_accepted_counterdict[dict_index] += observed_value_counts
                    ambiguous_spread_counterdict[dict_index] += observed_value_counts / float(matches.shape[0])
                    
                    # If there is no ambiguous matches, then add to ambiguous_ignored counter
                    if matches.shape[0] == 1:
                        ambiguous_ignored_counterdict[dict_index] += observed_value_counts
    
    # Helper function that converts defaultdict to series
    def create_counterseries(counterdict: DefaultDict[Tuple[str, Optional[str], Optional[str]], Union[int, float]]) -> pd.Series:
        counterseries: pd.Series = whitelist_guide_reporter_df.apply(lambda reporter: counterdict[tuple(reporter)], axis=1)
        counterseries.index = pd.MultiIndex.from_frame(whitelist_guide_reporter_df)
        return counterseries
    

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
def get_mismatchset_counterseries(observed_guide_reporter_umi_counts_inferred: DefaultDict[Tuple[str,Optional[str],Optional[str]], dict], attribute_name: str, whitelist_guide_reporter_df: pd.DataFrame, contains_umi: bool) -> SurrogateProtospacerMismatchSetWhitelistReporterCounterSeriesResults:
    #
    #   DEFINE THE DEFAULTDICTS FOR COUNTING
    #
    # MATCH counters
    ambiguous_ignored_umi_noncollapsed_match_counterdict : DefaultDict[Tuple[str, Optional[str], Optional[str]], int]  = defaultdict(int)
    ambiguous_ignored_umi_collapsed_match_counterdict : DefaultDict[Tuple[str, Optional[str], Optional[str]], int]  = defaultdict(int)
    ambiguous_ignored_match_counterdict : DefaultDict[Tuple[str, Optional[str], Optional[str]], int]  = defaultdict(int)

    ambiguous_accepted_umi_noncollapsed_match_counterdict : DefaultDict[Tuple[str, Optional[str], Optional[str]], int]  = defaultdict(int)
    ambiguous_accepted_umi_collapsed_match_counterdict : DefaultDict[Tuple[str, Optional[str], Optional[str]], int]  = defaultdict(int)
    ambiguous_accepted_match_counterdict : DefaultDict[Tuple[str, Optional[str], Optional[str]], int]  = defaultdict(int)

    ambiguous_spread_umi_noncollapsed_match_counterdict : DefaultDict[Tuple[str, Optional[str], Optional[str]], float]  = defaultdict(float)
    ambiguous_spread_umi_collapsed_match_counterdict : DefaultDict[Tuple[str, Optional[str], Optional[str]], float]  = defaultdict(float)
    ambiguous_spread_match_counterdict : DefaultDict[Tuple[str, Optional[str], Optional[str]], float]  = defaultdict(float)

    # MISMATCH counters (keys are PAIRS of indices, representing the protospacer and surrogate match separately)
    ambiguous_ignored_umi_noncollapsed_mismatch_counterdict : DefaultDict[Tuple[Tuple[str, Optional[str], Optional[str]], Tuple[str, Optional[str], Optional[str]]], int]  = defaultdict(int)
    ambiguous_ignored_umi_collapsed_mismatch_counterdict : DefaultDict[Tuple[Tuple[str, Optional[str], Optional[str]], Tuple[str, Optional[str], Optional[str]]], int]  = defaultdict(int)
    ambiguous_ignored_mismatch_counterdict : DefaultDict[Tuple[Tuple[str, Optional[str], Optional[str]], Tuple[str, Optional[str], Optional[str]]], int]  = defaultdict(int)

    ambiguous_accepted_umi_noncollapsed_mismatch_counterdict : DefaultDict[Tuple[Tuple[str, Optional[str], Optional[str]],Tuple[str, Optional[str], Optional[str]]], int]  = defaultdict(int)
    ambiguous_accepted_umi_collapsed_mismatch_counterdict : DefaultDict[Tuple[Tuple[str, Optional[str], Optional[str]],Tuple[str, Optional[str], Optional[str]]], int]  = defaultdict(int)
    ambiguous_accepted_mismatch_counterdict : DefaultDict[Tuple[Tuple[str, Optional[str], Optional[str]],Tuple[str, Optional[str], Optional[str]]], int]  = defaultdict(int)

    ambiguous_spread_umi_noncollapsed_mismatch_counterdict : DefaultDict[Tuple[Tuple[str, Optional[str], Optional[str]],Tuple[str, Optional[str], Optional[str]]], float]  = defaultdict(float)
    ambiguous_spread_umi_collapsed_mismatch_counterdict : DefaultDict[Tuple[Tuple[str, Optional[str], Optional[str]],Tuple[str, Optional[str], Optional[str]]], float]  = defaultdict(float)
    ambiguous_spread_mismatch_counterdict : DefaultDict[Tuple[Tuple[str, Optional[str], Optional[str]],Tuple[str, Optional[str], Optional[str]]], float]  = defaultdict(float)

    #
    # ITERATE THROUGH THE NON-ERROR INFERRED RESULTS (NOTE: If only one of the protospacer or surrogate matches but not the other, this is treated as an error, and will NOT be counted or considered (even in the single match series). For those counts, should just use the protospacer-only or surrogate-only match results)
    #
    for inferred_value_results in get_non_error_dict(observed_guide_reporter_umi_counts_inferred, attribute_name).values():
        #
        #   Get the relevant attributes
        #
        observed_value_counts: Union[CounterType[Optional[str]], int] = inferred_value_results.observed_value
        inferred_value_result: CompleteInferenceMatchResult =  inferred_value_results.inferred_value
        surrogate_protospacer_mismatch_single_inference_match_result : Optional[SurrogateProtospacerMismatchSingleInferenceMatchResult] = getattr(inferred_value_result, attribute_name)
        assert surrogate_protospacer_mismatch_single_inference_match_result is not None, "surrogate_protospacer_mismatch_single_inference_match_result should not be none since this is from the non error list. Developer error."
        
        mismatched: bool = surrogate_protospacer_mismatch_single_inference_match_result.value.mismatched
        surrogate_matches: pd.DataFrame = surrogate_protospacer_mismatch_single_inference_match_result.value.surrogate_matches
        protospacer_matches: pd.DataFrame = surrogate_protospacer_mismatch_single_inference_match_result.value.protospacer_matches
        protospacer_surrogate_matches: pd.DataFrame = surrogate_protospacer_mismatch_single_inference_match_result.value.protospacer_surrogate_matches



        if mismatched: # If mismatched, Add count to protospacer/surrogate pair. NOTE: Both the protospacer/surrogate matches should be gauranteed
            assert not surrogate_matches.empty, "Developer error: to be called a mismatch, there must be both separate protospacer and surrogate. No surrogate match (possible no protospacer match)"
            assert not protospacer_matches.empty, "Developer error: to be called a mismatch, there must be both separate protospacer and surrogate. No surrogate match (possible no protospacer match)"
            assert protospacer_surrogate_matches.empty, "Developer error: to be called a mismatch, matches dataframe must be empty."

            # NOTE: For the mismatches, we want to tally for the specific protospacer_match/surrogate_match pair, so we iterate through both match dataframes to tally
            for protospacer_matched_whitelist_reporter_series in protospacer_matches.iterrows():
                for surrogate_matched_whitelist_reporter_series in surrogate_matches.iterrows():
                    
                    # NOTE: For the counter_dict key, the order of the pairwise tuple is PROTOSPACER_MATCH, SURROGATE_MATCH.
                    dict_index = (tuple(protospacer_matched_whitelist_reporter_series[1]), tuple(surrogate_matched_whitelist_reporter_series[1]))
                    
                    if contains_umi:
                        assert isinstance(observed_value_counts, Counter), f"For UMI, expecting observed value is a Counter, but type is {type(observed_value_counts)}"
                        
                        ambiguous_accepted_umi_noncollapsed_mismatch_counterdict[dict_index] += sum(observed_value_counts.values())
                        ambiguous_accepted_umi_collapsed_mismatch_counterdict[dict_index] += len(observed_value_counts.values())

                        ambiguous_spread_umi_noncollapsed_mismatch_counterdict[dict_index] += sum(observed_value_counts.values()) / float(protospacer_matches.shape[0]*surrogate_matches.shape[0]) # NOTE: The denominator for the "spreaded count" is intuitively the number of pairs from the protospacer/surrogate-matches, so just the product
                        ambiguous_spread_umi_collapsed_mismatch_counterdict[dict_index] += len(observed_value_counts.values()) / float(protospacer_matches.shape[0]*surrogate_matches.shape[0])
                        
                        # If there is no ambiguous matches, then add to ambiguous_ignored counter
                        if (protospacer_matches.shape[0] == 1) and (surrogate_matches.shape[0] == 1):
                            ambiguous_ignored_umi_noncollapsed_mismatch_counterdict[dict_index] += sum(observed_value_counts.values())
                            ambiguous_ignored_umi_collapsed_mismatch_counterdict[dict_index] += len(observed_value_counts.values())
                    else:
                        assert isinstance(observed_value_counts, int), f"For non UMI, expecting observed value is an int, but type is {type(observed_value_counts)}"

                        ambiguous_accepted_mismatch_counterdict[dict_index] += observed_value_counts
                        ambiguous_spread_mismatch_counterdict[dict_index] += observed_value_counts / float(protospacer_matches.shape[0]*surrogate_matches.shape[0])
                        
                        # If there is no ambiguous matches, then add to ambiguous_ignored counter
                        if (protospacer_matches.shape[0] == 1) and (surrogate_matches.shape[0] == 1):
                            ambiguous_ignored_mismatch_counterdict[dict_index] += observed_value_counts


        else: # If matched, add count normally.(NOTE: As mentioned above, "If only one of the protospacer or surrogate matches but not the other, this is treated as an error, and will NOT be counted or considered (even in the single match series). For those counts, should just use the protospacer-only or surrogate-only match results")
            assert not protospacer_surrogate_matches.empty, f"mismatched==false, but the match dataframe is empty, unexpected paradox. Developer error"
            
            matches = protospacer_surrogate_matches
            # NOTE/TODO: FROM HERE IS THE SAME COUNTING LOGIC AS THE MATCH COUNTING FUNCTION - can modularize.
            for whitelist_reporter_series in matches.iterrows():
                dict_index = tuple(whitelist_reporter_series[1])
                if contains_umi:
                    assert isinstance(observed_value_counts, Counter), f"For UMI, expecting observed value is a Counter, but type is {type(observed_value_counts)}"
                    
                    ambiguous_accepted_umi_noncollapsed_match_counterdict[dict_index] += sum(observed_value_counts.values())
                    ambiguous_accepted_umi_collapsed_match_counterdict[dict_index] += len(observed_value_counts.values())

                    ambiguous_spread_umi_noncollapsed_match_counterdict[dict_index] += sum(observed_value_counts.values()) / float(matches.shape[0])
                    ambiguous_spread_umi_collapsed_match_counterdict[dict_index] += len(observed_value_counts.values()) / float(matches.shape[0])
                    
                    # If there is no ambiguous matches, then add to ambiguous_ignored counter
                    if matches.shape[0] == 1:
                        ambiguous_ignored_umi_noncollapsed_match_counterdict[dict_index] += sum(observed_value_counts.values())
                        ambiguous_ignored_umi_collapsed_match_counterdict[dict_index] += len(observed_value_counts.values())
                else:
                    assert isinstance(observed_value_counts, int), f"For non UMI, expecting observed value is an int, but type is {type(observed_value_counts)}"
                    ambiguous_accepted_match_counterdict[dict_index] += observed_value_counts
                    ambiguous_spread_match_counterdict[dict_index] += observed_value_counts / float(matches.shape[0])
                    
                    # If there is no ambiguous matches, then add to ambiguous_ignored counter
                    if matches.shape[0] == 1:
                        ambiguous_ignored_match_counterdict[dict_index] += observed_value_counts
    
    def create_match_counterseries(counterdict: DefaultDict[Tuple[str, Optional[str], Optional[str]], Union[int, float]]) -> pd.Series:
        counterseries: pd.Series = whitelist_guide_reporter_df.apply(lambda reporter: counterdict[tuple(reporter)], axis=1)
        counterseries.index = pd.MultiIndex.from_frame(whitelist_guide_reporter_df)
        return counterseries
    
    def create_mismatch_counterseries(counterdict: DefaultDict[Tuple[Tuple[str, Optional[str], Optional[str]], Tuple[str, Optional[str], Optional[str]]], Union[int, float]]) -> pd.Series:
        protospacer_match_suffix = "_ProtospacerMatch"
        surrogate_match_suffix = "_SurrogateMatch"
        whitelist_guide_reporter_df_product = whitelist_guide_reporter_df.merge(whitelist_guide_reporter_df, how='cross', suffixes=(protospacer_match_suffix, surrogate_match_suffix))
        counterseries: pd.Series = whitelist_guide_reporter_df_product.apply(lambda pairwise_reporter: counterdict[tuple(pairwise_reporter[[index for index in pairwise_reporter.index if index.endswith(protospacer_match_suffix)]]), tuple(pairwise_reporter[[index for index in pairwise_reporter.index if index.endswith(surrogate_match_suffix)]])], axis=1)
        counterseries.index = pd.MultiIndex.from_frame(whitelist_guide_reporter_df_product)
        return counterseries
    

    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results = SurrogateProtospacerMismatchSetWhitelistReporterCounterSeriesResults()
    
    # MATCH counters
    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_ignored_umi_noncollapsed_match_counterseries = ambiguous_ignored_umi_noncollapsed_match_counterdict
    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_ignored_umi_collapsed_match_counterseries = ambiguous_ignored_umi_collapsed_match_counterdict
    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_ignored_match_counterseries = ambiguous_ignored_match_counterdict

    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_accepted_umi_noncollapsed_match_counterseries = ambiguous_accepted_umi_noncollapsed_match_counterdict
    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_accepted_umi_collapsed_match_counterseries  = ambiguous_accepted_umi_collapsed_match_counterdict
    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_accepted_match_counterseries = ambiguous_accepted_match_counterdict

    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_spread_umi_noncollapsed_match_counterseries = ambiguous_spread_umi_noncollapsed_match_counterdict
    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_spread_umi_collapsed_match_counterseries = ambiguous_spread_umi_collapsed_match_counterdict
    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_spread_match_counterseries = ambiguous_spread_match_counterdict

    # MISMATCH counters (keys are PAIRS of indices, representing the protospacer and surrogate match separately)
    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_ignored_umi_noncollapsed_mismatch_counterseries = ambiguous_ignored_umi_noncollapsed_mismatch_counterdict
    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_ignored_umi_collapsed_mismatch_counterseries = ambiguous_ignored_umi_collapsed_mismatch_counterdict
    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_ignored_mismatch_counterseries = ambiguous_ignored_mismatch_counterdict

    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_accepted_umi_noncollapsed_mismatch_counterseries = ambiguous_accepted_umi_noncollapsed_mismatch_counterdict
    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_accepted_umi_collapsed_mismatch_counterseries = ambiguous_accepted_umi_collapsed_mismatch_counterdict
    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_accepted_mismatch_counterseries = ambiguous_accepted_mismatch_counterdict

    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_spread_umi_noncollapsed_mismatch_counterseries = ambiguous_spread_umi_noncollapsed_mismatch_counterdict
    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_spread_umi_collapsed_mismatch_counterseries = ambiguous_spread_umi_collapsed_mismatch_counterdict
    surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results.ambiguous_spread_mismatch_counterseries = ambiguous_spread_mismatch_counterdict

    return surrogate_protospacer_mismatch_set_whitelist_reporter_counter_series_results


#
# CALLED FUNCTION TO RETRIEVE ALL MATCHSET AND MISMATCHSET COUNTERSERIES
#
@typechecked
def get_counterseries_all_results(observed_guide_reporter_umi_counts_inferred: DefaultDict[Tuple[str,Optional[str],Optional[str]], dict], whitelist_guide_reporter_df: pd.DataFrame, contains_barcode: bool, contains_surrogate: bool,  contains_umi: bool) -> AllMatchSetWhitelistReporterCounterSeriesResults:
    all_match_set_whitelist_reporter_counter_series_results = AllMatchSetWhitelistReporterCounterSeriesResults()
    all_match_set_whitelist_reporter_counter_series_results.protospacer_match = get_matchset_counterseries(observed_guide_reporter_umi_counts_inferred, whitelist_guide_reporter_df, contains_umi, "protospacer_match")
    if contains_barcode:
        all_match_set_whitelist_reporter_counter_series_results.protospacer_match_barcode_match = get_matchset_counterseries(observed_guide_reporter_umi_counts_inferred, whitelist_guide_reporter_df, contains_umi, "protospacer_match_barcode_match")
        if contains_surrogate:
            all_match_set_whitelist_reporter_counter_series_results.protospacer_match_surrogate_match_barcode_match = get_matchset_counterseries(observed_guide_reporter_umi_counts_inferred, whitelist_guide_reporter_df, contains_umi, "protospacer_match_surrogate_match_barcode_match")
            
            all_match_set_whitelist_reporter_counter_series_results.protospacer_mismatch_surrogate_match_barcode_match = get_mismatchset_counterseries(observed_guide_reporter_umi_counts_inferred, whitelist_guide_reporter_df, contains_umi, "protospacer_mismatch_surrogate_match_barcode_match")
    if contains_surrogate:
        all_match_set_whitelist_reporter_counter_series_results.protospacer_match_surrogate_match = get_matchset_counterseries(observed_guide_reporter_umi_counts_inferred, whitelist_guide_reporter_df, contains_umi, "protospacer_match_surrogate_match")
        all_match_set_whitelist_reporter_counter_series_results.protospacer_mismatch_surrogate_match = get_mismatchset_counterseries(observed_guide_reporter_umi_counts_inferred, whitelist_guide_reporter_df, contains_umi, "protospacer_mismatch_surrogate_match")
    
    return all_match_set_whitelist_reporter_counter_series_results
