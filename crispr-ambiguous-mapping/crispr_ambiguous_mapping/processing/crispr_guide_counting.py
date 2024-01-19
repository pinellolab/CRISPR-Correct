from . import guide_raw_fastq_parsing
from . import reporter_tsv_parsing
from . import reporter_umitools_fastq_parsing
from . import sequence_encoding
import pandas as pd
import numpy as np
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from pandarallel import pandarallel
import matplotlib.pyplot as plt
from collections import Counter, defaultdict
from typeguard import typechecked
from datetime import date
import re
from multiprocessing import Pool
from functools import partial
from itertools import repeat
import gzip
import random
from enum import Enum
from typing import Callable
from typing import Union, List, Mapping, Tuple, Optional, Any, DefaultDict
from typing import Counter as CounterType
from concurrent.futures import ProcessPoolExecutor

from . import crispr_sequence_encoding
from . import crispr_guide_inference
from ..models.mapping_models import *
from ..models.editing_models import *
from crispr_editing_processing import check_match_result_non_error, get_non_error_dict


# TODO: There will probably be some type errors with the DefaultDict when testing on non UMI (since it requires CounterType), so make sure to test with different variations of inputs
@typechecked
def get_whitelist_reporter_counts_with_umi(observed_guide_reporter_umi_counts: DefaultDict[Tuple[str,Optional[str],Optional[str]], Union[int, CounterType[Optional[str]]]], whitelist_guide_reporter_df: pd.DataFrame, contains_surrogate:bool = False, contains_barcode:bool = False, contains_umi:bool = False, protospacer_hamming_threshold_strict: Optional[int] = 7, surrogate_hamming_threshold_strict: Optional[int] = 2, barcode_hamming_threshold_strict: Optional[int] = 2, cores: int=1):
    
    # Temporary bug fix. Pad sequences so they are all of same length - encoding numpy matrices requires consistent shape. Still pass the original guide table for selecting the matches.     
    def pad_series(series):
        max_surrogate_len = series.apply(len).max()
        return series.apply(lambda item: item.ljust(max_surrogate_len, 'X'))
    padded_whitelist_guide_reporter_df = whitelist_guide_reporter_df.apply(pad_series, axis=0)
    
    #
    # ENCODE THE WHITELISTED SEQUENCES INTO NUMPY MATRICES - THIS IS REQUIRED FOR HAMMING-BASED MAPPING
    #
    encoded_whitelist_protospacer_sequences_series = crispr_sequence_encoding.encode_guide_series(padded_whitelist_guide_reporter_df["protospacer"])
    if contains_surrogate:
        encoded_whitelist_surrogate_sequences_series = crispr_sequence_encoding.encode_guide_series(padded_whitelist_guide_reporter_df["surrogate"])
    if contains_barcode:
        encoded_whitelist_barcode_sequences_series = crispr_sequence_encoding.encode_guide_series(padded_whitelist_guide_reporter_df["barcode"])

    

    
    #
    #    SET THE PROTOSPACER HAMMING THRESHOLD
    #
    # NOTE: Could there be an even more dynamic threshold, that certain mapped guides can have a higher threshold (i.e. those with many editable bases) and other guides can have low hamming (i.e. those with not many editable bases)
    protospacer_hamming_threshold_dynamic = False
    if protospacer_hamming_threshold_strict is None:
        protospacer_hamming_threshold_dynamic = True
        # TODO: Pass in arguments to set this hamming threshold. 
        protospacer_hamming_threshold: int = crispr_guide_inference.determine_hamming_threshold(whitelist_guide_reporter_df["protospacer"], encoded_whitelist_protospacer_sequences_series, sample_count = 100, quantile = 0.05)
    else:
        protospacer_hamming_threshold: int = protospacer_hamming_threshold_strict
    print("Protospacer hamming threshold is " + str(protospacer_hamming_threshold))

    #
    #   SET THE SURROGATE HAMMING THRESHOLD
    #
    if contains_surrogate:
        surrogate_hamming_threshold_dynamic = False
        if surrogate_hamming_threshold_strict is None:
            surrogate_hamming_threshold_dynamic = True
            # TODO: Pass in arguments to set this hamming threshold. 
            surrogate_hamming_threshold: int = crispr_guide_inference.determine_hamming_threshold(whitelist_guide_reporter_df["surrogate"], encoded_whitelist_surrogate_sequences_series, sample_count = 100, quantile = 0.05)
            
        else:
            surrogate_hamming_threshold: int = surrogate_hamming_threshold_strict
        print("Surrogate hamming threshold is " + str(surrogate_hamming_threshold))

    #
    #   SET THE BARCODE HAMMING THRESHOLD
    #
    if contains_barcode:
        barcode_hamming_threshold_dynamic = False
        if barcode_hamming_threshold_strict is  None:
            barcode_hamming_threshold_dynamic = True
            barcode_hamming_threshold: int = crispr_guide_inference.determine_hamming_threshold(whitelist_guide_reporter_df["barcode"], encoded_whitelist_barcode_sequences_series, sample_count = 100, quantile = 0.05)
        else:
            barcode_hamming_threshold: int = barcode_hamming_threshold_strict
        print("Barcode hamming threshold is " + str(barcode_hamming_threshold))



    #
    #   FROM THE OBSERVED SEQUENCES, INFER THE WHITELIST SEQUENCE
    #
    print("Inferring the true guides from observed guides")

    infer_whitelist_sequence_p = partial(crispr_guide_inference.infer_whitelist_sequence,
            whitelist_guide_reporter_df=whitelist_guide_reporter_df,
            contains_surrogate=contains_surrogate,
            contains_barcode=contains_barcode,
            contains_umi=contains_umi,
            encoded_whitelist_protospacer_sequences_series=encoded_whitelist_protospacer_sequences_series,
            encoded_whitelist_surrogate_sequences_series=encoded_whitelist_surrogate_sequences_series,
            encoded_whitelist_barcode_sequences_series=encoded_whitelist_barcode_sequences_series,
            protospacer_hamming_threshold=protospacer_hamming_threshold, 
            surrogate_hamming_threshold=surrogate_hamming_threshold, 
            barcode_hamming_threshold=barcode_hamming_threshold)

    # Perform inference
    observed_guide_reporter_list = observed_guide_reporter_umi_counts.keys()
    inferred_true_reporter_sequences = None
    if cores > 1:
        print(f"Running inference parallelized with cores {cores}")
        with Pool(cores) as pool:
            inferred_true_reporter_sequences = pool.map(
            infer_whitelist_sequence_p,
            observed_guide_reporter_list
           )
    else:
        print("Running inference non-parallelized")
        inferred_true_reporter_sequences = [infer_whitelist_sequence_p(observed_guide_reporter) for observed_guide_reporter in observed_guide_reporter_list]
    
    # Map inference results to result object
    observed_guide_reporter_umi_counts_inferred: DefaultDict[Tuple[str,Optional[str],Optional[str]], dict] = defaultdict(dict)
    for observed_guide_reporter_key_index, observed_guide_reporter_key in enumerate(observed_guide_reporter_list):
        observed_guide_reporter_umi_counts_inferred[observed_guide_reporter_key] = InferenceResult(
            observed_value=observed_guide_reporter_umi_counts[observed_guide_reporter_key],
            inferred_value=inferred_true_reporter_sequences[observed_guide_reporter_key_index]
        )
    
    print("Completed inference")
    


    #
    #   GET THE WHITELIST COUNT PANDAS SERIES
    #
    # HELPER FUNCTION GETS COUNTS FOR THE THE MATCHES - defined in-function to reduce arguments being passed (NOTE: There is some duplicate code with mismatch counts function - keep in mind if making modifications)
    @typechecked
    def get_matchset_counterseries(attribute_name: str) -> MatchSetWhitelistReporterCounterSeriesResults: 
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
        for inferred_value_results in get_non_error_dict(attribute_name).values():
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
    def get_mismatchset_counterseries(attribute_name: str) -> SurrogateProtospacerMismatchSetWhitelistReporterCounterSeriesResults:
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
        for inferred_value_results in get_non_error_dict(attribute_name).values():
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
    #   Given th above functions, retrieve countseries results for each match type
    #
    all_match_set_whitelist_reporter_counter_series_results = AllMatchSetWhitelistReporterCounterSeriesResults()
    all_match_set_whitelist_reporter_counter_series_results.protospacer_match = get_matchset_counterseries("protospacer_match")
    if contains_barcode:
        all_match_set_whitelist_reporter_counter_series_results.protospacer_match_barcode_match = get_matchset_counterseries("protospacer_match_barcode_match")
        if contains_surrogate:
            all_match_set_whitelist_reporter_counter_series_results.protospacer_match_surrogate_match_barcode_match = get_matchset_counterseries("protospacer_match_surrogate_match_barcode_match")
            
            all_match_set_whitelist_reporter_counter_series_results.protospacer_mismatch_surrogate_match_barcode_match = get_mismatchset_counterseries("protospacer_mismatch_surrogate_match_barcode_match")
    if contains_surrogate:
        all_match_set_whitelist_reporter_counter_series_results.protospacer_match_surrogate_match = get_matchset_counterseries("protospacer_match_surrogate_match")
        all_match_set_whitelist_reporter_counter_series_results.protospacer_mismatch_surrogate_match = get_mismatchset_counterseries("protospacer_mismatch_surrogate_match")









    #
    # PERFORM THE QC
    #
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
        single_inference_quality_control_result.non_error_dict = get_non_error_dict(attribute_name)
        
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
    
    
    count_input = CountInput(whitelist_guide_reporter_df=whitelist_guide_reporter_df,
            contains_surrogate=contains_surrogate,
            contains_barcode=contains_barcode,
            contains_umi=contains_umi,
            protospacer_hamming_threshold_strict=protospacer_hamming_threshold_strict,
            surrogate_hamming_threshold_strict=surrogate_hamming_threshold_strict,
            barcode_hamming_threshold_strict=barcode_hamming_threshold_strict)
    
    return WhitelistReporterCountsResult(all_match_set_whitelist_reporter_counter_series_results=all_match_set_whitelist_reporter_counter_series_results, observed_guide_reporter_umi_counts_inferred=observed_guide_reporter_umi_counts_inferred, quality_control_result=quality_control_result)