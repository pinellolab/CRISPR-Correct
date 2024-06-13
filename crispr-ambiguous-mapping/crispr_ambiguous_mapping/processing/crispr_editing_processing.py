from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Optional, DefaultDict, Union, Tuple, List
import pandas as pd
from typing import Counter as CounterType
from ..models.types import * 
from ..models.mapping_models import CompleteInferenceMatchResult, MatchSetSingleInferenceMatchResult
from ..models.editing_models import (MatchSetWhitelistReporterObservedSequenceCounterSeriesResults, 
                                     MatchSetWhitelistReporterObservedSequenceMutationProfiles, 
                                     LinkedMutationCounters, 
                                     ObservedSequenceMutationProfile)

# Function that has conditional on whether a match result is not an errror
def check_match_result_non_error(match_result):
    return False if match_result is None else match_result.error is None # If match_result is None, treat as error. If match_result is not None, but error is None, then non_error

# Filter dict with observed sequence inference results for only those that do not contain any mapping errors
def get_non_error_dict(observed_guide_reporter_umi_counts_inferred: GeneralMappingInferenceDict, attribute_name: str) -> MatchSetWhitelistReporterObservedSequenceCounterSeriesResults:
    return {observed_guide_reporter_key: observed_guide_reporter_umi_counts_inferred_value for observed_guide_reporter_key, observed_guide_reporter_umi_counts_inferred_value in observed_guide_reporter_umi_counts_inferred.items() if check_match_result_non_error(getattr(observed_guide_reporter_umi_counts_inferred_value.inferred_value, attribute_name))}

#
# Given the datastructure containing the inference results "observed_guide_reporter_umi_counts_inferred", iterate through the entire datastructure to generate
# another datastructure that contains the observed alleles (protospacer/surrogate/barcode) for each whitelist reporter in either a dictionary or dataframe format.
# This is the foundational datastructure used for analyzing the mutations for each whitelist reporter.
#
def get_matchset_alleleseries(observed_guide_reporter_umi_counts_inferred: GeneralMappingInferenceDict, attribute_name: str, contains_surrogate: bool, contains_barcode: bool, contains_umi: bool): 
    #
    #   DEFINE THE DEFAULTDICTS FOR COUNTING
    #
    ambiguous_ignored_umi_noncollapsed_alleledict : GeneralAlleleDict  = defaultdict(lambda: defaultdict(int))
    ambiguous_ignored_umi_collapsed_alleledict : GeneralAlleleDict  = defaultdict(lambda: defaultdict(int))
    ambiguous_ignored_alleledict : GeneralAlleleDict  = defaultdict(lambda: defaultdict(int))

    ambiguous_accepted_umi_noncollapsed_alleledict : GeneralAlleleDict = defaultdict(lambda: defaultdict(int))
    ambiguous_accepted_umi_collapsed_alleledict : GeneralAlleleDict = defaultdict(lambda: defaultdict(int))
    ambiguous_accepted_alleledict : GeneralAlleleDict = defaultdict(lambda: defaultdict(int))

    ambiguous_spread_umi_noncollapsed_alleledict : GeneralAlleleDict = defaultdict(lambda: defaultdict(int))
    ambiguous_spread_umi_collapsed_alleledict : GeneralAlleleDict = defaultdict(lambda: defaultdict(int))
    ambiguous_spread_alleledict : GeneralAlleleDict = defaultdict(lambda: defaultdict(int))

    #
    # ITERATE THROUGH THE NON-ERROR INFERRED RESULTS AND FILL THE COUNTS
    #
    inferred_value_results: InferenceResult
    for observed_sequence, inferred_value_results in get_non_error_dict(observed_guide_reporter_umi_counts_inferred, attribute_name).items():
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
            for whitelist_reporter_series in matches.iterrows(): 
                # UMI-BASED COUNTING
                whitelist_sequence_index = tuple(whitelist_reporter_series[1])
                observed_sequence_index = tuple(observed_sequence)
                if contains_umi:
                    assert isinstance(observed_value_counts, Counter), f"For UMI, expecting observed value is a Counter, but type is {type(observed_value_counts)}"
                    ambiguous_accepted_umi_noncollapsed_alleledict[whitelist_sequence_index][observed_sequence_index] += sum(observed_value_counts.values())
                    ambiguous_accepted_umi_collapsed_alleledict[whitelist_sequence_index][observed_sequence_index] += len(observed_value_counts.values())

                    ambiguous_spread_umi_noncollapsed_alleledict[whitelist_sequence_index][observed_sequence_index] += sum(observed_value_counts.values()) / float(matches.shape[0])
                    ambiguous_spread_umi_collapsed_alleledict[whitelist_sequence_index][observed_sequence_index] += len(observed_value_counts.values()) / float(matches.shape[0])
                    
                    # If there is no ambiguous matches, then add to ambiguous_ignored counter
                    if matches.shape[0] == 1:
                        ambiguous_ignored_umi_noncollapsed_alleledict[whitelist_sequence_index][observed_sequence_index] += sum(observed_value_counts.values())
                        ambiguous_ignored_umi_collapsed_alleledict[whitelist_sequence_index][observed_sequence_index] += len(observed_value_counts.values())

                # STANDARD NON-UMI BASED COUNTING
                else:
                    assert isinstance(observed_value_counts, int), f"For non UMI, expecting observed value is an int, but type is {type(observed_value_counts)}"
                    ambiguous_accepted_alleledict[whitelist_sequence_index][observed_sequence_index] += observed_value_counts
                    ambiguous_spread_alleledict[whitelist_sequence_index][observed_sequence_index] += observed_value_counts / float(matches.shape[0])

                    # If there is no ambiguous matches, then add to ambiguous_ignored counter
                    if matches.shape[0] == 1:
                        ambiguous_ignored_alleledict[whitelist_sequence_index][observed_sequence_index] += observed_value_counts

                        
    # Helper function that converts defaultdict to series
    def create_dict_counterseries(alleledict: GeneralAlleleDict) -> GeneralAlleleCountSeriesDict: 
        return {whitelist_sequence_key: pd.Series(observed_sequence_counterdict) for whitelist_sequence_key, observed_sequence_counterdict in alleledict.items()}
    
    def create_df_from_dict_counterseries(alleledict_counterseries: GeneralAlleleCountSeriesDict, contains_surrogate: bool, contains_barcode: bool):
        whitelist_reporter_columns = ["protospacer"] # For dynamically creating immutable tuple from reporter sequences based on if provided
                
        if contains_surrogate:
            whitelist_reporter_columns.append("surrogate")
        if contains_barcode:
            whitelist_reporter_columns.append("barcode")


        if alleledict_counterseries is not None:
            observed_reporter_df_counts_list = []
            for whitelist_reporter_tuple, observed_reporter_series_counts in alleledict_counterseries.items():
                observed_reporter_df_counts = pd.DataFrame(observed_reporter_series_counts)
                observed_reporter_df_counts.columns = ["counts"]
                observed_reporter_df_counts_index = observed_reporter_df_counts.index.to_frame()
                observed_reporter_df_counts_index.columns = ["observed_" + label for label in whitelist_reporter_columns]
                observed_reporter_df_counts = pd.concat([observed_reporter_df_counts, observed_reporter_df_counts_index], axis=1)
                observed_reporter_df_counts[["true_" + label for label in whitelist_reporter_columns]] = whitelist_reporter_tuple 
                observed_reporter_df_counts_list.append(observed_reporter_df_counts)
            observed_reporter_df_counts_total = pd.concat(observed_reporter_df_counts_list)
            return observed_reporter_df_counts_total
        else:
            return None
    #
    #   CONVERT THE COUNT DICTS INTO PANDAS SERIES, since this is a more ideal structure.
    #
    match_set_whitelist_reporter_observed_sequence_counter_series_results = MatchSetWhitelistReporterObservedSequenceCounterSeriesResults()
    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_ignored_umi_noncollapsed_alleleseries_dict = create_dict_counterseries(ambiguous_ignored_umi_noncollapsed_alleledict)
    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_ignored_umi_collapsed_alleleseries_dict = create_dict_counterseries(ambiguous_ignored_umi_collapsed_alleledict)
    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_ignored_alleleseries_dict = create_dict_counterseries(ambiguous_ignored_alleledict)

    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_accepted_umi_noncollapsed_alleleseries_dict = create_dict_counterseries(ambiguous_accepted_umi_noncollapsed_alleledict)
    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_accepted_umi_collapsed_alleleseries_dict = create_dict_counterseries(ambiguous_accepted_umi_collapsed_alleledict)
    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_accepted_alleleseries_dict = create_dict_counterseries(ambiguous_accepted_alleledict)

    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_spread_umi_noncollapsed_alleleseries_dict = create_dict_counterseries(ambiguous_spread_umi_noncollapsed_alleledict)
    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_spread_umi_collapsed_alleleseries_dict = create_dict_counterseries(ambiguous_spread_umi_collapsed_alleledict)
    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_spread_alleleseries_dict = create_dict_counterseries(ambiguous_spread_alleledict)

    #
    #    Create the dataframes based on the countdicts
    #
    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_ignored_umi_noncollapsed_allele_df = create_df_from_dict_counterseries(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_ignored_umi_noncollapsed_alleleseries_dict, contains_surrogate=contains_surrogate, contains_barcode=contains_barcode)
    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_ignored_umi_collapsed_allele_df = create_df_from_dict_counterseries(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_ignored_umi_collapsed_alleleseries_dict, contains_surrogate=contains_surrogate, contains_barcode=contains_barcode)
    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_ignored_allele_df = create_df_from_dict_counterseries(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_ignored_alleleseries_dict, contains_surrogate=contains_surrogate, contains_barcode=contains_barcode)

    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_accepted_umi_noncollapsed_allele_df = create_df_from_dict_counterseries(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_accepted_umi_noncollapsed_alleleseries_dict, contains_surrogate=contains_surrogate, contains_barcode=contains_barcode)
    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_accepted_umi_collapsed_allele_df = create_df_from_dict_counterseries(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_accepted_umi_collapsed_alleleseries_dict, contains_surrogate=contains_surrogate, contains_barcode=contains_barcode)
    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_accepted_allele_df = create_df_from_dict_counterseries(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_accepted_alleleseries_dict, contains_surrogate=contains_surrogate, contains_barcode=contains_barcode)

    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_spread_umi_noncollapsed_allele_df = create_df_from_dict_counterseries(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_spread_umi_noncollapsed_alleleseries_dict, contains_surrogate=contains_surrogate, contains_barcode=contains_barcode)
    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_spread_umi_collapsed_allele_df = create_df_from_dict_counterseries(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_spread_umi_collapsed_alleleseries_dict, contains_surrogate=contains_surrogate, contains_barcode=contains_barcode)
    match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_spread_allele_df = create_df_from_dict_counterseries(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_spread_alleleseries_dict, contains_surrogate=contains_surrogate, contains_barcode=contains_barcode)
        
    return match_set_whitelist_reporter_observed_sequence_counter_series_results














from copy import deepcopy
import numpy as np
# TODO: Move to package

def pad_series(series, possible_max = 0):
    max_len = max(series.apply(len).max(), possible_max)
    return series.apply(lambda item: item.ljust(max_len, 'X')), max_len

def pad_sequence(sequence, sequence_max=0):
    return sequence.ljust(sequence_max, "X")

# Copying function from crispr-millipede package - don't want import entire package for a single function
def get_substitution_encoding(true_sequence, observed_sequence, skip_index=0):
    assert len(true_sequence) == len(observed_sequence), f"Sequences not same length. True={true_sequence}, Observed={observed_sequence}" # Ensure the aligned sequence (from allele table) is equal size to the reference sequence
    
    nucleotides = ["A","C","T","G","N","-", "X"] # List of possible nucleotides
    encodings_per_position = []
    mismatch_mappings_per_position = []
    for index in range(0, len(true_sequence)): # Iterate through each base and check for substitution
        # TODO Ensure sequences are uppercase
        nucleotides_mm = nucleotides[:]
        nucleotides_mm.remove(true_sequence[index])
        mm_encoding = pd.Series(np.repeat(0, len(nucleotides)-1))
        if observed_sequence[index] == true_sequence[index]: # If the observed sequence is same as reference
            pass
        else:
            mm_index = nucleotides_mm.index(observed_sequence[index])
            mm_encoding[mm_index] = 1
        mismatch_mappings_per_position.append(nucleotides_mm)
        encodings_per_position.append(mm_encoding)

    encodings_per_position_df = pd.DataFrame(encodings_per_position).T
    mismatch_mappings_per_position_df = pd.DataFrame(mismatch_mappings_per_position).T

    encodings_per_position_df.columns = list(true_sequence)
    mismatch_mappings_per_position_df.columns = list(true_sequence)

    mismatch_mappings_per_position_POS_list = np.arange(mismatch_mappings_per_position_df.shape[1]).repeat(mismatch_mappings_per_position_df.shape[0])
    mismatch_mappings_per_position_REF_list = np.asarray(list(true_sequence)).repeat(mismatch_mappings_per_position_df.shape[0]).astype(np.object_)
    mismatch_mappings_per_position_ALT_list = mismatch_mappings_per_position_df.T.values.flatten()
    mismatch_mappings_per_position_full_list = mismatch_mappings_per_position_POS_list.astype(np.str_).astype(object)+mismatch_mappings_per_position_REF_list + np.repeat(">", len(mismatch_mappings_per_position_REF_list)) + mismatch_mappings_per_position_ALT_list
    encodings_per_position_list = encodings_per_position_df.T.values.flatten()

    
    # Encodings per position DF, mismatch mappings per position DF, encodings per position flattened, mismatch mappings per position flattened, mismatch mapping position in flattened list, mismatch mapping ref in flattened list, mismatch mapping alt in flattened list, all substitutions made
    
    index = pd.MultiIndex.from_tuples(zip(mismatch_mappings_per_position_full_list, mismatch_mappings_per_position_POS_list, mismatch_mappings_per_position_REF_list, mismatch_mappings_per_position_ALT_list), names=["FullChange", "Position","Ref", "Alt"])
    
    assert len(encodings_per_position_list) == len(index), f"{len(encodings_per_position_list)} == {len(index)}"
    encodings_per_position_series = pd.Series(encodings_per_position_list, index = index, name="encoding")
    return encodings_per_position_series


def determine_mutations_in_sequence(true_sequence, observed_sequence):
    mutation_array = []
    lowest_len = min(len(true_sequence), len(observed_sequence))
    for i in range(lowest_len):
        if true_sequence[i] != observed_sequence[i]:
            preceding_nt_context = None
            succeeding_nt_context = None
            ref_nt = true_sequence[i]
            alt_nt = observed_sequence[i]
            if i == 0:
                preceding_nt_context = "_"
                succeeding_nt_context = true_sequence[i+1]
            elif i == len(true_sequence)-1:
                preceding_nt_context = true_sequence[i-1]
                succeeding_nt_context = "_"
            else:
                preceding_nt_context = true_sequence[i-1]
                succeeding_nt_context = true_sequence[i+1]
                
            mutation_series = pd.Series((preceding_nt_context, ref_nt, alt_nt,succeeding_nt_context, i), index=["preceding_nt_context", "ref_nt", "alt_nt", "succeeding_nt_context", "position"])
            mutation_array.append(mutation_series)
    observed_sequence_mutation_df = pd.DataFrame(mutation_array)
    return observed_sequence_mutation_df


def get_mutation_profile(match_set_whitelist_reporter_observed_sequence_counter_series_results: MatchSetWhitelistReporterObservedSequenceCounterSeriesResults, whitelist_reporter_df: pd.DataFrame, contains_surrogate: bool, contains_barcode: bool) -> MatchSetWhitelistReporterObservedSequenceMutationProfiles: 
    
    # Function to generate unlinked mutations for particular count type
    def generate_mutations_results(alleleseries: Optional[GeneralAlleleCountSeriesDict], whitelist_reporter_df: pd.DataFrame, contains_surrogate: bool, contains_barcode: bool) -> Optional[ObservedSequenceMutationProfile]:
        if alleleseries is not None:
            
            linked_mutations_whitelist_reporter_dict = {}
            all_observed_protospacer_unlinked_mutations_df_list = []
            if contains_surrogate:
                all_observed_surrogate_unlinked_mutations_df_list = []
            if contains_barcode:
                all_observed_barcode_unlinked_mutations_df_list = []

            for whitelist_reporter_sequence in whitelist_reporter_df.iterrows():
                #
                # GET OBSERVED SEQUENCE FRROM WHITELIST TUPLE KEY (the index of whitelist_reporter_df may have the tuple, though it may not be provided. I also want to ensure that the order of protospacer/surrogate/barcode is maintained if columns are swapped)
                #
                whitelist_sequence_pretuple_list = [] # For dynamically creating immutable tuple from reporter sequences based on if provided
                whitelist_protospacer_sequence = whitelist_reporter_sequence[1]["protospacer"][:]
                whitelist_sequence_pretuple_list.append(whitelist_protospacer_sequence)
                if contains_surrogate:
                    whitelist_surrogate_sequence = whitelist_reporter_sequence[1]["surrogate"][:]
                    whitelist_sequence_pretuple_list.append(whitelist_surrogate_sequence)
                if contains_barcode:
                    whitelist_barcode_sequence = whitelist_reporter_sequence[1]["barcode"][:]
                    whitelist_sequence_pretuple_list.append(whitelist_barcode_sequence)
                whitelist_reporter_tuple = tuple(whitelist_sequence_pretuple_list)
                
                linked_mutations_series_list = []
                try:
                    
                    # CALL THE DICT WITH WHITELIST TUPLE
                    observed_sequences_df = alleleseries[whitelist_reporter_tuple] # Get the observed sequences (DF) for the specific whitelist sequence
                    
                    
                    # PAD THE WHITELIST AND OBSERVED SEQUENCES TO ALLOW ENCODING OF EQUAL LENGTH SEQUENCES
                    whitelist_reporter_sequence_copy = deepcopy(whitelist_reporter_sequence) # Copy the sequences since it will be modified by padding
                    padded_observed_sequences_results: List[List[str], Optional[List[str]], Optional[List[str]]] = []
                    for level in range(observed_sequences_df.index.nlevels): # Iterate over sequence components (protospacer or surrogate or barcode)
                        series_padding_results = pad_series(observed_sequences_df.index.get_level_values(level).to_series(), possible_max=len(whitelist_reporter_sequence_copy[1].iloc[level])) # Pad the observed sequence series (protospacer or surrogate or barcode)
                        padded_observed_sequences_results.append(series_padding_results[0].to_list()) # Append the results to the result list
                        whitelist_reporter_sequence_copy[1].iloc[level] = pad_sequence(whitelist_reporter_sequence_copy[1].iloc[level], sequence_max=series_padding_results[1]) # Pad the result
                    padded_observed_sequence_results_tuples = [tuple(padded_observed_sequence_list) for padded_observed_sequence_list in zip(*padded_observed_sequences_results)]
                    observed_sequences_df.index = pd.MultiIndex.from_tuples(padded_observed_sequence_results_tuples)
                    
                    
                    for observed_sequences, count in observed_sequences_df.items(): # Iterate through each observed sequence to tally mutations
                        # Get protospacer unlinked mutations
                        observed_protospacer_sequence = observed_sequences[0]
                        observed_protospacer_unlinked_mutations_df = determine_mutations_in_sequence(true_sequence=whitelist_reporter_sequence_copy[1]["protospacer"], observed_sequence=observed_protospacer_sequence) # Get DF of mutations in 
                        observed_protospacer_unlinked_mutations_df["count"] = count # Add count
                        observed_protospacer_unlinked_mutations_df.loc[:, whitelist_reporter_sequence_copy[1].index.values] = whitelist_reporter_sequence_copy[1].values # Annotate the DF with the whitelist sequences
                        all_observed_protospacer_unlinked_mutations_df_list.append(observed_protospacer_unlinked_mutations_df) # Add to complete list
                        
                        # Get protospacer linked mutations 
                        observed_linked_mutations_series = get_substitution_encoding(true_sequence=whitelist_reporter_sequence_copy[1]["protospacer"], observed_sequence=observed_protospacer_sequence) # Get DF of mutations in 
                        observed_linked_mutations_series = pd.concat({'protospacer': observed_linked_mutations_series}, names=['SequenceType'])

                        # Get surrogate mutations (if surrogate provided)
                        if contains_surrogate:
                            observed_surrogate_sequence = observed_sequences[1]
                            observed_surrogate_unlinked_mutations_df = determine_mutations_in_sequence(true_sequence=whitelist_reporter_sequence_copy[1]["surrogate"], observed_sequence=observed_surrogate_sequence) # Get DF of mutations in 
                            observed_surrogate_unlinked_mutations_df["count"] = count # Add count
                            observed_surrogate_unlinked_mutations_df.loc[:, whitelist_reporter_sequence_copy[1].index.values] = whitelist_reporter_sequence_copy[1].values# Annotate the DF with the whitelist sequences
                            all_observed_surrogate_unlinked_mutations_df_list.append(observed_surrogate_unlinked_mutations_df)# Add to complete list
                            
                            observed_surrogate_linked_mutations_series = get_substitution_encoding(true_sequence=whitelist_reporter_sequence_copy[1]["surrogate"], observed_sequence=observed_surrogate_sequence) # Get DF of mutations in 
                            observed_surrogate_linked_mutations_series = pd.concat({'surrogate': observed_surrogate_linked_mutations_series}, names=['SequenceType'])
                            observed_linked_mutations_series = pd.concat([observed_linked_mutations_series,observed_surrogate_linked_mutations_series])
                            
                        
                        # Get barcode mutations (if barcode provided)
                        if contains_barcode:
                            observed_barcode_sequence = observed_sequences[2]
                            observed_barcode_unlinked_mutations_df = determine_mutations_in_sequence(true_sequence=whitelist_reporter_sequence_copy[1]["barcode"], observed_sequence=observed_barcode_sequence) # Get DF of mutations in 
                            observed_barcode_unlinked_mutations_df["count"] = count # Add count
                            observed_barcode_unlinked_mutations_df.loc[:, whitelist_reporter_sequence_copy[1].index.values] = whitelist_reporter_sequence_copy[1].values # Annotate the DF with the whitelist sequences
                            all_observed_barcode_unlinked_mutations_df_list.append(observed_barcode_unlinked_mutations_df) # Add to complete list
                            
                            observed_barcode_linked_mutations_series = get_substitution_encoding(true_sequence=whitelist_reporter_sequence_copy[1]["barcode"], observed_sequence=observed_barcode_sequence) # Get DF of mutations in 
                            observed_barcode_linked_mutations_series = pd.concat({'barcode': observed_barcode_linked_mutations_series}, names=['SequenceType'])
                            observed_linked_mutations_series = pd.concat([observed_linked_mutations_series,observed_barcode_linked_mutations_series])
                        
                        observed_linked_mutations_series["count"] = count # Add count
                        linked_mutations_series_list.append(observed_linked_mutations_series)
                    
                    # Add linked mutations DF to dict
                    linked_mutations_df = pd.concat(linked_mutations_series_list, axis=1).transpose()
                    linked_mutations_df.index = observed_sequences_df.index # Set the index of the encodings to be the observed sequences
                    linked_mutations_whitelist_reporter_dict.update({whitelist_reporter_tuple:linked_mutations_df})
                except KeyError as e:
                    print(f"No observed sequences found for: {whitelist_reporter_tuple}")
                
                
                
            all_observed_protospacer_unlinked_mutations_df = pd.concat(all_observed_protospacer_unlinked_mutations_df_list) # Concat list of all whitelist mutations to sequence dataframe
            observed_sequence_mutations = ObservedSequenceMutationProfile(all_observed_protospacer_unlinked_mutations_df=all_observed_protospacer_unlinked_mutations_df, linked_mutations_whitelist_reporter_dict=linked_mutations_whitelist_reporter_dict) # Initialize results object with groupby results
            
            # Add surrogate result to result object (if provided)
            if contains_surrogate:
                all_observed_surrogate_unlinked_mutations_df = pd.concat(all_observed_surrogate_unlinked_mutations_df_list)
                observed_sequence_mutations.all_observed_surrogate_unlinked_mutations_df=all_observed_surrogate_unlinked_mutations_df
                
            # Add barcode result to result object (if provided)
            if contains_barcode:
                all_observed_barcode_unlinked_mutations_df = pd.concat(all_observed_barcode_unlinked_mutations_df_list)
                observed_sequence_mutations.all_observed_barcode_unlinked_mutations_df=all_observed_barcode_unlinked_mutations_df
            
            return observed_sequence_mutations
        else:
            return None
    
    # Generate mutation counts for each count type
    mutations_results = MatchSetWhitelistReporterObservedSequenceMutationProfiles()
    print("Generating ambiguous_ignored_umi_noncollapsed_mutations")
    mutations_results.ambiguous_ignored_umi_noncollapsed_mutations = generate_mutations_results(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_ignored_umi_noncollapsed_alleleseries_dict, whitelist_reporter_df, contains_surrogate, contains_barcode)
    print("Generating ambiguous_ignored_umi_collapsed_mutations")
    mutations_results.ambiguous_ignored_umi_collapsed_mutations = generate_mutations_results(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_ignored_umi_collapsed_alleleseries_dict, whitelist_reporter_df, contains_surrogate, contains_barcode)
    print("Generating ambiguous_ignored_unlinked_mutations")
    mutations_results.ambiguous_ignored_mutations = generate_mutations_results(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_ignored_alleleseries_dict, whitelist_reporter_df, contains_surrogate, contains_barcode)

    print("Generating ambiguous_accepted_umi_noncollapsed_mutations")
    mutations_results.ambiguous_accepted_umi_noncollapsed_mutations = generate_mutations_results(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_accepted_umi_noncollapsed_alleleseries_dict, whitelist_reporter_df, contains_surrogate, contains_barcode)
    print("Generating ambiguous_accepted_umi_collapsed_mutations")
    mutations_results.ambiguous_accepted_umi_collapsed_mutations = generate_mutations_results(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_accepted_umi_collapsed_alleleseries_dict, whitelist_reporter_df, contains_surrogate, contains_barcode)
    print("Generating ambiguous_accepted_mutations")
    mutations_results.ambiguous_accepted_mutations = generate_mutations_results(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_accepted_alleleseries_dict, whitelist_reporter_df, contains_surrogate, contains_barcode)

    print("Generating ambiguous_spread_umi_noncollapsed_mutations")
    mutations_results.ambiguous_spread_umi_noncollapsed_mutations = generate_mutations_results(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_spread_umi_noncollapsed_alleleseries_dict, whitelist_reporter_df, contains_surrogate, contains_barcode)
    print("Generating ambiguous_spread_umi_collapsed_mutations")
    mutations_results.ambiguous_spread_umi_collapsed_mutations = generate_mutations_results(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_spread_umi_collapsed_alleleseries_dict, whitelist_reporter_df, contains_surrogate, contains_barcode)
    print("Generating ambiguous_spread_mutations")
    mutations_results.ambiguous_spread_mutations = generate_mutations_results(match_set_whitelist_reporter_observed_sequence_counter_series_results.ambiguous_spread_alleleseries_dict, whitelist_reporter_df, contains_surrogate, contains_barcode)
    
    return mutations_results

 
def tally_linked_mutation_count_per_sequence(mutations_results: MatchSetWhitelistReporterObservedSequenceMutationProfiles, contains_surrogate: bool, contains_barcode: bool, count_attribute_name: str="ambiguous_accepted_umi_noncollapsed_mutations") -> LinkedMutationCounters:
    protospacer_total_mutation_counter: CounterType = Counter()

    surrogate_total_mutation_counter: Optional[CounterType] = None
    barcode_total_mutation_counter: Optional[CounterType] = None
    if contains_surrogate:
        surrogate_total_mutation_counter = Counter()
    if contains_barcode:
        barcode_total_mutation_counter = Counter()

    # TODO: THIS IS HARDCODED TO ambiguous_accepted_umi_noncollapsed
        
    mutation_results_typed = getattr(mutations_results, count_attribute_name)
    if mutation_results_typed is not None:
        for linked_mutations_whitelist_reporter_df in mutation_results_typed.linked_mutations_whitelist_reporter_dict.values():
            accetable_SNP_columns = (linked_mutations_whitelist_reporter_df.columns.get_level_values("SequenceType") != "count") & (linked_mutations_whitelist_reporter_df.columns.get_level_values("Ref") != "X") & (linked_mutations_whitelist_reporter_df.columns.get_level_values("Alt") != "X")


            protospacer_columns = (linked_mutations_whitelist_reporter_df.columns.get_level_values("SequenceType") == "protospacer")
            protospacer_mutation_count_per_allele_series = linked_mutations_whitelist_reporter_df.loc[:, protospacer_columns & accetable_SNP_columns].sum(axis=1)
            for allele_index, mutation_count in enumerate(protospacer_mutation_count_per_allele_series):
                protospacer_total_mutation_counter[mutation_count] += linked_mutations_whitelist_reporter_df["count"][allele_index]

            if contains_surrogate:
                surrogate_columns = (linked_mutations_whitelist_reporter_df.columns.get_level_values("SequenceType") == "surrogate")
                surrogate_mutation_count_per_allele_series = linked_mutations_whitelist_reporter_df.loc[:, surrogate_columns & accetable_SNP_columns].sum(axis=1)
                for allele_index, mutation_count in enumerate(surrogate_mutation_count_per_allele_series):
                    surrogate_total_mutation_counter[mutation_count] += linked_mutations_whitelist_reporter_df["count"][allele_index]

            if contains_barcode:
                barcode_columns = (linked_mutations_whitelist_reporter_df.columns.get_level_values("SequenceType") == "barcode")
                barcode_mutation_count_per_allele_series = linked_mutations_whitelist_reporter_df.loc[:, barcode_columns & accetable_SNP_columns].sum(axis=1)
                for allele_index, mutation_count in enumerate(barcode_mutation_count_per_allele_series):
                    barcode_total_mutation_counter[mutation_count] += linked_mutations_whitelist_reporter_df["count"][allele_index]
    
    return LinkedMutationCounters(protospacer_total_mutation_counter=protospacer_total_mutation_counter,
                                surrogate_total_mutation_counter=surrogate_total_mutation_counter,
                                barcode_total_mutation_counter=barcode_total_mutation_counter)
