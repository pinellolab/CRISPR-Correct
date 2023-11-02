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
from typing import Union, List, Mapping, Tuple, Optional, Any
from concurrent.futures import ProcessPoolExecutor

from . import sequence_encoding

###
### GUIDE MAPPING MAIN FUNCTIONS
###
class GuideCountError(Enum):
    NO_MATCH = "No guide found within hamming distance"
    MULTIPLE_MATCH = "Multiple exact matches found for guide (likely a truncated guide read assuming guide series is unique)"
    NO_PROTOSPACER_WITH_SAME_LENGTH = "No whitelisted guides with same length as observed guide - maybe try enabling truncating whitelisted guides"
    NO_BARCODE_SAME_LENGTH = "No whitelisted barcode with the same length as the observed"
    NO_SURROGATE_SAME_LENGTH = "No whitelisted surrogate with the same length as the observed"
    NO_MATCH_PROTOSPACER_HAMMING_THRESHOLD = "No whitelisted guide that is below the provided guide hamming distance threshold"
    NO_MATCH_BARCODE_HAMMING_THRESHOLD = "No whitelisted barcode that is below the provided barcode hamming distance threshold"
    NO_MATCH_SURROGATE_HAMMING_THRESHOLD = "The inferred whitelisted surrogate does not match with the observed surrogate below the given hamming distance threshold"
    MULTIPLE_MATCH_EXACT = "Multiple exact guide matches, double check that there are no duplicates in your guide library (especially if truncation is enabled)"
    NO_PROTOSPACER_MATCH_MISSING_INFO = "The protospacer was not provided"
    NO_SURROGATE_MATCH_MISSING_INFO = "The surrogate was not provided"
    NO_BARCODE_MATCH_MISSING_INFO = "The barcode was not provided"
    NO_MATCH_OBSERVED_SURROGATE_SHORTER = "The observed surrogate is shorter than the inferred surrogate"
    
   
'''
    Given an observed, potentially self-edited guide (in 'row["observed_sequence"]'), try and find the true guide sequence from the whitelist ("guide_sequences_series") based on hamming distance
    # TODO (3/1/2023): This is the main bottleneck in terms of computational efficiency, see if it can be sped up
'''

@typechecked
def infer_true_guides(observed_guide_sequence: str, whitelist_guide_sequences_series: Union[List[str], pd.Series],
encoded_whitelist_guide_sequences_series, consider_truncated_sequences: bool = False, hamming_threshold: int = 3, verbose_result = False):
    #observed_guide_sequence = str(observed_guide_row["observed_sequence"])
    
    # If considering truncated guides, truncate all whitelist guides to the same size of the observed guide, else only consider whitelisted guides of the same length (likely all guides provided are same length of 20nt)
    if consider_truncated_sequences == True:
        whitelist_guide_sequences_series = whitelist_guide_sequences_series.apply(lambda guide: guide[0:len(observed_guide_sequence)])
    else:
        whitelist_guide_sequences_series = whitelist_guide_sequences_series[whitelist_guide_sequences_series.str.len() == len(observed_guide_sequence)]
        if len(whitelist_guide_sequences_series) == 0:
            return GuideCountError.NO_PROTOSPACER_WITH_SAME_LENGTH 
        

    # Determine if there are exact matches, hopefully just a single match
    whitelist_guide_sequences_series_match = whitelist_guide_sequences_series[whitelist_guide_sequences_series == observed_guide_sequence]
    
    # If there is a single exact match, great, no need for fancy mat
    if len(whitelist_guide_sequences_series_match) == 1: # Exact match found, return
        return tuple([str(whitelist_guide_sequences_series_match.index[0])])
    
    # If no matches, possible a self-edit or a sequencing error, try and find guide with lowest hamming distance
    elif len(whitelist_guide_sequences_series_match) == 0: # No match found, search based on hamming distance
        
        # Encode the whitelisted guides 
        # NOTE 20221202: Potentially improve efficiency by passing in the encoded guide series (assuming no truncation) so that this does not have to be ran on every guide
        #guide_sequences_series_encoded = encode_guide_series(whitelist_guide_sequences_series_match)
        
        # Encode the observed guide
        try:
            observed_guide_sequence_encoded = sequence_encoding.encode_DNA_base_vectorized(sequence_encoding.numpify_string_vectorized(observed_guide_sequence)) 
        except Exception as e:
            print(f"Error on observed guide sequence: {observed_guide_sequence}")
            raise e
        # Calculate the hamming distance of the guide with all whitelisted guides - vectorized operation
        observed_guide_sequence_dists = sequence_encoding.retrieve_hamming_distance_whitelist(observed_guide_sequence_encoded, encoded_whitelist_guide_sequences_series)
        
        # Get the minimum hamming distance calculated
        hamming_min = observed_guide_sequence_dists.min()
        
        # Get all whitelisted guides with the minimum hamming distance (could be multiple)
        guides_with_hamming_min = whitelist_guide_sequences_series[np.where(observed_guide_sequence_dists == hamming_min)[0]]
        
        # If the minimum hamming distance is greater than the specified threshold, then the guide is too ambigious to assign, so no match.
        if hamming_min >= hamming_threshold:
            if verbose_result:
                return {"Error": GuideCountError.NO_MATCH_PROTOSPACER_HAMMING_THRESHOLD, "hamming_min": hamming_min}
            else:
                return GuideCountError.NO_MATCH_PROTOSPACER_HAMMING_THRESHOLD
        
        # If there are multiple guides with the minimum hamming distance, then the guide is ambigious, so no mapping (due to multiple match)
        # NOTE (3/30/2023): This could potentially be an issue with igRNAs maybe?
        elif len(guides_with_hamming_min) > 1:
            if verbose_result:
                return {"Error": GuideCountError.MULTIPLE_MATCH,
                "exact_match": False, 
                "num_matches": len(guides_with_hamming_min),
                "matches": guides_with_hamming_min,
                "hamming_min": hamming_min}
            else:
                return GuideCountError.MULTIPLE_MATCH
        
        # Else if there is 1 guide with the match, then return the match
        else:
            return tuple([str(guides_with_hamming_min.index[0])])
    
    # Else if there are multiple exact match, which should never occur unless the whitelisted guide list is not unique, then return multiple match.
    else:
        
        if verbose_result:
            return {"Error": GuideCountError.MULTIPLE_MATCH_EXACT,
            "exact_match": True, 
            "num_matches": len(whitelist_guide_sequences_series_match),
            "matches": whitelist_guide_sequences_series_match,
            "hamming_min": 0}
        else:
            return GuideCountError.MULTIPLE_MATCH_EXACT
        #raise Exception("Multiple exact matches of the provided whitelisted guides - there are likely duplicates in the provided whitelist, please remove. Observed guide={}, guide matches={}".format(observed_guide_sequence, guide_sequences_series_match)) # NOTE 12/6/22: REMOVED THIS EXCEPTION - another reason is from truncated guides having multiple matches. In production code, just make sure to ensure that the whitelist is the set.

@typechecked
def infer_true_reporters(observed_reporter_sequences, whitelist_guide_reporter_df: pd.DataFrame,
encoded_whitelist_guide_sequences_series, encoded_whitelist_barcodes_series, surrogate_hamming_threshold: int = 10, barcode_hamming_threshold: int = 2, hamming_threshold: int = 7, verbose_result = False):
    for element in observed_reporter_sequences:
        if (element is None) or (element == "None") or (element == ""):
            if verbose_result:
                return {"Error": GuideCountError.NO_MATCH_MISSING_INFO, "message": "Observed protospacer/surrogate/barcode is None"}
            else:
                return GuideCountError.NO_MATCH_MISSING_INFO
    
    observed_reporter_sequences = pd.Series(observed_reporter_sequences, index=["protospacer", "surrogate", "barcode"])
    
    # Check if the observed protospacer length is the same length of the whitelisted guides
    if len(observed_reporter_sequences["protospacer"]) != encoded_whitelist_guide_sequences_series.shape[1]: 
        if verbose_result:
            return {"Error": GuideCountError.NO_PROTOSPACER_WITH_SAME_LENGTH, "message": f"Observed protospacer {observed_reporter_sequences['protospacer']} not of correct length: {encoded_whitelist_guide_sequences_series.shape[1]}"}
        else:
            return GuideCountError.NO_PROTOSPACER_WITH_SAME_LENGTH
    if len(observed_reporter_sequences["barcode"]) != encoded_whitelist_barcodes_series.shape[1]: 
        if verbose_result:
            return {"Error": GuideCountError.NO_BARCODE_SAME_LENGTH, "message": f"Observed barcode {observed_reporter_sequences['barcode']} not of correct length: {encoded_whitelist_barcodes_series.shape[1]}"}
        else:
            return GuideCountError.NO_BARCODE_SAME_LENGTH   
            
    # Determine if there are exact matches, hopefully just a single match
    whitelist_guide_reporter_df_match = whitelist_guide_reporter_df[(whitelist_guide_reporter_df["protospacer"]==observed_reporter_sequences["protospacer"]) & (whitelist_guide_reporter_df["surrogate"]==observed_reporter_sequences["surrogate"]) & (whitelist_guide_reporter_df["barcode"]==observed_reporter_sequences["barcode"])] #whitelist_guide_sequences_series == observed_guide_sequence
    
    # If there is a single exact match, great, no need for fancy math
    if whitelist_guide_reporter_df_match.shape[0] == 1: # Exact match found, return
        return tuple(whitelist_guide_reporter_df_match.iloc[0])
    
    # If no matches, possible a self-edit or a sequencing error, try and find guide with lowest hamming distance
    elif whitelist_guide_reporter_df_match.shape[0] == 0: # No match found, search based on hamming distance
        
        ###
        ### FILTER BY BARCODE (NOTE: 4/2/23: This is the main difference with the traditional filtering by protospacer)
        ### TODO 4/2/23: Assumes perfect barcode match, but in the future I can also select for barcodes of 1 hamming - 2 hamming if there is no matches with 0 hamming
        ###
        observed_barcode_encoded = sequence_encoding.encode_DNA_base_vectorized(sequence_encoding.numpify_string_vectorized(observed_reporter_sequences["barcode"]))  # Encode the observed barcode
        observed_barcode_dists = sequence_encoding.retrieve_hamming_distance_whitelist(observed_barcode_encoded, encoded_whitelist_barcodes_series) # Retrieve hamming distance with whitelist barcode
        barcode_hamming_min = observed_barcode_dists.min() # Get the barcode with the minimum  hamming distance
        if barcode_hamming_min >= barcode_hamming_threshold: # If the barcode surpasses the threshold, fail the read due to no barcode. 
            if verbose_result:
                return {"Error": GuideCountError.NO_MATCH_BARCODE_HAMMING_THRESHOLD, "barcode_hamming_min": barcode_hamming_min, "message": f"No barcode below threshold {barcode_hamming_threshold}"}
            else:
                return GuideCountError.NO_MATCH_BARCODE_HAMMING_THRESHOLD
            
        barcode_matches_indices = np.where(observed_barcode_dists == barcode_hamming_min)[0] # Get the indices of ALL barcode matches
        
        whitelist_guide_reporter_df_barcode = whitelist_guide_reporter_df.iloc[barcode_matches_indices] # Get the whitelist reporters with the matched barcode(s)
        encoded_whitelist_guide_sequences_series_barcode = encoded_whitelist_guide_sequences_series[barcode_matches_indices]
        
        ###
        ### Map the protospacer among those with the correct barcode
        ###
        observed_guide_sequence_encoded = sequence_encoding.encode_DNA_base_vectorized(sequence_encoding.numpify_string_vectorized(observed_reporter_sequences["protospacer"]))  # Encode the observed protospacer
        
        
        observed_guide_sequence_dists = sequence_encoding.retrieve_hamming_distance_whitelist(observed_guide_sequence_encoded, encoded_whitelist_guide_sequences_series_barcode) # Calculate the hamming distance of the guide with all whitelisted guides - vectorized operation
        
        
        hamming_min = observed_guide_sequence_dists.min() # Get the minimum hamming distance calculated
        
        reporters_with_hamming_min_df = whitelist_guide_reporter_df_barcode.iloc[np.where(observed_guide_sequence_dists == hamming_min)[0]] # Get all whitelisted guides with the minimum hamming distance (could be multiple)
        
        if hamming_min >= hamming_threshold: # If the minimum hamming distance is greater than the specified threshold, then the guide is too ambigious to assign, so no match.
            if verbose_result:
                return {"Error": GuideCountError.NO_MATCH_PROTOSPACER_HAMMING_THRESHOLD, "hamming_min": hamming_min}
            else:
                return GuideCountError.NO_MATCH_PROTOSPACER_HAMMING_THRESHOLD
        
        # If there are multiple guides with the minimum hamming distance, then the guide is ambigious, so no mapping (due to multiple match)
        # NOTE (3/30/2023): This could potentially be an issue with igRNAs maybe? (Response 6/30/23) No, the barcode should be able to distinguish between igRNAs
        elif reporters_with_hamming_min_df.shape[0] > 1:
            if verbose_result:
                return {"Error": GuideCountError.MULTIPLE_MATCH,
                "exact_match": False, 
                "num_matches": reporters_with_hamming_min_df.shape[0],
                "matches": reporters_with_hamming_min_df,
                "hamming_min": hamming_min}
            else:
                return GuideCountError.MULTIPLE_MATCH
        
        # Else if there is 1 guide with the match, then double check that the observed surrogate matches the mapped surrogate (or if it is due to recombination)
        else:
            inferred_reporter_sequences = tuple(reporters_with_hamming_min_df.iloc[0])
            inferred_surrogate_sequence = inferred_reporter_sequences[1]
            observed_surrogate_sequence = observed_reporter_sequences["surrogate"]
            
            if len(observed_surrogate_sequence) >= len(inferred_surrogate_sequence):
                observed_surrogate_sequence = observed_surrogate_sequence[-len(inferred_surrogate_sequence):] # Because their may be slippage of the polyT upstream of surrogate, slice relative to the downstream end.
                surrogate_hamming_distance = sequence_encoding.determine_hamming_distance_classic(inferred_surrogate_sequence, observed_surrogate_sequence)
                if surrogate_hamming_distance >= surrogate_hamming_threshold:
                    if verbose_result:
                        return {"Error": GuideCountError.NO_MATCH_SURROGATE_HAMMING_THRESHOLD, "surrogate_hamming_distance": surrogate_hamming_distance, "inferred_reporter_sequences": inferred_reporter_sequences}
                    else:
                        return GuideCountError.NO_MATCH_SURROGATE_HAMMING_THRESHOLD
                else:
                    return inferred_reporter_sequences
            else:
                if verbose_result:
                        return {"Error": GuideCountError.NO_MATCH_OBSERVED_SURROGATE_SHORTER, "inferred_reporter_sequences": inferred_reporter_sequences}
                else:
                    return GuideCountError.NO_MATCH_OBSERVED_SURROGATE_SHORTER
    
    # Else if there are multiple exact match, which should never occur unless the whitelisted guide list is not unique, then return multiple match.
    else:
        if verbose_result:
            return {"Error": GuideCountError.MULTIPLE_MATCH_EXACT,
            "exact_match": True, 
            "num_matches": whitelist_guide_reporter_df_match.shape[0],
            "matches": whitelist_guide_reporter_df_match,
            "hamming_min": 0}
        else:
            return GuideCountError.MULTIPLE_MATCH_EXACT
        #raise Exception("Multiple exact matches of the provided whitelisted guides - there are likely duplicates in the provided whitelist, please remove. Observed guide={}, guide matches={}".format(observed_guide_sequence, guide_sequences_series_match)) # NOTE 12/6/22: REMOVED THIS EXCEPTION - another reason is from truncated guides having multiple matches. In production code, just make sure to ensure that the whitelist is the set.

     
'''
    This performs simulation of determining how many mutations it takes for a guide to be ambigiously mapped based on hamming distance.

    This is performed on *sample_count* number of randomly selected guides, and the distribution of mutation count is determined. The 5% quantile is used as the threshold.
    
    This is useful in determing the ideal hamming distance threshold specific to a guide library
'''
@typechecked
def determine_hamming_threshold(whitelist_sequences_series: Union[List[str],pd.Series], encoded_whitelist_sequences_series, sample_count: int = 100, quantile: float = 0.05) -> int:
    #encoded_whitelist_guide_sequences_series = encode_guide_series(guide_sequences_series)
    
    mutation_count_until_nonunique = []
    for i in range(sample_count):
        # Sample a guide from whitelist
        sampled_guide = whitelist_sequences_series.sample()[0]
        
        # Generate order of base positions to mutate sequentially.
        guide_position_list = list(range(len(sampled_guide)))
        random.shuffle(guide_position_list)
        
        # Create temporary variable to represent mutated guide
        current_guide_sequence = sampled_guide
        
        # Iterate through positions to mutate
        for iteration, position in enumerate(guide_position_list):
            
            # Mutate the guide
            current_guide_sequence_separated = list(current_guide_sequence)
            guide_position_nt = current_guide_sequence_separated[position]
            nt_list = ["A", "C", "T", "G"]
            nt_list.remove(guide_position_nt.upper())
            new_nt = random.sample(nt_list, 1)[0]
            current_guide_sequence_separated[position] = new_nt
            current_guide_sequence = "".join(current_guide_sequence_separated)
            current_guide_sequence_encoded = sequence_encoding.encode_DNA_base_vectorized(sequence_encoding.numpify_string_vectorized(current_guide_sequence)) 
            
            # Calculate the hamming distance of the mutated guide to all other guides
            hamming_distances =  sequence_encoding.retrieve_hamming_distance_whitelist(current_guide_sequence_encoded, encoded_whitelist_sequences_series)
            if len(np.where(hamming_distances == hamming_distances.min())[0]) > 1: # TODO (3/30/23): Can also add to the conditional whether the minimum hamming distance guide is still the original guide - but it is probably rare cases where this would not be the case, so not too important to implement
                mutation_count_until_nonunique.append(iteration+1)
                break   
    mutation_count_until_nonunique = pd.Series(mutation_count_until_nonunique)
    
    # From the mutation count distribution, calculate the threshold based on the provided quantile.
    hamming_distance_threshold = mutation_count_until_nonunique.quantile(quantile)
    return int(hamming_distance_threshold) + 1 # Take ceil of result (to be conservative)


'''
    Take in input FASTQ filename, and a set of whitelisted guide sequences
'''
@typechecked
def get_whitelist_guide_counts(observed_guide_sequences_counts: Counter, whitelist_guide_sequences_series: pd.Series, hamming_threshold_strict: int = 3, hamming_threshold_dynamic: bool = False, cores: int=1):

    whitelist_guide_sequences_series.index = whitelist_guide_sequences_series.values
    # Create observed guide DF to contain all information
    # The "observed sequence" column represents all *unique* observed sequences that may have self-edits/errors. The "observed counts" column represents the count of each observed count.
    observed_guides_df = pd.DataFrame({"observed_sequence":[str(sequence) for sequence in observed_guide_sequences_counts.keys()], "observed_counts":observed_guide_sequences_counts.values()})
 
    # Get the hamming distance threshold. THe hamming distance must be below this threshold to assign an observed guide to a whitelist guide.
    encoded_whitelist_guide_sequences_series = sequence_encoding.encode_guide_series(whitelist_guide_sequences_series)
    if hamming_threshold_dynamic:
        hamming_threshold = int(determine_hamming_threshold(whitelist_guide_sequences_series, encoded_whitelist_guide_sequences_series, sample_count = 100, quantile = 0.05))
        print("Hamming threshold is " + str(hamming_threshold))
    else:
        hamming_threshold = hamming_threshold_strict
        
    # Infer whitelist guides from observed guides
    print("Inferring the true guides from observed guides")

    infer_true_guides_p = partial(infer_true_guides,
            whitelist_guide_sequences_series=whitelist_guide_sequences_series,
            encoded_whitelist_guide_sequences_series=encoded_whitelist_guide_sequences_series,
            hamming_threshold=hamming_threshold, verbose_result=False)

    inferred_true_guide_sequences = None
    
    if cores > 1:
        with Pool(cores) as pool:
            print("Inferencing parallel")
            inferred_true_guide_sequences = pool.map(
            infer_true_guides_p,
            observed_guides_df["observed_sequence"]
            )
    else:
        print("Inferencing sequentially")
        inferred_true_guide_sequences = [infer_true_guides_p(observed_sequence) for observed_sequence in observed_guides_df["observed_sequence"]]    
    
    print("Completed inference")

    observed_guides_df["inferred_guides"] = inferred_true_guide_sequences
    
    '''
        QC
    '''
    # QC: Determine the number of guides passed by determining that the result is not an error
    observed_guides_df_passed_inference = observed_guides_df[observed_guides_df["inferred_guides"].apply(lambda guide : type(guide) != GuideCountError)]
    # QC: Calculate number of guides that were unassigned
    observed_guides_df_no_match = observed_guides_df[observed_guides_df["inferred_guides"].apply(lambda guide : guide == GuideCountError.NO_MATCH)]
    # QC: Calculate number of guides with multiple inferred guides
    observed_guides_df_multiple_match = observed_guides_df[observed_guides_df["inferred_guides"].apply(lambda guide : guide == GuideCountError.MULTIPLE_MATCH)]
    # QC: Calculate percent mapped
    percent_mapped = observed_guides_df_passed_inference["observed_counts"].sum()/observed_guides_df["observed_counts"].sum()
    
    # Retrieve the observed sequences that were mapped and set the inferred guides
    inferred_guide_sequence_counter = Counter()
    for _, row in observed_guides_df_passed_inference.iterrows():
        inferred_guide_sequence_counter[row["inferred_guides"][0]] += row["observed_counts"]
    
    whitelist_guide_sequences_series_counts = whitelist_guide_sequences_series.apply(lambda guide: inferred_guide_sequence_counter[guide])
    whitelist_guide_sequences_series_counts.index = whitelist_guide_sequences_series
    
    multi_index = pd.MultiIndex.from_arrays([whitelist_guide_sequences_series], names=['protospacer'])
    whitelist_guide_sequences_series_counts.index = multi_index

    qc_dict = {"guide_sequences_unassigned_counts":observed_guides_df_no_match["observed_counts"].sum(), "guide_sequences_multiple_counts": observed_guides_df_multiple_match["observed_counts"].sum(), "total_guide_counts": observed_guides_df["observed_counts"].sum(), "percent_mapped": percent_mapped}
    
    return whitelist_guide_sequences_series_counts, observed_guides_df, qc_dict

@typechecked
def get_whitelist_reporter_counts(observed_guide_reporters_counts: Counter, whitelist_guide_reporter_df: pd.DataFrame, surrogate_hamming_threshold_strict: int = 2, barcode_hamming_threshold_strict: int = 2, hamming_threshold_strict: int = 7, hamming_threshold_dynamic: bool = False, cores: int=1):

    whitelist_guide_reporter_df.index = whitelist_guide_reporter_df["protospacer"]
    
    # Create observed guide DF to contain all information
    # The "observed sequence" column represents all *unique* observed sequences that may have self-edits/errors. The "observed counts" column represents the count of each observed count.
    observed_guides_df = pd.DataFrame({"observed_sequence":observed_guide_reporters_counts.keys(), "observed_counts":observed_guide_reporters_counts.values()})
 
    # Get the hamming distance threshold. THe hamming distance must be below this threshold to assign an observed guide to a whitelist guide.
    encoded_whitelist_guide_sequences_series = sequence_encoding.encode_guide_series(whitelist_guide_reporter_df["protospacer"])
    encoded_whitelist_barcodes_series = sequence_encoding.encode_guide_series(whitelist_guide_reporter_df["barcode"])
    
    if hamming_threshold_dynamic:
        hamming_threshold = determine_hamming_threshold(whitelist_guide_reporter_df["protospacer"], encoded_whitelist_guide_sequences_series, sample_count = 100, quantile = 0.05)
        print("Hamming threshold is " + str(hamming_threshold))
    else:
        hamming_threshold = hamming_threshold_strict
        
    # Infer whitelist guides from observed guides
    print("Inferring the true guides from observed guides")

    infer_true_reporters_p = partial(infer_true_reporters,
            whitelist_guide_reporter_df=whitelist_guide_reporter_df,
            encoded_whitelist_guide_sequences_series=encoded_whitelist_guide_sequences_series,
            encoded_whitelist_barcodes_series=encoded_whitelist_barcodes_series,
            hamming_threshold=hamming_threshold, barcode_hamming_threshold=barcode_hamming_threshold_strict, surrogate_hamming_threshold=surrogate_hamming_threshold_strict, verbose_result=False)

    inferred_true_reporter_sequences = None
    if cores > 1:
        with Pool(cores) as pool:
            inferred_true_reporter_sequences = pool.map(
            infer_true_reporters_p,
            observed_guides_df["observed_sequence"]
           )
    else:
        inferred_true_reporter_sequences = [infer_true_reporters_p(obs_reporter) for obs_reporter in observed_guides_df["observed_sequence"]]
    
    print("Completed inference")

    observed_guides_df["inferred_guides"] = inferred_true_reporter_sequences
    
    '''
        QC
    '''
    print("Retrieving QC tables")
    # QC: Determine the number of guides passed by determining that the result is not an error
    observed_guides_df_passed_inference = observed_guides_df[observed_guides_df["inferred_guides"].apply(lambda guide : type(guide) != GuideCountError)]
    # QC: Calculate number of guides that were unassigned
    observed_guides_df_no_match = observed_guides_df[observed_guides_df["inferred_guides"].apply(lambda guide : guide == GuideCountError.NO_MATCH)]
    # QC: Calculate number of guides with multiple inferred guides
    observed_guides_df_multiple_match = observed_guides_df[observed_guides_df["inferred_guides"].apply(lambda guide : guide == GuideCountError.MULTIPLE_MATCH)]
    # QC: Calculate percent mapped
    percent_mapped = observed_guides_df_passed_inference["observed_counts"].sum()/observed_guides_df["observed_counts"].sum()
    
    print("Get whitelist reporter counts")
    # Retrieve the observed sequences that were mapped and set the inferred guides
    inferred_guide_sequence_counter = observed_guides_df_passed_inference.groupby("inferred_guides")["observed_counts"].sum()

    whitelist_guide_reporter_counts = whitelist_guide_reporter_df.apply(lambda reporter: inferred_guide_sequence_counter.get(tuple(reporter), 0), axis=1)
    
    
    multi_index = pd.MultiIndex.from_arrays([whitelist_guide_reporter_df['protospacer'], whitelist_guide_reporter_df['surrogate'], whitelist_guide_reporter_df['barcode']], names=['protospacer', 'surrogate', 'barcode'])
    whitelist_guide_reporter_counts.index = multi_index
    

    qc_dict = {"guide_sequences_unassigned_counts":observed_guides_df_no_match["observed_counts"].sum(), "guide_sequences_multiple_counts": observed_guides_df_multiple_match["observed_counts"].sum(), "total_guide_counts": observed_guides_df["observed_counts"].sum(), "percent_mapped": percent_mapped}
    
    return whitelist_guide_reporter_counts, observed_guides_df, qc_dict



@typechecked
def get_whitelist_reporter_counts_with_umi(observed_guide_reporter_umi_counts: defaultdict(Tuple[str,Optional[str],Optional[str]], Counter[Optional[str]]), whitelist_guide_reporter_df: pd.DataFrame, contains_surrogate:bool = False, contains_barcode:bool = False, contains_umi:bool = False, protospacer_hamming_threshold_strict: Optional[int] = 7, surrogate_hamming_threshold_strict: Optional[int] = 2, barcode_hamming_threshold_strict: Optional[int] = 2, cores: int=1):

    whitelist_guide_reporter_df.index = whitelist_guide_reporter_df["protospacer"]
    observed_guide_reporter_sequences_df = pd.DataFrame("observed_sequence": observed_guide_reporter_umi_counts.keys())
    '''
        Commenting this out, since it may be redundant to just transfer the data to a dataframe with reduced information, when the data is already in a good structure.
    '''
    #observed_guides_df = None
    #if isinstance(observed_guide_reporter_umi_counts, defaultdict):
    #    total_counts = [sum(umi_counter) for umi_counter in observed_guide_reporter_umi_counts.values()]
    #    umi_counts = [len(umi_counter) for umi_counter in observed_guide_reporter_umi_counts.values()]
    #else:
    #    observed_guides_df = pd.DataFrame({"observed_sequence":observed_guide_reporter_umi_counts.keys(), 
    #                                       "observed_total_counts":observed_guide_reporter_umi_counts.values(),
    #                                       "observed_umi_counts": observed_guide_reporter_umi_counts.values()
    #                                       })
    
    # Create observed guide DF to contain all information
    # The "observed sequence" column represents all *unique* observed sequences that may have self-edits/errors. The "observed counts" column represents the count of each observed count.
    #observed_guides_df = pd.DataFrame({"observed_sequence":observed_guide_reporter_umi_counts.keys(), "observed_counts":observed_guide_reporter_umi_counts.values()})
 
    # Get the hamming distance threshold. The hamming distance must be below this threshold to assign an observed guide to a whitelist guide.

    encoded_whitelist_protospacer_sequences_series = sequence_encoding.encode_guide_series(whitelist_guide_reporter_df["protospacer"])
    if contains_surrogate:
        encoded_whitelist_surrogate_sequences_series = sequence_encoding.encode_guide_series(whitelist_guide_reporter_df["surrogate"])
    if contains_barcode:
        encoded_whitelist_barcodes_series = sequence_encoding.encode_guide_series(whitelist_guide_reporter_df["barcode"])
    
    protospacer_hamming_threshold_dynamic = False
    if protospacer_hamming_threshold_strict is  None:
        protospacer_hamming_threshold: int = determine_hamming_threshold(whitelist_guide_reporter_df["protospacer"], encoded_whitelist_protospacer_sequences_series, sample_count = 100, quantile = 0.05)
        protospacer_hamming_threshold_dynamic = True
    else:
        protospacer_hamming_threshold: int = protospacer_hamming_threshold_strict
    print("Protospacer hamming threshold is " + str(protospacer_hamming_threshold))

    if contains_surrogate:
        surrogate_hamming_threshold_dynamic = False
        if surrogate_hamming_threshold_strict is  None:
            surrogate_hamming_threshold: int = determine_hamming_threshold(whitelist_guide_reporter_df["surrogate"], encoded_whitelist_surrogate_sequences_series, sample_count = 100, quantile = 0.05)
            surrogate_hamming_threshold_dynamic = True
        else:
            surrogate_hamming_threshold: int = surrogate_hamming_threshold_strict
        print("Surrogate hamming threshold is " + str(surrogate_hamming_threshold))

    if contains_barcode:
        barcode_hamming_threshold_dynamic = False
        if barcode_hamming_threshold_strict is  None:
            barcode_hamming_threshold: int = determine_hamming_threshold(whitelist_guide_reporter_df["barcode"], encoded_whitelist_barcode_sequences_series, sample_count = 100, quantile = 0.05)
            barcode_hamming_threshold_dynamic = True
        else:
            barcode_hamming_threshold: int = barcode_hamming_threshold_strict
        print("Barcode hamming threshold is " + str(barcode_hamming_threshold))

    # Infer whitelist guides from observed guides
    print("Inferring the true guides from observed guides")
    infer_true_reporters_p = partial(infer_true_reporters,
            whitelist_guide_reporter_df=whitelist_guide_reporter_df,
            encoded_whitelist_guide_sequences_series=encoded_whitelist_protospacer_sequences_series,
            encoded_whitelist_barcodes_series=encoded_whitelist_barcodes_series,
            hamming_threshold=protospacer_hamming_threshold, barcode_hamming_threshold=barcode_hamming_threshold_strict, surrogate_hamming_threshold=surrogate_hamming_threshold_strict, verbose_result=False)

    inferred_true_reporter_sequences = None
    if cores > 1:
        with Pool(cores) as pool:
            inferred_true_reporter_sequences = pool.map(
            infer_true_reporters_p,
            observed_guides_df["observed_sequence"]
           )
    else:
        inferred_true_reporter_sequences = [infer_true_reporters_p(obs_reporter) for obs_reporter in observed_guides_df["observed_sequence"]]
    
    print("Completed inference")

    observed_guides_df["inferred_guides"] = inferred_true_reporter_sequences
    
    '''
        QC
    '''
    print("Retrieving QC tables")
    # QC: Determine the number of guides passed by determining that the result is not an error
    observed_guides_df_passed_inference = observed_guides_df[observed_guides_df["inferred_guides"].apply(lambda guide : type(guide) != GuideCountError)]
    # QC: Calculate number of guides that were unassigned
    observed_guides_df_no_match = observed_guides_df[observed_guides_df["inferred_guides"].apply(lambda guide : guide == GuideCountError.NO_MATCH)]
    # QC: Calculate number of guides with multiple inferred guides
    observed_guides_df_multiple_match = observed_guides_df[observed_guides_df["inferred_guides"].apply(lambda guide : guide == GuideCountError.MULTIPLE_MATCH)]
    # QC: Calculate percent mapped
    percent_mapped = observed_guides_df_passed_inference["observed_counts"].sum()/observed_guides_df["observed_counts"].sum()
    
    print("Get whitelist reporter counts")
    # Retrieve the observed sequences that were mapped and set the inferred guides
    inferred_guide_sequence_counter = observed_guides_df_passed_inference.groupby("inferred_guides")["observed_counts"].sum()

    whitelist_guide_reporter_counts = whitelist_guide_reporter_df.apply(lambda reporter: inferred_guide_sequence_counter.get(tuple(reporter), 0), axis=1)
    
    
    multi_index = pd.MultiIndex.from_arrays([whitelist_guide_reporter_df['protospacer'], whitelist_guide_reporter_df['surrogate'], whitelist_guide_reporter_df['barcode']], names=['protospacer', 'surrogate', 'barcode'])
    whitelist_guide_reporter_counts.index = multi_index
    

    qc_dict = {"guide_sequences_unassigned_counts":observed_guides_df_no_match["observed_counts"].sum(), "guide_sequences_multiple_counts": observed_guides_df_multiple_match["observed_counts"].sum(), "total_guide_counts": observed_guides_df["observed_counts"].sum(), "percent_mapped": percent_mapped}
    
    return whitelist_guide_reporter_counts, observed_guides_df, qc_dict



@typechecked
def infer_whitelist_sequence(observed_guide_reporter_sequence, whitelist_guide_reporter_df: pd.DataFrame, contains_surrogate:bool, contains_barcode:bool, contains_umi:bool, encoded_whitelist_protospacer_sequences_series: np.array, encoded_whitelist_surrogate_sequences_series: Optional[np.array] = None, encoded_whitelist_barcode_sequences_series: Optional[np.array] = None, count_duplicate_mappings: bool = False, protospacer_hamming_threshold: int = 7, surrogate_hamming_threshold: int = 10, barcode_hamming_threshold: int = 2, verbose_result = True):
    # Convert from tuple to labeled pandas series. TODO: May already be in this structure.
    observed_reporter_sequences_indices = ["protospacer"]
    if contains_surrogate:
        observed_reporter_sequences_indices.append("surrogate")
    if contains_barcode:
        observed_reporter_sequences_indices.append("barcode")
    observed_guide_reporter_sequence = pd.Series(observed_guide_reporter_sequence, index=observed_reporter_sequences_indices)
    
    # Check if the protospacer observed sequence is empty since this is a requirement for mapping.
    if (observed_guide_reporter_sequence["protospacer"] is None) or (observed_guide_reporter_sequence["protospacer"] == "None") or (observed_guide_reporter_sequence["protospacer"].strip() == ""):
        if verbose_result:
            return {"Error": GuideCountError.NO_PROTOSPACER_MATCH_MISSING_INFO, "message": f"Observed protospacer is None, Complete parsed sequence tuple: {observed_guide_reporter_sequence}"}
        else:
            return GuideCountError.NO_PROTOSPACER_MATCH_MISSING_INFO
    
     # Check if the observed protospacer length is the same length of the whitelisted guides
    if len(observed_guide_reporter_sequence["protospacer"]) != encoded_whitelist_protospacer_sequences_series.shape[1]: 
        if verbose_result:
            return {"Error": GuideCountError.NO_PROTOSPACER_WITH_SAME_LENGTH, "message": f"Observed protospacer {observed_guide_reporter_sequence['protospacer']} not of correct length: {encoded_whitelist_protospacer_sequences_series.shape[1]}"}
        else:
            return GuideCountError.NO_PROTOSPACER_WITH_SAME_LENGTH
        
        
    # COMMENTED FOR NOW - even if these error out, we still could map based on the guide. I want to do guide-based mapping and surrogate based mapping in one go.
    # if contains_surrogate:
    #     if (observed_guide_reporter_sequence["surrogate"] is None) or (observed_guide_reporter_sequence["surrogate"] == "None") or (observed_guide_reporter_sequence["surrogate"].strip() == ""):
    #         if verbose_result:
    #             return {"Error": GuideCountError.NO_SURROGATE_MATCH_MISSING_INFO, "message": f"Observed surrogate is None, Complete parsed sequence tuple: {observed_guide_reporter_sequence}"}
    #         else:
    #             return GuideCountError.NO_SURROGATE_MATCH_MISSING_INFO
        
    #     if len(observed_guide_reporter_sequence["surrogate"]) != encoded_whitelist_surrogate_sequences_series.shape[1]: 
    #         if verbose_result:
    #             return {"Error": GuideCountError.NO_SURROGATE_SAME_LENGTH, "message": f"Observed surrogate {observed_guide_reporter_sequence['surrogate']} not of correct length: {encoded_whitelist_surrogate_sequences_series.shape[1]}"}
    #         else:
    #             return GuideCountError.NO_SURROGATE_SAME_LENGTH   
            
    # if contains_barcode:
    #     if (observed_guide_reporter_sequence["barcode"] is None) or (observed_guide_reporter_sequence["barcode"] == "None") or (observed_guide_reporter_sequence["barcode"].strip() == ""):
    #         if verbose_result:
    #             return {"Error": GuideCountError.NO_BARCODE_MATCH_MISSING_INFO, "message": f"Observed barcode is None, Complete parsed sequence tuple: {observed_guide_reporter_sequence}"}
    #         else:
    #             return GuideCountError.NO_BARCODE_MATCH_MISSING_INFO
            
    #     if len(observed_guide_reporter_sequence["barcode"]) != encoded_whitelist_barcode_sequences_series.shape[1]: 
    #         if verbose_result:
    #             return {"Error": GuideCountError.NO_BARCODE_SAME_LENGTH, "message": f"Observed barcode {observed_guide_reporter_sequence['barcode']} not of correct length: {encoded_whitelist_barcode_sequences_series.shape[1]}"}
    #         else:
    #             return GuideCountError.NO_BARCODE_SAME_LENGTH   
            
    
    # PERFORM MAPPPING ON BOTH THE PROTOSPACER-ONLY, FULL, AND SURROGATE MISMAPPED-PROTOSPACER.

    # Determine if there are exact matches, hopefully just a single match
    # whitelist_protospacer_exact_matches = (whitelist_guide_reporter_df["protospacer"]==observed_guide_reporter_sequence["protospacer"])
    # if contains_surrogate:
    #     whitelist_surrogate_exact_matches = (whitelist_guide_reporter_df["surrogate"]==observed_guide_reporter_sequence["surrogate"])    
    # if contains_barcode:
    #     whitelist_barcode_exact_matches = (whitelist_guide_reporter_df["barcode"]==observed_guide_reporter_sequence["barcode"])    

    match_result = {
        "protospacer_match": None, # DONE
        "protospacer_match_surrogate_match_barcode_match": None, # DONE
        "protospacer_match_surrogate_match": None,
        "protospacer_match_barcode_match": None, # DONE
        "protospacer_mismatch_surrogate_match_barcode_match": None, # DONE
        "protospacer_mismatch_surrogate_match": None, # DONE
        
    }

    # NOTE: Probably won't do exact match lookup separately, will just rely on the hamming-based approach to find exact matches since it is fast.
    # def assign_exact_matches(match_key, whitelist_guide_reporter_df_reporter_match):
    #     if whitelist_guide_reporter_df_reporter_match.shape[0] == 1:
    #         match_result[match_key] = whitelist_guide_reporter_df_reporter_match
    #     elif whitelist_guide_reporter_df_reporter_match.shape[0] > 1:
    #         if count_duplicate_mappings:
    #             match_result[match_key] = whitelist_guide_reporter_df_reporter_match
    #         else:
    #             if verbose_result:
    #                 match_result[match_key] = {"Error": GuideCountError.MULTIPLE_MATCH_EXACT,
    #                 "exact_match": True, 
    #                 "num_matches": whitelist_guide_reporter_df_reporter_match.shape[0],
    #                 "matches": whitelist_guide_reporter_df_reporter_match,
    #                 "hamming_min": 0}
    #             else:
    #                 match_result[match_key] = GuideCountError.MULTIPLE_MATCH_EXACT


    # whitelist_guide_reporter_df_reporter_protospacer_match = whitelist_guide_reporter_df[whitelist_protospacer_exact_matches]
    # assign_exact_matches("protospacer_match", whitelist_guide_reporter_df_reporter_protospacer_match)
    # if contains_surrogate:
    #     if contains_barcode:
    #         whitelist_guide_reporter_df_reporter_protospacer_surrogate_barcode_match = whitelist_guide_reporter_df[whitelist_protospacer_exact_matches & whitelist_surrogate_exact_matches & whitelist_barcode_exact_matches]
    #         assign_exact_matches("protospacer_match_surrogate_match_barcode_match", whitelist_guide_reporter_df_reporter_protospacer_surrogate_barcode_match)

    #         whitelist_guide_reporter_df_reporter_nonprotospacer_surrogate_barcode_match = whitelist_guide_reporter_df[(~whitelist_protospacer_exact_matches) & whitelist_surrogate_exact_matches & whitelist_barcode_exact_matches]
    #         assign_exact_matches("protospacer_mismatch_surrogate_match_barcode_match", whitelist_guide_reporter_df_reporter_nonprotospacer_surrogate_barcode_match)
    #     else:
    #         whitelist_guide_reporter_df_reporter_protospacer_surrogate_match = whitelist_guide_reporter_df[whitelist_protospacer_exact_matches & whitelist_surrogate_exact_matches]
    #         assign_exact_matches("protospacer_match_surrogate_match", whitelist_guide_reporter_df_reporter_protospacer_surrogate_match)

    #         whitelist_guide_reporter_df_reporter_nonprotospacer_surrogate_match = whitelist_guide_reporter_df[(~whitelist_protospacer_exact_matches) & whitelist_surrogate_exact_matches]
    #         assign_exact_matches("protospacer_mismatch_surrogate_match_barcode_match", whitelist_guide_reporter_df_reporter_nonprotospacer_surrogate_match)
    # else:
    #     if contains_barcode:
    #         whitelist_guide_reporter_df_reporter_protospacer_barcode_match = whitelist_guide_reporter_df[whitelist_protospacer_exact_matches & whitelist_barcode_exact_matches]
    #         assign_exact_matches("protospacer_match_barcode_match", whitelist_guide_reporter_df_reporter_protospacer_barcode_match)

    


    # TODO: MAKE SURE TO CHECK THE LENGTHS OFF THE MATCHED BARCODES AND SURROGATES! ALSO, AS COMMENTED ABOVE, CHECK IF THEY ARE NONE/MISSING AS WELL! Might need a helper function to stay modular.
    # GET THE BARCODE-ONLY MATCHES
    if contains_barcode:
        observed_barcode_sequence_encoded = sequence_encoding.encode_DNA_base_vectorized(sequence_encoding.numpify_string_vectorized(observed_guide_reporter_sequence["barcode"]))  # Encode the observed barcode
        observed_barcode_sequence_dists = sequence_encoding.retrieve_hamming_distance_whitelist(observed_barcode_sequence_encoded, encoded_whitelist_barcode_sequences_series) # Retrieve hamming distance with whitelist barcode
        observed_barcode_sequence_dists_min = observed_barcode_sequence_dists.min() # Get the barcode with the minimum  hamming distance
        barcode_hamming_threshold_met = (observed_barcode_sequence_dists_min < barcode_hamming_threshold)
        if barcode_hamming_threshold_met: 
            barcode_matches_indices = np.where(observed_barcode_sequence_dists == observed_barcode_sequence_dists_min)[0] # Get the indices of ALL barcode matches
            whitelist_guide_reporter_df_barcode_match = whitelist_guide_reporter_df.iloc[barcode_matches_indices] # Get the whitelist reporters with the matched barcode(s)
            encoded_whitelist_protospacer_sequences_series_barcode_match = encoded_whitelist_protospacer_sequences_series[barcode_matches_indices] # Subset the protospacer encodings with the barcode matches for later
            if contains_surrogate:
                encoded_whitelist_surrogate_sequences_series_barcode_match = encoded_whitelist_surrogate_sequences_series[barcode_matches_indices] # Subset the surrogate encodings with the barcode matches for later
        else: # If the barcode surpasses the threshold, fail the read due to no barcode. 
            if verbose_result:
                error_result = {"Error": GuideCountError.NO_MATCH_BARCODE_HAMMING_THRESHOLD, "barcode_hamming_min": observed_barcode_sequence_dists_min, "message": f"No barcode below threshold {barcode_hamming_threshold}"}
            else:
                error_result = GuideCountError.NO_MATCH_BARCODE_HAMMING_THRESHOLD
            # Error, barcode hamming threshold not met. For any result requiring barcode, set error
            match_result["protospacer_match_surrogate_match_barcode_match"] = error_result
            match_result["protospacer_mismatch_surrogate_match_barcode_match"] = error_result
            match_result["protospacer_match_barcode_match"] = error_result
            

    
    # LEFTOFF (START HERE, THEN GO TO NEXT LEFTOFF): THE MISMATCH LOGIC IS DEFINITELY WRONG, IT WILL RETURN NO RESULTS AT ALL DUE TO HOW THE INDEXING WORKS (SUBSETTING-BASED), LOGIC FLOW MAY NOT BE AS CLEAR AS WHAT IS THE SUBSETTING DONE BLOW BELOW
    # PREPARE THE PROTOSPACER-ONLYMATCHES
    observed_protospacer_sequence_encoded = sequence_encoding.encode_DNA_base_vectorized(sequence_encoding.numpify_string_vectorized(observed_guide_reporter_sequence["protospacer"]))  # Encode the observed protospacer
    observed_protospacer_sequence_dists = sequence_encoding.retrieve_hamming_distance_whitelist(observed_protospacer_sequence_encoded, encoded_whitelist_protospacer_sequences_series) # Hamming distance among all whitelist protospacers
    observed_protospacer_sequence_dists_min = observed_protospacer_sequence_dists.min() # Get minimum hamming distance
    protospacer_hamming_threshold_met = (observed_protospacer_sequence_dists_min < protospacer_hamming_threshold)
    if protospacer_hamming_threshold_met: 
        # SET THE PROTOSPACER-ONLYMATCHES
        protospacer_matches_indices = np.where(observed_protospacer_sequence_dists == observed_protospacer_sequence_dists_min)[0]
        whitelist_guide_reporter_df_hamming_protospacer_match = whitelist_guide_reporter_df.iloc[protospacer_matches_indices] # Get all whitelisted guides with the minimum hamming distance (could be multiple)
        encoded_whitelist_protospacer_sequences_series_protospacer_match = encoded_whitelist_protospacer_sequences_series[protospacer_matches_indices] # Subset the protospacer encodings with the protospacer matches for later
        match_result["protospacer_match"] = whitelist_guide_reporter_df_hamming_protospacer_match
        
        # GET THE PROTOSPACER-NONMATCHES
        # protospacer_nonmatches_indices = np.where(observed_protospacer_sequence_dists != observed_protospacer_sequence_dists_min)[0]
        # whitelist_guide_reporter_df_hamming_protospacer_nonmatch = whitelist_guide_reporter_df.iloc[protospacer_nonmatches_indices] 
        # encoded_whitelist_protospacer_sequences_series_protospacer_nonmatch = encoded_whitelist_protospacer_sequences_series[protospacer_nonmatches_indices] # Subset the protospacer encodings with the protospacer matches for later
        if contains_barcode: 
            if barcode_hamming_threshold_met:
                
                # PREPARE PROTOSPACER-MATCH, BARCODE-MATCH
                observed_protospacer_sequence_dists_barcode_match = sequence_encoding.retrieve_hamming_distance_whitelist(observed_protospacer_sequence_encoded, encoded_whitelist_protospacer_sequences_series_barcode_match) # Hamming distance among barcode-match protospacers
                observed_protospacer_sequence_dists_barcode_match_min = observed_protospacer_sequence_dists_barcode_match.min() # Get minimum hamming distance
                barcode_match_protospacer_hamming_threshold_met = (observed_protospacer_sequence_dists_barcode_match_min < protospacer_hamming_threshold)
                if barcode_match_protospacer_hamming_threshold_met: 
                    
                    # SET PROTOSPACER-MATCH, BARCODE-MATCH
                    barcode_match_protospacer_match_indices = np.where(observed_protospacer_sequence_dists_barcode_match == observed_protospacer_sequence_dists_barcode_match_min)[0]
                    whitelist_guide_reporter_df_hamming_barcode_match_protospacer_match = whitelist_guide_reporter_df_barcode_match.iloc[barcode_match_protospacer_match_indices] # Get all whitelisted guides with the minimum hamming distance (could be multiple)
                    match_result["protospacer_match_barcode_match"] = whitelist_guide_reporter_df_hamming_barcode_match_protospacer_match

                    # GET PROTOSPACER-NONMATCH, BARCODE-MATCH
                    # barcode_match_protospacer_nonmatch_indices = np.where(observed_protospacer_sequence_dists_barcode_match != observed_protospacer_sequence_dists_barcode_match_min)[0]
                    # whitelist_guide_reporter_df_hamming_barcode_match_protospacer_nonmatch = whitelist_guide_reporter_df_barcode_match.iloc[barcode_match_protospacer_nonmatch_indices] # Get all whitelisted guides with the minimum hamming distance (could be multiple)
                    if contains_surrogate: 
                        # PREPARE PROTOSPACER-MATCH, BARCODE-MATCH, SURROGATE-MATCH
                        encoded_whitelist_surrogate_sequences_series_barcode_match_protospacer_match = encoded_whitelist_surrogate_sequences_series_barcode_match[barcode_match_protospacer_match_indices] # Subset the surrogate encodings with the protospacer and encoding matches for later
                        observed_surrogate_sequence_encoded = sequence_encoding.encode_DNA_base_vectorized(sequence_encoding.numpify_string_vectorized(observed_guide_reporter_sequence["surrogate"]))  # Encode the observed protospacer
                        observed_surrogate_sequence_dists_barcode_match_protospacer_match = sequence_encoding.retrieve_hamming_distance_whitelist(observed_surrogate_sequence_encoded, encoded_whitelist_surrogate_sequences_series_barcode_match_protospacer_match) # Hamming distance among all whitelist sub-selected surrogates
                        observed_surrogate_sequence_dists_barcode_match_protospacer_match_min = observed_surrogate_sequence_dists_barcode_match_protospacer_match.min()
                        barcode_match_protospacer_match_surrogate_hamming_threshold_met = (observed_surrogate_sequence_dists_barcode_match_protospacer_match_min < surrogate_hamming_threshold)
                        if barcode_match_protospacer_match_surrogate_hamming_threshold_met:
                            # SET PROTOSPACER-MATCH, BARCODE-MATCH, SURROGATE-MATCH
                            barcode_match_protospacer_match_surrogate_match_indices = np.where(observed_surrogate_sequence_dists_barcode_match_protospacer_match == observed_surrogate_sequence_dists_barcode_match_protospacer_match_min)[0]
                            whitelist_guide_reporter_df_hamming_barcode_match_protospacer_match_surrogate_match = whitelist_guide_reporter_df_hamming_barcode_match_protospacer_match.iloc[barcode_match_protospacer_match_surrogate_match_indices] # Get all whitelisted guides with the minimum hamming distance (could be multiple)
                            match_result["protospacer_match_barcode_match_surrogate_match"] = whitelist_guide_reporter_df_hamming_barcode_match_protospacer_match_surrogate_match
                            pass
                        else:
                            if verbose_result:
                                error_result = {"Error": GuideCountError.NO_MATCH_SURROGATE_HAMMING_THRESHOLD, "surrogate_hamming_min": observed_surrogate_sequence_dists_barcode_match_protospacer_match_min}
                            else:
                                error_result = GuideCountError.NO_MATCH_SURROGATE_HAMMING_THRESHOLD
                            match_result["protospacer_match_barcode_match_surrogate_match"] = error_result

                        # PREPARE PROTOSPACER-NONMATCH, BARCODE-MATCH, SURROGATE-MATCH
                        # encoded_whitelist_surrogate_sequences_series_barcode_match_protospacer_nonmatch = encoded_whitelist_surrogate_sequences_series_barcode_match[barcode_match_protospacer_nonmatch_indices] # Subset the surrogate encodings with the protospacer and encoding matches for later
                        # observed_surrogate_sequence_dists_barcode_match_protospacer_nonmatch = sequence_encoding.retrieve_hamming_distance_whitelist(observed_surrogate_sequence_encoded, encoded_whitelist_surrogate_sequences_series_barcode_match_protospacer_nonmatch) # Hamming distance among all whitelist sub-selected surrogates
                        # observed_surrogate_sequence_dists_barcode_match_protospacer_nonmatch_min = observed_surrogate_sequence_dists_barcode_match_protospacer_nonmatch.min()
                        # barcode_match_protospacer_nonmatch_surrogate_hamming_threshold_met = (observed_surrogate_sequence_dists_barcode_match_protospacer_nonmatch_min < surrogate_hamming_threshold)
                        # if barcode_match_protospacer_nonmatch_surrogate_hamming_threshold_met:
                        #     # SET PROTOSPACER-MATCH, BARCODE-MATCH, SURROGATE-MATCH
                        #     barcode_match_protospacer_nonmatch_surrogate_match_indices = np.where(observed_surrogate_sequence_dists_barcode_match_protospacer_nonmatch == observed_surrogate_sequence_dists_barcode_match_protospacer_nonmatch_min)[0]
                        #     whitelist_guide_reporter_df_hamming_barcode_match_protospacer_nonmatch_surrogate_match = whitelist_guide_reporter_df_hamming_barcode_match_protospacer_nonmatch.iloc[barcode_match_protospacer_nonmatch_surrogate_match_indices] # Get all whitelisted guides with the minimum hamming distance (could be multiple)
                        #     match_result["protospacer_mismatch_surrogate_match_barcode_match"] = whitelist_guide_reporter_df_hamming_barcode_match_protospacer_nonmatch_surrogate_match
                        # else:
                        #     if verbose_result:
                        #         error_result = {"Error": GuideCountError.NO_MATCH_SURROGATE_HAMMING_THRESHOLD, "surrogate_hamming_min": observed_surrogate_sequence_dists_barcode_match_protospacer_nonmatch_min}
                        #     else:
                        #         error_result = GuideCountError.NO_MATCH_SURROGATE_HAMMING_THRESHOLD
                        #     match_result["protospacer_mismatch_surrogate_match_barcode_match"] = error_result
                    
                else:
                    if verbose_result:
                        error_result = {"Error": GuideCountError.NO_MATCH_PROTOSPACER_HAMMING_THRESHOLD, "protospacer_hamming_min": observed_protospacer_sequence_dists_barcode_match_min}
                    else:
                        error_result = GuideCountError.NO_MATCH_PROTOSPACER_HAMMING_THRESHOLD
                    match_result["protospacer_match_barcode_match"] = error_result
        
        # PROTOSPACER-MATCH, SURROGATE-MATCH, NO BARCODE
        if contains_surrogate:
            encoded_whitelist_surrogate_sequences_series_protospacer_match = encoded_whitelist_surrogate_sequences_series[protospacer_matches_indices]
            
            observed_surrogate_sequence_encoded = sequence_encoding.encode_DNA_base_vectorized(sequence_encoding.numpify_string_vectorized(observed_guide_reporter_sequence["surrogate"]))  # Encode the observed protospacer
            observed_surrogate_sequence_dists_protospacer_match = sequence_encoding.retrieve_hamming_distance_whitelist(observed_surrogate_sequence_encoded, encoded_whitelist_surrogate_sequences_series_protospacer_match) # TODO: This should be on the protospacer-match array
            observed_surrogate_sequence_dists_protospacer_match_min = observed_surrogate_sequence_dists_protospacer_match.min()
            protospacer_match_surrogate_hamming_threshold_met = (observed_surrogate_sequence_dists_protospacer_match_min < surrogate_hamming_threshold)
            if protospacer_match_surrogate_hamming_threshold_met:
                protospacer_match_surrogate_match_indices = np.where(observed_surrogate_sequence_dists_protospacer_match == observed_surrogate_sequence_dists_protospacer_match_min)[0]
                whitelist_guide_reporter_df_hamming_protospacer_match_surrogate_match = whitelist_guide_reporter_df_hamming_protospacer_match.iloc[protospacer_match_surrogate_match_indices] # Get all whitelisted guides with the minimum hamming distance (could be multiple)
                match_result["protospacer_match_surrogate_match"] = whitelist_guide_reporter_df_hamming_protospacer_match_surrogate_match
            else:
                if verbose_result:
                    error_result = {"Error": GuideCountError.NO_MATCH_SURROGATE_HAMMING_THRESHOLD, "surrogate_hamming_min": observed_surrogate_sequence_dists_min}
                else:
                    error_result = GuideCountError.NO_MATCH_SURROGATE_HAMMING_THRESHOLD
                # Error, protospacer hamming threshold not met. For any result requiring protospacer, set error
                match_result["protospacer_match_surrogate_match"] = error_result







        # GET THE PROTOSPACER-SURROGATE MISMATCHES NOW (if surrogate is available)
        if contains_surrogate:
            if contains_barcode:
                if barcode_hamming_threshold_met:
                    #HEREEE
                    # PREPARE SURROGATE-MATCH, BARCODE-MATCH
                    observed_surrogate_sequence_encoded = sequence_encoding.encode_DNA_base_vectorized(sequence_encoding.numpify_string_vectorized(observed_guide_reporter_sequence["surrogate"]))  # Encode the observed protospacer
                    observed_surrogate_sequence_dists_barcode_match = sequence_encoding.retrieve_hamming_distance_whitelist(observed_surrogate_sequence_encoded, encoded_whitelist_surrogate_sequences_series_barcode_match) # Hamming distance among barcode-match protospacers
                    observed_surrogate_sequence_dists_barcode_match_min = observed_surrogate_sequence_dists_barcode_match.min() # Get minimum hamming distance
                    barcode_match_surrogate_hamming_threshold_met = (observed_surrogate_sequence_dists_barcode_match_min < surrogate_hamming_threshold)
                    if barcode_match_surrogate_hamming_threshold_met:
                        # SET SURROGATE-MATCH, BARCODE-MATCH
                        barcode_match_surrogate_match_indices = np.where(observed_surrogate_sequence_dists_barcode_match == observed_surrogate_sequence_dists_barcode_match_min)[0]
                        whitelist_guide_reporter_df_hamming_barcode_match_surrogate_match = whitelist_guide_reporter_df_barcode_match.iloc[barcode_match_surrogate_match_indices] # Get all whitelisted guides with the minimum hamming distance (could be multiple)

                        # See if there are identical matches between protospacer-only matches and surrogate-only matches
                        whitelist_guide_reporter_df_hamming_barcode_match_surrogate_match_protospacer_match = pd.merge(whitelist_guide_reporter_df_hamming_barcode_match_surrogate_match, whitelist_guide_reporter_df_hamming_protospacer_match, how='inner')
                        match_result["protospacer_mismatch_surrogate_match_barcode_match"] = {
                            "mismatched": whitelist_guide_reporter_df_hamming_barcode_match_surrogate_match_protospacer_match.empty, 
                            "surrogate_barcode_matches": whitelist_guide_reporter_df_hamming_barcode_match_surrogate_match_protospacer_match, 
                            "protospacer_matches": whitelist_guide_reporter_df_hamming_protospacer_match,
                            "protospacer_surrogate_barcode_matches": whitelist_guide_reporter_df_hamming_barcode_match_surrogate_match_protospacer_match
                            }
                    else:
                        if verbose_result:
                            error_result = {"Error": GuideCountError.NO_MATCH_SURROGATE_HAMMING_THRESHOLD, "surrogate_hamming_min": observed_surrogate_sequence_dists_barcode_match_min}
                        else:
                            error_result = GuideCountError.NO_MATCH_SURROGATE_HAMMING_THRESHOLD
                        # Error, protospacer hamming threshold not met. For any result requiring protospacer, set error
                        match_result["protospacer_mismatch_surrogate_match_barcode_match"] = error_result
            
            
            observed_surrogate_sequence_encoded = sequence_encoding.encode_DNA_base_vectorized(sequence_encoding.numpify_string_vectorized(observed_guide_reporter_sequence["surrogate"]))  # Encode the observed protospacer
            observed_surrogate_sequence_dists = sequence_encoding.retrieve_hamming_distance_whitelist(observed_surrogate_sequence_encoded, encoded_whitelist_surrogate_sequences_series) # Hamming distance among  all whitelist surrogates
            observed_surrogate_sequence_dists_min = observed_surrogate_sequence_dists.min()
            surrogate_hamming_threshold_met = (observed_surrogate_sequence_dists_min < surrogate_hamming_threshold)
            if surrogate_hamming_threshold_met:
                surrogate_match_indices = np.where(observed_surrogate_sequence_dists == observed_surrogate_sequence_dists_min)[0]
                whitelist_guide_reporter_df_hamming_surrogate_match = whitelist_guide_reporter_df.iloc[surrogate_match_indices] # Get all whitelisted guides with the minimum hamming distance (could be multiple)

                # See if there are identical matches between protospacer-only matches and surrogate-only matches
                whitelist_guide_reporter_df_hamming_surrogate_match_protospacer_match = pd.merge(whitelist_guide_reporter_df_hamming_surrogate_match, whitelist_guide_reporter_df_hamming_protospacer_match, how='inner')
                match_result["protospacer_mismatch_surrogate_match"] = {
                    "mismatched": whitelist_guide_reporter_df_hamming_surrogate_match_protospacer_match.empty, 
                    "surrogate_matches": whitelist_guide_reporter_df_hamming_surrogate_match, 
                    "protospacer_matches": whitelist_guide_reporter_df_hamming_protospacer_match,
                    "protospacer_surrogate_matches": whitelist_guide_reporter_df_hamming_surrogate_match_protospacer_match
                    }
            else:
                if verbose_result:
                    error_result = {"Error": GuideCountError.NO_MATCH_SURROGATE_HAMMING_THRESHOLD, "surrogate_hamming_min": observed_surrogate_sequence_dists_min}
                else:
                    error_result = GuideCountError.NO_MATCH_SURROGATE_HAMMING_THRESHOLD
                # Error, protospacer hamming threshold not met. For any result requiring protospacer, set error
                match_result["protospacer_mismatch_surrogate_match"] = error_result
            
    else: # If the barcode surpasses the threshold, fail the read due to no barcode. 
        if verbose_result:
            error_result = {"Error": GuideCountError.NO_MATCH_PROTOSPACER_HAMMING_THRESHOLD, "protospacer_hamming_min": observed_protospacer_sequence_dists_min}
        else:
            error_result = GuideCountError.NO_MATCH_PROTOSPACER_HAMMING_THRESHOLD
        # Error, protospacer hamming threshold not met. For any result requiring protospacer, set error
        match_result["protospacer_match"] = error_result
        match_result["protospacer_match_surrogate_match_barcode_match"] = error_result
        match_result["protospacer_match_surrogate_match"] = error_result
        match_result["protospacer_match_barcode_match"] = error_result
        match_result["protospacer_mismatch_surrogate_match_barcode_match"] = error_result
        match_result["protospacer_mismatch_surrogate_match"] = error_result

    

    # LEFTOFF HERE - just finished draft, commit, then cleanup code and fix errors.
    return match_result


###
### MAIN FUNCTIONS
###
'''
    Take in input FASTQ filename, and a set of whitelisted guide sequences
'''
# DEPRECATED use the all-in-one UMI-tools package
@typechecked
def get_whitelist_guide_counts_from_raw_fastq(whitelist_guide_sequences_series: pd.Series, fastq_fn: str, hamming_threshold_strict: int = 3, hamming_threshold_dynamic: bool = False, parse_left_flank: bool = True, parse_flank_sequence: Union[None, str] = None, cores: int=1):
    # Retrieve all observed guide sequences
    print("Retrieving FASTQ guide sequences and counting: " + fastq_fn)
    observed_guide_sequences_counts: Counter[str] = guide_raw_fastq_parsing.get_raw_fastq_observed_sequence_counts(fastq_fn, parse_left_flank=parse_left_flank, parse_flank_sequence=parse_flank_sequence, cores=cores)
    
    return get_whitelist_guide_counts(observed_guide_sequences_counts, whitelist_guide_sequences_series, hamming_threshold_strict, hamming_threshold_dynamic, cores)

'''
    Take in input FASTQ filename, and a set of whitelisted guide sequences
'''
@typechecked
def get_whitelist_reporter_counts_from_reporter_tsv(whitelist_guide_reporter_df: pd.DataFrame, reporter_tsv_fn: str, surrogate_hamming_threshold_strict: int = 10, barcode_hamming_threshold_strict: int = 2, hamming_threshold_strict: int = 7, hamming_threshold_dynamic: bool = False, cores: int=1):
    observed_guide_reporter_counts: Union[Counter[Tuple[str,str,str]], Counter[str]] = reporter_tsv_parsing.get_reporter_tsv_observed_sequence_counts(reporter_tsv_fn, include_surrogate = True, cores=cores)

    return get_whitelist_reporter_counts(observed_guide_reporters_counts=observed_guide_reporter_counts, whitelist_guide_reporter_df=whitelist_guide_reporter_df, surrogate_hamming_threshold_strict=surrogate_hamming_threshold_strict, barcode_hamming_threshold_strict=barcode_hamming_threshold_strict, hamming_threshold_strict=hamming_threshold_strict, hamming_threshold_dynamic=hamming_threshold_dynamic, cores=cores)


@typechecked
def get_whitelist_reporter_counts_from_umitools_output(whitelist_guide_reporter_df: pd.DataFrame, fastq_r1_fn: str, fastq_r2_fn: str, barcode_pattern_regex: Optional[str] = None, umi_pattern_regex: Optional[str] = None, surrogate_hamming_threshold_strict: int = 10, barcode_hamming_threshold_strict: int = 2, hamming_threshold_strict: int = 7, hamming_threshold_dynamic: bool = False, cores: int=1):
    observed_guide_reporter_umi_counts = reporter_tsv_parsing.get_umitools_observed_sequence_counts(r1_protospacer_fastq_file=fastq_r1_fn, r2_surrogate_fastq_file=fastq_r2_fn, barcode_pattern_regex=barcode_pattern_regex, umi_pattern_regex=umi_pattern_regex)
    if fastq_r2_fn is None: # ONLY R1
        if barcode_pattern_regex is None: # ONLY R1; NO BARCODE
            if umi_pattern_regex is None: # ONLY R1; NO BARCODE; NO UMI
                observed_guide_reporter_umi_counts: Counter[str] = observed_guide_reporter_umi_counts
                observed_guide_sequences_counts(observed_guide_reporter_umi_counts, whitelist_guide_reporter_df["protospacer"], hamming_threshold_strict, hamming_threshold_dynamic, cores)
            else: # ONLY R1; NO BARCODE; YES UMI
                observed_guide_reporter_umi_counts: defaultdict(str, Counter[str]) = observed_guide_reporter_umi_counts
                raise NotImplementedError("Not implemented")
        else: # ONLY R1; YES BARCODE
            if umi_pattern_regex is None: # ONLY R1; YES BARCODE; NO UMI
                observed_guide_reporter_umi_counts: Counter[Tuple[str, str]] = observed_guide_reporter_umi_counts
                raise NotImplementedError("Not implemented")
            else: # ONLY R1; YES BARCODE; YES UMI
                observed_guide_reporter_umi_counts: defaultdict(Tuple[str, str], Counter[str]) = observed_guide_reporter_umi_counts
                raise NotImplementedError("Not implemented")
    else: # YES R2
        if barcode_pattern_regex is None: # YES R2; NO BARCODE
            if umi_pattern_regex is None: # YES R2; NO BARCODE; NO UMI
                observed_guide_reporter_umi_counts: Counter[Tuple[str, str]] = observed_guide_reporter_umi_counts
                raise NotImplementedError("Not implemented")
            else: # YES R2; NO BARCODE; YES UMI
                observed_guide_reporter_umi_counts: defaultdict(Tuple[str, str], Counter[str]) = observed_guide_reporter_umi_counts
                raise NotImplementedError("Not implemented")
        else: # YES R2; YES BARCODE
            if umi_pattern_regex is None: # YES R2; YES BARCODE; NO UMI
                observed_guide_reporter_umi_counts: Counter[Tuple[str, str, str]] = observed_guide_reporter_umi_counts
                return get_whitelist_reporter_counts(observed_guide_reporters_counts=observed_guide_reporter_umi_counts, whitelist_guide_reporter_df=whitelist_guide_reporter_df, surrogate_hamming_threshold_strict=surrogate_hamming_threshold_strict, barcode_hamming_threshold_strict=barcode_hamming_threshold_strict, hamming_threshold_strict=hamming_threshold_strict, hamming_threshold_dynamic=hamming_threshold_dynamic, cores=cores)
            else: # YES R2; YES BARCODE; YES UMI
                observed_guide_reporter_umi_counts: defaultdict(Tuple[str, str, str], Counter[str]) = observed_guide_reporter_umi_counts
                raise NotImplementedError("Not implemented")



