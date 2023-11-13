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

from . import sequence_encoding
from .models import *

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
            return GuideCountErrorType.NO_PROTOSPACER_WITH_SAME_LENGTH 
        

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
                return {"Error": GuideCountErrorType.NO_MATCH_PROTOSPACER_HAMMING_THRESHOLD, "hamming_min": hamming_min}
            else:
                return GuideCountErrorType.NO_MATCH_PROTOSPACER_HAMMING_THRESHOLD
        
        # If there are multiple guides with the minimum hamming distance, then the guide is ambigious, so no mapping (due to multiple match)
        # NOTE (3/30/2023): This could potentially be an issue with igRNAs maybe?
        elif len(guides_with_hamming_min) > 1:
            if verbose_result:
                return {"Error": GuideCountErrorType.MULTIPLE_MATCH,
                "exact_match": False, 
                "num_matches": len(guides_with_hamming_min),
                "matches": guides_with_hamming_min,
                "hamming_min": hamming_min}
            else:
                return GuideCountErrorType.MULTIPLE_MATCH
        
        # Else if there is 1 guide with the match, then return the match
        else:
            return tuple([str(guides_with_hamming_min.index[0])])
    
    # Else if there are multiple exact match, which should never occur unless the whitelisted guide list is not unique, then return multiple match.
    else:
        
        if verbose_result:
            return {"Error": GuideCountErrorType.MULTIPLE_MATCH_EXACT,
            "exact_match": True, 
            "num_matches": len(whitelist_guide_sequences_series_match),
            "matches": whitelist_guide_sequences_series_match,
            "hamming_min": 0}
        else:
            return GuideCountErrorType.MULTIPLE_MATCH_EXACT
        #raise Exception("Multiple exact matches of the provided whitelisted guides - there are likely duplicates in the provided whitelist, please remove. Observed guide={}, guide matches={}".format(observed_guide_sequence, guide_sequences_series_match)) # NOTE 12/6/22: REMOVED THIS EXCEPTION - another reason is from truncated guides having multiple matches. In production code, just make sure to ensure that the whitelist is the set.

@typechecked
def infer_true_reporters(observed_reporter_sequences, whitelist_guide_reporter_df: pd.DataFrame,
encoded_whitelist_guide_sequences_series, encoded_whitelist_barcodes_series, surrogate_hamming_threshold: int = 10, barcode_hamming_threshold: int = 2, hamming_threshold: int = 7, verbose_result = False):
    for element in observed_reporter_sequences:
        if (element is None) or (element == "None") or (element == ""):
            if verbose_result:
                return {"Error": GuideCountErrorType.NO_MATCH_MISSING_INFO, "message": "Observed protospacer/surrogate/barcode is None"}
            else:
                return GuideCountErrorType.NO_MATCH_MISSING_INFO
    
    observed_reporter_sequences = pd.Series(observed_reporter_sequences, index=["protospacer", "surrogate", "barcode"])
    
    # Check if the observed protospacer length is the same length of the whitelisted guides
    if len(observed_reporter_sequences["protospacer"]) != encoded_whitelist_guide_sequences_series.shape[1]: 
        if verbose_result:
            return {"Error": GuideCountErrorType.NO_PROTOSPACER_WITH_SAME_LENGTH, "message": f"Observed protospacer {observed_reporter_sequences['protospacer']} not of correct length: {encoded_whitelist_guide_sequences_series.shape[1]}"}
        else:
            return GuideCountErrorType.NO_PROTOSPACER_WITH_SAME_LENGTH
    if len(observed_reporter_sequences["barcode"]) != encoded_whitelist_barcodes_series.shape[1]: 
        if verbose_result:
            return {"Error": GuideCountErrorType.NO_BARCODE_WITH_SAME_LENGTH, "message": f"Observed barcode {observed_reporter_sequences['barcode']} not of correct length: {encoded_whitelist_barcodes_series.shape[1]}"}
        else:
            return GuideCountErrorType.NO_BARCODE_WITH_SAME_LENGTH   
            
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
                return {"Error": GuideCountErrorType.NO_MATCH_BARCODE_HAMMING_THRESHOLD, "barcode_hamming_min": barcode_hamming_min, "message": f"No barcode below threshold {barcode_hamming_threshold}"}
            else:
                return GuideCountErrorType.NO_MATCH_BARCODE_HAMMING_THRESHOLD
            
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
                return {"Error": GuideCountErrorType.NO_MATCH_PROTOSPACER_HAMMING_THRESHOLD, "hamming_min": hamming_min}
            else:
                return GuideCountErrorType.NO_MATCH_PROTOSPACER_HAMMING_THRESHOLD
        
        # If there are multiple guides with the minimum hamming distance, then the guide is ambigious, so no mapping (due to multiple match)
        # NOTE (3/30/2023): This could potentially be an issue with igRNAs maybe? (Response 6/30/23) No, the barcode should be able to distinguish between igRNAs
        elif reporters_with_hamming_min_df.shape[0] > 1:
            if verbose_result:
                return {"Error": GuideCountErrorType.MULTIPLE_MATCH,
                "exact_match": False, 
                "num_matches": reporters_with_hamming_min_df.shape[0],
                "matches": reporters_with_hamming_min_df,
                "hamming_min": hamming_min}
            else:
                return GuideCountErrorType.MULTIPLE_MATCH
        
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
                        return {"Error": GuideCountErrorType.NO_MATCH_SURROGATE_HAMMING_THRESHOLD, "surrogate_hamming_distance": surrogate_hamming_distance, "inferred_reporter_sequences": inferred_reporter_sequences}
                    else:
                        return GuideCountErrorType.NO_MATCH_SURROGATE_HAMMING_THRESHOLD
                else:
                    return inferred_reporter_sequences
            else:
                if verbose_result:
                        return {"Error": GuideCountErrorType.NO_MATCH_OBSERVED_SURROGATE_SHORTER, "inferred_reporter_sequences": inferred_reporter_sequences}
                else:
                    return GuideCountErrorType.NO_MATCH_OBSERVED_SURROGATE_SHORTER
    
    # Else if there are multiple exact match, which should never occur unless the whitelisted guide list is not unique, then return multiple match.
    else:
        if verbose_result:
            return {"Error": GuideCountErrorType.MULTIPLE_MATCH_EXACT,
            "exact_match": True, 
            "num_matches": whitelist_guide_reporter_df_match.shape[0],
            "matches": whitelist_guide_reporter_df_match,
            "hamming_min": 0}
        else:
            return GuideCountErrorType.MULTIPLE_MATCH_EXACT
        #raise Exception("Multiple exact matches of the provided whitelisted guides - there are likely duplicates in the provided whitelist, please remove. Observed guide={}, guide matches={}".format(observed_guide_sequence, guide_sequences_series_match)) # NOTE 12/6/22: REMOVED THIS EXCEPTION - another reason is from truncated guides having multiple matches. In production code, just make sure to ensure that the whitelist is the set.


@typechecked
def infer_whitelist_sequence(observed_guide_reporter_sequence_input: Tuple[str, Optional[str], Optional[str]], 
        whitelist_guide_reporter_df: pd.DataFrame, 
        contains_surrogate:bool, 
        contains_barcode:bool, 
        contains_umi:bool, 
        encoded_whitelist_protospacer_sequences_series: np.ndarray, 
        encoded_whitelist_surrogate_sequences_series: Optional[np.ndarray] = None, 
        encoded_whitelist_barcode_sequences_series: Optional[np.ndarray] = None, 
        count_duplicate_mappings: bool = False, 
        protospacer_hamming_threshold: int = 7, 
        surrogate_hamming_threshold: int = 10, 
        barcode_hamming_threshold: int = 2):
    observed_reporter_sequences_indices = ["protospacer"]
    if contains_surrogate:
        observed_reporter_sequences_indices.append("surrogate")
    if contains_barcode:
        observed_reporter_sequences_indices.append("barcode")
    observed_guide_reporter_sequence = pd.Series(observed_guide_reporter_sequence_input, index=observed_reporter_sequences_indices)
    

    def validate_observed_sequence(sequence_name: str, encoded_whitelist_sequence_series: np.ndarray, missing_info_error: MissingInfoGuideCountError, unequal_length_error: UnequalLengthGuideCountError):
        if (observed_guide_reporter_sequence[sequence_name] is None) or (observed_guide_reporter_sequence[sequence_name] == "None") or (observed_guide_reporter_sequence[sequence_name].strip() == ""):
            return missing_info_error
        
        if len(observed_guide_reporter_sequence[sequence_name]) != encoded_whitelist_sequence_series.shape[1]: 
            return unequal_length_error   
        
        return None

    protospacer_error_result = validate_observed_sequence(sequence_name="protospacer", encoded_whitelist_sequence_series=encoded_whitelist_protospacer_sequences_series, missing_info_error=ProtospacerMissingInfoGuideCountError(sequence_value=observed_guide_reporter_sequence), unequal_length_error=ProtospacerUnequalLengthGuideCountError(sequence_length=len(observed_guide_reporter_sequence["protospacer"]), expected_length=encoded_whitelist_protospacer_sequences_series.shape[1]))
    surrogate_error_result = None
    if contains_surrogate:
        surrogate_error_result = validate_observed_sequence(sequence_name="surrogate", encoded_whitelist_sequence_series=encoded_whitelist_surrogate_sequences_series, missing_info_error=SurrogateMissingInfoGuideCountError(sequence_value=observed_guide_reporter_sequence), unequal_length_error=SurrogateUnequalLengthGuideCountError(sequence_length=len(observed_guide_reporter_sequence["surrogate"]), expected_length=encoded_whitelist_surrogate_sequences_series.shape[1]))
    barcode_error_result = None
    if contains_barcode:
        barcode_error_result = validate_observed_sequence(sequence_name="barcode", encoded_whitelist_sequence_series=encoded_whitelist_barcode_sequences_series, missing_info_error=BarcodeMissingInfoGuideCountError(sequence_value=observed_guide_reporter_sequence), unequal_length_error=BarcodeUnequalLengthGuideCountError(sequence_length=len(observed_guide_reporter_sequence["barcode"]), expected_length=encoded_whitelist_barcode_sequences_series.shape[1]))

    
    # PERFORM MAPPPING ON BOTH THE PROTOSPACER-ONLY, FULL, AND SURROGATE MISMAPPED-PROTOSPACER.
    complete_match_result = CompleteInferenceMatchResult()
    
    # GET THE BARCODE-ONLY MATCHES
    if contains_barcode:
        barcode_available = True
        if barcode_error_result is None:
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

                barcode_available = False
                # IMPORTANT CHECK TODO: I wonder if this will be overwritten during the flow. Perhaps set barcode_error_result to the error_result, and let the flow write the error in
                barcode_error_result = BarcodeHammingThresholdGuideCountError(hamming_min=observed_barcode_sequence_dists_min, hamming_threshold=barcode_hamming_threshold)
        else:
            barcode_available = False
            

    
    # PREPARE THE PROTOSPACER-ONLYMATCHES
    if protospacer_error_result is None:
        observed_protospacer_sequence_encoded = sequence_encoding.encode_DNA_base_vectorized(sequence_encoding.numpify_string_vectorized(observed_guide_reporter_sequence["protospacer"]))  # Encode the observed protospacer
        observed_protospacer_sequence_dists = sequence_encoding.retrieve_hamming_distance_whitelist(observed_protospacer_sequence_encoded, encoded_whitelist_protospacer_sequences_series) # Hamming distance among all whitelist protospacers
        observed_protospacer_sequence_dists_min = observed_protospacer_sequence_dists.min() # Get minimum hamming distance
        protospacer_hamming_threshold_met = (observed_protospacer_sequence_dists_min < protospacer_hamming_threshold)
        if protospacer_hamming_threshold_met: 
            # SET THE PROTOSPACER-ONLYMATCHES
            protospacer_matches_indices = np.where(observed_protospacer_sequence_dists == observed_protospacer_sequence_dists_min)[0]
            whitelist_guide_reporter_df_hamming_protospacer_match = whitelist_guide_reporter_df.iloc[protospacer_matches_indices] # Get all whitelisted guides with the minimum hamming distance (could be multiple)
            encoded_whitelist_protospacer_sequences_series_protospacer_match = encoded_whitelist_protospacer_sequences_series[protospacer_matches_indices] # Subset the protospacer encodings with the protospacer matches for later
            complete_match_result.protospacer_match = MatchSetSingleInferenceMatchResult(value=MatchSetSingleInferenceMatchResultValue(matches=whitelist_guide_reporter_df_hamming_protospacer_match))
            
            if contains_barcode: 
                if barcode_available:
                    
                    # PREPARE PROTOSPACER-MATCH, BARCODE-MATCH
                    observed_protospacer_sequence_dists_barcode_match = sequence_encoding.retrieve_hamming_distance_whitelist(observed_protospacer_sequence_encoded, encoded_whitelist_protospacer_sequences_series_barcode_match) # Hamming distance among barcode-match protospacers
                    observed_protospacer_sequence_dists_barcode_match_min = observed_protospacer_sequence_dists_barcode_match.min() # Get minimum hamming distance
                    barcode_match_protospacer_hamming_threshold_met = (observed_protospacer_sequence_dists_barcode_match_min < protospacer_hamming_threshold)
                    if barcode_match_protospacer_hamming_threshold_met: 
                        
                        # SET PROTOSPACER-MATCH, BARCODE-MATCH
                        barcode_match_protospacer_match_indices = np.where(observed_protospacer_sequence_dists_barcode_match == observed_protospacer_sequence_dists_barcode_match_min)[0]
                        whitelist_guide_reporter_df_hamming_barcode_match_protospacer_match = whitelist_guide_reporter_df_barcode_match.iloc[barcode_match_protospacer_match_indices] # Get all whitelisted guides with the minimum hamming distance (could be multiple)
                        complete_match_result.protospacer_match_barcode_match = MatchSetSingleInferenceMatchResult(value=MatchSetSingleInferenceMatchResultValue(matches=whitelist_guide_reporter_df_hamming_barcode_match_protospacer_match))
                        if contains_surrogate: 
                            if surrogate_error_result is None:
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
                                    complete_match_result.protospacer_match_surrogate_match_barcode_match = MatchSetSingleInferenceMatchResult(value=MatchSetSingleInferenceMatchResultValue(matches=whitelist_guide_reporter_df_hamming_barcode_match_protospacer_match_surrogate_match))
                                else:
                                    complete_match_result.protospacer_match_surrogate_match_barcode_match = MatchSetSingleInferenceMatchResult(
                                        error=SurrogateHammingThresholdGuideCountError(
                                            hamming_min=observed_surrogate_sequence_dists_barcode_match_protospacer_match_min, 
                                            hamming_threshold=surrogate_hamming_threshold))
                            else:
                                complete_match_result.protospacer_match_surrogate_match_barcode_match = MatchSetSingleInferenceMatchResult(error=surrogate_error_result)
                    else:
                        complete_match_result.protospacer_match_barcode_match = MatchSetSingleInferenceMatchResult(
                                        error=ProtospacerHammingThresholdGuideCountError(
                                            hamming_min=observed_protospacer_sequence_dists_barcode_match_min, 
                                            hamming_threshold=protospacer_hamming_threshold))
                else:
                    complete_match_result.protospacer_match_barcode_match = MatchSetSingleInferenceMatchResult(error=barcode_error_result)
                    complete_match_result.protospacer_match_surrogate_match_barcode_match = MatchSetSingleInferenceMatchResult(error=barcode_error_result)

            
            # PROTOSPACER-MATCH, SURROGATE-MATCH, NO BARCODE
            if contains_surrogate:
                if surrogate_error_result is None:
                    encoded_whitelist_surrogate_sequences_series_protospacer_match = encoded_whitelist_surrogate_sequences_series[protospacer_matches_indices]
                    
                    observed_surrogate_sequence_encoded = sequence_encoding.encode_DNA_base_vectorized(sequence_encoding.numpify_string_vectorized(observed_guide_reporter_sequence["surrogate"]))  # Encode the observed protospacer
                    observed_surrogate_sequence_dists_protospacer_match = sequence_encoding.retrieve_hamming_distance_whitelist(observed_surrogate_sequence_encoded, encoded_whitelist_surrogate_sequences_series_protospacer_match) # TODO: This should be on the protospacer-match array
                    observed_surrogate_sequence_dists_protospacer_match_min = observed_surrogate_sequence_dists_protospacer_match.min()
                    protospacer_match_surrogate_hamming_threshold_met = (observed_surrogate_sequence_dists_protospacer_match_min < surrogate_hamming_threshold)
                    if protospacer_match_surrogate_hamming_threshold_met:
                        protospacer_match_surrogate_match_indices = np.where(observed_surrogate_sequence_dists_protospacer_match == observed_surrogate_sequence_dists_protospacer_match_min)[0]
                        whitelist_guide_reporter_df_hamming_protospacer_match_surrogate_match = whitelist_guide_reporter_df_hamming_protospacer_match.iloc[protospacer_match_surrogate_match_indices] # Get all whitelisted guides with the minimum hamming distance (could be multiple)
                        complete_match_result.protospacer_match_surrogate_match = MatchSetSingleInferenceMatchResult(value=MatchSetSingleInferenceMatchResultValue(matches=whitelist_guide_reporter_df_hamming_protospacer_match_surrogate_match))
                    else:
                        complete_match_result.protospacer_match_surrogate_match = MatchSetSingleInferenceMatchResult(
                                        error=SurrogateHammingThresholdGuideCountError(
                                            hamming_min=observed_surrogate_sequence_dists_protospacer_match_min, 
                                            hamming_threshold=surrogate_hamming_threshold))
                else:
                    complete_match_result.protospacer_match_surrogate_match = MatchSetSingleInferenceMatchResult(error=surrogate_error_result)





            # GET THE PROTOSPACER-SURROGATE MISMATCHES NOW (if surrogate is available)
            if contains_surrogate:
                if surrogate_error_result is None:
                    if contains_barcode:
                        if barcode_available:
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
                                complete_match_result.protospacer_mismatch_surrogate_match_barcode_match = SurrogateProtospacerMismatchSingleInferenceMatchResult(value=SurrogateProtospacerMismatchSingleInferenceMatchResultValue(
                                    mismatched=whitelist_guide_reporter_df_hamming_barcode_match_surrogate_match_protospacer_match.empty, 
                                    surrogate_matches=whitelist_guide_reporter_df_hamming_barcode_match_surrogate_match_protospacer_match, 
                                    protospacer_matches=whitelist_guide_reporter_df_hamming_protospacer_match,
                                    protospacer_surrogate_matches=whitelist_guide_reporter_df_hamming_barcode_match_surrogate_match_protospacer_match
                                ))
                            else:
                                complete_match_result.protospacer_mismatch_surrogate_match_barcode_match = SurrogateProtospacerMismatchSingleInferenceMatchResult(
                                    error=SurrogateHammingThresholdGuideCountError(
                                        hamming_min=observed_surrogate_sequence_dists_barcode_match_min,
                                        hamming_threshold=surrogate_hamming_threshold
                                    ))
                        else:
                            complete_match_result.protospacer_mismatch_surrogate_match_barcode_match = SurrogateProtospacerMismatchSingleInferenceMatchResult(error=barcode_error_result)

                    observed_surrogate_sequence_encoded = sequence_encoding.encode_DNA_base_vectorized(sequence_encoding.numpify_string_vectorized(observed_guide_reporter_sequence["surrogate"]))  # Encode the observed protospacer
                    observed_surrogate_sequence_dists = sequence_encoding.retrieve_hamming_distance_whitelist(observed_surrogate_sequence_encoded, encoded_whitelist_surrogate_sequences_series) # Hamming distance among  all whitelist surrogates
                    observed_surrogate_sequence_dists_min = observed_surrogate_sequence_dists.min()
                    surrogate_hamming_threshold_met = (observed_surrogate_sequence_dists_min < surrogate_hamming_threshold)
                    if surrogate_hamming_threshold_met:
                        surrogate_match_indices = np.where(observed_surrogate_sequence_dists == observed_surrogate_sequence_dists_min)[0]
                        whitelist_guide_reporter_df_hamming_surrogate_match = whitelist_guide_reporter_df.iloc[surrogate_match_indices] # Get all whitelisted guides with the minimum hamming distance (could be multiple)

                        # See if there are identical matches between protospacer-only matches and surrogate-only matches
                        whitelist_guide_reporter_df_hamming_surrogate_match_protospacer_match = pd.merge(whitelist_guide_reporter_df_hamming_surrogate_match, whitelist_guide_reporter_df_hamming_protospacer_match, how='inner')
                        complete_match_result.protospacer_mismatch_surrogate_match = SurrogateProtospacerMismatchSingleInferenceMatchResult(value=SurrogateProtospacerMismatchSingleInferenceMatchResultValue(
                            mismatched=whitelist_guide_reporter_df_hamming_surrogate_match_protospacer_match.empty, 
                            surrogate_matches=whitelist_guide_reporter_df_hamming_surrogate_match, 
                            protospacer_matches=whitelist_guide_reporter_df_hamming_protospacer_match,
                            protospacer_surrogate_matches=whitelist_guide_reporter_df_hamming_surrogate_match_protospacer_match
                        ))
                    else:
                        complete_match_result.protospacer_mismatch_surrogate_match = SurrogateProtospacerMismatchSingleInferenceMatchResult(
                            error=SurrogateHammingThresholdGuideCountError(
                                hamming_min=observed_surrogate_sequence_dists_min,
                                hamming_threshold=surrogate_hamming_threshold
                            ))
                else:
                    complete_match_result.protospacer_mismatch_surrogate_match_barcode_match = SurrogateProtospacerMismatchSingleInferenceMatchResult(error=surrogate_error_result)
                    complete_match_result.protospacer_mismatch_surrogate_match = SurrogateProtospacerMismatchSingleInferenceMatchResult(error=surrogate_error_result)
                
        else: # If the barcode surpasses the threshold, fail the read due to no barcode. 
            error_result = ProtospacerHammingThresholdGuideCountError(
                hamming_min=observed_protospacer_sequence_dists_min,
                hamming_threshold=protospacer_hamming_threshold
            )
            complete_match_result.protospacer_match = MatchSetSingleInferenceMatchResult(error=error_result)
            complete_match_result.protospacer_match_surrogate_match_barcode_match = MatchSetSingleInferenceMatchResult(error=error_result)
            complete_match_result.protospacer_match_surrogate_match = MatchSetSingleInferenceMatchResult(error=error_result)
            complete_match_result.protospacer_match_barcode_match = MatchSetSingleInferenceMatchResult(error=error_result)
            complete_match_result.protospacer_mismatch_surrogate_match_barcode_match = SurrogateProtospacerMismatchSingleInferenceMatchResult(error=error_result)
            complete_match_result.protospacer_mismatch_surrogate_match = SurrogateProtospacerMismatchSingleInferenceMatchResult(error=error_result)
    else:
        complete_match_result.protospacer_match = MatchSetSingleInferenceMatchResult(error=protospacer_error_result)
        complete_match_result.protospacer_match_surrogate_match_barcode_match = MatchSetSingleInferenceMatchResult(error=protospacer_error_result)
        complete_match_result.protospacer_match_surrogate_match = MatchSetSingleInferenceMatchResult(error=protospacer_error_result)
        complete_match_result.protospacer_match_barcode_match = MatchSetSingleInferenceMatchResult(error=protospacer_error_result)
        complete_match_result.protospacer_mismatch_surrogate_match_barcode_match = SurrogateProtospacerMismatchSingleInferenceMatchResult(error=protospacer_error_result)
        complete_match_result.protospacer_mismatch_surrogate_match = SurrogateProtospacerMismatchSingleInferenceMatchResult(error=protospacer_error_result)
        
    
    return complete_match_result



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



# TODO: There will probably be some type errors with the DefaultDict when testing on non UMI (since it requires CounterType), so make sure to test with different variations of inputs
@typechecked
def get_whitelist_reporter_counts_with_umi(observed_guide_reporter_umi_counts: DefaultDict[Tuple[str,Optional[str],Optional[str]], Union[int, CounterType[Optional[str]]]], whitelist_guide_reporter_df: pd.DataFrame, contains_surrogate:bool = False, contains_barcode:bool = False, contains_umi:bool = False, protospacer_hamming_threshold_strict: Optional[int] = 7, surrogate_hamming_threshold_strict: Optional[int] = 2, barcode_hamming_threshold_strict: Optional[int] = 2, cores: int=1):
    encoded_whitelist_protospacer_sequences_series = sequence_encoding.encode_guide_series(whitelist_guide_reporter_df["protospacer"])
    if contains_surrogate:
        encoded_whitelist_surrogate_sequences_series = sequence_encoding.encode_guide_series(whitelist_guide_reporter_df["surrogate"])
    if contains_barcode:
        encoded_whitelist_barcode_sequences_series = sequence_encoding.encode_guide_series(whitelist_guide_reporter_df["barcode"])

    # TODO: Pass in arguments to set this hamming threshold. 
    # NOTE: Could there be an even more dynamic threshold, that certain mapped guides can have a higher threshold (i.e. those with many editable bases) and other guides can have low hamming (i.e. those with not many editable bases)
    protospacer_hamming_threshold_dynamic = False
    if protospacer_hamming_threshold_strict is None:
        protospacer_hamming_threshold_dynamic = True
        protospacer_hamming_threshold: int = determine_hamming_threshold(whitelist_guide_reporter_df["protospacer"], encoded_whitelist_protospacer_sequences_series, sample_count = 100, quantile = 0.05)
    else:
        protospacer_hamming_threshold: int = protospacer_hamming_threshold_strict
    print("Protospacer hamming threshold is " + str(protospacer_hamming_threshold))

    if contains_surrogate:
        surrogate_hamming_threshold_dynamic = False
        if surrogate_hamming_threshold_strict is  None:
            surrogate_hamming_threshold_dynamic = True
            surrogate_hamming_threshold: int = determine_hamming_threshold(whitelist_guide_reporter_df["surrogate"], encoded_whitelist_surrogate_sequences_series, sample_count = 100, quantile = 0.05)
            
        else:
            surrogate_hamming_threshold: int = surrogate_hamming_threshold_strict
        print("Surrogate hamming threshold is " + str(surrogate_hamming_threshold))

    if contains_barcode:
        barcode_hamming_threshold_dynamic = False
        if barcode_hamming_threshold_strict is  None:
            barcode_hamming_threshold_dynamic = True
            barcode_hamming_threshold: int = determine_hamming_threshold(whitelist_guide_reporter_df["barcode"], encoded_whitelist_barcode_sequences_series, sample_count = 100, quantile = 0.05)
        else:
            barcode_hamming_threshold: int = barcode_hamming_threshold_strict
        print("Barcode hamming threshold is " + str(barcode_hamming_threshold))



    # Infer whitelist guides from observed guides
    print("Inferring the true guides from observed guides")

    infer_whitelist_sequence_p = partial(infer_whitelist_sequence,
            whitelist_guide_reporter_df=whitelist_guide_reporter_df,
            contains_surrogate=contains_surrogate,
            contains_barcode=contains_barcode,
            contains_umi=contains_umi,
            encoded_whitelist_protospacer_sequences_series=encoded_whitelist_protospacer_sequences_series,
            encoded_whitelist_surrogate_sequences_series=encoded_whitelist_surrogate_sequences_series,
            encoded_whitelist_barcode_sequences_series=encoded_whitelist_barcode_sequences_series,
            count_duplicate_mappings=False,
            protospacer_hamming_threshold=protospacer_hamming_threshold, 
            surrogate_hamming_threshold=surrogate_hamming_threshold_strict, 
            barcode_hamming_threshold=barcode_hamming_threshold_strict)

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
        observed_guide_reporter_umi_counts_inferred[observed_guide_reporter_key] = {
            "observed_value": observed_guide_reporter_umi_counts[observed_guide_reporter_key],
            "inferred_value": inferred_true_reporter_sequences[observed_guide_reporter_key_index]
        }
    
    print("Completed inference")
    
    '''
        QC
    '''
    print("Retrieving QC tables")
    check_match_result_non_error = lambda match_result : False if match_result is None else match_result.error is None # If match_result is None, treat as error. If match_result is not None, but error is None, then non_error
    get_non_error_dict = lambda attribute_name : {observed_guide_reporter_key: observed_guide_reporter_umi_counts_inferred_value for observed_guide_reporter_key, observed_guide_reporter_umi_counts_inferred_value in observed_guide_reporter_umi_counts_inferred.items() if check_match_result_non_error(getattr(observed_guide_reporter_umi_counts_inferred_value["inferred_value"], attribute_name))}
    get_umi_noncollapsed_counts = lambda counts_inferred_dict : sum([sum(counts_inferred_value["observed_value"].values()) for counts_inferred_value in counts_inferred_dict.values()]) # NOTE: for UMI, observed_value is a Counter, so need to use .values()
    get_umi_collapsed_counts = lambda counts_inferred_dict : sum([len(counts_inferred_value["observed_value"].values()) for counts_inferred_value in counts_inferred_dict.values()]) # NOTE: for UMI, observed_value is a Counter, so need to use .values()
    get_counts = lambda counts_inferred_dict : sum([counts_inferred_value["observed_value"] for counts_inferred_value in counts_inferred_dict.values()]) # NOTE: for UMI, observed_value is an int, so NO NEED to use .values()
    
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
            error_result: Optional[GuideCountError] = getattr(observed_guide_reporter_umi_counts_inferred_value["inferred_value"], attribute_name).error if getattr(observed_guide_reporter_umi_counts_inferred_value["inferred_value"], attribute_name) is not None else GuideCountError(guide_count_error_type=GuideCountErrorType.NULL_MATCH_RESULT)
            if error_result is not None:
                if contains_umi:
                    guide_count_error_type_umi_noncollapsed_count[error_result.guide_count_error_type] += sum(observed_guide_reporter_umi_counts_inferred_value["observed_value"].values())
                    guide_count_error_type_umi_collapsed_count[error_result.guide_count_error_type] += len(observed_guide_reporter_umi_counts_inferred_value["observed_value"].values())
                else:
                    guide_count_error_type_count[error_result.guide_count_error_type] += observed_guide_reporter_umi_counts_inferred_value["observed_value"]

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
    
    '''
        Get finally count series (both UMI and non UMI version)
    '''
    @typechecked
    def get_matchset_counterseries(attribute_name: str) -> WhitelistReporterCounterSeriesResults:
        ambiguous_ignored_umi_noncollapsed_counterdict : DefaultDict[Tuple[str, Optional[str], Optional[str]], int]  = defaultdict(int)
        ambiguous_ignored_umi_collapsed_counterdict : DefaultDict[Tuple[str, Optional[str], Optional[str]], int]  = defaultdict(int)
        ambiguous_ignored_counterdict : DefaultDict[Tuple[str, Optional[str], Optional[str]], int]  = defaultdict(int)

        ambiguous_accepted_umi_noncollapsed_counterdict : DefaultDict[Tuple[str, Optional[str], Optional[str]], int]  = defaultdict(int)
        ambiguous_accepted_umi_collapsed_counterdict : DefaultDict[Tuple[str, Optional[str], Optional[str]], int]  = defaultdict(int)
        ambiguous_accepted_counterdict : DefaultDict[Tuple[str, Optional[str], Optional[str]], int]  = defaultdict(int)

        ambiguous_spread_umi_noncollapsed_counterdict : DefaultDict[Tuple[str, Optional[str], Optional[str]], float]  = defaultdict(float)
        ambiguous_spread_umi_collapsed_counterdict : DefaultDict[Tuple[str, Optional[str], Optional[str]], float]  = defaultdict(float)
        ambiguous_spread_counterdict : DefaultDict[Tuple[str, Optional[str], Optional[str]], float]  = defaultdict(float)

        for inferred_value_results in get_non_error_dict(attribute_name).values():
            observed_value_counts: Union[CounterType[Optional[str]], int] = inferred_value_results["observed_value"]
            inferred_value_result: CompleteInferenceMatchResult =  inferred_value_results["inferred_value"]
            match_set_single_inference_match_result : Optional[MatchSetSingleInferenceMatchResult] = getattr(inferred_value_result, attribute_name)
            assert match_set_single_inference_match_result is not None, "match_set_single_inference_match_result should not be none since this is from the non error list. Developer error."
            
            matches: pd.DataFrame = match_set_single_inference_match_result.value.matches
            if not matches.empty():

                for whitelist_reporter_series in matches.iterrows(): # LEFTOFF
                    if contains_umi:
                        assert isinstance(observed_value_counts, Counter), f"For UMI, expecting observed value is a Counter, but type is {type(observed_value_counts)}"
                        
                        ambiguous_accepted_umi_noncollapsed_counterdict[tuple(whitelist_reporter_series)] += sum(observed_value_counts.values())
                        ambiguous_accepted_umi_collapsed_counterdict[tuple(whitelist_reporter_series)] += len(observed_value_counts.values())

                        ambiguous_spread_umi_noncollapsed_counterdict[tuple(whitelist_reporter_series)] += sum(observed_value_counts.values()) / float(matches.shape[0])
                        ambiguous_spread_umi_collapsed_counterdict[tuple(whitelist_reporter_series)] += len(observed_value_counts.values()) / float(matches.shape[0])
                        
                        # If there is no ambiguous matches, then add to ambiguous_ignored counter
                        if matches.shape[0] == 1:
                            ambiguous_ignored_umi_noncollapsed_counterdict[tuple(whitelist_reporter_series)] += sum(observed_value_counts.values())
                            ambiguous_ignored_umi_collapsed_counterdict[tuple(whitelist_reporter_series)] += len(observed_value_counts.values())
                    else:
                        assert isinstance(observed_value_counts, int), f"For non UMI, expecting observed value is an int, but type is {type(observed_value_counts)}"
                        ambiguous_accepted_counterdict[tuple(whitelist_reporter_series)] += observed_value_counts
                        ambiguous_spread_counterdict[tuple(whitelist_reporter_series)] += observed_value_counts / float(matches.shape[0])
                        
                        # If there is no ambiguous matches, then add to ambiguous_ignored counter
                        if matches.shape[0] == 1:
                            ambiguous_ignored_counterdict[tuple(whitelist_reporter_series)] += observed_value_counts
        
        def create_counterseries(counterdict: DefaultDict[Tuple[str, Optional[str], Optional[str]], Union[int, float]]) -> pd.Series:
            counterseries: pd.Series = whitelist_guide_reporter_df.apply(lambda reporter: counterdict[tuple(reporter)])
            counterseries.index = pd.MultiIndex.from_frame(whitelist_guide_reporter_df)
            return counterseries
        
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

    @typechecked
    def get_mismatchset_counterseries(attribute_name: str) -> TODO:
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
        ambiguous_ignored_umi_noncollapsed_match_counterdict : DefaultDict[Tuple[Tuple[str, Optional[str], Optional[str]], Tuple[str, Optional[str], Optional[str]]], int]  = defaultdict(int)
        ambiguous_ignored_umi_collapsed_match_counterdict : DefaultDict[Tuple[Tuple[str, Optional[str], Optional[str]], Tuple[str, Optional[str], Optional[str]]], int]  = defaultdict(int)
        ambiguous_ignored_match_counterdict : DefaultDict[Tuple[Tuple[str, Optional[str], Optional[str]], Tuple[str, Optional[str], Optional[str]]], int]  = defaultdict(int)

        ambiguous_accepted_umi_noncollapsed_match_counterdict : DefaultDict[Tuple[Tuple[str, Optional[str], Optional[str]],Tuple[str, Optional[str], Optional[str]]], int]  = defaultdict(int)
        ambiguous_accepted_umi_collapsed_match_counterdict : DefaultDict[Tuple[Tuple[str, Optional[str], Optional[str]],Tuple[str, Optional[str], Optional[str]]], int]  = defaultdict(int)
        ambiguous_accepted_match_counterdict : DefaultDict[Tuple[Tuple[str, Optional[str], Optional[str]],Tuple[str, Optional[str], Optional[str]]], int]  = defaultdict(int)

        ambiguous_spread_umi_noncollapsed_match_counterdict : DefaultDict[Tuple[Tuple[str, Optional[str], Optional[str]],Tuple[str, Optional[str], Optional[str]]], float]  = defaultdict(float)
        ambiguous_spread_umi_collapsed_match_counterdict : DefaultDict[Tuple[Tuple[str, Optional[str], Optional[str]],Tuple[str, Optional[str], Optional[str]]], float]  = defaultdict(float)
        ambiguous_spread_match_counterdict : DefaultDict[Tuple[Tuple[str, Optional[str], Optional[str]],Tuple[str, Optional[str], Optional[str]]], float]  = defaultdict(float)

        for inferred_value_results in get_non_error_dict(attribute_name).values():
            observed_value_counts: Union[CounterType[Optional[str]], int] = inferred_value_results["observed_value"]
            inferred_value_result: CompleteInferenceMatchResult =  inferred_value_results["inferred_value"]
            surrogate_protospacer_mismatch_single_inference_match_result : Optional[SurrogateProtospacerMismatchSingleInferenceMatchResult] = getattr(inferred_value_result, attribute_name)
            assert surrogate_protospacer_mismatch_single_inference_match_result is not None, "surrogate_protospacer_mismatch_single_inference_match_result should not be none since this is from the non error list. Developer error."
            
            mismatched: bool = surrogate_protospacer_mismatch_single_inference_match_result.value.mismatched
            surrogate_matches: pd.DataFrame = surrogate_protospacer_mismatch_single_inference_match_result.value.mismatched
            protospacer_matches: pd.DataFrame = surrogate_protospacer_mismatch_single_inference_match_result.value.mismatched
            protospacer_surrogate_matches: pd.DataFrame = surrogate_protospacer_mismatch_single_inference_match_result.value.mismatched


            if mismatched: # If mismatched, Add count to protospacer/surrogate pair
                # TODO: Make sure to handle when there is a match in one sequence but not the other - will this be called as a mismatch or throw an error?
                # LEFTOFFF HERE: Count the pairs - but first check in the inference code on the above TODO.
                pass
            else: # If matched, add count normally.
                assert not protospacer_surrogate_matches.empty(), f"mismatched==false, but the match dataframe is empty, unexpected paradox. Developer error"
                
                #TODO: This is near duplicate code to get_matchset_counterseries - modularize! 
                matches = protospacer_surrogate_matches
                for whitelist_reporter_series in matches.iterrows(): # LEFTOFF
                    if contains_umi:
                        assert isinstance(observed_value_counts, Counter), f"For UMI, expecting observed value is a Counter, but type is {type(observed_value_counts)}"
                        
                        ambiguous_accepted_umi_noncollapsed_match_counterdict[tuple(whitelist_reporter_series)] += sum(observed_value_counts.values())
                        ambiguous_accepted_umi_collapsed_match_counterdict[tuple(whitelist_reporter_series)] += len(observed_value_counts.values())

                        ambiguous_spread_umi_noncollapsed_match_counterdict[tuple(whitelist_reporter_series)] += sum(observed_value_counts.values()) / float(matches.shape[0])
                        ambiguous_spread_umi_collapsed_match_counterdict[tuple(whitelist_reporter_series)] += len(observed_value_counts.values()) / float(matches.shape[0])
                        
                        # If there is no ambiguous matches, then add to ambiguous_ignored counter
                        if matches.shape[0] == 1:
                            ambiguous_ignored_umi_noncollapsed_match_counterdict[tuple(whitelist_reporter_series)] += sum(observed_value_counts.values())
                            ambiguous_ignored_umi_collapsed_match_counterdict[tuple(whitelist_reporter_series)] += len(observed_value_counts.values())
                    else:
                        assert isinstance(observed_value_counts, int), f"For non UMI, expecting observed value is an int, but type is {type(observed_value_counts)}"
                        ambiguous_accepted_match_counterdict[tuple(whitelist_reporter_series)] += observed_value_counts
                        ambiguous_spread_match_counterdict[tuple(whitelist_reporter_series)] += observed_value_counts / float(matches.shape[0])
                        
                        # If there is no ambiguous matches, then add to ambiguous_ignored counter
                        if matches.shape[0] == 1:
                            ambiguous_ignored_match_counterdict[tuple(whitelist_reporter_series)] += observed_value_counts
        


        def create_counterseries(counterdict: DefaultDict[Tuple[str, Optional[str], Optional[str]], Union[int, float]]) -> pd.Series:
            counterseries: pd.Series = whitelist_guide_reporter_df.apply(lambda reporter: counterdict[tuple(reporter)])
            counterseries.index = pd.MultiIndex.from_frame(whitelist_guide_reporter_df)
            return counterseries
        
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

    all_match_set_whitelist_reporter_counter_series_results = AllMatchSetWhitelistReporterCounterSeriesResults()
    all_match_set_whitelist_reporter_counter_series_results.protospacer_match = get_matchset_counterseries("protospacer_match")
    if contains_barcode:
        all_match_set_whitelist_reporter_counter_series_results.protospacer_match_barcode_match = get_matchset_counterseries("protospacer_match_barcode_match")
        if contains_surrogate:
            all_match_set_whitelist_reporter_counter_series_results.protospacer_match_surrogate_match_barcode_match = get_matchset_counterseries("protospacer_match_surrogate_match_barcode_match")
            
            #protospacer_mismatch_surrogate_match_barcode_match_nonerror_dict = get_non_error_dict("protospacer_mismatch_surrogate_match_barcode_match")
            #quality_control_result.protospacer_mismatch_surrogate_match_barcode_match = set_match_set_single_inference_quality_control_results("protospacer_mismatch_surrogate_match_barcode_match")
    if contains_surrogate:
        pass
        all_match_set_whitelist_reporter_counter_series_results.protospacer_match_surrogate_match = get_matchset_counterseries("protospacer_match_surrogate_match")
        #quality_control_result.protospacer_mismatch_surrogate_match = set_match_set_single_inference_quality_control_results("protospacer_mismatch_surrogate_match")


    return WhitelistReporterCountsResult(observed_guide_reporter_umi_counts_inferred=observed_guide_reporter_umi_counts_inferred, quality_control_result=quality_control_result)