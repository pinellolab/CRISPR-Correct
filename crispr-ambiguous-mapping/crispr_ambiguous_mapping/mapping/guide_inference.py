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
        protospacer_hamming_threshold: int = 7, 
        surrogate_hamming_threshold: int = 10, 
        barcode_hamming_threshold: int = 2):

    #
    #   FORMAT OBSERVED SEQUENCES AS PANDAS SERIES
    #
    observed_reporter_sequences_indices = ["protospacer"]
    if contains_surrogate:
        observed_reporter_sequences_indices.append("surrogate")
    if contains_barcode:
        observed_reporter_sequences_indices.append("barcode")
    observed_guide_reporter_sequence = pd.Series(observed_guide_reporter_sequence_input, index=observed_reporter_sequences_indices)
    

    #
    # QUICK VALIDATION OF THE OBSERVED SEQUENCES, SAVE RESULTS AS VARIABLES TO BE USED DOWNSTREAM
    #
    def validate_observed_sequence(sequence_name: str, encoded_whitelist_sequence_series: np.ndarray, missing_info_error: MissingInfoGuideCountError, insufficient_length_error: InsufficientLengthGuideCountError):
        if (observed_guide_reporter_sequence[sequence_name] is None) or (observed_guide_reporter_sequence[sequence_name] == "None") or (observed_guide_reporter_sequence[sequence_name].strip() == ""):
            return missing_info_error
        
        if len(observed_guide_reporter_sequence[sequence_name]) < insufficient_length_error.minimum_length: 
            return insufficient_length_error   
        
        return None
        
    protospacer_error_result = validate_observed_sequence(sequence_name="protospacer", encoded_whitelist_sequence_series=encoded_whitelist_protospacer_sequences_series, missing_info_error=ProtospacerMissingInfoGuideCountError(sequence_value=observed_guide_reporter_sequence), insufficient_length_error=ProtospacerInsufficientLengthGuideCountError(sequence_length=len(observed_guide_reporter_sequence["protospacer"]), minimum_length=whitelist_guide_reporter_df["protospacer"].apply(len).min()))
    surrogate_error_result = None
    if contains_surrogate:
        surrogate_error_result = validate_observed_sequence(sequence_name="surrogate", encoded_whitelist_sequence_series=encoded_whitelist_surrogate_sequences_series, missing_info_error=SurrogateMissingInfoGuideCountError(sequence_value=observed_guide_reporter_sequence), insufficient_length_error=SurrogateInsufficientLengthGuideCountError(sequence_length=len(observed_guide_reporter_sequence["surrogate"]), minimum_length=whitelist_guide_reporter_df["surrogate"].apply(len).min()))
    barcode_error_result = None
    if contains_barcode:
        barcode_error_result = validate_observed_sequence(sequence_name="barcode", encoded_whitelist_sequence_series=encoded_whitelist_barcode_sequences_series, missing_info_error=BarcodeMissingInfoGuideCountError(sequence_value=observed_guide_reporter_sequence), insufficient_length_error=BarcodeInsufficientLengthGuideCountError(sequence_length=len(observed_guide_reporter_sequence["barcode"]), minimum_length=whitelist_guide_reporter_df["barcode"].apply(len).min()))

    # CHANGE SHAPE OF ENCODING TO BE SAME LENGTH OF 
    encoded_whitelist_protospacer_sequences_series = encoded_whitelist_protospacer_sequences_series[:, :len(observed_guide_reporter_sequence["protospacer"]), :]
    encoded_whitelist_surrogate_sequences_series = encoded_whitelist_surrogate_sequences_series[:, :len(observed_guide_reporter_sequence["surrogate"]), :]
    encoded_whitelist_barcode_sequences_series = encoded_whitelist_barcode_sequences_series[:, :len(observed_guide_reporter_sequence["barcode"]), :]
    
    # PERFORM MAPPPING ON BOTH THE PROTOSPACER-ONLY, FULL, AND SURROGATE MISMAPPED-PROTOSPACER.
    complete_match_result = CompleteInferenceMatchResult()

    #
    # GET THE BARCODE-ONLY MATCHES
    # (Since barcode-matching will be done a lot, let's do first and save the match results)
    if contains_barcode: # IF BARCODE IS PARSED
        barcode_available = True # Let's save an indicator on if we successfully parsed barcode for convenience downstream
        if barcode_error_result is None: # IF NO BARCODE PARSING ISSUE

            #
            # FIND BARCODE MATCHES
            #
            observed_barcode_sequence_encoded = sequence_encoding.encode_DNA_base_vectorized(sequence_encoding.numpify_string_vectorized(observed_guide_reporter_sequence["barcode"]))  # Encode the observed barcode
            observed_barcode_sequence_dists = sequence_encoding.retrieve_hamming_distance_whitelist(observed_barcode_sequence_encoded, encoded_whitelist_barcode_sequences_series) # Retrieve hamming distance with whitelist barcode
            observed_barcode_sequence_dists_min = observed_barcode_sequence_dists.min() # Get the barcode with the minimum  hamming distance
            barcode_hamming_threshold_met = (observed_barcode_sequence_dists_min < barcode_hamming_threshold)

            barcode_matches_indices = np.where(observed_barcode_sequence_dists == observed_barcode_sequence_dists_min)[0] # Get the indices of ALL barcode matches
            whitelist_guide_reporter_df_barcode_match = whitelist_guide_reporter_df.iloc[barcode_matches_indices] # Get the whitelist reporters with the matched barcode(s)
            
            if barcode_hamming_threshold_met: # IF BARCODE MATCH
                encoded_whitelist_protospacer_sequences_series_barcode_match = encoded_whitelist_protospacer_sequences_series[barcode_matches_indices] # Subset the protospacer encodings with the barcode matches for later
                if contains_surrogate:
                    encoded_whitelist_surrogate_sequences_series_barcode_match = encoded_whitelist_surrogate_sequences_series[barcode_matches_indices] # Subset the surrogate encodings with the barcode matches for later
            else: # NO BARCODE MATCH, ERROR
                barcode_available = False
                barcode_error_result = BarcodeHammingThresholdGuideCountError(
                    hamming_min=observed_barcode_sequence_dists_min, 
                    hamming_threshold=barcode_hamming_threshold,
                    original_df=whitelist_guide_reporter_df,
                    hamming_min_match_df=whitelist_guide_reporter_df_barcode_match)
        else:
            barcode_available = False
            

    #
    # PREPARE THE PROTOSPACER-ONLYMATCHES
    #
    if protospacer_error_result is None:
        observed_protospacer_sequence_encoded = sequence_encoding.encode_DNA_base_vectorized(sequence_encoding.numpify_string_vectorized(observed_guide_reporter_sequence["protospacer"]))  # Encode the observed protospacer
        observed_protospacer_sequence_dists = sequence_encoding.retrieve_hamming_distance_whitelist(observed_protospacer_sequence_encoded, encoded_whitelist_protospacer_sequences_series) # Hamming distance among all whitelist protospacers
        observed_protospacer_sequence_dists_min = observed_protospacer_sequence_dists.min() # Get minimum hamming distance
        protospacer_hamming_threshold_met = (observed_protospacer_sequence_dists_min < protospacer_hamming_threshold)

        protospacer_matches_indices = np.where(observed_protospacer_sequence_dists == observed_protospacer_sequence_dists_min)[0]
        whitelist_guide_reporter_df_hamming_protospacer_match = whitelist_guide_reporter_df.iloc[protospacer_matches_indices] # Get all whitelisted guides with the minimum hamming distance (could be multiple)
        if protospacer_hamming_threshold_met: # PROTOSPACER MATCHED
            #
            # SET THE PROTOSPACER-ONLYMATCHES
            #
            encoded_whitelist_protospacer_sequences_series_protospacer_match = encoded_whitelist_protospacer_sequences_series[protospacer_matches_indices] # Subset the protospacer encodings with the protospacer matches for later
            complete_match_result.protospacer_match = MatchSetSingleInferenceMatchResult(value=MatchSetSingleInferenceMatchResultValue(matches=whitelist_guide_reporter_df_hamming_protospacer_match))
            
            if contains_barcode: # IF BARCODE IS PARSED, PROCEED TO BARCODE-MATCHING INFERENCE
                if barcode_available: # IF BARCODE IS MATCHED, PROCEED TO BARCODE-MATCHING INFERENCE
                    #
                    # PREPARE PROTOSPACER-MATCH, BARCODE-MATCH
                    #
                    observed_protospacer_sequence_dists_barcode_match = sequence_encoding.retrieve_hamming_distance_whitelist(observed_protospacer_sequence_encoded, encoded_whitelist_protospacer_sequences_series_barcode_match) # Hamming distance among barcode-match protospacers
                    observed_protospacer_sequence_dists_barcode_match_min = observed_protospacer_sequence_dists_barcode_match.min() # Get minimum hamming distance
                    barcode_match_protospacer_hamming_threshold_met = (observed_protospacer_sequence_dists_barcode_match_min < protospacer_hamming_threshold)

                    barcode_match_protospacer_match_indices = np.where(observed_protospacer_sequence_dists_barcode_match == observed_protospacer_sequence_dists_barcode_match_min)[0]
                    whitelist_guide_reporter_df_hamming_barcode_match_protospacer_match = whitelist_guide_reporter_df_barcode_match.iloc[barcode_match_protospacer_match_indices] # Get all whitelisted guides with the minimum hamming distance (could be multiple)
                    if barcode_match_protospacer_hamming_threshold_met:  # PROTOSPACER MATCHED WHEN BARCODE IS SPECIFIED
                        #
                        # SET PROTOSPACER-MATCH, BARCODE-MATCH
                        #
                        
                        complete_match_result.protospacer_match_barcode_match = MatchSetSingleInferenceMatchResult(value=MatchSetSingleInferenceMatchResultValue(matches=whitelist_guide_reporter_df_hamming_barcode_match_protospacer_match))
                        if contains_surrogate: # IF SURROGATE IS PARSED, PROCEED TO SURROGATE-MATCHING INFERENCE
                            if surrogate_error_result is None:
                                #
                                # PREPARE PROTOSPACER-MATCH, BARCODE-MATCH, SURROGATE-MATCH
                                #
                                encoded_whitelist_surrogate_sequences_series_barcode_match_protospacer_match = encoded_whitelist_surrogate_sequences_series_barcode_match[barcode_match_protospacer_match_indices] # Subset the surrogate encodings with the protospacer and encoding matches for later
                                observed_surrogate_sequence_encoded = sequence_encoding.encode_DNA_base_vectorized(sequence_encoding.numpify_string_vectorized(observed_guide_reporter_sequence["surrogate"]))  # Encode the observed protospacer
                                observed_surrogate_sequence_dists_barcode_match_protospacer_match = sequence_encoding.retrieve_hamming_distance_whitelist(observed_surrogate_sequence_encoded, encoded_whitelist_surrogate_sequences_series_barcode_match_protospacer_match) # Hamming distance among all whitelist sub-selected surrogates
                                observed_surrogate_sequence_dists_barcode_match_protospacer_match_min = observed_surrogate_sequence_dists_barcode_match_protospacer_match.min()
                                barcode_match_protospacer_match_surrogate_hamming_threshold_met = (observed_surrogate_sequence_dists_barcode_match_protospacer_match_min < surrogate_hamming_threshold)
                                
                                barcode_match_protospacer_match_surrogate_match_indices = np.where(observed_surrogate_sequence_dists_barcode_match_protospacer_match == observed_surrogate_sequence_dists_barcode_match_protospacer_match_min)[0]
                                whitelist_guide_reporter_df_hamming_barcode_match_protospacer_match_surrogate_match = whitelist_guide_reporter_df_hamming_barcode_match_protospacer_match.iloc[barcode_match_protospacer_match_surrogate_match_indices] # Get all whitelisted guides with the minimum hamming distance (could be multiple)
                                if barcode_match_protospacer_match_surrogate_hamming_threshold_met: # SURROGATE MATCHED
                                    #
                                    # SET PROTOSPACER-MATCH, BARCODE-MATCH, SURROGATE-MATCH
                                    #
                                    complete_match_result.protospacer_match_surrogate_match_barcode_match = MatchSetSingleInferenceMatchResult(value=MatchSetSingleInferenceMatchResultValue(matches=whitelist_guide_reporter_df_hamming_barcode_match_protospacer_match_surrogate_match))
                                else: # NO SURROGATE MATCH, ERROR
                                    complete_match_result.protospacer_match_surrogate_match_barcode_match = MatchSetSingleInferenceMatchResult(
                                        error=SurrogateHammingThresholdGuideCountError(
                                            hamming_min=observed_surrogate_sequence_dists_barcode_match_protospacer_match_min, 
                                            hamming_threshold=surrogate_hamming_threshold,
                                            original_df=whitelist_guide_reporter_df_hamming_barcode_match_protospacer_match,
                                            hamming_min_match_df=whitelist_guide_reporter_df_hamming_barcode_match_protospacer_match_surrogate_match,
                                            barcode_subsetted=True))
                            else: # PARSED SURROGATE HAS SOME UNEXPECTED ISSUE 
                                complete_match_result.protospacer_match_surrogate_match_barcode_match = MatchSetSingleInferenceMatchResult(error=surrogate_error_result)
                    else: # NO PROTOSPACER MATCH, ERROR
                        complete_match_result.protospacer_match_barcode_match = MatchSetSingleInferenceMatchResult(
                                        error=ProtospacerHammingThresholdGuideCountError(
                                            hamming_min=observed_protospacer_sequence_dists_barcode_match_min, 
                                            hamming_threshold=protospacer_hamming_threshold,
                                            original_df=whitelist_guide_reporter_df_barcode_match,
                                            hamming_min_match_df=whitelist_guide_reporter_df_hamming_barcode_match_protospacer_match,
                                            barcode_subsetted=True))
                else: # NO BARCODE MATCH, SET ERROR TO THE RESULTS REQUIREING A BARCODE MATCH
                    complete_match_result.protospacer_match_barcode_match = MatchSetSingleInferenceMatchResult(error=barcode_error_result)
                    complete_match_result.protospacer_match_surrogate_match_barcode_match = MatchSetSingleInferenceMatchResult(error=barcode_error_result)

            #
            # PROTOSPACER-MATCH, SURROGATE-MATCH, NO BARCODE
            #
            if contains_surrogate: # IF SURROGATE IS PARSED
                if surrogate_error_result is None: # AND NO SURROGATE PARSING ISSUE 

                    #
                    # GET SURROGATE-MATCHES ON THE PROTOSPACER-MATCHES
                    #
                    encoded_whitelist_surrogate_sequences_series_protospacer_match = encoded_whitelist_surrogate_sequences_series[protospacer_matches_indices]
                    observed_surrogate_sequence_encoded = sequence_encoding.encode_DNA_base_vectorized(sequence_encoding.numpify_string_vectorized(observed_guide_reporter_sequence["surrogate"]))  # Encode the observed protospacer
                    observed_surrogate_sequence_dists_protospacer_match = sequence_encoding.retrieve_hamming_distance_whitelist(observed_surrogate_sequence_encoded, encoded_whitelist_surrogate_sequences_series_protospacer_match) # TODO: This should be on the protospacer-match array
                    observed_surrogate_sequence_dists_protospacer_match_min = observed_surrogate_sequence_dists_protospacer_match.min()
                    protospacer_match_surrogate_hamming_threshold_met = (observed_surrogate_sequence_dists_protospacer_match_min < surrogate_hamming_threshold)

                    protospacer_match_surrogate_match_indices = np.where(observed_surrogate_sequence_dists_protospacer_match == observed_surrogate_sequence_dists_protospacer_match_min)[0]
                    whitelist_guide_reporter_df_hamming_protospacer_match_surrogate_match = whitelist_guide_reporter_df_hamming_protospacer_match.iloc[protospacer_match_surrogate_match_indices] # Get all whitelisted guides with the minimum hamming distance (could be multiple)
                    if protospacer_match_surrogate_hamming_threshold_met: # IF SURROGATE MATCH
                        complete_match_result.protospacer_match_surrogate_match = MatchSetSingleInferenceMatchResult(value=MatchSetSingleInferenceMatchResultValue(matches=whitelist_guide_reporter_df_hamming_protospacer_match_surrogate_match))
                    else: # NO SURROGATE MATCH, ERROR OUT
                        complete_match_result.protospacer_match_surrogate_match = MatchSetSingleInferenceMatchResult(
                                        error=SurrogateHammingThresholdGuideCountError(
                                            hamming_min=observed_surrogate_sequence_dists_protospacer_match_min, 
                                            hamming_threshold=surrogate_hamming_threshold,
                                            original_df=whitelist_guide_reporter_df_hamming_protospacer_match,
                                            hamming_min_match_df=whitelist_guide_reporter_df_hamming_protospacer_match_surrogate_match,
                                            barcode_subsetted=False))
                else: # SURROGATE PARSING ISSUE, ERROR
                    complete_match_result.protospacer_match_surrogate_match = MatchSetSingleInferenceMatchResult(error=surrogate_error_result)



            #
            # GET THE PROTOSPACER-SURROGATE MISMATCHES NOW (if surrogate is available)
            #
            if contains_surrogate: # IF SURROGATE IS PARSED
                if surrogate_error_result is None: # AND NO SURROGATE PARSING ISSUE
                    if contains_barcode: # IF BARCODE IS PARSED
                        if barcode_available: # AND NO BARCODE PARSING ISSUE
                            #
                            # PREPARE SURROGATE-MATCH, BARCODE-MATCH
                            #
                            observed_surrogate_sequence_encoded = sequence_encoding.encode_DNA_base_vectorized(sequence_encoding.numpify_string_vectorized(observed_guide_reporter_sequence["surrogate"]))  # Encode the observed protospacer
                            observed_surrogate_sequence_dists_barcode_match = sequence_encoding.retrieve_hamming_distance_whitelist(observed_surrogate_sequence_encoded, encoded_whitelist_surrogate_sequences_series_barcode_match) # Hamming distance among barcode-match protospacers
                            observed_surrogate_sequence_dists_barcode_match_min = observed_surrogate_sequence_dists_barcode_match.min() # Get minimum hamming distance
                            barcode_match_surrogate_hamming_threshold_met = (observed_surrogate_sequence_dists_barcode_match_min < surrogate_hamming_threshold)
                            
                             #
                            # SET SURROGATE-MATCH, BARCODE-MATCH
                            #
                            barcode_match_surrogate_match_indices = np.where(observed_surrogate_sequence_dists_barcode_match == observed_surrogate_sequence_dists_barcode_match_min)[0]
                            whitelist_guide_reporter_df_hamming_barcode_match_surrogate_match = whitelist_guide_reporter_df_barcode_match.iloc[barcode_match_surrogate_match_indices] # Get all whitelisted guides with the minimum hamming distance (could be multiple)

                            if barcode_match_surrogate_hamming_threshold_met: # If the surrogate matches are below the surrogate hamming threshold
                                # See if there are identical matches between protospacer-only matches and surrogate-only matches
                                whitelist_guide_reporter_df_hamming_barcode_match_surrogate_match_protospacer_match = pd.merge(whitelist_guide_reporter_df_hamming_barcode_match_surrogate_match, whitelist_guide_reporter_df_hamming_protospacer_match, how='inner')
                                
                                complete_match_result.protospacer_mismatch_surrogate_match_barcode_match = SurrogateProtospacerMismatchSingleInferenceMatchResult(value=SurrogateProtospacerMismatchSingleInferenceMatchResultValue(
                                    mismatched=whitelist_guide_reporter_df_hamming_barcode_match_surrogate_match_protospacer_match.empty, 
                                    surrogate_matches=whitelist_guide_reporter_df_hamming_barcode_match_surrogate_match, 
                                    protospacer_matches=whitelist_guide_reporter_df_hamming_protospacer_match,
                                    protospacer_surrogate_matches=whitelist_guide_reporter_df_hamming_barcode_match_surrogate_match_protospacer_match
                                ))
                            else: # NO SURROGATE MATCH
                                complete_match_result.protospacer_mismatch_surrogate_match_barcode_match = SurrogateProtospacerMismatchSingleInferenceMatchResult(
                                    error=SurrogateHammingThresholdGuideCountError(
                                        hamming_min=observed_surrogate_sequence_dists_barcode_match_min,
                                        hamming_threshold=surrogate_hamming_threshold,
                                        original_df=whitelist_guide_reporter_df_barcode_match,
                                        hamming_min_match_df=whitelist_guide_reporter_df_hamming_barcode_match_surrogate_match,
                                        barcode_subsetted=True,
                                    ))
                        else: # BARCODE PARSING ISSUE OR NO MATCH, ERROR OUT
                            complete_match_result.protospacer_mismatch_surrogate_match_barcode_match = SurrogateProtospacerMismatchSingleInferenceMatchResult(error=barcode_error_result)
                    #
                    # PREPARE SURROGATE MATCHES ON THE PROTOSPACER-MATCHES SEQUENCES
                    #
                    observed_surrogate_sequence_encoded = sequence_encoding.encode_DNA_base_vectorized(sequence_encoding.numpify_string_vectorized(observed_guide_reporter_sequence["surrogate"]))  # Encode the observed protospacer
                    observed_surrogate_sequence_dists = sequence_encoding.retrieve_hamming_distance_whitelist(observed_surrogate_sequence_encoded, encoded_whitelist_surrogate_sequences_series) # Hamming distance among  all whitelist surrogates
                    observed_surrogate_sequence_dists_min = observed_surrogate_sequence_dists.min()
                    surrogate_hamming_threshold_met = (observed_surrogate_sequence_dists_min < surrogate_hamming_threshold)

                    surrogate_match_indices = np.where(observed_surrogate_sequence_dists == observed_surrogate_sequence_dists_min)[0]
                    whitelist_guide_reporter_df_hamming_surrogate_match = whitelist_guide_reporter_df.iloc[surrogate_match_indices] # Get all whitelisted guides with the minimum hamming distance (could be multiple)
                    if surrogate_hamming_threshold_met: # IF SURROGATE MATCHED

                        # See if there are identical matches between protospacer-only matches and surrogate-only matches
                        whitelist_guide_reporter_df_hamming_surrogate_match_protospacer_match = pd.merge(whitelist_guide_reporter_df_hamming_surrogate_match, whitelist_guide_reporter_df_hamming_protospacer_match, how='inner')
                        complete_match_result.protospacer_mismatch_surrogate_match = SurrogateProtospacerMismatchSingleInferenceMatchResult(value=SurrogateProtospacerMismatchSingleInferenceMatchResultValue(
                            mismatched=whitelist_guide_reporter_df_hamming_surrogate_match_protospacer_match.empty, 
                            surrogate_matches=whitelist_guide_reporter_df_hamming_surrogate_match, 
                            protospacer_matches=whitelist_guide_reporter_df_hamming_protospacer_match,
                            protospacer_surrogate_matches=whitelist_guide_reporter_df_hamming_surrogate_match_protospacer_match
                        ))
                    else: # NO SURROGATE MATCH, ERROR
                        complete_match_result.protospacer_mismatch_surrogate_match = SurrogateProtospacerMismatchSingleInferenceMatchResult(
                            error=SurrogateHammingThresholdGuideCountError(
                                hamming_min=observed_surrogate_sequence_dists_min,
                                hamming_threshold=surrogate_hamming_threshold,
                                original_df=whitelist_guide_reporter_df,
                                hamming_min_match_df=whitelist_guide_reporter_df_hamming_surrogate_match,
                                barcode_subsetted=False
                            ))
                else: # SURROGATE PARSING ISSUE, ERROR
                    complete_match_result.protospacer_mismatch_surrogate_match_barcode_match = SurrogateProtospacerMismatchSingleInferenceMatchResult(error=surrogate_error_result)
                    complete_match_result.protospacer_mismatch_surrogate_match = SurrogateProtospacerMismatchSingleInferenceMatchResult(error=surrogate_error_result)
                
        else: # NO PROTOSPACER MATCH, ERROR
            error_result = ProtospacerHammingThresholdGuideCountError(
                hamming_min=observed_protospacer_sequence_dists_min,
                hamming_threshold=protospacer_hamming_threshold,
                original_df=whitelist_guide_reporter_df,
                hamming_min_match_df=whitelist_guide_reporter_df_hamming_protospacer_match,
                barcode_subsetted=False
            )
            complete_match_result.protospacer_match = MatchSetSingleInferenceMatchResult(error=error_result)
            complete_match_result.protospacer_match_surrogate_match_barcode_match = MatchSetSingleInferenceMatchResult(error=error_result)
            complete_match_result.protospacer_match_surrogate_match = MatchSetSingleInferenceMatchResult(error=error_result)
            complete_match_result.protospacer_match_barcode_match = MatchSetSingleInferenceMatchResult(error=error_result)
            complete_match_result.protospacer_mismatch_surrogate_match_barcode_match = SurrogateProtospacerMismatchSingleInferenceMatchResult(error=error_result)
            complete_match_result.protospacer_mismatch_surrogate_match = SurrogateProtospacerMismatchSingleInferenceMatchResult(error=error_result)
    else: # PARSED PROTOSPACER HAS SOME UNEXPECTED ISSUE, ERROR
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

