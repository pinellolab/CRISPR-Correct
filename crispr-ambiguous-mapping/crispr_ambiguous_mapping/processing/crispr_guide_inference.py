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
from ..models.mapping_models import CompleteInferenceMatchResult, MatchSetSingleInferenceMatchResult, MatchSetSingleInferenceMatchResultValue, SurrogateProtospacerMismatchSingleInferenceMatchResult, SurrogateProtospacerMismatchSingleInferenceMatchResultValue
from ..models.error_models import (GuideCountError, GuideCountErrorType,
                                   ProtospacerHammingThresholdGuideCountError, SurrogateHammingThresholdGuideCountError, BarcodeHammingThresholdGuideCountError,
                                   MissingInfoGuideCountError, InsufficientLengthGuideCountError,
                                   ProtospacerMissingInfoGuideCountError, SurrogateMissingInfoGuideCountError, BarcodeMissingInfoGuideCountError,
                                   ProtospacerInsufficientLengthGuideCountError, SurrogateInsufficientLengthGuideCountError, BarcodeInsufficientLengthGuideCountError)
import logging
_log = logging.getLogger(__name__)


@typechecked
def infer_whitelist_sequence(observed_guide_reporter_sequence_input: Union[str, Tuple[str], Tuple[str, str], Tuple[str,str, str]], 
        whitelist_guide_reporter_df: pd.DataFrame, 
        contains_guide_surrogate:bool, 
        contains_guide_barcode:bool, 
        contains_guide_umi:bool, 
        encoded_whitelist_protospacer_sequences_series: np.ndarray, 
        encoded_whitelist_surrogate_sequences_series: Optional[np.ndarray] = None, 
        encoded_whitelist_barcode_sequences_series: Optional[np.ndarray] = None, 
        protospacer_hamming_threshold: int = 7,
        surrogate_hamming_threshold: Optional[int] = 10,
        barcode_hamming_threshold: Optional[int] = 2,
        # PERF §3.6: accept precomputed library-column minimum lengths from the
        # caller so we don't redo `whitelist_df[col].apply(len).min()` for
        # every observed sequence (it's a run-wide constant).
        protospacer_min_len: Optional[int] = None,
        surrogate_min_len: Optional[int] = None,
        barcode_min_len: Optional[int] = None,
        ):

    #
    #   FORMAT OBSERVED SEQUENCES AS PANDAS SERIES
    #
    observed_reporter_sequences_indices = ["protospacer"]
    if contains_guide_surrogate:
        observed_reporter_sequences_indices.append("surrogate")
    if contains_guide_barcode:
        observed_reporter_sequences_indices.append("barcode")
    # PERF §3.5: replace per-observation `pd.Series(..., index=...)` construction
    # with a plain dict. Access pattern is by-name only (`seq["protospacer"]`
    # etc.) so Series infrastructure adds no value and constructs a new one per
    # observed read.
    observed_guide_reporter_sequence = dict(zip(observed_reporter_sequences_indices, observed_guide_reporter_sequence_input))
    

    #
    # QUICK VALIDATION OF THE OBSERVED SEQUENCES, SAVE RESULTS AS VARIABLES TO BE USED DOWNSTREAM
    #
    def validate_observed_sequence(sequence_name: str, encoded_whitelist_sequence_series: np.ndarray, missing_info_error: MissingInfoGuideCountError, insufficient_length_error: InsufficientLengthGuideCountError):
        if (observed_guide_reporter_sequence[sequence_name] is None) or (observed_guide_reporter_sequence[sequence_name] == "None") or (observed_guide_reporter_sequence[sequence_name].strip() == ""):
            return missing_info_error
        
        if len(observed_guide_reporter_sequence[sequence_name]) < insufficient_length_error.minimum_length: 
            return insufficient_length_error   
        
        return None
        
    # PERF §3.6: fall back to per-call recomputation only if the caller didn't
    # pre-compute the min (e.g. direct callers outside the main mapping path).
    if protospacer_min_len is None:
        protospacer_min_len = int(whitelist_guide_reporter_df["protospacer"].apply(len).min())
    if contains_guide_surrogate and surrogate_min_len is None:
        surrogate_min_len = int(whitelist_guide_reporter_df["surrogate"].apply(len).min())
    if contains_guide_barcode and barcode_min_len is None:
        barcode_min_len = int(whitelist_guide_reporter_df["barcode"].apply(len).min())

    protospacer_error_result = validate_observed_sequence(sequence_name="protospacer", encoded_whitelist_sequence_series=encoded_whitelist_protospacer_sequences_series, missing_info_error=ProtospacerMissingInfoGuideCountError(sequence_value=observed_guide_reporter_sequence), insufficient_length_error=ProtospacerInsufficientLengthGuideCountError(sequence_length=len(observed_guide_reporter_sequence["protospacer"]), minimum_length=protospacer_min_len))
    surrogate_error_result = None
    if contains_guide_surrogate:
        surrogate_error_result = validate_observed_sequence(sequence_name="surrogate", encoded_whitelist_sequence_series=encoded_whitelist_surrogate_sequences_series, missing_info_error=SurrogateMissingInfoGuideCountError(sequence_value=observed_guide_reporter_sequence), insufficient_length_error=SurrogateInsufficientLengthGuideCountError(sequence_length=len(observed_guide_reporter_sequence["surrogate"]), minimum_length=surrogate_min_len))
    barcode_error_result = None
    if contains_guide_barcode:
        barcode_error_result = validate_observed_sequence(sequence_name="barcode", encoded_whitelist_sequence_series=encoded_whitelist_barcode_sequences_series, missing_info_error=BarcodeMissingInfoGuideCountError(sequence_value=observed_guide_reporter_sequence), insufficient_length_error=BarcodeInsufficientLengthGuideCountError(sequence_length=len(observed_guide_reporter_sequence["barcode"]), minimum_length=barcode_min_len))

    # CHANGE SHAPE OF ENCODINGS TO BE SAME LENGTH
    if protospacer_error_result is None:
        observed_guide_reporter_sequence["protospacer"] = observed_guide_reporter_sequence["protospacer"][:encoded_whitelist_protospacer_sequences_series.shape[1]]
        encoded_whitelist_protospacer_sequences_series = encoded_whitelist_protospacer_sequences_series[:, :len(observed_guide_reporter_sequence["protospacer"]), :] # Change shape of whitelist library encodings
    if contains_guide_surrogate and (surrogate_error_result is None):
        # FIX §1.1: previously clamped observed surrogate to the *protospacer*
        # library length (20 bp) instead of the surrogate length (32 bp),
        # silently discarding 12 bp of every surrogate observation. The
        # consequence was that surrogate Hamming was computed only over
        # positions 0-19 (= 6 bp upstream + first 14 bp of the embedded
        # protospacer) instead of the full 32 bp. Use the correct tensor's
        # shape here and for the barcode clamp below.
        observed_guide_reporter_sequence["surrogate"] = observed_guide_reporter_sequence["surrogate"][:encoded_whitelist_surrogate_sequences_series.shape[1]]
        encoded_whitelist_surrogate_sequences_series = encoded_whitelist_surrogate_sequences_series[:, :len(observed_guide_reporter_sequence["surrogate"]), :]
    if contains_guide_barcode and (barcode_error_result is None):
        # FIX §1.1: same bug as above for the barcode clamp. In practice this
        # was a no-op because barcodes (4 bp) are always shorter than
        # protospacers (20 bp), but fix the cross-reference for consistency
        # and to avoid latent breakage if barcode length ever exceeds the
        # protospacer length.
        observed_guide_reporter_sequence["barcode"] = observed_guide_reporter_sequence["barcode"][:encoded_whitelist_barcode_sequences_series.shape[1]]
        encoded_whitelist_barcode_sequences_series = encoded_whitelist_barcode_sequences_series[:, :len(observed_guide_reporter_sequence["barcode"]), :]
    
    # PERFORM MAPPPING ON BOTH THE PROTOSPACER-ONLY, FULL, AND SURROGATE MISMAPPED-PROTOSPACER.
    complete_match_result = CompleteInferenceMatchResult()

    #
    # GET THE BARCODE-ONLY MATCHES
    # (Since barcode-matching will be done a lot, let's do first and save the match results)
    if contains_guide_barcode: # IF BARCODE IS PARSED
        barcode_available = True # Let's save an indicator on if we successfully parsed barcode for convenience downstream
        if barcode_error_result is None: # IF NO BARCODE PARSING ISSUE

            #
            # FIND BARCODE MATCHES
            #
            try:
                observed_barcode_sequence_encoded = crispr_sequence_encoding.encode_DNA_sequence_observed(observed_guide_reporter_sequence["barcode"])  # Encode the observed barcode
            
                observed_barcode_sequence_dists = crispr_sequence_encoding.retrieve_hamming_distance_whitelist(observed_barcode_sequence_encoded, encoded_whitelist_barcode_sequences_series) # Retrieve hamming distance with whitelist barcode
                observed_barcode_sequence_dists_min = observed_barcode_sequence_dists.min() # Get the barcode with the minimum  hamming distance
                barcode_hamming_threshold_met = (observed_barcode_sequence_dists_min < barcode_hamming_threshold)

                barcode_matches_indices = np.where(observed_barcode_sequence_dists == observed_barcode_sequence_dists_min)[0] # Get the indices of ALL barcode matches
                whitelist_guide_reporter_df_barcode_match = whitelist_guide_reporter_df.iloc[barcode_matches_indices] # Get the whitelist reporters with the matched barcode(s)
                
                if barcode_hamming_threshold_met: # IF BARCODE MATCH
                    # PERF §3.3: previously stashed subset-encoded tensors here
                    # for later Hamming re-computation. After the dedupe, all
                    # subset distances are index-gathered from the full
                    # distance vectors, so the subset tensors are no longer
                    # needed.
                    pass
                else: # NO BARCODE MATCH, ERROR
                    barcode_available = False
                    barcode_error_result = BarcodeHammingThresholdGuideCountError(
                        hamming_min=observed_barcode_sequence_dists_min,
                        hamming_threshold=barcode_hamming_threshold,
                        n_whitelist_candidates=len(whitelist_guide_reporter_df),
                        n_hamming_min_match=len(whitelist_guide_reporter_df_barcode_match))
            except Exception as e:
                barcode_available = False
                barcode_error_result = GuideCountError(guide_count_error_type=GuideCountErrorType.BARCODE_MISC, miscellaneous_info_dict={"message": f"Error encoding barcode sequence: {observed_guide_reporter_sequence['barcode']}; Exception message = {str(e)}"})
        else:
            barcode_available = False
            

    #
    # PERF §3.3: Encode the observed surrogate ONCE and compute its Hamming
    # distance against the FULL whitelist ONCE. All downstream branches
    # (protospacer-match subset, barcode-match subset, barcode+protospacer
    # double subset, whole-whitelist) are index-gathers into this single
    # vector — previously we re-encoded + re-ran Hamming up to 4× per observed
    # read. Hoisting lifts the encode/Hamming cost out of the inner conditional
    # chains without changing results (index-gather is mathematically
    # identical to subsetting the whitelist before Hamming).
    observed_surrogate_sequence_encoded = None
    observed_surrogate_sequence_dists_full = None
    encoding_surrogate_error_result_shared = None
    if contains_guide_surrogate and (surrogate_error_result is None):
        try:
            observed_surrogate_sequence_encoded = crispr_sequence_encoding.encode_DNA_sequence_observed(observed_guide_reporter_sequence["surrogate"])
            observed_surrogate_sequence_dists_full = crispr_sequence_encoding.retrieve_hamming_distance_whitelist(observed_surrogate_sequence_encoded, encoded_whitelist_surrogate_sequences_series)
        except Exception as e:
            encoding_surrogate_error_result_shared = GuideCountError(
                guide_count_error_type=GuideCountErrorType.SURROGATE_MISC,
                miscellaneous_info_dict={"message": f"Error encoding surrogate sequence: {observed_guide_reporter_sequence['surrogate']}; Exception message = {str(e)}"})

    #
    # PREPARE THE PROTOSPACER-ONLYMATCHES
    #
    if protospacer_error_result is None:
        encoding_protospacer_error_result = None
        try:
            observed_protospacer_sequence_encoded = crispr_sequence_encoding.encode_DNA_sequence_observed(observed_guide_reporter_sequence["protospacer"])  # Encode the observed protospacer
            observed_protospacer_sequence_dists = crispr_sequence_encoding.retrieve_hamming_distance_whitelist(observed_protospacer_sequence_encoded, encoded_whitelist_protospacer_sequences_series) # Hamming distance among all whitelist protospacers
            observed_protospacer_sequence_dists_min = observed_protospacer_sequence_dists.min() # Get minimum hamming distance
            protospacer_hamming_threshold_met = (observed_protospacer_sequence_dists_min < protospacer_hamming_threshold)

            protospacer_matches_indices = np.where(observed_protospacer_sequence_dists == observed_protospacer_sequence_dists_min)[0]
            whitelist_guide_reporter_df_hamming_protospacer_match = whitelist_guide_reporter_df.iloc[protospacer_matches_indices] # Get all whitelisted guides with the minimum hamming distance (could be multiple)
        except Exception as e:
            encoding_protospacer_error_result = GuideCountError(guide_count_error_type=GuideCountErrorType.PROTOSPACER_MISC, miscellaneous_info_dict={"message": f"Error encoding protospacer sequence: {observed_guide_reporter_sequence['protospacer']}; Exception message = {str(e)}"})
        
        if encoding_protospacer_error_result is not None:
            complete_match_result.protospacer_match = MatchSetSingleInferenceMatchResult(error=encoding_protospacer_error_result)
            complete_match_result.protospacer_match_surrogate_match_barcode_match = MatchSetSingleInferenceMatchResult(error=encoding_protospacer_error_result)
            complete_match_result.protospacer_match_surrogate_match = MatchSetSingleInferenceMatchResult(error=encoding_protospacer_error_result)
            complete_match_result.protospacer_match_barcode_match = MatchSetSingleInferenceMatchResult(error=encoding_protospacer_error_result)
            complete_match_result.protospacer_mismatch_surrogate_match_barcode_match = SurrogateProtospacerMismatchSingleInferenceMatchResult(error=encoding_protospacer_error_result)
            complete_match_result.protospacer_mismatch_surrogate_match = SurrogateProtospacerMismatchSingleInferenceMatchResult(error=encoding_protospacer_error_result)
        elif protospacer_hamming_threshold_met: # PROTOSPACER MATCHED
            #
            # SET THE PROTOSPACER-ONLYMATCHES
            #
            # PERF §3.3: subset tensor no longer needed (all Hamming subsets
            # are index-gathered from the full distance vectors).
            complete_match_result.protospacer_match = MatchSetSingleInferenceMatchResult(value=MatchSetSingleInferenceMatchResultValue(matches=whitelist_guide_reporter_df_hamming_protospacer_match))
            
            if contains_guide_barcode: # IF BARCODE IS PARSED, PROCEED TO BARCODE-MATCHING INFERENCE
                if barcode_available: # IF BARCODE IS MATCHED, PROCEED TO BARCODE-MATCHING INFERENCE
                    #
                    # PREPARE PROTOSPACER-MATCH, BARCODE-MATCH
                    #
                    # PERF §3.3: gather barcode-subset distances from the
                    # full-whitelist vector instead of re-running Hamming on
                    # the subset-encoded tensor — identical values.
                    observed_protospacer_sequence_dists_barcode_match = observed_protospacer_sequence_dists[barcode_matches_indices]
                    observed_protospacer_sequence_dists_barcode_match_min = observed_protospacer_sequence_dists_barcode_match.min() # Get minimum hamming distance
                    barcode_match_protospacer_hamming_threshold_met = (observed_protospacer_sequence_dists_barcode_match_min < protospacer_hamming_threshold)

                    barcode_match_protospacer_match_indices = np.where(observed_protospacer_sequence_dists_barcode_match == observed_protospacer_sequence_dists_barcode_match_min)[0]
                    whitelist_guide_reporter_df_hamming_barcode_match_protospacer_match = whitelist_guide_reporter_df_barcode_match.iloc[barcode_match_protospacer_match_indices] # Get all whitelisted guides with the minimum hamming distance (could be multiple)
                    if barcode_match_protospacer_hamming_threshold_met:  # PROTOSPACER MATCHED WHEN BARCODE IS SPECIFIED
                        #
                        # SET PROTOSPACER-MATCH, BARCODE-MATCH
                        #
                        
                        complete_match_result.protospacer_match_barcode_match = MatchSetSingleInferenceMatchResult(value=MatchSetSingleInferenceMatchResultValue(matches=whitelist_guide_reporter_df_hamming_barcode_match_protospacer_match))
                        if contains_guide_surrogate: # IF SURROGATE IS PARSED, PROCEED TO SURROGATE-MATCHING INFERENCE
                            if surrogate_error_result is None:
                                #
                                # PREPARE PROTOSPACER-MATCH, BARCODE-MATCH, SURROGATE-MATCH
                                #
                                # PERF §3.3: gather the double-subset surrogate
                                # distances from the pre-computed full-whitelist
                                # vector. barcode_matches_indices indexes into the
                                # full whitelist; barcode_match_protospacer_match_indices
                                # indexes within that subset — compose.
                                encoding_surrogate_error_result = encoding_surrogate_error_result_shared
                                if encoding_surrogate_error_result is None:
                                    combined_indices = barcode_matches_indices[barcode_match_protospacer_match_indices]
                                    observed_surrogate_sequence_dists_barcode_match_protospacer_match = observed_surrogate_sequence_dists_full[combined_indices]
                                    observed_surrogate_sequence_dists_barcode_match_protospacer_match_min = observed_surrogate_sequence_dists_barcode_match_protospacer_match.min()
                                    barcode_match_protospacer_match_surrogate_hamming_threshold_met = (observed_surrogate_sequence_dists_barcode_match_protospacer_match_min < surrogate_hamming_threshold)

                                    barcode_match_protospacer_match_surrogate_match_indices = np.where(observed_surrogate_sequence_dists_barcode_match_protospacer_match == observed_surrogate_sequence_dists_barcode_match_protospacer_match_min)[0]
                                    whitelist_guide_reporter_df_hamming_barcode_match_protospacer_match_surrogate_match = whitelist_guide_reporter_df_hamming_barcode_match_protospacer_match.iloc[barcode_match_protospacer_match_surrogate_match_indices]
                                
                                if encoding_surrogate_error_result is not None:
                                    complete_match_result.protospacer_match_surrogate_match_barcode_match = MatchSetSingleInferenceMatchResult(
                                        error=encoding_surrogate_error_result)
                                elif barcode_match_protospacer_match_surrogate_hamming_threshold_met: # SURROGATE MATCHED
                                    #
                                    # SET PROTOSPACER-MATCH, BARCODE-MATCH, SURROGATE-MATCH
                                    #
                                    complete_match_result.protospacer_match_surrogate_match_barcode_match = MatchSetSingleInferenceMatchResult(value=MatchSetSingleInferenceMatchResultValue(matches=whitelist_guide_reporter_df_hamming_barcode_match_protospacer_match_surrogate_match))
                                else: # NO SURROGATE MATCH, ERROR
                                    complete_match_result.protospacer_match_surrogate_match_barcode_match = MatchSetSingleInferenceMatchResult(
                                        error=SurrogateHammingThresholdGuideCountError(
                                            hamming_min=observed_surrogate_sequence_dists_barcode_match_protospacer_match_min, 
                                            hamming_threshold=surrogate_hamming_threshold,
                                            n_whitelist_candidates=len(whitelist_guide_reporter_df_hamming_barcode_match_protospacer_match),
                                            n_hamming_min_match=len(whitelist_guide_reporter_df_hamming_barcode_match_protospacer_match_surrogate_match),
                                            barcode_subsetted=True))
                            else: # PARSED SURROGATE HAS SOME UNEXPECTED ISSUE 
                                complete_match_result.protospacer_match_surrogate_match_barcode_match = MatchSetSingleInferenceMatchResult(error=surrogate_error_result)
                    else: # NO PROTOSPACER MATCH, ERROR
                        complete_match_result.protospacer_match_barcode_match = MatchSetSingleInferenceMatchResult(
                                        error=ProtospacerHammingThresholdGuideCountError(
                                            hamming_min=observed_protospacer_sequence_dists_barcode_match_min, 
                                            hamming_threshold=protospacer_hamming_threshold,
                                            n_whitelist_candidates=len(whitelist_guide_reporter_df_barcode_match),
                                            n_hamming_min_match=len(whitelist_guide_reporter_df_hamming_barcode_match_protospacer_match),
                                            barcode_subsetted=True))
                else: # NO BARCODE MATCH, SET ERROR TO THE RESULTS REQUIREING A BARCODE MATCH
                    complete_match_result.protospacer_match_barcode_match = MatchSetSingleInferenceMatchResult(error=barcode_error_result)
                    complete_match_result.protospacer_match_surrogate_match_barcode_match = MatchSetSingleInferenceMatchResult(error=barcode_error_result)

            #
            # PROTOSPACER-MATCH, SURROGATE-MATCH, NO BARCODE
            #
            if contains_guide_surrogate: # IF SURROGATE IS PARSED
                if surrogate_error_result is None: # AND NO SURROGATE PARSING ISSUE 

                    #
                    # GET SURROGATE-MATCHES ON THE PROTOSPACER-MATCHES
                    #
                    # PERF §3.3: gather the protospacer-subset surrogate
                    # distances from the pre-computed full-whitelist vector.
                    encoding_surrogate_error_result = encoding_surrogate_error_result_shared
                    if encoding_surrogate_error_result is None:
                        observed_surrogate_sequence_dists_protospacer_match = observed_surrogate_sequence_dists_full[protospacer_matches_indices]
                        observed_surrogate_sequence_dists_protospacer_match_min = observed_surrogate_sequence_dists_protospacer_match.min()
                        protospacer_match_surrogate_hamming_threshold_met = (observed_surrogate_sequence_dists_protospacer_match_min < surrogate_hamming_threshold)

                        protospacer_match_surrogate_match_indices = np.where(observed_surrogate_sequence_dists_protospacer_match == observed_surrogate_sequence_dists_protospacer_match_min)[0]
                        whitelist_guide_reporter_df_hamming_protospacer_match_surrogate_match = whitelist_guide_reporter_df_hamming_protospacer_match.iloc[protospacer_match_surrogate_match_indices]
                    
                    if encoding_surrogate_error_result is not None:
                        complete_match_result.protospacer_match_surrogate_match = MatchSetSingleInferenceMatchResult(
                                        error=encoding_surrogate_error_result)
                    elif protospacer_match_surrogate_hamming_threshold_met: # IF SURROGATE MATCH
                        complete_match_result.protospacer_match_surrogate_match = MatchSetSingleInferenceMatchResult(value=MatchSetSingleInferenceMatchResultValue(matches=whitelist_guide_reporter_df_hamming_protospacer_match_surrogate_match))
                    else: # NO SURROGATE MATCH, ERROR OUT
                        complete_match_result.protospacer_match_surrogate_match = MatchSetSingleInferenceMatchResult(
                                        error=SurrogateHammingThresholdGuideCountError(
                                            hamming_min=observed_surrogate_sequence_dists_protospacer_match_min, 
                                            hamming_threshold=surrogate_hamming_threshold,
                                            n_whitelist_candidates=len(whitelist_guide_reporter_df_hamming_protospacer_match),
                                            n_hamming_min_match=len(whitelist_guide_reporter_df_hamming_protospacer_match_surrogate_match),
                                            barcode_subsetted=False))
                else: # SURROGATE PARSING ISSUE, ERROR
                    complete_match_result.protospacer_match_surrogate_match = MatchSetSingleInferenceMatchResult(error=surrogate_error_result)



            #
            # GET THE PROTOSPACER-SURROGATE MISMATCHES NOW (if surrogate is available)
            #
            if contains_guide_surrogate: # IF SURROGATE IS PARSED
                if surrogate_error_result is None: # AND NO SURROGATE PARSING ISSUE
                    if contains_guide_barcode: # IF BARCODE IS PARSED
                        if barcode_available: # AND NO BARCODE PARSING ISSUE
                            #
                            # PREPARE SURROGATE-MATCH, BARCODE-MATCH
                            #
                            # PERF §3.3: gather barcode-subset surrogate dists
                            # from the pre-computed full-whitelist vector.
                            encoding_surrogate_error_result = encoding_surrogate_error_result_shared
                            if encoding_surrogate_error_result is None:
                                observed_surrogate_sequence_dists_barcode_match = observed_surrogate_sequence_dists_full[barcode_matches_indices]
                                observed_surrogate_sequence_dists_barcode_match_min = observed_surrogate_sequence_dists_barcode_match.min()
                                barcode_match_surrogate_hamming_threshold_met = (observed_surrogate_sequence_dists_barcode_match_min < surrogate_hamming_threshold)

                                barcode_match_surrogate_match_indices = np.where(observed_surrogate_sequence_dists_barcode_match == observed_surrogate_sequence_dists_barcode_match_min)[0]
                                whitelist_guide_reporter_df_hamming_barcode_match_surrogate_match = whitelist_guide_reporter_df_barcode_match.iloc[barcode_match_surrogate_match_indices]
                            
                            if encoding_surrogate_error_result is not None:
                                complete_match_result.protospacer_mismatch_surrogate_match_barcode_match = SurrogateProtospacerMismatchSingleInferenceMatchResult(
                                    error=encoding_surrogate_error_result)
                            if barcode_match_surrogate_hamming_threshold_met: # If the surrogate matches are below the surrogate hamming threshold
                                # PERF §3.4: replace pd.merge (heavy pandas overhead, per observed sequence) with a
                                # hash-based row intersection — the merge was an inner-join that only needed to
                                # identify shared rows between the two small DataFrames.
                                _proto_match_set = {tuple(r) for r in whitelist_guide_reporter_df_hamming_protospacer_match.itertuples(index=False, name=None)}
                                _mask = [tuple(r) in _proto_match_set for r in whitelist_guide_reporter_df_hamming_barcode_match_surrogate_match.itertuples(index=False, name=None)]
                                whitelist_guide_reporter_df_hamming_barcode_match_surrogate_match_protospacer_match = whitelist_guide_reporter_df_hamming_barcode_match_surrogate_match[_mask]
                                
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
                                        n_whitelist_candidates=len(whitelist_guide_reporter_df_barcode_match),
                                        n_hamming_min_match=len(whitelist_guide_reporter_df_hamming_barcode_match_surrogate_match),
                                        barcode_subsetted=True,
                                    ))
                        else: # BARCODE PARSING ISSUE OR NO MATCH, ERROR OUT
                            complete_match_result.protospacer_mismatch_surrogate_match_barcode_match = SurrogateProtospacerMismatchSingleInferenceMatchResult(error=barcode_error_result)
                    #
                    # PREPARE SURROGATE MATCHES ON THE PROTOSPACER-MATCHES SEQUENCES
                    #
                    # PERF §3.3: reuse the pre-computed full-whitelist
                    # surrogate distance vector — no per-observation re-encode.
                    encoding_surrogate_error_result = encoding_surrogate_error_result_shared
                    if encoding_surrogate_error_result is None:
                        observed_surrogate_sequence_dists = observed_surrogate_sequence_dists_full
                        observed_surrogate_sequence_dists_min = observed_surrogate_sequence_dists.min()
                        surrogate_hamming_threshold_met = (observed_surrogate_sequence_dists_min < surrogate_hamming_threshold)

                        surrogate_match_indices = np.where(observed_surrogate_sequence_dists == observed_surrogate_sequence_dists_min)[0]
                        whitelist_guide_reporter_df_hamming_surrogate_match = whitelist_guide_reporter_df.iloc[surrogate_match_indices]
                    
                    if encoding_surrogate_error_result is not None:
                        complete_match_result.protospacer_mismatch_surrogate_match = SurrogateProtospacerMismatchSingleInferenceMatchResult(
                            error=encoding_surrogate_error_result)
                    elif surrogate_hamming_threshold_met: # IF SURROGATE MATCHED

                        # PERF §3.4: same hash-intersection replacement for pd.merge.
                        _proto_match_set = {tuple(r) for r in whitelist_guide_reporter_df_hamming_protospacer_match.itertuples(index=False, name=None)}
                        _mask = [tuple(r) in _proto_match_set for r in whitelist_guide_reporter_df_hamming_surrogate_match.itertuples(index=False, name=None)]
                        whitelist_guide_reporter_df_hamming_surrogate_match_protospacer_match = whitelist_guide_reporter_df_hamming_surrogate_match[_mask]
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
                                n_whitelist_candidates=len(whitelist_guide_reporter_df),
                                n_hamming_min_match=len(whitelist_guide_reporter_df_hamming_surrogate_match),
                                barcode_subsetted=False
                            ))
                else: # SURROGATE PARSING ISSUE, ERROR
                    complete_match_result.protospacer_mismatch_surrogate_match_barcode_match = SurrogateProtospacerMismatchSingleInferenceMatchResult(error=surrogate_error_result)
                    complete_match_result.protospacer_mismatch_surrogate_match = SurrogateProtospacerMismatchSingleInferenceMatchResult(error=surrogate_error_result)
                
        else: # NO PROTOSPACER MATCH, ERROR
            error_result = ProtospacerHammingThresholdGuideCountError(
                hamming_min=observed_protospacer_sequence_dists_min,
                hamming_threshold=protospacer_hamming_threshold,
                n_whitelist_candidates=len(whitelist_guide_reporter_df),
                n_hamming_min_match=len(whitelist_guide_reporter_df_hamming_protospacer_match),
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
            try:
                current_guide_sequence_encoded = crispr_sequence_encoding.encode_DNA_sequence_observed(current_guide_sequence) 
            except Exception as e:
                _log.error(e)
                raise Exception(f"Error encoding surrogate sequence: {current_guide_sequence}")
            # Calculate the hamming distance of the mutated guide to all other guides
            hamming_distances =  crispr_sequence_encoding.retrieve_hamming_distance_whitelist(current_guide_sequence_encoded, encoded_whitelist_sequences_series)
            if len(np.where(hamming_distances == hamming_distances.min())[0]) > 1: # TODO (3/30/23): Can also add to the conditional whether the minimum hamming distance guide is still the original guide - but it is probably rare cases where this would not be the case, so not too important to implement
                mutation_count_until_nonunique.append(iteration+1)
                break   
    mutation_count_until_nonunique = pd.Series(mutation_count_until_nonunique)
    
    # From the mutation count distribution, calculate the threshold based on the provided quantile.
    hamming_distance_threshold = mutation_count_until_nonunique.quantile(quantile)
    return int(hamming_distance_threshold) + 1 # Take ceil of result (to be conservative)

