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
from datetime import datetime
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


# PERF §2.9: a module-level factory avoids `defaultdict(lambda: defaultdict(dict))`
# — the outer lambda closes over the inner lambda, so each missing key call
# constructs a fresh lambda stack. A named factory is equivalent but also
# pickleable (the previous lambda broke pickling if the outer dict was ever
# persisted directly).
def _inference_dict_factory():
    return defaultdict(dict)
from .crispr_count_processing import get_counterseries_all_results
from ..quality_control.crispr_mapping_quality_control import perform_counts_quality_control
from ..models.mapping_models import GeneralGuideCountType, GeneralMappingInferenceDict
from ..models.mapping_models import AllMatchSetWhitelistReporterCounterSeriesResults, WhitelistReporterCountsResult, InferenceResult, CountInput, QualityControlResult
import logging
_log = logging.getLogger(__name__)



# TODO: There will probably be some type errors with the DefaultDict when testing on non UMI (since it requires CounterType), so make sure to test with different variations of inputs
@typechecked
def get_whitelist_reporter_counts_with_umi(observed_guide_reporter_umi_counts: GeneralGuideCountType, 
                                           whitelist_guide_reporter_df: Optional[pd.DataFrame], 
                                           contains_guide_surrogate:bool = False, 
                                           contains_guide_barcode:bool = False, 
                                           contains_guide_umi:bool = False, 
                                           contains_sample_barcode:bool = False, 
                                           protospacer_hamming_threshold_strict: Optional[int] = 7, 
                                           surrogate_hamming_threshold_strict: Optional[int] = 10,
                                           guide_barcode_hamming_threshold_strict: Optional[int] = 2, 
                                           retain_inference_results: bool = False,
                                           cores: int=1) -> WhitelistReporterCountsResult:
    
    # §4.4: validate whitelist DF columns up front so misnamed columns fail
    # with a clear message instead of deep inside the encoding loop.
    if whitelist_guide_reporter_df is not None:
        _required_cols = {"protospacer"}
        if contains_guide_surrogate:
            _required_cols.add("surrogate")
        if contains_guide_barcode:
            _required_cols.add("barcode")
        _missing = _required_cols - set(whitelist_guide_reporter_df.columns)
        if _missing:
            raise ValueError(
                f"whitelist_guide_reporter_df is missing required columns: {sorted(_missing)}. "
                f"Expected {sorted(_required_cols)} given contains_guide_surrogate={contains_guide_surrogate}, "
                f"contains_guide_barcode={contains_guide_barcode}. Found columns: {list(whitelist_guide_reporter_df.columns)}."
            )

    # Generate whitelist dataframe based on all observed sequences if none provided
    if whitelist_guide_reporter_df is None:
        whitelist_dataframe_input = {}
        whitelist_dataframe_input["protospacer"] = [observed_sequence_tuple[0] for observed_sequence_tuple in observed_guide_reporter_umi_counts.keys()] # protospacer always in 0th index
        if contains_guide_surrogate:
            whitelist_dataframe_input["surrogate"] = [observed_sequence_tuple[1] for observed_sequence_tuple in observed_guide_reporter_umi_counts.keys()] # surrogate always in 1st index
            if contains_guide_barcode:
                whitelist_dataframe_input["barcode"] = [observed_sequence_tuple[2] for observed_sequence_tuple in observed_guide_reporter_umi_counts.keys()] # barcode in 2nd index if surrogate provided
        else:
            if contains_guide_barcode:
                whitelist_dataframe_input["barcode"] = [observed_sequence_tuple[1] for observed_sequence_tuple in observed_guide_reporter_umi_counts.keys()] # barcode in 1st index if surrogate not provided

        whitelist_guide_reporter_df = pd.DataFrame(whitelist_dataframe_input)


    # Strip all sequences:
    def strip_series(series):
        return series.apply(lambda item : item.rstrip())
    whitelist_guide_reporter_df = whitelist_guide_reporter_df.apply(strip_series, axis=0)

    # Temporary bug fix. Pad sequences so they are all of same length - encoding numpy matrices requires consistent shape. Still pass the original guide table for selecting the matches.     
    def pad_series(series):
        max_surrogate_len = series.apply(len).max()
        return series.apply(lambda item: item.ljust(max_surrogate_len, 'X'))
    padded_whitelist_guide_reporter_df = whitelist_guide_reporter_df.apply(pad_series, axis=0)
    

    #
    # ENCODE THE WHITELISTED SEQUENCES INTO NUMPY MATRICES - THIS IS REQUIRED FOR HAMMING-BASED MAPPING
    #
    encoded_whitelist_protospacer_sequences_series = crispr_sequence_encoding.encode_guide_series_whitelist(padded_whitelist_guide_reporter_df["protospacer"])
    encoded_whitelist_surrogate_sequences_series = None
    if contains_guide_surrogate:
        encoded_whitelist_surrogate_sequences_series = crispr_sequence_encoding.encode_guide_series_whitelist(padded_whitelist_guide_reporter_df["surrogate"])
    
    encoded_whitelist_barcode_sequences_series = None
    if contains_guide_barcode:
        encoded_whitelist_barcode_sequences_series = crispr_sequence_encoding.encode_guide_series_whitelist(padded_whitelist_guide_reporter_df["barcode"])

    

    
    #
    #    SET THE PROTOSPACER HAMMING THRESHOLD
    #
    # NOTE: Could there be an even more dynamic threshold, that certain mapped guides can have a higher threshold (i.e. those with many editable bases) and other guides can have low hamming (i.e. those with not many editable bases)
    protospacer_hamming_threshold_dynamic = False
    protospacer_hamming_threshold: int = protospacer_hamming_threshold_strict
    if protospacer_hamming_threshold_strict is None:
        protospacer_hamming_threshold_dynamic = True
        # TODO: Pass in arguments to set this hamming threshold. 
        protospacer_hamming_threshold: int = crispr_guide_inference.determine_hamming_threshold(whitelist_guide_reporter_df["protospacer"], encoded_whitelist_protospacer_sequences_series, sample_count = 100, quantile = 0.05)
        
    _log.info("Protospacer hamming threshold is " + str(protospacer_hamming_threshold))

    #
    #   SET THE SURROGATE HAMMING THRESHOLD
    #
    surrogate_hamming_threshold: Optional[int] = surrogate_hamming_threshold_strict
    if contains_guide_surrogate:
        surrogate_hamming_threshold_dynamic = False
        if surrogate_hamming_threshold_strict is None:
            surrogate_hamming_threshold_dynamic = True
            # TODO: Pass in arguments to set this hamming threshold. 
            surrogate_hamming_threshold: int = crispr_guide_inference.determine_hamming_threshold(whitelist_guide_reporter_df["surrogate"], encoded_whitelist_surrogate_sequences_series, sample_count = 100, quantile = 0.05)
        _log.info("Surrogate hamming threshold is " + str(surrogate_hamming_threshold))

    #
    #   SET THE BARCODE HAMMING THRESHOLD
    #
    guide_barcode_hamming_threshold: Optional[int] = guide_barcode_hamming_threshold_strict
    if contains_guide_barcode:
        guide_barcode_hamming_threshold_dynamic = False
        if guide_barcode_hamming_threshold_strict is  None:
            guide_barcode_hamming_threshold_dynamic = True
            guide_barcode_hamming_threshold: int = crispr_guide_inference.determine_hamming_threshold(whitelist_guide_reporter_df["barcode"], encoded_whitelist_barcode_sequences_series, sample_count = 100, quantile = 0.05)
        _log.info("Barcode hamming threshold is " + str(guide_barcode_hamming_threshold))



    #
    #   FROM THE OBSERVED SEQUENCES, INFER THE WHITELIST SEQUENCE
    #
    _log.info("Inferring the true guides from observed guides")

    # PERF §3.6: precompute the whitelist per-column min lengths once so
    # infer_whitelist_sequence doesn't recompute via
    # `whitelist_df[col].apply(len).min()` on every observed sequence.
    protospacer_min_len = int(whitelist_guide_reporter_df["protospacer"].apply(len).min())
    surrogate_min_len = int(whitelist_guide_reporter_df["surrogate"].apply(len).min()) if contains_guide_surrogate else None
    barcode_min_len = int(whitelist_guide_reporter_df["barcode"].apply(len).min()) if contains_guide_barcode else None

    infer_whitelist_sequence_p = partial(crispr_guide_inference.infer_whitelist_sequence,
            whitelist_guide_reporter_df=whitelist_guide_reporter_df,
            contains_guide_surrogate=contains_guide_surrogate,
            contains_guide_barcode=contains_guide_barcode,
            contains_guide_umi=contains_guide_umi,
            encoded_whitelist_protospacer_sequences_series=encoded_whitelist_protospacer_sequences_series,
            encoded_whitelist_surrogate_sequences_series=encoded_whitelist_surrogate_sequences_series,
            encoded_whitelist_barcode_sequences_series=encoded_whitelist_barcode_sequences_series,
            protospacer_hamming_threshold=protospacer_hamming_threshold,
            surrogate_hamming_threshold=surrogate_hamming_threshold,
            barcode_hamming_threshold=guide_barcode_hamming_threshold,
            protospacer_min_len=protospacer_min_len,
            surrogate_min_len=surrogate_min_len,
            barcode_min_len=barcode_min_len)

    # Perform inference: frpom the observed sequences (previously parsed), infer the true sequence from the whitelist DF.
    observed_guide_reporter_list = observed_guide_reporter_umi_counts.keys()
    inferred_true_reporter_sequences = None
    before_inference_time = datetime.now()
    if cores > 1:
        _log.info(f"Running inference parallelized on {len(observed_guide_reporter_list)} observed seqeunces with cores {cores}")
        # PERF §3.2: pool.map's default chunksize=1 means one IPC round-trip
        # per observed sequence, which dominates wall time for small
        # per-call work. Amortize by chunking; aim for ~4 chunks per core so
        # the pool stays busy while keeping per-chunk overhead low.
        chunksize = max(1, len(observed_guide_reporter_list) // (cores * 4))
        with Pool(cores) as pool:
            _imap = pool.imap(infer_whitelist_sequence_p, observed_guide_reporter_list, chunksize=chunksize)
            # §4.8: optional tqdm progress bar. tqdm is a soft dep; fall back
            # to the plain iterable when it's not installed.
            try:
                from tqdm.auto import tqdm as _tqdm
                _imap = _tqdm(_imap, total=len(observed_guide_reporter_list), desc="inference", disable=None)
            except ImportError:
                pass
            inferred_true_reporter_sequences = list(_imap)
    else:
        _log.info(f"Running inference on {len(observed_guide_reporter_list)} observed seqeunces non-parallelized")
        inferred_true_reporter_sequences = [infer_whitelist_sequence_p(observed_guide_reporter) for observed_guide_reporter in observed_guide_reporter_list]
    
    after_inference_time = datetime.now()
    _log.info(f"{(after_inference_time-before_inference_time).seconds} seconds for inference")


    _log.info(f"Mapping inference results of length {len(inferred_true_reporter_sequences)} to the result object")
    # Some organization: Map the inferred result of each observed sequence to a dict with the inferred result and correspoding count
    


    # NOTE 20251031: This may be able to be removed
    if contains_sample_barcode:
        observed_guide_reporter_umi_counts_inferred_all_samples: DefaultDict[str, GeneralMappingInferenceDict] = defaultdict(_inference_dict_factory)
        
        # Add all cell_barcodes
        for observed_guide_reporter_key_index, observed_guide_reporter_key in enumerate(observed_guide_reporter_list): # Iterate through each observed guide key
            observed_guide_reporter_cell_counts = observed_guide_reporter_umi_counts[observed_guide_reporter_key]
            observed_cell_barcodes = observed_guide_reporter_cell_counts.keys()
            for cell_barcode in observed_cell_barcodes:
                observed_guide_reporter_umi_counts_inferred_all_samples[cell_barcode][observed_guide_reporter_key] = InferenceResult(
                    observed_value=observed_guide_reporter_cell_counts[cell_barcode], # Add the count to the cell_barcode for the particular guide
                    inferred_value=inferred_true_reporter_sequences[observed_guide_reporter_key_index]
                )

        observed_guide_reporter_umi_counts_inferred_per_sample: GeneralMappingInferenceDict = defaultdict(dict)
    
        after_inference_processing_time = datetime.now()
        _log.info(f"{(after_inference_processing_time-after_inference_time).seconds} seconds for inference processing")
        _log.info("Completed inference")


        # GET THE MAPPED COUNT SERIES BASED ON THE INFERENCE RESULTS
        _log.info("Prepare the processed count series ")
        all_cell_barcodes: List[str] = list(observed_guide_reporter_umi_counts_inferred_all_samples.keys())
        
        all_match_set_whitelist_reporter_counter_series_results: AllMatchSetWhitelistReporterCounterSeriesResults
        quality_control_result: QualityControlResult
        
        all_match_set_whitelist_reporter_counter_series_results = get_counterseries_all_results(observed_guide_reporter_umi_counts_inferred_all_samples, whitelist_guide_reporter_df, contains_guide_barcode, contains_guide_surrogate, contains_guide_umi, contains_sample_barcode)
        quality_control_result: QualityControlResult = perform_counts_quality_control(observed_guide_reporter_umi_counts_inferred_all_samples, contains_guide_umi, contains_guide_surrogate, contains_guide_barcode, contains_sample_barcode)

        count_input= CountInput(whitelist_guide_reporter_df=whitelist_guide_reporter_df,
            contains_guide_surrogate=contains_guide_surrogate,
            contains_guide_barcode=contains_guide_barcode,
            contains_guide_umi=contains_guide_umi,
            contains_sample_barcode=contains_sample_barcode,
            protospacer_hamming_threshold_strict=protospacer_hamming_threshold,
            surrogate_hamming_threshold_strict=surrogate_hamming_threshold,
            guide_barcode_hamming_threshold_strict=guide_barcode_hamming_threshold)
        
        # MEM: by default drop the heavyweight per-observation inference dict —
        # it is only needed by allele/mutation post-processing, which users
        # opt into by passing retain_inference_results=True.
        return WhitelistReporterCountsResult(
            all_match_set_whitelist_reporter_counter_series_results=all_match_set_whitelist_reporter_counter_series_results,
            observed_guide_reporter_umi_counts_inferred=(
                observed_guide_reporter_umi_counts_inferred_all_samples if retain_inference_results else None
            ),
            quality_control_result=quality_control_result,
            count_input=count_input,
        )
    else:    
        
        observed_guide_reporter_umi_counts_inferred: GeneralMappingInferenceDict = defaultdict(dict)
        for observed_guide_reporter_key_index, observed_guide_reporter_key in enumerate(observed_guide_reporter_list): # Iterate through each observed guide key
            observed_guide_reporter_umi_counts_inferred[observed_guide_reporter_key] = InferenceResult(
                observed_value=observed_guide_reporter_umi_counts[observed_guide_reporter_key], # Store the observed 
                inferred_value=inferred_true_reporter_sequences[observed_guide_reporter_key_index]
            )
    
        after_inference_processing_time = datetime.now()
        _log.info(f"{(after_inference_processing_time-after_inference_time).seconds} seconds for inference processing")
        _log.info("Completed inference")


        # GET THE MAPPED COUNT SERIES BASED ON THE INFERENCE RESULTS
        _log.info("Prepare the processed count series ")
        # Count
        all_match_set_whitelist_reporter_counter_series_results = get_counterseries_all_results(observed_guide_reporter_umi_counts_inferred, whitelist_guide_reporter_df, contains_guide_barcode, contains_guide_surrogate, contains_guide_umi, contains_sample_barcode)

        after_counterseries_time = datetime.now()
        _log.info(f"{(after_counterseries_time-after_inference_processing_time).seconds} seconds for counter series generation")

        _log.info("Preparing quality control")
        quality_control_result: QualityControlResult = perform_counts_quality_control(observed_guide_reporter_umi_counts_inferred, contains_guide_umi, contains_guide_surrogate, contains_guide_barcode, contains_sample_barcode)
        
        after_qualitycontrol_time = datetime.now()
        _log.info(f"{(after_qualitycontrol_time-after_counterseries_time).seconds} seconds for quality control")

        count_input = CountInput(whitelist_guide_reporter_df=whitelist_guide_reporter_df,
                contains_guide_surrogate=contains_guide_surrogate,
                contains_guide_barcode=contains_guide_barcode,
                contains_guide_umi=contains_guide_umi,
                contains_sample_barcode=contains_sample_barcode,
                protospacer_hamming_threshold_strict=protospacer_hamming_threshold,
                surrogate_hamming_threshold_strict=surrogate_hamming_threshold,
                guide_barcode_hamming_threshold_strict=guide_barcode_hamming_threshold)
        
        # MEM: drop the inference dict by default — see note above.
        return WhitelistReporterCountsResult(
            all_match_set_whitelist_reporter_counter_series_results=all_match_set_whitelist_reporter_counter_series_results,
            observed_guide_reporter_umi_counts_inferred=(
                observed_guide_reporter_umi_counts_inferred if retain_inference_results else None
            ),
            quality_control_result=quality_control_result,
            count_input=count_input,
        )