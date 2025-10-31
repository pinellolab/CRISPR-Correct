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
from .crispr_count_processing import get_counterseries_all_results
from ..quality_control.crispr_mapping_quality_control import perform_counts_quality_control
from ..models.mapping_models import GeneralGuideCountType, GeneralMappingInferenceDict
from ..models.mapping_models import AllMatchSetWhitelistReporterCounterSeriesResults, WhitelistReporterCountsResult, SampleWhitelistReporterCountsResult, InferenceResult, CountInput


# TODO: There will probably be some type errors with the DefaultDict when testing on non UMI (since it requires CounterType), so make sure to test with different variations of inputs
@typechecked
def get_whitelist_reporter_counts_with_umi(observed_guide_reporter_umi_counts: GeneralGuideCountType, 
                                           whitelist_guide_reporter_df: Optional[pd.DataFrame], 
                                           contains_guide_surrogate:bool = False, 
                                           contains_guide_barcode:bool = False, 
                                           contains_guide_umi:bool = False, 
                                           contains_sample_barcode:bool = False, 
                                           protospacer_hamming_threshold_strict: Optional[int] = 7, 
                                           surrogate_hamming_threshold_strict: Optional[int] = 2, 
                                           guide_barcode_hamming_threshold_strict: Optional[int] = 2, 
                                           store_intermediates: bool = False,
                                           cores: int=1) -> Union[WhitelistReporterCountsResult, SampleWhitelistReporterCountsResult]:
    
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
        
    print("Protospacer hamming threshold is " + str(protospacer_hamming_threshold))

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
        print("Surrogate hamming threshold is " + str(surrogate_hamming_threshold))

    #
    #   SET THE BARCODE HAMMING THRESHOLD
    #
    guide_barcode_hamming_threshold: Optional[int] = guide_barcode_hamming_threshold_strict
    if contains_guide_barcode:
        guide_barcode_hamming_threshold_dynamic = False
        if guide_barcode_hamming_threshold_strict is  None:
            guide_barcode_hamming_threshold_dynamic = True
            guide_barcode_hamming_threshold: int = crispr_guide_inference.determine_hamming_threshold(whitelist_guide_reporter_df["barcode"], encoded_whitelist_barcode_sequences_series, sample_count = 100, quantile = 0.05)
        print("Barcode hamming threshold is " + str(guide_barcode_hamming_threshold))



    #
    #   FROM THE OBSERVED SEQUENCES, INFER THE WHITELIST SEQUENCE
    #
    print("Inferring the true guides from observed guides")

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
            store_intermediates=store_intermediates)

    # Perform inference: frpom the observed sequences (previously parsed), infer the true sequence from the whitelist DF.
    observed_guide_reporter_list = observed_guide_reporter_umi_counts.keys()
    inferred_true_reporter_sequences = None
    before_inference_time = datetime.now()
    if cores > 1:
        print(f"Running inference parallelized on {len(observed_guide_reporter_list)} observed seqeunces with cores {cores}")
        with Pool(cores) as pool:
            inferred_true_reporter_sequences = pool.map(
            infer_whitelist_sequence_p,
            observed_guide_reporter_list
           )
    else:
        print(f"Running inference on {len(observed_guide_reporter_list)} observed seqeunces non-parallelized")
        inferred_true_reporter_sequences = [infer_whitelist_sequence_p(observed_guide_reporter) for observed_guide_reporter in observed_guide_reporter_list]
    
    after_inference_time = datetime.now()
    print(f"{(after_inference_time-before_inference_time).seconds} seconds for inference")


    print(f"Mapping inference results of length {len(inferred_true_reporter_sequences)} to the result object")
    # Some organization: Map the inferred result of each observed sequence to a dict with the inferred result and correspoding count
    
    if contains_sample_barcode:
        observed_guide_reporter_umi_counts_inferred_all_samples: DefaultDict[str, GeneralMappingInferenceDict] = defaultdict(GeneralMappingInferenceDict)
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
        print(f"{(after_inference_processing_time-after_inference_time).seconds} seconds for inference processing")
        print("Completed inference")


        # GET THE MAPPED COUNT SERIES BASED ON THE INFERENCE RESULTS
        print("Prepare the processed count series ")
        all_cell_barcodes: List[str] = observed_guide_reporter_umi_counts_inferred_all_samples.keys()
        all_match_set_whitelist_reporter_counter_series_results_all_samples: DefaultDict[str, AllMatchSetWhitelistReporterCounterSeriesResults]
        quality_control_result_all_samples: DefaultDict[str, GeneralMappingInferenceDict]
        for cell_barcode in all_cell_barcodes:
            observed_guide_reporter_umi_counts_inferred_per_sample = observed_guide_reporter_umi_counts_inferred_all_samples[cell_barcode]
            all_match_set_whitelist_reporter_counter_series_results_per_sample = get_counterseries_all_results(observed_guide_reporter_umi_counts_inferred_per_sample, whitelist_guide_reporter_df, contains_guide_barcode, contains_guide_surrogate, contains_guide_umi)
            quality_control_result_per_sample = perform_counts_quality_control(observed_guide_reporter_umi_counts_inferred_per_sample, contains_guide_umi, contains_guide_surrogate, contains_guide_barcode)

            all_match_set_whitelist_reporter_counter_series_results_all_samples[cell_barcode] = all_match_set_whitelist_reporter_counter_series_results_per_sample
            quality_control_result_all_samples[cell_barcode] = quality_control_result_per_sample


        count_input= CountInput(whitelist_guide_reporter_df=whitelist_guide_reporter_df,
            contains_surrogate=contains_guide_surrogate,
            contains_guide_barcode=contains_guide_barcode,
            contains_guide_umi=contains_guide_umi,
            contains_sample_barcode=contains_sample_barcode,
            protospacer_hamming_threshold_strict=protospacer_hamming_threshold,
            surrogate_hamming_threshold_strict=surrogate_hamming_threshold,
            guide_barcode_hamming_threshold_strict=guide_barcode_hamming_threshold)
        
        return SampleWhitelistReporterCountsResult(all_match_set_whitelist_reporter_counter_series_results_all_samples=all_match_set_whitelist_reporter_counter_series_results_all_samples,
                                                    observed_guide_reporter_umi_counts_inferred_all_samples=observed_guide_reporter_umi_counts_inferred_all_samples, 
                                                    quality_control_result_all_samples=quality_control_result_all_samples, 
                                                    count_input=count_input)
    else:    
        
        observed_guide_reporter_umi_counts_inferred: GeneralMappingInferenceDict = defaultdict(dict)
        for observed_guide_reporter_key_index, observed_guide_reporter_key in enumerate(observed_guide_reporter_list): # Iterate through each observed guide key
            observed_guide_reporter_umi_counts_inferred[observed_guide_reporter_key] = InferenceResult(
                observed_value=observed_guide_reporter_umi_counts[observed_guide_reporter_key], # Store the observed 
                inferred_value=inferred_true_reporter_sequences[observed_guide_reporter_key_index]
            )
    
        after_inference_processing_time = datetime.now()
        print(f"{(after_inference_processing_time-after_inference_time).seconds} seconds for inference processing")
        print("Completed inference")


        # GET THE MAPPED COUNT SERIES BASED ON THE INFERENCE RESULTS
        print("Prepare the processed count series ")
        # Count
        all_match_set_whitelist_reporter_counter_series_results = get_counterseries_all_results(observed_guide_reporter_umi_counts_inferred, whitelist_guide_reporter_df, contains_guide_barcode, contains_guide_surrogate, contains_guide_umi)

        after_counterseries_time = datetime.now()
        print(f"{(after_counterseries_time-after_inference_processing_time).seconds} seconds for counter series generation")

        print("Preparing quality control")
        quality_control_result = perform_counts_quality_control(observed_guide_reporter_umi_counts_inferred, contains_guide_umi, contains_guide_surrogate, contains_guide_barcode)
        
        after_qualitycontrol_time = datetime.now()
        print(f"{(after_qualitycontrol_time-after_counterseries_time).seconds} seconds for quality control")

        count_input = CountInput(whitelist_guide_reporter_df=whitelist_guide_reporter_df,
                contains_surrogate=contains_guide_surrogate,
                contains_guide_barcode=contains_guide_barcode,
                contains_guide_umi=contains_guide_umi,
                protospacer_hamming_threshold_strict=protospacer_hamming_threshold,
                surrogate_hamming_threshold_strict=surrogate_hamming_threshold,
                guide_barcode_hamming_threshold_strict=guide_barcode_hamming_threshold)
        
        return WhitelistReporterCountsResult(all_match_set_whitelist_reporter_counter_series_results=all_match_set_whitelist_reporter_counter_series_results, observed_guide_reporter_umi_counts_inferred=observed_guide_reporter_umi_counts_inferred, quality_control_result=quality_control_result, count_input=count_input)