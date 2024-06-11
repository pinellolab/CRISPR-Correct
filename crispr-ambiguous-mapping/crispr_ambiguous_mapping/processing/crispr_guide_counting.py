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
from ..models.mapping_models import *
from ..models.editing_models import *
from ..models.types import *


# TODO: There will probably be some type errors with the DefaultDict when testing on non UMI (since it requires CounterType), so make sure to test with different variations of inputs
@typechecked
def get_whitelist_reporter_counts_with_umi(observed_guide_reporter_umi_counts: GeneralGuideCountType, 
                                           whitelist_guide_reporter_df: pd.DataFrame, 
                                           contains_surrogate:bool = False, 
                                           contains_barcode:bool = False, 
                                           contains_umi:bool = False, 
                                           protospacer_hamming_threshold_strict: Optional[int] = 7, 
                                           surrogate_hamming_threshold_strict: Optional[int] = 2, 
                                           barcode_hamming_threshold_strict: Optional[int] = 2, cores: int=1):
    
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
    observed_guide_reporter_umi_counts_inferred: DefaultDict[Tuple[str,Optional[str],Optional[str]], dict] = defaultdict(dict)
    for observed_guide_reporter_key_index, observed_guide_reporter_key in enumerate(observed_guide_reporter_list):
        observed_guide_reporter_umi_counts_inferred[observed_guide_reporter_key] = InferenceResult(
            observed_value=observed_guide_reporter_umi_counts[observed_guide_reporter_key],
            inferred_value=inferred_true_reporter_sequences[observed_guide_reporter_key_index]
        )
    
    after_inference_processing_time = datetime.now()
    print(f"{(after_inference_processing_time-after_inference_time).seconds} seconds for inference processing")
    print("Completed inference")

    # GET THE MAPPED COUNT SERIES BASED ON THE INFERENCE RESULTS
    print("Prepare the processed count series ")
    # Count
    all_match_set_whitelist_reporter_counter_series_results = get_counterseries_all_results(observed_guide_reporter_umi_counts_inferred, whitelist_guide_reporter_df, contains_barcode, contains_surrogate, contains_umi)

    after_counterseries_time = datetime.now()
    print(f"{(after_counterseries_time-after_inference_processing_time).seconds} seconds for counter series generation")

    print("Preparing quality control")
    quality_control_result = perform_counts_quality_control(observed_guide_reporter_umi_counts_inferred, contains_umi, contains_surrogate, contains_barcode)
    
    after_qualitycontrol_time = datetime.now()
    print(f"{(after_qualitycontrol_time-after_counterseries_time).seconds} seconds for quality control")

    count_input = CountInput(whitelist_guide_reporter_df=whitelist_guide_reporter_df,
            contains_surrogate=contains_surrogate,
            contains_barcode=contains_barcode,
            contains_umi=contains_umi,
            protospacer_hamming_threshold_strict=protospacer_hamming_threshold_strict,
            surrogate_hamming_threshold_strict=surrogate_hamming_threshold_strict,
            barcode_hamming_threshold_strict=barcode_hamming_threshold_strict)
    
    return WhitelistReporterCountsResult(all_match_set_whitelist_reporter_counter_series_results=all_match_set_whitelist_reporter_counter_series_results, observed_guide_reporter_umi_counts_inferred=observed_guide_reporter_umi_counts_inferred, quality_control_result=quality_control_result, count_input=count_input)