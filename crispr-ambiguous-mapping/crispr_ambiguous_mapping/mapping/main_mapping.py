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
from dataclasses import dataclass

from . import sequence_encoding

from . import guide_inference
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
    
    return guide_inference.get_whitelist_guide_counts(observed_guide_sequences_counts, whitelist_guide_sequences_series, hamming_threshold_strict, hamming_threshold_dynamic, cores)

'''
    Take in input FASTQ filename, and a set of whitelisted guide sequences
'''
@typechecked
def get_whitelist_reporter_counts_from_reporter_tsv(whitelist_guide_reporter_df: pd.DataFrame, reporter_tsv_fn: str, surrogate_hamming_threshold_strict: int = 10, barcode_hamming_threshold_strict: int = 2, hamming_threshold_strict: int = 7, hamming_threshold_dynamic: bool = False, cores: int=1):
    observed_guide_reporter_counts: Union[Counter[Tuple[str,str,str]], Counter[str]] = reporter_tsv_parsing.get_reporter_tsv_observed_sequence_counts(reporter_tsv_fn, include_surrogate = True, cores=cores)

    return guide_inference.get_whitelist_reporter_counts(observed_guide_reporters_counts=observed_guide_reporter_counts, whitelist_guide_reporter_df=whitelist_guide_reporter_df, surrogate_hamming_threshold_strict=surrogate_hamming_threshold_strict, barcode_hamming_threshold_strict=barcode_hamming_threshold_strict, hamming_threshold_strict=hamming_threshold_strict, hamming_threshold_dynamic=hamming_threshold_dynamic, cores=cores)


@typechecked
def get_whitelist_reporter_counts_from_umitools_output(whitelist_guide_reporter_df: pd.DataFrame, fastq_r1_fn: str, fastq_r2_fn: str, barcode_pattern_regex: Optional[str] = None, umi_pattern_regex: Optional[str] = None, surrogate_hamming_threshold_strict: Optional[int] = 10, barcode_hamming_threshold_strict: Optional[int] = 2, protospacer_hamming_threshold_strict: Optional[int] = 7, cores: int=1):
    #
    # Get counts of observed FASTQ sequences
    #
    observed_guide_reporter_umi_counts = reporter_umitools_fastq_parsing.get_umitools_observed_sequence_counts(r1_protospacer_fastq_file=fastq_r1_fn, r2_surrogate_fastq_file=fastq_r2_fn, barcode_pattern_regex=barcode_pattern_regex, umi_pattern_regex=umi_pattern_regex)

    #   
    # Map the observed sequences to the true whitelist sequence
    #
    get_whitelist_reporter_counts_with_umi_PARTIAL = partial(
        guide_inference.get_whitelist_reporter_counts_with_umi,
        observed_guide_reporter_umi_counts=observed_guide_reporter_umi_counts,
        whitelist_guide_reporter_df=whitelist_guide_reporter_df,
        contains_surrogate=False,
        contains_barcode=False,
        contains_umi=True,
        protospacer_hamming_threshold_strict=protospacer_hamming_threshold_strict,
        surrogate_hamming_threshold_strict=surrogate_hamming_threshold_strict,
        barcode_hamming_threshold_strict=barcode_hamming_threshold_strict,
        cores=cores
    )
    if fastq_r2_fn is None: # ONLY R1
        if barcode_pattern_regex is None: # ONLY R1; NO BARCODE
            if umi_pattern_regex is None: # ONLY R1; NO BARCODE; NO UMI
                # These lines are really added just to document the type... observed_guide_reporter_umi_counts is already added to the partial above anyways
                observed_guide_reporter_umi_counts: Counter[str] = observed_guide_reporter_umi_counts
                return get_whitelist_reporter_counts_with_umi_PARTIAL(contains_surrogate=False, contains_barcode=False, contains_umi = False)
            else: # ONLY R1; NO BARCODE; YES UMI
                observed_guide_reporter_umi_counts: defaultdict(str, Counter[str]) = observed_guide_reporter_umi_counts
                return get_whitelist_reporter_counts_with_umi_PARTIAL(contains_surrogate=False, contains_barcode=False, contains_umi = True)
        else: # ONLY R1; YES BARCODE
            if umi_pattern_regex is None: # ONLY R1; YES BARCODE; NO UMI
                observed_guide_reporter_umi_counts: Counter[Tuple[str, str]] = observed_guide_reporter_umi_counts
                return get_whitelist_reporter_counts_with_umi_PARTIAL(contains_surrogate=False, contains_barcode=True, contains_umi = False)
            else: # ONLY R1; YES BARCODE; YES UMI
                observed_guide_reporter_umi_counts: defaultdict(Tuple[str, str], Counter[str]) = observed_guide_reporter_umi_counts
                return get_whitelist_reporter_counts_with_umi_PARTIAL(contains_surrogate=False, contains_barcode=True, contains_umi = True)
    else: # YES R2
        if barcode_pattern_regex is None: # YES R2; NO BARCODE
            if umi_pattern_regex is None: # YES R2; NO BARCODE; NO UMI
                observed_guide_reporter_umi_counts: Counter[Tuple[str, str]] = observed_guide_reporter_umi_counts
                return get_whitelist_reporter_counts_with_umi_PARTIAL(contains_surrogate=True, contains_barcode=False, contains_umi = False)
            else: # YES R2; NO BARCODE; YES UMI
                observed_guide_reporter_umi_counts: defaultdict(Tuple[str, str], Counter[str]) = observed_guide_reporter_umi_counts
                return get_whitelist_reporter_counts_with_umi_PARTIAL(contains_surrogate=True, contains_barcode=False, contains_umi = True)
        else: # YES R2; YES BARCODE
            if umi_pattern_regex is None: # YES R2; YES BARCODE; NO UMI
                observed_guide_reporter_umi_counts: Counter[Tuple[str, str, str]] = observed_guide_reporter_umi_counts
                return get_whitelist_reporter_counts_with_umi_PARTIAL(contains_surrogate=True, contains_barcode=True, contains_umi = False)
            else: # YES R2; YES BARCODE; YES UMI
                observed_guide_reporter_umi_counts: defaultdict(Tuple[str, str, str], Counter[str]) = observed_guide_reporter_umi_counts
                return get_whitelist_reporter_counts_with_umi_PARTIAL(contains_surrogate=True, contains_barcode=True, contains_umi = True)



