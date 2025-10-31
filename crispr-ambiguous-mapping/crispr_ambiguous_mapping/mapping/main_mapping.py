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
from dataclasses import dataclass

from ..processing import crispr_sequence_encoding
from ..processing import crispr_guide_counting
from ..processing import crispr_sequence_encoding
from ..parsing import reporter_umitools_fastq_parsing, reporter_standard_fastq_parsing
from ..models.mapping_models import WhitelistReporterCountsResult, SampleWhitelistReporterCountsResult, GeneralGuideCountType

# 
# Deprecated 
#
@typechecked
def get_whitelist_reporter_counts_from_umitools_output(whitelist_guide_reporter_df: pd.DataFrame, 
                                                       fastq_r1_fn: str, 
                                                       fastq_r2_fn: str, 
                                                       barcode_pattern_regex: Optional[str] = None, 
                                                       umi_pattern_regex: Optional[str] = None, 
                                                       revcomp_protospacer: bool = False, 
                                                       revcomp_surrogate: bool = True, 
                                                       revcomp_barcode: bool = True, 
                                                       surrogate_hamming_threshold_strict: Optional[int] = 10, 
                                                       barcode_hamming_threshold_strict: Optional[int] = 2, 
                                                       protospacer_hamming_threshold_strict: Optional[int] = 7, 
                                                       cores: int=1) -> WhitelistReporterCountsResult:
    #
    # Get counts of observed FASTQ sequences
    #
    observed_guide_reporter_umi_counts = reporter_umitools_fastq_parsing.get_umitools_observed_sequence_counts(r1_protospacer_fastq_file=fastq_r1_fn, r2_surrogate_fastq_file=fastq_r2_fn, barcode_pattern_regex=barcode_pattern_regex, umi_pattern_regex=umi_pattern_regex, revcomp_protospacer = revcomp_protospacer, revcomp_surrogate = revcomp_surrogate, revcomp_barcode = revcomp_barcode)

    #   
    # Map the observed sequences to the true whitelist sequence
    #
    get_whitelist_reporter_counts_with_umi_PARTIAL = partial(
        crispr_guide_counting.get_whitelist_reporter_counts_with_umi,
        observed_guide_reporter_umi_counts=observed_guide_reporter_umi_counts,
        whitelist_guide_reporter_df=whitelist_guide_reporter_df,
        protospacer_hamming_threshold_strict=protospacer_hamming_threshold_strict,
        surrogate_hamming_threshold_strict=surrogate_hamming_threshold_strict,
        barcode_hamming_threshold_strict=barcode_hamming_threshold_strict,
        cores=cores
    )
    
    if fastq_r2_fn is None: # ONLY R1
        if barcode_pattern_regex is None: # ONLY R1; NO BARCODE
            if umi_pattern_regex is None: # ONLY R1; NO BARCODE; NO UMI
                # These lines are really added just to document the type... observed_guide_reporter_umi_counts is already added to the partial above anyways
                observed_guide_reporter_umi_counts: CounterType[str] = observed_guide_reporter_umi_counts
                return get_whitelist_reporter_counts_with_umi_PARTIAL(contains_surrogate=False, contains_barcode=False, contains_umi = False)
            else: # ONLY R1; NO BARCODE; YES UMI
                observed_guide_reporter_umi_counts: DefaultDict[str, CounterType[str]] = observed_guide_reporter_umi_counts
                return get_whitelist_reporter_counts_with_umi_PARTIAL(contains_surrogate=False, contains_barcode=False, contains_umi = True)
        else: # ONLY R1; YES BARCODE
            if umi_pattern_regex is None: # ONLY R1; YES BARCODE; NO UMI
                observed_guide_reporter_umi_counts: CounterType[Tuple[str, str]] = observed_guide_reporter_umi_counts
                return get_whitelist_reporter_counts_with_umi_PARTIAL(contains_surrogate=False, contains_barcode=True, contains_umi = False)
            else: # ONLY R1; YES BARCODE; YES UMI
                observed_guide_reporter_umi_counts: DefaultDict[Tuple[str, str], CounterType[str]] = observed_guide_reporter_umi_counts
                return get_whitelist_reporter_counts_with_umi_PARTIAL(contains_surrogate=False, contains_barcode=True, contains_umi = True)
    else: # YES R2
        if barcode_pattern_regex is None: # YES R2; NO BARCODE
            if umi_pattern_regex is None: # YES R2; NO BARCODE; NO UMI
                observed_guide_reporter_umi_counts: CounterType[Tuple[str, str]] = observed_guide_reporter_umi_counts
                return get_whitelist_reporter_counts_with_umi_PARTIAL(contains_surrogate=True, contains_barcode=False, contains_umi = False)
            else: # YES R2; NO BARCODE; YES UMI
                observed_guide_reporter_umi_counts: DefaultDict[Tuple[str, str], CounterType[str]] = observed_guide_reporter_umi_counts
                return get_whitelist_reporter_counts_with_umi_PARTIAL(contains_surrogate=True, contains_barcode=False, contains_umi = True)
        else: # YES R2; YES BARCODE
            if umi_pattern_regex is None: # YES R2; YES BARCODE; NO UMI
                observed_guide_reporter_umi_counts: CounterType[Tuple[str, str, str]] = observed_guide_reporter_umi_counts
                return get_whitelist_reporter_counts_with_umi_PARTIAL(contains_surrogate=True, contains_barcode=True, contains_umi = False)
            else: # YES R2; YES BARCODE; YES UMI
                observed_guide_reporter_umi_counts: DefaultDict[Tuple[str, str, str], CounterType[str]] = observed_guide_reporter_umi_counts
                return get_whitelist_reporter_counts_with_umi_PARTIAL(contains_surrogate=True, contains_barcode=True, contains_umi = True)


# TODO 11/20/2024: Implement UMI parsing via sequence or via R1 header
#@typechecked
def get_whitelist_reporter_counts_from_fastq(whitelist_guide_reporter_df: Optional[pd.DataFrame], 
                                                       fastq_r1_fns: List[str], 
                                                       fastq_r2_fns: Optional[List[str]] = None, 
                                                       
                                                       protospacer_pattern_regex: Optional[str] = None,
                                                       surrogate_pattern_regex: Optional[str] = None,
                                                       guide_barcode_pattern_regex: Optional[str] = None,
                                                       guide_umi_pattern_regex: Optional[str] = None,
                                                       sample_barcode_pattern_regex: Optional[str] = None,

                                                       protospacer_left_flank:Optional[str] = None,
                                                       protospacer_right_flank:Optional[str] = None,
                                                       protospacer_start_position:Optional[int] = None,
                                                       protospacer_end_position:Optional[int] = None,
                                                       protospacer_length: Optional[int] = None,

                                                       surrogate_left_flank:Optional[str] = None,
                                                       surrogate_right_flank:Optional[str] = None,
                                                       surrogate_start_position:Optional[int] = None,
                                                       surrogate_end_position:Optional[int] = None,
                                                       surrogate_length: Optional[int] = None,

                                                       guide_barcode_left_flank:Optional[str] = None,
                                                       guide_barcode_right_flank:Optional[str] = None,
                                                       guide_barcode_start_position:Optional[int] = None,
                                                       guide_barcode_end_position:Optional[int] = None,
                                                       guide_barcode_length: Optional[int] = None,

                                                       guide_umi_left_flank:Optional[str] = None,
                                                       guide_umi_right_flank:Optional[str] = None,
                                                       guide_umi_start_position:Optional[int] = None,
                                                       guide_umi_end_position:Optional[int] = None,
                                                       guide_umi_length: Optional[int] = None,

                                                       sample_barcode_left_flank:Optional[str] = None,
                                                       sample_barcode_right_flank:Optional[str] = None,
                                                       sample_barcode_start_position:Optional[int] = None,
                                                       sample_barcode_end_position:Optional[int] = None,
                                                       sample_barcode_length: Optional[int] = None,

                                                       is_protospacer_r1: Optional[bool] = None, 
                                                       is_surrogate_r1: Optional[bool] = None, 
                                                       is_guide_barcode_r1: Optional[bool] = None,
                                                       is_guide_umi_r1: Optional[bool] = None,
                                                       is_sample_barcode_r1: Optional[bool] = None,

                                                       is_protospacer_header: Optional[bool] = None, 
                                                       is_surrogate_header: Optional[bool] = None, 
                                                       is_guide_barcode_header: Optional[bool] = None,
                                                       is_guide_umi_header: Optional[bool] = None,
                                                       is_sample_barcode_header: Optional[bool] = None,
                                                    
                                                       revcomp_protospacer: Optional[bool] = None, 
                                                       revcomp_surrogate: Optional[bool] = None, 
                                                       revcomp_guide_barcode: Optional[bool] = None, 
                                                       revcomp_guide_umi: Optional[bool] = None,
                                                       revcomp_sample_barcode: Optional[bool] = None, 
                                                        
                                                       surrogate_hamming_threshold_strict: Optional[int] = None, 
                                                       guide_barcode_hamming_threshold_strict: Optional[int] = None, 
                                                       protospacer_hamming_threshold_strict: Optional[int] = None, 

                                                       store_intermediates: bool = False,
                                                       cores: int=1) -> Union[WhitelistReporterCountsResult, SampleWhitelistReporterCountsResult]:
    # Input parameter validation checks

    protospacer_pattern_regex = None if ((protospacer_pattern_regex is not None) and  (protospacer_pattern_regex.strip() == "")) else protospacer_pattern_regex
    surrogate_pattern_regex = None if ((surrogate_pattern_regex is not None) and (surrogate_pattern_regex.strip() == "")) else surrogate_pattern_regex
    guide_barcode_pattern_regex = None if ((guide_barcode_pattern_regex is not None) and  (guide_barcode_pattern_regex.strip() == "")) else guide_barcode_pattern_regex
    guide_umi_pattern_regex = None if ((guide_umi_pattern_regex is not None) and (guide_umi_pattern_regex.strip() == "")) else guide_umi_pattern_regex

    contains_surrogate = (surrogate_pattern_regex is not None) or (surrogate_left_flank is not None) or (surrogate_right_flank is not None) or (surrogate_start_position is not None) or (surrogate_end_position is not None) or (surrogate_length is not None)
    contains_guide_barcode = (guide_barcode_pattern_regex is not None) or (guide_barcode_left_flank is not None) or (guide_barcode_right_flank is not None) or (guide_barcode_start_position is not None) or (guide_barcode_end_position is not None) or (guide_barcode_length is not None)
    contains_guide_umi = (guide_umi_pattern_regex is not None) or (guide_umi_left_flank is not None) or (guide_umi_right_flank is not None) or (guide_umi_start_position is not None) or (guide_barcode_end_position is not None) or (guide_barcode_length is not None)
    contains_sample_barcode = (sample_barcode_pattern_regex is not None) or (sample_barcode_left_flank is not None) or (sample_barcode_right_flank is not None) or (sample_barcode_start_position is not None) or (sample_barcode_end_position is not None) or (sample_barcode_length is not None)
    
    def preprocess_sequence(sequence):
        if sequence is not None:
            sequence = sequence.upper().rstrip()
        return sequence
    
    protospacer_left_flank = preprocess_sequence(protospacer_left_flank)
    protospacer_right_flank = preprocess_sequence(protospacer_right_flank)
    surrogate_left_flank = preprocess_sequence(surrogate_left_flank)
    surrogate_right_flank = preprocess_sequence(surrogate_right_flank)
    guide_barcode_left_flank = preprocess_sequence(guide_barcode_left_flank)
    guide_barcode_right_flank = preprocess_sequence(guide_barcode_right_flank)
    guide_umi_left_flank = preprocess_sequence(guide_umi_left_flank)
    guide_umi_right_flank = preprocess_sequence(guide_umi_right_flank)
    sample_barcode_left_flank = preprocess_sequence(sample_barcode_left_flank)
    sample_barcode_right_flank = preprocess_sequence(sample_barcode_right_flank)

    print(f"Contains surrogate: {contains_surrogate}")
    print(f"Contains guide barcode: {contains_guide_barcode}")
    print(f"Contains guide UMI: {contains_guide_umi}")
    print(f"Contains sample barcode: {contains_sample_barcode}")
    #
    # Get counts of observed FASTQ sequences
    #
    observed_guide_reporter_umi_counts: GeneralGuideCountType = reporter_standard_fastq_parsing.get_standard_observed_sequence_counts(
                                            fastq_r1_fns=fastq_r1_fns, 
                                            fastq_r2_fns=fastq_r2_fns, 
                                            
                                            protospacer_pattern_regex=protospacer_pattern_regex,
                                            surrogate_pattern_regex=surrogate_pattern_regex,
                                            guide_barcode_pattern_regex=guide_barcode_pattern_regex,
                                            guide_umi_pattern_regex=guide_umi_pattern_regex,
                                            sample_barcode_pattern_regex=sample_barcode_pattern_regex,

                                            protospacer_left_flank=protospacer_left_flank,
                                            protospacer_right_flank=protospacer_right_flank,
                                            protospacer_start_position=protospacer_start_position,
                                            protospacer_end_position=protospacer_end_position,
                                            protospacer_length=protospacer_length,

                                            surrogate_left_flank=surrogate_left_flank,
                                            surrogate_right_flank=surrogate_right_flank,
                                            surrogate_start_position=surrogate_start_position,
                                            surrogate_end_position=surrogate_end_position,
                                            surrogate_length=surrogate_length,

                                            guide_barcode_left_flank=guide_barcode_left_flank,
                                            guide_barcode_right_flank=guide_barcode_right_flank,
                                            guide_barcode_start_position=guide_barcode_start_position,
                                            guide_barcode_end_position=guide_barcode_end_position,
                                            guide_barcode_length=guide_barcode_length,

                                            guide_umi_left_flank=guide_umi_left_flank,
                                            guide_umi_right_flank=guide_umi_right_flank,
                                            guide_umi_start_position=guide_umi_start_position,
                                            guide_umi_end_position=guide_umi_end_position,
                                            guide_umi_length=guide_umi_length,

                                            sample_barcode_left_flank=sample_barcode_left_flank,
                                            sample_barcode_right_flank=sample_barcode_right_flank,
                                            sample_barcode_start_position=sample_barcode_start_position,
                                            sample_barcode_end_position=sample_barcode_end_position,
                                            sample_barcode_length=sample_barcode_length,

                                            is_protospacer_r1=is_protospacer_r1, 
                                            is_surrogate_r1=is_surrogate_r1, 
                                            is_guide_barcode_r1=is_guide_barcode_r1,
                                            is_guide_umi_r1=is_guide_umi_r1,
                                            is_sample_barcode_r1=is_sample_barcode_r1,
                                            
                                            is_protospacer_header=is_protospacer_header, 
                                            is_surrogate_header=is_surrogate_header, 
                                            is_guide_barcode_header=is_guide_barcode_header,
                                            is_guide_umi_header=is_guide_umi_header,
                                            is_sample_barcode_header=is_sample_barcode_header,
                                        
                                            revcomp_protospacer=revcomp_protospacer, 
                                            revcomp_surrogate=revcomp_surrogate, 
                                            revcomp_guide_barcode=revcomp_guide_barcode,
                                            revcomp_guide_umi=revcomp_guide_umi,
                                            revcomp_sample_barcode=revcomp_sample_barcode,
                                            
                                            contains_guide_surrogate=contains_surrogate,
                                            contains_guide_barcode=contains_guide_barcode,
                                            contains_guide_umi=contains_guide_umi,
                                            contains_sample_barcode=contains_sample_barcode)

    print(f"Number of unique observed parsed sequences: {len(observed_guide_reporter_umi_counts.keys())}")

    #   
    # Map the observed sequences to the true whitelist sequence
    #
    return crispr_guide_counting.get_whitelist_reporter_counts_with_umi(
        observed_guide_reporter_umi_counts=observed_guide_reporter_umi_counts,
        whitelist_guide_reporter_df=whitelist_guide_reporter_df,
        contains_guide_surrogate=contains_surrogate, 
        contains_guide_barcode=contains_guide_barcode, 
        contains_guide_umi = contains_guide_umi,
        contains_sample_barcode=contains_sample_barcode, 
        protospacer_hamming_threshold_strict=protospacer_hamming_threshold_strict,
        surrogate_hamming_threshold_strict=surrogate_hamming_threshold_strict,
        guide_barcode_hamming_threshold_strict=guide_barcode_hamming_threshold_strict,
        store_intermediates=store_intermediates,
        cores=cores
    )

                

