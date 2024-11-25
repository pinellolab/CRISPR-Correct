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
from ..models.mapping_models import WhitelistReporterCountsResult


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
def get_whitelist_reporter_counts_from_fastq(whitelist_guide_reporter_df: pd.DataFrame, 
                                                       fastq_r1_fn: str, 
                                                       fastq_r2_fn: Optional[str] = None, 
                                                       
                                                       protospacer_pattern_regex: Optional[str] = None,
                                                       surrogate_pattern_regex: Optional[str] = None,
                                                       barcode_pattern_regex: Optional[str] = None,
                                                       umi_pattern_regex: Optional[str] = None,

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

                                                       barcode_left_flank:Optional[str] = None,
                                                       barcode_right_flank:Optional[str] = None,
                                                       barcode_start_position:Optional[int] = None,
                                                       barcode_end_position:Optional[int] = None,
                                                       barcode_length: Optional[int] = None,

                                                       umi_left_flank:Optional[str] = None,
                                                       umi_right_flank:Optional[str] = None,
                                                       umi_start_position:Optional[int] = None,
                                                       umi_end_position:Optional[int] = None,
                                                       umi_length: Optional[int] = None,

                                                       is_protospacer_r1: Optional[bool] = None, 
                                                       is_surrogate_r1: Optional[bool] = None, 
                                                       is_barcode_r1: Optional[bool] = None,
                                                       is_umi_r1: Optional[bool] = None,
                                                       is_protospacer_header: Optional[bool] = None, 
                                                       is_surrogate_header: Optional[bool] = None, 
                                                       is_barcode_header: Optional[bool] = None,
                                                       is_umi_header: Optional[bool] = None,
                                                    
                                                       revcomp_protospacer: Optional[bool] = None, 
                                                       revcomp_surrogate: Optional[bool] = None, 
                                                       revcomp_barcode: Optional[bool] = None, 
                                                       revcomp_umi: Optional[bool] = None,
                                                        
                                                       surrogate_hamming_threshold_strict: Optional[int] = None, 
                                                       barcode_hamming_threshold_strict: Optional[int] = None, 
                                                       protospacer_hamming_threshold_strict: Optional[int] = None, 
                                                       cores: int=1) -> WhitelistReporterCountsResult:
    # Input parameter validation checks

    protospacer_pattern_regex = None if ((protospacer_pattern_regex is not None) and  (protospacer_pattern_regex.strip() == "")) else protospacer_pattern_regex
    surrogate_pattern_regex = None if ((surrogate_pattern_regex is not None) and (surrogate_pattern_regex.strip() == "")) else surrogate_pattern_regex
    barcode_pattern_regex = None if ((barcode_pattern_regex is not None) and  (barcode_pattern_regex.strip() == "")) else barcode_pattern_regex
    umi_pattern_regex = None if ((umi_pattern_regex is not None) and (umi_pattern_regex.strip() == "")) else umi_pattern_regex

    contains_surrogate = (surrogate_pattern_regex is not None) or (surrogate_left_flank is not None) or (surrogate_right_flank is not None) or (surrogate_start_position is not None) or (surrogate_end_position is not None) or (surrogate_length is not None)
    contains_barcode = (barcode_pattern_regex is not None) or (barcode_left_flank is not None) or (barcode_right_flank is not None) or (barcode_start_position is not None) or (barcode_end_position is not None) or (barcode_length is not None)
    contains_umi = (umi_pattern_regex is not None) or (umi_left_flank is not None) or (umi_right_flank is not None) or (umi_start_position is not None) or (barcode_end_position is not None) or (barcode_length is not None)
    

    #
    # Get counts of observed FASTQ sequences
    #
    observed_guide_reporter_umi_counts = reporter_standard_fastq_parsing.get_standard_observed_sequence_counts(fastq_r1_fn=fastq_r1_fn, 
                                            fastq_r2_fn=fastq_r2_fn, 
                                            
                                            protospacer_pattern_regex=protospacer_pattern_regex,
                                            surrogate_pattern_regex=surrogate_pattern_regex,
                                            barcode_pattern_regex=barcode_pattern_regex,
                                            umi_pattern_regex=umi_pattern_regex,

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

                                            barcode_left_flank=barcode_left_flank,
                                            barcode_right_flank=barcode_right_flank,
                                            barcode_start_position=barcode_start_position,
                                            barcode_end_position=barcode_end_position,
                                            barcode_length=barcode_length,

                                            umi_left_flank=umi_left_flank,
                                            umi_right_flank=umi_right_flank,
                                            umi_start_position=umi_start_position,
                                            umi_end_position=umi_end_position,
                                            umi_length=umi_length,

                                            is_protospacer_r1=is_protospacer_r1, 
                                            is_surrogate_r1=is_surrogate_r1, 
                                            is_barcode_r1=is_barcode_r1,
                                            is_umi_r1=is_umi_r1,
                                            is_protospacer_header=is_protospacer_header, 
                                            is_surrogate_header=is_surrogate_header, 
                                            is_barcode_header=is_barcode_header,
                                            is_umi_header=is_umi_header,
                                        
                                            revcomp_protospacer=revcomp_protospacer, 
                                            revcomp_surrogate=revcomp_surrogate, 
                                            revcomp_barcode=revcomp_barcode,
                                            revcomp_umi=revcomp_umi,
                                            
                                            contains_surrogate=contains_surrogate,
                                            contains_barcode=contains_barcode,
                                            contains_umi=contains_umi)

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

    return get_whitelist_reporter_counts_with_umi_PARTIAL(contains_surrogate=contains_surrogate, contains_barcode=contains_barcode, contains_umi = contains_umi)
                

