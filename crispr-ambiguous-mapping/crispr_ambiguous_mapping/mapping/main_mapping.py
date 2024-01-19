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
from ..parsing import reporter_umitools_fastq_parsing


@typechecked
def get_whitelist_reporter_counts_from_umitools_output(whitelist_guide_reporter_df: pd.DataFrame, fastq_r1_fn: str, fastq_r2_fn: str, barcode_pattern_regex: Optional[str] = None, umi_pattern_regex: Optional[str] = None, revcomp_protospacer: bool = False, revcomp_surrogate: bool = True, revcomp_barcode: bool = True, surrogate_hamming_threshold_strict: Optional[int] = 10, barcode_hamming_threshold_strict: Optional[int] = 2, protospacer_hamming_threshold_strict: Optional[int] = 7, cores: int=1):
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



