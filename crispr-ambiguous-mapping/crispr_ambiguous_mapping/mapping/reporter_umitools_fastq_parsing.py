###
### GUIDE PARSING HELPER FUNCTIONS
###
from typeguard import typechecked
from typing import Union, Optional, List, Tuple, DefaultDict
from typing import Counter as CounterType
from collections import Counter
from typing import Callable
from functools import partial
import gzip
import re
from collections import Counter, defaultdict

# This is for grouping the FASTQ lines in 1 tuple.
def grouper(iterable, n=4):
    "s -> (s0,s1,...sn-1), (sn,sn+1,...s2n-1), (s2n,s2n+1,...s3n-1), ..."
    return zip(*[iter(iterable)]*n)

def perform_sequence_count(sequence_counter, r1_protospacer_fastq_filehandler, r2_surrogate_fastq_filehandler, barcode_pattern_regex: Optional[str] = None, i7_pattern_regex: Optional[str] = None):
    if r2_surrogate_fastq_filehandler is None:
        fastq_read_grouper = grouper(r1_protospacer_fastq_filehandler)   
    else:
        fastq_read_grouper = grouper(zip(r1_protospacer_fastq_filehandler, r2_surrogate_fastq_filehandler))

    for fastq_paired_read_group in fastq_read_grouper:
        if r2_surrogate_fastq_filehandler is None:
            r1_header_read = fastq_paired_read_group[0]
            r1_sequence_read = fastq_paired_read_group[1]
            protospacer=r1_sequence_read
            if barcode_pattern_regex:
                barcode = re.search(barcode_pattern_regex, r1_header_read)
                if i7_pattern_regex:
                    i7_index = re.search(i7_pattern_regex, r1_header_read)        
                else:
                    pass
            else:
                if i7_pattern_regex:
                    i7_index = re.search(i7_pattern_regex, r1_header_read)        
                else:
                    pass
        else:
            r1_header_read, _ = fastq_paired_read_group[0]
            r1_sequence_read, r2_sequence_read = fastq_paired_read_group[1]
            surrogate=r1_sequence_read
        
        

        if barcode_pattern_regex:
            barcode = re.search(barcode_pattern_regex, r1_header_read)
        if i7_pattern_regex:
            i7_index = re.search(i7_pattern_regex, r1_header_read)
        

        sequence_counter[(protospacer, surrogate, barcode)][i7_index] += 1

    return sequence_counter

@typechecked
def get_umitools_observed_sequence_counts(r1_protospacer_fastq_file: str, r2_surrogate_fastq_file: Optional[str] = None, barcode_pattern_regex: Optional[str] = None, umi_pattern_regex: Optional[str] = None) -> Union[CounterType[str], DefaultDict[str, CounterType[str]], CounterType[Tuple[str, str]], DefaultDict[Tuple[str, str], CounterType[str]], CounterType[Tuple[str, str, str]], DefaultDict[Tuple[str, str, str], CounterType[str]]]:
    def parse_fastq(r1_protospacer_fastq_filehandler, r2_surrogate_fastq_filehandler):
        if r2_surrogate_fastq_file is None: # ONLY R1
            fastq_single_read_group = grouper(r1_protospacer_fastq_filehandler)
            if barcode_pattern_regex is None: # ONLY R1; NO BARCODE
                if umi_pattern_regex is None: # ONLY R1; NO BARCODE; NO UMI
                    sequence_counter: CounterType[str] = Counter()
                    for fastq_paired_read_group in fastq_single_read_group:
                        r1_sequence_read = fastq_paired_read_group[1]
                        protospacer=r1_sequence_read
                        sequence_counter[protospacer] += 1
                else: # ONLY R1; NO BARCODE; YES UMI
                    sequence_counter: DefaultDict[str, CounterType[str]] = defaultdict(Counter)
                    for fastq_paired_read_group in fastq_single_read_group:
                        r1_header_read = fastq_paired_read_group[0]
                        r1_sequence_read = fastq_paired_read_group[1]
                        protospacer=r1_sequence_read
                        umi = re.search(umi_pattern_regex, r1_header_read)
                        sequence_counter[protospacer][umi] += 1
            else: # ONLY R1; YES BARCODE
                if umi_pattern_regex is None: # ONLY R1; YES BARCODE; NO UMI
                    sequence_counter: CounterType[Tuple[str, str]] = Counter()
                    for fastq_paired_read_group in fastq_single_read_group:
                        r1_sequence_read = fastq_paired_read_group[1]
                        protospacer=r1_sequence_read
                        barcode = re.search(barcode_pattern_regex, r1_header_read)
                        sequence_counter[(protospacer, barcode)] += 1
                else: # ONLY R1; YES BARCODE; YES UMI
                    sequence_counter: DefaultDict[Tuple[str, str], CounterType[str]] = defaultdict(Counter)
                    for fastq_paired_read_group in fastq_single_read_group:
                        r1_header_read = fastq_paired_read_group[0]
                        r1_sequence_read = fastq_paired_read_group[1]
                        protospacer=r1_sequence_read
                        barcode = re.search(barcode_pattern_regex, r1_header_read)
                        umi = re.search(umi_pattern_regex, r1_header_read)
                        sequence_counter[(barcode, protospacer)][umi] += 1
        else: # YES R2
            fastq_paired_read_group = grouper(zip(r1_protospacer_fastq_filehandler, r2_surrogate_fastq_filehandler))
            if barcode_pattern_regex is None: # YES R2; NO BARCODE
                if umi_pattern_regex is None: # YES R2; NO BARCODE; NO UMI
                    sequence_counter: CounterType[Tuple[str, str]] = Counter()
                    for fastq_paired_read_group in fastq_paired_read_group:
                        r1_sequence_read, r2_sequence_read = fastq_paired_read_group[1]
                        protospacer=r1_sequence_read
                        surrogate=r2_sequence_read
                        sequence_counter[(protospacer, surrogate)] += 1
                else: # YES R2; NO BARCODE; YES UMI
                    sequence_counter: DefaultDict[Tuple[str,str], CounterType[str]] = defaultdict(Counter)
                    for fastq_paired_read_group in fastq_paired_read_group:
                        r1_header_read, _ = fastq_paired_read_group[0]
                        r1_sequence_read, r2_sequence_read = fastq_paired_read_group[1]
                        
                        protospacer=r1_sequence_read
                        surrogate=r2_sequence_read
                        umi = re.search(umi_pattern_regex, r1_header_read)
                        
                        sequence_counter[(protospacer, surrogate)][umi] += 1
            else: # YES R2; YES BARCODE
                if umi_pattern_regex is None: # YES R2; YES BARCODE; NO UMI
                    sequence_counter: CounterType[Tuple[str, str, str]] = Counter()
                    for fastq_paired_read_group in fastq_paired_read_group:
                        r1_header_read, _ = fastq_paired_read_group[0]
                        r1_sequence_read, r2_sequence_read = fastq_paired_read_group[1]
                        
                        protospacer=r1_sequence_read
                        surrogate=r2_sequence_read
                        barcode = re.search(barcode_pattern_regex, r1_header_read)
                        umi = re.search(umi_pattern_regex, r1_header_read)
                        sequence_counter[(protospacer, surrogate, barcode)] += 1
                else: # YES R2; YES BARCODE; YES UMI
                    sequence_counter: DefaultDict[Tuple[str, str, str], CounterType[str]] = defaultdict(Counter)
                    for fastq_paired_read_group in fastq_paired_read_group:
                        r1_header_read, _ = fastq_paired_read_group[0]
                        r1_sequence_read, r2_sequence_read = fastq_paired_read_group[1]
                        
                        protospacer=r1_sequence_read
                        surrogate=r2_sequence_read
                        barcode = re.search(barcode_pattern_regex, r1_header_read)
                        umi = re.search(umi_pattern_regex, r1_header_read)
                        
                        sequence_counter[(protospacer, surrogate, barcode)][umi] += 1
        return sequence_counter
    
    # Open R1 file handler
    r1_protospacer_fastq_filehandler = None
    if r1_protospacer_fastq_file.endswith('.gz'):
        print(f"Opening FASTQ.gz file with gzip, filename={r1_protospacer_fastq_file}")
        r1_protospacer_fastq_filehandler = gzip.open(r1_protospacer_fastq_file, "rt", encoding="utf-8")
    else:
        print(f"Opening FASTQ file, filename={r1_protospacer_fastq_file}")
        r1_protospacer_fastq_filehandler = open(r1_protospacer_fastq_file, "r")
    
    # Open R2 file handler if provided
    r2_surrogate_fastq_filehandler = None
    if r2_surrogate_fastq_file is not None:
        if r2_surrogate_fastq_file.endswith('.gz'):
            print(f"Opening FASTQ.gz file with gzip, filename={r2_surrogate_fastq_file}")
            r2_surrogate_fastq_filehandler = gzip.open(r2_surrogate_fastq_file, "rt", encoding="utf-8")
        else:
            print(f"Opening FASTQ file, filename={r2_surrogate_fastq_file}")
            r2_surrogate_fastq_filehandler = open(r2_surrogate_fastq_file, "r")

    sequence_counter = parse_fastq(r1_protospacer_fastq_filehandler, r2_surrogate_fastq_filehandler)

    # Close the file handlers when done
    if r1_protospacer_fastq_filehandler is not None:
        r1_protospacer_fastq_filehandler.close()
    if r2_surrogate_fastq_filehandler is not None:
        r2_surrogate_fastq_filehandler.close()

    return sequence_counter