###
### GUIDE PARSING HELPER FUNCTIONS
###
from typeguard import typechecked
from typing import Union, Optional, List, Tuple, DefaultDict
from typing import Counter as CounterType
from ..models.types import ProtospacerCounter, ProtospacerDictUMICounter, ProtospacerBarcodeCounter, ProtospacerBarcodeDictUMICounter, ProtospacerSurrogateCounter, ProtospacerSurrogateDictUMICounter, ProtospacerSurrogateBarcodeCounter, ProtospacerSurrogateBarcodeDictUMICounter
from collections import Counter
from typing import Callable
from functools import partial
from datetime import datetime
import gzip
import re
from Bio.Seq import Seq
from collections import Counter, defaultdict

# This is for grouping the FASTQ lines in 1 tuple.
def grouper(iterable, n=4):
    "s -> (s0,s1,...sn-1), (sn,sn+1,...s2n-1), (s2n,s2n+1,...s3n-1), ..."
    return zip(*[iter(iterable)]*n)

# NOTE: Unsure if this function is used anymore
def perform_sequence_count(sequence_counter, r1_protospacer_fastq_filehandler, r2_surrogate_fastq_filehandler, barcode_pattern_regex: Optional[str] = None, umi_pattern_regex: Optional[str] = None):
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
                barcode_match = re.search(barcode_pattern_regex, r1_header_read)
                barcode = barcode_match.group(1)
                if umi_pattern_regex:
                    umi_match = re.search(umi_pattern_regex, r1_header_read)        
                    umi = umi_match.group(1)
                else:
                    pass
            else:
                if umi_pattern_regex:
                    umi = re.search(umi_pattern_regex, r1_header_read)        
                else:
                    pass
        else:
            r1_header_read, _ = fastq_paired_read_group[0]
            r1_sequence_read, r2_sequence_read = fastq_paired_read_group[1]
            surrogate=r1_sequence_read
        
        

        if barcode_pattern_regex:
            barcode = re.search(barcode_pattern_regex, r1_header_read)
        if umi_pattern_regex:
            umi = re.search(umi_pattern_regex, r1_header_read)
        

        sequence_counter[(protospacer, surrogate, barcode)][umi] += 1

    return sequence_counter

@typechecked
def get_umitools_observed_sequence_counts(r1_protospacer_fastq_file: str, r2_surrogate_fastq_file: Optional[str] = None, barcode_pattern_regex: Optional[str] = None, umi_pattern_regex: Optional[str] = None, revcomp_protospacer: bool = False, revcomp_surrogate: bool = True, revcomp_barcode: bool = True) -> Union[CounterType[str], DefaultDict[str, CounterType[str]], CounterType[Tuple[str, str]], DefaultDict[Tuple[str, str], CounterType[str]], CounterType[Tuple[str, str, str]], DefaultDict[Tuple[str, str, str], CounterType[str]]]:
    barcode_pattern_regex = None if ((barcode_pattern_regex is not None) and  (barcode_pattern_regex.strip() == "")) else barcode_pattern_regex
    umi_pattern_regex = None if ((umi_pattern_regex is not None) and (umi_pattern_regex.strip() == "")) else umi_pattern_regex

    print(f"Barcode pattern: {barcode_pattern_regex}")
    print(f"UMI pattern: {umi_pattern_regex}")

    revcomp = lambda sequence, do_revcomp: str(Seq(sequence).reverse_complement()) if do_revcomp else sequence

    def parse_fastq(r1_protospacer_fastq_filehandler, r2_surrogate_fastq_filehandler):
        if r2_surrogate_fastq_file is None: # ONLY R1
            fastq_single_read_group = grouper(r1_protospacer_fastq_filehandler)
            if barcode_pattern_regex is None: # ONLY R1; NO BARCODE
                if umi_pattern_regex is None: # ONLY R1; NO BARCODE; NO UMI
                    sequence_counter: ProtospacerCounter = Counter()
                    for fastq_paired_read_group in fastq_single_read_group:
                        r1_sequence_read = fastq_paired_read_group[1]
                        protospacer=revcomp(r1_sequence_read.strip(), revcomp_protospacer)
                        sequence_counter[protospacer] += 1
                else: # ONLY R1; NO BARCODE; YES UMI
                    sequence_counter: ProtospacerDictUMICounter = defaultdict(Counter)
                    for fastq_paired_read_group in fastq_single_read_group:
                        r1_header_read = fastq_paired_read_group[0]
                        r1_sequence_read = fastq_paired_read_group[1]
                        protospacer=revcomp(r1_sequence_read.strip(), revcomp_protospacer)
                        umi = re.search(umi_pattern_regex, r1_header_read).group(1).strip()
                        sequence_counter[protospacer][umi] += 1
            else: # ONLY R1; YES BARCODE
                if umi_pattern_regex is None: # ONLY R1; YES BARCODE; NO UMI
                    sequence_counter: ProtospacerBarcodeCounter = Counter()
                    for fastq_paired_read_group in fastq_single_read_group:
                        r1_sequence_read = fastq_paired_read_group[1]
                        protospacer=revcomp(r1_sequence_read.strip(), revcomp_protospacer)
                        barcode = revcomp(re.search(barcode_pattern_regex, r1_header_read).group(1).strip(), revcomp_barcode)
                        sequence_counter[(protospacer, barcode)] += 1
                else: # ONLY R1; YES BARCODE; YES UMI
                    sequence_counter: ProtospacerBarcodeDictUMICounter = defaultdict(Counter)
                    for fastq_paired_read_group in fastq_single_read_group:
                        r1_header_read = fastq_paired_read_group[0]
                        r1_sequence_read = fastq_paired_read_group[1]
                        protospacer=revcomp(r1_sequence_read.strip(), revcomp_protospacer)
                        barcode = revcomp(re.search(barcode_pattern_regex, r1_header_read).group(1).strip(), revcomp_barcode)
                        umi = re.search(umi_pattern_regex, r1_header_read).group(1).strip()
                        sequence_counter[(barcode, protospacer)][umi] += 1
        else: # YES R2
            fastq_paired_read_group = grouper(zip(r1_protospacer_fastq_filehandler, r2_surrogate_fastq_filehandler))
            if barcode_pattern_regex is None: # YES R2; NO BARCODE
                if umi_pattern_regex is None: # YES R2; NO BARCODE; NO UMI
                    sequence_counter: ProtospacerSurrogateCounter = Counter()
                    for fastq_paired_read_group in fastq_paired_read_group:
                        r1_sequence_read, r2_sequence_read = fastq_paired_read_group[1]
                        protospacer=revcomp(r1_sequence_read.strip(), revcomp_protospacer)
                        surrogate=revcomp(r2_sequence_read.strip(), revcomp_surrogate)
                        sequence_counter[(protospacer, surrogate)] += 1
                else: # YES R2; NO BARCODE; YES UMI
                    sequence_counter: ProtospacerSurrogateDictUMICounter = defaultdict(Counter)
                    for fastq_paired_read_group in fastq_paired_read_group:
                        r1_header_read, _ = fastq_paired_read_group[0]
                        r1_sequence_read, r2_sequence_read = fastq_paired_read_group[1]
                        
                        protospacer=revcomp(r1_sequence_read.strip(), revcomp_protospacer)
                        surrogate=revcomp(r2_sequence_read.strip(), revcomp_surrogate)
                        umi = re.search(umi_pattern_regex, r1_header_read).group(1).strip()
                        
                        sequence_counter[(protospacer, surrogate)][umi] += 1
            else: # YES R2; YES BARCODE
                if umi_pattern_regex is None: # YES R2; YES BARCODE; NO UMI
                    sequence_counter: ProtospacerSurrogateBarcodeCounter = Counter()
                    for fastq_paired_read_group in fastq_paired_read_group:
                        r1_header_read, _ = fastq_paired_read_group[0]
                        r1_sequence_read, r2_sequence_read = fastq_paired_read_group[1]
                        
                        protospacer=revcomp(r1_sequence_read.strip(), revcomp_protospacer)
                        surrogate=revcomp(r2_sequence_read.strip(), revcomp_surrogate)
                        barcode = revcomp(re.search(barcode_pattern_regex, r1_header_read).group(1).strip(), revcomp_barcode)
                        sequence_counter[(protospacer, surrogate, barcode)] += 1
                else: # YES R2; YES BARCODE; YES UMI
                    sequence_counter: ProtospacerSurrogateBarcodeDictUMICounter = defaultdict(Counter)
                    for fastq_paired_read_group in fastq_paired_read_group:
                        r1_header_read, _ = fastq_paired_read_group[0]
                        r1_sequence_read, r2_sequence_read = fastq_paired_read_group[1]
                        
                        protospacer=revcomp(r1_sequence_read.strip(), revcomp_protospacer)
                        surrogate=revcomp(r2_sequence_read.strip(), revcomp_surrogate)
                        barcode = revcomp(re.search(barcode_pattern_regex, r1_header_read).group(1).strip(), revcomp_barcode)
                        umi = re.search(umi_pattern_regex, r1_header_read).group(1).strip()
                        
                        sequence_counter[(protospacer, surrogate, barcode)][umi] += 1
        return sequence_counter
    
    before_file_loading_time = datetime.now()
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
    after_file_loading_time = datetime.now()
    print(f"{(after_file_loading_time-before_file_loading_time).seconds} seconds for file loading")

    sequence_counter = parse_fastq(r1_protospacer_fastq_filehandler, r2_surrogate_fastq_filehandler)
    after_parsing_time = datetime.now()
    print(f"{(after_parsing_time-after_file_loading_time).seconds} seconds for parsing")

    # Close the file handlers when done
    if r1_protospacer_fastq_filehandler is not None:
        r1_protospacer_fastq_filehandler.close()
    if r2_surrogate_fastq_filehandler is not None:
        r2_surrogate_fastq_filehandler.close()

    return sequence_counter