###
### GUIDE PARSING HELPER FUNCTIONS
###
from typeguard import typechecked
from typing import Union, Optional, List, Tuple, DefaultDict
from typing import Counter as CounterType
from ..models.mapping_models import ProtospacerCounter, ProtospacerDictUMICounter, ProtospacerBarcodeCounter, ProtospacerBarcodeDictUMICounter, ProtospacerSurrogateCounter, ProtospacerSurrogateDictUMICounter, ProtospacerSurrogateBarcodeCounter, ProtospacerSurrogateBarcodeDictUMICounter
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
def get_standard_observed_sequence_counts(  fastq_r1_fn: str, 
                                            fastq_r2_fn: Optional[str], 
                                            
                                            protospacer_pattern_regex: Optional[str],
                                            surrogate_pattern_regex: Optional[str],
                                            barcode_pattern_regex: Optional[str],
                                            umi_pattern_regex: Optional[str],

                                            protospacer_left_flank:Optional[str],
                                            protospacer_right_flank:Optional[str],
                                            protospacer_start_position:Optional[int],
                                            protospacer_end_position:Optional[int],
                                            protospacer_length: Optional[int],

                                            surrogate_left_flank:Optional[str],
                                            surrogate_right_flank:Optional[str],
                                            surrogate_start_position:Optional[int],
                                            surrogate_end_position:Optional[int],
                                            surrogate_length: Optional[int],

                                            barcode_left_flank:Optional[str],
                                            barcode_right_flank:Optional[str],
                                            barcode_start_position:Optional[int],
                                            barcode_end_position:Optional[int],
                                            barcode_length: Optional[int],

                                            umi_left_flank:Optional[str],
                                            umi_right_flank:Optional[str],
                                            umi_start_position:Optional[int],
                                            umi_end_position:Optional[int],
                                            umi_length: Optional[int],

                                            is_protospacer_r1: Optional[bool], 
                                            is_surrogate_r1: Optional[bool], 
                                            is_barcode_r1: Optional[bool],
                                            is_umi_r1: Optional[bool],
                                            is_protospacer_header: Optional[bool], 
                                            is_surrogate_header: Optional[bool], 
                                            is_barcode_header: Optional[bool],
                                            is_umi_header: Optional[bool],
                                        
                                            revcomp_protospacer: Optional[bool], 
                                            revcomp_surrogate: Optional[bool], 
                                            revcomp_barcode: Optional[bool],
                                            revcomp_umi: Optional[bool]) -> Union[CounterType[str], 
                                                                        DefaultDict[str, CounterType[str]], 
                                                                        CounterType[Tuple[str, str]], 
                                                                        DefaultDict[Tuple[str, str], CounterType[str]], 
                                                                        CounterType[Tuple[str, str, str]], 
                                                                        DefaultDict[Tuple[str, str, str], CounterType[str]]]:
    protospacer_pattern_regex = None if ((protospacer_pattern_regex is not None) and  (protospacer_pattern_regex.strip() == "")) else protospacer_pattern_regex
    surrogate_pattern_regex = None if ((surrogate_pattern_regex is not None) and (surrogate_pattern_regex.strip() == "")) else surrogate_pattern_regex
    barcode_pattern_regex = None if ((barcode_pattern_regex is not None) and  (barcode_pattern_regex.strip() == "")) else barcode_pattern_regex
    umi_pattern_regex = None if ((umi_pattern_regex is not None) and (umi_pattern_regex.strip() == "")) else umi_pattern_regex

    contains_surrogate = (surrogate_pattern_regex is not None) or (surrogate_left_flank is not None) or (surrogate_right_flank is not None) or (surrogate_start_position is not None) or (surrogate_end_position is not None) or (surrogate_length is not None)
    contains_barcode = (barcode_pattern_regex is not None) or (barcode_left_flank is not None) or (barcode_right_flank is not None) or (barcode_start_position is not None) or (barcode_end_position is not None) or (barcode_length is not None)
    contains_umi = (umi_pattern_regex is not None) or (umi_left_flank is not None) or (umi_right_flank is not None) or (umi_start_position is not None) or (barcode_end_position is not None) or (barcode_length is not None)
    
    revcomp = lambda sequence, do_revcomp: str(Seq(sequence).reverse_complement()) if do_revcomp else sequence

    def parse_sequence(
            template_sequence: str,
            sequence_pattern_regex: Optional[str],
            sequence_left_flank:Optional[str],
            sequence_right_flank:Optional[str],
            sequence_start_position:Optional[int],
            sequence_end_position:Optional[int],
            sequence_length: Optional[int],
            revcomp_sequence: Optional[bool],
            sequence_type: str):
        
        if sequence_pattern_regex is not None:
            """
                Parse entirely by regex
            """
            sequence = re.search(sequence_pattern_regex, template_sequence).group(1).strip()  
        else:
            """
                Get sequence position start
            """
            if sequence_left_flank is not None:
                """
                    Parse by left flank
                """
                position_start = template_sequence.find(sequence_left_flank) + len(sequence_left_flank)
            elif sequence_start_position is not None:
                """
                    Use provided position
                """
                position_start = sequence_start_position
            else:
                raise Exception(f"{sequence_type}_left_flank or {sequence_type}_start_position must be provided")
            



            """
                Get sequence position end
            """
            if sequence_right_flank is not None:
                """
                    Parse by left flank
                """
                position_end = template_sequence.find(sequence_right_flank)
            elif sequence_end_position is not None:
                """
                    Use provided position
                """
                position_end = sequence_end_position
            elif sequence_length is not None:
                """
                    Use provided sequence length
                """
                position_end = position_start + sequence_length
            else:
                raise Exception(f"{sequence_type}_right_flank or {sequence_type}_end_position or {sequence_type}_length must be provided")
            
            sequence = template_sequence[position_start:position_end]
        
        # Return sequence (reverse complement if necessary)
        assert revcomp_sequence is not None, f"revcomp_{sequence_type} must be provided, does the sequence need to be reverse complemented?"
        return revcomp(sequence, revcomp_sequence)

    protospacer_parse_sequence_partial = partial(parse_sequence, 
                                                 sequence_pattern_regex=protospacer_pattern_regex,
                                                 sequence_left_flank=protospacer_left_flank,
                                                 sequence_right_flank=protospacer_right_flank,
                                                 sequence_start_position=protospacer_start_position,
                                                 sequence_end_position=protospacer_end_position,
                                                 sequence_length=protospacer_length,
                                                 revcomp_sequence=revcomp_protospacer,
                                                 sequence_type="protospacer")
    surrogate_parse_sequence_partial = partial(parse_sequence, 
                                                 sequence_pattern_regex=surrogate_pattern_regex,
                                                 sequence_left_flank=surrogate_left_flank,
                                                 sequence_right_flank=surrogate_right_flank,
                                                 sequence_start_position=surrogate_start_position,
                                                 sequence_end_position=surrogate_end_position,
                                                 sequence_length=surrogate_length,
                                                 revcomp_sequence=revcomp_surrogate,
                                                 sequence_type="surrogate")
    barcode_parse_sequence_partial = partial(parse_sequence, 
                                                 sequence_pattern_regex=barcode_pattern_regex,
                                                 sequence_left_flank=barcode_left_flank,
                                                 sequence_right_flank=barcode_right_flank,
                                                 sequence_start_position=barcode_start_position,
                                                 sequence_end_position=barcode_end_position,
                                                 sequence_length=barcode_length,
                                                 revcomp_sequence=revcomp_barcode,
                                                 sequence_type="barcode")
    umi_parse_sequence_partial = partial(parse_sequence, 
                                                 sequence_pattern_regex=umi_pattern_regex,
                                                 sequence_left_flank=umi_left_flank,
                                                 sequence_right_flank=umi_right_flank,
                                                 sequence_start_position=umi_start_position,
                                                 sequence_end_position=umi_end_position,
                                                 sequence_length=umi_length,
                                                 revcomp_sequence=revcomp_umi,
                                                 sequence_type="umi")
    
    def get_sequence_from_read_group(fastq_read_group, parse_sequence_partial, is_r1:bool, is_header: bool, includes_r2: bool, sequence_type: str):
                if includes_r2:
                    r1_header_read, _ = fastq_read_group[0]
                    r1_sequence_read, r2_sequence_read = fastq_read_group[1]
                else:
                    r1_header_read = fastq_read_group[0]
                    r1_sequence_read = fastq_read_group[1]


                if is_r1 is True:
                    sequence = parse_sequence_partial(template_sequence=r1_sequence_read)
                else:
                    if is_header:
                        sequence = parse_sequence_partial(template_sequence=r1_header_read)
                    else:
                        if includes_r2:
                            sequence = parse_sequence_partial(template_sequence=r2_sequence_read)
                        else:
                            raise Exception(f"{sequence_type} cannot be parsed from R1 since is_{sequence_type}_r1 is False, is_{sequence_type}_header is False, and R2 is not provided")

                return sequence
    get_protospacer_from_read_group_partial = partial(get_sequence_from_read_group,
                                                      parse_sequence_partial=protospacer_parse_sequence_partial, 
                                                      is_r1=is_protospacer_r1, 
                                                      is_header=is_protospacer_header,
                                                      sequence_type="protospacer")
    get_surrogate_from_read_group_partial = partial(get_sequence_from_read_group,
                                                      parse_sequence_partial=surrogate_parse_sequence_partial, 
                                                      is_r1=is_surrogate_r1, 
                                                      is_header=is_surrogate_header,
                                                      sequence_type="surrogate")
    get_barcode_from_read_group_partial = partial(get_sequence_from_read_group,
                                                      parse_sequence_partial=barcode_parse_sequence_partial, 
                                                      is_r1=is_barcode_r1, 
                                                      is_header=is_barcode_header,
                                                      sequence_type="barcode")
    get_umi_from_read_group_partial = partial(get_sequence_from_read_group,
                                                      parse_sequence_partial=umi_parse_sequence_partial, 
                                                      is_r1=is_umi_r1, 
                                                      is_header=is_umi_header,
                                                      sequence_type="umi")
    
    def parse_fastq(fastq_r1_filehandler, fastq_r2_filehandler):
        if fastq_r2_fn is None: # ONLY R1
            fastq_single_read_grouper = grouper(fastq_r1_filehandler)
            if not contains_barcode: # ONLY R1; NO BARCODE
                if not contains_umi: # ONLY R1; NO BARCODE; NO UMI
                    if not contains_surrogate:
                        sequence_counter: ProtospacerCounter = Counter()
                        for fastq_single_read_group in fastq_single_read_grouper:
                            protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                            sequence_counter[protospacer] += 1
                    else:
                        # LEFTOFF HERE: Implement SURROGATE check then go to R2

                else: # ONLY R1; NO BARCODE; YES UMI
                    sequence_counter: ProtospacerDictUMICounter = defaultdict(Counter)
                    for fastq_single_read_group in fastq_single_read_grouper:
                        protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                        umi = get_umi_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                        sequence_counter[protospacer][umi] += 1
            else: # ONLY R1; YES BARCODE
                if not contains_umi: # ONLY R1; YES BARCODE; NO UMI
                    sequence_counter: ProtospacerBarcodeCounter = Counter()
                    for fastq_single_read_group in fastq_single_read_grouper:
                        protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                        barcode = get_barcode_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                        sequence_counter[(protospacer, barcode)] += 1
                else: # ONLY R1; YES BARCODE; YES UMI
                    sequence_counter: ProtospacerBarcodeDictUMICounter = defaultdict(Counter)
                    for fastq_single_read_group in fastq_single_read_grouper:
                        protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                        barcode = get_barcode_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                        umi = get_umi_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                        sequence_counter[(barcode, protospacer)][umi] += 1
        else: # YES R2
            fastq_paired_read_grouper = grouper(zip(r1_protospacer_fastq_filehandler, r2_surrogate_fastq_filehandler))
            if not contains_barcode: # YES R2; NO BARCODE
                if not contains_umi: # YES R2; NO BARCODE; NO UMI
                    sequence_counter: ProtospacerSurrogateCounter = Counter()
                    for fastq_paired_read_group in fastq_paired_read_grouper:
                        protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=True)
                        if contains_surrogate:
                            sequence_counter[(protospacer, surrogate)] += 1
                else: # YES R2; NO BARCODE; YES UMI
                    sequence_counter: ProtospacerSurrogateDictUMICounter = defaultdict(Counter)
                    for fastq_paired_read_group in fastq_paired_read_grouper:
                        r1_header_read, _ = fastq_paired_read_group[0]
                        r1_sequence_read, r2_sequence_read = fastq_paired_read_group[1]
                        
                        protospacer=revcomp(r1_sequence_read.strip(), revcomp_protospacer)
                        surrogate=revcomp(r2_sequence_read.strip(), revcomp_surrogate)
                        umi = re.search(umi_pattern_regex, r1_header_read).group(1).strip()
                        
                        sequence_counter[(protospacer, surrogate)][umi] += 1
            else: # YES R2; YES BARCODE
                if umi_pattern_regex is None: # YES R2; YES BARCODE; NO UMI
                    sequence_counter: ProtospacerSurrogateBarcodeCounter = Counter()
                    for fastq_paired_read_group in fastq_paired_read_grouper:
                        r1_header_read, _ = fastq_paired_read_group[0]
                        r1_sequence_read, r2_sequence_read = fastq_paired_read_group[1]
                        
                        protospacer=revcomp(r1_sequence_read.strip(), revcomp_protospacer)
                        surrogate=revcomp(r2_sequence_read.strip(), revcomp_surrogate)
                        barcode = revcomp(re.search(barcode_pattern_regex, r1_header_read).group(1).strip(), revcomp_barcode)
                        sequence_counter[(protospacer, surrogate, barcode)] += 1
                else: # YES R2; YES BARCODE; YES UMI
                    sequence_counter: ProtospacerSurrogateBarcodeDictUMICounter = defaultdict(Counter)
                    for fastq_paired_read_group in fastq_paired_read_grouper:
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