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
                                            revcomp_umi: Optional[bool],
                                            contains_surrogate: bool,
                                            contains_barcode: bool,
                                            contains_umi: bool) -> Union[CounterType[str], 
                                                                        DefaultDict[str, CounterType[str]], 
                                                                        CounterType[Tuple[str, str]], 
                                                                        DefaultDict[Tuple[str, str], CounterType[str]], 
                                                                        CounterType[Tuple[str, str, str]], 
                                                                        DefaultDict[Tuple[str, str, str], CounterType[str]]]:
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
            sequence = re.search(sequence_pattern_regex, template_sequence).group(1).rstrip()  
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
            
            sequence = template_sequence[position_start:position_end].rstrip() # rstrip is important, since need to remove newlines if sequence_length exceeds template_sequence length
        
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
            print("Only R1 FASTQ is provided, R2 NOT provided")
            fastq_single_read_grouper = grouper(fastq_r1_filehandler)
            if not contains_barcode: # ONLY R1; NO BARCODE
                if not contains_umi: # ONLY R1; NO BARCODE; NO UMI
                    if not contains_surrogate:
                        print(f"Performing protospacer-only parsing using function: {get_protospacer_from_read_group_partial}")
                        sequence_counter: ProtospacerCounter = Counter()
                        for fastq_single_read_group in fastq_single_read_grouper:
                            protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                            sequence_counter[(protospacer,)] += 1
                    else:
                        print(f"Performing protospacer+surrrogate parsing using protospacer function: {get_protospacer_from_read_group_partial} \nsurrogate function: {get_surrogate_from_read_group_partial}")
                        sequence_counter: ProtospacerSurrogateCounter = Counter()
                        for fastq_single_read_group in fastq_single_read_grouper:
                            protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                            surrogate = get_surrogate_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                            sequence_counter[(protospacer, surrogate)] += 1

                else: # ONLY R1; NO BARCODE; YES UMI
                    if not contains_surrogate:
                        print(f"Performing protospacer+UMI parsing using protospacer function: {get_protospacer_from_read_group_partial} \nUMI function: {get_umi_from_read_group_partial}")
                        sequence_counter: ProtospacerDictUMICounter = defaultdict(Counter)
                        for fastq_single_read_group in fastq_single_read_grouper:
                            protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                            umi = get_umi_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                            sequence_counter[(protospacer,)][umi] += 1
                    else:
                        print(f"Performing protospacer+surrogate+UMI parsing using protospacer function: {get_protospacer_from_read_group_partial} \nsurrogate function: {get_surrogate_from_read_group_partial} \nUMI function: {get_umi_from_read_group_partial}")
                        sequence_counter: ProtospacerSurrogateDictUMICounter = defaultdict(Counter)
                        for fastq_single_read_group in fastq_single_read_grouper:
                            protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                            surrogate = get_surrogate_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                            umi = get_umi_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                            sequence_counter[(protospacer, surrogate)][umi] += 1
            else: # ONLY R1; YES BARCODE
                if not contains_umi: # ONLY R1; YES BARCODE; NO UMI
                    if not contains_surrogate:
                        print(f"Performing protospacer+barcode parsing using protospacer function: {get_protospacer_from_read_group_partial} \nbarcode function: {get_barcode_from_read_group_partial}")
                        sequence_counter: ProtospacerBarcodeCounter = Counter()
                        for fastq_single_read_group in fastq_single_read_grouper:
                            protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                            barcode = get_barcode_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                            sequence_counter[(protospacer, barcode)] += 1
                    else:
                        print(f"Performing protospacer+surrogate+barcode parsing using protospacer function: {get_protospacer_from_read_group_partial} \nsurrogate function: {get_surrogate_from_read_group_partial} \nbarcode function: {get_barcode_from_read_group_partial}")
                        sequence_counter: ProtospacerSurrogateBarcodeCounter = Counter()
                        for fastq_single_read_group in fastq_single_read_grouper:
                            protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                            surrogate = get_surrogate_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                            barcode = get_barcode_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                            sequence_counter[(protospacer, surrogate, barcode)] += 1
                else: # ONLY R1; YES BARCODE; YES UMI
                    if not contains_surrogate:
                        print(f"Performing protospacer+barcode+UMI parsing using protospacer function: {get_protospacer_from_read_group_partial} \nbarcode function: {get_barcode_from_read_group_partial} \nUMI function: {get_umi_from_read_group_partial}")
                        sequence_counter: ProtospacerBarcodeDictUMICounter = defaultdict(Counter)
                        for fastq_single_read_group in fastq_single_read_grouper:
                            protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                            barcode = get_barcode_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                            umi = get_umi_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                            sequence_counter[(barcode, protospacer)][umi] += 1
                    else:
                        print(f"Performing protospacer+surrogate+barcode+UMI parsing using protospacer function: {get_protospacer_from_read_group_partial} \nsurrogate function: {get_surrogate_from_read_group_partial} \nbarcode function: {get_barcode_from_read_group_partial} \nUMI function: {get_umi_from_read_group_partial}")
                        sequence_counter: ProtospacerSurrogateBarcodeDictUMICounter = defaultdict(Counter)
                        for fastq_single_read_group in fastq_single_read_grouper:
                            protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                            surrogate = get_surrogate_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                            barcode = get_barcode_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                            umi = get_umi_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                            sequence_counter[(barcode, surrogate, protospacer)][umi] += 1
        else: # YES R2
            print("Both R1 FASTQ and R2 FASTQ is provided")
            fastq_paired_read_grouper = grouper(zip(fastq_r1_filehandler, fastq_r2_filehandler))
            if not contains_barcode: # YES R2; NO BARCODE
                if not contains_umi: # YES R2; NO BARCODE; NO UMI
                    if not contains_surrogate:
                        print(f"Performing protospacer-only parsing using function: {get_protospacer_from_read_group_partial}")
                        sequence_counter: ProtospacerCounter = Counter()
                        for fastq_paired_read_group in fastq_paired_read_grouper:
                            protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                            sequence_counter[(protospacer,)] += 1
                    else:
                        print(f"Performing protospacer+surrrogate parsing using protospacer function: {get_protospacer_from_read_group_partial} \nsurrogate function: {get_surrogate_from_read_group_partial}")
                        sequence_counter: ProtospacerSurrogateCounter = Counter()
                        for fastq_paired_read_group in fastq_paired_read_grouper:
                            protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                            surrogate = get_surrogate_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                            sequence_counter[(protospacer, surrogate)] += 1
                else: # YES R2; NO BARCODE; YES UMI
                    if not contains_surrogate:
                        print(f"Performing protospacer+UMI parsing using protospacer function: {get_protospacer_from_read_group_partial} \nUMI function: {get_umi_from_read_group_partial}")
                        sequence_counter: ProtospacerDictUMICounter = defaultdict(Counter)
                        for fastq_paired_read_group in fastq_paired_read_grouper:
                            protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                            umi = get_umi_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                            sequence_counter[(protospacer,)][umi] += 1
                    else:
                        print(f"Performing protospacer+surrogate+UMI parsing using protospacer function: {get_protospacer_from_read_group_partial} \nsurrogate function: {get_surrogate_from_read_group_partial} \nUMI function: {get_umi_from_read_group_partial}")
                        sequence_counter: ProtospacerSurrogateDictUMICounter = defaultdict(Counter)
                        for fastq_paired_read_group in fastq_paired_read_grouper:
                            protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                            surrogate = get_surrogate_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                            umi = get_umi_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                            sequence_counter[(protospacer, surrogate)][umi] += 1
            
            else: # YES R2; YES BARCODE
                if umi_pattern_regex is None: # YES R2; YES BARCODE; NO UMI
                    if not contains_surrogate:
                        print(f"Performing protospacer+barcode parsing using protospacer function: {get_protospacer_from_read_group_partial} \nbarcode function: {get_barcode_from_read_group_partial}")
                        sequence_counter: ProtospacerSurrogateBarcodeCounter = Counter()
                        for fastq_paired_read_group in fastq_paired_read_grouper:
                            protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                            barcode = get_barcode_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                            sequence_counter[(protospacer, barcode)] += 1
                    else:
                        print(f"Performing protospacer+surrogate+barcode parsing using protospacer function: {get_protospacer_from_read_group_partial} \nsurrogate function: {get_surrogate_from_read_group_partial} \nbarcode function: {get_barcode_from_read_group_partial}")
                        sequence_counter: ProtospacerSurrogateBarcodeCounter = Counter()
                        for fastq_paired_read_group in fastq_paired_read_grouper:
                            protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                            surrogate = get_surrogate_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                            barcode = get_barcode_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                            sequence_counter[(protospacer, surrogate, barcode)] += 1
                else: # YES R2; YES BARCODE; YES UMI
                    if not contains_surrogate:
                        print(f"Performing protospacer+barcode+UMI parsing using protospacer function: {get_protospacer_from_read_group_partial} \nbarcode function: {get_barcode_from_read_group_partial} \nUMI function: {get_umi_from_read_group_partial}")
                        sequence_counter: ProtospacerSurrogateBarcodeDictUMICounter = defaultdict(Counter)
                        for fastq_paired_read_group in fastq_paired_read_grouper:
                            protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                            barcode = get_barcode_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                            umi = get_umi_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                            sequence_counter[(protospacer, barcode)][umi] += 1
                    else:
                        print(f"Performing protospacer+surrogate+barcode+UMI parsing using protospacer function: {get_protospacer_from_read_group_partial} \nsurrogate function: {get_surrogate_from_read_group_partial} \nbarcode function: {get_barcode_from_read_group_partial} \nUMI function: {get_umi_from_read_group_partial}")
                        sequence_counter: ProtospacerSurrogateBarcodeDictUMICounter = defaultdict(Counter)
                        for fastq_paired_read_group in fastq_paired_read_grouper:
                            protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                            surrogate = get_surrogate_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                            barcode = get_barcode_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                            umi = get_umi_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                            sequence_counter[(protospacer, surrogate, barcode)][umi] += 1
        return sequence_counter
    
    before_file_loading_time = datetime.now()
    # Open R1 file handler
    fastq_r1_filehandler = None
    if fastq_r1_fn.endswith('.gz'):
        print(f"Opening FASTQ.gz file with gzip, filename={fastq_r1_fn}")
        fastq_r1_filehandler = gzip.open(fastq_r1_fn, "rt", encoding="utf-8")
    else:
        print(f"Opening FASTQ file, filename={fastq_r1_fn}")
        fastq_r1_filehandler = open(fastq_r1_fn, "r")
    
    # Open R2 file handler if provided
    fastq_r2_filehandler = None
    if fastq_r2_fn is not None:
        if fastq_r2_fn.endswith('.gz'):
            print(f"Opening FASTQ.gz file with gzip, filename={fastq_r2_fn}")
            fastq_r2_filehandler = gzip.open(fastq_r2_fn, "rt", encoding="utf-8")
        else:
            print(f"Opening FASTQ file, filename={fastq_r2_fn}")
            fastq_r2_filehandler = open(fastq_r2_fn, "r")
    after_file_loading_time = datetime.now()

    print(f"{(after_file_loading_time-before_file_loading_time).seconds} seconds for file loading")
    sequence_counter = parse_fastq(fastq_r1_filehandler, fastq_r2_filehandler)
    after_parsing_time = datetime.now()
    print(f"{(after_parsing_time-after_file_loading_time).seconds} seconds for parsing")

    # Close the file handlers when done
    if fastq_r1_filehandler is not None:
        fastq_r1_filehandler.close()
    if fastq_r2_filehandler is not None:
        fastq_r2_filehandler.close()

    return sequence_counter