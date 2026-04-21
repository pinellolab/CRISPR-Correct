###
### GUIDE PARSING HELPER FUNCTIONS
###
from typeguard import typechecked
from typing import Union, Optional, List, Tuple, DefaultDict
from typing import Counter as CounterType
from ..models.mapping_models import GeneralGuideCountType
from ..models.mapping_models import ProtospacerCounter, ProtospacerDictUMICounter, ProtospacerBarcodeCounter, ProtospacerBarcodeDictUMICounter, ProtospacerSurrogateCounter, ProtospacerSurrogateDictUMICounter, ProtospacerSurrogateBarcodeCounter, ProtospacerSurrogateBarcodeDictUMICounter
from ..models.mapping_models import SampleProtospacerCounter, SampleProtospacerDictUMICounter, SampleProtospacerBarcodeCounter, SampleProtospacerBarcodeDictUMICounter, SampleProtospacerSurrogateCounter, SampleProtospacerSurrogateDictUMICounter, SampleProtospacerSurrogateBarcodeCounter, SampleProtospacerSurrogateBarcodeDictUMICounter
from collections import Counter
from typing import Callable
from functools import partial
from datetime import datetime
import gzip
import re
from functools import lru_cache
from collections import Counter, defaultdict
import logging
_log = logging.getLogger(__name__)


# PERF §3.10: Bio.Seq.reverse_complement goes through a per-base Seq alphabet
# pipeline — ~20× slower than a translate-based revcomp on ACGTN strings. This
# table covers the IUPAC bases that actually appear in FASTQ reads.
_RCMAP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# PERF §3.9: Python's re module has a built-in 512-entry compiled-pattern
# cache, but it's keyed on the raw pattern string and only checked on every
# `re.search` call. A module-local lru_cache on `re.compile` is faster on the
# parsing hot path because the pattern is passed through as an argument every
# iteration.
@lru_cache(maxsize=64)
def _compile_cached(pattern: str) -> "re.Pattern":
    return re.compile(pattern)

# This is for grouping the FASTQ lines in 1 tuple.
def grouper(iterable, n=4):
    "s -> (s0,s1,...sn-1), (sn,sn+1,...s2n-1), (s2n,s2n+1,...s3n-1), ..."
    return zip(*[iter(iterable)]*n)

# TODO 20250908 multisample: Likely will change return descriptor to account for sample barcode
@typechecked
def get_standard_observed_sequence_counts(  fastq_r1_fns: List[str], 
                                            fastq_r2_fns: Optional[List[str]], 
                                            
                                            protospacer_pattern_regex: Optional[str],
                                            surrogate_pattern_regex: Optional[str],
                                            guide_barcode_pattern_regex: Optional[str],
                                            guide_umi_pattern_regex: Optional[str],
                                            sample_barcode_pattern_regex: Optional[str],

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

                                            guide_barcode_left_flank:Optional[str],
                                            guide_barcode_right_flank:Optional[str],
                                            guide_barcode_start_position:Optional[int],
                                            guide_barcode_end_position:Optional[int],
                                            guide_barcode_length: Optional[int],

                                            guide_umi_left_flank:Optional[str],
                                            guide_umi_right_flank:Optional[str],
                                            guide_umi_start_position:Optional[int],
                                            guide_umi_end_position:Optional[int],
                                            guide_umi_length: Optional[int],

                                            sample_barcode_left_flank:Optional[str],
                                            sample_barcode_right_flank:Optional[str],
                                            sample_barcode_start_position:Optional[int],
                                            sample_barcode_end_position:Optional[int],
                                            sample_barcode_length: Optional[int],


                                            is_protospacer_r1: Optional[bool], 
                                            is_surrogate_r1: Optional[bool], 
                                            is_guide_barcode_r1: Optional[bool],
                                            is_guide_umi_r1: Optional[bool],
                                            is_sample_barcode_r1: Optional[bool],

                                            is_protospacer_header: Optional[bool], 
                                            is_surrogate_header: Optional[bool], 
                                            is_guide_barcode_header: Optional[bool],
                                            is_guide_umi_header: Optional[bool],
                                            is_sample_barcode_header: Optional[bool],
                                        
                                            revcomp_protospacer: Optional[bool], 
                                            revcomp_surrogate: Optional[bool], 
                                            revcomp_guide_barcode: Optional[bool],
                                            revcomp_guide_umi: Optional[bool],
                                            revcomp_sample_barcode: Optional[bool],
                                            contains_guide_surrogate: bool,
                                            contains_guide_barcode: bool,
                                            contains_guide_umi: bool,
                                            contains_sample_barcode: bool,) -> GeneralGuideCountType:
    def revcomp(sequence, do_revcomp):
        return sequence.translate(_RCMAP)[::-1] if do_revcomp else sequence

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
            # FIX §1.3: previously `re.search(...).group(1)` raised
            # AttributeError on no-match (opaque to users); return None so the
            # caller treats it as a parse failure and moves on.
            _m = _compile_cached(sequence_pattern_regex).search(template_sequence)
            if _m is None:
                return None
            sequence = _m.group(1).rstrip()
        else:
            """
                Get sequence position start
            """
            if sequence_left_flank is not None:
                """
                    Parse by left flank
                """
                position_start_before_flank = template_sequence.find(sequence_left_flank)
                position_start = position_start_before_flank + len(sequence_left_flank)
                
                if position_start_before_flank == -1:
                    return None # Could not find left_flank, therefore unable to parse
                
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
                if position_end == -1:
                    return None # Could not find right flank, therefore unable to parse
                
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
            
            sequence = template_sequence[position_start:position_end].upper().rstrip() # rstrip is important, since need to remove newlines if sequence_length exceeds template_sequence length

            if (sequence_length is not None) and (len(sequence) != sequence_length):
                return None # The parsed sequence is not the expected length, so return None
        
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
    guide_barcode_parse_sequence_partial = partial(parse_sequence, 
                                                 sequence_pattern_regex=guide_barcode_pattern_regex,
                                                 sequence_left_flank=guide_barcode_left_flank,
                                                 sequence_right_flank=guide_barcode_right_flank,
                                                 sequence_start_position=guide_barcode_start_position,
                                                 sequence_end_position=guide_barcode_end_position,
                                                 sequence_length=guide_barcode_length,
                                                 revcomp_sequence=revcomp_guide_barcode,
                                                 sequence_type="guide_barcode")
    guide_umi_parse_sequence_partial = partial(parse_sequence, 
                                                 sequence_pattern_regex=guide_umi_pattern_regex,
                                                 sequence_left_flank=guide_umi_left_flank,
                                                 sequence_right_flank=guide_umi_right_flank,
                                                 sequence_start_position=guide_umi_start_position,
                                                 sequence_end_position=guide_umi_end_position,
                                                 sequence_length=guide_umi_length,
                                                 revcomp_sequence=revcomp_guide_umi,
                                                 sequence_type="guide_umi")
    
    sample_barcode_parse_sequence_partial = partial(parse_sequence, 
                                                 sequence_pattern_regex=sample_barcode_pattern_regex,
                                                 sequence_left_flank=sample_barcode_left_flank,
                                                 sequence_right_flank=sample_barcode_right_flank,
                                                 sequence_start_position=sample_barcode_start_position,
                                                 sequence_end_position=sample_barcode_end_position,
                                                 sequence_length=sample_barcode_length,
                                                 revcomp_sequence=revcomp_sample_barcode,
                                                 sequence_type="sample_barcode")
    # TODO 09082025 multisample - CREATE sample_barcode parse sequence function.

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
    get_guide_barcode_from_read_group_partial = partial(get_sequence_from_read_group,
                                                      parse_sequence_partial=guide_barcode_parse_sequence_partial, 
                                                      is_r1=is_guide_barcode_r1, 
                                                      is_header=is_guide_barcode_header,
                                                      sequence_type="guide_barcode")
    get_guide_umi_from_read_group_partial = partial(get_sequence_from_read_group,
                                                      parse_sequence_partial=guide_umi_parse_sequence_partial, 
                                                      is_r1=is_guide_umi_r1, 
                                                      is_header=is_guide_umi_header,
                                                      sequence_type="guide_umi")
    get_sample_barcode_from_read_group_partial = partial(get_sequence_from_read_group,
                                                      parse_sequence_partial=sample_barcode_parse_sequence_partial, 
                                                      is_r1=is_sample_barcode_r1, 
                                                      is_header=is_sample_barcode_header,
                                                      sequence_type="sample_barcode")
    
    def parse_fastq(fastq_r1_filehandler, fastq_r2_filehandler, sequence_counter: Optional[GeneralGuideCountType] = None):
        NestedDict = lambda: defaultdict(lambda: Counter())

        # FIX: the original line referenced `fastq_r2_fn` (singular), which is a
        # name left over from the master API and only assigned inside the
        # caller's fastq_index loop when fastq_r2_fns is not None — raising
        # NameError when the caller passes fastq_r2_fns=None (surrogate-off
        # mode). Use the file handler object instead, which is always defined.
        if fastq_r2_filehandler is None: # ONLY R1
            
            _log.info("Only R1 FASTQ is provided, R2 NOT provided")
            fastq_single_read_grouper = grouper(fastq_r1_filehandler)
            
            # CONTAINS SAMPLE BARCODE
            if contains_sample_barcode:
                if not contains_guide_barcode: # ONLY R1; NO BARCODE
                    if not contains_guide_umi: # ONLY R1; NO BARCODE; NO UMI
                        if not contains_guide_surrogate:
                            _log.info(f"Performing protospacer-only parsing using function: {get_protospacer_from_read_group_partial}")
                            if sequence_counter is None:
                                sequence_counter: SampleProtospacerCounter = defaultdict(Counter)
                            for fastq_single_read_group in fastq_single_read_grouper:
                                sample_barcode = get_sample_barcode_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                if (protospacer is not None) and (sample_barcode is not None):
                                    sequence_counter[(protospacer,)][sample_barcode] += 1
                        else:
                            _log.info(f"Performing protospacer+surrrogate parsing using protospacer function: {get_protospacer_from_read_group_partial} \nsurrogate function: {get_surrogate_from_read_group_partial}")
                            if sequence_counter is None:
                                sequence_counter: SampleProtospacerSurrogateCounter = defaultdict(Counter)
                            for fastq_single_read_group in fastq_single_read_grouper:
                                sample_barcode = get_sample_barcode_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                surrogate = get_surrogate_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                if (protospacer is not None) and (surrogate is not None) and (sample_barcode is not None):
                                    sequence_counter[(protospacer, surrogate)][sample_barcode] += 1

                    else: # ONLY R1; NO BARCODE; YES UMI
                        if not contains_guide_surrogate:
                            _log.info(f"Performing protospacer+UMI parsing using protospacer function: {get_protospacer_from_read_group_partial} \nUMI function: {get_guide_umi_from_read_group_partial}")
                            if sequence_counter is None:
                                sequence_counter: SampleProtospacerDictUMICounter = defaultdict(NestedDict)
                            for fastq_single_read_group in fastq_single_read_grouper:
                                sample_barcode = get_sample_barcode_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                umi = get_guide_umi_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                if (protospacer is not None) and (umi is not None) and (sample_barcode is not None):
                                    sequence_counter[(protospacer,)][sample_barcode][umi] += 1
                        else:
                            _log.info(f"Performing protospacer+surrogate+UMI parsing using protospacer function: {get_protospacer_from_read_group_partial} \nsurrogate function: {get_surrogate_from_read_group_partial} \nUMI function: {get_guide_umi_from_read_group_partial}")
                            if sequence_counter is None:
                                sequence_counter: SampleProtospacerSurrogateDictUMICounter = defaultdict(NestedDict)
                            for fastq_single_read_group in fastq_single_read_grouper:
                                sample_barcode = get_sample_barcode_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                surrogate = get_surrogate_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                umi = get_guide_umi_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                if (protospacer is not None) and (surrogate is not None) and (umi is not None) and (sample_barcode is not None):
                                    sequence_counter[(protospacer, surrogate)][sample_barcode][umi] += 1

                else: # ONLY R1; YES BARCODE
                    if not contains_guide_umi: # ONLY R1; YES BARCODE; NO UMI
                        if not contains_guide_surrogate:
                            _log.info(f"Performing protospacer+barcode parsing using protospacer function: {get_protospacer_from_read_group_partial} \nbarcode function: {get_guide_barcode_from_read_group_partial}")
                            if sequence_counter is None:
                                sequence_counter: SampleProtospacerBarcodeCounter = defaultdict(Counter)
                            for fastq_single_read_group in fastq_single_read_grouper:
                                sample_barcode = get_sample_barcode_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                barcode = get_guide_barcode_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                if (protospacer is not None) and (barcode is not None) and (sample_barcode is not None):
                                    sequence_counter[(protospacer, barcode)][sample_barcode] += 1
                        else:
                            _log.info(f"Performing protospacer+surrogate+barcode parsing using protospacer function: {get_protospacer_from_read_group_partial} \nsurrogate function: {get_surrogate_from_read_group_partial} \nbarcode function: {get_guide_barcode_from_read_group_partial}")
                            if sequence_counter is None:
                                sequence_counter: SampleProtospacerSurrogateBarcodeCounter = defaultdict(Counter)
                            for fastq_single_read_group in fastq_single_read_grouper:
                                sample_barcode = get_sample_barcode_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                surrogate = get_surrogate_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                barcode = get_guide_barcode_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                if (protospacer is not None) and (surrogate is not None) and (barcode is not None) and (sample_barcode is not None):
                                    sequence_counter[(protospacer, surrogate, barcode)][sample_barcode] += 1
                    else: # ONLY R1; YES BARCODE; YES UMI
                        if not contains_guide_surrogate:
                            _log.info(f"Performing protospacer+barcode+UMI parsing using protospacer function: {get_protospacer_from_read_group_partial} \nbarcode function: {get_guide_barcode_from_read_group_partial} \nUMI function: {get_guide_umi_from_read_group_partial}")
                            if sequence_counter is None:
                                sequence_counter: SampleProtospacerBarcodeDictUMICounter = defaultdict(NestedDict)
                            for fastq_single_read_group in fastq_single_read_grouper:
                                sample_barcode = get_sample_barcode_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                barcode = get_guide_barcode_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                umi = get_guide_umi_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                if (protospacer is not None) and (barcode is not None) and (umi is not None) and (sample_barcode is not None):
                                    sequence_counter[(protospacer, barcode)][sample_barcode][umi] += 1
                        else:
                            _log.info(f"Performing protospacer+surrogate+barcode+UMI parsing using protospacer function: {get_protospacer_from_read_group_partial} \nsurrogate function: {get_surrogate_from_read_group_partial} \nbarcode function: {get_guide_barcode_from_read_group_partial} \nUMI function: {get_guide_umi_from_read_group_partial}")
                            if sequence_counter is None:
                                sequence_counter: SampleProtospacerSurrogateBarcodeDictUMICounter = defaultdict(NestedDict)
                            for fastq_single_read_group in fastq_single_read_grouper:
                                sample_barcode = get_sample_barcode_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                surrogate = get_surrogate_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                barcode = get_guide_barcode_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                umi = get_guide_umi_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                if (protospacer is not None) and (surrogate is not None) and (barcode is not None) and (umi is not None) and (sample_barcode is not None):
                                    sequence_counter[(protospacer, surrogate,barcode)][sample_barcode][umi] += 1
        

            # (Original: no sample barcode)
            elif not contains_sample_barcode:
                if not contains_guide_barcode: # ONLY R1; NO BARCODE
                    if not contains_guide_umi: # ONLY R1; NO BARCODE; NO UMI
                        if not contains_guide_surrogate:
                            _log.info(f"Performing protospacer-only parsing using function: {get_protospacer_from_read_group_partial}")
                            if sequence_counter is None:
                                sequence_counter: ProtospacerCounter = Counter()
                            for fastq_single_read_group in fastq_single_read_grouper:
                                protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                if protospacer is not None:
                                    sequence_counter[(protospacer,)] += 1
                        else:
                            _log.info(f"Performing protospacer+surrrogate parsing using protospacer function: {get_protospacer_from_read_group_partial} \nsurrogate function: {get_surrogate_from_read_group_partial}")
                            if sequence_counter is None:
                                sequence_counter: ProtospacerSurrogateCounter = Counter()
                            for fastq_single_read_group in fastq_single_read_grouper:
                                protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                surrogate = get_surrogate_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                if (protospacer is not None) and (surrogate is not None):
                                    sequence_counter[(protospacer, surrogate)] += 1

                    else: # ONLY R1; NO BARCODE; YES UMI
                        if not contains_guide_surrogate:
                            _log.info(f"Performing protospacer+UMI parsing using protospacer function: {get_protospacer_from_read_group_partial} \nUMI function: {get_guide_umi_from_read_group_partial}")
                            if sequence_counter is None:
                                sequence_counter: ProtospacerDictUMICounter = defaultdict(Counter)
                            for fastq_single_read_group in fastq_single_read_grouper:
                                protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                umi = get_guide_umi_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                if (protospacer is not None) and (umi is not None):
                                    sequence_counter[(protospacer,)][umi] += 1
                        else:
                            _log.info(f"Performing protospacer+surrogate+UMI parsing using protospacer function: {get_protospacer_from_read_group_partial} \nsurrogate function: {get_surrogate_from_read_group_partial} \nUMI function: {get_guide_umi_from_read_group_partial}")
                            if sequence_counter is None:
                                sequence_counter: ProtospacerSurrogateDictUMICounter = defaultdict(Counter)
                            for fastq_single_read_group in fastq_single_read_grouper:
                                protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                surrogate = get_surrogate_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                umi = get_guide_umi_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                if (protospacer is not None) and (surrogate is not None) and (umi is not None):
                                    sequence_counter[(protospacer, surrogate)][umi] += 1
                else: # ONLY R1; YES BARCODE
                    if not contains_guide_umi: # ONLY R1; YES BARCODE; NO UMI
                        if not contains_guide_surrogate:
                            _log.info(f"Performing protospacer+barcode parsing using protospacer function: {get_protospacer_from_read_group_partial} \nbarcode function: {get_guide_barcode_from_read_group_partial}")
                            if sequence_counter is None:
                                sequence_counter: ProtospacerBarcodeCounter = Counter()
                            for fastq_single_read_group in fastq_single_read_grouper:
                                protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                barcode = get_guide_barcode_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                if (protospacer is not None) and (barcode is not None):
                                    sequence_counter[(protospacer, barcode)] += 1
                        else:
                            _log.info(f"Performing protospacer+surrogate+barcode parsing using protospacer function: {get_protospacer_from_read_group_partial} \nsurrogate function: {get_surrogate_from_read_group_partial} \nbarcode function: {get_guide_barcode_from_read_group_partial}")
                            if sequence_counter is None:
                                sequence_counter: ProtospacerSurrogateBarcodeCounter = Counter()
                            for fastq_single_read_group in fastq_single_read_grouper:
                                protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                surrogate = get_surrogate_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                barcode = get_guide_barcode_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                if (protospacer is not None) and (surrogate is not None) and (barcode is not None):
                                    sequence_counter[(protospacer, surrogate, barcode)] += 1
                    else: # ONLY R1; YES BARCODE; YES UMI
                        if not contains_guide_surrogate:
                            _log.info(f"Performing protospacer+barcode+UMI parsing using protospacer function: {get_protospacer_from_read_group_partial} \nbarcode function: {get_guide_barcode_from_read_group_partial} \nUMI function: {get_guide_umi_from_read_group_partial}")
                            if sequence_counter is None:
                                sequence_counter: ProtospacerBarcodeDictUMICounter = defaultdict(Counter)
                            for fastq_single_read_group in fastq_single_read_grouper:
                                protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                barcode = get_guide_barcode_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                umi = get_guide_umi_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                if (protospacer is not None) and (barcode is not None) and (umi is not None):
                                    sequence_counter[(protospacer, barcode)][umi] += 1
                        else:
                            _log.info(f"Performing protospacer+surrogate+barcode+UMI parsing using protospacer function: {get_protospacer_from_read_group_partial} \nsurrogate function: {get_surrogate_from_read_group_partial} \nbarcode function: {get_guide_barcode_from_read_group_partial} \nUMI function: {get_guide_umi_from_read_group_partial}")
                            if sequence_counter is None:
                                sequence_counter: ProtospacerSurrogateBarcodeDictUMICounter = defaultdict(Counter)
                            for fastq_single_read_group in fastq_single_read_grouper:
                                protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                surrogate = get_surrogate_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                barcode = get_guide_barcode_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                umi = get_guide_umi_from_read_group_partial(fastq_read_group=fastq_single_read_group, includes_r2=False)
                                if (protospacer is not None) and (surrogate is not None) and (barcode is not None) and (umi is not None):
                                    sequence_counter[(protospacer, surrogate,barcode)][umi] += 1
        
        else: # YES R2
        
            _log.info("Both R1 FASTQ and R2 FASTQ is provided")
            fastq_paired_read_grouper = grouper(zip(fastq_r1_filehandler, fastq_r2_filehandler))
            
            if contains_sample_barcode:
                if not contains_guide_barcode: # YES R2; NO BARCODE
                    if not contains_guide_umi: # YES R2; NO BARCODE; NO UMI
                        if not contains_guide_surrogate:
                            _log.info(f"Performing protospacer-only parsing using function: {get_protospacer_from_read_group_partial}")
                            if sequence_counter is None:
                                sequence_counter: SampleProtospacerCounter = defaultdict(Counter)
                            for fastq_paired_read_group in fastq_paired_read_grouper:
                                sample_barcode = get_sample_barcode_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                if (protospacer is not None) and (sample_barcode is not None):
                                    sequence_counter[(protospacer,)][sample_barcode] += 1
                        else:
                            _log.info(f"Performing protospacer+surrogate parsing using protospacer function: {get_protospacer_from_read_group_partial} \nsurrogate function: {get_surrogate_from_read_group_partial}")
                            if sequence_counter is None:
                                sequence_counter: SampleProtospacerSurrogateCounter = defaultdict(Counter)
                            for fastq_paired_read_group in fastq_paired_read_grouper:
                                sample_barcode = get_sample_barcode_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                surrogate = get_surrogate_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                if (protospacer is not None) and (surrogate is not None) and (sample_barcode is not None):
                                    sequence_counter[(protospacer, surrogate)][sample_barcode] += 1
                    else: # YES R2; NO BARCODE; YES UMI
                        if not contains_guide_surrogate:
                            _log.info(f"Performing protospacer+UMI parsing using protospacer function: {get_protospacer_from_read_group_partial} \nUMI function: {get_guide_umi_from_read_group_partial}")
                            if sequence_counter is None:
                                sequence_counter: SampleProtospacerDictUMICounter = defaultdict(NestedDict)
                            for fastq_paired_read_group in fastq_paired_read_grouper:
                                sample_barcode = get_sample_barcode_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                umi = get_guide_umi_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                if (protospacer is not None) and (umi is not None) and (sample_barcode is not None):
                                    sequence_counter[(protospacer,)][sample_barcode][umi] += 1
                        else:
                            _log.info(f"Performing protospacer+surrogate+UMI parsing using protospacer function: {get_protospacer_from_read_group_partial} \nsurrogate function: {get_surrogate_from_read_group_partial} \nUMI function: {get_guide_umi_from_read_group_partial}")
                            if sequence_counter is None:
                                sequence_counter: SampleProtospacerSurrogateDictUMICounter = defaultdict(NestedDict)
                            for fastq_paired_read_group in fastq_paired_read_grouper:
                                sample_barcode = get_sample_barcode_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                surrogate = get_surrogate_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                umi = get_guide_umi_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                if (protospacer is not None) and (surrogate is not None) and (umi is not None) and (sample_barcode is not None):
                                    sequence_counter[(protospacer, surrogate)][sample_barcode][umi] += 1
                
                else: # YES R2; YES BARCODE
                    if guide_umi_pattern_regex is None: # YES R2; YES BARCODE; NO UMI
                        if not contains_guide_surrogate:
                            _log.info(f"Performing protospacer+barcode parsing using protospacer function: {get_protospacer_from_read_group_partial} \nbarcode function: {get_guide_barcode_from_read_group_partial}")
                            if sequence_counter is None:
                                sequence_counter: SampleProtospacerSurrogateBarcodeCounter = defaultdict(Counter)
                            for fastq_paired_read_group in fastq_paired_read_grouper:
                                sample_barcode = get_sample_barcode_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                barcode = get_guide_barcode_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                if (protospacer is not None) and (barcode is not None) and (sample_barcode is not None):
                                    sequence_counter[(protospacer, barcode)][sample_barcode] += 1
                        else:
                            _log.info(f"Performing protospacer+surrogate+barcode parsing using protospacer function: {get_protospacer_from_read_group_partial} \nsurrogate function: {get_surrogate_from_read_group_partial} \nbarcode function: {get_guide_barcode_from_read_group_partial}")
                            if sequence_counter is None:
                                sequence_counter: SampleProtospacerSurrogateBarcodeCounter = defaultdict(Counter)
                            for fastq_paired_read_group in fastq_paired_read_grouper:
                                sample_barcode = get_sample_barcode_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                surrogate = get_surrogate_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                barcode = get_guide_barcode_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                if (protospacer is not None) and (surrogate is not None) and (barcode is not None) and (sample_barcode is not None):
                                    sequence_counter[(protospacer, surrogate, barcode)][sample_barcode] += 1
                    else: # YES R2; YES BARCODE; YES UMI
                        if not contains_guide_surrogate:
                            _log.info(f"Performing protospacer+barcode+UMI parsing using protospacer function: {get_protospacer_from_read_group_partial} \nbarcode function: {get_guide_barcode_from_read_group_partial} \nUMI function: {get_guide_umi_from_read_group_partial}")
                            if sequence_counter is None:
                                sequence_counter: SampleProtospacerSurrogateBarcodeDictUMICounter = defaultdict(NestedDict)
                            for fastq_paired_read_group in fastq_paired_read_grouper:
                                sample_barcode = get_sample_barcode_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                barcode = get_guide_barcode_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                umi = get_guide_umi_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                if (protospacer is not None) and (barcode is not None) and (umi is not None) and (sample_barcode is not None):
                                    sequence_counter[(protospacer, barcode)][sample_barcode][umi] += 1
                        else:
                            _log.info(f"Performing protospacer+surrogate+barcode+UMI parsing using protospacer function: {get_protospacer_from_read_group_partial} \nsurrogate function: {get_surrogate_from_read_group_partial} \nbarcode function: {get_guide_barcode_from_read_group_partial} \nUMI function: {get_guide_umi_from_read_group_partial}")
                            if sequence_counter is None:
                                sequence_counter: SampleProtospacerSurrogateBarcodeDictUMICounter = defaultdict(NestedDict)
                            for fastq_paired_read_group in fastq_paired_read_grouper:
                                sample_barcode = get_sample_barcode_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                surrogate = get_surrogate_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                barcode = get_guide_barcode_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                umi = get_guide_umi_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                if (protospacer is not None) and (surrogate is not None) and (barcode is not None) and (umi is not None) and (sample_barcode is not None):
                                    sequence_counter[(protospacer, surrogate, barcode)][sample_barcode][umi] += 1


            # (Original: no sample barcode)
            if not contains_sample_barcode:
                if not contains_guide_barcode: # YES R2; NO BARCODE
                    if not contains_guide_umi: # YES R2; NO BARCODE; NO UMI
                        if not contains_guide_surrogate:
                            _log.info(f"Performing protospacer-only parsing using function: {get_protospacer_from_read_group_partial}")
                            if sequence_counter is None:
                                sequence_counter: ProtospacerCounter = Counter()
                            for fastq_paired_read_group in fastq_paired_read_grouper:
                                protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                if protospacer is not None:
                                    sequence_counter[(protospacer,)] += 1
                        else:
                            _log.info(f"Performing protospacer+surrrogate parsing using protospacer function: {get_protospacer_from_read_group_partial} \nsurrogate function: {get_surrogate_from_read_group_partial}")
                            if sequence_counter is None:
                                sequence_counter: ProtospacerSurrogateCounter = Counter()
                            for fastq_paired_read_group in fastq_paired_read_grouper:
                                protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                surrogate = get_surrogate_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                if (protospacer is not None) and (surrogate is not None):
                                    sequence_counter[(protospacer, surrogate)] += 1
                    else: # YES R2; NO BARCODE; YES UMI
                        if not contains_guide_surrogate:
                            _log.info(f"Performing protospacer+UMI parsing using protospacer function: {get_protospacer_from_read_group_partial} \nUMI function: {get_guide_umi_from_read_group_partial}")
                            if sequence_counter is None:
                                sequence_counter: ProtospacerDictUMICounter = defaultdict(Counter)
                            for fastq_paired_read_group in fastq_paired_read_grouper:
                                protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                umi = get_guide_umi_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                if (protospacer is not None) and (umi is not None):
                                    sequence_counter[(protospacer,)][umi] += 1
                        else:
                            _log.info(f"Performing protospacer+surrogate+UMI parsing using protospacer function: {get_protospacer_from_read_group_partial} \nsurrogate function: {get_surrogate_from_read_group_partial} \nUMI function: {get_guide_umi_from_read_group_partial}")
                            if sequence_counter is None:
                                sequence_counter: ProtospacerSurrogateDictUMICounter = defaultdict(Counter)
                            for fastq_paired_read_group in fastq_paired_read_grouper:
                                protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                surrogate = get_surrogate_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                umi = get_guide_umi_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                if (protospacer is not None) and (surrogate is not None) and (umi is not None):
                                    sequence_counter[(protospacer, surrogate)][umi] += 1
                
                else: # YES R2; YES BARCODE
                    if guide_umi_pattern_regex is None: # YES R2; YES BARCODE; NO UMI
                        if not contains_guide_surrogate:
                            _log.info(f"Performing protospacer+barcode parsing using protospacer function: {get_protospacer_from_read_group_partial} \nbarcode function: {get_guide_barcode_from_read_group_partial}")
                            if sequence_counter is None:
                                sequence_counter: ProtospacerSurrogateBarcodeCounter = Counter()
                            for fastq_paired_read_group in fastq_paired_read_grouper:
                                protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                barcode = get_guide_barcode_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                if (protospacer is not None) and (barcode is not None):
                                    sequence_counter[(protospacer, barcode)] += 1
                        else:
                            _log.info(f"Performing protospacer+surrogate+barcode parsing using protospacer function: {get_protospacer_from_read_group_partial} \nsurrogate function: {get_surrogate_from_read_group_partial} \nbarcode function: {get_guide_barcode_from_read_group_partial}")
                            if sequence_counter is None:
                                sequence_counter: ProtospacerSurrogateBarcodeCounter = Counter()
                            for fastq_paired_read_group in fastq_paired_read_grouper:
                                protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                surrogate = get_surrogate_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                barcode = get_guide_barcode_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                if (protospacer is not None) and (surrogate is not None) and (barcode is not None):
                                    sequence_counter[(protospacer, surrogate, barcode)] += 1
                    else: # YES R2; YES BARCODE; YES UMI
                        if not contains_guide_surrogate:
                            _log.info(f"Performing protospacer+barcode+UMI parsing using protospacer function: {get_protospacer_from_read_group_partial} \nbarcode function: {get_guide_barcode_from_read_group_partial} \nUMI function: {get_guide_umi_from_read_group_partial}")
                            if sequence_counter is None:
                                sequence_counter: ProtospacerSurrogateBarcodeDictUMICounter = defaultdict(Counter)
                            for fastq_paired_read_group in fastq_paired_read_grouper:
                                protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                barcode = get_guide_barcode_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                umi = get_guide_umi_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                if (protospacer is not None) and (barcode is not None) and (umi is not None):
                                    sequence_counter[(protospacer, barcode)][umi] += 1
                        else:
                            _log.info(f"Performing protospacer+surrogate+barcode+UMI parsing using protospacer function: {get_protospacer_from_read_group_partial} \nsurrogate function: {get_surrogate_from_read_group_partial} \nbarcode function: {get_guide_barcode_from_read_group_partial} \nUMI function: {get_guide_umi_from_read_group_partial}")
                            if sequence_counter is None:
                                sequence_counter: ProtospacerSurrogateBarcodeDictUMICounter = defaultdict(Counter)
                            for fastq_paired_read_group in fastq_paired_read_grouper:
                                protospacer = get_protospacer_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                surrogate = get_surrogate_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                barcode = get_guide_barcode_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                umi = get_guide_umi_from_read_group_partial(fastq_read_group=fastq_paired_read_group, includes_r2=True)
                                if (protospacer is not None) and (surrogate is not None) and (barcode is not None) and (umi is not None):
                                    sequence_counter[(protospacer, surrogate, barcode)][umi] += 1
        
        return sequence_counter
    


    # Iterate through FASTQs
    sequence_counter: Optional[GeneralGuideCountType] = None
    for fastq_index in range(len(fastq_r1_fns)):
        _log.info(f"Processing FASTQ {fastq_index} of {len(fastq_r1_fns)}")

        before_file_loading_time = datetime.now()
        fastq_r1_fn = fastq_r1_fns[fastq_index]
        
        # Open R1 file handler
        fastq_r1_filehandler = None
        if fastq_r1_fn.endswith('.gz'):
            _log.info(f"Opening FASTQ.gz file with gzip, filename={fastq_r1_fn}")
            fastq_r1_filehandler = gzip.open(fastq_r1_fn, "rt", encoding="utf-8")
        else:
            _log.info(f"Opening FASTQ file, filename={fastq_r1_fn}")
            fastq_r1_filehandler = open(fastq_r1_fn, "r")
        
        # Open R2 file handler if provided
        fastq_r2_filehandler = None
        if fastq_r2_fns is not None:
            fastq_r2_fn = fastq_r2_fns[fastq_index]

            if fastq_r2_fn.endswith('.gz'):
                _log.info(f"Opening FASTQ.gz file with gzip, filename={fastq_r2_fn}")
                fastq_r2_filehandler = gzip.open(fastq_r2_fn, "rt", encoding="utf-8")
            else:
                _log.info(f"Opening FASTQ file, filename={fastq_r2_fn}")
                fastq_r2_filehandler = open(fastq_r2_fn, "r")
        after_file_loading_time = datetime.now()

        _log.info(f"{(after_file_loading_time-before_file_loading_time).seconds} seconds for file loading")
        sequence_counter = parse_fastq(fastq_r1_filehandler, fastq_r2_filehandler, sequence_counter)
        after_parsing_time = datetime.now()
        _log.info(f"{(after_parsing_time-after_file_loading_time).seconds} seconds for parsing")

        # Close the file handlers when done
        if fastq_r1_filehandler is not None:
            fastq_r1_filehandler.close()
        if fastq_r2_filehandler is not None:
            fastq_r2_filehandler.close()

    return sequence_counter