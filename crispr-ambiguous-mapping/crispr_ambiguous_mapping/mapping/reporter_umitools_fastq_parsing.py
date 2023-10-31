###
### GUIDE PARSING HELPER FUNCTIONS
###
from typeguard import typechecked
from typing import Union, Optional, List
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

@typechecked
def retrieve_fastq_guide_sequences(r1_protospacer_fastq_file: str, r2_surrogate_fastq_file: Optional[str] = None, barcode_pattern_regex: Optional[str] = None, i7_pattern_regex: Optional[str] = None):
    def parse_fastq(r1_protospacer_fastq_filehandler, r2_surrogate_fastq_filehandler):
        sequence_counter_dict = defaultdict(Counter)

        # If the R2 is also provided
        if r2_surrogate_fastq_file is not None:
            
            # Iterate through each read group (see grouper function above)
            for fastq_paired_read_group in grouper(zip(r1_protospacer_fastq_filehandler, r2_surrogate_fastq_filehandler)):
                r1_header_read, r2_header_read = fastq_paired_read_group[0]
                r1_sequence_read, r2_sequence_read = fastq_paired_read_group[1]
                
                protospacer=r1_sequence_read
                surrogate=r2_sequence_read
                barcode = None
                i7_index = None

                # Perform parsing
                if barcode_pattern_regex:
                    barcode = re.search(barcode_pattern_regex, r1_header_read)
                if i7_pattern_regex:
                    i7_index = re.search(i7_pattern_regex, r1_header_read)
                
                # Count
                sequence_counter_dict[(protospacer, surrogate, barcode)][i7_index] += 1
        else: # If only the R1 is provided
            for fastq_single_read_group in grouper(r1_protospacer_fastq_filehandler):
                r1_header_read = fastq_paired_read_group[0]
                r1_sequence_read = fastq_paired_read_group[1]
                
                protospacer=r1_sequence_read
                surrogate=None
                barcode = None
                i7_index = None

                if barcode_pattern_regex:
                    barcode = re.search(barcode_pattern_regex, r1_header_read)
                if i7_pattern_regex:
                    i7_index = re.search(i7_pattern_regex, r1_header_read)
                
                sequence_counter_dict[(protospacer, surrogate, barcode)][i7_index] += 1
        return sequence_counter_dict
    
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

    sequence_counter_dict = parse_fastq(r1_protospacer_fastq_filehandler, r2_surrogate_fastq_filehandler)

    # Close the file handlers when done
    if r1_protospacer_fastq_filehandler is not None:
        r1_protospacer_fastq_filehandler.close()
    if r2_surrogate_fastq_filehandler is not None:
        r2_surrogate_fastq_filehandler.close()

    return sequence_counter_dict