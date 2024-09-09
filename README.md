# CRISPR-Correct
**CRISPR-Correct** is a Python package for performing guide RNA mapping given the raw FASTQs and guide library dataframe. Specifically, <ins>CRISPR-Correct is able to handle imperfect mapping </ins> (due to self-editing using SpRY-based base-editors or sequencing errors) by mapping observed protospacer sequences to the guide RNA with the closest hamming distance; CRISPR-Correct is unable to handle in-dels in the protospacer.

For several large samples, CRISPR-Correct can also be ran on Broad Institute's Terra Platform. The workflow file is located at: https://portal.firecloud.org/?return=terra#methods/pinellolab/CrisprSelfEditMappingOrchestratorWorkflowSampleEntity/2

## Installation
CRISPR-Correct can be installed from PyPi `pip install crispr-ambiguous-mapping`. CRISPR-Correct requires Python versions >=3.8,<3.12.

## Prepare inputs to guide mapping

You will need a TSV file containing the list of guides for mapping. The table will contain three columns: protospacer (required), surrogate (optional), barcde (optional). The lengths of the protospacer, surrogate, and barcode can be any length, however, all values of each column must be the same length.

Example "guide_library.tsv"
```
protospacer	surrogate	barcode
TGTCGTGAGGTAGCTACGAC	CAGCAATGTCGTGAGGTAGCTACGACTTGTCA	GCTC
AGTCGTAGCTACCTCACGAC	ATGACAAGTCGTAGCTACCTCACGACATTGCT	GTTG
CCTAGTGGTTATTCGATGTC	AGGTTACCTAGTGGTTATTCGATGTCTCAGAA	CGAA
```

Additionally, you will need the demultiplexed R1 (required) and R2 (optional) FASTQ files, where the R1 contains the protospacer sequence and the R2 contains the surrogate sequence. The barcode and UMI sequence will be parsed from the FASTQ read headers, although you will need to prepare the regex strings for parsing, i.e.:

Example header:

`@lh00134:140:225VLGLT3:7:1101:1028:1080_ANGC_GGCA 1:N:0:GAAATAAG+ACGTCCTG`

Example barcode regex extracting "GGCA":

`BARCODE_REGEX = r"_([^_ ]+)[\s+]"`

Example UMI regex extracting "GAAATAAG":

`UMI_REGEX = r":([^+:]{6})(.{2})\+"`

You will also need to provide the function the hamming distance threshold for guide mapping, i.e. if there is no guide with a minimum hamming distance as the threshold, the sequence will be thrown out.

For a protospacer length of 20, surrogate length of 32, and barcode length of 4, the defaults used are:
```
PROTOSPACER_HAMMING_THRESHOLD = 7
SURROGATE_HAMMING_THRESHOLD = 10 
BARCODE_HAMMING_THRESHOLD = 2
```

## Running the guide mapping

Import the library "crispr_ambiguous_mapping" and pandas, load in the guide library TSV, and pass in the inputs arguments (along with the number of cores for parallelization).

You can then run the guide mapping function:

```
import crispr_ambiguous_mapping
import pandas as pd
        
GUIDE_LIBRARY_DATAFRAME = pd.read_table("{GUIDE_LIBRARY_TSV}")

count_result = crispr_ambiguous_mapping.mapping.get_whitelist_reporter_counts_from_umitools_output(
            whitelist_guide_reporter_df=GUIDE_LIBRARY_DATAFRAME, 
            fastq_r1_fn='{READ1_FASTQ}', 
            fastq_r2_fn='{READ1_FASTQ}',
            barcode_pattern_regex={BARCODE_REGEX},
            umi_pattern_regex={UMI_REGEX},
            surrogate_hamming_threshold_strict={SURROGATE_HAMMING_THRESHOLD},
            barcode_hamming_threshold_strict={BARCODE_HAMMING_THRESHOLD},
            protospacer_hamming_threshold_strict={PROTOSPACER_HAMMING_THRESHOLD},
            cores={CPUS})
```

You can retrieve the count results at:
`count_result.all_match_set_whitelist_reporter_counter_series_results`


