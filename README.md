![PyPI - Version](https://img.shields.io/pypi/v/crispr-ambiguous-mapping)

# CRISPR-Correct
**CRISPR-Correct** was developed by the *Pinello Lab* as an easy-to-use Python package for performing guide RNA mapping given the raw FASTQs and guide library dataframe. Specifically, <ins>CRISPR-Correct is able to handle imperfect mapping </ins> (either due to self-editing using SpRY-based base-editors or sequencing errors) by mapping observed protospacer sequences to the guide RNA with the closest hamming distance; CRISPR-Correct is **unable** to handle in-dels in the protospacer (consider using another *Pinello Lab* tool [CRISPR BEAN]()). If you aren't expecting self-editing or sequencing error, this tool will still work though it may be easier and faster to use another guide mapping tool such as *Pinello Lab* tool [CRISPR SURF](https://github.com/pinellolab/CRISPR-SURF/).

<ins>CRISPR-Correct is also able to handle guide RNA sensor constructs and UMIs</ins> (see **Figure 1** below). Specifically, the mapping is performed on each protospacer/surrogate/barcode permutation, and the editing outcomes at the protospacer and surrogate is characterized. You will just need to provide regex strings for extracting the protospacer/surrogate/barcode sequences from the read.

<img src="https://github.com/user-attachments/assets/e79e2dc2-ae96-4328-a1b8-06c9ec7a36ae" alt="CRISPR-CLEAR framework" width="300"></img>

<em>**Figure 1.** Schematic of guide RNA sensor construct. Typically, the expressed guide RNA will edit both the endogenous target site and the surrogate target site. Typically, a short barcode will be added to distinguish between any similar protospacer/surrogate sequences in the guide RNA library to facilitate guide RNA mapping. In this design, paired-end sequencing is performed to sequence the protospacer (in the R1 read) and the surrogate and barcode (in the R2 read).</em>

If you are mapping several large samples that would take too long to run on a personal computer, CRISPR-Correct can also be ran on [Broad Institute's Terra Platform](https://terra.bio/)! The workflow file is located at the [Terra Firecloud repository](https://portal.firecloud.org/?return=terra#methods/pinellolab/CrisprSelfEditMappingOrchestratorWorkflowSampleEntity/2).

### Installation
CRISPR-Correct can be easily installed from PyPi `pip install crispr-ambiguous-mapping`, which should only take a couple minutes. CRISPR-Correct requires **Python versions >=3.8,<3.12** which can be installed from the [Python download page](https://www.python.org/downloads/). 

### System Requirements
CRISPR-Correct can run on [any operating system where Python versions >=3.8,<3.12 can be installed](https://www.python.org/downloads/operating-systems/). To speed up model performance, CRISPR-Correct can utilize multiple CPUs (for multi-threading) and is highly recommended especially for large samples. Depending on how large the sample is, it is suggested that there is sufficient RAM to allow in-memory storage of the guide RNA mapping results. Additionally, disk I/O speed can **substantially** increase mapping performance (such as using a solid state drive over a hard disk drive).

## Prepare inputs to guide mapping
Prior to running the tool, you need to spend time preparing the relevant inputs: 1) the pre-processed R1 (and R2) demultiplexed FASTQs, 2) guide library table, 3) regex strings to parse barcode and UMI sequences from FASTQ header if provided.

### R1 and R2 demultiplexed FASTQs
You will need to provide demultiplexed R1 (and R2) FASTQs where the R1 FASTQ sequence **only** contains the protospacer sequence and (if sequenced) the R2 FASTQ sequece **only** contains the surrogate sequence if included in the guide RNA design. The FASTQ header should contain the barcode and UMI (if included in the guide RNA design), which will be parsed by the input regex strings (discussed in later section). You can use any tool to accomplish ths, you can use [**UMITools**](https://github.com/CGATOxford/UMI-tools) to parse the protospacer/surrogate/barcode/UMI sequences and especially any in-read index sequence for demultiplexing. We have then used [**BBMap** ](https://github.com/BioInfoTools/BBMap/blob/master/sh/demuxbyname.sh) to demultiplex the UMITools-processed FASTQs, which allowed us to demultiplex by the in-read index containing the the FASTQ header. We have provided a [Terra Firecloud method](https://portal.firecloud.org/?return=terra#methods/pinellolab/CrisprMillipedeGuideDemultiplex/5) for performing both UMITools and BBMap especially for large samples, however this method may still be useful to check out if running on a personal computer. See [the example pipeline JSON file](https://github.com/pinellolab/CRISPR-Correct/blob/master/examples/demultiplex_input.json) for an example input file for the demultiplexing and pre-processing pipeline.

### Barcode and UMI regex strings
The barcode and UMI sequence will be parsed from the FASTQ read headers, although you will need to prepare the regex strings for parsing, i.e.:

Example header:

`@lh00134:140:225VLGLT3:7:1101:1028:1080_ANGC_GGCA 1:N:0:GAAATAAG+ACGTCCTG`

Example barcode regex extracting "GGCA" barcode sequence:

`BARCODE_REGEX = r"_([^_ ]+)[\s+]"`

Example UMI regex extracting "GAAATAAG" UMI sequence:

`UMI_REGEX = r":([^+:]{6})(.{2})\+"`

### Guide library table

You will need a TSV file containing the list of guides for mapping. The table will contain three columns: protospacer (required), surrogate (optional), barcode (optional). The lengths of the protospacer, surrogate, and barcode can each be any length, however, all values of each column must be the same length. The guide mapping will directly map the R1 protospacer sequence to the protospacer column, the R2 surrogate sequence to the surrogate column, and the header barcode sequence to the barcode column.

Example "guide_library.tsv"
```
protospacer	surrogate	barcode
TGTCGTGAGGTAGCTACGAC	CAGCAATGTCGTGAGGTAGCTACGACTTGTCA	GCTC
AGTCGTAGCTACCTCACGAC	ATGACAAGTCGTAGCTACCTCACGACATTGCT	GTTG
CCTAGTGGTTATTCGATGTC	AGGTTACCTAGTGGTTATTCGATGTCTCAGAA	CGAA
...
```

### Other inputs

You will also need to provide the function the hamming distance threshold for guide mapping, i.e. if there is no guide with a minimum hamming distance as the threshold, the sequence will be thrown out.

For a protospacer length of 20, surrogate length of 32, and barcode length of 4, the defaults used are:
```
PROTOSPACER_HAMMING_THRESHOLD = 7
SURROGATE_HAMMING_THRESHOLD = 10 
BARCODE_HAMMING_THRESHOLD = 2
```
## Running the guide mapping

### How is mapping performed?
For each sample, the tool will iterate through the R1 (and R2) sequences and extract the protospacer/surrogate/barcode/UMI sequences based on the provided regex strings. To prevent redundant mapping of the same sequences, the tool will create a Python *Counter* object to store unique sequences as the key and their read counts as the value. For each unique observed sequence, the hamming distance for the protospacer/surrogate/barcode sequence will be calculated for all sequences in the provided guide library in an efficient vectorized manner. The sequence with the minimum hamming distance is the mapped sequence (unless the minimum hamming distance is above the thresholds provided). If there are multiple sequences with the same minimum hamming distance (ambiguous mapping), the mapping count is either assigned to all ambiguous mappings, spread across the ambiguous mappings, or counted to neither of the ambiguous mappings.  The protospacer/surrogate/barcode sequences are mapped separately to diagnose any mismatching between protospacer and surrogate sequences (due to possible recombination effects). 

### Run the mapping function

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

### See the mapping outputs
You can retrieve the count results at:
`count_result.all_match_set_whitelist_reporter_counter_series_results`


