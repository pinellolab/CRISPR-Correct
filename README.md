![PyPI - Version](https://img.shields.io/pypi/v/crispr-ambiguous-mapping)

# CRISPR-Correct
**CRISPR-Correct** was developed by the *Pinello Lab* as an easy-to-use Python package for performing guide RNA mapping given the raw FASTQs and guide library dataframe. Specifically, <ins>CRISPR-Correct is able to handle imperfect mapping </ins> (either due to self-editing using SpRY-based base-editors or sequencing errors) by mapping observed protospacer sequences to the guide RNA with the closest hamming distance; CRISPR-Correct is **unable** to handle in-dels in the protospacer (consider using another *Pinello Lab* tool [CRISPR BEAN]()). If you aren't expecting self-editing or sequencing error, this tool will still work though it may be easier and faster to use another guide mapping tool such as *Pinello Lab* tool [CRISPR SURF](https://github.com/pinellolab/CRISPR-SURF/).

<ins>CRISPR-Correct is also able to handle guide RNA sensor constructs and UMIs</ins> (see **Figure 1** below). Specifically, the mapping is performed on each protospacer/surrogate/barcode permutation, and the editing outcomes at the protospacer and surrogate is characterized. You will just need to provide regex strings for extracting the protospacer/surrogate/barcode sequences from the read.

<img src="https://github.com/user-attachments/assets/e79e2dc2-ae96-4328-a1b8-06c9ec7a36ae" alt="CRISPR-CLEAR framework" width="300"></img>

<em>**Figure 1.** Schematic of guide RNA sensor construct. Typically, the expressed guide RNA will edit both the endogenous target site and the surrogate target site. Typically, a short barcode will be added to distinguish between any similar protospacer/surrogate sequences in the guide RNA library to facilitate guide RNA mapping. In this design, paired-end sequencing is performed to sequence the protospacer (in the R1 read) and the surrogate and barcode (in the R2 read).</em>

If you are mapping several large samples that would take too long to run on a personal computer, CRISPR-Correct can also be ran on [Broad Institute's Terra Platform](https://terra.bio/)! The workflow file is located at the [Terra Firecloud repository](https://portal.firecloud.org/?return=terra#methods/pinellolab/CrisprSelfEditMappingOrchestratorWorkflowSampleEntity/2).

### Installation
CRISPR-Correct can be easily installed from PyPi `pip install crispr-ambiguous-mapping==0.0.177`, which should only take a couple minutes. CRISPR-Correct requires **Python versions >=3.8,<3.12** which can be installed from the [Python download page](https://www.python.org/downloads/). 

### System Requirements
CRISPR-Correct can run on [any operating system where Python versions >=3.8,<3.12 can be installed](https://www.python.org/downloads/operating-systems/). To speed up model performance, CRISPR-Correct can utilize multiple CPUs (for multi-threading) and is highly recommended especially for large samples. Depending on how large the sample is, it is suggested that there is sufficient RAM to allow in-memory storage of the guide RNA mapping results. Additionally, disk I/O speed can **substantially** increase mapping performance (such as using a solid state drive over a hard disk drive).

## Prepare inputs to guide mapping
Prior to running the tool, you need to spend time preparing the relevant inputs: 1) the R1 (and R2) demultiplexed FASTQs, 2) guide library table, 3) specifications for parsing the protospacer (and surrogate, barcode, and UMI if provided).

### R1 and R2 demultiplexed FASTQs
You will need to provide demultiplexed R1 (and R2) FASTQs. 

**Tip:** If you need to parse the sequence before demultiplexing (i.e. in-read index), you could use [**UMITools**](https://github.com/CGATOxford/UMI-tools) to parse the protospacer/surrogate/barcode/UMI sequences and especially any in-read index sequence for demultiplexing. We have then used [**BBMap** ](https://github.com/BioInfoTools/BBMap/blob/master/sh/demuxbyname.sh) to demultiplex the UMITools-processed FASTQs, which allowed us to demultiplex by the in-read index containing the the FASTQ header. We have provided a [Terra Firecloud method](https://portal.firecloud.org/?return=terra#methods/pinellolab/CrisprMillipedeGuideDemultiplex/5) for performing both UMITools and BBMap especially for large samples, however this method may still be useful to check out if running on a personal computer. See [the example pipeline JSON file](https://github.com/pinellolab/CRISPR-Correct/blob/master/examples/demultiplex_input.json) for an example input file for the demultiplexing and pre-processing pipeline.

### Barcode and UMI regex strings
If the barcode and UMI sequence will be parsed from the FASTQ read headers (especially if using UMITools to pre-parse), the best (but not only way) to parse is by providing the regex strings for parsing, i.e.:

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

For a protospacer length of 20, surrogate length of 32, and barcode length of 4, we suggest the following thresholds:
```
PROTOSPACER_HAMMING_THRESHOLD = 7
SURROGATE_HAMMING_THRESHOLD = 10 
BARCODE_HAMMING_THRESHOLD = 2
```

**Tip**: These thresholds are based on an arbitrary maximum number of edited bases given a canonical editing window of a base-editor. As long as the numbers are not too low (i.e. 2-4 where highly-edited sequences are filtered out, or too high i.e. 13+ where random sequences get mapped), then the thresholds should be sufficient and don't need to be optimized.

## Running the guide mapping

### How is mapping performed?
For each sample, the tool will iterate through the R1 (and R2) sequences and extract the protospacer/surrogate/barcode/UMI sequences based on the provided parsing specification. To prevent redundant mapping of the same sequences, the tool will create a Python *Counter* object to store unique sequences as the key and their read counts as the value. For each unique observed sequence, the hamming distance for the protospacer/surrogate/barcode sequence will be calculated for all sequences in the provided guide library in an efficient vectorized manner. The sequence with the minimum hamming distance is the mapped sequence (unless the minimum hamming distance is above the thresholds provided). If there are multiple sequences with the same minimum hamming distance (ambiguous mapping), the mapping count is either assigned to all ambiguous mappings, spread across the ambiguous mappings, or counted to neither of the ambiguous mappings. The protospacer/surrogate/barcode sequences are mapped separately to diagnose any mismatching between protospacer and surrogate sequences (due to possible recombination effects). **The tool will not be effective at mapping sequences with in-del edits since it is based on base-per-base hamming distance.**

### Run the mapping function

Import the library "crispr_ambiguous_mapping" and pandas, load in the guide library TSV, and pass in the inputs arguments (along with the number of cores for parallelization).

Below is the complete function definition of the primary mapping function:

```
def get_whitelist_reporter_counts_from_fastq(whitelist_guide_reporter_df: pd.DataFrame, 
                                               fastq_r1_fn: str, 
                                               fastq_r2_fn: Optional[str], 
                                               
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
                                               cores: int=1)
```

- Provide your guide library table as `whitelist_guide_reporter_df`
- Provide your R1 (can be .gz or fastq.gz) as `fastq_r1_fn`, and your optional R2 as `fastq_r1_fn` if available.
- Provide the parsing specifications for each of your sequences: you can provide a regex string as `{sequence_type}_pattern_regex`, or you can provide the 0-based starting index of the sequence as `{sequence_type}_start_position` or a string to the left of the sequence as `{sequence_type}_left_flank` (i.e. "CACCG" for the protospacer_left_flank), then you can provide the 0-based end index of the sequence as `{sequence_type}_end_position` or a string to the right of the sequence as `{sequence_type}_right_flank` (i.e. "GTTTTA" for the protospacer_right_flank) or the size of the sequence as `{sequence_type}_length` (i.e. 20 for the protospacer_length). These specifications will depend on the guide design, sequencing strategy,and any pre-processing done, therefore you should look at your FASTQs by eye to determine the best approach.
- You will need to specify where the sequence resides and needs to be parsed from. If it is in the R1 sequence, set `is_{sequence_type}_r1` to True. If it is in the R1 header sequence, set `is_{sequence_type}_header` to True. If it is in the R2 sequence and the `fastq_r2_fn` is provided, set `is_{sequence_type}_r1` and `is_{sequence_type}_header` to False.
- If the sequence needs to be reverse complemented (i.e. if it is on the R2 read), then set `revcomp_{sequence_type}` to True, else set to False.
- Provide the thresholds for each sequence as `{sequence_type}_hamming_threshold_strict`
- Set `cores` to the number of CPUs on your system that you want to use for parallelization, i.e. if on an 8-core computer, you can set `cores` to 8 max.
  
Here is an example run (will likely be differen to your configuration):

```
import crispr_ambiguous_mapping
import pandas as pd

# Provide guide library as structure in section above
GUIDE_LIBRARY_DATAFRAME = pd.read_table("{GUIDE_LIBRARY_TSV}")

# Provide regex to parse barcode and UMI as specified in above section
barcode_pattern = r"_([^_ ]+)[\s+]"
umi_pattern = r":([^+:]{6})(.{2})\+"

result = crispr_ambiguous_mapping.mapping.get_whitelist_reporter_counts_from_fastq(
       whitelist_guide_reporter_df=whitelist_guide_reporter_df, 
       fastq_r1_fn=read1_fn, 
       fastq_r2_fn=read2_fn, 

       # Protospacer and surrogate parsing don't via positional indices 
       #protospacer_pattern_regex = None,
       #surrogate_pattern_regex = None,

       barcode_pattern_regex=barcode_pattern,
       umi_pattern_regex=umi_pattern,

       # Protospacer parsing done via first 20 bases 
       #protospacer_left_flank = None,
       #protospacer_right_flank = None,
       protospacer_start_position = 0,
       #protospacer_end_position = None,
       protospacer_length = 20,

       # Protospacer parsing done via first 32 bases 
       #surrogate_left_flank = None,
       #surrogate_right_flank = None,
       surrogate_start_position = 0,
       #surrogate_end_position = None,
       surrogate_length = 32,


       # Barcode parsing done via regex, below arguments not needed 
       #barcode_left_flank = None,
       #barcode_right_flank = None,
       #barcode_start_position = None,
       #barcode_end_position = None,
       #barcode_length = None,

       # UMI parsing done via regex, below arguments not needed
       #umi_left_flank = None,
       #umi_right_flank = None,
       #umi_start_position = None,
       #umi_end_position = None,
       #umi_length = None,

       # Protospacer parsed from R1, surrogate from R2, barcode and UMI from header
       is_protospacer_r1 = True, 
       is_surrogate_r1 = False, 
       is_barcode_r1 = False,
       is_umi_r1 = False,
       is_protospacer_header = False, 
       is_surrogate_header = False, 
       is_barcode_header = True,
       is_umi_header = True,

       # Both surrogate and barcode came from R2 read, so both will be reverse complemented 
       revcomp_protospacer = False, 
       revcomp_surrogate = True, 
       revcomp_barcode = True, 
       revcomp_umi = False,

       # Provide thresholds as specified in above section 
       surrogate_hamming_threshold_strict=10, 
       barcode_hamming_threshold_strict=2, 
       protospacer_hamming_threshold_strict=7,

       # Set cores to max number of CPUs on my computer 
       cores=8)
```

### See the mapping outputs

You can save the mapping outputs to your system as below:
```
sample_name = "SAMPLE"

match_set_whitelist_reporter_observed_sequence_counter_series_results = crispr_ambiguous_mapping.processing.get_matchset_alleleseries(result.observed_guide_reporter_umi_counts_inferred, "protospacer_match_surrogate_match_barcode_match", contains_surrogate=result.count_input.contains_surrogate, contains_barcode=result.count_input.contains_barcode, contains_umi=result.count_input.contains_umi) 
mutations_results = crispr_ambiguous_mapping.processing.get_mutation_profile(match_set_whitelist_reporter_observed_sequence_counter_series_results, whitelist_reporter_df=whitelist_guide_reporter_df, contains_surrogate=result.count_input.contains_surrogate, contains_barcode=result.count_input.contains_barcode) 
linked_mutation_counters = crispr_ambiguous_mapping.processing.tally_linked_mutation_count_per_sequence(mutations_results=mutations_results, contains_surrogate = result.count_input.contains_surrogate, contains_barcode = result.count_input.contains_barcode)
crispr_ambiguous_mapping.visualization.plot_mutation_count_histogram(linked_mutation_counters.protospacer_total_mutation_counter, filename=f"{sample_name}_protospacer_total_mutation_histogram.png")
crispr_ambiguous_mapping.visualization.plot_mutation_count_histogram(linked_mutation_counters.surrogate_total_mutation_counter, filename=f"{sample_name}_surrogate_total_mutation_histogram.png")
crispr_ambiguous_mapping.visualization.plot_mutation_count_histogram(linked_mutation_counters.barcode_total_mutation_counter, filename=f"{sample_name}_barcode_total_mutation_histogram.png")

with open(f"{sample_name}_protospacer_editing_efficiency.txt", "w") as text_file:
        print(crispr_ambiguous_mapping.utility.calculate_average_editing_frequency(linked_mutation_counters.protospacer_total_mutation_counter), file=text_file)
with open(f"{sample_name}_surrogate_editing_efficiency.txt", "w") as text_file:
        print(crispr_ambiguous_mapping.utility.calculate_average_editing_frequency(linked_mutation_counters.surrogate_total_mutation_counter), file=text_file)
with open(f"{sample_name}_barcode_editing_efficiency.txt", "w") as text_file:
        print(crispr_ambiguous_mapping.utility.calculate_average_editing_frequency(linked_mutation_counters.barcode_total_mutation_counter), file=text_file)

crispr_ambiguous_mapping.visualization.plot_trinucleotide_mutational_signature(mutations_results=mutations_results, count_attribute_name="ambiguous_accepted_umi_noncollapsed_mutations", unlinked_mutation_attribute_name = "all_observed_surrogate_unlinked_mutations_df", label=f'{sample_name}', filename=f"{sample_name}_surrogate_trinucleotide_mutational_signature.png")
crispr_ambiguous_mapping.visualization.plot_positional_mutational_signature(mutations_results=mutations_results, count_attribute_name="ambiguous_accepted_umi_noncollapsed_mutations", unlinked_mutation_attribute_name = "all_observed_surrogate_unlinked_mutations_df", label=f'{sample_name}', min_position = 6, max_position=20, filename=f"{sample_name}_surrogate_trinucleotide_positional_signature.png")

crispr_ambiguous_mapping.utility.save_or_load_pickle("./", f"{sample_name}_match_set_whitelist_reporter_observed_sequence_counter_series_results", py_object = match_set_whitelist_reporter_observed_sequence_counter_series_results, date_string="")
crispr_ambiguous_mapping.utility.save_or_load_pickle("./", f"{sample_name}_mutations_results", py_object = mutations_results, date_string="")
crispr_ambiguous_mapping.utility.save_or_load_pickle("./", f"{sample_name}_linked_mutation_counters", py_object = linked_mutation_counters, date_string="")
crispr_ambiguous_mapping.utility.save_or_load_pickle("./", f"{sample_name}_whitelist_guide_reporter_df", py_object = whitelist_guide_reporter_df, date_string="")


# Store the complete count result object. This will be a very large object
crispr_ambiguous_mapping.utility.save_or_load_pickle("./", f"{sample_name}_result", py_object = result, date_string="")

# Store the components of the result object, so that the user can load the information as needed
crispr_ambiguous_mapping.utility.save_or_load_pickle("./", f"{sample_name}_count_series_result", py_object = result.all_match_set_whitelist_reporter_counter_series_results, date_string="")
crispr_ambiguous_mapping.utility.save_or_load_pickle("./", f"{sample_name}_observed_guide_reporter_umi_counts_inferred", py_object = result.observed_guide_reporter_umi_counts_inferred, date_string="")
crispr_ambiguous_mapping.utility.save_or_load_pickle("./", f"{sample_name}_quality_control_result", py_object = result.quality_control_result, date_string="")
crispr_ambiguous_mapping.utility.save_or_load_pickle("./", f"{sample_name}_count_input", py_object = result.count_input, date_string="")
```

You can retrieve the count results directly at:
`result.all_match_set_whitelist_reporter_counter_series_results`


