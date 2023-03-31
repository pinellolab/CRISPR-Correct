import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from pandarallel import pandarallel
import matplotlib.pyplot as plt
from collections import Counter
import pickle
from datetime import date
from typeguard import typechecked
from os import listdir
from os.path import isfile, join
import re
from multiprocessing import Pool
from functools import partial
from itertools import repeat
import gzip
import random
from enum import Enum
from typing import Callable
from typing import Union, List, Mapping, Tuple, Optional


###
### UTILITY FUNCTIONS
###

def save_or_load_pickle(directory, label, py_object = None, date_string = None):
    '''Save a pickle for caching that is notated by the date'''
    
    if date_string == None:
        today = date.today()
        date_string = str(today.year) + ("0" + str(today.month) if today.month < 10 else str(today.month)) + str(today.day)
    
    filename = directory + label + "_" + date_string + '.pickle'
    print(filename)
    if py_object == None:
        with open(filename, 'rb') as handle:
            py_object = pickle.load(handle)
            return py_object
    else:
        with open(filename, 'wb') as handle:
            pickle.dump(py_object, handle, protocol=pickle.HIGHEST_PROTOCOL)

def display_all_pickle_versions(directory, label):
    '''Retrieve all pickles with a label, specifically to identify versions available'''
    return [f for f in listdir(directory) if isfile(join(directory, f)) and label == f[:len(label)]]


###
### GUIDE MAPPING HELPER FUNCTIONS
###

# Vector encoding of each base, since we want to vectorized
non_ambigious_encoding_dict = dict({
    "A": np.asarray([1,0,0,0]),
    "C": np.asarray([0,1,0,0]),
    "G": np.asarray([0,0,1,0]),
    "T": np.asarray([0,0,0,1]),
    "U": np.asarray([0,0,0,1])
})

'''
 Helper function for retrieving the encoding for ambigious bases based on IUPAC codes
'''
def ambiguity_encoder(bases):
    return np.logical_or.reduce([non_ambigious_encoding_dict[base] for base in bases]).astype(int)

'''
    Final dictionary for getting encoding of each IUPAC base
'''
full_encoding_dict = dict({
    "A": non_ambigious_encoding_dict["A"],
    "C": non_ambigious_encoding_dict["C"],
    "G": non_ambigious_encoding_dict["G"],
    "T": non_ambigious_encoding_dict["T"], 
    "R": ambiguity_encoder(["A", "G"]),
    "Y": ambiguity_encoder(["C", "T"]),
    "S": ambiguity_encoder(["G", "C"]),
    "W": ambiguity_encoder(["A", "T"]),
    "K": ambiguity_encoder(["G", "T"]),
    "M": ambiguity_encoder(["A", "C"]),
    "B": ambiguity_encoder(["C", "G", "T"]),
    "D": ambiguity_encoder(["A", "G", "T"]),
    "H": ambiguity_encoder(["A", "C", "T"]),
    "V": ambiguity_encoder(["A", "C", "G"]),
    "N": ambiguity_encoder(["A", "C", "G", "T"]),
})

'''
    Main function to encode a single base
'''
def encode_DNA_base(char):
    return full_encoding_dict[char]
encode_DNA_base_vectorized = np.vectorize(encode_DNA_base, signature='()->(n)') # Vectorized function for a string (i.e. gRNA)

'''
    Function for converting string (i.e. gRNA) into a np array of chars  - may be deprecated (NOTE 20221202)
'''
def numpify_string(string):
    return np.array(list(string), dtype=str)
numpify_string_vectorized = np.vectorize(numpify_string, signature='()->(n)') # Vectorize the function

def encode_guide_series(guide_series) -> np.array:
    guide_numpy = guide_series.to_numpy(dtype=object)
    guide_numpy = guide_numpy.astype(str)
    guide_numpy_char = np.array(list(map(list, guide_numpy))) # Map into a list of list of characters
    guide_numpy_encoding = encode_DNA_base_vectorized(guide_numpy_char)
    return guide_numpy_encoding


def retrieve_hamming_distance_whitelist(target_guide_encoded, whitelist_guide_encoded):
    '''
        This takes a encoded guide sequence and a list of encoded whitelisted guides and matrix computes the hamming distance of the 
        encoded guide across all whitelisted guides in a single operation
        
        (target_guide_encoded*whitelist_guide_encoded[:, np.newaxis]).sum(axis=3) # Determines 
        
    '''
    return ((target_guide_encoded*whitelist_guide_encoded[:, np.newaxis]).sum(axis=3)^1).sum(axis=2).flatten()

###
### GUIDE PARSING HELPER FUNCTIONS
###

'''
    Strategies to parse read
'''
@typechecked
def parse_read_positional(read_sequence: Union[str, Seq], position_start: int, position_end: int) -> Union[str, Seq]:  
    return read_sequence[position_start:position_end]

@typechecked
def parse_read_left_flank(read_sequence: Union[str, Seq], left_flank:Union[str, Seq], guide_sequence_length:int) -> Union[str, Seq]: 
    position_start = read_sequence.find(left_flank) + len(left_flank)
    return read_sequence[position_start:position_start+guide_sequence_length]

@typechecked
def parse_read_right_flank(read_sequence: Union[str, Seq], right_flank:Union[str, Seq], guide_sequence_length:int) -> Union[str, Seq]:
    position_end = read_sequence.find(right_flank)
    return read_sequence[position_end-guide_sequence_length:position_end]


'''
    Extract the guide sequence from the read provided
'''
@typechecked
def parse_guide_sequence(read_sequence: Union[str, Seq], parser_function: Callable) -> Union[str, Seq]:
    read_guide_sequence = parser_function(read_sequence)
    return read_guide_sequence


'''
    Iterate over all the reads in the FASTQ (parallelized) and retrieve the observed guide sequence
'''
@typechecked
def retrieve_fastq_guide_sequences(fastq_file: str, cores: int=1) -> Union[List[str], List[Seq]]:
    parse_read_left_flank_p = partial(parse_read_left_flank, left_flank="CACCG", guide_sequence_length=20)
    parse_guide_sequence_p = partial(parse_guide_sequence, parser_function=parse_read_left_flank_p)
    
    with gzip.open(fastq_file, "rt", encoding="utf-8") as handle, Pool(cores) as pool:
        fastq_guide_sequences = pool.map(
        parse_guide_sequence_p,
        (seq.seq for seq in SeqIO.parse(handle, 'fastq')),
        chunksize=2000,
        )

    return fastq_guide_sequences


###
### GUIDE MAPPING MAIN FUNCTIONS
###
class GuideCountError(Enum):
    NO_MATCH = "No guide found within hamming distance"
    MULTIPLE_MATCH = "Multiple exact matches found for guide (likely a truncated guide read assuming guide series is unique)"
    NO_GUIDE_WITH_SAME_LENGTH = "No whitelisted guides with same length as observed guide - maybe try enabling truncating whitelisted guides"
   
'''
    Given an observed, potentially self-edited guide (in 'row["observed_sequence"]'), try and find the true guide sequence from the whitelist ("guide_sequences_series") based on hamming distance
    # TODO (3/1/2023): This is the main bottleneck in terms of computational efficiency, see if it can be sped up
'''

@typechecked
def infer_true_guides(observed_guide_sequence: str, whitelist_guide_sequences_series: Union[List[str], pd.Series],
encoded_whitelist_guide_sequences_series: np.array, consider_truncated_sequences: bool = False, hamming_threshold: int = 3):
    #observed_guide_sequence = str(observed_guide_row["observed_sequence"])
    
    # If considering truncated guides, truncate all whitelist guides to the same size of the observed guide, else only consider whitelisted guides of the same length (likely all guides provided are same length of 20nt)
    if consider_truncated_sequences == True:
        whitelist_guide_sequences_series = whitelist_guide_sequences_series.apply(lambda guide: guide[0:len(observed_guide_sequence)])
    else:
        whitelist_guide_sequences_series = whitelist_guide_sequences_series[whitelist_guide_sequences_series.str.len() == len(observed_guide_sequence)]
        if len(whitelist_guide_sequences_series) == 0:
            return GuideCountError.NO_GUIDE_WITH_SAME_LENGTH 
        

    # Determine if there are exact matches, hopefully just a single match
    whitelist_guide_sequences_series_match = whitelist_guide_sequences_series[whitelist_guide_sequences_series == observed_guide_sequence]
    
    # If there is a single exact match, great, no need for fancy mat
    if len(whitelist_guide_sequences_series_match) == 1: # Exact match found, return
        return whitelist_guide_sequences_series_match.index[0]
    
    # If no matches, possible a self-edit or a sequencing error, try and find guide with lowest hamming distance
    elif len(whitelist_guide_sequences_series_match) == 0: # No match found, search based on hamming distance
        
        # Encode the whitelisted guides 
        # NOTE 20221202: Potentially improve efficiency by passing in the encoded guide series (assuming no truncation) so that this does not have to be ran on every guide
        #guide_sequences_series_encoded = encode_guide_series(whitelist_guide_sequences_series_match)
        
        # Encode the observed guide
        observed_guide_sequence_encoded = encode_DNA_base_vectorized(numpify_string_vectorized(observed_guide_sequence)) 
        
        # Calculate the hamming distance of the guide with all whitelisted guides - vectorized operation
        observed_guide_sequence_dists = retrieve_hamming_distance_whitelist(observed_guide_sequence_encoded, encoded_whitelist_guide_sequences_series)
        
        # Get the minimum hamming distance calculated
        hamming_min = observed_guide_sequence_dists.min()
        
        # Get all whitelisted guides with the minimum hamming distance (could be multiple)
        guides_with_hamming_min = whitelist_guide_sequences_series[np.where(observed_guide_sequence_dists == hamming_min)[0]]
        
        # If the minimum hamming distance is greater than the specified threshold, then the guide is too ambigious to assign, so no match.
        if hamming_min >= hamming_threshold:
            return GuideCountError.NO_MATCH
        
        # If there are multiple guides with the minimum hamming distance, then the guide is ambigious, so no mapping (due to multiple match)
        # NOTE (3/30/2023): This could potentially be an issue with igRNAs maybe?
        elif len(guides_with_hamming_min) > 1:
            return GuideCountError.MULTIPLE_MATCH
        
        # Else if there is 1 guide with the match, then return the match
        else:
            return guides_with_hamming_min.index[0]
    
    # Else if there are multiple exact match, which should never occur unless the whitelisted guide list is not unique, then return multiple match.
    else:
        return GuideCountError.MULTIPLE_MATCH
        #raise Exception("Multiple exact matches of the provided whitelisted guides - there are likely duplicates in the provided whitelist, please remove. Observed guide={}, guide matches={}".format(observed_guide_sequence, guide_sequences_series_match)) # NOTE 12/6/22: REMOVED THIS EXCEPTION - another reason is from truncated guides having multiple matches. In production code, just make sure to ensure that the whitelist is the set.

'''
    This performs simulation of determining how many mutations it takes for a guide to be ambigiously mapped based on hamming distance.

    This is performed on *sample_count* number of randomly selected guides, and the distribution of mutation count is determined. The 5% quantile is used as the threshold.
    
    This is useful in determing the ideal hamming distance threshold specific to a guide library
'''
@typechecked
def determine_hamming_threshold(whitelist_guide_sequences_series: Union[List[str],pd.Series], encoded_whitelist_guide_sequences_series: np.array, sample_count: int = 100, quantile: float = 0.05) -> float:
    #encoded_whitelist_guide_sequences_series = encode_guide_series(guide_sequences_series)
    
    mutation_count_until_nonunique = []
    for i in range(sample_count):
        # Sample a guide from whitelist
        sampled_guide = whitelist_guide_sequences_series.sample()[0]
        
        # Generate order of base positions to mutate sequentially.
        guide_position_list = list(range(len(sampled_guide)))
        random.shuffle(guide_position_list)
        
        # Create temporary variable to represent mutated guide
        current_guide_sequence = sampled_guide
        
        # Iterate through positions to mutate
        for iteration, position in enumerate(guide_position_list):
            
            # Mutate the guide
            current_guide_sequence_separated = list(current_guide_sequence)
            guide_position_nt = current_guide_sequence_separated[position]
            nt_list = ["A", "C", "T", "G"]
            nt_list.remove(guide_position_nt.upper())
            new_nt = random.sample(nt_list, 1)[0]
            current_guide_sequence_separated[position] = new_nt
            current_guide_sequence = "".join(current_guide_sequence_separated)
            current_guide_sequence_encoded = encode_DNA_base_vectorized(numpify_string_vectorized(current_guide_sequence)) 
            
            # Calculate the hamming distance of the mutated guide to all other guides
            hamming_distances =  retrieve_hamming_distance_whitelist(current_guide_sequence_encoded, encoded_whitelist_guide_sequences_series)
            if len(np.where(hamming_distances == hamming_distances.min())[0]) > 1: # TODO (3/30/23): Can also add to the conditional whether the minimum hamming distance guide is still the original guide - but it is probably rare cases where this would not be the case, so not too important to implement
                mutation_count_until_nonunique.append(iteration+1)
                break   
    mutation_count_until_nonunique = pd.Series(mutation_count_until_nonunique)
    
    # From the mutation count distribution, calculate the threshold based on the provided quantile.
    hamming_distance_threshold = mutation_count_until_nonunique.quantile(quantile)
    return hamming_distance_threshold


'''
    Take in input FASTQ filename, and a set of whitelisted guide sequences
'''
@typechecked
def get_guide_counts(whitelist_guide_sequences_series: pd.Series, fastq_fn: str, hamming_threshold_strict: int = 3, hamming_threshold_dynamic: bool = False, cores: int=1):
    # Retrieve all observed guide sequences
    print("Retrieving FASTQ guide sequences and counting: " + fastq_fn)
    observed_guide_raw_sequences = retrieve_fastq_guide_sequences(fastq_fn, cores=cores)
    
    # Count each unique observed guide sequence
    observed_guide_sequences_counts = Counter(observed_guide_raw_sequences)
    
    # Create observed guide DF to contain all information
    # The "observed sequence" column represents all *unique* observed sequences that may have self-edits/errors. The "observed counts" column represents the count of each observed count.
    observed_guides_df = pd.DataFrame({"observed_sequence":[str(sequence) for sequence in observed_guide_sequences_counts.keys()], "observed_counts":observed_guide_sequences_counts.values()})
 
    # Get the hamming distance threshold. THe hamming distance must be below this threshold to assign an observed guide to a whitelist guide.
    encoded_whitelist_guide_sequences_series = encode_guide_series(whitelist_guide_sequences_series)
    if hamming_threshold_dynamic:
        hamming_threshold = int(determine_hamming_threshold(whitelist_guide_sequences_series, encoded_whitelist_guide_sequences_series, sample_count = 100, quantile = 0.05))
        print("Hamming threshold is " + str(hamming_threshold))
    else:
        hamming_threshold = hamming_threshold_strict
        
    # Infer whitelist guides from observed guides
    print("Inferring the true guides from observed guides")

    inferred_true_guide_sequences = None
    with Pool(cores) as pool:
        inferred_true_guide_sequences = pool.map(
        infer_true_guides,
        observed_guides_df["observed_sequence"],
        chunksize=len(observed_guides_df["observed_sequence"])/cores,
        )
        
    
    observed_guides_df["inferred_guides"] = inferred_true_guide_sequences
    
    '''
        QC
    '''
    # QC: Determine the number of guides passed by determining that the result is not an error
    observed_guides_df_passed_inference = observed_guides_df[observed_guides_df["inferred_guides"].apply(lambda guide : type(guide) != GuideCountError)]
    # QC: Calculate number of guides that were unassigned
    observed_guides_df_no_match = observed_guides_df[observed_guides_df["inferred_guides"].apply(lambda guide : guide == GuideCountError.NO_MATCH)]
    # QC: Calculate number of guides with multiple inferred guides
    observed_guides_df_multiple_match = observed_guides_df[observed_guides_df["inferred_guides"].apply(lambda guide : guide == GuideCountError.MULTIPLE_MATCH)]
    # QC: Calculate percent mapped
    percent_mapped = observed_guides_df_passed_inference["observed_counts"].sum()/observed_guides_df["observed_counts"].sum()
    
    # Retrieve the observed sequences that were mapped and set the inferred guides
    inferred_guide_sequence_counter = Counter()
    for _, row in observed_guides_df_passed_inference.iterrows():
        inferred_guide_sequence_counter[row["inferred_guides"]] += row["observed_counts"]
    
    whitelist_guide_sequences_series_counts = whitelist_guide_sequences_series.apply(lambda guide: inferred_guide_sequence_counter[guide])
    whitelist_guide_sequences_series_counts.index = whitelist_guide_sequences_series
    
    qc_dict = {"guide_sequences_unassigned_counts":observed_guides_df_no_match.sum(), "guide_sequences_multiple_counts": observed_guides_df_multiple_match.sum(), "total_guide_counts": observed_guides_df["observed_counts"].sum(), "percent_mapped": percent_mapped}
    
    return observed_guide_raw_sequences, whitelist_guide_sequences_series_counts, observed_guides_df, qc_dict

