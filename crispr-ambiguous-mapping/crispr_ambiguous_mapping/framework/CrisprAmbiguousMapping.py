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
from typing import Union, List, Mapping, Tuple, Optional, Any


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
def retrieve_fastq_guide_sequences(fastq_file: str, parse_left_flank: bool = True, parse_flank_sequence: Union[None, str] = None, cores: int=1) -> Union[List[str], List[Seq]]:
    parse_guide_sequence_p = None
    if parse_left_flank:
        if parse_flank_sequence is None:
            print("No flank sequence passed. Setting left-flank default sequence to CACCG assuming U6 G+N20 guide")
            flank_sequence = "CACCG"
        else:
            flank_sequence = parse_flank_sequence

        parse_read_left_flank_p = partial(parse_read_left_flank, left_flank=flank_sequence, guide_sequence_length=20)
        parse_guide_sequence_p = partial(parse_guide_sequence, parser_function=parse_read_left_flank_p)
    else:
        if parse_flank_sequence is None:
            print("No flank sequence passed. Setting right-flank default sequence to GTTTT. If you are using sgRNA(F+E) design, flank sequence may need to be changed")
            flank_sequence = "GTTTT"
        else:
            flank_sequence = parse_flank_sequence
        
        parse_read_right_flank_p = partial(parse_read_right_flank, right_flank=flank_sequence, guide_sequence_length=20)
        parse_guide_sequence_p = partial(parse_guide_sequence, parser_function=parse_read_right_flank_p)

    if fastq_file.endswith('.gz'):
        with gzip.open(fastq_file, "rt", encoding="utf-8") as handle, Pool(cores) as pool:
            fastq_guide_sequences = pool.map(
                parse_guide_sequence_p,
                (seq.seq for seq in SeqIO.parse(handle, 'fastq')),
                chunksize=2000,
            )
    else:
        with Pool(cores) as pool:
            fastq_guide_sequences = pool.map(
                parse_guide_sequence_p,
                (seq.seq for seq in SeqIO.parse(fastq_file, 'fastq')),
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
encoded_whitelist_guide_sequences_series, consider_truncated_sequences: bool = False, hamming_threshold: int = 3, verbose_result = False):
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
            if verbose_result:
                return {"Error": GuideCountError.NO_MATCH, "hamming_min": hamming_min}
            else:
                return GuideCountError.NO_MATCH
        
        # If there are multiple guides with the minimum hamming distance, then the guide is ambigious, so no mapping (due to multiple match)
        # NOTE (3/30/2023): This could potentially be an issue with igRNAs maybe?
        elif len(guides_with_hamming_min) > 1:
            if verbose_result:
                return {"Error": GuideCountError.MULTIPLE_MATCH,
                "exact_match": False, 
                "num_matches": len(guides_with_hamming_min),
                "matches": guides_with_hamming_min,
                "hamming_min": hamming_min}
            else:
                return GuideCountError.MULTIPLE_MATCH
        
        # Else if there is 1 guide with the match, then return the match
        else:
            return guides_with_hamming_min.index[0]
    
    # Else if there are multiple exact match, which should never occur unless the whitelisted guide list is not unique, then return multiple match.
    else:
        
        if verbose_result:
            return {"Error": GuideCountError.MULTIPLE_MATCH,
            "exact_match": True, 
            "num_matches": len(whitelist_guide_sequences_series_match),
            "matches": whitelist_guide_sequences_series_match,
            "hamming_min": 0}
        else:
            return GuideCountError.MULTIPLE_MATCH
        #raise Exception("Multiple exact matches of the provided whitelisted guides - there are likely duplicates in the provided whitelist, please remove. Observed guide={}, guide matches={}".format(observed_guide_sequence, guide_sequences_series_match)) # NOTE 12/6/22: REMOVED THIS EXCEPTION - another reason is from truncated guides having multiple matches. In production code, just make sure to ensure that the whitelist is the set.

@typechecked
def infer_true_reporters(observed_reporter_sequences, whitelist_guide_reporter_df: pd.DataFrame,
encoded_whitelist_guide_sequences_series, encoded_whitelist_barcodes_series, barcode_hamming_threshold: int = 2, hamming_threshold: int = 4, verbose_result = False):
    for element in observed_reporter_sequences:
        if (element is None) or (element == "None") or (element == ""):
            if verbose_result:
                return {"Error": GuideCountError.NO_MATCH, "message": "Observed protospacer/surrogate/barcode is None"}
            else:
                return GuideCountError.NO_MATCH
    
    observed_reporter_sequences = pd.Series(observed_reporter_sequences, index=["protospacer", "surrogate", "barcode"])
    
    # Check if the observed protospacer length is the same length of the whitelisted guides
    if len(observed_reporter_sequences["protospacer"]) != encoded_whitelist_guide_sequences_series.shape[1]: 
        if verbose_result:
            return {"Error": GuideCountError.NO_MATCH, "message": f"Observed protospacer {observed_reporter_sequences['protospacer']} not of correct length: {encoded_whitelist_guide_sequences_series.shape[1]}"}
        else:
            return GuideCountError.NO_MATCH
    if len(observed_reporter_sequences["barcode"]) != encoded_whitelist_barcodes_series.shape[1]: 
        if verbose_result:
            return {"Error": GuideCountError.NO_MATCH, "message": f"Observed barcode {observed_reporter_sequences['barcode']} not of correct length: {encoded_whitelist_barcodes_series.shape[1]}"}
        else:
            return GuideCountError.NO_MATCH   
            
    # Determine if there are exact matches, hopefully just a single match
    whitelist_guide_reporter_df_match = whitelist_guide_reporter_df[(whitelist_guide_reporter_df["protospacer"]==observed_reporter_sequences["protospacer"]) & (whitelist_guide_reporter_df["surrogate"]==observed_reporter_sequences["surrogate"]) & (whitelist_guide_reporter_df["barcode"]==observed_reporter_sequences["barcode"])] #whitelist_guide_sequences_series == observed_guide_sequence
    
    # If there is a single exact match, great, no need for fancy mat
    if whitelist_guide_reporter_df_match.shape[0] == 1: # Exact match found, return
        return tuple(whitelist_guide_reporter_df_match.iloc[0])
    
    # If no matches, possible a self-edit or a sequencing error, try and find guide with lowest hamming distance
    elif whitelist_guide_reporter_df_match.shape[0] == 0: # No match found, search based on hamming distance
        
        
        ###
        ### FILTER BY BARCODE (NOTE: 4/2/23: This is the main difference with the traditional filtering by protospacer)
        ### TODO 4/2/23: Assumes perfect barcode match, but in the future I can also select for barcodes of 1 hamming - 2 hamming if there is no matches with 0 hamming
        ###
        observed_barcode_encoded = encode_DNA_base_vectorized(numpify_string_vectorized(observed_reporter_sequences["barcode"])) 
        observed_barcode_dists = retrieve_hamming_distance_whitelist(observed_barcode_encoded, encoded_whitelist_barcodes_series)
        barcode_hamming_min = observed_barcode_dists.min()
        if barcode_hamming_min >= barcode_hamming_threshold:
            if verbose_result:
                return {"Error": GuideCountError.NO_MATCH, "barcode_hamming_min": barcode_hamming_min, "message": f"No barcode below threshold {barcode_hamming_threshold}"}
            else:
                return GuideCountError.NO_MATCH
            
        barcode_matches_indices = np.where(observed_barcode_dists == barcode_hamming_min)[0]
        
        whitelist_guide_reporter_df_barcode = whitelist_guide_reporter_df.iloc[barcode_matches_indices]
        encoded_whitelist_guide_sequences_series_barcode = encoded_whitelist_guide_sequences_series[barcode_matches_indices]
        ###
        ### Now do similar as before
        ###

        
        # Encode the observed guide
        observed_guide_sequence_encoded = encode_DNA_base_vectorized(numpify_string_vectorized(observed_reporter_sequences["protospacer"])) 
        
        # Calculate the hamming distance of the guide with all whitelisted guides - vectorized operation
        observed_guide_sequence_dists = retrieve_hamming_distance_whitelist(observed_guide_sequence_encoded, encoded_whitelist_guide_sequences_series_barcode)
        
        # Get the minimum hamming distance calculated
        hamming_min = observed_guide_sequence_dists.min()
        
        # Get all whitelisted guides with the minimum hamming distance (could be multiple)
        reporters_with_hamming_min_df = whitelist_guide_reporter_df_barcode.iloc[np.where(observed_guide_sequence_dists == hamming_min)[0]]
        
        # If the minimum hamming distance is greater than the specified threshold, then the guide is too ambigious to assign, so no match.
        if hamming_min >= hamming_threshold:
            if verbose_result:
                return {"Error": GuideCountError.NO_MATCH, "hamming_min": hamming_min}
            else:
                return GuideCountError.NO_MATCH
        
        # If there are multiple guides with the minimum hamming distance, then the guide is ambigious, so no mapping (due to multiple match)
        # NOTE (3/30/2023): This could potentially be an issue with igRNAs maybe?
        elif reporters_with_hamming_min_df.shape[0] > 1:
            if verbose_result:
                return {"Error": GuideCountError.MULTIPLE_MATCH,
                "exact_match": False, 
                "num_matches": reporters_with_hamming_min_df.shape[0],
                "matches": reporters_with_hamming_min_df,
                "hamming_min": hamming_min}
            else:
                return GuideCountError.MULTIPLE_MATCH
        
        # Else if there is 1 guide with the match, then return the match
        else:
            return tuple(reporters_with_hamming_min_df.iloc[0])
    
    # Else if there are multiple exact match, which should never occur unless the whitelisted guide list is not unique, then return multiple match.
    else:
        
        if verbose_result:
            return {"Error": GuideCountError.MULTIPLE_MATCH,
            "exact_match": True, 
            "num_matches": whitelist_guide_reporter_df_match.shape[0],
            "matches": whitelist_guide_reporter_df_match,
            "hamming_min": 0}
        else:
            return GuideCountError.MULTIPLE_MATCH
        #raise Exception("Multiple exact matches of the provided whitelisted guides - there are likely duplicates in the provided whitelist, please remove. Observed guide={}, guide matches={}".format(observed_guide_sequence, guide_sequences_series_match)) # NOTE 12/6/22: REMOVED THIS EXCEPTION - another reason is from truncated guides having multiple matches. In production code, just make sure to ensure that the whitelist is the set.

     
'''
    This performs simulation of determining how many mutations it takes for a guide to be ambigiously mapped based on hamming distance.

    This is performed on *sample_count* number of randomly selected guides, and the distribution of mutation count is determined. The 5% quantile is used as the threshold.
    
    This is useful in determing the ideal hamming distance threshold specific to a guide library
'''
@typechecked
def determine_hamming_threshold(whitelist_guide_sequences_series: Union[List[str],pd.Series], encoded_whitelist_guide_sequences_series, sample_count: int = 100, quantile: float = 0.05) -> float:
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
def get_guide_counts(observed_guide_sequences_counts: Counter, whitelist_guide_sequences_series: pd.Series, hamming_threshold_strict: int = 3, hamming_threshold_dynamic: bool = False, cores: int=1):

    whitelist_guide_sequences_series.index = whitelist_guide_sequences_series.values
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

    infer_true_guides_p = partial(infer_true_guides,
            whitelist_guide_sequences_series=whitelist_guide_sequences_series,
            encoded_whitelist_guide_sequences_series=encoded_whitelist_guide_sequences_series,
            hamming_threshold=hamming_threshold, verbose_result=False)

    inferred_true_guide_sequences = None
    with Pool(cores) as pool:
        inferred_true_guide_sequences = pool.map(
        infer_true_guides_p,
        observed_guides_df["observed_sequence"]
        )
        
    print("Completed inference")

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
    
    qc_dict = {"guide_sequences_unassigned_counts":observed_guides_df_no_match["observed_counts"].sum(), "guide_sequences_multiple_counts": observed_guides_df_multiple_match["observed_counts"].sum(), "total_guide_counts": observed_guides_df["observed_counts"].sum(), "percent_mapped": percent_mapped}
    
    return whitelist_guide_sequences_series_counts, observed_guides_df, qc_dict

@typechecked
def get_reporter_counts(observed_guide_reporters_counts: Counter, whitelist_guide_reporter_df: pd.DataFrame, barcode_hamming_threshold_strict: int = 2, hamming_threshold_strict: int = 4, hamming_threshold_dynamic: bool = False, cores: int=1):

    whitelist_guide_reporter_df.index = whitelist_guide_reporter_df["protospacer"]
    
    # Create observed guide DF to contain all information
    # The "observed sequence" column represents all *unique* observed sequences that may have self-edits/errors. The "observed counts" column represents the count of each observed count.
    observed_guides_df = pd.DataFrame({"observed_sequence":observed_guide_reporters_counts.keys(), "observed_counts":observed_guide_reporters_counts.values()})
 
    # Get the hamming distance threshold. THe hamming distance must be below this threshold to assign an observed guide to a whitelist guide.
    encoded_whitelist_guide_sequences_series = encode_guide_series(whitelist_guide_reporter_df["protospacer"])
    encoded_whitelist_barcodes_series = encode_guide_series(whitelist_guide_reporter_df["barcode"])
    
    if hamming_threshold_dynamic:
        hamming_threshold = int(determine_hamming_threshold(whitelist_guide_reporter_df["protospacer"], encoded_whitelist_guide_sequences_series, sample_count = 100, quantile = 0.05))
        print("Hamming threshold is " + str(hamming_threshold))
    else:
        hamming_threshold = hamming_threshold_strict
        
    # Infer whitelist guides from observed guides
    print("Inferring the true guides from observed guides")

    infer_true_reporters_p = partial(infer_true_reporters,
            whitelist_guide_reporter_df=whitelist_guide_reporter_df,
            encoded_whitelist_guide_sequences_series=encoded_whitelist_guide_sequences_series,
            encoded_whitelist_barcodes_series=encoded_whitelist_barcodes_series,
            hamming_threshold=hamming_threshold, barcode_hamming_threshold=barcode_hamming_threshold_strict, verbose_result=False)

    inferred_true_reporter_sequences = None
    if cores > 1:
        with Pool(cores) as pool:
            inferred_true_reporter_sequences = pool.map(
            infer_true_reporters_p,
            observed_guides_df["observed_sequence"]
           )
    else:
        inferred_true_reporter_sequences = [infer_true_reporters_p(obs_reporter) for obs_reporter in observed_guides_df["observed_sequence"]]
    
    print("Completed inference")

    observed_guides_df["inferred_guides"] = inferred_true_reporter_sequences
    
    '''
        QC
    '''
    print("Retrieving QC tables")
    # QC: Determine the number of guides passed by determining that the result is not an error
    observed_guides_df_passed_inference = observed_guides_df[observed_guides_df["inferred_guides"].apply(lambda guide : type(guide) != GuideCountError)]
    # QC: Calculate number of guides that were unassigned
    observed_guides_df_no_match = observed_guides_df[observed_guides_df["inferred_guides"].apply(lambda guide : guide == GuideCountError.NO_MATCH)]
    # QC: Calculate number of guides with multiple inferred guides
    observed_guides_df_multiple_match = observed_guides_df[observed_guides_df["inferred_guides"].apply(lambda guide : guide == GuideCountError.MULTIPLE_MATCH)]
    # QC: Calculate percent mapped
    percent_mapped = observed_guides_df_passed_inference["observed_counts"].sum()/observed_guides_df["observed_counts"].sum()
    
    print("Get whitelist reporter counts")
    # Retrieve the observed sequences that were mapped and set the inferred guides
    inferred_guide_sequence_counter = observed_guides_df_passed_inference.groupby("inferred_guides")["observed_counts"].sum()

    whitelist_guide_reporter_counts = whitelist_guide_reporter_df.apply(lambda reporter: inferred_guide_sequence_counter.get(tuple(reporter), 0), axis=1)
    
    
    multi_index = pd.MultiIndex.from_arrays([whitelist_guide_reporter_df['protospacer'], whitelist_guide_reporter_df['surrogate'], whitelist_guide_reporter_df['barcode']], names=['protospacer', 'surrogate', 'barcode'])
    whitelist_guide_reporter_counts.index = multi_index
    

    qc_dict = {"guide_sequences_unassigned_counts":observed_guides_df_no_match["observed_counts"].sum(), "guide_sequences_multiple_counts": observed_guides_df_multiple_match["observed_counts"].sum(), "total_guide_counts": observed_guides_df["observed_counts"].sum(), "percent_mapped": percent_mapped}
    
    return whitelist_guide_reporter_counts, observed_guides_df, qc_dict



###
### Read in sequences from a pre-parsed reporter file
### 
import csv
import multiprocessing as mp
from collections import Counter
from functools import reduce

def gen_chunks(reader, chunksize=1000):
    """
    Chunk generator. Take a CSV `reader` and yield
    `chunksize` sized slices.
    Source: https://gist.github.com/miku/820490
    """
    chunk = []
    for index, line in enumerate(reader):
        if index % chunksize == 0 and index > 0:
            yield chunk
            chunk = []
        chunk.append(line)
    yield chunk

def process(tsv_chunk, include_surrogate = False):
    local_counter = Counter()
    
    for row in tsv_chunk.copy():
        if include_surrogate:
            local_counter[(row[1], row[2], row[3])] += 1
        else:
            local_counter[row[1]] += 1
    
    return local_counter


@typechecked
def map_sample_protospacers(parsing_demult_handler, include_surrogate = False, cores=1):
    tsv_reader = csv.reader(parsing_demult_handler, delimiter='\t')  # change delimiter for normal csv files
    header = next(tsv_reader)
      
    read_chunks = gen_chunks(tsv_reader, chunksize=10000)

    process_func = partial(process, include_surrogate=include_surrogate)

    combined_counter = None
    if cores > 1:
        with mp.Pool(processes=cores) as pool:
            chunk_counters = pool.map(process_func, read_chunks)
            combined_counter = reduce(lambda x, y: x + y, chunk_counters)
    else:
        chunk_counters = [process_func(read_chunk) for read_chunk in read_chunks]
        combined_counter = reduce(lambda x, y: x + y, chunk_counters)

    return combined_counter 


###
### MAIN FUNCTIONS
###

'''
    Take in input FASTQ filename, and a set of whitelisted guide sequences
'''
@typechecked
def get_guide_counts_from_fastq(whitelist_guide_sequences_series: pd.Series, fastq_fn: str, hamming_threshold_strict: int = 3, hamming_threshold_dynamic: bool = False, parse_left_flank: bool = True, parse_flank_sequence: Union[None, str] = None, cores: int=1):
    # Retrieve all observed guide sequences
    print("Retrieving FASTQ guide sequences and counting: " + fastq_fn)
    observed_guide_raw_sequences = retrieve_fastq_guide_sequences(fastq_fn, parse_left_flank=parse_left_flank, parse_flank_sequence=parse_flank_sequence, cores=cores)
    
    # Count each unique observed guide sequence
    observed_guide_sequences_counts = Counter(observed_guide_raw_sequences)
    
    return get_guide_counts(observed_guide_sequences_counts, whitelist_guide_sequences_series, hamming_threshold_strict, hamming_threshold_dynamic, cores)


'''
    Take in input FASTQ filename, and a set of whitelisted guide sequences
'''
@typechecked
def get_guide_counts_from_reporter_tsv(whitelist_guide_sequences_series: pd.Series, reporter_tsv_fn: str, hamming_threshold_strict: int = 3, hamming_threshold_dynamic: bool = False, cores: int=1):
    combined_counter = None
    with open(reporter_tsv_fn, "r", newline='') as reporter_tsv_handler:
        combined_counter = map_sample_protospacers(reporter_tsv_handler, cores)

    return get_guide_counts(combined_counter, whitelist_guide_sequences_series, hamming_threshold_strict, hamming_threshold_dynamic, cores)

@typechecked
def get_reporter_counts_from_reporter_tsv(whitelist_guide_reporter_df: pd.DataFrame, reporter_tsv_fn: str, barcode_hamming_threshold_strict: int = 2, hamming_threshold_strict: int = 3, hamming_threshold_dynamic: bool = False, cores: int=1):
    
    combined_counter = None
    with open(reporter_tsv_fn, "r", newline='') as reporter_tsv_handler:
        combined_counter = map_sample_protospacers(reporter_tsv_handler, include_surrogate = True, cores=cores)

    return get_reporter_counts(observed_guide_reporters_counts=combined_counter, whitelist_guide_reporter_df=whitelist_guide_reporter_df, barcode_hamming_threshold_strict=barcode_hamming_threshold_strict, hamming_threshold_strict=hamming_threshold_strict, hamming_threshold_dynamic=hamming_threshold_dynamic, cores=cores)