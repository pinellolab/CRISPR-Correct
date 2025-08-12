import numpy as np

# Vector encoding of each base, since we want to vectorized.
# This is for encoding the whitelist library
non_ambigious_encoding_dict_whitelist = dict({
    "A": np.asarray([1,0,0,0]),
    "C": np.asarray([0,1,0,0]),
    "G": np.asarray([0,0,1,0]),
    "T": np.asarray([0,0,0,1]),
    "U": np.asarray([0,0,0,1]),
    "-": np.asarray([0,0,0,0]), # NOTE: This representes a deletion token
    "X": np.asarray([0,0,0,0]) # NOTE: This represents a padding token for sequences. Set to some high number, so that the calculated hamming distance is very high and won't be selected (which will happen if the observed sequence is greater than the whitelist sequence) # Correct 20250106 should be set to 0, so that it will just be a hamming increment of 1)
})

# Encoding of observed sequences
# Set mismatches to 3 so that after bitflip, hamming penalty is 2
non_ambigious_encoding_dict_observed = dict({
    "A": np.asarray([1,0,0,0]),
    "C": np.asarray([0,1,0,0]),
    "G": np.asarray([0,0,1,0]),
    "T": np.asarray([0,0,0,1]),
    "U": np.asarray([0,0,0,1]),
    "-": np.asarray([0,0,0,0]), # NOTE: This representes a deletion token
    "X": np.asarray([0,0,0,0]) # NOTE: This represents a padding token for sequences. Set to some high number, so that the calculated hamming distance is very high and won't be selected (which will happen if the observed sequence is greater than the whitelist sequence) # Correct 20250106 should be set to 0, so that it will just be a hamming increment of 1)
})

'''
 Helper function for retrieving the encoding for ambigious bases based on IUPAC codes
'''
def ambiguity_encoder(bases):
    return np.logical_or.reduce([non_ambigious_encoding_dict_whitelist[base] for base in bases]).astype(int)

'''
    Final dictionary for getting encoding of each IUPAC base
'''
# Final encoding of whitelist
full_encoding_dict_whitelist = dict({
    "A": non_ambigious_encoding_dict_whitelist["A"],
    "C": non_ambigious_encoding_dict_whitelist["C"],
    "G": non_ambigious_encoding_dict_whitelist["G"],
    "T": non_ambigious_encoding_dict_whitelist["T"], 
    "X": non_ambigious_encoding_dict_whitelist["X"], # NOTE: This represents a padding token for sequences
    "-": non_ambigious_encoding_dict_whitelist["-"], # NOTE: This represents a padding token for sequences 
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

full_encoding_dict_observed = dict({
    "A": non_ambigious_encoding_dict_observed["A"],
    "C": non_ambigious_encoding_dict_observed["C"],
    "G": non_ambigious_encoding_dict_observed["G"],
    "T": non_ambigious_encoding_dict_observed["T"], 
    "X": non_ambigious_encoding_dict_observed["X"], # NOTE: This represents a padding token for sequences
    "-": non_ambigious_encoding_dict_observed["-"], # NOTE: This represents a padding token for sequences 
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
def encode_DNA_base_whitelist(char):
    return full_encoding_dict_whitelist[char]
encode_DNA_base_whitelist_vectorized = np.vectorize(encode_DNA_base_whitelist, signature='()->(n)') # Vectorized function for a string (i.e. gRNA)

def encode_DNA_base_observed(char):
    return full_encoding_dict_observed[char]
encode_DNA_base_observed_vectorized = np.vectorize(encode_DNA_base_observed, signature='()->(n)') # Vectorized function for a string (i.e. gRNA)

'''
    Function for converting string (i.e. gRNA) into a np array of chars  - may be deprecated (NOTE 20221202)
'''
def numpify_string(string):
    return np.array(list(string), dtype=str)
numpify_string_vectorized = np.vectorize(numpify_string, signature='()->(n)') # Vectorize the function

def encode_guide_series_whitelist(guide_series) -> np.array:
    guide_numpy = guide_series.to_numpy(dtype=object)
    guide_numpy = guide_numpy.astype(str)
    guide_numpy_char = np.array(list(map(list, guide_numpy))) # Map into a list of list of characters
    guide_numpy_encoding = encode_DNA_base_whitelist_vectorized(guide_numpy_char)
    return guide_numpy_encoding

def encode_guide_series_observed(guide_series) -> np.array:
    guide_numpy = guide_series.to_numpy(dtype=object)
    guide_numpy = guide_numpy.astype(str)
    guide_numpy_char = np.array(list(map(list, guide_numpy))) # Map into a list of list of characters
    guide_numpy_encoding = encode_DNA_base_observed_vectorized(guide_numpy_char)
    return guide_numpy_encoding

def retrieve_hamming_distance_whitelist(observed_guide_encoded, whitelist_guide_encoded):
    '''
        This takes a encoded guide sequence and a list of encoded whitelisted guides and matrix computes the hamming distance of the 
        encoded guide across all whitelisted guides in a single operation.
        
        
        (observed_guide_encoded*whitelist_guide_encoded[:, np.newaxis]).clip(0, 2).sum(axis=3) # Determines 

        # observed_guide_encoded*whitelist_guide_encoded[:, np.newaxis] performs multiplication of the target encoding with every whitelist encoding
        #
        
    '''
    return ((observed_guide_encoded*whitelist_guide_encoded[:, np.newaxis]).sum(axis=3)^1).sum(axis=2).flatten()
