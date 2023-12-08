###
### GUIDE MAPPING HELPER FUNCTIONS
###
import numpy as np

# Vector encoding of each base, since we want to vectorized
non_ambigious_encoding_dict = dict({
    "A": np.asarray([1,0,0,0]),
    "C": np.asarray([0,1,0,0]),
    "G": np.asarray([0,0,1,0]),
    "T": np.asarray([0,0,0,1]),
    "U": np.asarray([0,0,0,1]),
    "X": np.asarray([9999,9999,9999,9999]) # NOTE: This represents a padding token for sequences. Set to some high number, so that the calculated hamming distance is very high and won't be selected (which will happen if the observed sequence is greater than the whitelist sequence)
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
    "X": non_ambigious_encoding_dict["X"], # NOTE: This represents a padding token for sequences 
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

# Deprecated
def determine_hamming_distance_classic(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError(f"Sequences must have equal length. {seq1} and {seq2}")
    
    distance = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            distance += 1
    
    return distance
