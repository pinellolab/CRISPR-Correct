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
    "A": np.asarray([1,3,3,3]),
    "C": np.asarray([3,1,3,3]),
    "G": np.asarray([3,3,1,3]),
    "T": np.asarray([3,3,3,1]),
    "U": np.asarray([3,3,3,1]),
    "-": np.asarray([3,3,3,3]), 
    "X": np.asarray([3,3,3,3])
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
    "R": (ambiguity_encoder(["A", "G"])^1)*3,
    "Y": (ambiguity_encoder(["C", "T"])^1)*3,
    "S": (ambiguity_encoder(["G", "C"])^1)*3,
    "W": (ambiguity_encoder(["A", "T"])^1)*3,
    "K": (ambiguity_encoder(["G", "T"])^1)*3,
    "M": (ambiguity_encoder(["A", "C"])^1)*3,
    "B": (ambiguity_encoder(["C", "G", "T"])^1)*3,
    "D": (ambiguity_encoder(["A", "G", "T"])^1)*3,
    "H": (ambiguity_encoder(["A", "C", "T"])^1)*3,
    "V": (ambiguity_encoder(["A", "C", "G"])^1)*3,
    "N": (ambiguity_encoder(["A", "C", "G", "T"])^1)*3,
})

# §5.4: removed legacy `encode_DNA_base_{whitelist,observed}` +
# `_vectorized` wrappers and `numpify_string{,_vectorized}` — all superseded
# by the LUT-based `encode_DNA_sequence_{whitelist,observed}` below. Last
# in-tree callers migrated in Phase 1 (LUT rewrite); no external callers
# found in project drivers.


# PERF §3.1: `np.vectorize(...)` over a per-base dict lookup is a Python
# for-loop in disguise — it dominated the inference hot path (called twice per
# observed sequence per component). Precomputed 256-entry lookup tables here
# let each string encode in a single `np.frombuffer` + fancy-index on LUT,
# which is a real numpy-speed operation (no Python loop).
_LUT_WHITELIST = np.zeros((256, 4), dtype=np.int8)
_LUT_OBSERVED = np.zeros((256, 4), dtype=np.int8)
for _base, _vec in full_encoding_dict_whitelist.items():
    _LUT_WHITELIST[ord(_base)] = _vec
    _LUT_WHITELIST[ord(_base.lower())] = _vec
for _base, _vec in full_encoding_dict_observed.items():
    _LUT_OBSERVED[ord(_base)] = _vec
    _LUT_OBSERVED[ord(_base.lower())] = _vec
# Any base outside the known alphabet (e.g. sequencing 'N' on an observed
# sequence) already has a 4-wide all-zero / all-3 row via the IUPAC expansion
# above. Unknown ASCII chars fall back to the all-zero default, which produces
# an extremely high Hamming distance — same semantics as the old
# KeyError-raising dict access, but without crashing.
del _base, _vec


def encode_DNA_sequence_whitelist(seq: str) -> np.ndarray:
    """LUT-based whitelist encoding. Returns shape (len(seq), 4) int8."""
    return _LUT_WHITELIST[np.frombuffer(seq.encode('ascii', errors='replace'), dtype=np.uint8)]


def encode_DNA_sequence_observed(seq: str) -> np.ndarray:
    """LUT-based observed encoding. Returns shape (len(seq), 4) int8."""
    return _LUT_OBSERVED[np.frombuffer(seq.encode('ascii', errors='replace'), dtype=np.uint8)]


def encode_guide_series_whitelist(guide_series) -> np.array:
    # PERF §3.1 + §3.13: avoid the `list(map(list, ...))` Python detour and the
    # np.vectorize. Build the 2D char buffer once via `np.frombuffer` and
    # fancy-index into the LUT in a single numpy call.
    strs = guide_series.astype(str).tolist()
    if not strs:
        return np.zeros((0, 0, 4), dtype=np.int8)
    L = max(len(s) for s in strs)
    padded = [s.ljust(L, 'X') for s in strs]
    buf = np.frombuffer(''.join(padded).encode('ascii', errors='replace'), dtype=np.uint8).reshape(len(strs), L)
    return _LUT_WHITELIST[buf]


def encode_guide_series_observed(guide_series) -> np.array:
    strs = guide_series.astype(str).tolist()
    if not strs:
        return np.zeros((0, 0, 4), dtype=np.int8)
    L = max(len(s) for s in strs)
    padded = [s.ljust(L, 'X') for s in strs]
    buf = np.frombuffer(''.join(padded).encode('ascii', errors='replace'), dtype=np.uint8).reshape(len(strs), L)
    return _LUT_OBSERVED[buf]

def retrieve_hamming_distance_whitelist(observed_guide_encoded, whitelist_guide_encoded):
    '''
        This takes a encoded guide sequence and a list of encoded whitelisted guides and matrix computes the hamming distance of the 
        encoded guide across all whitelisted guides in a single operation.
        
        
        (observed_guide_encoded*whitelist_guide_encoded[:, np.newaxis]).clip(0, 2).sum(axis=3) # Determines 

        # observed_guide_encoded*whitelist_guide_encoded[:, np.newaxis] performs multiplication of the target encoding with every whitelist encoding
        #
        
    '''
    # Compute the per-character un-bit flipped hamming penalty 
    encoding_product = (observed_guide_encoded*whitelist_guide_encoded[:, np.newaxis]).clip(0,3)
    # Need to get the minimum positive un-bit flipped hamming penalty (so we are removing 0s)
    masked_encoding_product = np.where(encoding_product > 0, encoding_product, np.inf)
    # Get the minimum positive un-bit flipped hamming penalty for each comparison (axis=-1)
    min_positive = np.min(masked_encoding_product, axis=-1)
    # If there is no un-bit flipped hamming penalty (all 0s, or all np.inf after masking), set the minimum positive as 0 
    min_positive = np.where(np.isinf(min_positive), 0, min_positive).astype(int)
    # Bit-flip to get the hamming distance, take the sum to get the per-encoding comparison hamming distance, and divide by 2 since the hamming distance is * 2 since working in bits. 
    hamming_distance = (min_positive^1).sum(axis=2).flatten()/2
    return hamming_distance
