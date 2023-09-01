'''
    NOTE 8/28/223: Most of the code is located in Juypter notebook http://localhost:30111/notebooks/PROJECTS/2023_03_BB_SensorDemultiplexAndAnalysis/20230324_NovaSeq_Demultiplex.ipynb in anchor 20230828-millipede-encoding-anchor-new
'''
from typing import Tuple, Union, Mapping
from multiprocessing import Pool
from functools import partial
import pandas as pd
import numpy as np
import os
from functools import reduce
import pickle
from collections import defaultdict
from dataclasses import dataclass
from .CrisprAmbiguousMapping import save_or_load_pickle

def determine_mutations_in_sequence(true_sequence, observed_sequence):
    mutation_array = []
    lowest_len = min(len(true_sequence), len(observed_sequence))
    for i in range(lowest_len):
        if true_sequence[i] != observed_sequence[i]:
            preceding_nt_context = None
            succeeding_nt_context = None
            ref_nt = true_sequence[i]
            alt_nt = observed_sequence[i]
            if i == 0:
                preceding_nt_context = "_"
                succeeding_nt_context = true_sequence[i+1]
            elif i == len(true_sequence)-1:
                preceding_nt_context = true_sequence[i-1]
                succeeding_nt_context = "_"
            else:
                preceding_nt_context = true_sequence[i-1]
                succeeding_nt_context = true_sequence[i+1]
                
            mutation_series = pd.Series((preceding_nt_context, ref_nt, alt_nt,succeeding_nt_context, i), index=["preceding_nt_context", "ref_nt", "alt_nt", "succeeding_nt_context", "position"])
            mutation_array.append(mutation_series)
    observed_sequence_mutation_df = pd.DataFrame(mutation_array)
    return observed_sequence_mutation_df


# TODO 4/3/23: All this could probably be sped up through vectorization
def get_mutational_profile(observed_reporter_series, whitelist_reporter_tuple: Tuple[str, str, str]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    observed_reporter_seq = observed_reporter_series["observed_sequence"]
    
    true_guide_seq = whitelist_reporter_tuple[0]
    true_surrogate_seq = whitelist_reporter_tuple[1]
    
    observed_guide_seq = observed_reporter_seq[0]
    observed_surrogate_seq = observed_reporter_seq[1]
    observed_surrogate_seq = observed_surrogate_seq[-len(true_surrogate_seq)-1:]
    
    observed_guide_mutation_df = determine_mutations_in_sequence(true_sequence=true_guide_seq, observed_sequence=observed_guide_seq)
    observed_surrogate_mutation_df = determine_mutations_in_sequence(true_sequence=true_surrogate_seq, observed_sequence=observed_surrogate_seq)

    return (observed_guide_mutation_df, observed_surrogate_mutation_df)


# Iterate through each whitelist reporter

'''
    Get the protospacer and surrogate mutation list for each whitelist entry. 
'''
def get_observed_mutation_per_whitelist_reporter(row, observed_guides_df):
    whitelist_reporter_tuple = tuple(row.values())
    
    observed_reporter_df = observed_guides_df[observed_guides_df["inferred_guides"]==whitelist_reporter_tuple]

    all_observed_guide_mutation_df_list = []
    all_observed_surrogate_mutation_df_list = []
    all_observed_reporter_series_mutation_count_series_list = []
    
    for _, observed_reporter_series in observed_reporter_df.iterrows():
        observed_count = observed_reporter_series["observed_counts"]
        observed_guide_mutation_df, observed_surrogate_mutation_df = get_mutational_profile(observed_reporter_series, whitelist_reporter_tuple)

        observed_guide_mutation_df["observed_count"] = observed_count
        observed_surrogate_mutation_df["observed_count"] = observed_count

        observed_guide_mutation_df["true_protospacer"] = row["protospacer"]
        observed_surrogate_mutation_df["true_protospacer"] = row["protospacer"]
        observed_surrogate_mutation_df["true_surrogate"] = row["surrogate"]

        # Get count of mutations per guide/surrogate
        
        observed_reporter_series_mutation_count_series = observed_reporter_series.copy()
        observed_reporter_series_mutation_count_series["guide_mutation_count"] = observed_guide_mutation_df.shape[0]
        observed_reporter_series_mutation_count_series["surrogate_mutation_count"] = observed_surrogate_mutation_df.shape[0]
        
        all_observed_guide_mutation_df_list.append(observed_guide_mutation_df)
        all_observed_surrogate_mutation_df_list.append(observed_surrogate_mutation_df)
        all_observed_reporter_series_mutation_count_series_list.append(observed_reporter_series_mutation_count_series)

    
    all_observed_guide_mutation_df = pd.concat(all_observed_guide_mutation_df_list) if len(all_observed_guide_mutation_df_list) != 0 else pd.DataFrame()
    all_observed_surrogate_mutation_df = pd.concat(all_observed_surrogate_mutation_df_list) if len(all_observed_surrogate_mutation_df_list) != 0 else pd.DataFrame()
    all_observed_reporter_series_mutation_count_df = pd.DataFrame(all_observed_reporter_series_mutation_count_series_list) if len(all_observed_reporter_series_mutation_count_series_list) != 0 else pd.DataFrame()

    return {
        "all_observed_guide_mutation_df": all_observed_guide_mutation_df,
        "all_observed_surrogate_mutation_df": all_observed_surrogate_mutation_df,
        "all_observed_reporter_series_mutation_count_df": all_observed_reporter_series_mutation_count_df
    }


# define a function to apply multiprocessing
'''
    Given the count result, generate the dataframes summarizing the mutations in the observed sequences
'''
def get_observed_mutation_dfs(count_result: Tuple[pd.Series, pd.DataFrame, Mapping], cores=1):
    whitelist_guide_reporter_counts, observed_guides_df, qc_dict = count_result
    whitelist_guide_reporter_df = whitelist_guide_reporter_counts.index.to_frame().reset_index(drop=True)

    # define the number of processes to use

    # create a Pool object with the number of processes defined above
    
    results = None
    
    get_observed_mutation_per_whitelist_reporter_p = partial(get_observed_mutation_per_whitelist_reporter, observed_guides_df=observed_guides_df)
    if cores == 1:
        results = [get_observed_mutation_per_whitelist_reporter_p(record) for record in whitelist_guide_reporter_df.to_dict("records")]
    else:
        print(f"Running parallelized with {cores} cores")
        with Pool(processes=cores) as pool:
            # apply the process_row function to each row of the whitelist_guide_reporter_df DataFrame using map()
            results = pool.map(get_observed_mutation_per_whitelist_reporter_p, whitelist_guide_reporter_df.to_dict("records"))
    
    # concatenate all the resulting DataFrames
    all_observed_guide_mutation_df = pd.concat([r["all_observed_guide_mutation_df"] for r in results])
    all_observed_surrogate_mutation_df = pd.concat([r["all_observed_surrogate_mutation_df"] for r in results])
    all_observed_reporter_series_mutation_count_df = pd.concat([r["all_observed_reporter_series_mutation_count_df"] for r in results], axis=0)

    return {
        "all_observed_guide_mutation_df": all_observed_guide_mutation_df,
        "all_observed_surrogate_mutation_df": all_observed_surrogate_mutation_df,
        "all_observed_reporter_series_mutation_count_df": all_observed_reporter_series_mutation_count_df
    }


'''
function copied from "anchor 20230703-millipede-encoding-anchor" in notebook "2021_11_BB_Shared_Tiling_Screen_Analysis/20220504_Davide_pilot_BE_analysis/STEP%202-%20Prepare%20base-editing%20allele%20encoding%20.ipynb"

TODO: Once we build the package for crispr-millipede, replace with using the function from crispr-millipede. There are a few differences, specifically had to change np.object and np.str to np.object_ and np.str_ probably due to numpy version differences. Also added N and changed assertion statement

'''
def get_substitution_encoding(aligned_sequence, original_seq, skip_index=0):
    assert len(aligned_sequence) == len(original_seq) # Ensure the aligned sequence (from allele table) is equal size to the reference sequence
    
    
    nucleotides = ["A","C","T","G","N","-"] # List of possible nucleotides
    encodings_per_position = []
    mismatch_mappings_per_position = []
    for index in range(0, len(original_seq)): # Iterate through each base and check for substitution
        # TODO Ensure sequences are uppercase
        nucleotides_mm = nucleotides[:]
        nucleotides_mm.remove(original_seq[index])
        mm_encoding = pd.Series(np.repeat(0, len(nucleotides)-1))
        if aligned_sequence[index] == original_seq[index]: # If the aligned sequence is same as reference
            pass
        else:
            mm_index = nucleotides_mm.index(aligned_sequence[index])
            mm_encoding[mm_index] = 1
        mismatch_mappings_per_position.append(nucleotides_mm)
        encodings_per_position.append(mm_encoding)

    encodings_per_position_df = pd.DataFrame(encodings_per_position).T
    mismatch_mappings_per_position_df = pd.DataFrame(mismatch_mappings_per_position).T

    encodings_per_position_df.columns = list(original_seq)
    mismatch_mappings_per_position_df.columns = list(original_seq)

    mismatch_mappings_per_position_POS_list = np.arange(mismatch_mappings_per_position_df.shape[1]).repeat(mismatch_mappings_per_position_df.shape[0])
    mismatch_mappings_per_position_REF_list = np.asarray(list(original_seq)).repeat(mismatch_mappings_per_position_df.shape[0]).astype(np.object_)
    mismatch_mappings_per_position_ALT_list = mismatch_mappings_per_position_df.T.values.flatten()
    mismatch_mappings_per_position_full_list = mismatch_mappings_per_position_POS_list.astype(np.str_).astype(object)+mismatch_mappings_per_position_REF_list + np.repeat(">", len(mismatch_mappings_per_position_REF_list)) + mismatch_mappings_per_position_ALT_list
    encodings_per_position_list = encodings_per_position_df.T.values.flatten()

    
    # Encodings per position DF, mismatch mappings per position DF, encodings per position flattened, mismatch mappings per position flattened, mismatch mapping position in flattened list, mismatch mapping ref in flattened list, mismatch mapping alt in flattened list, all substitutions made
    
    index = pd.MultiIndex.from_tuples(zip(mismatch_mappings_per_position_full_list, mismatch_mappings_per_position_POS_list, mismatch_mappings_per_position_REF_list, mismatch_mappings_per_position_ALT_list), names=["FullChange", "Position","Ref", "Alt"])
    
    assert len(encodings_per_position_list) == len(index), f"{len(encodings_per_position_list)} == {len(index)}"
    encodings_per_position_series = pd.Series(encodings_per_position_list, index = index, name="encoding")
    return encodings_per_position_series


def process_encoding(encoding_set):
    for encoding_df in encoding_set:
        encoding_df.columns = encoding_df.columns.get_level_values("FullChange")
        
        
def add_read_column(original_dfs, encoded_dfs, suffix):
    for i, original_dfs_rep in enumerate(original_dfs):
        encoded_dfs[i]["#Reads{}".format(suffix)] = original_dfs_rep["#Reads"]
        
        
def collapse_encodings(encoded_dfs):
    encoded_dfs_collapsed = []
    for encoded_df_rep in encoded_dfs:
        feature_colnames = [name for name in list(encoded_df_rep.columns) if "#Reads" not in name]
        encoded_dfs_collapsed.append(encoded_df_rep.groupby(feature_colnames, as_index=True).sum().reset_index())
    return encoded_dfs_collapsed


def merge_conditions_by_rep(first_encodings_collapsed, second_encodings_collapsed, third_encodings_collapsed):
    assert len(first_encodings_collapsed) == len(second_encodings_collapsed) == len(third_encodings_collapsed)
    encoded_dfs_merged = []
    for rep_i in range(len(first_encodings_collapsed)):
        feature_colnames |= [name for name in list(first_encodings_collapsed[rep_i].columns) if "#Reads" not in name]
        samples = [first_encodings_collapsed[rep_i], second_encodings_collapsed[rep_i], third_encodings_collapsed[rep_i]]
        df_encoding_rep1 = reduce(lambda  left,right: pd.merge(left,right,on=feature_colnames,
                                                how='outer'), samples).fillna(0)
        
        encoded_dfs_merged.append(df_encoding_rep1)
    return encoded_dfs_merged


@dataclass
class GuideReporterEncodingStrategy:
    output_dir: str
    date_string: str


def merge_and_encode_observed_surrogate(args):
    true_surrogate_sequence, condition_dataframes = args
    dataframes = [dataframe for condition, dataframe in condition_dataframes]
    
    if len(dataframes) > 1:
        merged_observed_surrogate_sequence_df = reduce(lambda left, right: pd.merge(left, right, how="outer", on=["observed_protospacer", "observed_surrogate"]), dataframes)
    else:
        merged_observed_surrogate_sequence_df = pd.DataFrame(dataframes[0])
    
    merged_observed_surrogate_sequence_df = merged_observed_surrogate_sequence_df.fillna(0.0)

    merged_observed_surrogate_sequence_df["observed_surrogate"] = merged_observed_surrogate_sequence_df["observed_surrogate"].apply(lambda sequence: sequence[-len(true_surrogate_sequence[1]):])
    merged_observed_surrogate_sequence_df = merged_observed_surrogate_sequence_df.groupby(["observed_protospacer", "observed_surrogate"]).sum()
    merged_observed_surrogate_sequence_df["edited_protospacer"] = merged_observed_surrogate_sequence_df.index.get_level_values('observed_protospacer') != true_surrogate_sequence[0]
    merged_observed_surrogate_sequence_df["edited_surrogate"] = merged_observed_surrogate_sequence_df.index.get_level_values('observed_surrogate') != true_surrogate_sequence[1]
    merged_observed_surrogate_sequence_df = merged_observed_surrogate_sequence_df.sort_values(f"observed_counts_{condition_dataframes[0][0]}", ascending=False)
    
    merged_observed_surrogate_sequence_surrogate_encoding_df = None
    if len(merged_observed_surrogate_sequence_df.index) > 0:
        merged_observed_surrogate_sequence_surrogate_encoding_df = pd.concat([get_substitution_encoding(observed_surrogate, true_surrogate_sequence[1]) for observed_surrogate in merged_observed_surrogate_sequence_df.index.get_level_values('observed_surrogate')], axis=1).transpose()
        nt_columns = merged_observed_surrogate_sequence_surrogate_encoding_df.columns
        merged_observed_surrogate_sequence_surrogate_encoding_df.index = merged_observed_surrogate_sequence_df.index.get_level_values('observed_surrogate')

        for colname in merged_observed_surrogate_sequence_df.columns:
            merged_observed_surrogate_sequence_surrogate_encoding_df[colname] = merged_observed_surrogate_sequence_df[colname].values
        
        all_columns = merged_observed_surrogate_sequence_surrogate_encoding_df.columns
        # Can only groupby on single index, so temporarily setting column as FullChange index
        merged_observed_surrogate_sequence_surrogate_encoding_df.columns = list(merged_observed_surrogate_sequence_surrogate_encoding_df.columns.get_level_values("FullChange"))
        merged_observed_surrogate_sequence_surrogate_encoding_df_groups = merged_observed_surrogate_sequence_surrogate_encoding_df.groupby(list(nt_columns.get_level_values("FullChange")), as_index=False)
        group_index = [group.index[0] for _, group in merged_observed_surrogate_sequence_surrogate_encoding_df_groups]
        # TODO: The edited_surrogate and edited_protospacer are edited as well
        merged_observed_surrogate_sequence_surrogate_encoding_df = pd.concat([group.sum() for _, group in merged_observed_surrogate_sequence_surrogate_encoding_df_groups], axis=1).transpose()
        merged_observed_surrogate_sequence_surrogate_encoding_df.index = group_index
        merged_observed_surrogate_sequence_surrogate_encoding_df.columns = all_columns

    merged_observed_surrogate_sequence_protospacer_encoding_df = None
    if len(merged_observed_surrogate_sequence_df.index) > 0:
        merged_observed_surrogate_sequence_protospacer_encoding_df = pd.concat([get_substitution_encoding(observed_protospacer, true_surrogate_sequence[0]) for observed_protospacer in merged_observed_surrogate_sequence_df.index.get_level_values('observed_protospacer')], axis=1).transpose()
        nt_columns = merged_observed_surrogate_sequence_protospacer_encoding_df.columns
        merged_observed_surrogate_sequence_protospacer_encoding_df.index = merged_observed_surrogate_sequence_df.index.get_level_values('observed_protospacer')

        for colname in merged_observed_surrogate_sequence_df.columns:
            merged_observed_surrogate_sequence_protospacer_encoding_df[colname] = merged_observed_surrogate_sequence_df[colname].values
        
        all_columns = merged_observed_surrogate_sequence_protospacer_encoding_df.columns
        # Can only groupby on single index, so temporarily setting column as FullChange index
        merged_observed_surrogate_sequence_protospacer_encoding_df.columns = list(merged_observed_surrogate_sequence_protospacer_encoding_df.columns.get_level_values("FullChange"))
        merged_observed_surrogate_sequence_protospacer_encoding_df_groups = merged_observed_surrogate_sequence_protospacer_encoding_df.groupby(list(nt_columns.get_level_values("FullChange")), as_index=False)
        group_index = [group.index[0] for _, group in merged_observed_surrogate_sequence_protospacer_encoding_df_groups]
        merged_observed_surrogate_sequence_protospacer_encoding_df = pd.concat([group.sum() for _, group in merged_observed_surrogate_sequence_protospacer_encoding_df_groups], axis=1).transpose()
        merged_observed_surrogate_sequence_protospacer_encoding_df.index = group_index
        merged_observed_surrogate_sequence_protospacer_encoding_df.columns = all_columns

        
        
    return true_surrogate_sequence, merged_observed_surrogate_sequence_df, merged_observed_surrogate_sequence_surrogate_encoding_df, merged_observed_surrogate_sequence_protospacer_encoding_df

    
def process_tiling_reporter_counts(sample_pooling_guidelibrary_df: pd.DataFrame, date_string: str, screen_column_group: str, replicate_column: Union[int, str], condition_column, guide_reporter_encoding_strategy, cores=1):
    # Modify below:
    surrogate_output_sub_dir = f"{guide_reporter_encoding_strategy.output_dir}"

    if not os.path.exists(surrogate_output_sub_dir):
        os.makedirs(surrogate_output_sub_dir)
        
    merged_observed_surrogate_sequence_df_dict_screen_reps_library = {}
    merged_observed_surrogate_sequence_surrogate_encoding_df_dict_screen_reps_library = {}
    merged_observed_surrogate_sequence_protospacer_encoding_df_dict_screen_reps_library = {}
    for guide_library_id, library_group in sample_pooling_guidelibrary_df.groupby("guide_library_id"):
        print(f"Processing guide library {guide_library_id}")
        
        merged_observed_surrogate_sequence_df_dict_screen_reps = {}
        merged_observed_surrogate_sequence_surrogate_encoding_df_dict_screen_reps = {}
        merged_observed_surrogate_sequence_protospacer_encoding_df_dict_screen_reps = {}
            
        # Iterate through each screen group
        for screen_name, screen_group in library_group.groupby(screen_column_group):
            print(f"Processing screen {screen_name}")

            merged_observed_surrogate_sequence_df_dict_screen_rep = {}
            merged_observed_surrogate_sequence_surrogate_encoding_df_dict_screen_rep = {}
            merged_observed_surrogate_sequence_protospacer_encoding_df_dict_screen_rep = {}
            for replicate, replicate_group in screen_group.groupby(replicate_column):
                print(f"Processing replicate {replicate}")

                condition_observed_surrogate_sequence_df_dict = defaultdict(list)

                # Iterate through each screen condition
                conditions = list(set(replicate_group[condition_column].astype(str)))
                for condition in conditions:
                    print(f"Processing condition {condition}")
                    # Get the barcode for the screen condition sample
                    
                    # Get the metadata needed to construct the filename of the count matrix file
                    i5_index = replicate_group.loc[replicate_group[condition_column].astype(str) == condition, "i5_index"].iloc[0]
                    sample_sub_dir = sample_indices_df_subset.loc[sample_indices_df_subset["i5_index"] == i5_index,"sub_library_name"].iloc[0]
                    U6PE1_barcode = replicate_group.loc[replicate_group[condition_column].astype(str) == condition, "U6PE1_Barcode"].iloc[0]
                    with open(f"{guide_reporter_counting_strategy.output_dir}/{sample_sub_dir}/{sample_sub_dir}_#{U6PE1_barcode}_count_{date_string}.pickle", 'rb') as handle:

                        count_result = pickle.load(handle)

                        whitelist_guide_reporter_counts, observed_guides_df, qc_dict = count_result
                        whitelist_guide_reporter_list = list(whitelist_guide_reporter_counts.index)

                        # Get the surrogate count df for the screen condition sample
                        for true_surrogate_sequence in whitelist_guide_reporter_list:
                            observed_surrogate_info_sequence_df = observed_guides_df[observed_guides_df["inferred_guides"] == true_surrogate_sequence].loc[:,["observed_sequence", "observed_counts"]]
                            observed_surrogate_sequence_df = pd.DataFrame({"observed_protospacer": observed_surrogate_info_sequence_df["observed_sequence"].apply(lambda info: info[0]), "observed_surrogate": observed_surrogate_info_sequence_df["observed_sequence"].apply(lambda info: info[1]), f"observed_counts_{condition}": observed_surrogate_info_sequence_df["observed_counts"]})
                            condition_observed_surrogate_sequence_df_dict[true_surrogate_sequence].append((condition, observed_surrogate_sequence_df))

                # TODO: For cores == 1, using regular list comprehension. Not sure how to deal with the two dicts that need updating, perhaps modify function to return if cores==1
                # Create the encoding for each whitelist surrogate (parallelized)
                print(f"Performing encoding of {replicate} with cores {cores}")
                if cores > 1:
                    with Pool(processes=cores) as pool:
                        # Prepare the arguments for the parallel processing
                        args_list = [(true_surrogate_sequence, condition_dataframes)
                                     for true_surrogate_sequence, condition_dataframes in condition_observed_surrogate_sequence_df_dict.items()]

                        # Use the Pool's map function to process the loop iterations in parallel
                        results = pool.map(merge_and_encode_observed_surrogate, args_list)
                else:
                    # Prepare the arguments for the parallel processing
                    args_list = [(true_surrogate_sequence, condition_dataframes)
                                 for true_surrogate_sequence, condition_dataframes in condition_observed_surrogate_sequence_df_dict.items()]

                    # Use the Pool's map function to process the loop iterations in parallel
                    results = [merge_and_encode_observed_surrogate(args) for args in args_list]
                print(f"Completed encoding of {replicate}. Processing results") 
                # Parse the parallelized result into the dictionaries
                merged_observed_surrogate_sequence_df_dict = {} # Contains counts per condition samples
                merged_observed_surrogate_sequence_surrogate_encoding_df_dict = {} # Contain encodings per condition samples
                merged_observed_surrogate_sequence_protospacer_encoding_df_dict = {} # Contain encodings per condition samples
                for true_surrogate_sequence, merged_observed_surrogate_sequence_df, merged_observed_surrogate_sequence_surrogate_encoding_df, merged_observed_surrogate_sequence_protospacer_encoding_df in results:
                    merged_observed_surrogate_sequence_df_dict[true_surrogate_sequence] = merged_observed_surrogate_sequence_df
                    merged_observed_surrogate_sequence_surrogate_encoding_df_dict[true_surrogate_sequence] = merged_observed_surrogate_sequence_surrogate_encoding_df
                    merged_observed_surrogate_sequence_protospacer_encoding_df_dict[true_surrogate_sequence] = merged_observed_surrogate_sequence_protospacer_encoding_df

                # Add dictionaries to main list
                merged_observed_surrogate_sequence_df_dict_screen_rep[replicate] = merged_observed_surrogate_sequence_df_dict
                merged_observed_surrogate_sequence_surrogate_encoding_df_dict_screen_rep[replicate] = merged_observed_surrogate_sequence_surrogate_encoding_df_dict
                merged_observed_surrogate_sequence_protospacer_encoding_df_dict_screen_rep[replicate] = merged_observed_surrogate_sequence_protospacer_encoding_df_dict

            merged_observed_surrogate_sequence_df_dict_screen_reps[screen_name] = merged_observed_surrogate_sequence_df_dict_screen_rep
            merged_observed_surrogate_sequence_surrogate_encoding_df_dict_screen_reps[screen_name] = merged_observed_surrogate_sequence_surrogate_encoding_df_dict_screen_rep
            merged_observed_surrogate_sequence_protospacer_encoding_df_dict_screen_reps[screen_name] = merged_observed_surrogate_sequence_protospacer_encoding_df_dict_screen_rep
        
        print(f"Completed guide library encoding: {guide_library_id}. Now saving.")
        save_or_load_pickle(f"{surrogate_output_sub_dir}/", f"merged_observed_surrogate_sequence_df_dict_screen_reps_{guide_library_id}", py_object=merged_observed_surrogate_sequence_df_dict_screen_reps, date_string=guide_reporter_encoding_strategy.date_string)
        save_or_load_pickle(f"{surrogate_output_sub_dir}/", f"merged_observed_surrogate_sequence_surrogate_encoding_df_dict_screen_reps_{guide_library_id}", py_object=merged_observed_surrogate_sequence_surrogate_encoding_df_dict_screen_reps, date_string=guide_reporter_encoding_strategy.date_string)
        save_or_load_pickle(f"{surrogate_output_sub_dir}/", f"merged_observed_surrogate_sequence_protospacer_encoding_df_dict_screen_reps_{guide_library_id}", py_object=merged_observed_surrogate_sequence_protospacer_encoding_df_dict_screen_reps, date_string=guide_reporter_encoding_strategy.date_string)
    
        merged_observed_surrogate_sequence_df_dict_screen_reps_library[guide_library_id] = merged_observed_surrogate_sequence_df_dict_screen_reps
        merged_observed_surrogate_sequence_surrogate_encoding_df_dict_screen_reps_library[guide_library_id] = merged_observed_surrogate_sequence_surrogate_encoding_df_dict_screen_reps
        merged_observed_surrogate_sequence_protospacer_encoding_df_dict_screen_reps_library[guide_library_id] = merged_observed_surrogate_sequence_protospacer_encoding_df_dict_screen_reps

    print("Completed all encoding, saving results")
    save_or_load_pickle(f"{surrogate_output_sub_dir}/", "merged_observed_surrogate_sequence_df_dict_screen_reps_library", py_object=merged_observed_surrogate_sequence_df_dict_screen_reps_library, date_string=guide_reporter_encoding_strategy.date_string)
    save_or_load_pickle(f"{surrogate_output_sub_dir}/", "merged_observed_surrogate_sequence_surrogate_encoding_df_dict_screen_reps_library", py_object=merged_observed_surrogate_sequence_surrogate_encoding_df_dict_screen_reps_library, date_string=guide_reporter_encoding_strategy.date_string)
    save_or_load_pickle(f"{surrogate_output_sub_dir}/", "merged_observed_surrogate_sequence_protospacer_encoding_df_dict_screen_reps_library", py_object=merged_observed_surrogate_sequence_protospacer_encoding_df_dict_screen_reps_library, date_string=guide_reporter_encoding_strategy.date_string)

    return merged_observed_surrogate_sequence_df_dict_screen_reps_library, merged_observed_surrogate_sequence_surrogate_encoding_df_dict_screen_reps_library, merged_observed_surrogate_sequence_protospacer_encoding_df_dict_screen_reps_library
