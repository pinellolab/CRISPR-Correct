from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd

from typing import List, Tuple, Optional

def plot_mutation_count_histogram(mutation_counter, title="Bar Plot of Counts Sorted by Index", filename: Optional[str] = None):
    sorted_mutation_counter = sorted(mutation_counter.items())

    # Extract keys (indexes) and values from the sorted Counter
    indexes, counts = zip(*sorted_mutation_counter)

    # Create a bar plot
    fig, ax = plt.subplots(1, figsize=(5, 5))
    ax.bar(indexes, counts, color='blue')

    # Add labels and title
    ax.set_xlabel('Mutations Per Allele')
    ax.set_ylabel('Read Count')
    fig.suptitle(title)

    if filename is None:
        plt.show()
    else:
        fig.savefig(filename)



def plot_trinucleotide_mutational_signature(unlinked_mutations_df: pd.DataFrame, label: str, possible_refs: List[str] = ["A","C","G","T"], possible_alts: List[str] = ["A","C","G","T","X","N","-","X"], empty_context_character: str = "_", filename: Optional[str] = None):
    fig, ax = plt.subplots(1, figsize=(15, 5*1))

    if unlinked_mutations_df.shape[0] > 0: # If there are mutations
        
        '''
            Pre-processing mutation data
        '''
        # Prepare unique trinucleotide contexts based on provided REF and ALT sequences
        unique_trinucleotide_contexts_series_list: List[pd.Series] = []
        for ref_nt in possible_refs:
            for alt_nt in possible_alts:
                if ref_nt != alt_nt:
                    for preceding_ref_nt in possible_refs + [empty_context_character]:
                        for succeeding_ref_nt in possible_refs + [empty_context_character]:
                            unique_trinucleotide_contexts_series_list.append(pd.Series([preceding_ref_nt, ref_nt, alt_nt, succeeding_ref_nt], index=["preceding_nt_context", "ref_nt", "alt_nt", "succeeding_nt_context"]))
        unique_trinucleotide_contexts_df = pd.DataFrame(unique_trinucleotide_contexts_series_list)

        # Calculate the counts for each trinucleotide context
        unique_trinucleotide_contexts_counts_series_list: List[pd.Series] = []
        for mononucleotide_context_indices, mononucleotide_context_group in unique_trinucleotide_contexts_df.groupby(["ref_nt", "alt_nt"]):
            for row_index, row in mononucleotide_context_group.iterrows():
                row["count"] = unlinked_mutations_df.loc[(unlinked_mutations_df["preceding_nt_context"] == row["preceding_nt_context"]) & (unlinked_mutations_df["ref_nt"] == row["ref_nt"]) & (unlinked_mutations_df["alt_nt"] == row["alt_nt"]) & (unlinked_mutations_df["succeeding_nt_context"] == row["succeeding_nt_context"]), "count"].sum()
                unique_trinucleotide_contexts_counts_series_list.append(row)
        unique_trinucleotide_contexts_counts_df = pd.DataFrame(unique_trinucleotide_contexts_counts_series_list)

        '''
            Prepare constants for plotting
        '''
        # Define constants normalization and axis setting
        unlinked_mutations_df_groupby_trinucleotide_context_summed = unlinked_mutations_df.groupby(["preceding_nt_context", "ref_nt", "alt_nt", "succeeding_nt_context"]).sum("count")
        total_mutations = unlinked_mutations_df_groupby_trinucleotide_context_summed["count"].sum()
        max_context_frequency = (np.round(max(unlinked_mutations_df_groupby_trinucleotide_context_summed["count"]/total_mutations) * 100)) + 1 # Get the max context frequency to set plot y-axis limit
        
        '''
            Generation annotations of each group and their positions.
        '''
        current_position = 0
        variant_class_positions: List[Tuple[str, int]] = []
        
        # TODO: Do I need to sort by preceeding and succeeding context to ensure order is correct? Check later by printing
        
        # Collect the variant position labels
        unique_trinucleotide_contexts_counts_df_groupby_mononucleotide_context = unique_trinucleotide_contexts_counts_df.groupby(["ref_nt", "alt_nt"])
        for mononucleotide_context_indices, mononucleotide_context_group in unique_trinucleotide_contexts_counts_df_groupby_mononucleotide_context:
            position = current_position + (mononucleotide_context_group.shape[0]/2)
            current_position = current_position + mononucleotide_context_group.shape[0]
            variant_class_positions.append((f"{mononucleotide_context_indices[0]}>{mononucleotide_context_indices[1]}", position)) # Add the label 
        
        '''
            Perform plotting
        '''
        color_palette = sns.color_palette("Set2",unique_trinucleotide_contexts_counts_df_groupby_mononucleotide_context.ngroups)
        bar_colors = np.asarray([color_palette[group_index] for group_index, (_, mononucleotide_context_group) in enumerate(unique_trinucleotide_contexts_counts_df_groupby_mononucleotide_context) for _ in range(mononucleotide_context_group.shape[0])])
        variant_frequencies = pd.concat([(mononucleotide_context_group["count"]/total_mutations)*100 for mononucleotide_context_indices, mononucleotide_context_group in unique_trinucleotide_contexts_counts_df_groupby_mononucleotide_context])
        variant_frequencies.plot.bar(x='lab', y='val', rot=0, color = bar_colors, ax=ax)
        ax.set_ylim(0, max_context_frequency)
        ax.set_xticks([position for variant_class, position in variant_class_positions], [variant_class for variant_class, _ in variant_class_positions])
        ax.set_ylabel("Proportion of mutations (%)")
        ax.set_xlabel("Variant Type")
        ax.set_title(label)
    
    if filename is None:
        plt.show()    
    else:
        fig.savefig(filename)
    
def plot_positional_mutational_signature(unlinked_mutations_df: pd.DataFrame, label: str, possible_refs: List[str] = ["A","C","G","T"], possible_alts: List[str] = ["A","C","G","T","X","N","-","X"], empty_context_character: str = "_", min_position: int = 0, max_position: Optional[int] = None, filename: Optional[str] = None):
    fig, ax = plt.subplots(1, figsize=(15, 5*1))

    if unlinked_mutations_df.shape[0] > 0: # If there are mutations
        '''
            Pre-processing mutation data
        '''
        max_position = int(unlinked_mutations_df.position.quantile(q=0.95)) if max_position is None else max_position
        unique_positional_contexts_series_list: List[pd.Series] = []
        for ref_nt in possible_refs:
            for alt_nt in possible_alts:
                if ref_nt != alt_nt:
                    for position in range(0, max_position):
                        unique_positional_contexts_series_list.append(pd.Series([ref_nt, alt_nt, position], index=["ref_nt", "alt_nt", "position"]))
        unique_positional_contexts_df = pd.DataFrame(unique_positional_contexts_series_list)

        unique_positional_contexts_counts_series_list: List[pd.Series] = []
        for mononucleotide_context_indices, mononucleotide_context_group in unique_positional_contexts_df.groupby(["ref_nt", "alt_nt"]):
            for row_index, row in mononucleotide_context_group.iterrows():
                row["count"] = unlinked_mutations_df.loc[(unlinked_mutations_df["ref_nt"] == row["ref_nt"]) & (unlinked_mutations_df["alt_nt"] == row["alt_nt"]) & (unlinked_mutations_df["position"] == row["position"]), "count"].sum()
                unique_positional_contexts_counts_series_list.append(row)
        unique_positional_contexts_counts_df = pd.DataFrame(unique_positional_contexts_counts_series_list)

        '''
            Prepare constants for plotting
        '''
        # Define constants normalization and axis setting
        unlinked_mutations_df_groupby_positional_context_summed = unlinked_mutations_df.groupby(["ref_nt", "alt_nt", "position"]).sum("count")
        total_mutations = unlinked_mutations_df_groupby_positional_context_summed["count"].sum()
        max_context_frequency = (np.round(max(unlinked_mutations_df_groupby_positional_context_summed["count"]/total_mutations) * 100)) + 1 # Get the max context frequency to set plot y-axis limit
        
        '''
            Generation annotations of each group and their positions.
        '''
        current_position = 0
        variant_class_positions: List[Tuple[str, int]] = []
        
        # TODO: Do I need to sort by preceeding and succeeding context to ensure order is correct? Check later by printing
        
        # Collect the variant position labels
        unique_positional_contexts_counts_df_groupby_mononucleotide_context = unique_positional_contexts_counts_df.groupby(["ref_nt", "alt_nt"])
        for mononucleotide_context_indices, mononucleotide_context_group in unique_positional_contexts_counts_df_groupby_mononucleotide_context:
            position = current_position + (mononucleotide_context_group.shape[0]/2)
            current_position = current_position + mononucleotide_context_group.shape[0]
            variant_class_positions.append((f"{mononucleotide_context_indices[0]}>{mononucleotide_context_indices[1]}", position)) # Add the label 
        
        '''
            Perform plotting
        '''
        color_palette = sns.color_palette("Set2",unique_positional_contexts_counts_df_groupby_mononucleotide_context.ngroups)
        bar_colors = np.asarray([color_palette[group_index] for group_index, (_, mononucleotide_context_group) in enumerate(unique_positional_contexts_counts_df_groupby_mononucleotide_context) for _ in range(mononucleotide_context_group.shape[0])])
        variant_frequencies = pd.concat([(mononucleotide_context_group["count"]/total_mutations)*100 for mononucleotide_context_indices, mononucleotide_context_group in unique_positional_contexts_counts_df_groupby_mononucleotide_context])
        variant_frequencies.plot.bar(x='lab', y='val', rot=0, color = bar_colors, ax=ax)
        ax.set_ylim(0, max_context_frequency)
        ax.set_xticks([position for variant_class, position in variant_class_positions], [variant_class for variant_class, _ in variant_class_positions])
        ax.set_ylabel("Proportion of mutations (%)")
        ax.set_xlabel("Variant Type")
        ax.set_title(label)
    
    if filename is None:
        plt.show()    
    else:
        fig.savefig(filename)