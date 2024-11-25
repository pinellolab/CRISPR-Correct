from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import anndata as ad
from ..models.editing_models import MatchSetWhitelistReporterObservedSequenceMutationProfiles
from Bio.Seq import Seq
from typing import List, Tuple, Optional, Dict
import matplotlib.patches as mpatches
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Rectangle


def show_reporter_encoding_frequency(encoding_df, sequence_type, count_column_name="count"):
    encoding_df_count = encoding_df.loc[:, count_column_name]
    encoding_df_type = encoding_df.loc[:, encoding_df.columns.get_level_values("SequenceType") == sequence_type]
    encoding_frequency_series = encoding_df_type.mul(encoding_df_count, axis=0).sum(axis=0) / encoding_df_count.sum() * 100
    return encoding_frequency_series

def plot_editing_frequency_heatmap(replicate_scores, protospacer_start_coordinate:int, protospacer_size:int = 20, title: str=""):
    # Group by Position and sum the frequencies
    frequency_dfs = []
    for replicate_score in replicate_scores:
        frequency_df = pd.DataFrame(replicate_score, columns=["Frequency"]).groupby(['Position', 'Ref'])["Frequency"].sum().reset_index()
        frequency_dfs.append(frequency_df)

    # Create a single-row DataFrame for the heatmap
    heatmap_data = pd.DataFrame(columns=range(frequency_dfs[0].shape[0]))
    for rep_i, frequency_df in enumerate(frequency_dfs):
        heatmap_data.loc[rep_i] = [0.0] * frequency_df.shape[0]  # Initialize with zeros

    # Add the editing frequencies to the heatmap input dataframe
    for index, _ in frequency_dfs[0].iterrows():
        position = frequency_dfs[0].loc[index]['Position']
        for rep_i, frequency_df in enumerate(frequency_dfs):
            heatmap_data.at[rep_i, position] = frequency_df.loc[index]['Frequency']

    # Prepare annotations (only show greater than 5% editing)
    new_annot = heatmap_data.copy()
    new_annot[new_annot <= 5] = ""

    # Custom function to format numbers to ".2f" and leave strings unchanged
    def format_numbers(val):
        if isinstance(val, (int, float)):
            return f'{val:.1f}%'
        return val

    # Apply the custom function to each element in the DataFrame
    new_annot = new_annot.applymap(format_numbers)

    # Plot the heatmap using seaborn
    plt.figure(figsize=(35, 1.5))
    heatmap = sns.heatmap(heatmap_data, cmap='Reds', annot=new_annot, fmt='', cbar=True, vmin=0, vmax=100,
                          linecolor='grey', linewidths=1, cbar_kws={"shrink": 0.7, "pad": 0.01} )
    plt.xlabel('Position')
    plt.ylabel('Replicate')
    plt.title('Editing Frequencies Heatmap')

    # Set reference bases as labels
    bold_font = FontProperties(weight='bold')
    heatmap.set_xticklabels(frequency_df["Ref"], fontproperties=bold_font)
    heatmap.tick_params(axis='x', labelsize=16)

    heatmap.set_yticklabels(heatmap.get_yticklabels(), rotation=0)

    rect = Rectangle((protospacer_start_coordinate, 0.00), protospacer_size, 1, fill=False, edgecolor='black', linewidth=3)

    heatmap.add_patch(rect)
    sns.despine()
    plt.title(title)
    plt.show()

def plot_lfc_scatter_plot(anndata: ad.AnnData, 
                          reference_sequence_dict: Dict[str, Seq],
                          groupby_columns_list: List[str], 
                          enriched_population_label: str, 
                          baseline_population_label: str, 
                          population_column_name:str, 
                          count_size_population_label: Optional[str] = None,
                          guide_type_column_name:str = "type",
                          guide_targeting_type_label: str = "targeting",
                          is_igRNA_column_name: Optional[str] = "is_igRNA", 
                          strand_column_name: Optional[str] = "strand", 
                          coordinate_column_name: str = "coordinate", # Set to hg38_editsite_coordinate
                          chromosome_column_name: str = "chromosome",
                          score_layer_name: str ="lfc",
                          count_size_max_threshold: int = 300,
                          figure_width: int = 100,
                          figure_height: int = 6,
                          count_to_scale_factor: int = 0.5,
                          perform_assertion: bool = True
                          ):
    for group in anndata.obs.groupby(groupby_columns_list):
        enriched_pop_sample_df = group[1][group[1][population_column_name] == enriched_population_label]
        if perform_assertion is True:
            assert enriched_pop_sample_df.shape[0] == 1, f"Multiple samples matching enriched condition in group {group[1]}"
        enriched_pop_sample = enriched_pop_sample_df
        
        baseline_pop_sample_df = group[1][group[1][population_column_name] == baseline_population_label]
        if perform_assertion is True:
            assert baseline_pop_sample_df.shape[0] == 1, f"Multiple samples matching enriched condition in group {group[1]}"
        baseline_pop_sample = baseline_pop_sample_df

        if count_size_population_label is not None:
            count_size_pop_sample_df = group[1][group[1][population_column_name] == count_size_population_label]
            if perform_assertion is True:
                assert count_size_pop_sample_df.shape[0] == 1, f"Multiple samples matching condition to calculate count size in group {group[1]}"
            count_size_pop_sample = count_size_pop_sample_df
        
        enriched_anndata = anndata[enriched_pop_sample.index, :] # DEVELOPER NOTE: Will thow error of indices of samples are changed from ints
        baseline_anndata = anndata[baseline_pop_sample.index, :] # DEVELOPER NOTE: Will thow error of indices of samples are changed from ints
        if count_size_population_label is not None:
            count_size_anndata = anndata[count_size_pop_sample.index, :] # DEVELOPER NOTE: Will thow error of indices of samples are changed from ints
        
        
        enriched_targeting_anndata = enriched_anndata[:, enriched_anndata.var[guide_type_column_name] == guide_targeting_type_label]
        baseline_targeting_anndata = baseline_anndata[:, baseline_anndata.var[guide_type_column_name] == guide_targeting_type_label]
        if count_size_population_label is not None:
            count_size_targeting_anndata = count_size_anndata[:, count_size_anndata.var[guide_type_column_name] == guide_targeting_type_label]

        # Get igRNAs and pgRNtargeting_anndataAs separately
        is_igRNA = None
        is_pgRNA = None
        if is_igRNA_column_name is not None:
            is_igRNA = enriched_targeting_anndata.var[is_igRNA_column_name]
            is_pgRNA = ~is_igRNA
        else:
            # If igRNA not provided, set all guides as pgRNA by default
            is_igRNA = pd.Series([False] * len(enriched_targeting_anndata), index=enriched_targeting_anndata.index)
            is_pgRNA = pd.Series([True] * len(enriched_targeting_anndata), index=enriched_targeting_anndata.index)

        positive_strand = None
        negative_strand = None
        if strand_column_name is not None:
            positive_strand = enriched_targeting_anndata.var[strand_column_name] == "+"
            negative_strand = ~positive_strand
        else:
            # If strand not provided, set all guides as positive strand by default
            negative_strand = pd.Series([False] * len(enriched_targeting_anndata), index=enriched_targeting_anndata.index)
            positive_strand = pd.Series([True] * len(enriched_targeting_anndata), index=enriched_targeting_anndata.index)

        # Get the guide coordinates
        coordinates = enriched_targeting_anndata.var[coordinate_column_name]

        chromosome = enriched_targeting_anndata.var[chromosome_column_name][0]
        ids_pgRNA = enriched_targeting_anndata[:, is_pgRNA].var.index
        positive_strand_pgRNA = positive_strand[is_pgRNA] 

        coordinates_pgRNA = coordinates[is_pgRNA.values]
        scores_pgRNA = enriched_targeting_anndata[:, is_pgRNA].layers[score_layer_name].mean(axis=0).flatten()
        if count_size_population_label is not None:
            counts_pgRNA = count_size_targeting_anndata[:, is_pgRNA].X[0]
        else:
            counts_pgRNA = enriched_targeting_anndata[:, is_pgRNA].X[0] + baseline_targeting_anndata[:, is_pgRNA].X[0]
        counts_pgRNA[counts_pgRNA>count_size_max_threshold] = count_size_max_threshold

        coordinates_igRNA = coordinates[is_igRNA]
        ids_igRNA = enriched_targeting_anndata[:, is_igRNA].var.index
        scores_igRNA = enriched_targeting_anndata[:, is_igRNA].layers[score_layer_name].mean(axis=0).flatten()
        if count_size_population_label is not None:
            counts_pgRNA = count_size_targeting_anndata[:, is_igRNA].X[0]
        else:
            counts_igRNA = enriched_targeting_anndata[:, is_igRNA].X[0] + baseline_targeting_anndata[:, is_igRNA].X[0]
        counts_igRNA[counts_igRNA>count_size_max_threshold] = count_size_max_threshold
        
        positive_strand_igRNA = positive_strand[is_igRNA] 

        fig, axes = plt.subplots(1, 1, figsize=(figure_width, figure_height), sharex=True)
        edge_colors = ['black', 'orange']
        
        scatter = axes.scatter(coordinates_pgRNA, scores_pgRNA, alpha=0.7, c="blue", edgecolors=np.where(positive_strand_pgRNA, edge_colors[0], edge_colors[1]), label="pgRNA", s=counts_pgRNA*count_to_scale_factor)
        if is_igRNA_column_name is not None:
            axes.scatter(coordinates_igRNA, scores_igRNA, alpha=0.3, c="black", edgecolors=np.where(positive_strand_pgRNA, edge_colors[0], edge_colors[1]), label="igRNA", s=counts_igRNA*count_to_scale_factor)

        ylim = axes.get_ylim()
        axes.set_ylim(ylim[0], ylim[1]+1)
        xlim = axes.get_xlim()
        for coordinate in range(int(xlim[0])+2,int(xlim[1])-1):
            axes.text(coordinate-1, ylim[0]+0.1, reference_sequence_dict[chromosome][coordinate], fontsize=10, alpha=1)

        # Create legends for point sizes, edge colors, and dot colors
        if is_igRNA_column_name is not None:
            dotcolor_legend = axes.legend(title="Type", loc="upper right")
            plt.gca().add_artist(dotcolor_legend)
        if strand_column_name:
            edgecolor_patches = [mpatches.Patch(color=color, label=["+", "-"][i]) for i, color in enumerate(edge_colors)]
            edgecolor_legend = plt.legend(handles=edgecolor_patches, title="Strand (Edge Color)", loc="upper right", bbox_to_anchor=(1.0, 0.83))
            plt.gca().add_artist(edgecolor_legend)

        if count_size_population_label is not None:
            size_title = f"{count_size_population_label} Count"
        else:
            size_title = "Total Count"
        size_legend = plt.legend(*scatter.legend_elements(prop="sizes", num=5, func=lambda x: x * 1/count_to_scale_factor), title=size_title, loc="upper right", bbox_to_anchor=(1.0, 0.65))
        plt.gca().add_artist(size_legend)
        
        fig.suptitle("_".join([str(term) for term in group[0]]))
        plt.show() 


def plot_mutation_count_histogram(mutation_counter, title="Bar Plot of Counts Sorted by Index", filename: Optional[str] = None):
    fig, ax = plt.subplots(1, figsize=(5, 5))

    if sum(mutation_counter) > 0: # If there exists mutation counts
        sorted_mutation_counter = sorted(mutation_counter.items())

        # Extract keys (indexes) and values from the sorted Counter
        indexes, counts = zip(*sorted_mutation_counter)

        # Create a bar plot
        
        ax.bar(indexes, counts, color='blue')

        # Add labels and title
        ax.set_xlabel('Mutations Per Allele')
        ax.set_ylabel('Read Count')
    
    fig.suptitle(title)
    if filename is None:
        plt.show()
    else:
        fig.savefig(filename)



def plot_trinucleotide_mutational_signature(mutations_results: MatchSetWhitelistReporterObservedSequenceMutationProfiles, label: str, count_attribute_name="ambiguous_accepted_umi_noncollapsed_mutations", unlinked_mutation_attribute_name = "all_observed_surrogate_unlinked_mutations_df", possible_refs: List[str] = ["A","C","G","T"], possible_alts: List[str] = ["A","C","G","T","X","N","-","X"], empty_context_character: str = "_", filename: Optional[str] = None):
    fig, ax = plt.subplots(1, figsize=(15, 5*1))

    mutation_results_typed = getattr(mutations_results, count_attribute_name)
    if mutation_results_typed is not None:
        unlinked_mutations_df: pd.DataFrame = getattr(mutation_results_typed, unlinked_mutation_attribute_name)

        if (unlinked_mutations_df is not None) and (unlinked_mutations_df.shape[0] > 0): # If there are mutations
            
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
    
def plot_positional_mutational_signature(mutations_results: MatchSetWhitelistReporterObservedSequenceMutationProfiles, label: str, count_attribute_name="ambiguous_accepted_umi_noncollapsed_mutations", unlinked_mutation_attribute_name = "all_observed_surrogate_unlinked_mutations_df", possible_refs: List[str] = ["A","C","G","T"], possible_alts: List[str] = ["A","C","G","T","X","N","-","X"], empty_context_character: str = "_", min_position: int = 0, max_position: Optional[int] = None, filename: Optional[str] = None):
    fig, ax = plt.subplots(1, figsize=(15, 5*1))

    mutation_results_typed = getattr(mutations_results, count_attribute_name)
    if mutation_results_typed is not None:
        unlinked_mutations_df: pd.DataFrame = getattr(mutation_results_typed, unlinked_mutation_attribute_name)

        if (unlinked_mutations_df is not None) and (unlinked_mutations_df.shape[0] > 0): # If there are mutations
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