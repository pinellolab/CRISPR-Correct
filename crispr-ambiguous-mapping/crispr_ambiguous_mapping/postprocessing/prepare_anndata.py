import pandas as pd
import anndata as ad
from typing import Optional, List
import numpy as np

def initialize_anndata(sample_info_df: pd.DataFrame, 
                       guide_library_df: pd.DataFrame,
                       count_df: Optional[pd.DataFrame] = None, count_series_list: Optional[List[pd.Series]] = None) -> ad.AnnData:
    # Prepare sample_info_df
    sample_info_df = sample_info_df.reset_index(drop=True)

    # Prepare count_df
    if count_df is None:
        if count_series_list is not None:
            count_df = pd.DataFrame(count_series_list).transpose()
            count_df.columns = sample_info_df.index
        else:
            raise Exception("Must provide processed count_df or unprocessed count_series_list")
    
    # Set guide_library_df index
    guide_columns = ['protospacer', 'surrogate', 'barcode']
    available_columns = [col for col in guide_columns if col in guide_library_df.columns]
    guide_library_df.index = guide_library_df["protospacer"] + "_" +  guide_library_df["surrogate"] + "_" + guide_library_df["barcode"]
    guide_library_df.index = guide_library_df[available_columns].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)

    # Set count_df indices
    count_df_index = count_df.index.to_frame()
    count_df.index = count_df_index[available_columns].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
    count_df = count_df.loc[guide_library_df.index] # Filter count_df from available guides
    count_df.columns = sample_info_df.index

    # Initialize AnnData
    count_anndata = ad.AnnData(var=guide_library_df, obs=sample_info_df, X=count_df.transpose())

    return count_anndata

def add_lfc_layer(anndata: ad.AnnData, 
                  groupby_columns_list: List[str], 
                  enriched_population_label: str, 
                  baseline_population_label: str, 
                  population_column_name: str) -> ad.AnnData:
    lfc_dict = {} 
    for group in anndata.obs.groupby(groupby_columns_list):
        enriched_pop_sample_df = group[1][group[1][population_column_name] == enriched_population_label]
        assert enriched_pop_sample_df.shape[0] == 1, f"Multiple samples matching enriched condition in group {group[1]}"
        enriched_pop_sample = enriched_pop_sample_df.iloc[0,:]
        
        baseline_pop_sample_df = group[1][group[1][population_column_name] == baseline_population_label]
        assert baseline_pop_sample_df.shape[0] == 1, f"Multiple samples matching enriched condition in group {group[1]}"
        baseline_pop_sample = baseline_pop_sample_df.iloc[0,:]

        enriched_anndata = anndata[int(enriched_pop_sample.name), :] # DEVELOPER NOTE: Will thow error of indices of samples are changed from ints
        baseline_anndata = anndata[int(baseline_pop_sample.name), :] # DEVELOPER NOTE: Will thow error of indices of samples are changed from ints

        enriched_X = enriched_anndata.X[0]
        baseline_X = baseline_anndata.X[0]

        enriched_X_safetargeting_total = enriched_anndata[:, enriched_anndata.var["type"] == "safe-targeting"].X[0].sum()
        baseline_X_safetargeting_total = baseline_anndata[:, baseline_anndata.var["type"] == "safe-targeting"].X[0].sum()

        baseline_normalization_constant = enriched_X_safetargeting_total / baseline_X_safetargeting_total

        baseline_X_normalized = baseline_X * baseline_normalization_constant

        enriched_lfc = np.log(np.divide(enriched_X + 1, baseline_X_normalized + 1))
        baseline_lfc = np.log(np.divide(baseline_X_normalized + 1, enriched_X + 1))

        lfc_dict.update({int(enriched_pop_sample.name): enriched_lfc})
        lfc_dict.update({int(baseline_pop_sample.name): baseline_lfc})
    
    # NOTE: LFC layer for all samples must be added at once
    anndata.layers["lfc"] = np.asarray([lfc_dict[idx] for idx in range(0, anndata.shape[0])])
    return anndata