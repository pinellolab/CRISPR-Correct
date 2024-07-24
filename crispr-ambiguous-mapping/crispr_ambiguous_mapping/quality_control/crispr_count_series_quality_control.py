from ..models.mapping_models import MatchSetWhitelistReporterCounterSeriesResults
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from typing import List

def calculate_gini_index(count_series: pd.Series, eps=1e-8):
    '''
        Code from https://gist.github.com/fclesio/dfde0c512a0af55ef3c3dbdbec4bd445
    '''
    # Convert series to numpy array
    arr = count_series.to_numpy()

    # All values are treated equally, arrays must be 1d and > 0:
    arr = np.abs(arr).flatten() + eps

    # Values must be sorted:
    arr = np.sort(arr)

    # Index per array element:
    index = np.arange(1, arr.shape[0]+1)

    # Number of array elements:
    N = arr.shape[0]

    # Gini coefficient:
    return(np.sum((2*index - N - 1)*arr))/(N*np.sum(arr))

def log_count_series_warnings(match_type_string: str,
         count_type_string: str,
         protospacer_match_gt_given_match_type_total: Optional[int],
         given_match_type_vs_protospacer_match_total_ratio: Optional[float],
         given_match_type_vs_protospacer_match_total_ratio_threshold: Optional[float],
         zerocount: int,
         gini: float,
         std: float,
         gini_index_threshold: float,
         coverage_lower_quantile: float,
         lower_quantile_coverage_percentage: float,
         lower_quantile_coverage_threshold: float,
         coverage_mean: float,
         mean_coverage_threshold: int,
         umidups_mean: Optional[float],
         umidups_lower_threshold: Optional[float],
         umidups_upper_threshold: Optional[float]):
    
    if protospacer_match_gt_given_match_type_total is not None:
        # Check that protospacer-only count is greater than protospacer+(surrogate)+(barcode) count
        if protospacer_match_gt_given_match_type_total > 0:
            print(f"WARNING for {match_type_string}; {count_type_string}: There are {protospacer_match_gt_given_match_type_total} unexpected sequences greater than the protospacer-match. This may be due to the match type having improved mapping performance, though these sequences many be worth investigating manually.")
    if (given_match_type_vs_protospacer_match_total_ratio is not None) and (given_match_type_vs_protospacer_match_total_ratio_threshold is not None):
        # Check that the protospacer-only to protospacer+(surrogate)+(barcode) count proportion is not too low
        if given_match_type_vs_protospacer_match_total_ratio <= given_match_type_vs_protospacer_match_total_ratio_threshold:
            print(f"WARNING for {match_type_string}; {count_type_string}: The match type has a much less proportion of coverage (of {given_match_type_vs_protospacer_match_total_ratio}) compared to the protospacer-match. This may be due to high recombination or mapping issues with the given match type. Manually check if this persists across all match types and samples.")
    if zerocount > 0: 
        # Check that there are no guides with a count of zero.
        print(f"WARNING for {match_type_string}; {count_type_string}: There are {zerocount} guides with count of 0. This may be due to low coverage or a mapping issue. Recommended to investigate coverage issue or investigate these sequences manually for mapping issues, especially if recurrent across different samples (especially high-coverage plasmid samples).")
    if gini >= gini_index_threshold:
        # Check that the Gini index is not too high
        print(f"WARNING for {match_type_string}; {count_type_string}: Guide counts have a high gini index of {gini} greather than threshold {gini_index_threshold} with standard deviation {std}. There may be many guides with low coverage that may be underpowered for statistical testing, even if the mean coverage is high. This may be due to a bottleneck during the screen. Recommended to compare the gini index with the plasmid sample to see how much the non-uniformity increases during the screen.")
    if coverage_mean <= mean_coverage_threshold:
        # Check that the mean coverage is not too low
        print(f"WARNING: Mean coverage {coverage_mean} is below specified threshold {mean_coverage_threshold}.")
    if coverage_lower_quantile <= lower_quantile_coverage_threshold:
        # Check that the lower quantile coverage is not too low
        print(f"WARNING for {match_type_string}; {count_type_string}: Lower {lower_quantile_coverage_percentage} quantile of coverage {coverage_lower_quantile} is below specified threshold {lower_quantile_coverage_threshold}. These guides may be underpowered for statistical testing unless the effect sizes are large.")
    
    if (umidups_mean is not None) and (umidups_lower_threshold is not None):
        # Check that the UMI duplication rate is not too low
        if umidups_mean <= umidups_threshold:
            print(f"WARNING for {match_type_string}; {count_type_string}: Average UMI duplication rate of {umidups_mean} is below specified threshold of {umidups_lower_threshold}, refer to the UMI duplication histogram plotted. Your sample may be undersequenced as low average UMI duplication suggests spare material unsequenced. If guide coverage is low, it is suggested to sequence more material to increase coverage. It is also suggested to use the UMI-noncollapsed count for downstream analysis")
            
            # If the guide coverage is also low, recommend sequencing more
            if (coverage_mean <= mean_coverage_threshold) or (coverage_lower_quantile <= lower_quantile_coverage_threshold):
                print(f"WARNING for {match_type_string}; {count_type_string}: Since both the average UMI duplication rate ({umidups_mean}) and guide coverage of {coverage_mean} (with lower quantile coverage of {coverage_lower_quantile}) is low, it is suggested to sequence more material to increase guide coverage.")

    if (umidups_mean is not None) and (umidups_upper_threshold is not None):
        # Check that the UMI duplication rate is not too high.
        if umidups_mean >= umidups_threshold:
            print(f"NOTE for {match_type_string}; {count_type_string}: Average UMI duplication rate of {umidups_mean} is higher than specified threshold of {umidups_upper_threshold}, refer to the UMI duplication histogram plotted. Your sample may be oversequenced as high average UMI duplication suggests material repetitively sequence post-PCR amplification. Suggested to use the UMI-collapsed count for downstream analysis.")

            # If the guide coverage is also low, recommend performing another screen with higher cell coverage to gain sufficient power
            if (coverage_mean <= mean_coverage_threshold) or (coverage_lower_quantile <= lower_quantile_coverage_threshold):
                print(f"WARNING for {match_type_string}; {count_type_string}: Since the average UMI duplication rate ({umidups_mean}) is high but the guide coverage of {coverage_mean} (with lower quantile coverage of {coverage_lower_quantile}) is low, sequencing more material may only increase UMI duplication but not the UMI-collapsed guide coverage. It is suggested to perform another screen with increased cell coverage to ensure sufficient power.")
    

def count_series_result_quality_control_stats(count_series_result: MatchSetWhitelistReporterCounterSeriesResults,
                                              contains_umi: bool, 
                                              contains_surrogate: bool, 
                                              contains_barcode: bool, 
                                              display_visualizations: bool = True,
                                              gini_index_threshold: float = 0.4,
                                              mean_coverage_threshold: float = 100,
                                              lower_quantile_coverage_percentage: float = 0.1,
                                              lower_quantile_coverage_threshold: float = 60,
                                              given_match_type_vs_protospacer_match_total_ratio_threshold= 0.3,
                                              umidups_lower_threshold: int = 2,
                                              umidups_upper_threshold: int = 4):

    
    count_series_result_quality_control_stats_dict = dict()
    if contains_umi:

        gini_index_bar_plot_labels: List[str] = []
        collapsed_gini_index_bar_plot_values: List[float] = []
        noncollapsed_gini_index_bar_plot_values: List[float] = []
        
        # PROTOSPACER - COLLAPSED
        protospacer_match_collapsed_zerocount = sum(count_series_result.protospacer_match.ambiguous_accepted_umi_collapsed_counterseries == 0)
        protospacer_match_collapsed_totalcount = len(count_series_result.protospacer_match.ambiguous_accepted_umi_collapsed_counterseries)
        protospacer_match_collapsed_gini = calculate_gini_index(count_series_result.protospacer_match.ambiguous_accepted_umi_collapsed_counterseries)
        protospacer_match_collapsed_coverage_mean = np.mean(count_series_result.protospacer_match.ambiguous_accepted_umi_collapsed_counterseries)
        protospacer_match_collapsed_coverage_std = np.std(count_series_result.protospacer_match.ambiguous_accepted_umi_collapsed_counterseries)
        protospacer_match_collapsed_coverage_lower_quantile = count_series_result.protospacer_match.ambiguous_accepted_umi_collapsed_counterseries.quantile(q=low_quantile_coverage_percentage)
        
        count_series_result_quality_control_stats_dict.update({
            "protospacer_match_collapsed_totalcount": protospacer_match_collapsed_totalcount,
            "protospacer_match_collapsed_zerocount": protospacer_match_collapsed_zerocount,
            "protospacer_match_collapsed_gini": protospacer_match_collapsed_gini,
            "protospacer_match_collapsed_coverage_mean": protospacer_match_collapsed_coverage_mean,
            "protospacer_match_collapsed_coverage_std": protospacer_match_collapsed_coverage_std,
            f"protospacer_match_collapsed_coverage_lower_quantile_{lower_quantile_coverage_percentage}": protospacer_match_collapsed_coverage_lower_quantile
        })

        if display_visualizations:
            # Add values for plots
            gini_index_bar_plot_labels.append("PROTOSPACER")
            collapsed_gini_index_bar_plot_values.append(protospacer_match_collapsed_gini)

            # Log warnings from count statistics
            log_count_series_warnings(match_type_string="protospacer-match",
                                      count_type_string="UMI-collapsed",
                                      protospacer_match_gt_given_match_type_total=None,
                                      given_match_type_vs_protospacer_match_total_ratio=None,
                                      given_match_type_vs_protospacer_match_total_ratio_threshold=given_match_type_vs_protospacer_match_total_ratio_threshold,
                                      zerocount=protospacer_match_collapsed_zerocount,
                                      gini=protospacer_match_collapsed_gini,
                                      std=protospacer_match_collapsed_coverage_std,
                                      gini_index_threshold=gini_index_threshold,
                                      coverage_lower_quantile=protospacer_match_collapsed_coverage_lower_quantile,
                                      lower_quantile_coverage_percentage=lower_quantile_coverage_percentage,
                                      lower_quantile_coverage_threshold=lower_quantile_coverage_threshold,
                                      coverage_mean=protospacer_match_collapsed_coverage_mean,
                                      mean_coverage_threshold=mean_coverage_threshold)

        # PROTOSPACER - NONCOLLAPSED
        protospacer_match_noncollapsed_zerocount = sum(count_series_result.protospacer_match.ambiguous_accepted_umi_noncollapsed_counterseries == 0)
        protospacer_match_noncollapsed_totalcount = len(count_series_result.protospacer_match.ambiguous_accepted_umi_noncollapsed_counterseries)
        protospacer_match_noncollapsed_gini = calculate_gini_index(count_series_result.protospacer_match.ambiguous_accepted_umi_noncollapsed_counterseries)
        protospacer_match_noncollapsed_coverage_mean = np.mean(count_series_result.protospacer_match.ambiguous_accepted_umi_noncollapsed_counterseries)
        protospacer_match_noncollapsed_coverage_std = np.std(count_series_result.protospacer_match.ambiguous_accepted_umi_noncollapsed_counterseries)
        protospacer_match_noncollapsed_coverage_lower_quantile = count_series_result.protospacer_match.ambiguous_accepted_umi_noncollapsed_counterseries.quantile(q=low_quantile_coverage_percentage)

        protospacer_match_umidups = count_series_result.protospacer_match.ambiguous_accepted_umi_noncollapsed_counterseries / count_series_result.protospacer_match.ambiguous_accepted_umi_collapsed_counterseries
        protospacer_match_umidups_mean = protospacer_match_umidups.mean()

        count_series_result_quality_control_stats_dict.update({
            "protospacer_match_noncollapsed_totalcount": protospacer_match_noncollapsed_totalcount,
            "protospacer_match_noncollapsed_zerocount": protospacer_match_noncollapsed_zerocount,
            "protospacer_match_noncollapsed_gini": protospacer_match_noncollapsed_gini,
            "protospacer_match_noncollapsed_coverage_mean": protospacer_match_noncollapsed_coverage_mean,
            "protospacer_match_noncollapsed_coverage_std": protospacer_match_noncollapsed_coverage_std,
            f"protospacer_match_noncollapsed_coverage_lower_quantile_{lower_quantile_coverage_percentage}": protospacer_match_noncollapsed_coverage_lower_quantile,
            "protospacer_match_umidups_mean": protospacer_match_umidups_mean
        })

        if display_visualizations:
            # Add values for plots
            noncollapsed_gini_index_bar_plot_values.append(protospacer_match_noncollapsed_gini)

            # Log warnings from count statistics
            log_count_series_warnings(match_type_string="protospacer-match",
                                      count_type_string="UMI-noncollapsed",
                                      protospacer_match_gt_given_match_type_total=None,
                                      given_match_type_vs_protospacer_match_total_ratio=None,
                                      given_match_type_vs_protospacer_match_total_ratio_threshold=given_match_type_vs_protospacer_match_total_ratio_threshold,
                                      zerocount=protospacer_match_noncollapsed_zerocount,
                                      gini=protospacer_match_noncollapsed_gini,
                                      std=protospacer_match_noncollapsed_coverage_std,
                                      gini_index_threshold=gini_index_threshold,
                                      coverage_lower_quantile=protospacer_match_collapsed_coverage_lower_quantile,
                                      lower_quantile_coverage_percentage=lower_quantile_coverage_percentage,
                                      lower_quantile_coverage_threshold=lower_quantile_coverage_threshold,
                                      coverage_mean=protospacer_match_noncollapsed_coverage_mean,
                                      mean_coverage_threshold=mean_coverage_threshold,
                                      umidups_mean=protospacer_match_umidups_mean,
                                      umidups_lower_threshold=umidups_lower_threshold,
                                      umidups_upper_threshold=umidups_upper_threshold)


        if contains_surrogate:
                # PROTOSPACER+SURROGATE vs. PROTOSPACER - COLLAPSED
                protospacer_match_gt_protospacer_surrogate_match_collapsed_total = sum(count_series_result.protospacer_match.ambiguous_accepted_umi_collapsed_counterseries >= count_series_result.protospacer_match_surrogate_match.ambiguous_accepted_umi_collapsed_counterseries)
                protospacer_surrogate_match_vs_protospacer_match_collapsed_total_ratio_threshold = sum(count_series_result.protospacer_match_surrogate_match.ambiguous_accepted_umi_collapsed_counterseries) / sum(count_series_result.protospacer_match.ambiguous_accepted_umi_collapsed_counterseries)
                protospacer_surrogate_match_collapsed_zerocount = sum(count_series_result.protospacer_match_surrogate_match.ambiguous_accepted_umi_collapsed_counterseries == 0)
                protospacer_surrogate_match_collapsed_totalcount = len(count_series_result.protospacer_match_surrogate_match.ambiguous_accepted_umi_collapsed_counterseries)
                protospacer_surrogate_match_collapsed_gini = calculate_gini_index(count_series_result.protospacer_match_surrogate_match.ambiguous_accepted_umi_collapsed_counterseries)
                protospacer_surrogate_match_collapsed_coverage_mean = np.mean(count_series_result.protospacer_match_surrogate_match.ambiguous_accepted_umi_collapsed_counterseries)
                protospacer_surrogate_match_collapsed_coverage_std = np.std(count_series_result.protospacer_match_surrogate_match.ambiguous_accepted_umi_collapsed_counterseries)
                protospacer_surrogate_match_collapsed_coverage_lower_quantile = count_series_result.protospacer_match_surrogate_match.ambiguous_accepted_umi_collapsed_counterseries.quantile(q=low_quantile_coverage_percentage)

                count_series_result_quality_control_stats_dict.update({
                    "protospacer_surrogate_match_collapsed_totalcount": protospacer_surrogate_match_collapsed_totalcount,
                    "protospacer_match_gt_protospacer_surrogate_match_collapsed_total": protospacer_match_gt_protospacer_surrogate_match_collapsed_total,
                    "protospacer_surrogate_match_vs_protospacer_match_collapsed_total_ratio_threshold": protospacer_surrogate_match_vs_protospacer_match_collapsed_total_ratio_threshold,
                    "protospacer_surrogate_match_collapsed_zerocount": protospacer_surrogate_match_collapsed_zerocount,
                    "protospacer_surrogate_match_collapsed_gini": protospacer_surrogate_match_collapsed_gini,
                    "protospacer_surrogate_match_collapsed_coverage_mean": protospacer_surrogate_match_collapsed_coverage_mean,
                    "protospacer_surrogate_match_collapsed_coverage_std": protospacer_surrogate_match_collapsed_coverage_std,
                    f"protospacer_surrogate_match_collapsed_coverage_lower_quantile_{lower_quantile_coverage_percentage}": protospacer_surrogate_match_collapsed_coverage_lower_quantile
                })
                if display_visualizations:
                    gini_index_bar_plot_labels.append("PROTOSPACER+SURROGATE")
                    collapsed_gini_index_bar_plot_values.append(protospacer_surrogate_match_collapsed_gini)

                    # Log warnings from count statistics
                    log_count_series_warnings(match_type_string="protospacer+surrogate-match",
                                      count_type_string="UMI-collapsed",
                                      protospacer_match_gt_given_match_type_total=protospacer_match_gt_protospacer_surrogate_match_collapsed_total,
                                      given_match_type_vs_protospacer_match_total_ratio=protospacer_surrogate_match_vs_protospacer_match_collapsed_total_ratio_threshold,
                                      given_match_type_vs_protospacer_match_total_ratio_threshold=given_match_type_vs_protospacer_match_total_ratio_threshold,
                                      zerocount=protospacer_surrogate_match_collapsed_zerocount,
                                      gini=protospacer_surrogate_match_collapsed_gini,
                                      std=protospacer_surrogate_match_collapsed_coverage_std,
                                      gini_index_threshold=gini_index_threshold,
                                      coverage_lower_quantile=protospacer_surrogate_match_collapsed_coverage_lower_quantile,
                                      lower_quantile_coverage_percentage=lower_quantile_coverage_percentage,
                                      lower_quantile_coverage_threshold=lower_quantile_coverage_threshold,
                                      coverage_mean=protospacer_surrogate_match_collapsed_coverage_mean,
                                      mean_coverage_threshold=mean_coverage_threshold)

                # PROTOSPACER+SURROGATE vs. PROTOSPACER - NONCOLLAPSED
                protospacer_match_gt_protospacer_surrogate_match_noncollapsed_total = sum(count_series_result.protospacer_match.ambiguous_accepted_umi_noncollapsed_counterseries >= count_series_result.protospacer_match_surrogate_match.ambiguous_accepted_umi_noncollapsed_counterseries)
                protospacer_surrogate_match_vs_protospacer_match_noncollapsed_total_ratio_threshold = sum(count_series_result.protospacer_match_surrogate_match.ambiguous_accepted_umi_noncollapsed_counterseries) / sum(count_series_result.protospacer_match.ambiguous_accepted_umi_noncollapsed_counterseries)
                protospacer_surrogate_match_noncollapsed_zerocount = sum(count_series_result.protospacer_match_surrogate_match.ambiguous_accepted_umi_noncollapsed_counterseries == 0)
                protospacer_surrogate_match_noncollapsed_totalcount = len(count_series_result.protospacer_match_surrogate_match.ambiguous_accepted_umi_noncollapsed_counterseries)
                protospacer_surrogate_match_noncollapsed_gini = calculate_gini_index(count_series_result.protospacer_match_surrogate_match.ambiguous_accepted_umi_noncollapsed_counterseries)
                protospacer_surrogate_match_noncollapsed_coverage_mean = np.mean(count_series_result.protospacer_match_surrogate_match.ambiguous_accepted_umi_noncollapsed_counterseries)
                protospacer_surrogate_match_noncollapsed_coverage_std = np.std(count_series_result.protospacer_match_surrogate_match.ambiguous_accepted_umi_noncollapsed_counterseries)
                protospacer_surrogate_match_noncollapsed_coverage_lower_quantile = count_series_result.protospacer_match_surrogate_match.ambiguous_accepted_umi_noncollapsed_counterseries.quantile(q=low_quantile_coverage_percentage)

                protospacer_surrogate_match_umidups = count_series_result.protospacer_match_surrogate_match.ambiguous_accepted_umi_noncollapsed_counterseries / count_series_result.protospacer_match_surrogate_match.ambiguous_accepted_umi_collapsed_counterseries
                protospacer_surrogate_match_umidups_mean = protospacer_surrogate_match_umidups.mean()

                count_series_result_quality_control_stats_dict.update({
                    "protospacer_surrogate_match_noncollapsed_totalcount": protospacer_surrogate_match_noncollapsed_totalcount,
                    "protospacer_match_gt_protospacer_surrogate_match_noncollapsed_total": protospacer_match_gt_protospacer_surrogate_match_noncollapsed_total,
                    "protospacer_surrogate_match_vs_protospacer_match_noncollapsed_total_ratio_threshold": protospacer_surrogate_match_vs_protospacer_match_noncollapsed_total_ratio_threshold,
                    "protospacer_surrogate_match_noncollapsed_zerocount": protospacer_surrogate_match_noncollapsed_zerocount,
                    "protospacer_surrogate_match_noncollapsed_gini": protospacer_surrogate_match_noncollapsed_gini,
                    "protospacer_surrogate_match_noncollapsed_coverage_mean": protospacer_surrogate_match_noncollapsed_coverage_mean,
                    "protospacer_surrogate_match_noncollapsed_coverage_std": protospacer_surrogate_match_noncollapsed_coverage_std,
                    f"protospacer_surrogate_match_noncollapsed_coverage_lower_quantile_{lower_quantile_coverage_percentage}": protospacer_surrogate_match_noncollapsed_coverage_lower_quantile,
                    "protospacer_surrogate_match_umidups_mean": protospacer_surrogate_match_umidups_mean
                })
                
                if display_visualizations:
                    noncollapsed_gini_index_bar_plot_values.append(protospacer_surrogate_match_noncollapsed_gini)

                    # Log warnings from count statistics
                    log_count_series_warnings(match_type_string="protospacer+surrogate-match",
                                      count_type_string="UMI-noncollapsed",
                                      protospacer_match_gt_given_match_type_total=protospacer_match_gt_protospacer_surrogate_match_noncollapsed_total,
                                      given_match_type_vs_protospacer_match_total_ratio=protospacer_surrogate_match_vs_protospacer_match_noncollapsed_total_ratio_threshold,
                                      given_match_type_vs_protospacer_match_total_ratio_threshold=given_match_type_vs_protospacer_match_total_ratio_threshold,
                                      zerocount=protospacer_surrogate_match_noncollapsed_zerocount,
                                      gini=protospacer_surrogate_match_noncollapsed_gini,
                                      std=protospacer_surrogate_match_noncollapsed_coverage_std,
                                      gini_index_threshold=gini_index_threshold,
                                      coverage_lower_quantile=protospacer_surrogate_match_noncollapsed_coverage_lower_quantile,
                                      lower_quantile_coverage_percentage=lower_quantile_coverage_percentage,
                                      lower_quantile_coverage_threshold=lower_quantile_coverage_threshold,
                                      coverage_mean=protospacer_surrogate_match_noncollapsed_coverage_mean,
                                      mean_coverage_threshold=mean_coverage_threshold,
                                      umidups_mean=protospacer_surrogate_match_umidups_mean,
                                      umidups_lower_threshold=umidups_lower_threshold,
                                      umidups_upper_threshold=umidups_upper_threshold
                                      )
                

                if contains_barcode:
                        # PROTOSPACER+SURROGATE+BARCODE vs. PROTOSPACER - COLLAPSED
                        protospacer_match_gt_protospacer_surrogate_barcode_match_collapsed_total = sum(count_series_result.protospacer_match.ambiguous_accepted_umi_collapsed_counterseries >= count_series_result.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_umi_collapsed_counterseries)
                        protospacer_surrogate_barcode_match_vs_protospacer_match_collapsed_total_ratio_threshold = sum(count_series_result.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_umi_collapsed_counterseries) / sum(count_series_result.protospacer_match.ambiguous_accepted_umi_collapsed_counterseries)
                        protospacer_surrogate_barcode_match_collapsed_zerocount = sum(count_series_result.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_umi_collapsed_counterseries == 0)
                        protospacer_surrogate_barcode_match_collapsed_totalcount = len(count_series_result.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_umi_collapsed_counterseries)
                        protospacer_surrogate_barcode_match_collapsed_gini = calculate_gini_index(count_series_result.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_umi_collapsed_counterseries)
                        protospacer_surrogate_barcode_match_collapsed_coverage_mean = np.mean(count_series_result.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_umi_collapsed_counterseries)
                        protospacer_surrogate_barcode_match_collapsed_coverage_std = np.std(count_series_result.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_umi_collapsed_counterseries)
                        protospacer_surrogate_barcode_match_collapsed_coverage_lower_quantile = count_series_result.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_umi_collapsed_counterseries.quantile(q=low_quantile_coverage_percentage)

                        count_series_result_quality_control_stats_dict.update({
                            "protospacer_surrogate_barcode_match_collapsed_totalcount": protospacer_surrogate_barcode_match_collapsed_totalcount,
                            "protospacer_match_gt_protospacer_surrogate_barcode_match_collapsed_total": protospacer_match_gt_protospacer_surrogate_barcode_match_collapsed_total,
                            "protospacer_surrogate_barcode_match_vs_protospacer_match_collapsed_total_ratio_threshold": protospacer_surrogate_barcode_match_vs_protospacer_match_collapsed_total_ratio_threshold,
                            "protospacer_surrogate_barcode_match_collapsed_zerocount": protospacer_surrogate_barcode_match_collapsed_zerocount,
                            "protospacer_surrogate_barcode_match_collapsed_gini": protospacer_surrogate_barcode_match_collapsed_gini,
                            "protospacer_surrogate_barcode_match_collapsed_coverage_mean": protospacer_surrogate_barcode_match_collapsed_coverage_mean,
                            "protospacer_surrogate_barcode_match_collapsed_coverage_std": protospacer_surrogate_barcode_match_collapsed_coverage_std,
                            f"protospacer_surrogate_barcode_match_collapsed_coverage_lower_quantile_{lower_quantile_coverage_percentage}": protospacer_surrogate_barcode_match_collapsed_coverage_lower_quantile
                        })
                        
                        if display_visualizations:
                            gini_index_bar_plot_labels.append("PROTOSPACER+SURROGATE+BARCODE")
                            collapsed_gini_index_bar_plot_values.append(protospacer_surrogate_barcode_match_collapsed_gini)

                            # Log warnings from count statistics
                            log_count_series_warnings(match_type_string="protospacer+surrogate+barcode-match",
                                            count_type_string="UMI-collapsed",
                                            protospacer_match_gt_given_match_type_total=protospacer_match_gt_protospacer_surrogate_barcode_match_collapsed_total,
                                            given_match_type_vs_protospacer_match_total_ratio=protospacer_surrogate_barcode_match_vs_protospacer_match_collapsed_total_ratio_threshold,
                                            given_match_type_vs_protospacer_match_total_ratio_threshold=given_match_type_vs_protospacer_match_total_ratio_threshold,
                                            zerocount=protospacer_surrogate_barcode_match_collapsed_zerocount,
                                            gini=protospacer_surrogate_barcode_match_collapsed_gini,
                                            std=protospacer_surrogate_barcode_match_collapsed_coverage_std,
                                            gini_index_threshold=gini_index_threshold,
                                            coverage_lower_quantile=protospacer_surrogate_barcode_match_collapsed_coverage_lower_quantile,
                                            lower_quantile_coverage_percentage=lower_quantile_coverage_percentage,
                                            lower_quantile_coverage_threshold=lower_quantile_coverage_threshold,
                                            coverage_mean=protospacer_surrogate_barcode_match_collapsed_coverage_mean,
                                            mean_coverage_threshold=mean_coverage_threshold)

                        # PROTOSPACER+SURROGATE+BARCODE vs. PROTOSPACER - NONCOLLAPSED
                        protospacer_match_gt_protospacer_surrogate_barcode_match_noncollapsed_total = sum(count_series_result.protospacer_match.ambiguous_accepted_umi_noncollapsed_counterseries >= count_series_result.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_umi_noncollapsed_counterseries)
                        protospacer_surrogate_barcode_match_vs_protospacer_match_noncollapsed_total_ratio_threshold = sum(count_series_result.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_umi_noncollapsed_counterseries) / sum(count_series_result.protospacer_match.ambiguous_accepted_umi_noncollapsed_counterseries)
                        protospacer_surrogate_barcode_match_noncollapsed_zerocount = sum(count_series_result.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_umi_noncollapsed_counterseries == 0)
                        protospacer_surrogate_barcode_match_noncollapsed_totalcount = len(count_series_result.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_umi_noncollapsed_counterseries)
                        protospacer_surrogate_barcode_match_noncollapsed_gini = calculate_gini_index(count_series_result.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_umi_noncollapsed_counterseries)
                        protospacer_surrogate_barcode_match_noncollapsed_coverage_mean = np.mean(count_series_result.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_umi_noncollapsed_counterseries)
                        protospacer_surrogate_barcode_match_noncollapsed_coverage_std = np.std(count_series_result.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_umi_noncollapsed_counterseries)
                        protospacer_surrogate_barcode_match_noncollapsed_coverage_lower_quantile = count_series_result.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_umi_noncollapsed_counterseries.quantile(q=low_quantile_coverage_percentage)
                        
                        protospacer_surrogate_barcode_match_umidups = count_series_result.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_umi_noncollapsed_counterseries / count_series_result.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_umi_collapsed_counterseries
                        protospacer_surrogate_barcode_match_umidups_mean = protospacer_surrogate_barcode_match_umidups.mean()

                        count_series_result_quality_control_stats_dict.update({
                            "protospacer_surrogate_barcode_match_noncollapsed_totalcount": protospacer_surrogate_barcode_match_noncollapsed_totalcount,
                            "protospacer_match_gt_protospacer_surrogate_barcode_match_noncollapsed_total": protospacer_match_gt_protospacer_surrogate_barcode_match_noncollapsed_total,
                            "protospacer_surrogate_barcode_match_noncollapsed_zerocount": protospacer_surrogate_barcode_match_noncollapsed_zerocount,
                            "protospacer_surrogate_barcode_match_noncollapsed_gini": protospacer_surrogate_barcode_match_noncollapsed_gini,
                            "protospacer_surrogate_barcode_match_noncollapsed_coverage_mean": protospacer_surrogate_barcode_match_noncollapsed_coverage_mean,
                            "protospacer_surrogate_barcode_match_noncollapsed_coverage_std": protospacer_surrogate_barcode_match_noncollapsed_coverage_std,
                            f"protospacer_surrogate_barcode_match_noncollapsed_coverage_lower_quantile_{lower_quantile_coverage_percentage}": protospacer_surrogate_barcode_match_noncollapsed_coverage_lower_quantile,
                            "protospacer_surrogate_barcode_match_umidups": protospacer_surrogate_barcode_match_umidups
                        })

                        if display_visualizations:
                            noncollapsed_gini_index_bar_plot_values.append(protospacer_surrogate_barcode_match_noncollapsed_gini)

                            # Log warnings from count statistics
                            log_count_series_warnings(match_type_string="protospacer+surrogate+barcode-match",
                                            count_type_string="UMI-noncollapsed",
                                            protospacer_match_gt_given_match_type_total=,
                                            given_match_type_vs_protospacer_match_total_ratio=,
                                            given_match_type_vs_protospacer_match_total_ratio_threshold=given_match_type_vs_protospacer_match_total_ratio_threshold,
                                            zerocount=protospacer_surrogate_barcode_match_noncollapsed_zerocount,
                                            gini=protospacer_surrogate_barcode_match_noncollapsed_gini,
                                            std=protospacer_surrogate_barcode_match_noncollapsed_coverage_std,
                                            gini_index_threshold=protospacer_surrogate_barcode_match_noncollapsed_gini,
                                            coverage_lower_quantile=protospacer_surrogate_barcode_match_collapsed_coverage_lower_quantile,
                                            lower_quantile_coverage_percentage=lower_quantile_coverage_percentage,
                                            lower_quantile_coverage_threshold=lower_quantile_coverage_threshold,
                                            coverage_mean=protospacer_surrogate_barcode_match_noncollapsed_coverage_mean,
                                            mean_coverage_threshold=mean_coverage_threshold,
                                            umidups_mean=protospacer_surrogate_barcode_match_umidups,
                                            umidups_lower_threshold=umidups_lower_threshold,
                                            umidups_upper_threshold=umidups_upper_threshold)

        else:
            if contains_barcode:
                    # PROTOSPACER+BARCODE vs. PROTOSPACER - COLLAPSED
                    protospacer_match_gt_protospacer_barcode_match_collapsed_total = sum(count_series_result.protospacer_match.ambiguous_accepted_umi_collapsed_counterseries >= count_series_result.protospacer_match_barcode_match.ambiguous_accepted_umi_collapsed_counterseries)
                    protospacer_barcode_match_vs_protospacer_match_collapsed_total_ratio_threshold = sum(count_series_result.protospacer_match_barcode_match.ambiguous_accepted_umi_collapsed_counterseries) / sum(count_series_result.protospacer_match.ambiguous_accepted_umi_collapsed_counterseries)
                    protospacer_barcode_match_collapsed_zerocount = sum(count_series_result.protospacer_match_barcode_match.ambiguous_accepted_umi_collapsed_counterseries == 0)
                    protospacer_barcode_match_collapsed_totalcount = len(count_series_result.protospacer_match_barcode_match.ambiguous_accepted_umi_collapsed_counterseries)
                    protospacer_barcode_match_collapsed_gini = calculate_gini_index(count_series_result.protospacer_match_barcode_match.ambiguous_accepted_umi_collapsed_counterseries)
                    protospacer_barcode_match_collapsed_coverage_mean = np.mean(count_series_result.protospacer_match_barcode_match.ambiguous_accepted_umi_collapsed_counterseries)
                    protospacer_barcode_match_collapsed_coverage_std = np.std(count_series_result.protospacer_match_barcode_match.ambiguous_accepted_umi_collapsed_counterseries)
                    protospacer_barcode_match_collapsed_coverage_lower_quantile = count_series_result.protospacer_match_barcode_match.ambiguous_accepted_umi_collapsed_counterseries.quantile(q=low_quantile_coverage_percentage)

                    count_series_result_quality_control_stats_dict.update({
                        "protospacer_barcode_match_collapsed_totalcount": protospacer_barcode_match_collapsed_totalcount,
                        "protospacer_match_gt_protospacer_barcode_match_collapsed_total": protospacer_match_gt_protospacer_barcode_match_collapsed_total,
                        "protospacer_barcode_match_vs_protospacer_match_collapsed_total_ratio_threshold": protospacer_barcode_match_vs_protospacer_match_collapsed_total_ratio_threshold,
                        "protospacer_barcode_match_collapsed_zerocount": protospacer_barcode_match_collapsed_zerocount,
                        "protospacer_barcode_match_collapsed_gini": protospacer_barcode_match_collapsed_gini,
                        "protospacer_barcode_match_collapsed_coverage_mean": protospacer_barcode_match_collapsed_coverage_mean,
                        "protospacer_barcode_match_collapsed_coverage_std": protospacer_barcode_match_collapsed_coverage_std,
                        f"protospacer_barcode_match_collapsed_coverage_lower_quantile_{lower_quantile_coverage_percentage}": protospacer_barcode_match_collapsed_coverage_lower_quantile
                    })
                    
                    if display_visualizations:
                        gini_index_bar_plot_labels.append("PROTOSPACER+BARCODE")
                        collapsed_gini_index_bar_plot_values.append(protospacer_barcode_match_collapsed_gini)

                        # Log warnings from count statistics
                        log_count_series_warnings(match_type_string="protospacer+barcode-match",
                                        count_type_string="UMI-collapsed",
                                        protospacer_match_gt_given_match_type_total=protospacer_match_gt_protospacer_barcode_match_collapsed_total,
                                        given_match_type_vs_protospacer_match_total_ratio=protospacer_barcode_match_vs_protospacer_match_collapsed_total_ratio_threshold,
                                        given_match_type_vs_protospacer_match_total_ratio_threshold=given_match_type_vs_protospacer_match_total_ratio_threshold,
                                        zerocount=protospacer_barcode_match_collapsed_zerocount,
                                        gini=protospacer_barcode_match_collapsed_gini,
                                        std=protospacer_barcode_match_collapsed_coverage_std,
                                        gini_index_threshold=protospacer_barcode_match_collapsed_gini,
                                        coverage_lower_quantile=protospacer_barcode_match_collapsed_coverage_lower_quantile,
                                        lower_quantile_coverage_percentage=lower_quantile_coverage_percentage,
                                        lower_quantile_coverage_threshold=lower_quantile_coverage_threshold,
                                        coverage_mean=protospacer_barcode_match_collapsed_coverage_mean,
                                        mean_coverage_threshold=mean_coverage_threshold)



                    # PROTOSPACER+BARCODE vs. PROTOSPACER - NONCOLLAPSED
                    protospacer_match_gt_protospacer_barcode_match_noncollapsed_total = sum(count_series_result.protospacer_match.ambiguous_accepted_umi_noncollapsed_counterseries >= count_series_result.protospacer_match_barcode_match.ambiguous_accepted_umi_noncollapsed_counterseries)
                    protospacer_barcode_match_vs_protospacer_match_noncollapsed_total_ratio_threshold = sum(count_series_result.protospacer_match_barcode_match.ambiguous_accepted_umi_noncollapsed_counterseries) / sum(count_series_result.protospacer_match.ambiguous_accepted_umi_noncollapsed_counterseries)
                    protospacer_barcode_match_noncollapsed_zerocount = sum(count_series_result.protospacer_match_barcode_match.ambiguous_accepted_umi_noncollapsed_counterseries == 0)
                    protospacer_barcode_match_noncollapsed_totalcount = len(count_series_result.protospacer_match_barcode_match.ambiguous_accepted_umi_noncollapsed_counterseries)
                    protospacer_barcode_match_noncollapsed_gini = calculate_gini_index(count_series_result.protospacer_match_barcode_match.ambiguous_accepted_umi_noncollapsed_counterseries)
                    protospacer_barcode_match_noncollapsed_coverage_mean = np.mean(count_series_result.protospacer_match_barcode_match.ambiguous_accepted_umi_noncollapsed_counterseries)
                    protospacer_barcode_match_noncollapsed_coverage_std = np.std(count_series_result.protospacer_match_barcode_match.ambiguous_accepted_umi_noncollapsed_counterseries)
                    protospacer_barcode_match_noncollapsed_coverage_lower_quantile = count_series_result.protospacer_match_barcode_match.ambiguous_accepted_umi_noncollapsed_counterseries.quantile(q=low_quantile_coverage_percentage)

                    protospacer_barcode_match_umidups = count_series_result.protospacer_match_barcode_match.ambiguous_accepted_umi_noncollapsed_counterseries / count_series_result.protospacer_match_barcode_match.ambiguous_accepted_umi_collapsed_counterseries
                    protospacer_barcode_match_umidups_mean = protospacer_barcode_match_umidups.mean()

                    count_series_result_quality_control_stats_dict.update({
                        "protospacer_barcode_match_noncollapsed_totalcount": protospacer_barcode_match_noncollapsed_totalcount,
                        "protospacer_match_gt_protospacer_barcode_match_noncollapsed_total": protospacer_match_gt_protospacer_barcode_match_noncollapsed_total,
                        "protospacer_barcode_match_vs_protospacer_match_noncollapsed_total_ratio_threshold": protospacer_barcode_match_vs_protospacer_match_noncollapsed_total_ratio_threshold,
                        "protospacer_barcode_match_noncollapsed_zerocount": protospacer_barcode_match_noncollapsed_zerocount,
                        "protospacer_barcode_match_noncollapsed_gini": protospacer_barcode_match_noncollapsed_gini,
                        "protospacer_barcode_match_noncollapsed_coverage_mean": protospacer_barcode_match_noncollapsed_coverage_mean,
                        "protospacer_barcode_match_noncollapsed_coverage_std": protospacer_barcode_match_noncollapsed_coverage_std,
                        f"protospacer_barcode_match_noncollapsed_coverage_lower_quantile_{lower_quantile_coverage_percentage}": protospacer_barcode_match_noncollapsed_coverage_lower_quantile,
                        "protospacer_barcode_match_umidups": protospacer_barcode_match_umidups
                    })

                    if display_visualizations:
                        noncollapsed_gini_index_bar_plot_values.append(protospacer_barcode_match_noncollapsed_gini)

                        # Log warnings from count statistics
                        log_count_series_warnings(match_type_string="protospacer+barcode-match",
                                        count_type_string="UMI-collapsed",
                                        protospacer_match_gt_given_match_type_total=protospacer_match_gt_protospacer_barcode_match_noncollapsed_total,
                                        given_match_type_vs_protospacer_match_total_ratio=protospacer_barcode_match_vs_protospacer_match_noncollapsed_total_ratio_threshold,
                                        given_match_type_vs_protospacer_match_total_ratio_threshold=given_match_type_vs_protospacer_match_total_ratio_threshold,
                                        zerocount=protospacer_barcode_match_noncollapsed_zerocount,
                                        gini=protospacer_barcode_match_noncollapsed_gini,
                                        std=protospacer_barcode_match_noncollapsed_coverage_std,
                                        gini_index_threshold=gini_index_threshold,
                                        coverage_lower_quantile=protospacer_barcode_match_noncollapsed_coverage_lower_quantile,
                                        lower_quantile_coverage_percentage=lower_quantile_coverage_percentage,
                                        lower_quantile_coverage_threshold=lower_quantile_coverage_threshold,
                                        coverage_mean=protospacer_barcode_match_noncollapsed_coverage_mean,
                                        mean_coverage_threshold=mean_coverage_threshold,
                                        umidups_mean==protospacer_barcode_match_umidups,
                                        umidups_lower_threshold=umidups_lower_threshold,
                                        umidups_upper_threshold=umidups_upper_threshold
                                        )

        # PLOT GINI INDEX:
        if display_visualizations:
            gini_index_bar_plot_input_df = pd.DataFrame({
                "UMI Collapsed": collapsed_gini_index_bar_plot_values,
                "UMI Non-collapsed": noncollapsed_gini_index_bar_plot_values}, index=gini_index_bar_plot_labels
                )
            gini_index_bar_plot_input_df.plot.bar(rot=0)
            plt.title("Gini-Index")
            plt.show()
    else:
        gini_index_bar_plot_labels: List[str] = []
        collapsed_gini_index_bar_plot_values: List[float] = []

        # PROTOSPACER+SURROGATE vs. PROTOSPACER
        protospacer_match_noumi_totalcount = len(count_series_result.protospacer_match.ambiguous_accepted_counterseries)
        protospacer_match_noumi_zerocount = sum(count_series_result.protospacer_match.ambiguous_accepted_counterseries == 0)
        protospacer_match_noumi_gini = calculate_gini_index(count_series_result.protospacer_match.ambiguous_accepted_counterseries)
        protospacer_match_noumi_coverage_mean = np.mean(count_series_result.protospacer_match.ambiguous_accepted_counterseries)
        protospacer_match_noumi_coverage_std = np.std(count_series_result.protospacer_match.ambiguous_accepted_counterseries)
        protospacer_match_noumi_coverage_lower_quantile = count_series_result.protospacer_match.ambiguous_accepted_counterseries.quantile(q=low_quantile_coverage_percentage)

        count_series_result_quality_control_stats_dict.update({
            "protospacer_match_noumi_totalcount": protospacer_match_noumi_totalcount,
            "protospacer_match_noumi_zerocount": protospacer_match_noumi_zerocount,
            "protospacer_match_noumi_gini": protospacer_match_noumi_gini,
            "protospacer_match_noumi_coverage_mean": protospacer_match_noumi_coverage_mean,
            "protospacer_match_noumi_coverage_std": protospacer_match_noumi_coverage_std,
            f"protospacer_match_noumi_coverage_lower_quantile_{lower_quantile_coverage_percentage}": protospacer_match_noumi_coverage_lower_quantile
        })
        if display_visualizations:
            # Log warnings from count statistics
            log_count_series_warnings(match_type_string="protospacer-match",
                            count_type_string="No-UMI",
                            protospacer_match_gt_given_match_type_total=None,
                            given_match_type_vs_protospacer_match_total_ratio=None,
                            given_match_type_vs_protospacer_match_total_ratio_threshold=given_match_type_vs_protospacer_match_total_ratio_threshold,
                            zerocount=protospacer_match_noumi_zerocount,
                            gini=protospacer_match_noumi_gini,
                            std=protospacer_match_noumi_coverage_std,
                            gini_index_threshold=gini_index_threshold,
                            coverage_lower_quantile=protospacer_match_noumi_coverage_lower_quantile,
                            lower_quantile_coverage_percentage=lower_quantile_coverage_percentage,
                            lower_quantile_coverage_threshold=lower_quantile_coverage_threshold,
                            coverage_mean=protospacer_match_noumi_coverage_mean,
                            mean_coverage_threshold=mean_coverage_threshold)
        if contains_surrogate:
            # PROTOSPACER+SURROGATE vs. PROTOSPACER
            protospacer_match_gt_protospacer_surrogate_match_noumi_total = sum(count_series_result.protospacer_match.ambiguous_accepted_counterseries >= count_series_result.protospacer_match_surrogate_match.ambiguous_accepted_counterseries)
            protospacer_surrogate_match_vs_protospacer_match_noumi_total_ratio_threshold = sum(count_series_result.protospacer_match_surrogate_match.ambiguous_accepted_counterseries) / sum(count_series_result.protospacer_match.ambiguous_accepted_counterseries)
            protospacer_surrogate_match_noumi_totalcount = len(count_series_result.protospacer_match_surrogate_match.ambiguous_accepted_counterseries)
            protospacer_surrogate_match_noumi_zerocount = sum(count_series_result.protospacer_match_surrogate_match.ambiguous_accepted_counterseries == 0)
            protospacer_surrogate_match_noumi_gini = calculate_gini_index(count_series_result.protospacer_match_surrogate_match.ambiguous_accepted_counterseries)
            protospacer_surrogate_match_noumi_coverage_mean = np.mean(count_series_result.protospacer_match_surrogate_match.ambiguous_accepted_counterseries)
            protospacer_surrogate_match_noumi_coverage_std = np.std(count_series_result.protospacer_match_surrogate_match.ambiguous_accepted_counterseries)
            protospacer_surrogate_match_noumi_coverage_lower_quantile = count_series_result.protospacer_match_surrogate_match.ambiguous_accepted_counterseries.quantile(q=low_quantile_coverage_percentage)


            count_series_result_quality_control_stats_dict.update({
                "protospacer_surrogate_match_noumi_totalcount": protospacer_surrogate_match_noumi_totalcount,
                "protospacer_match_gt_protospacer_surrogate_match_noumi_total": protospacer_match_gt_protospacer_surrogate_match_noumi_total,
                "protospacer_surrogate_match_vs_protospacer_match_noumi_total_ratio_threshold": protospacer_surrogate_match_vs_protospacer_match_noumi_total_ratio_threshold,
                "protospacer_surrogate_match_noumi_zerocount": protospacer_surrogate_match_noumi_zerocount,
                "protospacer_surrogate_match_noumi_gini": protospacer_surrogate_match_noumi_gini,
                "protospacer_surrogate_match_noumi_coverage_mean": protospacer_surrogate_match_noumi_coverage_mean,
                "protospacer_surrogate_match_noumi_coverage_std": protospacer_surrogate_match_noumi_coverage_std,
                f"protospacer_surrogate_match_noumi_coverage_lower_quantile_{lower_quantile_coverage_percentage}": protospacer_surrogate_match_noumi_coverage_lower_quantile
            })
            if display_visualizations:
                # Log warnings from count statistics
                log_count_series_warnings(match_type_string="protospacer+surrogate-match",
                                count_type_string="No-UMI",
                                protospacer_match_gt_given_match_type_total=protospacer_match_gt_protospacer_surrogate_match_noumi_total,
                                given_match_type_vs_protospacer_match_total_ratio=protospacer_surrogate_match_vs_protospacer_match_noumi_total_ratio_threshold,
                                given_match_type_vs_protospacer_match_total_ratio_threshold=given_match_type_vs_protospacer_match_total_ratio_threshold,
                                zerocount=protospacer_surrogate_match_noumi_zerocount,
                                gini=protospacer_surrogate_match_noumi_gini,
                                std=protospacer_surrogate_match_noumi_coverage_std,
                                gini_index_threshold=gini_index_threshold,
                                coverage_lower_quantile=protospacer_surrogate_match_noumi_coverage_lower_quantile,
                                lower_quantile_coverage_percentage=lower_quantile_coverage_percentage,
                                lower_quantile_coverage_threshold=lower_quantile_coverage_threshold,
                                coverage_mean=protospacer_surrogate_match_noumi_coverage_mean,
                                mean_coverage_threshold=mean_coverage_threshold)

            if contains_barcode:
                # PROTOSPACER+SURROGATE+BARCODE vs. PROTOSPACER
                protospacer_match_gt_protospacer_surrogate_barcode_match_noumi_total = sum(count_series_result.protospacer_match.ambiguous_accepted_counterseries >= count_series_result.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_counterseries)
                protospacer_surrogate_barcode_match_vs_protospacer_match_noumi_total_ratio_threshold = sum(count_series_result.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_counterseries) / sum(count_series_result.protospacer_match.ambiguous_accepted_counterseries)
                protospacer_surrogate_barcode_match_noumi_totalcount = len(count_series_result.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_counterseries)
                protospacer_surrogate_barcode_match_noumi_zerocount = sum(count_series_result.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_counterseries == 0)
                protospacer_surrogate_barcode_match_noumi_gini = calculate_gini_index(count_series_result.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_counterseries)
                protospacer_surrogate_barcode_match_noumi_coverage_mean = np.mean(count_series_result.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_counterseries)
                protospacer_surrogate_barcode_match_noumi_coverage_std = np.std(count_series_result.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_counterseries)
                protospacer_surrogate_barcode_match_noumi_coverage_lower_quantile = count_series_result.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_counterseries.quantile(q=low_quantile_coverage_percentage)

                count_series_result_quality_control_stats_dict.update({
                    "protospacer_surrogate_barcode_match_noumi_totalcount": protospacer_surrogate_barcode_match_noumi_totalcount,
                    "protospacer_match_gt_protospacer_surrogate_barcode_match_noumi_total": protospacer_match_gt_protospacer_surrogate_barcode_match_noumi_total,
                    "protospacer_surrogate_barcode_match_vs_protospacer_match_noumi_total_ratio_threshold": protospacer_surrogate_barcode_match_vs_protospacer_match_noumi_total_ratio_threshold,
                    "protospacer_surrogate_barcode_match_noumi_zerocount": protospacer_surrogate_barcode_match_noumi_zerocount,
                    "protospacer_surrogate_barcode_match_noumi_gini": protospacer_surrogate_barcode_match_noumi_gini,
                    "protospacer_surrogate_barcode_match_noumi_coverage_mean": protospacer_surrogate_barcode_match_noumi_coverage_mean,
                    "protospacer_surrogate_barcode_match_noumi_coverage_std": protospacer_surrogate_barcode_match_noumi_coverage_std,
                    f"protospacer_surrogate_barcode_match_noumi_coverage_lower_quantile_{lower_quantile_coverage_percentage}": protospacer_surrogate_barcode_match_noumi_coverage_lower_quantile
                })

                if display_visualizations:
                    # Log warnings from count statistics
                    log_count_series_warnings(match_type_string="protospacer+surrogate+barcode-match",
                                    count_type_string="No-UMI",
                                    protospacer_match_gt_given_match_type_total=protospacer_match_gt_protospacer_surrogate_barcode_match_noumi_total,
                                    given_match_type_vs_protospacer_match_total_ratio=protospacer_surrogate_barcode_match_vs_protospacer_match_noumi_total_ratio_threshold,
                                    given_match_type_vs_protospacer_match_total_ratio_threshold=given_match_type_vs_protospacer_match_total_ratio_threshold,
                                    zerocount=protospacer_surrogate_match_noumi_zerocount,
                                    gini=protospacer_surrogate_match_noumi_gini,
                                    std=protospacer_surrogate_match_noumi_coverage_std,
                                    gini_index_threshold=gini_index_threshold,
                                    coverage_lower_quantile=protospacer_surrogate_match_noumi_coverage_lower_quantile,
                                    lower_quantile_coverage_percentage=lower_quantile_coverage_percentage,
                                    lower_quantile_coverage_threshold=lower_quantile_coverage_threshold,
                                    coverage_mean=protospacer_surrogate_match_noumi_coverage_mean,
                                    mean_coverage_threshold=mean_coverage_threshold)

        else:
            if contains_barcode:
                # PROTOSPACER+BARCODE vs. PROTOSPACER
                protospacer_match_gt_protospacer_barcode_match_noumi_total = sum(count_series_result.protospacer_match.ambiguous_accepted_counterseries >= count_series_result.protospacer_match_barcode_match.ambiguous_accepted_counterseries)
                protospacer_barcode_match_vs_protospacer_match_noumi_total_ratio_threshold = sum(count_series_result.protospacer_match_barcode_match.ambiguous_accepted_counterseries) / sum(count_series_result.protospacer_match.ambiguous_accepted_counterseries)
                protospacer_barcode_match_noumi_totalcount = len(count_series_result.protospacer_match_barcode_match.ambiguous_accepted_counterseries)
                protospacer_barcode_match_noumi_zerocount = sum(count_series_result.protospacer_match_barcode_match.ambiguous_accepted_counterseries == 0)
                protospacer_barcode_match_noumi_gini = calculate_gini_index(count_series_result.protospacer_match_barcode_match.ambiguous_accepted_counterseries)
                protospacer_barcode_match_noumi_coverage_mean = np.mean(count_series_result.protospacer_match_barcode_match.ambiguous_accepted_counterseries)
                protospacer_barcode_match_noumi_coverage_std = np.std(count_series_result.protospacer_match_barcode_match.ambiguous_accepted_counterseries)
                protospacer_barcode_match_noumi_coverage_lower_quantile = count_series_result.protospacer_match_barcode_match.ambiguous_accepted_counterseries.quantile(q=low_quantile_coverage_percentage)

                count_series_result_quality_control_stats_dict.update({
                    "protospacer_barcode_match_noumi_totalcount", protospacer_barcode_match_noumi_totalcount,
                    "protospacer_match_gt_protospacer_barcode_match_noumi_total": protospacer_match_gt_protospacer_barcode_match_noumi_total,
                    "protospacer_barcode_match_vs_protospacer_match_noumi_total_ratio_threshold": protospacer_barcode_match_vs_protospacer_match_noumi_total_ratio_threshold,
                    "protospacer_barcode_match_noumi_zerocount": protospacer_barcode_match_noumi_zerocount,
                    "protospacer_barcode_match_noumi_gini": protospacer_barcode_match_noumi_gini,
                    "protospacer_barcode_match_noumi_coverage_mean": protospacer_barcode_match_noumi_coverage_mean,
                    "protospacer_barcode_match_noumi_coverage_std": protospacer_barcode_match_noumi_coverage_std,
                    f"protospacer_barcode_match_noumi_coverage_lower_quantile_{lower_quantile_coverage_percentage}": protospacer_barcode_match_noumi_coverage_lower_quantile
                })

                if display_visualizations:
                    # Log warnings from count statistics
                    log_count_series_warnings(match_type_string="protospacer+barcode-match",
                                    count_type_string="No-UMI",
                                    protospacer_match_gt_given_match_type_total=protospacer_match_gt_protospacer_barcode_match_noumi_total,
                                    given_match_type_vs_protospacer_match_total_ratio=protospacer_barcode_match_vs_protospacer_match_noumi_total_ratio_threshold,
                                    given_match_type_vs_protospacer_match_total_ratio_threshold=given_match_type_vs_protospacer_match_total_ratio_threshold,
                                    zerocount=protospacer_barcode_match_noumi_zerocount,
                                    gini=protospacer_barcode_match_noumi_gini,
                                    std=protospacer_barcode_match_noumi_coverage_std,
                                    gini_index_threshold=gini_index_threshold,
                                    coverage_lower_quantile=protospacer_barcode_match_noumi_coverage_lower_quantile,
                                    lower_quantile_coverage_percentage=lower_quantile_coverage_percentage,
                                    lower_quantile_coverage_threshold=lower_quantile_coverage_threshold,
                                    coverage_mean=protospacer_barcode_match_noumi_coverage_mean,
                                    mean_coverage_threshold=mean_coverage_threshold)
    return count_series_result_quality_control_stats_dict
    

def plot_match_type_comparisons(count_series_result: MatchSetWhitelistReporterCounterSeriesResults, contains_umi: bool, contains_surrogate: bool, contains_barcode: bool, collapsed:bool=False, alpha: float=0.3):
    if contains_umi:
        if contains_surrogate:
            if collapsed:
                # PROTOSPACER+SURROGATE vs. PROTOSPACER - COLLAPSED
                plt.scatter(count_series_result.protospacer_match_surrogate_match.ambiguous_accepted_umi_collapsed_counterseries,
                        count_series_result.protospacer_match.ambiguous_accepted_umi_collapsed_counterseries,
                    alpha=alpha)
                plt.xlabel("Protospacer+Surrogate Match")
                plt.ylabel("Protospacer Match")
                plt.title("Comparison of COLLAPSED Protospacer+Surrogate Match vs. Protospacer-only match")
                plt.show()
            else:
                # PROTOSPACER+SURROGATE vs. PROTOSPACER - NONCOLLAPSED
                plt.scatter(count_series_result.protospacer_match_surrogate_match.ambiguous_accepted_umi_noncollapsed_counterseries,
                        count_series_result.protospacer_match.ambiguous_accepted_umi_noncollapsed_counterseries,
                    alpha=alpha)
                plt.xlabel("Protospacer+Surrogate Match")
                plt.ylabel("Protospacer Match")
                plt.title("Comparison of NON-COLLAPSED Protospacer+Surrogate Match vs. Protospacer-only match")
                plt.show()

            if contains_barcode:
                if collapsed:
                    # PROTOSPACER+SURROGATE+BARCODE vs. PROTOSPACER - COLLAPSED
                    plt.scatter(count_series_result.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_umi_collapsed_counterseries,
                        count_series_result.protospacer_match.ambiguous_accepted_umi_collapsed_counterseries,
                    alpha=alpha)
                    plt.xlabel("Protospacer+Surrogate+Barcode Match")
                    plt.ylabel("Protospacer Match")
                    plt.title("Comparison of COLLAPSED Protospacer+Surrogate+Barcode Match vs. Protospacer-only match")
                    plt.show()
                else:
                    # PROTOSPACER+SURROGATE+BARCODE vs. PROTOSPACER - NONCOLLAPSED
                    plt.scatter(count_series_result.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_umi_noncollapsed_counterseries,
                        count_series_result.protospacer_match.ambiguous_accepted_umi_noncollapsed_counterseries,
                    alpha=alpha)
                    plt.xlabel("Protospacer+Surrogate+Barcode Match")
                    plt.ylabel("Protospacer Match")
                    plt.title("Comparison of NON-COLLAPSED Protospacer+Surrogate+Barcode Match vs. Protospacer-only match")
                    plt.show()
        else:
            if contains_barcode:
                if collapsed:
                    # PROTOSPACER+BARCODE vs. PROTOSPACER - COLLAPSED
                    plt.scatter(count_series_result.protospacer_match_barcode_match.ambiguous_accepted_umi_collapsed_counterseries,
                        count_series_result.protospacer_match.ambiguous_accepted_umi_collapsed_counterseries,
                    alpha=alpha)
                    plt.xlabel("Protospacer+Barcode Match")
                    plt.ylabel("Protospacer Match")
                    plt.title("Comparison of COLLAPSED Protospacer+Barcode Match vs. Protospacer-only match")
                    plt.show()
                else:
                    # PROTOSPACER+BARCODE vs. PROTOSPACER - NONCOLLAPSED
                    plt.scatter(count_series_result.protospacer_match_barcode_match.ambiguous_accepted_umi_noncollapsed_counterseries,
                        count_series_result.protospacer_match.ambiguous_accepted_umi_noncollapsed_counterseries,
                    alpha=alpha)
                    plt.xlabel("Protospacer+Barcode Match")
                    plt.ylabel("Protospacer Match")
                    plt.title("Comparison of NON-COLLAPSED Protospacer+Barcode Match vs. Protospacer-only match")
                    plt.show()
    else:

        if contains_surrogate:
            # PROTOSPACER+SURROGATE vs. PROTOSPACER
            plt.scatter(count_series_result.protospacer_match_surrogate_match.ambiguous_accepted_counterseries,
                    count_series_result.protospacer_match.ambiguous_accepted_counterseries,
                alpha=alpha)
            plt.xlabel("Protospacer+Surrogate Match")
            plt.ylabel("Protospacer Match")
            plt.title("Comparison of Protospacer+Surrogate Match vs. Protospacer-only match")
            plt.show()
            
            if contains_barcode:

                # PROTOSPACER+SURROGATE+BARCODE vs. PROTOSPACER
                plt.scatter(count_series_result.protospacer_match_surrogate_match_barcode_match.ambiguous_accepted_counterseries,
                    count_series_result.protospacer_match.ambiguous_accepted_counterseries,
                alpha=alpha)
                plt.xlabel("Protospacer+Surrogate+Barcode Match")
                plt.ylabel("Protospacer Match")
                plt.title("Comparison of Protospacer+Surrogate+Barcode Match vs. Protospacer-only match")
                plt.show()
        else:
            if contains_barcode:

                # PROTOSPACER+BARCODE vs. PROTOSPACER
                plt.scatter(count_series_result.protospacer_match_barcode_match.ambiguous_accepted_counterseries,
                    count_series_result.protospacer_match.ambiguous_accepted_counterseries,
                alpha=alpha)
                plt.xlabel("Protospacer+Barcode Match")
                plt.ylabel("Protospacer Match")
                plt.title("Comparison of Protospacer+Barcode Match vs. Protospacer-only match")
                plt.show()
