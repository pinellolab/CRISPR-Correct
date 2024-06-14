from dataclasses import dataclass
from typing import Counter as CounterType
from typing import Tuple, Optional, DefaultDict, Dict
from .mapping_models import GeneralAlleleCountSeriesDict
import pandas as pd

@dataclass
class MatchSetWhitelistReporterObservedSequenceCounterSeriesResults:
    
    # Storing as a dictionary 
    ambiguous_ignored_umi_noncollapsed_alleleseries_dict : Optional[GeneralAlleleCountSeriesDict] = None
    ambiguous_ignored_umi_collapsed_alleleseries_dict : Optional[GeneralAlleleCountSeriesDict] = None
    ambiguous_ignored_alleleseries_dict : Optional[GeneralAlleleCountSeriesDict] = None

    ambiguous_accepted_umi_noncollapsed_alleleseries_dict : Optional[GeneralAlleleCountSeriesDict] = None
    ambiguous_accepted_umi_collapsed_alleleseries_dict : Optional[GeneralAlleleCountSeriesDict] = None
    ambiguous_accepted_alleleseries_dict : Optional[GeneralAlleleCountSeriesDict] = None

    ambiguous_spread_umi_noncollapsed_alleleseries_dict : Optional[GeneralAlleleCountSeriesDict] = None
    ambiguous_spread_umi_collapsed_alleleseries_dict : Optional[GeneralAlleleCountSeriesDict] = None
    ambiguous_spread_alleleseries_dict : Optional[GeneralAlleleCountSeriesDict] = None
        
    # Storing as a dataframe
    ambiguous_ignored_umi_noncollapsed_allele_df : Optional[pd.DataFrame] = None
    ambiguous_ignored_umi_collapsed_allele_df : Optional[pd.DataFrame] = None
    ambiguous_ignored_allele_df : Optional[pd.DataFrame] = None

    ambiguous_accepted_umi_noncollapsed_allele_df : Optional[pd.DataFrame] = None
    ambiguous_accepted_umi_collapsed_allele_df : Optional[pd.DataFrame] = None
    ambiguous_accepted_allele_df : Optional[pd.DataFrame] = None

    ambiguous_spread_umi_noncollapsed_allele_df : Optional[pd.DataFrame] = None
    ambiguous_spread_umi_collapsed_allele_df : Optional[pd.DataFrame] = None
    ambiguous_spread_allele_df : Optional[pd.DataFrame] = None

    # This ensures that any empty series are kept at None
    def __setattr__(self, name, value):
        if isinstance(value, dict):
            if len(value) == 0: # If no items in dictionary, just return None
                super().__setattr__(name, None)
            else:
                super().__setattr__(name, value)
        elif isinstance(value, pd.DataFrame):
            if value.empty: # If no items in dataframe, just return None
                super().__setattr__(name, None)
            else:
                super().__setattr__(name, value)
        else:
            # Set the attribute
            super().__setattr__(name, value)

@dataclass
class ObservedSequenceMutationProfile:
    linked_mutations_whitelist_reporter_dict: Dict[Tuple[str, Optional[str], Optional[str]], pd.DataFrame] = None
    all_observed_protospacer_unlinked_mutations_df: pd.DataFrame = None
    all_observed_surrogate_unlinked_mutations_df: Optional[pd.DataFrame] = None
    all_observed_barcode_unlinked_mutations_df: Optional[pd.DataFrame] = None
    
@dataclass
class MatchSetWhitelistReporterObservedSequenceMutationProfiles:
    ambiguous_ignored_umi_noncollapsed_mutations : Optional[ObservedSequenceMutationProfile] = None
    ambiguous_ignored_umi_collapsed_mutations : Optional[ObservedSequenceMutationProfile] = None
    ambiguous_ignored_mutations : Optional[ObservedSequenceMutationProfile] = None

    ambiguous_accepted_umi_noncollapsed_mutations : Optional[ObservedSequenceMutationProfile] = None
    ambiguous_accepted_umi_collapsed_mutations : Optional[ObservedSequenceMutationProfile] = None
    ambiguous_accepted_mutations : Optional[ObservedSequenceMutationProfile] = None

    ambiguous_spread_umi_noncollapsed_mutations : Optional[ObservedSequenceMutationProfile] = None
    ambiguous_spread_umi_collapsed_mutations : Optional[ObservedSequenceMutationProfile] = None
    ambiguous_spread_mutations : Optional[ObservedSequenceMutationProfile] = None

@dataclass
class LinkedMutationCounters:
    protospacer_total_mutation_counter: CounterType
    surrogate_total_mutation_counter: Optional[CounterType]
    barcode_total_mutation_counter: Optional[CounterType]
    