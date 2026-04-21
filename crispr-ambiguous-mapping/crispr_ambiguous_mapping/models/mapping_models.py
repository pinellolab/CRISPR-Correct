from dataclasses import dataclass
from enum import Enum
from typing import Union, List, Mapping, Tuple, Optional, Any, DefaultDict, Dict
from typing import Counter as CounterType

from collections import Counter, defaultdict
import pandas as pd

from .error_models import GuideCountError
from .quality_control_models import QualityControlResult


class MatchTier(str, Enum):
    """Type-safe names for the six mapping tiers exposed by `CompleteInferenceMatchResult`.

    §4.2 / §7.5: previously users passed magic strings like
    ``"protospacer_match_surrogate_match_barcode_match"`` to post-processing
    functions. The enum subclasses ``str`` so existing magic-string callers
    keep working (``MatchTier.PM_SM_BM == "protospacer_match_surrogate_match_barcode_match"``
    is True); new code can use the typed members for IDE autocompletion and
    static checking.
    """
    PM = "protospacer_match"
    PM_SM = "protospacer_match_surrogate_match"
    PM_BM = "protospacer_match_barcode_match"
    PM_SM_BM = "protospacer_match_surrogate_match_barcode_match"
    PM_MISMATCH_SM = "protospacer_mismatch_surrogate_match"
    PM_MISMATCH_SM_BM = "protospacer_mismatch_surrogate_match_barcode_match"

from typing import Union, List, Mapping, Tuple, Optional, Any, DefaultDict, Dict
from typing import Counter as CounterType
import pandas as pd

# Provides a mapping between the whitelist sequence and the list of observed sequences. Each observed sequence has the corresponding count Union[int, Dict [ str, int ]] (either int if no UMI, or Dict[str, int] for UMI-collapsed and non-collpased count)
WhitelistReporterObservedSequenceMapping = DefaultDict[  Tuple[str, Optional[str], Optional[str]],    List[   Tuple[ Tuple[str, Optional[str], Optional[str]], Union[int, Dict [ str, int ]], int, List[Tuple[str, Optional[str], Optional[str]]] ]    ]]

### Sequence Count Result Objects
## No cell barcode
# For non UMI counts
ProtospacerCounter = CounterType[str]
ProtospacerSurrogateCounter = CounterType[Tuple[str, str]]
ProtospacerBarcodeCounter = CounterType[Tuple[str, str]]
ProtospacerSurrogateBarcodeCounter = CounterType[Tuple[str, str, str]]

# For UMI counts
ProtospacerDictUMICounter = DefaultDict[str, CounterType[str]]
ProtospacerSurrogateDictUMICounter = DefaultDict[Tuple[str, str], CounterType[str]]
ProtospacerBarcodeDictUMICounter = DefaultDict[Tuple[str, str], CounterType[str]]
ProtospacerSurrogateBarcodeDictUMICounter = DefaultDict[Tuple[str, str, str], CounterType[str]]


## With sample barcode
# For non UMI counts
SampleProtospacerCounter = DefaultDict[str, CounterType[str]]
SampleProtospacerSurrogateCounter = DefaultDict[Tuple[str, str], CounterType[str]]
SampleProtospacerBarcodeCounter = DefaultDict[Tuple[str, str], CounterType[str]]
SampleProtospacerSurrogateBarcodeCounter = DefaultDict[Tuple[str, str, str], CounterType[str]]

# For UMI counts
SampleProtospacerDictUMICounter = DefaultDict[str, DefaultDict[str, CounterType[str]]]
SampleProtospacerSurrogateDictUMICounter = DefaultDict[Tuple[str, str], DefaultDict[str, CounterType[str]]]
SampleProtospacerBarcodeDictUMICounter = DefaultDict[Tuple[str, str], DefaultDict[str, CounterType[str]]]
SampleProtospacerSurrogateBarcodeDictUMICounter = DefaultDict[Tuple[str, str, str], DefaultDict[str, CounterType[str]]]

GeneralGuideCountType = Union[ProtospacerCounter, 
                       ProtospacerSurrogateCounter, 
                       ProtospacerBarcodeCounter, 
                       ProtospacerSurrogateBarcodeCounter,
                       ProtospacerDictUMICounter,
                       ProtospacerSurrogateDictUMICounter,
                       ProtospacerBarcodeDictUMICounter, 
                       ProtospacerSurrogateBarcodeDictUMICounter,
                       SampleProtospacerCounter, 
                       SampleProtospacerSurrogateCounter, 
                       SampleProtospacerBarcodeCounter, 
                       SampleProtospacerSurrogateBarcodeCounter,
                       SampleProtospacerDictUMICounter,
                       SampleProtospacerSurrogateDictUMICounter,
                       SampleProtospacerBarcodeDictUMICounter, 
                       SampleProtospacerSurrogateBarcodeDictUMICounter]

# Mapping Count Dict Object

ProtospacerSurrogateBarcodeMatchCountDict = DefaultDict[Tuple[str, str, str], Union[int, float]]
ProtospacerSurrogateMatchCountDict = DefaultDict[Tuple[str, str], Union[int, float]]
ProtospacerBarcodeMatchCountDict = DefaultDict[Tuple[str, str], Union[int, float]]
ProtospacerMatchCountDict = DefaultDict[str, Union[int, float]]
GeneralMatchCountDict = Union[ProtospacerSurrogateBarcodeMatchCountDict, 
      ProtospacerSurrogateMatchCountDict, 
      ProtospacerBarcodeMatchCountDict, 
      ProtospacerMatchCountDict]

ProtospacerSurrogateBarcodeMismatchCountDict = DefaultDict[Tuple[Tuple[str, str, str], Tuple[str, str, str]], Union[int, float]]
ProtospacerSurrogateMismatchCountDict = DefaultDict[Tuple[Tuple[str, str], Tuple[str, str]], Union[int, float]]
ProtospacerBarcodeMismatchCountDict = DefaultDict[Tuple[Tuple[str, str], Tuple[str, str]], Union[int, float]]
ProtospacerMismatchCountDict = DefaultDict[Tuple[str, str], Union[int, float]]

GeneralMismatchCountDict = Union[ProtospacerSurrogateBarcodeMismatchCountDict, 
      ProtospacerSurrogateMismatchCountDict, 
      ProtospacerBarcodeMismatchCountDict, 
      ProtospacerMismatchCountDict]

# Allele nested dict Object (Key of first dict is inferred, key of second dict is observed, value of second dict is count)
ProtospacerSurrogateBarcodeAlleleDict = DefaultDict[Tuple[str, str, str], DefaultDict[Tuple[str, str, str], Union[int, float]]]
ProtospacerSurrogateAlleleDict = DefaultDict[Tuple[str, str], DefaultDict[Tuple[str, str], Union[int, float]]]
ProtospacerBarcodeAlleleDict = DefaultDict[Tuple[str, str], DefaultDict[Tuple[str, str], Union[int, float]]]
ProtospacerAlleleDict = DefaultDict[str, DefaultDict[str, Union[int, float]]]

GeneralAlleleDict = Union[ProtospacerSurrogateBarcodeAlleleDict,
      ProtospacerSurrogateAlleleDict,
      ProtospacerBarcodeAlleleDict,
      ProtospacerAlleleDict]




# Allele Count Series Dict
ProtospacerSurrogateBarcodeAlleleCountSeriesDict = DefaultDict[Tuple[str, str, str], pd.Series]
ProtospacerSurrogateAlleleCountSeriesDict = DefaultDict[Tuple[str, str], pd.Series]
ProtospacerBarcodeAlleleCountSeriesDict = DefaultDict[Tuple[str, str], pd.Series]
ProtospacerAlleleCountSeriesDict = DefaultDict[str, pd.Series]

GeneralAlleleCountSeriesDict = Union[ProtospacerSurrogateBarcodeAlleleCountSeriesDict,
                                     ProtospacerSurrogateAlleleCountSeriesDict,
                                     ProtospacerBarcodeAlleleCountSeriesDict,
                                     ProtospacerAlleleCountSeriesDict]







@dataclass
class SingleInferenceMatchResultValue:
    pass

@dataclass
class MatchSetSingleInferenceMatchResultValue(SingleInferenceMatchResultValue):
    matches: Optional[pd.DataFrame] = None

@dataclass
class SurrogateProtospacerMismatchSingleInferenceMatchResultValue(SingleInferenceMatchResultValue):
    mismatched: bool
    surrogate_matches: pd.DataFrame
    protospacer_matches: pd.DataFrame 
    protospacer_surrogate_matches: pd.DataFrame

@dataclass
class SingleInferenceMatchResult:
    error: Optional[GuideCountError] = None
    value: Optional[SingleInferenceMatchResultValue] = None

@dataclass
class MatchSetSingleInferenceMatchResult(SingleInferenceMatchResult):
    value: Optional[MatchSetSingleInferenceMatchResultValue] = None

@dataclass
class SurrogateProtospacerMismatchSingleInferenceMatchResult(SingleInferenceMatchResult):
    value: Optional[SurrogateProtospacerMismatchSingleInferenceMatchResultValue] = None

@dataclass
class CompleteInferenceMatchResult:
    protospacer_match: Optional[MatchSetSingleInferenceMatchResult] = None
    protospacer_match_surrogate_match_barcode_match: Optional[MatchSetSingleInferenceMatchResult] = None
    protospacer_match_surrogate_match: Optional[MatchSetSingleInferenceMatchResult] = None
    protospacer_match_barcode_match: Optional[MatchSetSingleInferenceMatchResult] =None
    protospacer_mismatch_surrogate_match_barcode_match: Optional[SurrogateProtospacerMismatchSingleInferenceMatchResult] = None
    protospacer_mismatch_surrogate_match: Optional[SurrogateProtospacerMismatchSingleInferenceMatchResult] = None

@dataclass
class InferenceResult:
    observed_value: Union[int, CounterType[Optional[str]]] # Count
    inferred_value: CompleteInferenceMatchResult


# Inference Result Object
ProtospacerSurrogateBarcodeMappingInferenceDict = DefaultDict[Tuple[str,str,str], Dict[Tuple[str,str,str], InferenceResult]]
ProtospacerSurrogateMappingInferenceDict = DefaultDict[Tuple[str,str], Dict[Tuple[str,str], InferenceResult]]
ProtospacerBarcodeMappingInferenceDict = DefaultDict[Tuple[str,str], Dict[Tuple[str,str], InferenceResult]]
ProtospacerMappingInferenceDict = DefaultDict[str, Dict[str, InferenceResult]]

GeneralMappingInferenceDict = Union[ProtospacerSurrogateBarcodeMappingInferenceDict,
      ProtospacerSurrogateMappingInferenceDict,
      ProtospacerBarcodeMappingInferenceDict,
      ProtospacerMappingInferenceDict]


@dataclass
class MatchSetWhitelistReporterCounterSeriesResults:
    ambiguous_ignored_umi_noncollapsed_counterseries : Optional[pd.Series] = None
    ambiguous_ignored_umi_collapsed_counterseries : Optional[pd.Series] = None
    ambiguous_ignored_counterseries : Optional[pd.Series] = None

    ambiguous_accepted_umi_noncollapsed_counterseries : Optional[pd.Series] = None
    ambiguous_accepted_umi_collapsed_counterseries : Optional[pd.Series] = None
    ambiguous_accepted_counterseries : Optional[pd.Series] = None

    ambiguous_spread_umi_noncollapsed_counterseries : Optional[pd.Series] = None
    ambiguous_spread_umi_collapsed_counterseries : Optional[pd.Series] = None
    ambiguous_spread_counterseries : Optional[pd.Series] = None

    # This ensures that any empty series are kept at None
    def __setattr__(self, name, value):
        if isinstance(value, pd.Series):
            if value.sum() == 0: # If sum of pandas series is 0, just return None
                super().__setattr__(name, None)
            else:
                super().__setattr__(name, value)
        else:
            # Set the attribute
            super().__setattr__(name, value)

@dataclass
class SurrogateProtospacerMismatchSetWhitelistReporterCounterSeriesResults:
    
    # MATCH counters
    ambiguous_ignored_umi_noncollapsed_match_counterseries : Optional[pd.Series]  = None
    ambiguous_ignored_umi_collapsed_match_counterseries : Optional[pd.Series]  = None
    ambiguous_ignored_match_counterseries : Optional[pd.Series] = None

    ambiguous_accepted_umi_noncollapsed_match_counterseries : Optional[pd.Series] = None
    ambiguous_accepted_umi_collapsed_match_counterseries : Optional[pd.Series] = None
    ambiguous_accepted_match_counterseries : Optional[pd.Series] = None

    ambiguous_spread_umi_noncollapsed_match_counterseries : Optional[pd.Series] = None
    ambiguous_spread_umi_collapsed_match_counterseries : Optional[pd.Series] = None
    ambiguous_spread_match_counterseries : Optional[pd.Series] = None

    # MISMATCH counters (keys are PAIRS of indices, representing the protospacer and surrogate match separately), [str, Optional[str], Optional[str], str, Optional[str], Optional[str]], first three are Protospacer_Match, last three is Surrogate_Match
    ambiguous_ignored_umi_noncollapsed_mismatch_counterseries : Optional[pd.Series] = None
    ambiguous_ignored_umi_collapsed_mismatch_counterseries : Optional[pd.Series] = None
    ambiguous_ignored_mismatch_counterseries : Optional[pd.Series] = None

    ambiguous_accepted_umi_noncollapsed_mismatch_counterseries : Optional[pd.Series] = None
    ambiguous_accepted_umi_collapsed_mismatch_counterseries : Optional[pd.Series] = None
    ambiguous_accepted_mismatch_counterseries : Optional[pd.Series] = None

    ambiguous_spread_umi_noncollapsed_mismatch_counterseries : Optional[pd.Series] = None
    ambiguous_spread_umi_collapsed_mismatch_counterseries : Optional[pd.Series] = None
    ambiguous_spread_mismatch_counterseries : Optional[pd.Series] = None

    # This ensures that any empty series are kept at None
    def __setattr__(self, name, value):
        if isinstance(value, pd.Series):
            if value.sum() == 0: # If sum of pandas series is 0, just return None
                super().__setattr__(name, None)
            else:
                super().__setattr__(name, value)
        else:
            # Set the attribute
            super().__setattr__(name, value)

@dataclass
class AllMatchSetWhitelistReporterCounterSeriesResults:
    protospacer_match: Optional[MatchSetWhitelistReporterCounterSeriesResults] = None
    protospacer_match_surrogate_match_barcode_match: Optional[MatchSetWhitelistReporterCounterSeriesResults] = None
    protospacer_match_surrogate_match: Optional[MatchSetWhitelistReporterCounterSeriesResults] = None
    protospacer_match_barcode_match: Optional[MatchSetWhitelistReporterCounterSeriesResults] = None
    protospacer_mismatch_surrogate_match_barcode_match: Optional[SurrogateProtospacerMismatchSetWhitelistReporterCounterSeriesResults] = None
    protospacer_mismatch_surrogate_match: Optional[SurrogateProtospacerMismatchSetWhitelistReporterCounterSeriesResults] = None

@dataclass
class CountInput:
    whitelist_guide_reporter_df: pd.DataFrame
    contains_surrogate:bool
    contains_guide_barcode:bool
    contains_guide_umi:bool
    contains_sample_barcode:bool
    protospacer_hamming_threshold_strict: Optional[int]
    surrogate_hamming_threshold_strict: Optional[int]
    guide_barcode_hamming_threshold_strict: Optional[int]

@dataclass
class WhitelistReporterCountsResult:
    all_match_set_whitelist_reporter_counter_series_results: AllMatchSetWhitelistReporterCounterSeriesResults
    # `observed_guide_reporter_umi_counts_inferred` is None by default on new
    # runs (pass retain_inference_results=True to the mapping call to keep it).
    # Only allele / mutation post-processing reads this dict; the slim default
    # drops ~77% of the result pickle size.
    observed_guide_reporter_umi_counts_inferred: Optional[Union[GeneralMappingInferenceDict, DefaultDict[str, GeneralMappingInferenceDict]]]
    quality_control_result: QualityControlResult
    count_input: CountInput

