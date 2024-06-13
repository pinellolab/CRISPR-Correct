from typing import Union, List, Mapping, Tuple, Optional, Any, DefaultDict, Dict
from typing import Counter as CounterType
from .mapping_models import InferenceResult
import pandas as pd


# Sequence Count Result Objects
ProtospacerCounter = CounterType[str]
ProtospacerSurrogateCounter = CounterType[Tuple[str, str]]
ProtospacerBarcodeCounter = CounterType[Tuple[str, str]]
ProtospacerSurrogateBarcodeCounter = CounterType[Tuple[str, str, str]]

ProtospacerDictUMICounter = DefaultDict[str, CounterType[str]]
ProtospacerSurrogateDictUMICounter = DefaultDict[Tuple[str, str], CounterType[str]]
ProtospacerBarcodeDictUMICounter = DefaultDict[Tuple[str, str], CounterType[str]]
ProtospacerSurrogateBarcodeDictUMICounter = DefaultDict[Tuple[str, str, str], CounterType[str]]

GeneralGuideCountType = Union[ProtospacerCounter, 
                       ProtospacerSurrogateCounter, 
                       ProtospacerBarcodeCounter, 
                       ProtospacerSurrogateBarcodeCounter,
                       ProtospacerDictUMICounter,
                       ProtospacerSurrogateDictUMICounter,
                       ProtospacerBarcodeDictUMICounter, 
                       ProtospacerSurrogateBarcodeDictUMICounter]





# Inference Result Object
ProtospacerSurrogateBarcodeMappingInferenceDict = DefaultDict[Tuple[str,str,str], Dict[InferenceResult]]
ProtospacerSurrogateMappingInferenceDict = DefaultDict[Tuple[str,str], Dict[InferenceResult]]
ProtospacerBarcodeMappingInferenceDict = DefaultDict[Tuple[str,str], Dict[InferenceResult]]
ProtospacerMappingInferenceDict = DefaultDict[str, Dict[InferenceResult]]

GeneralMappingInferenceDict = Union[ProtospacerSurrogateBarcodeMappingInferenceDict,
      ProtospacerSurrogateMappingInferenceDict,
      ProtospacerBarcodeMappingInferenceDict,
      ProtospacerMappingInferenceDict]


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




