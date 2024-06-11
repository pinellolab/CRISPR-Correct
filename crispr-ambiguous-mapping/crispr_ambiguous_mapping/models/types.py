from typing import Union, List, Mapping, Tuple, Optional, Any, DefaultDict
from typing import Counter as CounterType

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