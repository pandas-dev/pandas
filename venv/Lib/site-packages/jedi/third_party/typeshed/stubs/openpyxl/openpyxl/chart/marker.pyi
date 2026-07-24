from _typeshed import ConvertibleToFloat, ConvertibleToInt, Unused
from typing import ClassVar, Literal
from typing_extensions import TypeAlias

from openpyxl.chart.picture import PictureOptions
from openpyxl.chart.shapes import GraphicalProperties
from openpyxl.descriptors.base import Alias, Typed, _ConvertibleToBool
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.nested import NestedBool, NestedInteger, NestedMinMax, NestedNoneSet, _NestedNoneSetParam
from openpyxl.descriptors.serialisable import Serialisable

from ..xml._functions_overloads import _HasTagAndGet

_MarkerSymbol: TypeAlias = Literal[
    "circle", "dash", "diamond", "dot", "picture", "plus", "square", "star", "triangle", "x", "auto"
]

class Marker(Serialisable):
    tagname: ClassVar[str]
    symbol: NestedNoneSet[_MarkerSymbol]
    size: NestedMinMax[float, Literal[True]]
    spPr: Typed[GraphicalProperties, Literal[True]]
    graphicalProperties: Alias
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        symbol: _NestedNoneSetParam[_MarkerSymbol] = None,
        size: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None = None,
        spPr: GraphicalProperties | None = None,
        extLst: Unused = None,
    ) -> None: ...

class DataPoint(Serialisable):
    tagname: ClassVar[str]
    idx: NestedInteger[Literal[False]]
    invertIfNegative: NestedBool[Literal[True]]
    marker: Typed[Marker, Literal[True]]
    bubble3D: NestedBool[Literal[True]]
    explosion: NestedInteger[Literal[True]]
    spPr: Typed[GraphicalProperties, Literal[True]]
    graphicalProperties: Alias
    pictureOptions: Typed[PictureOptions, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        idx: _HasTagAndGet[ConvertibleToInt] | ConvertibleToInt,
        invertIfNegative: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        marker: Marker | None = None,
        bubble3D: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        explosion: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        spPr: GraphicalProperties | None = None,
        pictureOptions: PictureOptions | None = None,
        extLst: Unused = None,
    ) -> None: ...
