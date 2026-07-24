from _typeshed import ConvertibleToFloat, ConvertibleToInt, Incomplete, Unused
from typing import ClassVar, Literal
from typing_extensions import TypeAlias

from openpyxl.descriptors.base import Alias, Integer, MinMax, NoneSet, Typed, _ConvertibleToBool
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.nested import EmptyTag, NestedInteger, NestedNoneSet, _NestedNoneSetParam
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.drawing.colors import ColorChoice, ColorChoiceDescriptor
from openpyxl.drawing.fill import GradientFillProperties, PatternFillProperties

from ..xml._functions_overloads import _HasTagAndGet

_LineEndPropertiesType: TypeAlias = Literal["none", "triangle", "stealth", "diamond", "oval", "arrow"]
_LineEndPropertiesWLen: TypeAlias = Literal["sm", "med", "lg"]
_LinePropertiesCap: TypeAlias = Literal["rnd", "sq", "flat"]
_LinePropertiesCmpd: TypeAlias = Literal["sng", "dbl", "thickThin", "thinThick", "tri"]
_LinePropertiesAlgn: TypeAlias = Literal["ctr", "in"]
_LinePropertiesPrstDash: TypeAlias = Literal[
    "solid", "dot", "dash", "lgDash", "dashDot", "lgDashDot", "lgDashDotDot", "sysDash", "sysDot", "sysDashDot", "sysDashDotDot"
]

class LineEndProperties(Serialisable):
    tagname: ClassVar[str]
    namespace: ClassVar[str]
    type: NoneSet[_LineEndPropertiesType]
    w: NoneSet[_LineEndPropertiesWLen]
    len: NoneSet[_LineEndPropertiesWLen]
    def __init__(
        self,
        type: _LineEndPropertiesType | Literal["none"] | None = None,
        w: _LineEndPropertiesWLen | Literal["none"] | None = None,
        len: _LineEndPropertiesWLen | Literal["none"] | None = None,
    ) -> None: ...

class DashStop(Serialisable):
    tagname: ClassVar[str]
    namespace: ClassVar[str]
    d: Integer[Literal[False]]
    length: Alias
    sp: Integer[Literal[False]]
    space: Alias
    def __init__(self, d: ConvertibleToInt = 0, sp: ConvertibleToInt = 0) -> None: ...

class DashStopList(Serialisable):
    ds: Incomplete
    def __init__(self, ds=None) -> None: ...

class LineProperties(Serialisable):
    tagname: ClassVar[str]
    namespace: ClassVar[str]
    w: MinMax[float, Literal[True]]
    width: Alias
    cap: NoneSet[_LinePropertiesCap]
    cmpd: NoneSet[_LinePropertiesCmpd]
    algn: NoneSet[_LinePropertiesAlgn]
    noFill: EmptyTag[Literal[False]]
    solidFill: ColorChoiceDescriptor
    gradFill: Typed[GradientFillProperties, Literal[True]]
    pattFill: Typed[PatternFillProperties, Literal[True]]
    prstDash: NestedNoneSet[_LinePropertiesPrstDash]
    dashStyle: Alias
    custDash: Typed[DashStop, Literal[True]]
    round: EmptyTag[Literal[False]]
    bevel: EmptyTag[Literal[False]]
    miter: NestedInteger[Literal[True]]
    headEnd: Typed[LineEndProperties, Literal[True]]
    tailEnd: Typed[LineEndProperties, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        w: ConvertibleToFloat | None = None,
        cap: _LinePropertiesCap | Literal["none"] | None = None,
        cmpd: _LinePropertiesCmpd | Literal["none"] | None = None,
        algn: _LinePropertiesAlgn | Literal["none"] | None = None,
        noFill: _HasTagAndGet[_ConvertibleToBool] | _ConvertibleToBool = None,
        solidFill: str | ColorChoice | None = None,
        gradFill: GradientFillProperties | None = None,
        pattFill: PatternFillProperties | None = None,
        prstDash: _NestedNoneSetParam[_LinePropertiesPrstDash] = None,
        custDash: DashStop | None = None,
        round: _HasTagAndGet[_ConvertibleToBool] | _ConvertibleToBool = None,
        bevel: _HasTagAndGet[_ConvertibleToBool] | _ConvertibleToBool = None,
        miter: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        headEnd: LineEndProperties | None = None,
        tailEnd: LineEndProperties | None = None,
        extLst: Unused = None,
    ) -> None: ...
