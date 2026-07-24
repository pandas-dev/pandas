from _typeshed import ConvertibleToFloat, ConvertibleToInt, Incomplete
from typing import ClassVar, Literal
from typing_extensions import TypeAlias

from openpyxl.descriptors.base import Alias, Bool, Integer, MinMax, NoneSet, Set, Typed, _ConvertibleToBool
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.nested import NestedNoneSet, NestedValue, _NestedNoneSetParam
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.drawing.colors import ColorChoice, HSLColor, RGBPercent as _RGBPercent, SchemeColor, SystemColor, _PresetColors
from openpyxl.drawing.effect import (
    AlphaBiLevelEffect,
    AlphaCeilingEffect,
    AlphaFloorEffect,
    AlphaInverseEffect,
    AlphaModulateEffect,
    AlphaModulateFixedEffect,
    AlphaReplaceEffect,
    BiLevelEffect,
    BlurEffect,
    ColorChangeEffect,
    ColorReplaceEffect,
    DuotoneEffect,
    FillOverlayEffect,
    GrayscaleEffect,
    HSLEffect,
    LuminanceEffect,
    TintEffect,
)

from ..xml._functions_overloads import _HasTagAndGet

_PatternFillPropertiesPrst: TypeAlias = Literal[
    "pct5",
    "pct10",
    "pct20",
    "pct25",
    "pct30",
    "pct40",
    "pct50",
    "pct60",
    "pct70",
    "pct75",
    "pct80",
    "pct90",
    "horz",
    "vert",
    "ltHorz",
    "ltVert",
    "dkHorz",
    "dkVert",
    "narHorz",
    "narVert",
    "dashHorz",
    "dashVert",
    "cross",
    "dnDiag",
    "upDiag",
    "ltDnDiag",
    "ltUpDiag",
    "dkDnDiag",
    "dkUpDiag",
    "wdDnDiag",
    "wdUpDiag",
    "dashDnDiag",
    "dashUpDiag",
    "diagCross",
    "smCheck",
    "lgCheck",
    "smGrid",
    "lgGrid",
    "dotGrid",
    "smConfetti",
    "lgConfetti",
    "horzBrick",
    "diagBrick",
    "solidDmnd",
    "openDmnd",
    "dotDmnd",
    "plaid",
    "sphere",
    "weave",
    "divot",
    "shingle",
    "wave",
    "trellis",
    "zigZag",
]
_PropertiesFlip: TypeAlias = Literal["x", "y", "xy"]
_TileInfoPropertiesAlgn: TypeAlias = Literal["tl", "t", "tr", "l", "ctr", "r", "bl", "b", "br"]
_BlipCstate: TypeAlias = Literal["email", "screen", "print", "hqprint"]
_PathShadePropertiesPath: TypeAlias = Literal["shape", "circle", "rect"]

class PatternFillProperties(Serialisable):
    tagname: ClassVar[str]
    namespace: ClassVar[str]
    prst: NoneSet[_PatternFillPropertiesPrst]
    preset: Alias
    fgClr: Typed[ColorChoice, Literal[True]]
    foreground: Alias
    bgClr: Typed[ColorChoice, Literal[True]]
    background: Alias
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        prst: _PatternFillPropertiesPrst | Literal["none"] | None = None,
        fgClr: ColorChoice | None = None,
        bgClr: ColorChoice | None = None,
    ) -> None: ...

class RelativeRect(Serialisable):
    tagname: ClassVar[str]
    namespace: ClassVar[str]
    l: Incomplete
    left: Alias
    t: Incomplete
    top: Alias
    r: Incomplete
    right: Alias
    b: Incomplete
    bottom: Alias
    def __init__(self, l=None, t=None, r=None, b=None) -> None: ...

class StretchInfoProperties(Serialisable):
    tagname: ClassVar[str]
    namespace: ClassVar[str]
    fillRect: Typed[RelativeRect, Literal[True]]
    def __init__(self, fillRect: RelativeRect = ...) -> None: ...

class GradientStop(Serialisable):
    tagname: ClassVar[str]
    namespace: ClassVar[str]
    pos: MinMax[float, Literal[True]]
    scrgbClr: Typed[_RGBPercent, Literal[True]]
    RGBPercent: Alias
    srgbClr: NestedValue[_RGBPercent, Literal[True]]
    RGB: Alias
    hslClr: Typed[HSLColor, Literal[True]]
    sysClr: Typed[SystemColor, Literal[True]]
    schemeClr: Typed[SchemeColor, Literal[True]]
    prstClr: NestedNoneSet[_PresetColors]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        pos: ConvertibleToFloat | None = None,
        scrgbClr: _RGBPercent | None = None,
        srgbClr: _HasTagAndGet[_RGBPercent | None] | _RGBPercent | None = None,
        hslClr: HSLColor | None = None,
        sysClr: SystemColor | None = None,
        schemeClr: SchemeColor | None = None,
        prstClr: _NestedNoneSetParam[_PresetColors] = None,
    ) -> None: ...

class LinearShadeProperties(Serialisable):
    tagname: ClassVar[str]
    namespace: ClassVar[str]
    ang: Integer[Literal[False]]
    scaled: Bool[Literal[True]]
    def __init__(self, ang: ConvertibleToInt, scaled: _ConvertibleToBool | None = None) -> None: ...

class PathShadeProperties(Serialisable):
    tagname: ClassVar[str]
    namespace: ClassVar[str]
    path: Set[_PathShadePropertiesPath]
    fillToRect: Typed[RelativeRect, Literal[True]]
    def __init__(self, path: _PathShadePropertiesPath, fillToRect: RelativeRect | None = None) -> None: ...

class GradientFillProperties(Serialisable):
    tagname: ClassVar[str]
    namespace: ClassVar[str]
    flip: NoneSet[_PropertiesFlip]
    rotWithShape: Bool[Literal[True]]
    gsLst: Incomplete
    stop_list: Alias
    lin: Typed[LinearShadeProperties, Literal[True]]
    linear: Alias
    path: Typed[PathShadeProperties, Literal[True]]
    tileRect: Typed[RelativeRect, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        flip: _PropertiesFlip | Literal["none"] | None = None,
        rotWithShape: _ConvertibleToBool | None = None,
        gsLst=(),
        lin: LinearShadeProperties | None = None,
        path: PathShadeProperties | None = None,
        tileRect: RelativeRect | None = None,
    ) -> None: ...

class SolidColorFillProperties(Serialisable):
    tagname: ClassVar[str]
    scrgbClr: Typed[_RGBPercent, Literal[True]]
    RGBPercent: Alias
    srgbClr: NestedValue[_RGBPercent, Literal[True]]
    RGB: Alias
    hslClr: Typed[HSLColor, Literal[True]]
    sysClr: Typed[SystemColor, Literal[True]]
    schemeClr: Typed[SchemeColor, Literal[True]]
    prstClr: NestedNoneSet[_PresetColors]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        scrgbClr: _RGBPercent | None = None,
        srgbClr: _HasTagAndGet[_RGBPercent | None] | _RGBPercent | None = None,
        hslClr: HSLColor | None = None,
        sysClr: SystemColor | None = None,
        schemeClr: SchemeColor | None = None,
        prstClr: _NestedNoneSetParam[_PresetColors] = None,
    ) -> None: ...

class Blip(Serialisable):
    tagname: ClassVar[str]
    namespace: ClassVar[str]
    cstate: NoneSet[_BlipCstate]
    embed: Incomplete
    link: Incomplete
    noGrp: Bool[Literal[True]]
    noSelect: Bool[Literal[True]]
    noRot: Bool[Literal[True]]
    noChangeAspect: Bool[Literal[True]]
    noMove: Bool[Literal[True]]
    noResize: Bool[Literal[True]]
    noEditPoints: Bool[Literal[True]]
    noAdjustHandles: Bool[Literal[True]]
    noChangeArrowheads: Bool[Literal[True]]
    noChangeShapeType: Bool[Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    alphaBiLevel: Typed[AlphaBiLevelEffect, Literal[True]]
    alphaCeiling: Typed[AlphaCeilingEffect, Literal[True]]
    alphaFloor: Typed[AlphaFloorEffect, Literal[True]]
    alphaInv: Typed[AlphaInverseEffect, Literal[True]]
    alphaMod: Typed[AlphaModulateEffect, Literal[True]]
    alphaModFix: Typed[AlphaModulateFixedEffect, Literal[True]]
    alphaRepl: Typed[AlphaReplaceEffect, Literal[True]]
    biLevel: Typed[BiLevelEffect, Literal[True]]
    blur: Typed[BlurEffect, Literal[True]]
    clrChange: Typed[ColorChangeEffect, Literal[True]]
    clrRepl: Typed[ColorReplaceEffect, Literal[True]]
    duotone: Typed[DuotoneEffect, Literal[True]]
    fillOverlay: Typed[FillOverlayEffect, Literal[True]]
    grayscl: Typed[GrayscaleEffect, Literal[True]]
    hsl: Typed[HSLEffect, Literal[True]]
    lum: Typed[LuminanceEffect, Literal[True]]
    tint: Typed[TintEffect, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        cstate: _BlipCstate | Literal["none"] | None = None,
        embed=None,
        link=None,
        noGrp: _ConvertibleToBool | None = None,
        noSelect: _ConvertibleToBool | None = None,
        noRot: _ConvertibleToBool | None = None,
        noChangeAspect: _ConvertibleToBool | None = None,
        noMove: _ConvertibleToBool | None = None,
        noResize: _ConvertibleToBool | None = None,
        noEditPoints: _ConvertibleToBool | None = None,
        noAdjustHandles: _ConvertibleToBool | None = None,
        noChangeArrowheads: _ConvertibleToBool | None = None,
        noChangeShapeType: _ConvertibleToBool | None = None,
        extLst: ExtensionList | None = None,
        alphaBiLevel: AlphaBiLevelEffect | None = None,
        alphaCeiling: AlphaCeilingEffect | None = None,
        alphaFloor: AlphaFloorEffect | None = None,
        alphaInv: AlphaInverseEffect | None = None,
        alphaMod: AlphaModulateEffect | None = None,
        alphaModFix: AlphaModulateFixedEffect | None = None,
        alphaRepl: AlphaReplaceEffect | None = None,
        biLevel: BiLevelEffect | None = None,
        blur: BlurEffect | None = None,
        clrChange: ColorChangeEffect | None = None,
        clrRepl: ColorReplaceEffect | None = None,
        duotone: DuotoneEffect | None = None,
        fillOverlay: FillOverlayEffect | None = None,
        grayscl: GrayscaleEffect | None = None,
        hsl: HSLEffect | None = None,
        lum: LuminanceEffect | None = None,
        tint: TintEffect | None = None,
    ) -> None: ...

class TileInfoProperties(Serialisable):
    tx: Integer[Literal[True]]
    ty: Integer[Literal[True]]
    sx: Integer[Literal[True]]
    sy: Integer[Literal[True]]
    flip: NoneSet[_PropertiesFlip]
    algn: Set[_TileInfoPropertiesAlgn]
    def __init__(
        self,
        tx: ConvertibleToInt | None = None,
        ty: ConvertibleToInt | None = None,
        sx: ConvertibleToInt | None = None,
        sy: ConvertibleToInt | None = None,
        flip: _PropertiesFlip | Literal["none"] | None = None,
        *,
        algn: _TileInfoPropertiesAlgn,
    ) -> None: ...

class BlipFillProperties(Serialisable):
    tagname: ClassVar[str]
    dpi: Integer[Literal[True]]
    rotWithShape: Bool[Literal[True]]
    blip: Typed[Blip, Literal[True]]
    srcRect: Typed[RelativeRect, Literal[True]]
    tile: Typed[TileInfoProperties, Literal[True]]
    stretch: Typed[StretchInfoProperties, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        dpi: ConvertibleToInt | None = None,
        rotWithShape: _ConvertibleToBool | None = None,
        blip: Blip | None = None,
        tile: TileInfoProperties | None = None,
        stretch: StretchInfoProperties = ...,
        srcRect: RelativeRect | None = None,
    ) -> None: ...
