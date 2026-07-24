from _typeshed import ConvertibleToFloat, ConvertibleToInt, Incomplete
from typing import ClassVar, Final, Literal, overload
from typing_extensions import TypeAlias

from openpyxl.descriptors import Strict, Typed
from openpyxl.descriptors.base import Alias, Integer, MinMax, Set, _ConvertibleToBool
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.nested import EmptyTag, NestedInteger, NestedNoneSet, NestedValue, _NestedNoneSetParam
from openpyxl.descriptors.serialisable import Serialisable

from ..xml._functions_overloads import _HasTagAndGet

_ColorSetType: TypeAlias = Literal[
    "dk1", "lt1", "dk2", "lt2", "accent1", "accent2", "accent3", "accent4", "accent5", "accent6", "hlink", "folHlink"
]
_SystemColorVal: TypeAlias = Literal[
    "scrollBar",
    "background",
    "activeCaption",
    "inactiveCaption",
    "menu",
    "window",
    "windowFrame",
    "menuText",
    "windowText",
    "captionText",
    "activeBorder",
    "inactiveBorder",
    "appWorkspace",
    "highlight",
    "highlightText",
    "btnFace",
    "btnShadow",
    "grayText",
    "btnText",
    "inactiveCaptionText",
    "btnHighlight",
    "3dDkShadow",
    "3dLight",
    "infoText",
    "infoBk",
    "hotLight",
    "gradientActiveCaption",
    "gradientInactiveCaption",
    "menuHighlight",
    "menuBar",
]
_SchemeColors: TypeAlias = Literal[
    "bg1",
    "tx1",
    "bg2",
    "tx2",
    "accent1",
    "accent2",
    "accent3",
    "accent4",
    "accent5",
    "accent6",
    "hlink",
    "folHlink",
    "phClr",
    "dk1",
    "lt1",
    "dk2",
    "lt2",
]
_PresetColors: TypeAlias = Literal[
    "aliceBlue",
    "antiqueWhite",
    "aqua",
    "aquamarine",
    "azure",
    "beige",
    "bisque",
    "black",
    "blanchedAlmond",
    "blue",
    "blueViolet",
    "brown",
    "burlyWood",
    "cadetBlue",
    "chartreuse",
    "chocolate",
    "coral",
    "cornflowerBlue",
    "cornsilk",
    "crimson",
    "cyan",
    "darkBlue",
    "darkCyan",
    "darkGoldenrod",
    "darkGray",
    "darkGrey",
    "darkGreen",
    "darkKhaki",
    "darkMagenta",
    "darkOliveGreen",
    "darkOrange",
    "darkOrchid",
    "darkRed",
    "darkSalmon",
    "darkSeaGreen",
    "darkSlateBlue",
    "darkSlateGray",
    "darkSlateGrey",
    "darkTurquoise",
    "darkViolet",
    "dkBlue",
    "dkCyan",
    "dkGoldenrod",
    "dkGray",
    "dkGrey",
    "dkGreen",
    "dkKhaki",
    "dkMagenta",
    "dkOliveGreen",
    "dkOrange",
    "dkOrchid",
    "dkRed",
    "dkSalmon",
    "dkSeaGreen",
    "dkSlateBlue",
    "dkSlateGray",
    "dkSlateGrey",
    "dkTurquoise",
    "dkViolet",
    "deepPink",
    "deepSkyBlue",
    "dimGray",
    "dimGrey",
    "dodgerBlue",
    "firebrick",
    "floralWhite",
    "forestGreen",
    "fuchsia",
    "gainsboro",
    "ghostWhite",
    "gold",
    "goldenrod",
    "gray",
    "grey",
    "green",
    "greenYellow",
    "honeydew",
    "hotPink",
    "indianRed",
    "indigo",
    "ivory",
    "khaki",
    "lavender",
    "lavenderBlush",
    "lawnGreen",
    "lemonChiffon",
    "lightBlue",
    "lightCoral",
    "lightCyan",
    "lightGoldenrodYellow",
    "lightGray",
    "lightGrey",
    "lightGreen",
    "lightPink",
    "lightSalmon",
    "lightSeaGreen",
    "lightSkyBlue",
    "lightSlateGray",
    "lightSlateGrey",
    "lightSteelBlue",
    "lightYellow",
    "ltBlue",
    "ltCoral",
    "ltCyan",
    "ltGoldenrodYellow",
    "ltGray",
    "ltGrey",
    "ltGreen",
    "ltPink",
    "ltSalmon",
    "ltSeaGreen",
    "ltSkyBlue",
    "ltSlateGray",
    "ltSlateGrey",
    "ltSteelBlue",
    "ltYellow",
    "lime",
    "limeGreen",
    "linen",
    "magenta",
    "maroon",
    "medAquamarine",
    "medBlue",
    "medOrchid",
    "medPurple",
    "medSeaGreen",
    "medSlateBlue",
    "medSpringGreen",
    "medTurquoise",
    "medVioletRed",
    "mediumAquamarine",
    "mediumBlue",
    "mediumOrchid",
    "mediumPurple",
    "mediumSeaGreen",
    "mediumSlateBlue",
    "mediumSpringGreen",
    "mediumTurquoise",
    "mediumVioletRed",
    "midnightBlue",
    "mintCream",
    "mistyRose",
    "moccasin",
    "navajoWhite",
    "navy",
    "oldLace",
    "olive",
    "oliveDrab",
    "orange",
    "orangeRed",
    "orchid",
    "paleGoldenrod",
    "paleGreen",
    "paleTurquoise",
    "paleVioletRed",
    "papayaWhip",
    "peachPuff",
    "peru",
    "pink",
    "plum",
    "powderBlue",
    "purple",
    "red",
    "rosyBrown",
    "royalBlue",
    "saddleBrown",
    "salmon",
    "sandyBrown",
    "seaGreen",
    "seaShell",
    "sienna",
    "silver",
    "skyBlue",
    "slateBlue",
    "slateGray",
    "slateGrey",
    "snow",
    "springGreen",
    "steelBlue",
    "tan",
    "teal",
    "thistle",
    "tomato",
    "turquoise",
    "violet",
    "wheat",
    "white",
    "whiteSmoke",
    "yellow",
    "yellowGreen",
]

PRESET_COLORS: Final[list[_PresetColors]]
SCHEME_COLORS: Final[list[_SchemeColors]]

class Transform(Serialisable): ...

class SystemColor(Serialisable):
    tagname: ClassVar[str]
    namespace: ClassVar[str]
    tint: NestedInteger[Literal[True]]
    shade: NestedInteger[Literal[True]]
    comp: Typed[Transform, Literal[True]]
    inv: Typed[Transform, Literal[True]]
    gray: Typed[Transform, Literal[True]]
    alpha: NestedInteger[Literal[True]]
    alphaOff: NestedInteger[Literal[True]]
    alphaMod: NestedInteger[Literal[True]]
    hue: NestedInteger[Literal[True]]
    hueOff: NestedInteger[Literal[True]]
    hueMod: NestedInteger[Literal[True]]
    sat: NestedInteger[Literal[True]]
    satOff: NestedInteger[Literal[True]]
    satMod: NestedInteger[Literal[True]]
    lum: NestedInteger[Literal[True]]
    lumOff: NestedInteger[Literal[True]]
    lumMod: NestedInteger[Literal[True]]
    red: NestedInteger[Literal[True]]
    redOff: NestedInteger[Literal[True]]
    redMod: NestedInteger[Literal[True]]
    green: NestedInteger[Literal[True]]
    greenOff: NestedInteger[Literal[True]]
    greenMod: NestedInteger[Literal[True]]
    blue: NestedInteger[Literal[True]]
    blueOff: NestedInteger[Literal[True]]
    blueMod: NestedInteger[Literal[True]]
    gamma: Typed[Transform, Literal[True]]
    invGamma: Typed[Transform, Literal[True]]
    val: Set[_SystemColorVal]
    lastClr: Incomplete
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        val: _SystemColorVal = "windowText",
        lastClr=None,
        tint: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        shade: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        comp: Transform | None = None,
        inv: Transform | None = None,
        gray: Transform | None = None,
        alpha: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        alphaOff: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        alphaMod: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        hue: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        hueOff: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        hueMod: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        sat: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        satOff: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        satMod: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        lum: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        lumOff: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        lumMod: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        red: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        redOff: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        redMod: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        green: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        greenOff: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        greenMod: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        blue: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        blueOff: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        blueMod: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        gamma: Transform | None = None,
        invGamma: Transform | None = None,
    ) -> None: ...

class HSLColor(Serialisable):
    tagname: ClassVar[str]
    hue: Integer[Literal[False]]
    sat: MinMax[float, Literal[False]]
    lum: MinMax[float, Literal[False]]
    def __init__(self, hue: ConvertibleToInt, sat: ConvertibleToFloat, lum: ConvertibleToFloat) -> None: ...

class RGBPercent(Serialisable):
    tagname: ClassVar[str]
    r: MinMax[float, Literal[False]]
    g: MinMax[float, Literal[False]]
    b: MinMax[float, Literal[False]]
    def __init__(self, r: ConvertibleToFloat, g: ConvertibleToFloat, b: ConvertibleToFloat) -> None: ...

_RGBPercent: TypeAlias = RGBPercent

class SchemeColor(Serialisable):
    tagname: ClassVar[str]
    namespace: ClassVar[str]
    tint: NestedInteger[Literal[True]]
    shade: NestedInteger[Literal[True]]
    comp: EmptyTag[Literal[True]]
    inv: NestedInteger[Literal[True]]
    gray: NestedInteger[Literal[True]]
    alpha: NestedInteger[Literal[True]]
    alphaOff: NestedInteger[Literal[True]]
    alphaMod: NestedInteger[Literal[True]]
    hue: NestedInteger[Literal[True]]
    hueOff: NestedInteger[Literal[True]]
    hueMod: NestedInteger[Literal[True]]
    sat: NestedInteger[Literal[True]]
    satOff: NestedInteger[Literal[True]]
    satMod: NestedInteger[Literal[True]]
    lum: NestedInteger[Literal[True]]
    lumOff: NestedInteger[Literal[True]]
    lumMod: NestedInteger[Literal[True]]
    red: NestedInteger[Literal[True]]
    redOff: NestedInteger[Literal[True]]
    redMod: NestedInteger[Literal[True]]
    green: NestedInteger[Literal[True]]
    greenOff: NestedInteger[Literal[True]]
    greenMod: NestedInteger[Literal[True]]
    blue: NestedInteger[Literal[True]]
    blueOff: NestedInteger[Literal[True]]
    blueMod: NestedInteger[Literal[True]]
    gamma: EmptyTag[Literal[True]]
    invGamma: EmptyTag[Literal[True]]
    val: Set[_SchemeColors]
    __elements__: ClassVar[tuple[str, ...]]
    @overload
    def __init__(
        self,
        tint: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        shade: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        comp: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        inv: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        gray: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        alpha: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        alphaOff: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        alphaMod: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        hue: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        hueOff: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        hueMod: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        sat: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        satOff: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        satMod: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        lum: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        lumOff: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        lumMod: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        red: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        redOff: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        redMod: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        green: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        greenOff: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        greenMod: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        blue: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        blueOff: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        blueMod: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        gamma: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        invGamma: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        *,
        val: _SchemeColors,
    ) -> None: ...
    @overload
    def __init__(
        self,
        tint: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None,
        shade: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None,
        comp: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None,
        inv: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None,
        gray: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None,
        alpha: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None,
        alphaOff: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None,
        alphaMod: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None,
        hue: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None,
        hueOff: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None,
        hueMod: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None,
        sat: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None,
        satOff: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None,
        satMod: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None,
        lum: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None,
        lumOff: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None,
        lumMod: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None,
        red: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None,
        redOff: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None,
        redMod: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None,
        green: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None,
        greenOff: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None,
        greenMod: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None,
        blue: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None,
        blueOff: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None,
        blueMod: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None,
        gamma: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None,
        invGamma: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None,
        val: _SchemeColors,
    ) -> None: ...

class ColorChoice(Serialisable):
    tagname: ClassVar[str]
    namespace: ClassVar[str]
    scrgbClr: Typed[_RGBPercent, Literal[True]]
    RGBPercent: Alias
    srgbClr: NestedValue[_RGBPercent, Literal[True]]
    RGB: Alias
    hslClr: Typed[HSLColor, Literal[True]]
    sysClr: Typed[SystemColor, Literal[True]]
    schemeClr: Typed[SystemColor, Literal[True]]
    prstClr: NestedNoneSet[_PresetColors]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        scrgbClr: _RGBPercent | None = None,
        srgbClr: _HasTagAndGet[_RGBPercent | None] | _RGBPercent | None = None,
        hslClr: HSLColor | None = None,
        sysClr: SystemColor | None = None,
        schemeClr: SystemColor | None = None,
        prstClr: _NestedNoneSetParam[_PresetColors] = None,
    ) -> None: ...

_COLOR_SET: Final[tuple[_ColorSetType, ...]]

class ColorMapping(Serialisable):
    tagname: ClassVar[str]
    bg1: Set[_ColorSetType]
    tx1: Set[_ColorSetType]
    bg2: Set[_ColorSetType]
    tx2: Set[_ColorSetType]
    accent1: Set[_ColorSetType]
    accent2: Set[_ColorSetType]
    accent3: Set[_ColorSetType]
    accent4: Set[_ColorSetType]
    accent5: Set[_ColorSetType]
    accent6: Set[_ColorSetType]
    hlink: Set[_ColorSetType]
    folHlink: Set[_ColorSetType]
    extLst: Typed[ExtensionList, Literal[True]]
    def __init__(
        self,
        bg1: str = "lt1",
        tx1: str = "dk1",
        bg2: str = "lt2",
        tx2: str = "dk2",
        accent1: str = "accent1",
        accent2: str = "accent2",
        accent3: str = "accent3",
        accent4: str = "accent4",
        accent5: str = "accent5",
        accent6: str = "accent6",
        hlink: str = "hlink",
        folHlink: str = "folHlink",
        extLst: ExtensionList | None = None,
    ) -> None: ...

class ColorChoiceDescriptor(Typed[ColorChoice, Literal[True]]):
    expected_type: type[ColorChoice]
    allow_none: Literal[True]
    def __init__(self, name: str | None = None) -> None: ...
    def __set__(self, instance: Serialisable | Strict, value: str | ColorChoice | None) -> None: ...
