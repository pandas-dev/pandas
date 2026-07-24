from collections.abc import Iterable, Iterator
from typing import Final, Literal, TypeVar, overload, type_check_only
from typing_extensions import Self, TypeAlias

_ColorT = TypeVar("_ColorT", bound=Color)
# NOTE: Reportlab is very inconsistent and sometimes uses the interpretation
#       used in reportlab.pdfgen.textobject instead, so we pick a different name
_ConvertibleToColor: TypeAlias = Color | list[float] | tuple[float, float, float, float] | tuple[float, float, float] | str | int

__version__: Final[str]

class Color:
    red: float
    green: float
    blue: float
    alpha: float
    def __init__(self, red: float = 0, green: float = 0, blue: float = 0, alpha: float = 1) -> None: ...
    @property
    def __key__(self) -> tuple[float, ...]: ...
    def __hash__(self) -> int: ...
    def __comparable__(self, other: object) -> bool: ...
    def __lt__(self, other: object) -> bool: ...
    def __eq__(self, other: object) -> bool: ...
    def __le__(self, other: object) -> bool: ...
    def __gt__(self, other: object) -> bool: ...
    def __ge__(self, other: object) -> bool: ...
    def rgb(self) -> tuple[float, float, float]: ...
    def rgba(self) -> tuple[float, float, float, float]: ...
    def bitmap_rgb(self) -> tuple[int, int, int]: ...
    def bitmap_rgba(self) -> tuple[int, int, int, int]: ...
    def hexval(self) -> str: ...
    def hexvala(self) -> str: ...
    def int_rgb(self) -> int: ...
    def int_rgba(self) -> int: ...
    def int_argb(self) -> int: ...
    @property
    def cKwds(self) -> Iterator[tuple[str, int]]: ...
    # NOTE: Possible arguments depend on __init__, so this violates LSP
    #       For now we just leave it unchecked
    def clone(self, **kwds) -> Self: ...
    @property
    def normalizedAlpha(self) -> float: ...

def opaqueColor(c: object) -> bool: ...

class CMYKColor(Color):
    cyan: float
    magenta: float
    yellow: float
    black: float
    spotName: str | None
    density: float
    knockout: bool | None
    alpha: float
    def __init__(
        self,
        cyan: float = 0,
        magenta: float = 0,
        yellow: float = 0,
        black: float = 0,
        spotName: str | None = None,
        density: float = 1,
        knockout: bool | None = None,
        alpha: float = 1,
    ) -> None: ...
    def fader(self, n: int, reverse: bool = False) -> list[Self]: ...
    def cmyk(self) -> tuple[float, float, float, float]: ...
    def cmyka(self) -> tuple[float, float, float, float, float]: ...

class PCMYKColor(CMYKColor):
    def __init__(
        self,
        cyan: float,
        magenta: float,
        yellow: float,
        black: float,
        density: float = 100,
        spotName: str | None = None,
        knockout: bool | None = None,
        alpha: float = 100,
    ) -> None: ...

class CMYKColorSep(CMYKColor):
    def __init__(
        self,
        cyan: float = 0,
        magenta: float = 0,
        yellow: float = 0,
        black: float = 0,
        spotName: str | None = None,
        density: float = 1,
        alpha: float = 1,
    ) -> None: ...

class PCMYKColorSep(PCMYKColor, CMYKColorSep):
    def __init__(
        self,
        cyan: float = 0,
        magenta: float = 0,
        yellow: float = 0,
        black: float = 0,
        spotName: str | None = None,
        density: float = 100,
        alpha: float = 100,
    ) -> None: ...

def cmyk2rgb(cmyk: tuple[float, float, float, float], density: float = 1) -> tuple[float, float, float]: ...
def rgb2cmyk(r: float, g: float, b: float) -> tuple[float, float, float, float]: ...
def color2bw(colorRGB: Color) -> Color: ...
def HexColor(val: str | int, htmlOnly: bool = False, hasAlpha: bool = False) -> Color: ...
def linearlyInterpolatedColor(c0: _ColorT, c1: _ColorT, x0: float, x1: float, x: float) -> _ColorT: ...
@overload
def obj_R_G_B(
    c: Color | list[float] | tuple[float, float, float, float] | tuple[float, float, float],
) -> tuple[float, float, float]: ...
@overload
def obj_R_G_B(c: None) -> None: ...

transparent: Color
ReportLabBlueOLD: Color
ReportLabBlue: Color
ReportLabBluePCMYK: Color
ReportLabLightBlue: Color
ReportLabFidBlue: Color
ReportLabFidRed: Color
ReportLabGreen: Color
ReportLabLightGreen: Color
aliceblue: Color
antiquewhite: Color
aqua: Color
aquamarine: Color
azure: Color
beige: Color
bisque: Color
black: Color
blanchedalmond: Color
blue: Color
blueviolet: Color
brown: Color
burlywood: Color
cadetblue: Color
chartreuse: Color
chocolate: Color
coral: Color
cornflowerblue: Color
cornflower: Color
cornsilk: Color
crimson: Color
cyan: Color
darkblue: Color
darkcyan: Color
darkgoldenrod: Color
darkgray: Color
darkgrey: Color
darkgreen: Color
darkkhaki: Color
darkmagenta: Color
darkolivegreen: Color
darkorange: Color
darkorchid: Color
darkred: Color
darksalmon: Color
darkseagreen: Color
darkslateblue: Color
darkslategray: Color
darkslategrey: Color
darkturquoise: Color
darkviolet: Color
deeppink: Color
deepskyblue: Color
dimgray: Color
dimgrey: Color
dodgerblue: Color
firebrick: Color
floralwhite: Color
forestgreen: Color
fuchsia: Color
gainsboro: Color
ghostwhite: Color
gold: Color
goldenrod: Color
gray: Color
grey: Color
green: Color
greenyellow: Color
honeydew: Color
hotpink: Color
indianred: Color
indigo: Color
ivory: Color
khaki: Color
lavender: Color
lavenderblush: Color
lawngreen: Color
lemonchiffon: Color
lightblue: Color
lightcoral: Color
lightcyan: Color
lightgoldenrodyellow: Color
lightgreen: Color
lightgrey: Color
lightpink: Color
lightsalmon: Color
lightseagreen: Color
lightskyblue: Color
lightslategray: Color
lightslategrey: Color
lightsteelblue: Color
lightyellow: Color
lime: Color
limegreen: Color
linen: Color
magenta: Color
maroon: Color
mediumaquamarine: Color
mediumblue: Color
mediumorchid: Color
mediumpurple: Color
mediumseagreen: Color
mediumslateblue: Color
mediumspringgreen: Color
mediumturquoise: Color
mediumvioletred: Color
midnightblue: Color
mintcream: Color
mistyrose: Color
moccasin: Color
navajowhite: Color
navy: Color
oldlace: Color
olive: Color
olivedrab: Color
orange: Color
orangered: Color
orchid: Color
palegoldenrod: Color
palegreen: Color
paleturquoise: Color
palevioletred: Color
papayawhip: Color
peachpuff: Color
peru: Color
pink: Color
plum: Color
powderblue: Color
purple: Color
red: Color
rosybrown: Color
royalblue: Color
saddlebrown: Color
salmon: Color
sandybrown: Color
seagreen: Color
seashell: Color
sienna: Color
silver: Color
skyblue: Color
slateblue: Color
slategray: Color
slategrey: Color
snow: Color
springgreen: Color
steelblue: Color
tan: Color
teal: Color
thistle: Color
tomato: Color
turquoise: Color
violet: Color
wheat: Color
white: Color
whitesmoke: Color
yellow: Color
yellowgreen: Color
fidblue: Color
fidred: Color
fidlightblue: Color
ColorType: type[Color]

def colorDistance(col1: Color, col2: Color) -> float: ...
def cmykDistance(col1: Color, col2: Color) -> float: ...
def getAllNamedColors() -> dict[str, Color]: ...
@overload
def describe(aColor: Color, mode: Literal[0] = 0) -> None: ...
@overload
def describe(aColor: Color, mode: Literal[1]) -> str: ...
@overload
def describe(aColor: Color, mode: Literal[2]) -> tuple[str, float]: ...
def hue2rgb(m1: float, m2: float, h: float) -> float: ...
def hsl2rgb(h: float, s: float, l: float) -> tuple[float, float, float]: ...
@type_check_only
class _cssParse:
    def pcVal(self, v: str) -> float: ...
    def rgbPcVal(self, v: str) -> float: ...
    def rgbVal(self, v: str) -> float: ...
    def hueVal(self, v: str) -> float: ...
    def alphaVal(self, v: str, c: float = 1, n: str = "alpha") -> float: ...
    s: str
    def __call__(self, s: str) -> Color: ...

cssParse: _cssParse

@type_check_only
class _toColor:
    extraColorsNS: dict[str, Color]
    def __init__(self) -> None: ...
    def setExtraColorsNameSpace(self, NS: dict[str, Color]) -> None: ...
    def __call__(self, arg: _ConvertibleToColor, default: Color | None = None) -> Color: ...

toColor: _toColor

@overload
def toColorOrNone(arg: None, default: Color | None) -> None: ...
@overload
def toColorOrNone(arg: _ConvertibleToColor, default: Color | None = None) -> Color: ...
def setColors(**kw: _ConvertibleToColor) -> None: ...
def Whiter(c: _ColorT, f: float) -> _ColorT: ...
def Blacker(c: _ColorT, f: float) -> _ColorT: ...
def fade(aSpotColor: CMYKColor, percentages: Iterable[float]) -> list[CMYKColor]: ...
