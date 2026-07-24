from _typeshed import ConvertibleToFloat, ConvertibleToInt, Incomplete, Unused
from typing import ClassVar, Literal, overload
from typing_extensions import TypeAlias

from openpyxl.descriptors.base import Alias, Bool, Float, Integer, MinMax, NoneSet, Set, String, Typed, _ConvertibleToBool
from openpyxl.descriptors.excel import Coordinate, ExtensionList, Percentage
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.styles.colors import Color

_BevelPrst: TypeAlias = Literal[
    "relaxedInset",
    "circle",
    "slope",
    "cross",
    "angle",
    "softRound",
    "convex",
    "coolSlant",
    "divot",
    "riblet",
    "hardEdge",
    "artDeco",
]
_Shape3DPrstMaterial: TypeAlias = Literal[
    "legacyMatte",
    "legacyPlastic",
    "legacyMetal",
    "legacyWireframe",
    "matte",
    "plastic",
    "metal",
    "warmMatte",
    "translucentPowder",
    "powder",
    "dkEdge",
    "softEdge",
    "clear",
    "flat",
    "softmetal",
]
_Path2DFill: TypeAlias = Literal["norm", "lighten", "lightenLess", "darken", "darkenLess"]
_FontReferenceIdx: TypeAlias = Literal["major", "minor"]
_CameraPrst: TypeAlias = Literal[
    "legacyObliqueTopLeft",
    "legacyObliqueTop",
    "legacyObliqueTopRight",
    "legacyObliqueLeft",
    "legacyObliqueFront",
    "legacyObliqueRight",
    "legacyObliqueBottomLeft",
    "legacyObliqueBottom",
    "legacyObliqueBottomRight",
    "legacyPerspectiveTopLeft",
    "legacyPerspectiveTop",
    "legacyPerspectiveTopRight",
    "legacyPerspectiveLeft",
    "legacyPerspectiveFront",
    "legacyPerspectiveRight",
    "legacyPerspectiveBottomLeft",
    "legacyPerspectiveBottom",
    "legacyPerspectiveBottomRight",
    "orthographicFront",
    "isometricTopUp",
    "isometricTopDown",
    "isometricBottomUp",
    "isometricBottomDown",
    "isometricLeftUp",
    "isometricLeftDown",
    "isometricRightUp",
    "isometricRightDown",
    "isometricOffAxis1Left",
    "isometricOffAxis1Right",
    "isometricOffAxis1Top",
    "isometricOffAxis2Left",
    "isometricOffAxis2Right",
    "isometricOffAxis2Top",
    "isometricOffAxis3Left",
    "isometricOffAxis3Right",
    "isometricOffAxis3Bottom",
    "isometricOffAxis4Left",
    "isometricOffAxis4Right",
    "isometricOffAxis4Bottom",
    "obliqueTopLeft",
    "obliqueTop",
    "obliqueTopRight",
    "obliqueLeft",
    "obliqueRight",
    "obliqueBottomLeft",
    "obliqueBottom",
    "obliqueBottomRight",
    "perspectiveFront",
    "perspectiveLeft",
    "perspectiveRight",
    "perspectiveAbove",
    "perspectiveBelow",
    "perspectiveAboveLeftFacing",
    "perspectiveAboveRightFacing",
    "perspectiveContrastingLeftFacing",
    "perspectiveContrastingRightFacing",
    "perspectiveHeroicLeftFacing",
    "perspectiveHeroicRightFacing",
    "perspectiveHeroicExtremeLeftFacing",
    "perspectiveHeroicExtremeRightFacing",
    "perspectiveRelaxed",
    "perspectiveRelaxedModerately",
]
_LightRigRig: TypeAlias = Literal[
    "legacyFlat1",
    "legacyFlat2",
    "legacyFlat3",
    "legacyFlat4",
    "legacyNormal1",
    "legacyNormal2",
    "legacyNormal3",
    "legacyNormal4",
    "legacyHarsh1",
    "legacyHarsh2",
    "legacyHarsh3",
    "legacyHarsh4",
    "threePt",
    "balanced",
    "soft",
    "harsh",
    "flood",
    "contrasting",
    "morning",
    "sunrise",
    "sunset",
    "chilly",
    "freezing",
    "flat",
    "twoPt",
    "glow",
    "brightRoom",
]
_LightRigDir: TypeAlias = Literal["tl", "t", "tr", "l", "r", "bl", "b", "br"]
_PresetGeometry2DPrst: TypeAlias = Literal[
    "line",
    "lineInv",
    "triangle",
    "rtTriangle",
    "rect",
    "diamond",
    "parallelogram",
    "trapezoid",
    "nonIsoscelesTrapezoid",
    "pentagon",
    "hexagon",
    "heptagon",
    "octagon",
    "decagon",
    "dodecagon",
    "star4",
    "star5",
    "star6",
    "star7",
    "star8",
    "star10",
    "star12",
    "star16",
    "star24",
    "star32",
    "roundRect",
    "round1Rect",
    "round2SameRect",
    "round2DiagRect",
    "snipRoundRect",
    "snip1Rect",
    "snip2SameRect",
    "snip2DiagRect",
    "plaque",
    "ellipse",
    "teardrop",
    "homePlate",
    "chevron",
    "pieWedge",
    "pie",
    "blockArc",
    "donut",
    "noSmoking",
    "rightArrow",
    "leftArrow",
    "upArrow",
    "downArrow",
    "stripedRightArrow",
    "notchedRightArrow",
    "bentUpArrow",
    "leftRightArrow",
    "upDownArrow",
    "leftUpArrow",
    "leftRightUpArrow",
    "quadArrow",
    "leftArrowCallout",
    "rightArrowCallout",
    "upArrowCallout",
    "downArrowCallout",
    "leftRightArrowCallout",
    "upDownArrowCallout",
    "quadArrowCallout",
    "bentArrow",
    "uturnArrow",
    "circularArrow",
    "leftCircularArrow",
    "leftRightCircularArrow",
    "curvedRightArrow",
    "curvedLeftArrow",
    "curvedUpArrow",
    "curvedDownArrow",
    "swooshArrow",
    "cube",
    "can",
    "lightningBolt",
    "heart",
    "sun",
    "moon",
    "smileyFace",
    "irregularSeal1",
    "irregularSeal2",
    "foldedCorner",
    "bevel",
    "frame",
    "halfFrame",
    "corner",
    "diagStripe",
    "chord",
    "arc",
    "leftBracket",
    "rightBracket",
    "leftBrace",
    "rightBrace",
    "bracketPair",
    "bracePair",
    "straightConnector1",
    "bentConnector2",
    "bentConnector3",
    "bentConnector4",
    "bentConnector5",
    "curvedConnector2",
    "curvedConnector3",
    "curvedConnector4",
    "curvedConnector5",
    "callout1",
    "callout2",
    "callout3",
    "accentCallout1",
    "accentCallout2",
    "accentCallout3",
    "borderCallout1",
    "borderCallout2",
    "borderCallout3",
    "accentBorderCallout1",
    "accentBorderCallout2",
    "accentBorderCallout3",
    "wedgeRectCallout",
    "wedgeRoundRectCallout",
    "wedgeEllipseCallout",
    "cloudCallout",
    "cloud",
    "ribbon",
    "ribbon2",
    "ellipseRibbon",
    "ellipseRibbon2",
    "leftRightRibbon",
    "verticalScroll",
    "horizontalScroll",
    "wave",
    "doubleWave",
    "plus",
    "flowChartProcess",
    "flowChartDecision",
    "flowChartInputOutput",
    "flowChartPredefinedProcess",
    "flowChartInternalStorage",
    "flowChartDocument",
    "flowChartMultidocument",
    "flowChartTerminator",
    "flowChartPreparation",
    "flowChartManualInput",
    "flowChartManualOperation",
    "flowChartConnector",
    "flowChartPunchedCard",
    "flowChartPunchedTape",
    "flowChartSummingJunction",
    "flowChartOr",
    "flowChartCollate",
    "flowChartSort",
    "flowChartExtract",
    "flowChartMerge",
    "flowChartOfflineStorage",
    "flowChartOnlineStorage",
    "flowChartMagneticTape",
    "flowChartMagneticDisk",
    "flowChartMagneticDrum",
    "flowChartDisplay",
    "flowChartDelay",
    "flowChartAlternateProcess",
    "flowChartOffpageConnector",
    "actionButtonBlank",
    "actionButtonHome",
    "actionButtonHelp",
    "actionButtonInformation",
    "actionButtonForwardNext",
    "actionButtonBackPrevious",
    "actionButtonEnd",
    "actionButtonBeginning",
    "actionButtonReturn",
    "actionButtonDocument",
    "actionButtonSound",
    "actionButtonMovie",
    "gear6",
    "gear9",
    "funnel",
    "mathPlus",
    "mathMinus",
    "mathMultiply",
    "mathDivide",
    "mathEqual",
    "mathNotEqual",
    "cornerTabs",
    "squareTabs",
    "plaqueTabs",
    "chartX",
    "chartStar",
    "chartPlus",
]

class Point2D(Serialisable):
    tagname: ClassVar[str]
    namespace: ClassVar[str]
    x: Incomplete
    y: Incomplete
    def __init__(self, x=None, y=None) -> None: ...

class PositiveSize2D(Serialisable):
    tagname: ClassVar[str]
    namespace: ClassVar[str]
    cx: Integer[Literal[False]]
    width: Alias
    cy: Integer[Literal[False]]
    height: Alias
    def __init__(self, cx: ConvertibleToInt, cy: ConvertibleToInt) -> None: ...

class Transform2D(Serialisable):
    tagname: ClassVar[str]
    namespace: ClassVar[str]
    rot: Integer[Literal[True]]
    flipH: Bool[Literal[True]]
    flipV: Bool[Literal[True]]
    off: Typed[Point2D, Literal[True]]
    ext: Typed[PositiveSize2D, Literal[True]]
    chOff: Typed[Point2D, Literal[True]]
    chExt: Typed[PositiveSize2D, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        rot: ConvertibleToInt | None = None,
        flipH: _ConvertibleToBool | None = None,
        flipV: _ConvertibleToBool | None = None,
        off: Point2D | None = None,
        ext: PositiveSize2D | None = None,
        chOff: Point2D | None = None,
        chExt: PositiveSize2D | None = None,
    ) -> None: ...

class GroupTransform2D(Serialisable):
    tagname: ClassVar[str]
    namespace: ClassVar[str]
    rot: Integer[Literal[True]]
    flipH: Bool[Literal[True]]
    flipV: Bool[Literal[True]]
    off = Typed(expected_type=Point2D, allow_none=True)
    ext = Typed(expected_type=PositiveSize2D, allow_none=True)
    chOff = Typed(expected_type=Point2D, allow_none=True)
    chExt = Typed(expected_type=PositiveSize2D, allow_none=True)
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        rot: ConvertibleToInt | None = 0,
        flipH: _ConvertibleToBool | None = None,
        flipV: _ConvertibleToBool | None = None,
        off: Point2D | None = None,
        ext: PositiveSize2D | None = None,
        chOff: Point2D | None = None,
        chExt: PositiveSize2D | None = None,
    ) -> None: ...

class SphereCoords(Serialisable):
    tagname: ClassVar[str]
    lat: Integer[Literal[False]]
    lon: Integer[Literal[False]]
    rev: Integer[Literal[False]]
    def __init__(self, lat: ConvertibleToInt, lon: ConvertibleToInt, rev: ConvertibleToInt) -> None: ...

class Camera(Serialisable):
    tagname: ClassVar[str]
    prst: Set[_CameraPrst]
    fov: Integer[Literal[True]]
    zoom: Typed[Percentage, Literal[True]]
    rot: Typed[SphereCoords, Literal[True]]
    def __init__(
        self,
        prst: _CameraPrst,
        fov: ConvertibleToInt | None = None,
        zoom: Percentage | None = None,
        rot: SphereCoords | None = None,
    ) -> None: ...

class LightRig(Serialisable):
    tagname: ClassVar[str]
    rig: Set[_LightRigRig]
    dir: Set[_LightRigDir]
    rot: Typed[SphereCoords, Literal[True]]
    def __init__(self, rig: _LightRigRig, dir: _LightRigDir, rot: SphereCoords | None = None) -> None: ...

class Vector3D(Serialisable):
    tagname: ClassVar[str]
    dx: Integer[Literal[False]]
    dy: Integer[Literal[False]]
    dz: Integer[Literal[False]]
    def __init__(self, dx: ConvertibleToInt, dy: ConvertibleToInt, dz: ConvertibleToInt) -> None: ...

class Point3D(Serialisable):
    tagname: ClassVar[str]
    x: Integer[Literal[False]]
    y: Integer[Literal[False]]
    z: Integer[Literal[False]]
    def __init__(self, x: ConvertibleToInt, y: ConvertibleToInt, z: ConvertibleToInt) -> None: ...

class Backdrop(Serialisable):
    anchor: Typed[Point3D, Literal[False]]
    norm: Typed[Vector3D, Literal[False]]
    up: Typed[Vector3D, Literal[False]]
    extLst: Typed[ExtensionList, Literal[True]]
    def __init__(self, anchor: Point3D, norm: Vector3D, up: Vector3D, extLst: ExtensionList | None = None) -> None: ...

class Scene3D(Serialisable):
    camera: Typed[Camera, Literal[False]]
    lightRig: Typed[LightRig, Literal[False]]
    backdrop: Typed[Backdrop, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    def __init__(
        self, camera: Camera, lightRig: LightRig, backdrop: Backdrop | None = None, extLst: ExtensionList | None = None
    ) -> None: ...

class Bevel(Serialisable):
    tagname: ClassVar[str]
    w: Integer[Literal[False]]
    h: Integer[Literal[False]]
    prst: NoneSet[_BevelPrst]
    def __init__(self, w: ConvertibleToInt, h: ConvertibleToInt, prst: _BevelPrst | Literal["none"] | None = None) -> None: ...

class Shape3D(Serialisable):
    namespace: ClassVar[str]
    z: Typed[Coordinate[bool], Literal[True]]
    extrusionH: Integer[Literal[True]]
    contourW: Integer[Literal[True]]
    prstMaterial: NoneSet[_Shape3DPrstMaterial]
    bevelT: Typed[Bevel, Literal[True]]
    bevelB: Typed[Bevel, Literal[True]]
    extrusionClr: Typed[Color, Literal[True]]
    contourClr: Typed[Color, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    def __init__(
        self,
        z: Coordinate[bool] | None = None,
        extrusionH: ConvertibleToInt | None = None,
        contourW: ConvertibleToInt | None = None,
        prstMaterial: _Shape3DPrstMaterial | Literal["none"] | None = None,
        bevelT: Bevel | None = None,
        bevelB: Bevel | None = None,
        extrusionClr: Color | None = None,
        contourClr: Color | None = None,
        extLst: ExtensionList | None = None,
    ) -> None: ...

class Path2D(Serialisable):
    w: Float[Literal[False]]
    h: Float[Literal[False]]
    fill: NoneSet[_Path2DFill]
    stroke: Bool[Literal[True]]
    extrusionOk: Bool[Literal[True]]
    def __init__(
        self,
        w: ConvertibleToFloat,
        h: ConvertibleToFloat,
        fill: _Path2DFill | Literal["none"] | None = None,
        stroke: _ConvertibleToBool | None = None,
        extrusionOk: _ConvertibleToBool | None = None,
    ) -> None: ...

class Path2DList(Serialisable):
    path: Typed[Path2D, Literal[True]]
    def __init__(self, path: Path2D | None = None) -> None: ...

class GeomRect(Serialisable):
    l: Incomplete
    t: Incomplete
    r: Incomplete
    b: Incomplete
    def __init__(self, l=None, t=None, r=None, b=None) -> None: ...

class AdjPoint2D(Serialisable):
    x: Incomplete
    y: Incomplete
    def __init__(self, x=None, y=None) -> None: ...

class ConnectionSite(Serialisable):
    ang: MinMax[float, Literal[False]]
    pos: Typed[AdjPoint2D, Literal[False]]
    def __init__(self, ang: ConvertibleToFloat, pos: AdjPoint2D) -> None: ...

class ConnectionSiteList(Serialisable):
    cxn: Typed[ConnectionSite, Literal[True]]
    def __init__(self, cxn: ConnectionSite | None = None) -> None: ...

class AdjustHandleList(Serialisable): ...

class GeomGuide(Serialisable):
    name: String[Literal[False]]
    fmla: String[Literal[False]]
    def __init__(self, name: str, fmla: str) -> None: ...

class GeomGuideList(Serialisable):
    gd: Typed[GeomGuide, Literal[True]]
    def __init__(self, gd: GeomGuide | None = None) -> None: ...

class CustomGeometry2D(Serialisable):
    avLst: Typed[GeomGuideList, Literal[True]]
    gdLst: Typed[GeomGuideList, Literal[True]]
    ahLst: Typed[AdjustHandleList, Literal[True]]
    cxnLst: Typed[ConnectionSiteList, Literal[True]]
    pathLst: Typed[Path2DList, Literal[False]]
    rect: GeomRect | None
    @overload
    def __init__(
        self,
        avLst: GeomGuideList | None = None,
        gdLst: GeomGuideList | None = None,
        ahLst: AdjustHandleList | None = None,
        cxnLst: ConnectionSiteList | None = None,
        rect: Unused = None,
        *,
        pathLst: Path2DList,
    ) -> None: ...
    @overload
    def __init__(
        self,
        avLst: GeomGuideList | None,
        gdLst: GeomGuideList | None,
        ahLst: AdjustHandleList | None,
        cxnLst: ConnectionSiteList | None,
        rect: Unused,
        pathLst: Path2DList,
    ) -> None: ...

class PresetGeometry2D(Serialisable):
    namespace: ClassVar[str]
    prst: Set[_PresetGeometry2DPrst]
    avLst: Typed[GeomGuideList, Literal[True]]
    def __init__(self, prst: _PresetGeometry2DPrst, avLst: GeomGuideList | None = None) -> None: ...

class FontReference(Serialisable):
    idx: NoneSet[_FontReferenceIdx]
    def __init__(self, idx: _FontReferenceIdx | Literal["none"] | None = None) -> None: ...

class StyleMatrixReference(Serialisable):
    idx: Integer[Literal[False]]
    def __init__(self, idx: ConvertibleToInt) -> None: ...

class ShapeStyle(Serialisable):
    lnRef: Typed[StyleMatrixReference, Literal[False]]
    fillRef: Typed[StyleMatrixReference, Literal[False]]
    effectRef: Typed[StyleMatrixReference, Literal[False]]
    fontRef: Typed[FontReference, Literal[False]]
    def __init__(
        self, lnRef: StyleMatrixReference, fillRef: StyleMatrixReference, effectRef: StyleMatrixReference, fontRef: FontReference
    ) -> None: ...
