from _typeshed import ConvertibleToFloat, ConvertibleToInt, Incomplete, Unused
from typing import ClassVar, Literal
from typing_extensions import TypeAlias

from openpyxl.descriptors.base import Alias, Bool, Integer, MinMax, NoneSet, Set, String, Typed, _ConvertibleToBool
from openpyxl.descriptors.excel import Coordinate, ExtensionList
from openpyxl.descriptors.nested import EmptyTag, NestedBool, NestedInteger, NestedText, NestedValue
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.drawing.colors import ColorChoice, ColorChoiceDescriptor
from openpyxl.drawing.effect import Color, EffectContainer, EffectList
from openpyxl.drawing.fill import Blip, BlipFillProperties, GradientFillProperties, PatternFillProperties
from openpyxl.drawing.geometry import Scene3D
from openpyxl.drawing.line import LineProperties

from ..xml._functions_overloads import _HasTagAndGet

_CharacterPropertiesU: TypeAlias = Literal[
    "words",
    "sng",
    "dbl",
    "heavy",
    "dotted",
    "dottedHeavy",
    "dash",
    "dashHeavy",
    "dashLong",
    "dashLongHeavy",
    "dotDash",
    "dotDashHeavy",
    "dotDotDash",
    "dotDotDashHeavy",
    "wavy",
    "wavyHeavy",
    "wavyDbl",
]
_CharacterPropertiesStrike: TypeAlias = Literal["noStrike", "sngStrike", "dblStrike"]
_CharacterPropertiesCap: TypeAlias = Literal["small", "all"]
_ParagraphPropertiesAlgn: TypeAlias = Literal["l", "ctr", "r", "just", "justLow", "dist", "thaiDist"]
_ParagraphPropertiesFontAlgn: TypeAlias = Literal["auto", "t", "ctr", "base", "b"]
_RichTextPropertiesVertOverflow: TypeAlias = Literal["overflow", "ellipsis", "clip"]
_RichTextPropertiesHorzOverflow: TypeAlias = Literal["overflow", "clip"]
_RichTextPropertiesVert: TypeAlias = Literal[
    "horz", "vert", "vert270", "wordArtVert", "eaVert", "mongolianVert", "wordArtVertRtl"
]
_RichTextPropertiesWrap: TypeAlias = Literal["none", "square"]
_RichTextPropertiesAnchor: TypeAlias = Literal["t", "ctr", "b", "just", "dist"]
_AutonumberBulletType: TypeAlias = Literal[
    "alphaLcParenBoth",
    "alphaUcParenBoth",
    "alphaLcParenR",
    "alphaUcParenR",
    "alphaLcPeriod",
    "alphaUcPeriod",
    "arabicParenBoth",
    "arabicParenR",
    "arabicPeriod",
    "arabicPlain",
    "romanLcParenBoth",
    "romanUcParenBoth",
    "romanLcParenR",
    "romanUcParenR",
    "romanLcPeriod",
    "romanUcPeriod",
    "circleNumDbPlain",
    "circleNumWdBlackPlain",
    "circleNumWdWhitePlain",
    "arabicDbPeriod",
    "arabicDbPlain",
    "ea1ChsPeriod",
    "ea1ChsPlain",
    "ea1ChtPeriod",
    "ea1ChtPlain",
    "ea1JpnChsDbPeriod",
    "ea1JpnKorPlain",
    "ea1JpnKorPeriod",
    "arabic1Minus",
    "arabic2Minus",
    "hebrew2Minus",
    "thaiAlphaPeriod",
    "thaiAlphaParenR",
    "thaiAlphaParenBoth",
    "thaiNumPeriod",
    "thaiNumParenR",
    "thaiNumParenBoth",
    "hindiAlphaPeriod",
    "hindiNumPeriod",
    "hindiNumParenR",
    "hindiAlpha1Period",
]
_TabStopAlgn: TypeAlias = Literal["l", "ctr", "r", "dec"]
_PresetTextShapePrst: TypeAlias = Literal[
    "textNoShape",
    "textPlain",
    "textStop",
    "textTriangle",
    "textTriangleInverted",
    "textChevron",
    "textChevronInverted",
    "textRingInside",
    "textRingOutside",
    "textArchUp",
    "textArchDown",
    "textCircle",
    "textButton",
    "textArchUpPour",
    "textArchDownPour",
    "textCirclePour",
    "textButtonPour",
    "textCurveUp",
    "textCurveDown",
    "textCanUp",
    "textCanDown",
    "textWave1",
    "textWave2",
    "textDoubleWave1",
    "textWave4",
    "textInflate",
    "textDeflate",
    "textInflateBottom",
    "textDeflateBottom",
    "textInflateTop",
    "textDeflateTop",
    "textDeflateInflate",
    "textDeflateInflateDeflate",
    "textFadeRight",
    "textFadeLeft",
    "textFadeUp",
    "textFadeDown",
    "textSlantUp",
    "textSlantDown",
    "textCascadeUp",
    "textCascadeDown",
]

class EmbeddedWAVAudioFile(Serialisable):
    name: String[Literal[True]]
    def __init__(self, name: str | None = None) -> None: ...

class Hyperlink(Serialisable):
    tagname: ClassVar[str]
    namespace: ClassVar[str]
    invalidUrl: String[Literal[True]]
    action: String[Literal[True]]
    tgtFrame: String[Literal[True]]
    tooltip: String[Literal[True]]
    history: Bool[Literal[True]]
    highlightClick: Bool[Literal[True]]
    endSnd: Bool[Literal[True]]
    snd: Typed[EmbeddedWAVAudioFile, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    id: Incomplete
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        invalidUrl: str | None = None,
        action: str | None = None,
        tgtFrame: str | None = None,
        tooltip: str | None = None,
        history: _ConvertibleToBool | None = None,
        highlightClick: _ConvertibleToBool | None = None,
        endSnd: _ConvertibleToBool | None = None,
        snd: EmbeddedWAVAudioFile | None = None,
        extLst: ExtensionList | None = None,
        id=None,
    ) -> None: ...

class Font(Serialisable):
    tagname: ClassVar[str]
    namespace: ClassVar[str]
    typeface: String[Literal[False]]
    panose: Incomplete
    pitchFamily: MinMax[float, Literal[True]]
    charset: Integer[Literal[True]]
    def __init__(
        self, typeface: str, panose=None, pitchFamily: ConvertibleToFloat | None = None, charset: ConvertibleToInt | None = None
    ) -> None: ...

class CharacterProperties(Serialisable):
    tagname: ClassVar[str]
    namespace: ClassVar[str]
    kumimoji: Bool[Literal[True]]
    lang: String[Literal[True]]
    altLang: String[Literal[True]]
    sz: MinMax[float, Literal[True]]
    b: Bool[Literal[True]]
    i: Bool[Literal[True]]
    u: NoneSet[_CharacterPropertiesU]
    strike: NoneSet[_CharacterPropertiesStrike]
    kern: Integer[Literal[True]]
    cap: NoneSet[_CharacterPropertiesCap]
    spc: Integer[Literal[True]]
    normalizeH: Bool[Literal[True]]
    baseline: Integer[Literal[True]]
    noProof: Bool[Literal[True]]
    dirty: Bool[Literal[True]]
    err: Bool[Literal[True]]
    smtClean: Bool[Literal[True]]
    smtId: Integer[Literal[True]]
    bmk: String[Literal[True]]
    ln: Typed[LineProperties, Literal[True]]
    highlight: Typed[Color, Literal[True]]
    latin: Typed[Font, Literal[True]]
    ea: Typed[Font, Literal[True]]
    cs: Typed[Font, Literal[True]]
    sym: Typed[Font, Literal[True]]
    hlinkClick: Typed[Hyperlink, Literal[True]]
    hlinkMouseOver: Typed[Hyperlink, Literal[True]]
    rtl: NestedBool[Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    noFill: EmptyTag[Literal[False]]
    solidFill: ColorChoiceDescriptor
    gradFill: Typed[GradientFillProperties, Literal[True]]
    blipFill: Typed[BlipFillProperties, Literal[True]]
    pattFill: Typed[PatternFillProperties, Literal[True]]
    grpFill: EmptyTag[Literal[False]]
    effectLst: Typed[EffectList, Literal[True]]
    effectDag: Typed[EffectContainer, Literal[True]]
    uLnTx: EmptyTag[Literal[False]]
    uLn: Typed[LineProperties, Literal[True]]
    uFillTx: EmptyTag[Literal[False]]
    uFill: EmptyTag[Literal[False]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        kumimoji: _ConvertibleToBool | None = None,
        lang: str | None = None,
        altLang: str | None = None,
        sz: ConvertibleToFloat | None = None,
        b: _ConvertibleToBool | None = None,
        i: _ConvertibleToBool | None = None,
        u: _CharacterPropertiesU | Literal["none"] | None = None,
        strike: _CharacterPropertiesStrike | Literal["none"] | None = None,
        kern: ConvertibleToInt | None = None,
        cap: _CharacterPropertiesCap | Literal["none"] | None = None,
        spc: ConvertibleToInt | None = None,
        normalizeH: _ConvertibleToBool | None = None,
        baseline: ConvertibleToInt | None = None,
        noProof: _ConvertibleToBool | None = None,
        dirty: _ConvertibleToBool | None = None,
        err: _ConvertibleToBool | None = None,
        smtClean: _ConvertibleToBool | None = None,
        smtId: ConvertibleToInt | None = None,
        bmk: str | None = None,
        ln: LineProperties | None = None,
        highlight: Color | None = None,
        latin: Font | None = None,
        ea: Font | None = None,
        cs: Font | None = None,
        sym: Font | None = None,
        hlinkClick: Hyperlink | None = None,
        hlinkMouseOver: Hyperlink | None = None,
        rtl: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        extLst: Unused = None,
        noFill: _HasTagAndGet[_ConvertibleToBool] | _ConvertibleToBool = None,
        solidFill: str | ColorChoice | None = None,
        gradFill: GradientFillProperties | None = None,
        blipFill: BlipFillProperties | None = None,
        pattFill: PatternFillProperties | None = None,
        grpFill: _HasTagAndGet[_ConvertibleToBool] | _ConvertibleToBool = None,
        effectLst: EffectList | None = None,
        effectDag: EffectContainer | None = None,
        uLnTx: _HasTagAndGet[_ConvertibleToBool] | _ConvertibleToBool = None,
        uLn: LineProperties | None = None,
        uFillTx: _HasTagAndGet[_ConvertibleToBool] | _ConvertibleToBool = None,
        uFill: _HasTagAndGet[_ConvertibleToBool] | _ConvertibleToBool = None,
    ) -> None: ...

class TabStop(Serialisable):
    pos: Typed[Coordinate[bool], Literal[True]]
    algn: Typed[Set[_TabStopAlgn], Literal[False]]
    def __init__(self, pos: Coordinate[bool] | None = None, algn: Set[_TabStopAlgn] | None = None) -> None: ...

class TabStopList(Serialisable):
    tab: Typed[TabStop, Literal[True]]
    def __init__(self, tab=None) -> None: ...

class Spacing(Serialisable):
    spcPct: NestedInteger[Literal[True]]
    spcPts: NestedInteger[Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        spcPct: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        spcPts: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
    ) -> None: ...

class AutonumberBullet(Serialisable):
    type: Set[_AutonumberBulletType]
    startAt: Integer[Literal[False]]
    def __init__(self, type: _AutonumberBulletType, startAt: ConvertibleToInt) -> None: ...

class ParagraphProperties(Serialisable):
    tagname: ClassVar[str]
    namespace: ClassVar[str]
    marL: Integer[Literal[True]]
    marR: Integer[Literal[True]]
    lvl: Integer[Literal[True]]
    indent: Integer[Literal[True]]
    algn: NoneSet[_ParagraphPropertiesAlgn]
    defTabSz: Integer[Literal[True]]
    rtl: Bool[Literal[True]]
    eaLnBrk: Bool[Literal[True]]
    fontAlgn: NoneSet[_ParagraphPropertiesFontAlgn]
    latinLnBrk: Bool[Literal[True]]
    hangingPunct: Bool[Literal[True]]
    lnSpc: Typed[Spacing, Literal[True]]
    spcBef: Typed[Spacing, Literal[True]]
    spcAft: Typed[Spacing, Literal[True]]
    tabLst: Typed[TabStopList, Literal[True]]
    defRPr: Typed[CharacterProperties, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    buClrTx: EmptyTag[Literal[False]]
    buClr: Typed[Color, Literal[True]]
    buSzTx: EmptyTag[Literal[False]]
    buSzPct: NestedInteger[Literal[True]]
    buSzPts: NestedInteger[Literal[True]]
    buFontTx: EmptyTag[Literal[False]]
    buFont: Typed[Font, Literal[True]]
    buNone: EmptyTag[Literal[False]]
    buAutoNum: EmptyTag[Literal[False]]
    buChar: NestedValue[str, Literal[True]]
    buBlip: NestedValue[Blip, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        marL: ConvertibleToInt | None = None,
        marR: ConvertibleToInt | None = None,
        lvl: ConvertibleToInt | None = None,
        indent: ConvertibleToInt | None = None,
        algn: _ParagraphPropertiesAlgn | Literal["none"] | None = None,
        defTabSz: ConvertibleToInt | None = None,
        rtl: _ConvertibleToBool | None = None,
        eaLnBrk: _ConvertibleToBool | None = None,
        fontAlgn: _ParagraphPropertiesFontAlgn | Literal["none"] | None = None,
        latinLnBrk: _ConvertibleToBool | None = None,
        hangingPunct: _ConvertibleToBool | None = None,
        lnSpc: Spacing | None = None,
        spcBef: Spacing | None = None,
        spcAft: Spacing | None = None,
        tabLst: TabStopList | None = None,
        defRPr: CharacterProperties | None = None,
        extLst: ExtensionList | None = None,
        buClrTx: _HasTagAndGet[_ConvertibleToBool] | _ConvertibleToBool = None,
        buClr: Color | None = None,
        buSzTx: _HasTagAndGet[_ConvertibleToBool] | _ConvertibleToBool = None,
        buSzPct: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        buSzPts: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        buFontTx: _HasTagAndGet[_ConvertibleToBool] | _ConvertibleToBool = None,
        buFont: Font | None = None,
        buNone: _HasTagAndGet[_ConvertibleToBool] | _ConvertibleToBool = None,
        buAutoNum: _HasTagAndGet[_ConvertibleToBool] | _ConvertibleToBool = None,
        buChar: object = None,
        buBlip: object = None,
    ) -> None: ...

class ListStyle(Serialisable):
    tagname: ClassVar[str]
    namespace: ClassVar[str]
    defPPr: Typed[ParagraphProperties, Literal[True]]
    lvl1pPr: Typed[ParagraphProperties, Literal[True]]
    lvl2pPr: Typed[ParagraphProperties, Literal[True]]
    lvl3pPr: Typed[ParagraphProperties, Literal[True]]
    lvl4pPr: Typed[ParagraphProperties, Literal[True]]
    lvl5pPr: Typed[ParagraphProperties, Literal[True]]
    lvl6pPr: Typed[ParagraphProperties, Literal[True]]
    lvl7pPr: Typed[ParagraphProperties, Literal[True]]
    lvl8pPr: Typed[ParagraphProperties, Literal[True]]
    lvl9pPr: Typed[ParagraphProperties, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        defPPr: ParagraphProperties | None = None,
        lvl1pPr: ParagraphProperties | None = None,
        lvl2pPr: ParagraphProperties | None = None,
        lvl3pPr: ParagraphProperties | None = None,
        lvl4pPr: ParagraphProperties | None = None,
        lvl5pPr: ParagraphProperties | None = None,
        lvl6pPr: ParagraphProperties | None = None,
        lvl7pPr: ParagraphProperties | None = None,
        lvl8pPr: ParagraphProperties | None = None,
        lvl9pPr: ParagraphProperties | None = None,
        extLst: ParagraphProperties | None = None,
    ) -> None: ...

class RegularTextRun(Serialisable):
    tagname: ClassVar[str]
    namespace: ClassVar[str]
    rPr: Typed[CharacterProperties, Literal[True]]
    properties: Alias
    t: NestedText[str, Literal[False]]
    value: Alias
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, rPr: CharacterProperties | None = None, t: object = "") -> None: ...

class LineBreak(Serialisable):
    tagname: ClassVar[str]
    namespace: ClassVar[str]
    rPr: Typed[CharacterProperties, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, rPr: CharacterProperties | None = None) -> None: ...

class TextField(Serialisable):
    id: String[Literal[False]]
    type: String[Literal[True]]
    rPr: Typed[CharacterProperties, Literal[True]]
    pPr: Typed[CharacterProperties, Literal[True]]
    t: String[Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        id: str,
        type: str | None = None,
        rPr: CharacterProperties | None = None,
        pPr: CharacterProperties | None = None,
        t: str | None = None,
    ) -> None: ...

class Paragraph(Serialisable):
    tagname: ClassVar[str]
    namespace: ClassVar[str]
    pPr: Typed[ParagraphProperties, Literal[True]]
    properties: Alias
    endParaRPr: Typed[CharacterProperties, Literal[True]]
    r: Incomplete
    text: Alias
    br: Typed[LineBreak, Literal[True]]
    fld: Typed[TextField, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        pPr: ParagraphProperties | None = None,
        endParaRPr: CharacterProperties | None = None,
        r=None,
        br: LineBreak | None = None,
        fld: TextField | None = None,
    ) -> None: ...

class GeomGuide(Serialisable):
    name: String[Literal[False]]
    fmla: String[Literal[False]]
    def __init__(self, name: str, fmla: str) -> None: ...

class GeomGuideList(Serialisable):
    gd: Incomplete
    def __init__(self, gd=None) -> None: ...

class PresetTextShape(Serialisable):
    prst: Typed[Set[_PresetTextShapePrst], Literal[False]]
    avLst: Typed[GeomGuideList, Literal[True]]
    def __init__(self, prst: Set[_PresetTextShapePrst], avLst: GeomGuideList | None = None) -> None: ...

class TextNormalAutofit(Serialisable):
    fontScale: Integer[Literal[False]]
    lnSpcReduction: Integer[Literal[False]]
    def __init__(self, fontScale: ConvertibleToInt, lnSpcReduction: ConvertibleToInt) -> None: ...

class RichTextProperties(Serialisable):
    tagname: ClassVar[str]
    namespace: ClassVar[str]
    rot: Integer[Literal[True]]
    spcFirstLastPara: Bool[Literal[True]]
    vertOverflow: NoneSet[_RichTextPropertiesVertOverflow]
    horzOverflow: NoneSet[_RichTextPropertiesHorzOverflow]
    vert: NoneSet[_RichTextPropertiesVert]
    wrap: NoneSet[_RichTextPropertiesWrap]
    lIns: Integer[Literal[True]]
    tIns: Integer[Literal[True]]
    rIns: Integer[Literal[True]]
    bIns: Integer[Literal[True]]
    numCol: Integer[Literal[True]]
    spcCol: Integer[Literal[True]]
    rtlCol: Bool[Literal[True]]
    fromWordArt: Bool[Literal[True]]
    anchor: NoneSet[_RichTextPropertiesAnchor]
    anchorCtr: Bool[Literal[True]]
    forceAA: Bool[Literal[True]]
    upright: Bool[Literal[True]]
    compatLnSpc: Bool[Literal[True]]
    prstTxWarp: Typed[PresetTextShape, Literal[True]]
    scene3d: Typed[Scene3D, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    noAutofit: EmptyTag[Literal[False]]
    normAutofit: EmptyTag[Literal[False]]
    spAutoFit: EmptyTag[Literal[False]]
    flatTx: NestedInteger[Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        rot: ConvertibleToInt | None = None,
        spcFirstLastPara: _ConvertibleToBool | None = None,
        vertOverflow: _RichTextPropertiesVertOverflow | Literal["none"] | None = None,
        horzOverflow: _RichTextPropertiesHorzOverflow | Literal["none"] | None = None,
        vert: _RichTextPropertiesVert | Literal["none"] | None = None,
        wrap: _RichTextPropertiesWrap | Literal["none"] | None = None,
        lIns: ConvertibleToInt | None = None,
        tIns: ConvertibleToInt | None = None,
        rIns: ConvertibleToInt | None = None,
        bIns: ConvertibleToInt | None = None,
        numCol: ConvertibleToInt | None = None,
        spcCol: ConvertibleToInt | None = None,
        rtlCol: _ConvertibleToBool | None = None,
        fromWordArt: _ConvertibleToBool | None = None,
        anchor: _RichTextPropertiesAnchor | Literal["none"] | None = None,
        anchorCtr: _ConvertibleToBool | None = None,
        forceAA: _ConvertibleToBool | None = None,
        upright: _ConvertibleToBool | None = None,
        compatLnSpc: _ConvertibleToBool | None = None,
        prstTxWarp: PresetTextShape | None = None,
        scene3d: Scene3D | None = None,
        extLst: Unused = None,
        noAutofit: _HasTagAndGet[_ConvertibleToBool] | _ConvertibleToBool = None,
        normAutofit: _HasTagAndGet[_ConvertibleToBool] | _ConvertibleToBool = None,
        spAutoFit: _HasTagAndGet[_ConvertibleToBool] | _ConvertibleToBool = None,
        flatTx: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
    ) -> None: ...
