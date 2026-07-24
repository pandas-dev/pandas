from _typeshed import Incomplete
from weakref import ReferenceType

from reportlab.pdfbase.pdfdoc import PDFDictionary, PDFObject, PDFStream

__all__ = ("AcroForm",)

visibilities: dict[str, int]
orientations: dict[str, list[Incomplete]]
fieldFlagValues: dict[str, int]
annotationFlagValues: dict[str, int]

def bsPDF(borderWidth: int, borderStyle: str, dashLen) -> PDFDictionary: ...
def escPDF(s) -> str: ...
def makeFlags(s: int | str, d: dict[str, int] = ...) -> int: ...

class PDFFromString(PDFObject):
    def __init__(self, s: str | bytes) -> None: ...
    def format(self, document) -> bytes: ...

class RadioGroup(PDFObject):
    TU: Incomplete
    Ff: int
    kids: list[Incomplete]
    T: Incomplete
    V: Incomplete
    def __init__(self, name, tooltip: str = "", fieldFlags: str = "noToggleToOff required radio") -> None: ...
    def format(self, doc) -> bytes: ...

class AcroForm(PDFObject):
    formFontNames: dict[str, str]
    referenceMap: dict[Incomplete, Incomplete]
    fonts: dict[str, str]
    fields: list[Incomplete]
    sigFlags: Incomplete
    extras: dict[Incomplete, Incomplete]
    def __init__(self, canv, **kwds) -> None: ...
    @property
    def canv(self) -> ReferenceType[Incomplete]: ...
    def fontRef(self, f) -> str: ...
    def format(self, doc) -> bytes: ...
    def colorTuple(self, c): ...
    def streamFillColor(self, c) -> str: ...
    def streamStrokeColor(self, c) -> str: ...
    def checkboxAP(
        self,
        key,
        value,
        buttonStyle: str = "circle",
        shape: str = "square",
        fillColor=None,
        borderColor=None,
        textColor=None,
        borderWidth: int = 1,
        borderStyle: str = "solid",
        size: int = 20,
        dashLen: int = 3,
    ) -> PDFStream: ...
    @staticmethod
    def circleArcStream(size, r, arcs=(0, 1, 2, 3), rotated: bool = False) -> str: ...
    def zdMark(self, c, size, ds, iFontName) -> str: ...
    def getRef(self, obj): ...
    def getRefStr(self, obj) -> str: ...
    @staticmethod
    def stdColors(t, b, f) -> tuple[Incomplete, Incomplete, Incomplete]: ...
    @staticmethod
    def varyColors(key, t, b, f) -> tuple[Incomplete, Incomplete, Incomplete]: ...
    def checkForceBorder(
        self, x, y, width, height, forceBorder, shape, borderStyle, borderWidth, borderColor, fillColor
    ) -> None: ...
    def checkbox(
        self,
        checked: bool = False,
        buttonStyle: str = "check",
        shape: str = "square",
        fillColor=None,
        borderColor=None,
        textColor=None,
        borderWidth: int = 1,
        borderStyle: str = "solid",
        size: int = 20,
        x: int = 0,
        y: int = 0,
        tooltip=None,
        name=None,
        annotationFlags: str = "print",
        fieldFlags: str = "required",
        forceBorder: bool = False,
        relative: bool = False,
        dashLen: int = 3,
    ) -> None: ...
    def radio(
        self,
        value=None,
        selected: bool = False,
        buttonStyle: str = "circle",
        shape: str = "circle",
        fillColor=None,
        borderColor=None,
        textColor=None,
        borderWidth: int = 1,
        borderStyle: str = "solid",
        size: int = 20,
        x: int = 0,
        y: int = 0,
        tooltip=None,
        name=None,
        annotationFlags: str = "print",
        fieldFlags: str = "noToggleToOff required radio",
        forceBorder: bool = False,
        relative: bool = False,
        dashLen: int = 3,
    ) -> None: ...
    def makeStream(self, width, height, stream, **D) -> PDFStream: ...
    def txAP(
        self,
        key,
        value,
        iFontName,
        rFontName,
        fontSize,
        shape: str = "square",
        fillColor=None,
        borderColor=None,
        textColor=None,
        borderWidth: int = 1,
        borderStyle: str = "solid",
        width: int = 120,
        height: int = 36,
        dashLen: int = 3,
        wkind: str = "textfield",
        labels=[],
        I=[],
        sel_bg: str = "0.600006 0.756866 0.854904 rg",
        sel_fg: str = "0 g",
    ) -> PDFStream: ...
    def makeFont(self, fontName: str | None) -> tuple[str, str]: ...
    def textfield(
        self,
        value: str = "",
        fillColor=None,
        borderColor=None,
        textColor=None,
        borderWidth: int = 1,
        borderStyle: str = "solid",
        width: int = 120,
        height: int = 36,
        x: int = 0,
        y: int = 0,
        tooltip=None,
        name=None,
        annotationFlags: str = "print",
        fieldFlags: str = "",
        forceBorder: bool = False,
        relative: bool = False,
        maxlen: int = 100,
        fontName: str | None = None,
        fontSize=None,
        dashLen: int = 3,
    ) -> None: ...
    def listbox(
        self,
        value: str = "",
        fillColor=None,
        borderColor=None,
        textColor=None,
        borderWidth: int = 1,
        borderStyle: str = "solid",
        width: int = 120,
        height: int = 36,
        x: int = 0,
        y: int = 0,
        tooltip=None,
        name=None,
        annotationFlags: str = "print",
        fieldFlags: str = "",
        forceBorder: bool = False,
        relative: bool = False,
        fontName: str | None = None,
        fontSize=None,
        dashLen: int = 3,
        maxlen=None,
        options=[],
    ) -> None: ...
    def choice(
        self,
        value: str = "",
        fillColor=None,
        borderColor=None,
        textColor=None,
        borderWidth: int = 1,
        borderStyle: str = "solid",
        width: int = 120,
        height: int = 36,
        x: int = 0,
        y: int = 0,
        tooltip=None,
        name=None,
        annotationFlags: str = "print",
        fieldFlags: str = "combo",
        forceBorder: bool = False,
        relative: bool = False,
        fontName: str | None = None,
        fontSize=None,
        dashLen: int = 3,
        maxlen=None,
        options=[],
    ) -> None: ...
    def checkboxRelative(
        self,
        *,
        checked: bool = False,
        buttonStyle: str = "check",
        shape: str = "square",
        fillColor=None,
        borderColor=None,
        textColor=None,
        borderWidth: int = 1,
        borderStyle: str = "solid",
        size: int = 20,
        x: int = 0,
        y: int = 0,
        tooltip=None,
        name=None,
        annotationFlags: str = "print",
        fieldFlags: str = "required",
        forceBorder: bool = False,
        dashLen: int = 3,
    ) -> None: ...
    def radioRelative(
        self,
        *,
        value=None,
        selected: bool = False,
        buttonStyle: str = "circle",
        shape: str = "circle",
        fillColor=None,
        borderColor=None,
        textColor=None,
        borderWidth: int = 1,
        borderStyle: str = "solid",
        size: int = 20,
        x: int = 0,
        y: int = 0,
        tooltip=None,
        name=None,
        annotationFlags: str = "print",
        fieldFlags: str = "noToggleToOff required radio",
        forceBorder: bool = False,
        dashLen: int = 3,
    ) -> None: ...
    def textfieldRelative(
        self,
        *,
        value: str = "",
        fillColor=None,
        borderColor=None,
        textColor=None,
        borderWidth: int = 1,
        borderStyle: str = "solid",
        width: int = 120,
        height: int = 36,
        x: int = 0,
        y: int = 0,
        tooltip=None,
        name=None,
        annotationFlags: str = "print",
        fieldFlags: str = "",
        forceBorder: bool = False,
        maxlen: int = 100,
        fontName: str | None = None,
        fontSize=None,
        dashLen: int = 3,
    ) -> None: ...
    def listboxRelative(
        self,
        *,
        value: str = "",
        fillColor=None,
        borderColor=None,
        textColor=None,
        borderWidth: int = 1,
        borderStyle: str = "solid",
        width: int = 120,
        height: int = 36,
        x: int = 0,
        y: int = 0,
        tooltip=None,
        name=None,
        annotationFlags: str = "print",
        fieldFlags: str = "",
        forceBorder: bool = False,
        maxlen: int = 100,
        fontName: str | None = None,
        fontSize=None,
        dashLen: int = 3,
    ) -> None: ...
    def choiceRelative(
        self,
        *,
        value: str = "",
        fillColor=None,
        borderColor=None,
        textColor=None,
        borderWidth: int = 1,
        borderStyle: str = "solid",
        width: int = 120,
        height: int = 36,
        x: int = 0,
        y: int = 0,
        tooltip=None,
        name=None,
        annotationFlags: str = "print",
        fieldFlags: str = "",
        forceBorder: bool = False,
        maxlen: int = 100,
        fontName: str | None = None,
        fontSize=None,
        dashLen: int = 3,
    ) -> None: ...
    @property
    def encRefStr(self) -> str: ...

class CBMark:
    opNames: list[str]
    opCount: tuple[int, ...]
    ops: Incomplete
    xmin: Incomplete
    ymin: Incomplete
    xmax: Incomplete
    ymax: Incomplete
    points: Incomplete
    slack: Incomplete
    def __init__(self, ops, points, bounds, slack: float = 0.05) -> None: ...
    def scaledRender(self, size, ds: int = 0) -> str: ...
