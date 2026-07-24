from _typeshed import Unused
from typing import ClassVar, Literal

from openpyxl.chart.shapes import GraphicalProperties
from openpyxl.descriptors.base import Alias, Bool, String, Typed, _ConvertibleToBool
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.drawing.fill import BlipFillProperties
from openpyxl.drawing.geometry import ShapeStyle
from openpyxl.drawing.properties import NonVisualDrawingProps

class PictureLocking(Serialisable):
    tagname: ClassVar[str]
    namespace: ClassVar[str]
    noCrop: Bool[Literal[True]]
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
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        noCrop: _ConvertibleToBool | None = None,
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
        extLst: Unused = None,
    ) -> None: ...

class NonVisualPictureProperties(Serialisable):
    tagname: ClassVar[str]
    preferRelativeResize: Bool[Literal[True]]
    picLocks: Typed[PictureLocking, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, preferRelativeResize: _ConvertibleToBool | None = None, picLocks=None, extLst: Unused = None) -> None: ...

class PictureNonVisual(Serialisable):
    tagname: ClassVar[str]
    cNvPr: Typed[NonVisualDrawingProps, Literal[False]]
    cNvPicPr: Typed[NonVisualPictureProperties, Literal[False]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self, cNvPr: NonVisualDrawingProps | None = None, cNvPicPr: NonVisualPictureProperties | None = None
    ) -> None: ...

class PictureFrame(Serialisable):
    tagname: ClassVar[str]
    macro: String[Literal[True]]
    fPublished: Bool[Literal[True]]
    nvPicPr: Typed[PictureNonVisual, Literal[False]]
    blipFill: Typed[BlipFillProperties, Literal[False]]
    spPr: Typed[GraphicalProperties, Literal[False]]
    graphicalProperties: Alias
    style: Typed[ShapeStyle, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        macro: str | None = None,
        fPublished: _ConvertibleToBool | None = None,
        nvPicPr: PictureNonVisual | None = None,
        blipFill: BlipFillProperties | None = None,
        spPr: GraphicalProperties | None = None,
        style: ShapeStyle | None = None,
    ) -> None: ...
