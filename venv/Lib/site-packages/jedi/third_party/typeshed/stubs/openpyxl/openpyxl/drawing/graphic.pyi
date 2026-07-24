from typing import ClassVar, Literal

from openpyxl.descriptors.base import Alias, Bool, String, Typed, _ConvertibleToBool
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.drawing.picture import PictureFrame
from openpyxl.drawing.properties import GroupShapeProperties, NonVisualGroupShape
from openpyxl.drawing.relation import ChartRelation
from openpyxl.drawing.xdr import XDRTransform2D

class GraphicFrameLocking(Serialisable):
    noGrp: Bool[Literal[True]]
    noDrilldown: Bool[Literal[True]]
    noSelect: Bool[Literal[True]]
    noChangeAspect: Bool[Literal[True]]
    noMove: Bool[Literal[True]]
    noResize: Bool[Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    def __init__(
        self,
        noGrp: _ConvertibleToBool | None = None,
        noDrilldown: _ConvertibleToBool | None = None,
        noSelect: _ConvertibleToBool | None = None,
        noChangeAspect: _ConvertibleToBool | None = None,
        noMove: _ConvertibleToBool | None = None,
        noResize: _ConvertibleToBool | None = None,
        extLst: ExtensionList | None = None,
    ) -> None: ...

class NonVisualGraphicFrameProperties(Serialisable):
    tagname: ClassVar[str]
    graphicFrameLocks: Typed[GraphicFrameLocking, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    def __init__(self, graphicFrameLocks: GraphicFrameLocking | None = None, extLst: ExtensionList | None = None) -> None: ...

class NonVisualGraphicFrame(Serialisable):
    tagname: ClassVar[str]
    cNvPr: Typed[ExtensionList, Literal[False]]
    cNvGraphicFramePr: Typed[ExtensionList, Literal[False]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, cNvPr=None, cNvGraphicFramePr=None) -> None: ...

class GraphicData(Serialisable):
    tagname: ClassVar[str]
    namespace: ClassVar[str]
    uri: String[Literal[False]]
    chart: Typed[ChartRelation, Literal[True]]
    def __init__(self, uri: str = ..., chart: ChartRelation | None = None) -> None: ...

class GraphicObject(Serialisable):
    tagname: ClassVar[str]
    namespace: ClassVar[str]
    graphicData: Typed[GraphicData, Literal[False]]
    def __init__(self, graphicData: GraphicData | None = None) -> None: ...

class GraphicFrame(Serialisable):
    tagname: ClassVar[str]
    nvGraphicFramePr: Typed[NonVisualGraphicFrame, Literal[False]]
    xfrm: Typed[XDRTransform2D, Literal[False]]
    graphic: Typed[GraphicObject, Literal[False]]
    macro: String[Literal[True]]
    fPublished: Bool[Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        nvGraphicFramePr: NonVisualGraphicFrame | None = None,
        xfrm: XDRTransform2D | None = None,
        graphic: GraphicObject | None = None,
        macro: str | None = None,
        fPublished: _ConvertibleToBool | None = None,
    ) -> None: ...

class GroupShape(Serialisable):
    nvGrpSpPr: Typed[NonVisualGroupShape, Literal[False]]
    nonVisualProperties: Alias
    grpSpPr: Typed[GroupShapeProperties, Literal[False]]
    visualProperties: Alias
    pic: Typed[PictureFrame, Literal[True]]
    # Source incorrectly uses a list here instead of a tuple
    __elements__: ClassVar[list[str]]  # type: ignore[assignment]
    def __init__(
        self, nvGrpSpPr: NonVisualGroupShape, grpSpPr: GroupShapeProperties, pic: PictureFrame | None = None
    ) -> None: ...
