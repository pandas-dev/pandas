from _typeshed import ConvertibleToInt
from typing import ClassVar, Literal, overload

from openpyxl.chart.shapes import GraphicalProperties
from openpyxl.chart.text import RichText
from openpyxl.descriptors import Typed
from openpyxl.descriptors.base import Alias, Bool, Integer, String, _ConvertibleToBool
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.drawing.geometry import ShapeStyle
from openpyxl.drawing.properties import NonVisualDrawingProps, NonVisualDrawingShapeProps

class Connection(Serialisable):
    id: Integer[Literal[False]]
    idx: Integer[Literal[False]]
    def __init__(self, id: ConvertibleToInt, idx: ConvertibleToInt) -> None: ...

class ConnectorLocking(Serialisable):
    extLst: Typed[ExtensionList, Literal[True]]
    def __init__(self, extLst: ExtensionList | None = None) -> None: ...

class NonVisualConnectorProperties(Serialisable):
    cxnSpLocks: Typed[ConnectorLocking, Literal[True]]
    stCxn: Typed[Connection, Literal[True]]
    endCxn: Typed[Connection, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    def __init__(
        self,
        cxnSpLocks: ConnectorLocking | None = None,
        stCxn: Connection | None = None,
        endCxn: Connection | None = None,
        extLst: ExtensionList | None = None,
    ) -> None: ...

class ConnectorNonVisual(Serialisable):
    cNvPr: Typed[NonVisualDrawingProps, Literal[False]]
    cNvCxnSpPr: Typed[NonVisualConnectorProperties, Literal[False]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, cNvPr: NonVisualDrawingProps, cNvCxnSpPr: NonVisualConnectorProperties) -> None: ...

class ConnectorShape(Serialisable):
    tagname: ClassVar[str]
    nvCxnSpPr: Typed[ConnectorNonVisual, Literal[False]]
    spPr: Typed[GraphicalProperties, Literal[False]]
    style: Typed[ShapeStyle, Literal[True]]
    macro: String[Literal[True]]
    fPublished: Bool[Literal[True]]
    def __init__(
        self,
        nvCxnSpPr: ConnectorNonVisual,
        spPr: GraphicalProperties,
        style: ShapeStyle | None = None,
        macro: str | None = None,
        fPublished: _ConvertibleToBool | None = None,
    ) -> None: ...

class ShapeMeta(Serialisable):
    tagname: ClassVar[str]
    cNvPr: Typed[NonVisualDrawingProps, Literal[False]]
    cNvSpPr: Typed[NonVisualDrawingShapeProps, Literal[False]]
    def __init__(self, cNvPr: NonVisualDrawingProps, cNvSpPr: NonVisualDrawingShapeProps) -> None: ...

class Shape(Serialisable):
    macro: String[Literal[True]]
    textlink: String[Literal[True]]
    fPublished: Bool[Literal[True]]
    fLocksText: Bool[Literal[True]]
    nvSpPr: Typed[ShapeMeta, Literal[True]]
    meta: Alias
    spPr: Typed[GraphicalProperties, Literal[False]]
    graphicalProperties: Alias
    style: Typed[ShapeStyle, Literal[True]]
    txBody: Typed[RichText, Literal[True]]
    @overload
    def __init__(
        self,
        macro: str | None = None,
        textlink: str | None = None,
        fPublished: _ConvertibleToBool | None = None,
        fLocksText: _ConvertibleToBool | None = None,
        nvSpPr: ShapeMeta | None = None,
        *,
        spPr: GraphicalProperties,
        style: ShapeStyle | None = None,
        txBody: RichText | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        macro: str | None,
        textlink: str | None,
        fPublished: _ConvertibleToBool | None,
        fLocksText: _ConvertibleToBool | None,
        nvSpPr: ShapeMeta | None,
        spPr: GraphicalProperties,
        style: ShapeStyle | None = None,
        txBody: RichText | None = None,
    ) -> None: ...
