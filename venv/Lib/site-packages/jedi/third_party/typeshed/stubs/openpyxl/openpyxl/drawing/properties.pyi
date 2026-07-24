from _typeshed import Incomplete, Unused
from typing import ClassVar, Literal, overload
from typing_extensions import TypeAlias

from openpyxl.descriptors.base import Bool, NoneSet, String, Typed, _ConvertibleToBool
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.drawing.geometry import GroupTransform2D, Scene3D
from openpyxl.drawing.text import Hyperlink

_GroupShapePropertiesBwMode: TypeAlias = Literal[
    "clr", "auto", "gray", "ltGray", "invGray", "grayWhite", "blackGray", "blackWhite", "black", "white", "hidden"
]

class GroupShapeProperties(Serialisable):
    tagname: ClassVar[str]
    bwMode: NoneSet[_GroupShapePropertiesBwMode]
    xfrm: Typed[GroupTransform2D, Literal[True]]
    scene3d: Typed[Scene3D, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    def __init__(
        self,
        bwMode: _GroupShapePropertiesBwMode | Literal["none"] | None = None,
        xfrm: GroupTransform2D | None = None,
        scene3d: Scene3D | None = None,
        extLst: ExtensionList | None = None,
    ) -> None: ...

class GroupLocking(Serialisable):
    tagname: ClassVar[str]
    namespace: ClassVar[str]
    noGrp: Bool[Literal[True]]
    noUngrp: Bool[Literal[True]]
    noSelect: Bool[Literal[True]]
    noRot: Bool[Literal[True]]
    noChangeAspect: Bool[Literal[True]]
    noMove: Bool[Literal[True]]
    noResize: Bool[Literal[True]]
    noChangeArrowheads: Bool[Literal[True]]
    noEditPoints: Bool[Literal[True]]
    noAdjustHandles: Bool[Literal[True]]
    noChangeShapeType: Bool[Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        noGrp: _ConvertibleToBool | None = None,
        noUngrp: _ConvertibleToBool | None = None,
        noSelect: _ConvertibleToBool | None = None,
        noRot: _ConvertibleToBool | None = None,
        noChangeAspect: _ConvertibleToBool | None = None,
        noChangeArrowheads: _ConvertibleToBool | None = None,
        noMove: _ConvertibleToBool | None = None,
        noResize: _ConvertibleToBool | None = None,
        noEditPoints: _ConvertibleToBool | None = None,
        noAdjustHandles: _ConvertibleToBool | None = None,
        noChangeShapeType: _ConvertibleToBool | None = None,
        extLst: Unused = None,
    ) -> None: ...

class NonVisualGroupDrawingShapeProps(Serialisable):
    tagname: ClassVar[str]
    grpSpLocks: Typed[GroupLocking, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, grpSpLocks=None, extLst: Unused = None) -> None: ...

class NonVisualDrawingShapeProps(Serialisable):
    tagname: ClassVar[str]
    spLocks: Typed[GroupLocking, Literal[True]]
    txBax: Bool[Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    txBox: Incomplete
    def __init__(self, spLocks=None, txBox: _ConvertibleToBool | None = None, extLst: Unused = None) -> None: ...

class NonVisualDrawingProps(Serialisable):
    tagname: ClassVar[str]
    id: Incomplete
    name: String[Literal[False]]
    descr: String[Literal[True]]
    hidden: Bool[Literal[True]]
    title: String[Literal[True]]
    hlinkClick: Typed[Hyperlink, Literal[True]]
    hlinkHover: Typed[Hyperlink, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    # Source incorrectly uses a list here instead of a tuple
    __elements__: ClassVar[list[str]]  # type: ignore[assignment]
    @overload
    def __init__(
        self,
        id=None,
        *,
        name: str,
        descr: str | None = None,
        hidden: _ConvertibleToBool | None = None,
        title: str | None = None,
        hlinkClick: Hyperlink | None = None,
        hlinkHover: Hyperlink | None = None,
        extLst: ExtensionList | None = None,
    ) -> None: ...
    @overload
    def __init__(
        self,
        id: Incomplete | None,
        name: str,
        descr: str | None = None,
        hidden: _ConvertibleToBool | None = None,
        title: str | None = None,
        hlinkClick: Hyperlink | None = None,
        hlinkHover: Hyperlink | None = None,
        extLst: ExtensionList | None = None,
    ) -> None: ...

class NonVisualGroupShape(Serialisable):
    tagname: ClassVar[str]
    cNvPr: Typed[NonVisualDrawingProps, Literal[False]]
    cNvGrpSpPr: Typed[NonVisualGroupDrawingShapeProps, Literal[False]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, cNvPr: NonVisualDrawingProps, cNvGrpSpPr: NonVisualGroupDrawingShapeProps) -> None: ...
