from _typeshed import Incomplete, Unused
from typing import ClassVar, Literal
from typing_extensions import TypeAlias

from openpyxl.descriptors.base import Alias, NoneSet, Typed, _ConvertibleToBool
from openpyxl.descriptors.nested import EmptyTag
from openpyxl.descriptors.serialisable import Serialisable
from openpyxl.drawing.colors import ColorChoice, ColorChoiceDescriptor
from openpyxl.drawing.fill import GradientFillProperties, PatternFillProperties
from openpyxl.drawing.geometry import CustomGeometry2D, PresetGeometry2D, Scene3D, Shape3D, Transform2D
from openpyxl.drawing.line import LineProperties

from ..xml._functions_overloads import _HasTagAndGet

_GraphicalPropertiesBwMode: TypeAlias = Literal[
    "clr", "auto", "gray", "ltGray", "invGray", "grayWhite", "blackGray", "blackWhite", "black", "white", "hidden"
]

class GraphicalProperties(Serialisable):
    tagname: ClassVar[str]
    bwMode: NoneSet[_GraphicalPropertiesBwMode]
    xfrm: Typed[Transform2D, Literal[True]]
    transform: Alias
    custGeom: Typed[CustomGeometry2D, Literal[True]]
    prstGeom: Typed[PresetGeometry2D, Literal[True]]
    noFill: EmptyTag[Literal[False]]
    solidFill: ColorChoiceDescriptor
    gradFill: Typed[GradientFillProperties, Literal[True]]
    pattFill: Typed[PatternFillProperties, Literal[True]]
    ln: Typed[LineProperties, Literal[True]]
    line: Alias
    scene3d: Typed[Scene3D, Literal[True]]
    sp3d: Typed[Shape3D, Literal[True]]
    shape3D: Alias
    extLst: Typed[Incomplete, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        bwMode: _GraphicalPropertiesBwMode | Literal["none"] | None = None,
        xfrm: Transform2D | None = None,
        noFill: _HasTagAndGet[_ConvertibleToBool] | _ConvertibleToBool = None,
        solidFill: str | ColorChoice | None = None,
        gradFill: GradientFillProperties | None = None,
        pattFill: PatternFillProperties | None = None,
        ln=None,
        scene3d: Scene3D | None = None,
        custGeom: CustomGeometry2D | None = None,
        prstGeom: PresetGeometry2D | None = None,
        sp3d: Shape3D | None = None,
        extLst: Unused = None,
    ) -> None: ...
