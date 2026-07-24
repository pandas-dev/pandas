from _typeshed import ConvertibleToFloat, ConvertibleToInt, Unused
from typing import ClassVar, Literal

from openpyxl.chart.picture import PictureOptions
from openpyxl.chart.shapes import GraphicalProperties
from openpyxl.descriptors.base import Alias, Typed, _ConvertibleToBool
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.nested import NestedBool, NestedInteger, NestedMinMax
from openpyxl.descriptors.serialisable import Serialisable

from ..xml._functions_overloads import _HasTagAndGet

class View3D(Serialisable):
    tagname: ClassVar[str]
    rotX: NestedMinMax[float, Literal[True]]
    x_rotation: Alias
    hPercent: NestedMinMax[float, Literal[True]]
    height_percent: Alias
    rotY: NestedInteger[Literal[True]]
    y_rotation: Alias
    depthPercent: NestedInteger[Literal[True]]
    rAngAx: NestedBool[Literal[True]]
    right_angle_axes: Alias
    perspective: NestedInteger[Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        rotX: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None = 15,
        hPercent: _HasTagAndGet[ConvertibleToFloat | None] | ConvertibleToFloat | None = None,
        rotY: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = 20,
        depthPercent: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        rAngAx: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = True,
        perspective: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        extLst: Unused = None,
    ) -> None: ...

class Surface(Serialisable):
    tagname: ClassVar[str]
    thickness: NestedInteger[Literal[True]]
    spPr: Typed[GraphicalProperties, Literal[True]]
    graphicalProperties: Alias
    pictureOptions: Typed[PictureOptions, Literal[True]]
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        thickness: _HasTagAndGet[ConvertibleToInt | None] | ConvertibleToInt | None = None,
        spPr: GraphicalProperties | None = None,
        pictureOptions: PictureOptions | None = None,
        extLst: Unused = None,
    ) -> None: ...

class _3DBase(Serialisable):
    tagname: ClassVar[str]
    view3D: Typed[View3D, Literal[True]]
    floor: Typed[Surface, Literal[True]]
    sideWall: Typed[Surface, Literal[True]]
    backWall: Typed[Surface, Literal[True]]
    def __init__(
        self,
        view3D: View3D | None = None,
        floor: Surface | None = None,
        sideWall: Surface | None = None,
        backWall: Surface | None = None,
    ) -> None: ...
