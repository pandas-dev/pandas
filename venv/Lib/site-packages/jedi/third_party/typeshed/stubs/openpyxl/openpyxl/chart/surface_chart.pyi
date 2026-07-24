from _typeshed import ConvertibleToInt, Incomplete
from typing import ClassVar, Literal

from openpyxl.chart.axis import NumericAxis, SeriesAxis, TextAxis
from openpyxl.chart.shapes import GraphicalProperties
from openpyxl.descriptors.base import Alias, Typed, _ConvertibleToBool
from openpyxl.descriptors.excel import ExtensionList
from openpyxl.descriptors.nested import NestedBool, NestedInteger
from openpyxl.descriptors.serialisable import Serialisable

from ..xml._functions_overloads import _HasTagAndGet
from ._3d import _3DBase
from ._chart import ChartBase

class BandFormat(Serialisable):
    tagname: ClassVar[str]
    idx: NestedInteger[Literal[False]]
    spPr: Typed[GraphicalProperties, Literal[True]]
    graphicalProperties: Alias
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self, idx: _HasTagAndGet[ConvertibleToInt] | ConvertibleToInt = 0, spPr: GraphicalProperties | None = None
    ) -> None: ...

class BandFormatList(Serialisable):
    tagname: ClassVar[str]
    bandFmt: Incomplete
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, bandFmt=()) -> None: ...

class _SurfaceChartBase(ChartBase):
    wireframe: NestedBool[Literal[True]]
    ser: Incomplete
    bandFmts: Typed[BandFormatList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(
        self,
        wireframe: _HasTagAndGet[_ConvertibleToBool | None] | _ConvertibleToBool | None = None,
        ser=(),
        bandFmts: BandFormatList | None = None,
        **kw,
    ) -> None: ...

class SurfaceChart3D(_SurfaceChartBase, _3DBase):
    tagname: ClassVar[str]
    # Same as parent
    # wireframe = _SurfaceChartBase.wireframe
    # ser = _SurfaceChartBase.ser
    # bandFmts = _SurfaceChartBase.bandFmts
    extLst: Typed[ExtensionList, Literal[True]]
    x_axis: Typed[TextAxis, Literal[False]]
    y_axis: Typed[NumericAxis, Literal[False]]
    z_axis: Typed[SeriesAxis, Literal[False]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, **kw) -> None: ...

class SurfaceChart(SurfaceChart3D):
    tagname: ClassVar[str]
    # Same as parent and grandparent
    # wireframe = _SurfaceChartBase.wireframe
    # ser = _SurfaceChartBase.ser
    # bandFmts = _SurfaceChartBase.bandFmts
    extLst: Typed[ExtensionList, Literal[True]]
    __elements__: ClassVar[tuple[str, ...]]
    def __init__(self, **kw) -> None: ...
