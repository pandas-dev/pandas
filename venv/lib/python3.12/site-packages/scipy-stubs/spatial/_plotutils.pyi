from typing import Any, TypeAlias, TypedDict
from typing_extensions import Unpack

from ._qhull import ConvexHull, Delaunay, Voronoi

__all__ = ["convex_hull_plot_2d", "delaunay_plot_2d", "voronoi_plot_2d"]

_Axes: TypeAlias = object
_Figure: TypeAlias = Any

class _VoronoiPlotKwargs(TypedDict, total=False):
    show_points: bool
    show_vertices: bool
    line_colors: str
    line_width: float
    line_alpha: float
    point_size: float

def delaunay_plot_2d(tri: Delaunay, ax: _Axes | None = None) -> _Figure: ...
def convex_hull_plot_2d(hull: ConvexHull, ax: _Axes | None = None) -> _Figure: ...
def voronoi_plot_2d(vor: Voronoi, ax: _Axes | None = None, **kwargs: Unpack[_VoronoiPlotKwargs]) -> _Figure: ...
