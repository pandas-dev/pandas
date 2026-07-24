from typing import Any, Literal, overload

from matplotlib.axes import Axes  # type: ignore[import-not-found]
from matplotlib.lines import Line2D  # type: ignore[import-not-found]
from matplotlib.patches import PathPatch  # type: ignore[import-not-found]
from matplotlib.typing import ColorType  # type: ignore[import-not-found]

from .geometry import LinearRing, LineString, MultiLineString, MultiPolygon, Polygon
from .lib import Geometry

def patch_from_polygon(polygon: Polygon | MultiPolygon, **kwargs: Any) -> PathPatch: ...
@overload
def plot_polygon(
    polygon: Polygon | MultiPolygon,
    ax: Axes | None = None,
    add_points: Literal[True] = True,
    color: ColorType | None = None,
    facecolor: ColorType | None = None,
    edgecolor: ColorType | None = None,
    linewidth: float | None = None,
    **kwargs: Any,
) -> tuple[PathPatch, Line2D]: ...
@overload
def plot_polygon(
    polygon: Polygon | MultiPolygon,
    ax: Axes | None = None,
    *,
    add_points: Literal[False],
    color: ColorType | None = None,
    facecolor: ColorType | None = None,
    edgecolor: ColorType | None = None,
    linewidth: float | None = None,
    **kwargs: Any,
) -> PathPatch: ...
@overload
def plot_polygon(
    polygon: Polygon | MultiPolygon,
    ax: Axes | None,
    add_points: Literal[False],
    color: ColorType | None = None,
    facecolor: ColorType | None = None,
    edgecolor: ColorType | None = None,
    linewidth: float | None = None,
    **kwargs: Any,
) -> PathPatch: ...
@overload
def plot_line(
    line: LineString | LinearRing | MultiLineString,
    ax: Axes | None = None,
    add_points: Literal[True] = True,
    color: ColorType | None = None,
    linewidth: float = 2,
    **kwargs: Any,
) -> tuple[PathPatch, Line2D]: ...
@overload
def plot_line(
    line: LineString | LinearRing | MultiLineString,
    ax: Axes | None = None,
    *,
    add_points: Literal[False],
    color: ColorType | None = None,
    linewidth: float = 2,
    **kwargs: Any,
) -> PathPatch: ...
@overload
def plot_line(
    line: LineString | LinearRing | MultiLineString,
    ax: Axes | None,
    add_points: Literal[False],
    color: ColorType | None = None,
    linewidth: float = 2,
    **kwargs: Any,
) -> PathPatch: ...
def plot_points(
    geom: Geometry, ax: Axes | None = None, color: ColorType | None = None, marker: str = "o", **kwargs: Any
) -> Line2D: ...
