from _typeshed import Incomplete, SupportsItems
from collections.abc import Callable, Collection, Hashable, Iterable, Mapping, Sequence
from typing import Any, Generic, Literal, TypedDict, TypeVar, overload, type_check_only
from typing_extensions import TypeAlias, Unpack

import numpy as np
from matplotlib.axes import Axes  # type: ignore[import-not-found]
from matplotlib.collections import LineCollection, PathCollection  # type: ignore[import-not-found]
from matplotlib.colors import Colormap  # type: ignore[import-not-found]
from matplotlib.patches import FancyArrowPatch  # type: ignore[import-not-found]
from matplotlib.text import Text  # type: ignore[import-not-found]
from matplotlib.typing import ColorType  # type: ignore[import-not-found]
from networkx._typing import Array2D
from networkx.classes.digraph import DiGraph
from networkx.classes.graph import Graph, _Node

__all__ = [
    "display",
    "apply_matplotlib_colors",
    "draw",
    "draw_networkx",
    "draw_networkx_nodes",
    "draw_networkx_edges",
    "draw_networkx_labels",
    "draw_networkx_edge_labels",
    "draw_bipartite",
    "draw_circular",
    "draw_kamada_kawai",
    "draw_random",
    "draw_spectral",
    "draw_spring",
    "draw_planar",
    "draw_shell",
    "draw_forceatlas2",
]

_G = TypeVar("_G", bound=Graph[Any])

# types from matplotlib
_FontSize: TypeAlias = Literal["xx-small", "x-small", "small", "medium", "large", "x-large", "xx-large"] | float
_FontWeight: TypeAlias = (
    Literal[
        "ultralight",
        "light",
        "normal",
        "regular",
        "book",
        "medium",
        "roman",
        "semibold",
        "demibold",
        "demi",
        "bold",
        "heavy",
        "extra bold",
        "black",
    ]
    | int
)
_HAlign: TypeAlias = Literal["left", "center", "right"]
_VAlign: TypeAlias = Literal["baseline", "bottom", "center", "center_baseline", "top"]

@type_check_only
class _DrawNetworkxKwds(TypedDict, Generic[_Node], total=False):
    # draw nodes keywords; keep in sync with draw_networkx_nodes
    node_color: ColorType | Collection[ColorType] | Collection[float]
    cmap: str | Colormap | None
    vmin: float | None
    vmax: float | None
    linewidths: float | Collection[float] | None
    edgecolors: Literal["face", "none"] | ColorType | Collection[ColorType] | Collection[float] | None
    margins: float | tuple[float, float] | None
    # draw edges keywords; keep in sync with draw_networkx_edges
    edgelist: Collection[_Node | Hashable] | None  # (u, v, k) for multigraphs and (u, v) for simple graphs
    width: float | Collection[float]
    edge_color: ColorType | Collection[ColorType]
    style: str | Collection[str]
    arrowstyle: str | Collection[str] | None
    arrowsize: float | list[int] | list[float]
    edge_cmap: Colormap | None
    edge_vmin: float | None
    edge_vmax: float | None
    connectionstyle: str | Iterable[str]
    min_source_margin: int | Collection[int]
    min_target_margin: int | Collection[int]
    # draw labels keywords; keep in sync with draw_networkx_labels
    labels: Mapping[_Node, object] | None
    font_size: _FontSize | Mapping[_Node, _FontSize]
    font_color: ColorType | Mapping[_Node, Colormap]
    font_family: str | Mapping[_Node, str]
    font_weight: _FontWeight | Mapping[_Node, _FontWeight]
    bbox: dict[str, Any] | None
    horizontalalignment: _HAlign
    verticalalignment: _VAlign
    clip_on: bool
    # common keywords
    nodelist: Sequence[_Node] | None
    node_size: float | Collection[float]
    node_shape: str
    alpha: float | Collection[float] | None
    label: str | None
    hide_ticks: bool

def apply_matplotlib_colors(
    G: Graph[_Node],
    src_attr: str,
    dest_attr: str,
    map: str | Colormap,
    vmin: float | None = None,
    vmax: float | None = None,
    nodes: bool = True,
) -> None: ...

class CurvedArrowTextBase:
    arrow: FancyArrowPatch
    label_pos: float
    labels_horizontal: bool
    ax: Axes
    x: Incomplete
    y: Incomplete
    angle: Incomplete
    def __init__(
        self,
        arrow: FancyArrowPatch,
        *args,
        label_pos: float = 0.5,
        labels_horizontal: bool = False,
        ax: Axes | None = None,
        **kwargs,
    ) -> None: ...
    def draw(self, renderer) -> None: ...

def display(
    G: _G,
    canvas: Axes | None = None,
    *,
    pos: str | Callable[[_G], Mapping[_Node, Collection[float]]] = ...,
    node_visible: str | bool = ...,
    node_color: str = ...,
    node_size: str | float = ...,
    node_label: str | bool = ...,
    node_shape: str = ...,
    node_alpha: str = ...,
    node_border_width: str = ...,
    node_border_color: str = ...,
    edge_visible: str | bool = ...,
    edge_width: str | int = ...,
    edge_color: str | ColorType = ...,
    edge_label: str = ...,
    edge_style: str = ...,
    edge_alpha: str | float = ...,
    arrowstyle: str = ...,
    arrowsize: str | int = ...,
    edge_curvature: str = ...,
    edge_source_margin: str | int = ...,
    edge_target_margin: str | int = ...,
    hide_ticks: bool = True,
) -> _G: ...
def draw(
    G: Graph[_Node],
    pos: Mapping[_Node, Collection[float]] | None = None,
    ax: Axes | None = None,
    *,
    with_labels: bool = ...,  # default depends on whether a label argument is passed
    **kwds: Unpack[_DrawNetworkxKwds[_Node]],
) -> None: ...
def draw_networkx(
    G: Graph[_Node],
    pos: Mapping[_Node, Collection[float]] | None = None,
    arrows: bool | None = None,
    with_labels: bool = True,
    *,
    ax: Axes | None = None,
    **kwds: Unpack[_DrawNetworkxKwds[_Node]],
) -> None: ...
def draw_networkx_nodes(  # keep in sync with _DrawNetworkxKwds above
    G: Graph[_Node],
    pos: Mapping[_Node, Collection[float]],
    nodelist: Collection[_Node] | None = None,
    node_size: float | Collection[float] = 300,
    node_color: ColorType | Collection[ColorType] | Collection[float] = "#1f78b4",
    node_shape: str = "o",
    alpha: float | Collection[float] | None = None,
    cmap: str | Colormap | None = None,
    vmin: float | None = None,
    vmax: float | None = None,
    ax: Axes | None = None,
    linewidths: float | Collection[float] | None = None,
    edgecolors: Literal["face", "none"] | ColorType | Collection[ColorType] | Collection[float] | None = None,
    label: str | None = None,
    margins: float | tuple[float, float] | None = None,
    hide_ticks: bool = True,
) -> PathCollection: ...
@overload  # arrows=None -> LineCollection if G is undirected, list[FancyArrowPatch] if G is directed
def draw_networkx_edges(  # keep in sync with _DrawNetworkxKwds above
    G: Graph[_Node],
    pos: Mapping[_Node, Collection[float]],
    edgelist: Collection[_Node | Hashable] | None = None,  # (u, v, k) for multigraphs and (u, v) for simple graphs
    width: float | Collection[float] = 1.0,
    edge_color: ColorType | Collection[ColorType] = "k",
    style: str | Collection[str] = "solid",
    alpha: float | Collection[float] | None = None,
    arrowstyle: str | Collection[str] | None = None,
    arrowsize: float | list[int] | list[float] = 10,  # documented as int, mpl accepts float
    edge_cmap: Colormap | None = None,
    edge_vmin: float | None = None,
    edge_vmax: float | None = None,
    ax: Axes | None = None,
    arrows: None = None,
    label: str | None = None,  # documented as str, mpl accepts any object as it calls str on it
    node_size: float | Collection[float] = 300,
    nodelist: Sequence[_Node] | None = None,
    node_shape: str = "o",
    connectionstyle: str | Iterable[str] = "arc3",
    min_source_margin: int | Collection[int] = 0,  # documented as int, mpl accepts float
    min_target_margin: int | Collection[int] = 0,  # documented as int, mpl accepts float
    hide_ticks: bool = True,
) -> LineCollection | list[FancyArrowPatch]: ...
@overload  # directed graph and arrows=None -> list[FancyArrowPatch]
def draw_networkx_edges(
    G: DiGraph[_Node],
    pos: Mapping[_Node, Collection[float]],
    edgelist: Collection[_Node | Hashable] | None = None,  # (u, v, k) for multigraphs and (u, v) for simple graphs
    width: float | Collection[float] = 1.0,
    edge_color: ColorType | Collection[ColorType] = "k",
    style: str | Collection[str] = "solid",
    alpha: float | Collection[float] | None = None,
    arrowstyle: str | Collection[str] | None = None,
    arrowsize: float | list[int] | list[float] = 10,  # documented as int, mpl accepts float
    edge_cmap: Colormap | None = None,
    edge_vmin: float | None = None,
    edge_vmax: float | None = None,
    ax: Axes | None = None,
    arrows: None = None,
    label: str | None = None,  # documented as str, mpl accepts any object as it calls str on it
    node_size: float | Collection[float] = 300,
    nodelist: Sequence[_Node] | None = None,
    node_shape: str = "o",
    connectionstyle: str | Iterable[str] = "arc3",
    min_source_margin: int | Collection[int] = 0,  # documented as int, mpl accepts float
    min_target_margin: int | Collection[int] = 0,  # documented as int, mpl accepts float
    hide_ticks: bool = True,
) -> list[FancyArrowPatch]: ...
@overload  # arrows=True -> list[FancyArrowPatch]
def draw_networkx_edges(
    G: Graph[_Node],
    pos: Mapping[_Node, Collection[float]],
    edgelist: Collection[_Node | Hashable] | None = None,  # (u, v, k) for multigraphs and (u, v) for simple graphs
    width: float | Collection[float] = 1.0,
    edge_color: ColorType | Collection[ColorType] = "k",
    style: str | Collection[str] = "solid",
    alpha: float | Collection[float] | None = None,
    arrowstyle: str | Collection[str] | None = None,
    arrowsize: float | list[int] | list[float] = 10,  # documented as int, mpl accepts float
    edge_cmap: Colormap | None = None,
    edge_vmin: float | None = None,
    edge_vmax: float | None = None,
    ax: Axes | None = None,
    *,
    arrows: Literal[True],
    label: str | None = None,  # documented as str, mpl accepts any object as it calls str on it
    node_size: float | Collection[float] = 300,
    nodelist: Sequence[_Node] | None = None,
    node_shape: str = "o",
    connectionstyle: str | Iterable[str] = "arc3",
    min_source_margin: int | Collection[int] = 0,  # documented as int, mpl accepts float
    min_target_margin: int | Collection[int] = 0,  # documented as int, mpl accepts float
    hide_ticks: bool = True,
) -> list[FancyArrowPatch]: ...
@overload  # arrows=False -> LineCollection
def draw_networkx_edges(
    G: Graph[_Node],
    pos: Mapping[_Node, Collection[float]],
    edgelist: Collection[_Node | Hashable] | None = None,  # (u, v, k) for multigraphs and (u, v) for simple graphs
    width: float | Collection[float] = 1.0,
    edge_color: ColorType | Collection[ColorType] = "k",
    style: str | Collection[str] = "solid",
    alpha: float | Collection[float] | None = None,
    *,
    edge_cmap: Colormap | None = None,
    edge_vmin: float | None = None,
    edge_vmax: float | None = None,
    ax: Axes | None = None,
    arrows: Literal[False],
    label: str | None = None,  # documented as str, mpl accepts any object as it calls str on it
    node_size: float | Collection[float] = 300,
    nodelist: Sequence[_Node] | None = None,
    node_shape: str = "o",
    hide_ticks: bool = True,
) -> LineCollection: ...
def draw_networkx_labels(  # keep in sync with _DrawNetworkxKwds above
    G: Graph[_Node],
    pos: Mapping[_Node, Collection[float]],
    labels: Mapping[_Node, object] | None = None,  # labels are explicitly converted to str
    font_size: _FontSize | Mapping[_Node, _FontSize] = 12,
    font_color: ColorType | Mapping[_Node, Colormap] = "k",
    font_family: str | Mapping[_Node, str] = "sans-serif",
    font_weight: _FontWeight | Mapping[_Node, _FontWeight] = "normal",
    alpha: float | Mapping[_Node, float] | None = None,
    bbox: dict[str, Any] | None = None,  # Any comes from mpl
    horizontalalignment: _HAlign = "center",  # doc is wrong, doesn't really accept array
    verticalalignment: _VAlign = "center",  # doc is wrong, doesn't really accept array
    ax: Axes | None = None,
    clip_on: bool = True,
    hide_ticks: bool = True,
) -> dict[_Node, Text]: ...
def draw_networkx_edge_labels(
    # TODO: find a way to have a covariant list for params annotated with `something | list[Incomplete]`
    G: Graph[_Node],
    pos: Mapping[_Node, Collection[float]],
    edge_labels: (
        SupportsItems[
            Collection[_Node | Hashable],  # (u, v, k) for multigraphs and (u, v) for simple graphs
            object,  # labels are explicitly converted to str and nx internally passes non-str
        ]
        | None
    ) = None,
    label_pos: float | list[Incomplete] = 0.5,
    font_size: _FontSize | list[Incomplete] = 10,
    font_color: ColorType | list[Incomplete] = "k",
    font_family: str = "sans-serif",
    font_weight: _FontWeight | list[Incomplete] = "normal",
    alpha: float | list[Incomplete] | None = None,
    bbox: dict[str, Any] | None = None,  # Any comes from mpl
    horizontalalignment: _HAlign | list[Incomplete] = "center",
    verticalalignment: _VAlign | list[Incomplete] = "center",
    ax: Axes | None = None,
    rotate: bool | list[bool] = True,
    clip_on: bool = True,
    node_size: float | Collection[float] = 300,
    nodelist: Sequence[_Node] | None = None,
    connectionstyle: str | Iterable[str] = "arc3",
    hide_ticks: bool = True,
) -> dict[tuple[_Node, _Node] | tuple[_Node, _Node, Any], Text]: ...  # Any is for multigraph key
def draw_bipartite(
    G: Graph[_Node], *, ax: Axes | None = None, with_labels: bool = ..., **kwargs: Unpack[_DrawNetworkxKwds[_Node]]
) -> None: ...
def draw_circular(
    G: Graph[_Node], *, ax: Axes | None = None, with_labels: bool = ..., **kwargs: Unpack[_DrawNetworkxKwds[_Node]]
) -> None: ...
def draw_kamada_kawai(
    G: Graph[_Node], *, ax: Axes | None = None, with_labels: bool = ..., **kwargs: Unpack[_DrawNetworkxKwds[_Node]]
) -> None: ...
def draw_random(
    G: Graph[_Node], *, ax: Axes | None = None, with_labels: bool = ..., **kwargs: Unpack[_DrawNetworkxKwds[_Node]]
) -> None: ...
def draw_spectral(
    G: Graph[_Node], *, ax: Axes | None = None, with_labels: bool = ..., **kwargs: Unpack[_DrawNetworkxKwds[_Node]]
) -> None: ...
def draw_spring(
    G: Graph[_Node], *, ax: Axes | None = None, with_labels: bool = ..., **kwargs: Unpack[_DrawNetworkxKwds[_Node]]
) -> None: ...
def draw_shell(
    G: Graph[_Node],
    nlist: Collection[Collection[_Node]] | None = None,
    *,
    ax: Axes | None = None,
    with_labels: bool = ...,
    **kwargs: Unpack[_DrawNetworkxKwds[_Node]],
) -> None: ...
def draw_planar(
    G: Graph[_Node], *, ax: Axes | None = None, with_labels: bool = ..., **kwargs: Unpack[_DrawNetworkxKwds[_Node]]
) -> None: ...
def draw_forceatlas2(
    G: Graph[_Node], *, ax: Axes | None = None, with_labels: bool = ..., **kwargs: Unpack[_DrawNetworkxKwds[_Node]]
) -> None: ...
def apply_alpha(
    colors: ColorType | Collection[ColorType] | Collection[float],
    alpha: float | Collection[float],
    elem_list: Collection[object],  # nx objects (nodes, edges, labels) but its content is not used!
    cmap: str | Colormap | None = None,
    vmin: float | None = None,
    vmax: float | None = None,
) -> Array2D[np.float64]: ...
