import os
from _typeshed import Incomplete
from collections.abc import Callable, Generator, Iterable, Mapping
from typing import IO, Any, Literal, TypeVar
from typing_extensions import Concatenate, ParamSpec, Self, TypeAlias, deprecated

import numpy as np
from matplotlib.artist import Artist
from matplotlib.axes import Axes
from matplotlib.backend_bases import MouseEvent, RendererBase
from matplotlib.colors import Colormap
from matplotlib.figure import Figure
from matplotlib.font_manager import FontProperties
from matplotlib.gridspec import SubplotSpec
from matplotlib.legend import Legend
from matplotlib.patches import Patch
from matplotlib.path import Path as mpl_Path
from matplotlib.patheffects import AbstractPathEffect
from matplotlib.scale import ScaleBase
from matplotlib.text import Text
from matplotlib.transforms import Bbox, BboxBase, Transform, TransformedPath
from matplotlib.typing import ColorType, LineStyleType, MarkerType
from numpy.typing import ArrayLike, NDArray
from pandas import DataFrame, Series

from ._core.typing import ColumnName, DataSource, NormSpec, SupportsDataFrame
from .palettes import _RGBColorPalette
from .utils import _DataSourceWideForm, _Palette, _Vector

__all__ = ["FacetGrid", "PairGrid", "JointGrid", "pairplot", "jointplot"]

_P = ParamSpec("_P")
_R = TypeVar("_R")

_LiteralFont: TypeAlias = Literal["xx-small", "x-small", "small", "medium", "large", "x-large", "xx-large"]

class _BaseGrid:
    def set(
        self,
        *,
        # Keywords follow `matplotlib.axes.Axes.set`. Each keyword <KW> corresponds to a `set_<KW>` method
        adjustable: Literal["box", "datalim"] = ...,
        agg_filter: Callable[[ArrayLike, float], tuple[NDArray[np.floating[Any]], float, float]] | None = ...,
        alpha: float | None = ...,
        anchor: str | tuple[float, float] = ...,
        animated: bool = ...,
        aspect: float | Literal["auto", "equal"] = ...,
        autoscale_on: bool = ...,
        autoscalex_on: bool = ...,
        autoscaley_on: bool = ...,
        axes_locator: Callable[[Axes, RendererBase], Bbox] = ...,
        axisbelow: bool | Literal["line"] = ...,
        box_aspect: float | None = ...,
        clip_box: BboxBase | None = ...,
        clip_on: bool = ...,
        clip_path: Patch | mpl_Path | TransformedPath | None = ...,
        facecolor: ColorType | None = ...,
        frame_on: bool = ...,
        gid: str | None = ...,
        in_layout: bool = ...,
        label: object = ...,
        mouseover: bool = ...,
        navigate: bool = ...,
        path_effects: list[AbstractPathEffect] = ...,
        picker: bool | float | Callable[[Artist, MouseEvent], tuple[bool, dict[Any, Any]]] | None = ...,
        position: Bbox | tuple[float, float, float, float] = ...,
        prop_cycle=...,  # TODO: use cycler.Cycler when cycler gets typed
        rasterization_zorder: float | None = ...,
        rasterized: bool = ...,
        sketch_params: float | None = ...,
        snap: bool | None = ...,
        subplotspec: SubplotSpec = ...,
        title: str = ...,
        transform: Transform | None = ...,
        url: str | None = ...,
        visible: bool = ...,
        xbound: float | None | tuple[float | None, float | None] = ...,
        xlabel: str = ...,
        xlim: float | None | tuple[float | None, float | None] = ...,
        xmargin: float = ...,
        xscale: str | ScaleBase = ...,
        xticklabels: Iterable[str | Text] = ...,
        xticks: ArrayLike = ...,
        ybound: float | None | tuple[float | None, float | None] = ...,
        ylabel: str = ...,
        ylim: float | None | tuple[float | None, float | None] = ...,
        ymargin: float = ...,
        yscale: str | ScaleBase = ...,
        yticklabels: Iterable[str | Text] = ...,
        yticks: ArrayLike = ...,
        zorder: float = ...,
        **kwargs: Any,
    ) -> Self: ...
    @property
    @deprecated("Attribute `fig` is deprecated in favor of `figure`")
    def fig(self) -> Figure: ...
    @property
    def figure(self) -> Figure: ...
    def apply(self, func: Callable[Concatenate[Self, _P], object], *args: _P.args, **kwargs: _P.kwargs) -> Self: ...
    def pipe(self, func: Callable[Concatenate[Self, _P], _R], *args: _P.args, **kwargs: _P.kwargs) -> _R: ...
    def savefig(
        self,
        # Signature follows `matplotlib.figure.Figure.savefig`
        fname: str | os.PathLike[Any] | IO[Any],
        *,
        transparent: bool | None = None,
        dpi: float | Literal["figure"] | None = 96,
        facecolor: ColorType | Literal["auto"] | None = "auto",
        edgecolor: ColorType | Literal["auto"] | None = "auto",
        orientation: Literal["landscape", "portrait"] = "portrait",
        format: str | None = None,
        bbox_inches: Literal["tight"] | Bbox | None = "tight",
        pad_inches: float | Literal["layout"] | None = None,
        backend: str | None = None,
        **kwargs: Any,
    ) -> None: ...

class Grid(_BaseGrid):
    def __init__(self) -> None: ...
    def tight_layout(
        self,
        *,
        # Keywords follow `matplotlib.figure.Figure.tight_layout`
        pad: float = 1.08,
        h_pad: float | None = None,
        w_pad: float | None = None,
        rect: tuple[float, float, float, float] | None = None,
    ) -> Self: ...
    def add_legend(
        self,
        # Cannot use precise key type with union for legend_data because of invariant Mapping keys
        legend_data: Mapping[Any, Artist] | None = None,
        title: str | None = None,
        label_order: list[str] | None = None,
        adjust_subtitles: bool = False,
        *,
        # Keywords follow `matplotlib.legend.Legend`
        loc: str | int | tuple[float, float] | None = None,
        numpoints: int | None = None,
        markerscale: float | None = None,
        markerfirst: bool = True,
        reverse: bool = False,
        scatterpoints: int | None = None,
        scatteryoffsets: Iterable[float] | None = None,
        prop: FontProperties | dict[str, Any] | None = None,
        fontsize: int | _LiteralFont | None = None,
        labelcolor: str | Iterable[str] | None = None,
        borderpad: float | None = None,
        labelspacing: float | None = None,
        handlelength: float | None = None,
        handleheight: float | None = None,
        handletextpad: float | None = None,
        borderaxespad: float | None = None,
        columnspacing: float | None = None,
        ncols: int = 1,
        mode: Literal["expand"] | None = None,
        fancybox: bool | None = None,
        shadow: bool | dict[str, int] | dict[str, float] | None = None,
        title_fontsize: int | _LiteralFont | None = None,
        framealpha: float | None = None,
        edgecolor: ColorType | None = None,
        facecolor: ColorType | None = None,
        bbox_to_anchor: BboxBase | tuple[float, float] | tuple[float, float, float, float] | None = None,
        bbox_transform: Transform | None = None,
        frameon: bool | None = None,
        handler_map: None = None,
        title_fontproperties: FontProperties | None = None,
        alignment: Literal["center", "left", "right"] = "center",
        ncol: int = 1,
        draggable: bool = False,
    ) -> Self: ...
    @property
    def legend(self) -> Legend | None: ...
    def tick_params(
        self,
        axis: Literal["x", "y", "both"] = "both",
        *,
        # Keywords follow `matplotlib.axes.Axes.tick_params`
        which: Literal["major", "minor", "both"] = "major",
        reset: bool = False,
        direction: Literal["in", "out", "inout"] = ...,
        length: float = ...,
        width: float = ...,
        color: ColorType = ...,
        pad: float = ...,
        labelsize: float | str = ...,
        labelcolor: ColorType = ...,
        labelfontfamily: str = ...,
        colors: ColorType = ...,
        zorder: float = ...,
        bottom: bool = ...,
        top: bool = ...,
        left: bool = ...,
        right: bool = ...,
        labelbottom: bool = ...,
        labeltop: bool = ...,
        labelleft: bool = ...,
        labelright: bool = ...,
        labelrotation: float = ...,
        grid_color: ColorType = ...,
        grid_alpha: float = ...,
        grid_linewidth: float = ...,
        grid_linestyle: str = ...,
        **kwargs: Any,
    ) -> Self: ...

class FacetGrid(Grid):
    data: DataFrame
    row_names: list[Any]
    col_names: list[Any]
    hue_names: list[Any] | None
    hue_kws: dict[str, Any]
    def __init__(
        self,
        data: DataFrame | SupportsDataFrame,
        *,
        row: str | None = None,
        col: str | None = None,
        hue: str | None = None,
        col_wrap: int | None = None,
        sharex: bool | Literal["col", "row"] = True,
        sharey: bool | Literal["col", "row"] = True,
        height: float = 3,
        aspect: float = 1,
        palette: _Palette | None = None,
        row_order: Iterable[Any] | None = None,
        col_order: Iterable[Any] | None = None,
        hue_order: Iterable[Any] | None = None,
        hue_kws: dict[str, Any] | None = None,
        dropna: bool = False,
        legend_out: bool = True,
        despine: bool = True,
        margin_titles: bool = False,
        xlim: tuple[float, float] | None = None,
        ylim: tuple[float, float] | None = None,
        subplot_kws: dict[str, Any] | None = None,
        gridspec_kws: dict[str, Any] | None = None,
    ) -> None: ...
    def facet_data(self) -> Generator[tuple[tuple[int, int, int], DataFrame]]: ...
    def map(self, func: Callable[..., object], *args: str, **kwargs: Any) -> Self: ...
    def map_dataframe(self, func: Callable[..., object], *args: str, **kwargs: Any) -> Self: ...
    def facet_axis(self, row_i: int, col_j: int, modify_state: bool = True) -> Axes: ...
    # `despine` should be kept roughly in line with `seaborn.utils.despine`
    def despine(
        self,
        *,
        ax: Axes | None = None,
        top: bool = True,
        right: bool = True,
        left: bool = False,
        bottom: bool = False,
        offset: int | Mapping[str, int] | None = None,
        trim: bool = False,
    ) -> Self: ...
    def set_axis_labels(
        self, x_var: str | None = None, y_var: str | None = None, clear_inner: bool = True, **kwargs: Any
    ) -> Self: ...
    def set_xlabels(self, label: str | None = None, clear_inner: bool = True, **kwargs: Any) -> Self: ...
    def set_ylabels(self, label: str | None = None, clear_inner: bool = True, **kwargs: Any) -> Self: ...
    def set_xticklabels(self, labels: Iterable[str | Text] | None = None, step: int | None = None, **kwargs: Any) -> Self: ...
    def set_yticklabels(self, labels: Iterable[str | Text] | None = None, **kwargs: Any) -> Self: ...
    def set_titles(
        self, template: str | None = None, row_template: str | None = None, col_template: str | None = None, **kwargs: Any
    ) -> Self: ...
    def refline(
        self,
        *,
        x: float | None = None,
        y: float | None = None,
        color: ColorType = ".5",
        linestyle: LineStyleType = "--",
        **line_kws: Any,
    ) -> Self: ...
    @property
    def axes(self) -> NDArray[Incomplete]: ...  # array of `Axes`
    @property
    def ax(self) -> Axes: ...
    @property
    def axes_dict(self) -> dict[Any, Axes]: ...

class PairGrid(Grid):
    x_vars: list[str]
    y_vars: list[str]
    square_grid: bool
    axes: NDArray[Incomplete]  # two-dimensional array of `Axes`
    data: DataFrame
    diag_sharey: bool
    diag_vars: list[str] | None
    diag_axes: list[Axes] | None
    hue_names: list[str]
    hue_vals: Series[Any]
    hue_kws: dict[str, Any]
    palette: _RGBColorPalette
    def __init__(
        self,
        data: DataFrame | SupportsDataFrame,
        *,
        hue: str | None = None,
        vars: Iterable[str] | None = None,
        x_vars: Iterable[str] | str | None = None,
        y_vars: Iterable[str] | str | None = None,
        hue_order: Iterable[str] | None = None,
        palette: _Palette | None = None,
        hue_kws: dict[str, Any] | None = None,
        corner: bool = False,
        diag_sharey: bool = True,
        height: float = 2.5,
        aspect: float = 1,
        layout_pad: float = 0.5,
        despine: bool = True,
        dropna: bool = False,
    ) -> None: ...
    def map(self, func: Callable[..., object], **kwargs: Any) -> Self: ...
    def map_lower(self, func: Callable[..., object], **kwargs: Any) -> Self: ...
    def map_upper(self, func: Callable[..., object], **kwargs: Any) -> Self: ...
    def map_offdiag(self, func: Callable[..., object], **kwargs: Any) -> Self: ...
    def map_diag(self, func: Callable[..., object], **kwargs: Any) -> Self: ...

class JointGrid(_BaseGrid):
    ax_joint: Axes
    ax_marg_x: Axes
    ax_marg_y: Axes
    x: Series[Any]
    y: Series[Any]
    hue: Series[Any]
    def __init__(
        self,
        data: DataSource | _DataSourceWideForm | None = None,
        *,
        x: ColumnName | _Vector | None = None,
        y: ColumnName | _Vector | None = None,
        hue: ColumnName | _Vector | None = None,
        height: float = 6,
        ratio: float = 5,
        space: float = 0.2,
        palette: _Palette | Colormap | None = None,
        hue_order: Iterable[ColumnName] | None = None,
        hue_norm: NormSpec = None,
        dropna: bool = False,
        xlim: float | tuple[float, float] | None = None,
        ylim: float | tuple[float, float] | None = None,
        marginal_ticks: bool = False,
    ) -> None: ...
    def plot(self, joint_func: Callable[..., object], marginal_func: Callable[..., object], **kwargs: Any) -> Self: ...
    def plot_joint(self, func: Callable[..., object], **kwargs: Any) -> Self: ...
    def plot_marginals(self, func: Callable[..., object], **kwargs: Any) -> Self: ...
    def refline(
        self,
        *,
        x: float | None = None,
        y: float | None = None,
        joint: bool = True,
        marginal: bool = True,
        color: ColorType = ".5",
        linestyle: LineStyleType = "--",
        **line_kws: Any,
    ) -> Self: ...
    def set_axis_labels(self, xlabel: str = "", ylabel: str = "", **kwargs: Any) -> Self: ...

def pairplot(
    data: DataFrame,
    *,
    hue: str | None = None,
    hue_order: Iterable[str] | None = None,
    palette: _Palette | None = None,
    vars: Iterable[str] | None = None,
    x_vars: Iterable[str] | str | None = None,
    y_vars: Iterable[str] | str | None = None,
    kind: Literal["scatter", "kde", "hist", "reg"] = "scatter",
    diag_kind: Literal["auto", "hist", "kde"] | None = "auto",
    markers: MarkerType | list[MarkerType] | None = None,
    height: float = 2.5,
    aspect: float = 1,
    corner: bool = False,
    dropna: bool = False,
    plot_kws: dict[str, Any] | None = None,
    diag_kws: dict[str, Any] | None = None,
    grid_kws: dict[str, Any] | None = None,
    size: float | None = None,  # deprecated
) -> PairGrid: ...
def jointplot(
    data: DataSource | _DataSourceWideForm | None = None,
    *,
    x: ColumnName | _Vector | None = None,
    y: ColumnName | _Vector | None = None,
    hue: ColumnName | _Vector | None = None,
    kind: Literal["scatter", "kde", "hist", "hex", "reg", "resid"] = "scatter",
    height: float = 6,
    ratio: float = 5,
    space: float = 0.2,
    dropna: bool = False,
    xlim: float | tuple[float, float] | None = None,
    ylim: float | tuple[float, float] | None = None,
    color: ColorType | None = None,
    palette: _Palette | Colormap | None = None,
    hue_order: Iterable[ColumnName] | None = None,
    hue_norm: NormSpec = None,
    marginal_ticks: bool = False,
    joint_kws: dict[str, Any] | None = None,
    marginal_kws: dict[str, Any] | None = None,
    **kwargs: Any,
) -> JointGrid: ...
