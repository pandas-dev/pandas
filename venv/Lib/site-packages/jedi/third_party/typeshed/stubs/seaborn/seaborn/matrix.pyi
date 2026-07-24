from _typeshed import Incomplete
from collections.abc import Hashable, Iterable, Mapping, Sequence
from typing import Literal, TypedDict, type_check_only
from typing_extensions import Self, TypeAlias

import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.colors import Colormap, ListedColormap, Normalize
from matplotlib.gridspec import GridSpec
from matplotlib.typing import ColorType
from numpy._typing import _ArrayLikeInt_co
from numpy.typing import ArrayLike, NDArray

from .axisgrid import Grid

# pandas._typing.ListLikeU is partially Unknown
_ListLikeU: TypeAlias = Sequence[Incomplete] | NDArray[Incomplete] | pd.Series[Incomplete] | pd.Index[Incomplete]
_ConvertibleToDataFrame: TypeAlias = (
    _ListLikeU
    | pd.DataFrame
    | dict[Incomplete, Incomplete]
    | Iterable[_ListLikeU | tuple[Hashable, _ListLikeU] | dict[Incomplete, Incomplete]]
    | None
)
_FlatOrNestedSequenceOfColors: TypeAlias = (
    Sequence[ColorType]
    | Sequence[Iterable[ColorType]]
    | NDArray[Incomplete]
    | pd.Index[Incomplete]
    | pd.Series[Incomplete]
    | pd.DataFrame
)

__all__ = ["heatmap", "clustermap"]

def heatmap(
    data: pd.DataFrame | ArrayLike,
    *,
    vmin: float | None = None,
    vmax: float | None = None,
    cmap: str | list[ColorType] | Colormap | None = None,
    center: float | None = None,
    robust: bool = False,
    annot: bool | ArrayLike | None = None,
    fmt: str = ".2g",
    annot_kws: dict[str, Incomplete] | None = None,
    linewidths: float = 0,
    linecolor: ColorType = "white",
    cbar: bool = True,
    cbar_kws: dict[str, Incomplete] | None = None,
    cbar_ax: Axes | None = None,
    square: bool = False,
    xticklabels: Literal["auto"] | bool | int | Sequence[str] = "auto",
    yticklabels: Literal["auto"] | bool | int | Sequence[str] = "auto",
    mask: NDArray[np.bool_] | pd.DataFrame | None = None,
    ax: Axes | None = None,
    # Kwargs below passed to matplotlib.axes.Axes.pcolormesh
    alpha: float | None = None,
    norm: str | Normalize | None = None,
    shading: Literal["flat", "nearest", "gouraud", "auto"] | None = None,
    antialiased: bool = False,
    **kwargs,
) -> Axes: ...
@type_check_only
class _Dendogram(TypedDict):
    icoord: list[list[float]]
    dcoord: list[list[float]]
    ivl: list[str]
    leaves: list[int]
    color_list: list[str]
    leaves_color_list: list[str]

class _DendrogramPlotter:
    axis: int
    array: NDArray[np.floating]
    data: pd.DataFrame
    shape: tuple[int, int]
    metric: str
    method: str
    label: bool
    rotate: bool
    linkage: NDArray[np.floating]
    dendrogram: _Dendogram
    xticks: list[float] | NDArray[np.floating]
    yticks: list[float] | NDArray[np.floating]
    xticklabels: list[str]
    yticklabels: list[str]
    ylabel: str
    xlabel: str
    dependent_coord: list[list[float]]
    independent_coord: list[list[float]]
    def __init__(
        self,
        data: pd.DataFrame,
        linkage: NDArray[np.floating] | None,
        metric: str,
        method: str,
        axis: int,
        label: bool,
        rotate: bool,
    ) -> None: ...
    @property
    def calculated_linkage(self) -> NDArray[np.float64]: ...
    def calculate_dendrogram(self) -> _Dendogram: ...
    @property
    def reordered_ind(self) -> list[int]: ...
    def plot(self, ax: Axes, tree_kws: dict[str, Incomplete]) -> Self: ...

def dendrogram(
    data: pd.DataFrame,
    *,
    linkage: NDArray[np.floating] | None = None,
    axis: int = 1,
    label: bool = True,
    metric: str = "euclidean",
    method: str = "average",
    rotate: bool = False,
    tree_kws: dict[str, Incomplete] | None = None,
    ax: Axes | None = None,
) -> _DendrogramPlotter: ...

class ClusterGrid(Grid):
    data: pd.DataFrame
    data2d: pd.DataFrame
    mask: pd.DataFrame
    row_colors: list[list[tuple[float, float, float]]] | None
    row_color_labels: list[str] | None
    col_colors: list[list[tuple[float, float, float]]] | None
    col_color_labels: list[str] | None
    gs: GridSpec
    ax_row_dendrogram: Axes
    ax_col_dendrogram: Axes
    ax_row_colors: Axes | None
    ax_col_colors: Axes | None
    ax_heatmap: Axes
    ax_cbar: Axes | None
    cax: Axes | None
    cbar_pos: tuple[float, float, float, float] | None
    dendrogram_row: _DendrogramPlotter | None
    dendrogram_col: _DendrogramPlotter | None
    def __init__(
        self,
        data: _ConvertibleToDataFrame,
        pivot_kws: Mapping[str, Incomplete] | None = None,
        z_score: int | None = None,
        standard_scale: int | None = None,
        figsize: tuple[float, float] | None = None,
        row_colors: _FlatOrNestedSequenceOfColors | None = None,
        col_colors: _FlatOrNestedSequenceOfColors | None = None,
        mask: NDArray[np.bool_] | pd.DataFrame | None = None,
        dendrogram_ratio: float | tuple[float, float] | None = None,
        colors_ratio: float | tuple[float, float] | None = None,
        cbar_pos: tuple[float, float, float, float] | None = None,
    ) -> None: ...
    def format_data(
        self,
        data: pd.DataFrame,
        pivot_kws: Mapping[str, Incomplete] | None,
        z_score: int | None = None,
        standard_scale: int | None = None,
    ) -> pd.DataFrame: ...
    @staticmethod
    def z_score(data2d: pd.DataFrame, axis: int = 1) -> pd.DataFrame: ...
    @staticmethod
    def standard_scale(data2d: pd.DataFrame, axis: int = 1) -> pd.DataFrame: ...
    def dim_ratios(self, colors: ArrayLike | None, dendrogram_ratio: float, colors_ratio: float) -> list[float]: ...
    @staticmethod
    def color_list_to_matrix_and_cmap(
        colors: _FlatOrNestedSequenceOfColors, ind: _ArrayLikeInt_co, axis: int = 0
    ) -> tuple[NDArray[np.int_], ListedColormap]: ...
    def plot_dendrograms(
        self,
        row_cluster: bool,
        col_cluster: bool,
        metric: str,
        method: str,
        row_linkage: NDArray[np.floating] | None,
        col_linkage: NDArray[np.floating] | None,
        tree_kws: dict[str, Incomplete] | None,
    ) -> None: ...
    def plot_colors(self, xind: _ArrayLikeInt_co, yind: _ArrayLikeInt_co, **kws) -> None: ...
    def plot_matrix(self, colorbar_kws: dict[str, Incomplete], xind: _ArrayLikeInt_co, yind: _ArrayLikeInt_co, **kws) -> None: ...
    def plot(
        self,
        metric: str,
        method: str,
        colorbar_kws: dict[str, Incomplete] | None,
        row_cluster: bool,
        col_cluster: bool,
        row_linkage: NDArray[np.floating] | None,
        col_linkage: NDArray[np.floating] | None,
        tree_kws: dict[str, Incomplete] | None,
        **kws,
    ) -> Self: ...

def clustermap(
    data: _ConvertibleToDataFrame,
    *,
    pivot_kws: dict[str, Incomplete] | None = None,
    method: str = "average",
    metric: str = "euclidean",
    z_score: int | None = None,
    standard_scale: int | None = None,
    figsize: tuple[float, float] | None = (10, 10),
    cbar_kws: dict[str, Incomplete] | None = None,
    row_cluster: bool = True,
    col_cluster: bool = True,
    row_linkage: NDArray[np.floating] | None = None,
    col_linkage: NDArray[np.floating] | None = None,
    row_colors: _FlatOrNestedSequenceOfColors | None = None,
    col_colors: _FlatOrNestedSequenceOfColors | None = None,
    mask: NDArray[np.bool_] | pd.DataFrame | None = None,
    dendrogram_ratio: float | tuple[float, float] = 0.2,
    colors_ratio: float | tuple[float, float] = 0.03,
    cbar_pos: tuple[float, float, float, float] | None = (0.02, 0.8, 0.05, 0.18),
    tree_kws: dict[str, Incomplete] | None = None,
    **kwargs,
) -> ClusterGrid: ...
