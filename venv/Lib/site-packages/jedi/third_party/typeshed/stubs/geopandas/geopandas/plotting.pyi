from _typeshed import Incomplete
from collections.abc import Collection, Hashable, Iterable, Mapping, Sequence
from typing import Literal, overload
from typing_extensions import TypeAlias

import numpy as np
import pandas as pd
from matplotlib.axes import Axes  # type: ignore[import-not-found]
from matplotlib.colors import Colormap, Normalize  # type: ignore[import-not-found]
from matplotlib.typing import ColorType  # type: ignore[import-not-found]
from numpy.typing import ArrayLike, NDArray
from pandas.plotting import PlotAccessor

from .geodataframe import GeoDataFrame
from .geoseries import GeoSeries

_ColorOrColors: TypeAlias = ColorType | Sequence[ColorType] | ArrayLike

def plot_series(
    s: GeoSeries,
    cmap: str | Colormap | None = None,
    color: _ColorOrColors | None = None,
    ax: Axes | None = None,
    figsize: tuple[float, float] | None = None,
    aspect: Literal["auto", "equal"] | float | None = "auto",
    autolim: bool = True,
    *,
    # Extracted from `**style_kwds`
    vmin: float = ...,
    vmax: float = ...,
    facecolor: _ColorOrColors | None = None,
    norm: Normalize | None = None,
    **style_kwds,
) -> Axes: ...

# IMPORTANT: keep roughly in sync with `GeoplotAccessor` methods below
def plot_dataframe(
    df: GeoDataFrame,
    column: Hashable | None = None,
    cmap: str | Colormap | None = None,
    color: _ColorOrColors | None = None,
    ax: Axes | None = None,
    cax: Axes | None = None,
    categorical: bool = False,
    legend: bool = False,
    scheme: str | None = None,
    k: int = 5,
    vmin: float | None = None,
    vmax: float | None = None,
    markersize: str | float | Iterable[float] | ArrayLike | None = None,
    figsize: tuple[float, float] | None = None,
    legend_kwds: dict[str, Incomplete] | None = None,
    categories: Iterable[Hashable] | None = None,
    classification_kwds: dict[str, Incomplete] | None = None,
    missing_kwds: dict[str, Incomplete] | None = None,
    aspect: Literal["auto", "equal"] | float | None = "auto",
    autolim: bool = True,
    *,
    # Extracted from `**style_kwds`
    norm: Normalize | None = None,
    alpha: float = 1,
    facecolor: _ColorOrColors | None = None,
    edgecolor: _ColorOrColors | None = None,
    linewidth: float = ...,
    label: str = "NaN",
    **style_kwds,
) -> Axes: ...

# IMPORTANT: keep roughly in sync with `plot_dataframe`
class GeoplotAccessor(PlotAccessor):
    # The first 3 overloads of calls are from pandas, the last overload is geopandas specific
    @overload  # type: ignore[override]
    def __call__(
        self,
        x: Hashable = ...,
        y: Hashable | Sequence[Hashable] = ...,
        *,
        kind: Literal["line", "bar", "barh", "hist", "box", "kde", "density", "area", "pie", "scatter", "hexbin"],
        ax: Axes | None = None,
        subplots: Literal[False] = False,
        sharex: bool | None = None,
        sharey: bool | None = None,
        layout: tuple[int, int] | None = None,
        figsize: tuple[float, float] | None = None,
        use_index: bool = True,
        title: str | None = None,
        grid: bool | None = None,
        legend: bool | Literal["reverse"] = True,
        style: str | Sequence[str] | Mapping[Incomplete, str] | None = None,
        logx: bool | Literal["sym"] = False,
        logy: bool | Literal["sym"] = False,
        loglog: bool | Literal["sym"] = False,
        xticks: Sequence[float] | None = None,
        yticks: Sequence[float] | None = None,
        xlim: tuple[float, float] | list[float] | None = None,
        ylim: tuple[float, float] | list[float] | None = None,
        xlabel: str | None = None,
        ylabel: str | None = None,
        rot: float | None = None,
        fontsize: float | None = None,
        cmap: str | Colormap | None = None,  # also accepts `colormap` but plot_dataframe uses `cmap`
        colorbar: bool | None = None,
        position: float = 0.5,
        table: bool | pd.Series[Incomplete] | pd.DataFrame = False,
        yerr: pd.DataFrame | pd.Series[float] | ArrayLike | Mapping[Incomplete, ArrayLike] | str = ...,
        xerr: pd.DataFrame | pd.Series[float] | ArrayLike | Mapping[Incomplete, ArrayLike] | str = ...,
        stacked: bool = ...,  # default value depends on kind
        secondary_y: bool | Sequence[Hashable] = False,
        mark_right: bool = True,
        include_bool: bool = False,
        backend: str | None = None,
        **kwargs,
    ) -> Axes: ...
    @overload
    def __call__(
        self,
        x: Hashable = ...,
        y: Hashable | Sequence[Hashable] = ...,
        *,
        kind: Literal["line", "bar", "barh", "hist", "kde", "density", "area", "pie", "scatter", "hexbin"],
        ax: Sequence[Axes] | None = None,
        subplots: Literal[True] | Iterable[Iterable[Hashable]],
        sharex: bool | None = None,
        sharey: bool | None = None,
        layout: tuple[int, int] | None = None,
        figsize: tuple[float, float] | None = None,
        use_index: bool = True,
        title: str | Collection[str] | None = None,
        grid: bool | None = None,
        legend: bool | Literal["reverse"] = True,
        style: str | Sequence[str] | Mapping[Incomplete, str] | None = None,
        logx: bool | Literal["sym"] = False,
        logy: bool | Literal["sym"] = False,
        loglog: bool | Literal["sym"] = False,
        xticks: Sequence[float] | None = None,
        yticks: Sequence[float] | None = None,
        xlim: tuple[float, float] | list[float] | None = None,
        ylim: tuple[float, float] | list[float] | None = None,
        xlabel: str | None = None,
        ylabel: str | None = None,
        rot: float | None = None,
        fontsize: float | None = None,
        cmap: str | Colormap | None = None,  # also accepts `colormap` but plot_dataframe uses `cmap`
        colorbar: bool | None = None,
        position: float = 0.5,
        table: bool | pd.Series[Incomplete] | pd.DataFrame = False,
        yerr: pd.DataFrame | pd.Series[float] | ArrayLike | Mapping[Incomplete, ArrayLike] | str = ...,
        xerr: pd.DataFrame | pd.Series[float] | ArrayLike | Mapping[Incomplete, ArrayLike] | str = ...,
        stacked: bool = ...,  # default value depends on kind
        secondary_y: bool | Sequence[Hashable] = False,
        mark_right: bool = True,
        include_bool: bool = False,
        backend: str | None = None,
        **kwargs,
    ) -> NDArray[np.object_]: ...  # should be NDArray[Axes] but it is not supported
    @overload
    def __call__(
        self,
        x: Hashable = ...,
        y: Hashable | Sequence[Hashable] = ...,
        *,
        kind: Literal["box"],
        ax: Sequence[Axes] | None = None,
        subplots: Literal[True] | Iterable[Iterable[Hashable]],
        sharex: bool | None = None,
        sharey: bool | None = None,
        layout: tuple[int, int] | None = None,
        figsize: tuple[float, float] | None = None,
        use_index: bool = True,
        title: str | Collection[str] | None = None,
        grid: bool | None = None,
        legend: bool | Literal["reverse"] = True,
        style: str | Sequence[str] | Mapping[Incomplete, str] | None = None,
        logx: bool | Literal["sym"] = False,
        logy: bool | Literal["sym"] = False,
        loglog: bool | Literal["sym"] = False,
        xticks: Sequence[float] | None = None,
        yticks: Sequence[float] | None = None,
        xlim: tuple[float, float] | list[float] | None = None,
        ylim: tuple[float, float] | list[float] | None = None,
        xlabel: str | None = None,
        ylabel: str | None = None,
        rot: float | None = None,
        fontsize: float | None = None,
        cmap: str | Colormap | None = None,  # also accepts `colormap` but plot_dataframe uses `cmap`
        colorbar: bool | None = None,
        position: float = 0.5,
        table: bool | pd.Series[Incomplete] | pd.DataFrame = False,
        yerr: pd.DataFrame | pd.Series[float] | ArrayLike | Mapping[Incomplete, ArrayLike] | str = ...,
        xerr: pd.DataFrame | pd.Series[float] | ArrayLike | Mapping[Incomplete, ArrayLike] | str = ...,
        stacked: bool = ...,  # default value depends on kind
        secondary_y: bool | Sequence[Hashable] = False,
        mark_right: bool = True,
        include_bool: bool = False,
        backend: str | None = None,
        **kwargs,
    ) -> pd.Series[Axes]: ...  # type: ignore[type-var] # pyright: ignore[reportInvalidTypeArguments]
    @overload
    def __call__(
        self,
        column: Hashable | pd.Series | pd.Index | NDArray | None = None,
        cmap: str | Colormap | None = None,
        color: _ColorOrColors | None = None,
        ax: Axes | None = None,
        cax: Axes | None = None,
        categorical: bool = False,
        legend: bool = False,
        scheme: str | None = None,
        k: int = 5,
        vmin: float | None = None,
        vmax: float | None = None,
        markersize: str | float | Iterable[float] | ArrayLike | None = None,
        figsize: tuple[float, float] | None = None,
        legend_kwds: dict[str, Incomplete] | None = None,
        categories: Iterable[Hashable] | None = None,
        classification_kwds: dict[str, Incomplete] | None = None,
        missing_kwds: dict[str, Incomplete] | None = None,
        aspect: Literal["auto", "equal"] | float | None = "auto",
        *,
        kind: Literal["geo"] = "geo",
        # Extracted from `**style_kwds`
        norm: Normalize | None = None,
        alpha: float = 1,
        facecolor: _ColorOrColors | None = None,
        edgecolor: _ColorOrColors | None = None,
        linewidth: float = ...,
        label: str = "NaN",
        **style_kwds,
    ) -> Axes: ...
    def geo(
        self,
        column: Hashable | pd.Series | pd.Index | NDArray | None = None,
        cmap: str | Colormap | None = None,
        color: _ColorOrColors | None = None,
        ax: Axes | None = None,
        cax: Axes | None = None,
        categorical: bool = False,
        legend: bool = False,
        scheme: str | None = None,
        k: int = 5,
        vmin: float | None = None,
        vmax: float | None = None,
        markersize: str | float | Iterable[float] | ArrayLike | None = None,
        figsize: tuple[float, float] | None = None,
        legend_kwds: dict[str, Incomplete] | None = None,
        categories: Iterable[Hashable] | None = None,
        classification_kwds: dict[str, Incomplete] | None = None,
        missing_kwds: dict[str, Incomplete] | None = None,
        aspect: Literal["auto", "equal"] | float | None = "auto",
        *,
        # Extracted from `**style_kwds`
        norm: Normalize | None = None,
        alpha: float = 1,
        facecolor: _ColorOrColors | None = None,
        edgecolor: _ColorOrColors | None = None,
        linewidth: float = ...,
        label: str = "NaN",
        **style_kwds,
    ) -> Axes: ...
