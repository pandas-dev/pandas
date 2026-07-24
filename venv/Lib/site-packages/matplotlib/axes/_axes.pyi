from matplotlib.axes._base import _AxesBase
from matplotlib.axes._secondary_axes import SecondaryAxis

from matplotlib.artist import Artist
from matplotlib.collections import (
    Collection,
    FillBetweenPolyCollection,
    LineCollection,
    PathCollection,
    PolyCollection,
    EventCollection,
    QuadMesh,
)
from matplotlib.colorizer import Colorizer
from matplotlib.colors import Colormap, Normalize
from matplotlib.container import (
    BarContainer, PieContainer, ErrorbarContainer, StemContainer)
from matplotlib.contour import ContourSet, QuadContourSet
from matplotlib.image import AxesImage, PcolorImage
from matplotlib.inset import InsetIndicator
from matplotlib.legend import Legend
from matplotlib.legend_handler import HandlerBase
from matplotlib.lines import Line2D, AxLine
from matplotlib.mlab import GaussianKDE
from matplotlib.patches import Rectangle, FancyArrow, Polygon, StepPatch
from matplotlib.quiver import Quiver, QuiverKey, Barbs
from matplotlib.text import Annotation, Text
from matplotlib.transforms import Transform
from matplotlib.typing import CoordsType
import matplotlib.tri as mtri
import matplotlib.table as mtable
import matplotlib.stackplot as mstack
import matplotlib.streamplot as mstream

import PIL.Image
from collections.abc import Callable, Iterable, Sequence
from typing import Any, Literal, overload
import numpy as np
from numpy.typing import ArrayLike
from matplotlib.typing import (
    ColorType, DataParamType, MarkerType, LegendLocType, LineStyleType)
import pandas as pd


class _GroupedBarReturn:
    bar_containers: list[BarContainer]
    def __init__(self, bar_containers: list[BarContainer]) -> None: ...
    def remove(self) -> None: ...

class Axes(_AxesBase):
    def get_title(self, loc: Literal["left", "center", "right"] = ...) -> str: ...
    def set_title(
        self,
        label: str,
        fontdict: dict[str, Any] | None = ...,
        loc: Literal["left", "center", "right"] | None = ...,
        pad: float | None = ...,
        *,
        y: float | None = ...,
        **kwargs
    ) -> Text: ...
    def get_legend_handles_labels(
        self, legend_handler_map: dict[type, HandlerBase] | None = ...
    ) -> tuple[list[Artist], list[Any]]: ...
    legend_: Legend | None

    @overload
    def legend(self) -> Legend: ...
    @overload
    def legend(self, handles: Iterable[Artist | tuple[Artist, ...]], labels: Iterable[str],
               *, loc: LegendLocType | None = ..., **kwargs) -> Legend: ...
    @overload
    def legend(self, *, handles: Iterable[Artist | tuple[Artist, ...]],
               loc: LegendLocType | None = ..., **kwargs) -> Legend: ...
    @overload
    def legend(self, labels: Iterable[str],
               *, loc: LegendLocType | None = ..., **kwargs) -> Legend: ...
    @overload
    def legend(self, *, loc: LegendLocType | None = ..., **kwargs) -> Legend: ...

    def inset_axes(
        self,
        bounds: tuple[float, float, float, float],
        *,
        transform: Transform | None = ...,
        zorder: float = ...,
        **kwargs
    ) -> Axes: ...
    def indicate_inset(
        self,
        bounds: tuple[float, float, float, float] | None = ...,
        inset_ax: Axes | None = ...,
        *,
        transform: Transform | None = ...,
        facecolor: ColorType = ...,
        edgecolor: ColorType = ...,
        alpha: float = ...,
        zorder: float | None = ...,
        **kwargs
    ) -> InsetIndicator: ...
    def indicate_inset_zoom(self, inset_ax: Axes, **kwargs) -> InsetIndicator: ...
    def secondary_xaxis(
        self,
        location: Literal["top", "bottom"] | float,
        functions: tuple[
            Callable[[ArrayLike], ArrayLike], Callable[[ArrayLike], ArrayLike]
        ]
        | Transform
        | None = ...,
        *,
        transform: Transform | None = ...,
        **kwargs
    ) -> SecondaryAxis: ...
    def secondary_yaxis(
        self,
        location: Literal["left", "right"] | float,
        functions: tuple[
            Callable[[ArrayLike], ArrayLike], Callable[[ArrayLike], ArrayLike]
        ]
        | Transform
        | None = ...,
        *,
        transform: Transform | None = ...,
        **kwargs
    ) -> SecondaryAxis: ...
    def text(
        self,
        x: float,
        y: float,
        s: str,
        fontdict: dict[str, Any] | None = ...,
        **kwargs
    ) -> Text: ...
    def annotate(
        self,
        text: str,
        xy: tuple[Any, Any],
        xytext: tuple[Any, Any] | None = ...,
        xycoords: CoordsType = ...,
        textcoords: CoordsType | None = ...,
        arrowprops: dict[str, Any] | None = ...,
        annotation_clip: bool | None = ...,
        **kwargs
    ) -> Annotation: ...
    def axhline(
        self, y: float = ..., xmin: float = ..., xmax: float = ..., **kwargs
    ) -> Line2D: ...
    def axvline(
        self, x: float = ..., ymin: float = ..., ymax: float = ..., **kwargs
    ) -> Line2D: ...

    # TODO: Could separate the xy2 and slope signatures
    def axline(
        self,
        xy1: tuple[float, float],
        xy2: tuple[float, float] | None = ...,
        *,
        slope: float | None = ...,
        **kwargs
    ) -> AxLine: ...
    def axhspan(
        self, ymin: float, ymax: float, xmin: float = ..., xmax: float = ..., **kwargs
    ) -> Rectangle: ...
    def axvspan(
        self, xmin: float, xmax: float, ymin: float = ..., ymax: float = ..., **kwargs
    ) -> Rectangle: ...
    def hlines(
        self,
        y: float | ArrayLike,
        xmin: float | ArrayLike,
        xmax: float | ArrayLike,
        colors: ColorType | Sequence[ColorType] | None = ...,
        linestyles: LineStyleType = ...,
        *,
        label: str = ...,
        data: DataParamType = ...,
        **kwargs
    ) -> LineCollection: ...
    def vlines(
        self,
        x: float | ArrayLike,
        ymin: float | ArrayLike,
        ymax: float | ArrayLike,
        colors: ColorType | Sequence[ColorType] | None = ...,
        linestyles: LineStyleType = ...,
        *,
        label: str = ...,
        data: DataParamType = ...,
        **kwargs
    ) -> LineCollection: ...
    def eventplot(
        self,
        positions: ArrayLike | Sequence[ArrayLike],
        *,
        orientation: Literal["horizontal", "vertical"] = ...,
        lineoffsets: float | Sequence[float] = ...,
        linelengths: float | Sequence[float] = ...,
        linewidths: float | Sequence[float] | None = ...,
        colors: ColorType | Sequence[ColorType] | None = ...,
        alpha: float | Sequence[float] | None = ...,
        linestyles: LineStyleType | Sequence[LineStyleType] = ...,
        data: DataParamType = ...,
        **kwargs
    ) -> EventCollection: ...
    def plot(
        self,
        *args: float | ArrayLike | str,
        scalex: bool = ...,
        scaley: bool = ...,
        data: DataParamType = ...,
        **kwargs
    ) -> list[Line2D]: ...
    def loglog(self, *args, **kwargs) -> list[Line2D]: ...
    def semilogx(self, *args, **kwargs) -> list[Line2D]: ...
    def semilogy(self, *args, **kwargs) -> list[Line2D]: ...
    def acorr(
        self, x: ArrayLike, *, data: DataParamType = ..., **kwargs
    ) -> tuple[np.ndarray, np.ndarray, LineCollection | Line2D, Line2D | None]: ...
    def xcorr(
        self,
        x: ArrayLike,
        y: ArrayLike,
        *,
        normed: bool = ...,
        detrend: Callable[[ArrayLike], ArrayLike] = ...,
        usevlines: bool = ...,
        maxlags: int = ...,
        data: DataParamType = ...,
        **kwargs
    ) -> tuple[np.ndarray, np.ndarray, LineCollection | Line2D, Line2D | None]: ...
    def step(
        self,
        x: ArrayLike,
        y: ArrayLike,
        *args,
        where: Literal["pre", "post", "mid"] = ...,
        data: DataParamType = ...,
        **kwargs
    ) -> list[Line2D]: ...
    def bar(
        self,
        x: float | ArrayLike,
        height: float | ArrayLike,
        width: float | ArrayLike = ...,
        bottom: float | ArrayLike | None = ...,
        *,
        align: Literal["center", "edge"] = ...,
        data: DataParamType = ...,
        **kwargs
    ) -> BarContainer: ...
    def barh(
        self,
        y: float | ArrayLike,
        width: float | ArrayLike,
        height: float | ArrayLike = ...,
        left: float | ArrayLike | None = ...,
        *,
        align: Literal["center", "edge"] = ...,
        data: DataParamType = ...,
        **kwargs
    ) -> BarContainer: ...
    def bar_label(
        self,
        container: BarContainer,
        labels: ArrayLike | None = ...,
        *,
        fmt: str | Callable[[float], str] = ...,
        label_type: Literal["center", "edge"] = ...,
        padding: float | ArrayLike = ...,
        **kwargs
    ) -> list[Annotation]: ...
    def broken_barh(
        self,
        xranges: Sequence[tuple[float, float]],
        yrange: tuple[float, float],
        align: Literal["bottom", "center", "top"] = ...,
        *,
        data: DataParamType = ...,
        **kwargs
    ) -> PolyCollection: ...
    def grouped_bar(
        self,
        heights: Sequence[ArrayLike] | dict[str, ArrayLike] | np.ndarray | pd.DataFrame,
        *,
        positions: ArrayLike | None = ...,
        tick_labels: Sequence[str] | None = ...,
        labels: Sequence[str] | None = ...,
        group_spacing: float | None = ...,
        bar_spacing: float | None = ...,
        orientation: Literal["vertical", "horizontal"] = ...,
        colors: Iterable[ColorType] | None = ...,
        **kwargs
    ) -> list[BarContainer]: ...
    def stem(
        self,
        *args: ArrayLike | str,
        linefmt: str | None = ...,
        markerfmt: str | None = ...,
        basefmt: str | None = ...,
        bottom: float = ...,
        label: str | None = ...,
        orientation: Literal["vertical", "horizontal"] = ...,
        data: DataParamType = ...,
    ) -> StemContainer: ...

    # TODO: data kwarg preprocessor?
    def pie(
        self,
        x: ArrayLike,
        *,
        explode: ArrayLike | None = ...,
        labels: Sequence[str] | None = ...,
        colors: ColorType | Sequence[ColorType] | None = ...,
        autopct: str | Callable[[float], str] | None = ...,
        pctdistance: float = ...,
        shadow: bool = ...,
        labeldistance: float | None = ...,
        startangle: float = ...,
        radius: float = ...,
        counterclock: bool = ...,
        wedgeprops: dict[str, Any] | None = ...,
        textprops: dict[str, Any] | None = ...,
        center: tuple[float, float] = ...,
        frame: bool = ...,
        rotatelabels: bool = ...,
        normalize: bool = ...,
        hatch: str | Sequence[str] | None = ...,
        data: DataParamType = ...,
    ) -> PieContainer: ...
    def pie_label(
        self,
        container: PieContainer,
        /,
        labels: str | Sequence[str],
        *,
        distance: float = ...,
        textprops: dict | None = ...,
        rotate: bool = ...,
        alignment: str = ...,
    ) -> list[Text]: ...
    def errorbar(
        self,
        x: float | ArrayLike,
        y: float | ArrayLike,
        yerr: float | ArrayLike | None = ...,
        xerr: float | ArrayLike | None = ...,
        fmt: str = ...,
        *,
        ecolor: ColorType | None = ...,
        elinewidth: float | None = ...,
        elinestyle: LineStyleType | None = ...,
        capsize: float | None = ...,
        barsabove: bool = ...,
        lolims: bool | ArrayLike = ...,
        uplims: bool | ArrayLike = ...,
        xlolims: bool | ArrayLike = ...,
        xuplims: bool | ArrayLike = ...,
        errorevery: int | tuple[int, int] = ...,
        capthick: float | None = ...,
        data: DataParamType = ...,
        **kwargs
    ) -> ErrorbarContainer: ...
    def boxplot(
        self,
        x: ArrayLike | Sequence[ArrayLike],
        *,
        notch: bool | None = ...,
        sym: str | None = ...,
        vert: bool | None = ...,
        orientation: Literal["vertical", "horizontal"] = ...,
        whis: float | tuple[float, float] | None = ...,
        positions: ArrayLike | None = ...,
        widths: float | ArrayLike | None = ...,
        patch_artist: bool | None = ...,
        bootstrap: int | None = ...,
        usermedians: ArrayLike | None = ...,
        conf_intervals: ArrayLike | None = ...,
        meanline: bool | None = ...,
        showmeans: bool | None = ...,
        showcaps: bool | None = ...,
        showbox: bool | None = ...,
        showfliers: bool | None = ...,
        boxprops: dict[str, Any] | None = ...,
        tick_labels: Sequence[str] | None = ...,
        flierprops: dict[str, Any] | None = ...,
        medianprops: dict[str, Any] | None = ...,
        meanprops: dict[str, Any] | None = ...,
        capprops: dict[str, Any] | None = ...,
        whiskerprops: dict[str, Any] | None = ...,
        manage_ticks: bool = ...,
        autorange: bool = ...,
        zorder: float | None = ...,
        capwidths: float | ArrayLike | None = ...,
        label: Sequence[str] | None = ...,
        data: DataParamType = ...,
    ) -> dict[str, Any]: ...
    def bxp(
        self,
        bxpstats: Sequence[dict[str, Any]],
        positions: ArrayLike | None = ...,
        *,
        widths: float | ArrayLike | None = ...,
        vert: bool | None = ...,
        orientation: Literal["vertical", "horizontal"] = ...,
        patch_artist: bool = ...,
        shownotches: bool = ...,
        showmeans: bool = ...,
        showcaps: bool = ...,
        showbox: bool = ...,
        showfliers: bool = ...,
        boxprops: dict[str, Any] | None = ...,
        whiskerprops: dict[str, Any] | None = ...,
        flierprops: dict[str, Any] | None = ...,
        medianprops: dict[str, Any] | None = ...,
        capprops: dict[str, Any] | None = ...,
        meanprops: dict[str, Any] | None = ...,
        meanline: bool = ...,
        manage_ticks: bool = ...,
        zorder: float | None = ...,
        capwidths: float | ArrayLike | None = ...,
        label: Sequence[str] | None = ...,
    ) -> dict[str, Any]: ...
    def scatter(
        self,
        x: float | ArrayLike,
        y: float | ArrayLike,
        s: float | ArrayLike | None = ...,
        c: ArrayLike | Sequence[ColorType] | ColorType | None = ...,
        *,
        marker: MarkerType | None = ...,
        cmap: str | Colormap | None = ...,
        norm: str | Normalize | None = ...,
        vmin: float | None = ...,
        vmax: float | None = ...,
        alpha: float | None = ...,
        linewidths: float | Sequence[float] | None = ...,
        edgecolors: Literal["face", "none"] | ColorType | Sequence[ColorType] | None = ...,
        colorizer: Colorizer | None = ...,
        plotnonfinite: bool = ...,
        data: DataParamType = ...,
        **kwargs
    ) -> PathCollection: ...
    def hexbin(
        self,
        x: ArrayLike,
        y: ArrayLike,
        C: ArrayLike | None = ...,
        *,
        gridsize: int | tuple[int, int] = ...,
        bins: Literal["log"] | int | Sequence[float] | None = ...,
        xscale: Literal["linear", "log"] = ...,
        yscale: Literal["linear", "log"] = ...,
        extent: tuple[float, float, float, float] | None = ...,
        cmap: str | Colormap | None = ...,
        norm: str | Normalize | None = ...,
        vmin: float | None = ...,
        vmax: float | None = ...,
        alpha: float | None = ...,
        linewidths: float | None = ...,
        edgecolors: Literal["face", "none"] | ColorType = ...,
        reduce_C_function: Callable[[np.ndarray | list[float]], float] = ...,
        mincnt: int | None = ...,
        marginals: bool = ...,
        colorizer: Colorizer | None = ...,
        data: DataParamType = ...,
        **kwargs
    ) -> PolyCollection: ...
    def arrow(
        self, x: float, y: float, dx: float, dy: float, **kwargs
    ) -> FancyArrow: ...
    def quiverkey(
        self, Q: Quiver, X: float, Y: float, U: float, label: str, **kwargs
    ) -> QuiverKey: ...
    def quiver(self, *args, data: DataParamType = ..., **kwargs) -> Quiver: ...
    def barbs(self, *args, data: DataParamType = ..., **kwargs) -> Barbs: ...
    def fill(self, *args, data: DataParamType = ..., **kwargs) -> list[Polygon]: ...
    def fill_between(
        self,
        x: ArrayLike,
        y1: ArrayLike | float,
        y2: ArrayLike | float = ...,
        where: ArrayLike | None = ...,
        interpolate: bool = ...,
        step: Literal["pre", "post", "mid"] | None = ...,
        *,
        data: DataParamType = ...,
        **kwargs
    ) -> FillBetweenPolyCollection: ...
    def fill_betweenx(
        self,
        y: ArrayLike,
        x1: ArrayLike | float,
        x2: ArrayLike | float = ...,
        where: ArrayLike | None = ...,
        step: Literal["pre", "post", "mid"] | None = ...,
        interpolate: bool = ...,
        *,
        data: DataParamType = ...,
        **kwargs
    ) -> FillBetweenPolyCollection: ...
    def imshow(
        self,
        X: ArrayLike | PIL.Image.Image,
        cmap: str | Colormap | None = ...,
        norm: str | Normalize | None = ...,
        *,
        aspect: Literal["equal", "auto"] | float | None = ...,
        interpolation: str | None = ...,
        alpha: float | ArrayLike | None = ...,
        vmin: float | None = ...,
        vmax: float | None = ...,
        colorizer: Colorizer | None = ...,
        origin: Literal["upper", "lower"] | None = ...,
        extent: tuple[float, float, float, float] | None = ...,
        interpolation_stage: Literal["data", "rgba", "auto"] | None = ...,
        filternorm: bool = ...,
        filterrad: float = ...,
        resample: bool | None = ...,
        url: str | None = ...,
        data: DataParamType = ...,
        **kwargs
    ) -> AxesImage: ...
    def pcolor(
        self,
        *args: ArrayLike,
        shading: Literal["flat", "nearest", "auto"] | None = ...,
        alpha: float | None = ...,
        norm: str | Normalize | None = ...,
        cmap: str | Colormap | None = ...,
        vmin: float | None = ...,
        vmax: float | None = ...,
        colorizer: Colorizer | None = ...,
        data: DataParamType = ...,
        **kwargs
    ) -> Collection: ...
    def pcolormesh(
        self,
        *args: ArrayLike,
        alpha: float | None = ...,
        norm: str | Normalize | None = ...,
        cmap: str | Colormap | None = ...,
        vmin: float | None = ...,
        vmax: float | None = ...,
        colorizer: Colorizer | None = ...,
        shading: Literal["flat", "nearest", "gouraud", "auto"] | None = ...,
        antialiased: bool = ...,
        data: DataParamType = ...,
        **kwargs
    ) -> QuadMesh: ...
    def pcolorfast(
        self,
        *args: ArrayLike | tuple[float, float],
        alpha: float | None = ...,
        norm: str | Normalize | None = ...,
        cmap: str | Colormap | None = ...,
        vmin: float | None = ...,
        vmax: float | None = ...,
        colorizer: Colorizer | None = ...,
        data: DataParamType = ...,
        **kwargs
    ) -> AxesImage | PcolorImage | QuadMesh: ...
    def contour(self, *args, data: DataParamType = ..., **kwargs) -> QuadContourSet: ...
    def contourf(self, *args, data: DataParamType = ..., **kwargs) -> QuadContourSet: ...
    def clabel(
        self, CS: ContourSet, levels: ArrayLike | None = ..., **kwargs
    ) -> list[Text]: ...
    def hist(
        self,
        x: ArrayLike | Sequence[ArrayLike],
        bins: int | Sequence[float] | str | None = ...,
        *,
        range: tuple[float, float] | None = ...,
        density: bool = ...,
        weights: ArrayLike | None = ...,
        cumulative: bool | float = ...,
        bottom: ArrayLike | float | None = ...,
        histtype: Literal["bar", "barstacked", "step", "stepfilled"] = ...,
        align: Literal["left", "mid", "right"] = ...,
        orientation: Literal["vertical", "horizontal"] = ...,
        rwidth: float | None = ...,
        log: bool = ...,
        color: ColorType | Sequence[ColorType] | None = ...,
        label: str | Sequence[str] | None = ...,
        stacked: bool = ...,
        data: DataParamType = ...,
        **kwargs
    ) -> tuple[
        np.ndarray | list[np.ndarray],
        np.ndarray,
        BarContainer | Polygon | list[BarContainer | Polygon],
    ]: ...
    def stairs(
        self,
        values: ArrayLike,
        edges: ArrayLike | None = ...,
        *,
        orientation: Literal["vertical", "horizontal"] = ...,
        baseline: float | ArrayLike | None = ...,
        fill: bool = ...,
        data: DataParamType = ...,
        **kwargs
    ) -> StepPatch: ...
    def hist2d(
        self,
        x: ArrayLike,
        y: ArrayLike,
        bins: None
        | int
        | tuple[int, int]
        | ArrayLike
        | tuple[ArrayLike, ArrayLike] = ...,
        *,
        range: ArrayLike | None = ...,
        density: bool = ...,
        weights: ArrayLike | None = ...,
        cmin: float | None = ...,
        cmax: float | None = ...,
        data: DataParamType = ...,
        **kwargs
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, QuadMesh]: ...
    def ecdf(
        self,
        x: ArrayLike,
        weights: ArrayLike | None = ...,
        *,
        complementary: bool=...,
        orientation: Literal["vertical", "horizontal"]=...,
        compress: bool=...,
        data: DataParamType = ...,
        **kwargs
    ) -> Line2D: ...
    def psd(
        self,
        x: ArrayLike,
        *,
        NFFT: int | None = ...,
        Fs: float | None = ...,
        Fc: int | None = ...,
        detrend: Literal["none", "mean", "linear"]
        | Callable[[ArrayLike], ArrayLike]
        | None = ...,
        window: Callable[[ArrayLike], ArrayLike] | ArrayLike | None = ...,
        noverlap: int | None = ...,
        pad_to: int | None = ...,
        sides: Literal["default", "onesided", "twosided"] | None = ...,
        scale_by_freq: bool | None = ...,
        return_line: bool | None = ...,
        data: DataParamType = ...,
        **kwargs
    ) -> tuple[np.ndarray, np.ndarray] | tuple[np.ndarray, np.ndarray, Line2D]: ...
    def csd(
        self,
        x: ArrayLike,
        y: ArrayLike,
        *,
        NFFT: int | None = ...,
        Fs: float | None = ...,
        Fc: int | None = ...,
        detrend: Literal["none", "mean", "linear"]
        | Callable[[ArrayLike], ArrayLike]
        | None = ...,
        window: Callable[[ArrayLike], ArrayLike] | ArrayLike | None = ...,
        noverlap: int | None = ...,
        pad_to: int | None = ...,
        sides: Literal["default", "onesided", "twosided"] | None = ...,
        scale_by_freq: bool | None = ...,
        return_line: bool | None = ...,
        data: DataParamType = ...,
        **kwargs
    ) -> tuple[np.ndarray, np.ndarray] | tuple[np.ndarray, np.ndarray, Line2D]: ...
    def magnitude_spectrum(
        self,
        x: ArrayLike,
        *,
        Fs: float | None = ...,
        Fc: int | None = ...,
        window: Callable[[ArrayLike], ArrayLike] | ArrayLike | None = ...,
        pad_to: int | None = ...,
        sides: Literal["default", "onesided", "twosided"] | None = ...,
        scale: Literal["default", "linear", "dB"] | None = ...,
        data: DataParamType = ...,
        **kwargs
    ) -> tuple[np.ndarray, np.ndarray, Line2D]: ...
    def angle_spectrum(
        self,
        x: ArrayLike,
        *,
        Fs: float | None = ...,
        Fc: int | None = ...,
        window: Callable[[ArrayLike], ArrayLike] | ArrayLike | None = ...,
        pad_to: int | None = ...,
        sides: Literal["default", "onesided", "twosided"] | None = ...,
        data: DataParamType = ...,
        **kwargs
    ) -> tuple[np.ndarray, np.ndarray, Line2D]: ...
    def phase_spectrum(
        self,
        x: ArrayLike,
        *,
        Fs: float | None = ...,
        Fc: int | None = ...,
        window: Callable[[ArrayLike], ArrayLike] | ArrayLike | None = ...,
        pad_to: int | None = ...,
        sides: Literal["default", "onesided", "twosided"] | None = ...,
        data: DataParamType = ...,
        **kwargs
    ) -> tuple[np.ndarray, np.ndarray, Line2D]: ...
    def cohere(
        self,
        x: ArrayLike,
        y: ArrayLike,
        *,
        NFFT: int = ...,
        Fs: float = ...,
        Fc: int = ...,
        detrend: Literal["none", "mean", "linear"]
        | Callable[[ArrayLike], ArrayLike] = ...,
        window: Callable[[ArrayLike], ArrayLike] | ArrayLike = ...,
        noverlap: int = ...,
        pad_to: int | None = ...,
        sides: Literal["default", "onesided", "twosided"] = ...,
        scale_by_freq: bool | None = ...,
        data: DataParamType = ...,
        **kwargs
    ) -> tuple[np.ndarray, np.ndarray]: ...
    def specgram(
        self,
        x: ArrayLike,
        *,
        NFFT: int | None = ...,
        Fs: float | None = ...,
        Fc: int | None = ...,
        detrend: Literal["none", "mean", "linear"]
        | Callable[[ArrayLike], ArrayLike]
        | None = ...,
        window: Callable[[ArrayLike], ArrayLike] | ArrayLike | None = ...,
        noverlap: int | None = ...,
        cmap: str | Colormap | None = ...,
        xextent: tuple[float, float] | None = ...,
        pad_to: int | None = ...,
        sides: Literal["default", "onesided", "twosided"] | None = ...,
        scale_by_freq: bool | None = ...,
        mode: Literal["default", "psd", "magnitude", "angle", "phase"] | None = ...,
        scale: Literal["default", "linear", "dB"] | None = ...,
        vmin: float | None = ...,
        vmax: float | None = ...,
        data: DataParamType = ...,
        **kwargs
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, AxesImage]: ...
    def spy(
        self,
        Z: ArrayLike,
        *,
        precision: float | Literal["present"] = ...,
        marker: str | None = ...,
        markersize: float | None = ...,
        aspect: Literal["equal", "auto"] | float | None = ...,
        origin: Literal["upper", "lower"] = ...,
        **kwargs
    ) -> AxesImage: ...
    def matshow(self, Z: ArrayLike, **kwargs) -> AxesImage: ...
    def violinplot(
        self,
        dataset: ArrayLike | Sequence[ArrayLike],
        positions: ArrayLike | None = ...,
        *,
        vert: bool | None = ...,
        orientation: Literal["vertical", "horizontal"] = ...,
        widths: float | ArrayLike = ...,
        showmeans: bool = ...,
        showextrema: bool = ...,
        showmedians: bool = ...,
        quantiles: Sequence[float | Sequence[float]] | None = ...,
        points: int = ...,
        bw_method: Literal["scott", "silverman"]
        | float
        | Callable[[GaussianKDE], float]
        | None = ...,
        side: Literal["both", "low", "high"] = ...,
        facecolor: Sequence[ColorType] | ColorType | None = ...,
        linecolor: Sequence[ColorType] | ColorType | None = ...,
        data: DataParamType = ...,
    ) -> dict[str, Collection]: ...
    def violin(
        self,
        vpstats: Sequence[dict[str, Any]],
        positions: ArrayLike | None = ...,
        *,
        vert: bool | None = ...,
        orientation: Literal["vertical", "horizontal"] = ...,
        widths: float | ArrayLike = ...,
        showmeans: bool = ...,
        showextrema: bool = ...,
        showmedians: bool = ...,
        side: Literal["both", "low", "high"] = ...,
        facecolor: Sequence[ColorType] | ColorType | None = ...,
        linecolor: Sequence[ColorType] | ColorType | None = ...,
    ) -> dict[str, Collection]: ...

    table = mtable.table
    stackplot = mstack.stackplot
    streamplot = mstream.streamplot
    tricontour = mtri.tricontour
    tricontourf = mtri.tricontourf
    tripcolor = mtri.tripcolor
    triplot = mtri.triplot
