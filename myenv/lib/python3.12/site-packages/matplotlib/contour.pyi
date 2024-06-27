import matplotlib.cm as cm
from matplotlib.artist import Artist
from matplotlib.axes import Axes
from matplotlib.collections import Collection, PathCollection
from matplotlib.colors import Colormap, Normalize
from matplotlib.font_manager import FontProperties
from matplotlib.path import Path
from matplotlib.patches import Patch
from matplotlib.text import Text
from matplotlib.transforms import Transform, TransformedPatchPath, TransformedPath
from matplotlib.ticker import Locator, Formatter

from numpy.typing import ArrayLike
import numpy as np
from collections.abc import Callable, Iterable, Sequence
from typing import Literal
from .typing import ColorType



class ContourLabeler:
    labelFmt: str | Formatter | Callable[[float], str] | dict[float, str]
    labelManual: bool | Iterable[tuple[float, float]]
    rightside_up: bool
    labelLevelList: list[float]
    labelIndiceList: list[int]
    labelMappable: cm.ScalarMappable
    labelCValueList: list[ColorType]
    labelXYs: list[tuple[float, float]]
    def clabel(
        self,
        levels: ArrayLike | None = ...,
        *,
        fontsize: str | float | None = ...,
        inline: bool = ...,
        inline_spacing: float = ...,
        fmt: str | Formatter | Callable[[float], str] | dict[float, str] | None = ...,
        colors: ColorType | Sequence[ColorType] | None = ...,
        use_clabeltext: bool = ...,
        manual: bool | Iterable[tuple[float, float]] = ...,
        rightside_up: bool = ...,
        zorder: float | None = ...
    ) -> list[Text]: ...
    def print_label(self, linecontour: ArrayLike, labelwidth: float) -> bool: ...
    def too_close(self, x: float, y: float, lw: float) -> bool: ...
    def get_text(
        self,
        lev: float,
        fmt: str | Formatter | Callable[[float], str] | dict[float, str],
    ) -> str: ...
    def locate_label(
        self, linecontour: ArrayLike, labelwidth: float
    ) -> tuple[float, float, float]: ...
    def calc_label_rot_and_inline(
        self,
        slc: ArrayLike,
        ind: int,
        lw: float,
        lc: ArrayLike | None = ...,
        spacing: int = ...,
    ) -> tuple[float, list[ArrayLike]]: ...
    def add_label(
        self, x: float, y: float, rotation: float, lev: float, cvalue: ColorType
    ) -> None: ...
    def add_label_clabeltext(
        self, x: float, y: float, rotation: float, lev: float, cvalue: ColorType
    ) -> None: ...
    def add_label_near(
        self,
        x: float,
        y: float,
        inline: bool = ...,
        inline_spacing: int = ...,
        transform: Transform | Literal[False] | None = ...,
    ) -> None: ...
    def pop_label(self, index: int = ...) -> None: ...
    def labels(self, inline: bool, inline_spacing: int) -> None: ...
    def remove(self) -> None: ...

class ContourSet(ContourLabeler, Collection):
    axes: Axes
    levels: Iterable[float]
    filled: bool
    linewidths: float | ArrayLike | None
    hatches: Iterable[str | None]
    origin: Literal["upper", "lower", "image"] | None
    extent: tuple[float, float, float, float] | None
    colors: ColorType | Sequence[ColorType]
    extend: Literal["neither", "both", "min", "max"]
    nchunk: int
    locator: Locator | None
    logscale: bool
    negative_linestyles: None | Literal[
        "solid", "dashed", "dashdot", "dotted"
    ] | Iterable[Literal["solid", "dashed", "dashdot", "dotted"]]
    clip_path: Patch | Path | TransformedPath | TransformedPatchPath | None
    labelTexts: list[Text]
    labelCValues: list[ColorType]
    @property
    def tcolors(self) -> list[tuple[tuple[float, float, float, float]]]: ...

    # only for not filled
    @property
    def tlinewidths(self) -> list[tuple[float]]: ...

    @property
    def allkinds(self) -> list[list[np.ndarray | None]]: ...
    @property
    def allsegs(self) -> list[list[np.ndarray]]: ...
    @property
    def alpha(self) -> float | None: ...
    @property
    def antialiased(self) -> bool: ...
    @antialiased.setter
    def antialiased(self, aa: bool | Sequence[bool]) -> None: ...
    @property
    def collections(self) -> list[PathCollection]: ...
    @property
    def linestyles(self) -> (
        None |
        Literal["solid", "dashed", "dashdot", "dotted"] |
        Iterable[Literal["solid", "dashed", "dashdot", "dotted"]]
    ): ...

    def __init__(
        self,
        ax: Axes,
        *args,
        levels: Iterable[float] | None = ...,
        filled: bool = ...,
        linewidths: float | ArrayLike | None = ...,
        linestyles: Literal["solid", "dashed", "dashdot", "dotted"]
        | Iterable[Literal["solid", "dashed", "dashdot", "dotted"]]
        | None = ...,
        hatches: Iterable[str | None] = ...,
        alpha: float | None = ...,
        origin: Literal["upper", "lower", "image"] | None = ...,
        extent: tuple[float, float, float, float] | None = ...,
        cmap: str | Colormap | None = ...,
        colors: ColorType | Sequence[ColorType] | None = ...,
        norm: str | Normalize | None = ...,
        vmin: float | None = ...,
        vmax: float | None = ...,
        extend: Literal["neither", "both", "min", "max"] = ...,
        antialiased: bool | None = ...,
        nchunk: int = ...,
        locator: Locator | None = ...,
        transform: Transform | None = ...,
        negative_linestyles: Literal["solid", "dashed", "dashdot", "dotted"]
        | Iterable[Literal["solid", "dashed", "dashdot", "dotted"]]
        | None = ...,
        clip_path: Patch | Path | TransformedPath | TransformedPatchPath | None = ...,
        **kwargs
    ) -> None: ...
    def legend_elements(
        self, variable_name: str = ..., str_format: Callable[[float], str] = ...
    ) -> tuple[list[Artist], list[str]]: ...
    def find_nearest_contour(
        self, x: float, y: float, indices: Iterable[int] | None = ..., pixel: bool = ...
    ) -> tuple[int, int, int, float, float, float]: ...

class QuadContourSet(ContourSet): ...
