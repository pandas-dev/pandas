import inspect
import os
from _typeshed import Incomplete, SupportsKeysAndGetItem
from collections.abc import Callable, Generator
from contextlib import contextmanager
from typing import IO, Any, Literal, NoReturn, TypedDict, TypeVar
from typing_extensions import Never, Self

import matplotlib as mpl
from matplotlib.artist import Artist
from matplotlib.axes import Axes
from matplotlib.figure import Figure, SubFigure
from matplotlib.transforms import Bbox
from matplotlib.typing import ColorType
from seaborn._core.data import PlotData
from seaborn._core.moves import Move
from seaborn._core.scales import Scale
from seaborn._core.typing import DataSource, Default, OrderSpec, VariableSpec, VariableSpecList
from seaborn._marks.base import Mark
from seaborn._stats.base import Stat

_ClsT = TypeVar("_ClsT", bound=type[Any])

default: Default

class Layer(TypedDict, total=False):
    mark: Mark
    stat: Stat | None
    move: Move | list[Move] | None
    data: PlotData
    source: DataSource
    vars: dict[str, VariableSpec]
    orient: str
    legend: bool
    label: str | None

class FacetSpec(TypedDict, total=False):
    variables: dict[str, VariableSpec]
    structure: dict[str, list[str]]
    wrap: int | None

class PairSpec(TypedDict, total=False):
    variables: dict[str, VariableSpec]
    structure: dict[str, list[str]]
    cross: bool
    wrap: int | None

@contextmanager
def theme_context(params: dict[str, Any]) -> Generator[None]: ...
def build_plot_signature(cls: _ClsT) -> _ClsT: ...  # -> _ClsT & "__signature__ protocol"

class ThemeConfig(mpl.RcParams):
    THEME_GROUPS: list[str]
    def __init__(self) -> None: ...
    def reset(self) -> None: ...
    def update(self, other: SupportsKeysAndGetItem[Incomplete, Incomplete] | None = None, /, **kwds) -> None: ...  # type: ignore[override]

class DisplayConfig(TypedDict):
    format: Literal["png", "svg"]
    scaling: float
    hidpi: bool

class PlotConfig:
    def __init__(self) -> None: ...
    @property
    def theme(self) -> dict[str, Any]: ...
    @property
    def display(self) -> DisplayConfig: ...

@build_plot_signature
class Plot:
    __signature__: inspect.Signature
    config: PlotConfig
    def __init__(self, *args: DataSource | VariableSpec, data: DataSource = None, **variables: VariableSpec) -> None: ...
    def __add__(self, other: Never) -> NoReturn: ...
    def on(self, target: Axes | SubFigure | Figure) -> Plot: ...
    def add(
        self,
        mark: Mark,
        *transforms: Stat | Move,
        orient: str | None = None,
        legend: bool = True,
        label: str | None = None,
        data: DataSource = None,
        **variables: VariableSpec,
    ) -> Plot: ...
    def pair(
        self, x: VariableSpecList = None, y: VariableSpecList = None, wrap: int | None = None, cross: bool = True
    ) -> Plot: ...
    def facet(
        self,
        col: VariableSpec = None,
        row: VariableSpec = None,
        order: OrderSpec | dict[str, OrderSpec] = None,
        wrap: int | None = None,
    ) -> Plot: ...
    def scale(self, **scales: Scale) -> Plot: ...
    def share(self, **shares: bool | str) -> Plot: ...
    def limit(self, **limits: tuple[Any, Any]) -> Plot: ...
    def label(self, *, title: str | None = None, legend: str | None = None, **variables: str | Callable[[str], str]) -> Plot: ...
    def layout(
        self,
        *,
        size: tuple[float, float] | Default = ...,
        engine: str | None | Default = ...,
        extent: tuple[float, float, float, float] | Default = ...,
    ) -> Plot: ...
    def theme(self, config: dict[str, Any], /) -> Plot: ...
    # Same signature as Plotter.save
    def save(
        self,
        loc: str | os.PathLike[Any] | IO[Any],
        *,
        transparent: bool | None = None,
        dpi: float | None = 96,
        facecolor: ColorType | Literal["auto"] | None = ...,
        edgecolor: ColorType | Literal["auto"] | None = ...,
        orientation: str = "protrait",
        format: str | None = None,
        bbox_inches: Literal["tight"] | Bbox | None = None,
        pad_inches: float | None = None,
        bbox_extra_artists: list[Artist] | None = None,
        backend: str | None = None,
        **kwargs: Any,
    ) -> Self: ...
    # Same signature as Plotter.show
    def show(self, *, block: bool | None = None) -> None: ...
    def plot(self, pyplot: bool = False) -> Plotter: ...

class Plotter:
    def __init__(self, pyplot: bool, theme: dict[str, Any]) -> None: ...
    def save(
        self,
        loc: str | os.PathLike[Any] | IO[Any],
        *,
        # From matplotlib.figure.Figure.savefig
        transparent: bool | None = None,
        # keyword-only arguments below are the same as matplotlib.backend_bases.FigureCanvasBase.print_figure
        # but with different defaults
        dpi: float | None = 96,
        facecolor: ColorType | Literal["auto"] | None = ...,
        edgecolor: ColorType | Literal["auto"] | None = ...,
        orientation: str = "protrait",
        format: str | None = None,
        bbox_inches: Literal["tight"] | Bbox | None = None,
        pad_inches: float | None = None,
        bbox_extra_artists: list[Artist] | None = None,
        backend: str | None = None,
        # Further **kwargs can truly be anything from an overridden Canvas method that is still passed down
        **kwargs: Any,
    ) -> Self: ...
    # Same as matplotlib.backend_bases._Bases. No other backend override show
    def show(self, *, block: bool | None = None) -> None: ...
