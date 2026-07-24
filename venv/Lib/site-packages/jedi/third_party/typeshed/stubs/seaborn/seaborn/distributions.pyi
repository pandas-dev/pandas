from collections.abc import Iterable
from typing import Any, Literal, Protocol, TypeVar, type_check_only
from typing_extensions import TypeAlias, deprecated

from matplotlib.axes import Axes
from matplotlib.colors import Colormap
from matplotlib.typing import ColorType
from numpy.typing import ArrayLike

from ._core.typing import ColumnName, DataSource, NormSpec
from .axisgrid import FacetGrid
from .external.kde import _BwMethodType
from .utils import _DataSourceWideForm, _LogScale, _Palette, _Vector

__all__ = ["displot", "histplot", "kdeplot", "ecdfplot", "rugplot", "distplot"]

_T = TypeVar("_T")
_OneOrPair: TypeAlias = _T | tuple[_T, _T]

@type_check_only
class _Fit(Protocol):
    def fit(self, a: ArrayLike) -> tuple[ArrayLike, ...]: ...
    def pdf(self, x: ArrayLike, *params: ArrayLike) -> ArrayLike: ...

def histplot(
    data: DataSource | _DataSourceWideForm | None = None,
    *,
    x: ColumnName | _Vector | None = None,
    y: ColumnName | _Vector | None = None,
    hue: ColumnName | _Vector | None = None,
    weights: ColumnName | _Vector | None = None,
    stat: str = "count",
    bins: _OneOrPair[str | int | ArrayLike] = "auto",
    binwidth: float | tuple[float, float] | None = None,
    binrange: _OneOrPair[tuple[float, float]] | None = None,
    discrete: bool | None = None,
    cumulative: bool = False,
    common_bins: bool = True,
    common_norm: bool = True,
    multiple: Literal["layer", "dodge", "stack", "fill"] = "layer",
    element: Literal["bars", "step", "poly"] = "bars",
    fill: bool = True,
    shrink: float = 1,
    kde: bool = False,
    kde_kws: dict[str, Any] | None = None,
    line_kws: dict[str, Any] | None = None,
    thresh: float | None = 0,
    pthresh: float | None = None,
    pmax: float | None = None,
    cbar: bool = False,
    cbar_ax: Axes | None = None,
    cbar_kws: dict[str, Any] | None = None,
    palette: _Palette | Colormap | None = None,
    hue_order: Iterable[ColumnName] | None = None,
    hue_norm: NormSpec = None,
    color: ColorType | None = None,
    log_scale: _LogScale | None = None,
    legend: bool = True,
    ax: Axes | None = None,
    **kwargs: Any,
) -> Axes: ...
def kdeplot(
    data: DataSource | _DataSourceWideForm | None = None,
    *,
    x: ColumnName | _Vector | None = None,
    y: ColumnName | _Vector | None = None,
    hue: ColumnName | _Vector | None = None,
    weights: ColumnName | _Vector | None = None,
    palette: _Palette | Colormap | None = None,
    hue_order: Iterable[ColumnName] | None = None,
    hue_norm: NormSpec = None,
    color: ColorType | None = None,
    fill: bool | None = None,
    multiple: Literal["layer", "stack", "fill"] = "layer",
    common_norm: bool = True,
    common_grid: bool = False,
    cumulative: bool = False,
    bw_method: _BwMethodType = "scott",
    bw_adjust: float = 1,
    warn_singular: bool = True,
    log_scale: _LogScale | None = None,
    levels: int | Iterable[float] = 10,
    thresh: float = 0.05,
    gridsize: int = 200,
    cut: float = 3,
    clip: _OneOrPair[tuple[float | None, float | None]] | None = None,
    legend: bool = True,
    cbar: bool = False,
    cbar_ax: Axes | None = None,
    cbar_kws: dict[str, Any] | None = None,
    ax: Axes | None = None,
    **kwargs: Any,
) -> Axes: ...
def ecdfplot(
    data: DataSource | _DataSourceWideForm | None = None,
    *,
    x: ColumnName | _Vector | None = None,
    y: ColumnName | _Vector | None = None,
    hue: ColumnName | _Vector | None = None,
    weights: ColumnName | _Vector | None = None,
    stat: Literal["proportion", "percent", "count"] = "proportion",
    complementary: bool = False,
    palette: _Palette | Colormap | None = None,
    hue_order: Iterable[ColumnName] | None = None,
    hue_norm: NormSpec = None,
    log_scale: _LogScale | None = None,
    legend: bool = True,
    ax: Axes | None = None,
    **kwargs: Any,
) -> Axes: ...
def rugplot(
    data: DataSource | _DataSourceWideForm | None = None,
    *,
    x: ColumnName | _Vector | None = None,
    y: ColumnName | _Vector | None = None,
    hue: ColumnName | _Vector | None = None,
    height: float = 0.025,
    expand_margins: bool = True,
    palette: _Palette | Colormap | None = None,
    hue_order: Iterable[ColumnName] | None = None,
    hue_norm: NormSpec = None,
    legend: bool = True,
    ax: Axes | None = None,
    **kwargs: Any,
) -> Axes: ...
def displot(
    data: DataSource | _DataSourceWideForm | None = None,
    *,
    x: ColumnName | _Vector | None = None,
    y: ColumnName | _Vector | None = None,
    hue: ColumnName | _Vector | None = None,
    row: ColumnName | _Vector | None = None,
    col: ColumnName | _Vector | None = None,
    weights: ColumnName | _Vector | None = None,
    kind: Literal["hist", "kde", "ecdf"] = "hist",
    rug: bool = False,
    rug_kws: dict[str, Any] | None = None,
    log_scale: _LogScale | None = None,
    legend: bool = True,
    palette: _Palette | Colormap | None = None,
    hue_order: Iterable[ColumnName] | None = None,
    hue_norm: NormSpec = None,
    color: ColorType | None = None,
    col_wrap: int | None = None,
    row_order: Iterable[ColumnName] | None = None,
    col_order: Iterable[ColumnName] | None = None,
    height: float = 5,
    aspect: float = 1,
    facet_kws: dict[str, Any] | None = None,
    **kwargs: Any,
) -> FacetGrid: ...
@deprecated("Function `distplot` is deprecated and will be removed in seaborn v0.14.0")
def distplot(
    a: ArrayLike | None = None,
    bins: ArrayLike | None = None,
    hist: bool = True,
    kde: bool = True,
    rug: bool = False,
    fit: _Fit | None = None,
    hist_kws: dict[str, Any] | None = None,
    kde_kws: dict[str, Any] | None = None,
    rug_kws: dict[str, Any] | None = None,
    fit_kws: dict[str, Any] | None = None,
    color: ColorType | None = None,
    vertical: bool = False,
    norm_hist: bool = False,
    axlabel: str | Literal[False] | None = None,
    label: str | None = None,
    ax: Axes | None = None,
    x: ArrayLike | None = None,
) -> Axes: ...
