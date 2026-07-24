from collections.abc import Iterable, Mapping, Sequence
from typing import Any, Literal
from typing_extensions import TypeAlias

from matplotlib.axes import Axes
from matplotlib.colors import Colormap
from matplotlib.typing import MarkerType

from ._core.typing import ColumnName, DataSource, NormSpec
from .axisgrid import FacetGrid
from .utils import _DataSourceWideForm, _ErrorBar, _Estimator, _Legend, _Palette, _Seed, _Vector

__all__ = ["relplot", "scatterplot", "lineplot"]

_Sizes: TypeAlias = list[int] | list[float] | dict[str, int] | dict[str, float] | tuple[float, float]
_DashType: TypeAlias = tuple[None, None] | Sequence[float]  # See matplotlib.lines.Line2D.set_dashes
# "dashes" and "markers" require dict but we use mapping to avoid long unions because dict is invariant in its value type
_Dashes: TypeAlias = bool | Sequence[_DashType] | Mapping[Any, _DashType]
_Markers: TypeAlias = bool | Sequence[MarkerType] | Mapping[Any, MarkerType]

def lineplot(
    data: DataSource | _DataSourceWideForm | None = None,
    *,
    x: ColumnName | _Vector | None = None,
    y: ColumnName | _Vector | None = None,
    hue: ColumnName | _Vector | None = None,
    size: ColumnName | _Vector | None = None,
    style: ColumnName | _Vector | None = None,
    units: ColumnName | _Vector | None = None,
    weights: ColumnName | _Vector | None = None,
    palette: _Palette | Colormap | None = None,
    hue_order: Iterable[ColumnName] | None = None,
    hue_norm: NormSpec = None,
    sizes: _Sizes | None = None,
    size_order: Iterable[ColumnName] | None = None,
    size_norm: NormSpec = None,
    dashes: _Dashes | None = True,
    markers: _Markers | None = None,
    style_order: Iterable[ColumnName] | None = None,
    estimator: _Estimator | None = "mean",
    errorbar: _ErrorBar | None = ("ci", 95),
    n_boot: int = 1000,
    seed: _Seed | None = None,
    orient: Literal["x", "y"] = "x",
    sort: bool = True,
    err_style: Literal["band", "bars"] = "band",
    err_kws: dict[str, Any] | None = None,
    legend: _Legend = "auto",
    ci: str | int | None = "deprecated",  # deprecated
    ax: Axes | None = None,
    **kwargs: Any,
) -> Axes: ...
def scatterplot(
    data: DataSource | _DataSourceWideForm | None = None,
    *,
    x: ColumnName | _Vector | None = None,
    y: ColumnName | _Vector | None = None,
    hue: ColumnName | _Vector | None = None,
    size: ColumnName | _Vector | None = None,
    style: ColumnName | _Vector | None = None,
    palette: _Palette | Colormap | None = None,
    hue_order: Iterable[ColumnName] | None = None,
    hue_norm: NormSpec = None,
    sizes: _Sizes | None = None,
    size_order: Iterable[ColumnName] | None = None,
    size_norm: NormSpec = None,
    markers: _Markers | None = True,
    style_order: Iterable[ColumnName] | None = None,
    legend: _Legend = "auto",
    ax: Axes | None = None,
    **kwargs: Any,
) -> Axes: ...
def relplot(
    data: DataSource | _DataSourceWideForm | None = None,
    *,
    x: ColumnName | _Vector | None = None,
    y: ColumnName | _Vector | None = None,
    hue: ColumnName | _Vector | None = None,
    size: ColumnName | _Vector | None = None,
    style: ColumnName | _Vector | None = None,
    units: ColumnName | _Vector | None = None,
    weights: ColumnName | _Vector | None = None,
    row: ColumnName | _Vector | None = None,
    col: ColumnName | _Vector | None = None,
    col_wrap: int | None = None,
    row_order: Iterable[ColumnName] | None = None,
    col_order: Iterable[ColumnName] | None = None,
    palette: _Palette | Colormap | None = None,
    hue_order: Iterable[ColumnName] | None = None,
    hue_norm: NormSpec = None,
    sizes: _Sizes | None = None,
    size_order: Iterable[ColumnName] | None = None,
    size_norm: NormSpec = None,
    markers: _Markers | None = None,
    dashes: _Dashes | None = None,
    style_order: Iterable[ColumnName] | None = None,
    legend: _Legend = "auto",
    kind: Literal["scatter", "line"] = "scatter",
    height: float = 5,
    aspect: float = 1,
    facet_kws: dict[str, Any] | None = None,
    **kwargs: Any,
) -> FacetGrid: ...
