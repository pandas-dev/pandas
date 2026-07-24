from _typeshed import Incomplete
from collections.abc import Callable, Iterable
from typing import Any, Literal, overload
from typing_extensions import TypeAlias

import pandas as pd
from matplotlib.axes import Axes
from matplotlib.typing import ColorType
from numpy.typing import NDArray

from .axisgrid import FacetGrid
from .utils import _Palette, _Seed

__all__ = ["lmplot", "regplot", "residplot"]

_Vector: TypeAlias = list[Incomplete] | pd.Series[Incomplete] | pd.Index[Incomplete] | NDArray[Incomplete]

def lmplot(
    data: pd.DataFrame,
    *,
    x: str | None = None,
    y: str | None = None,
    hue: str | None = None,
    col: str | None = None,
    row: str | None = None,
    palette: _Palette | None = None,
    col_wrap: int | None = None,
    height: float = 5,
    aspect: float = 1,
    markers: str = "o",
    sharex: bool | Literal["col", "row"] | None = None,  # deprecated
    sharey: bool | Literal["col", "row"] | None = None,  # deprecated
    hue_order: Iterable[str] | None = None,
    col_order: Iterable[str] | None = None,
    row_order: Iterable[str] | None = None,
    legend: bool = True,
    legend_out: bool | None = None,  # deprecated
    x_estimator: Callable[[Incomplete], Incomplete] | None = None,
    x_bins: int | _Vector | None = None,
    x_ci: Literal["ci", "sd"] | int | None = "ci",
    scatter: bool = True,
    fit_reg: bool = True,
    ci: int | None = 95,
    n_boot: int = 1000,
    units: str | None = None,
    seed: _Seed | None = None,
    order: int = 1,
    logistic: bool = False,
    lowess: bool = False,
    robust: bool = False,
    logx: bool = False,
    x_partial: str | None = None,
    y_partial: str | None = None,
    truncate: bool = True,
    x_jitter: float | None = None,
    y_jitter: float | None = None,
    scatter_kws: dict[str, Any] | None = None,
    line_kws: dict[str, Any] | None = None,
    facet_kws: dict[str, Any] | None = None,
) -> FacetGrid: ...
@overload
def regplot(
    data: None = None,
    *,
    x: _Vector | None = None,
    y: _Vector | None = None,
    x_estimator: Callable[[Incomplete], Incomplete] | None = None,
    x_bins: int | _Vector | None = None,
    x_ci: Literal["ci", "sd"] | int | None = "ci",
    scatter: bool = True,
    fit_reg: bool = True,
    ci: int | None = 95,
    n_boot: int = 1000,
    units: _Vector | None = None,
    seed: _Seed | None = None,
    order: int = 1,
    logistic: bool = False,
    lowess: bool = False,
    robust: bool = False,
    logx: bool = False,
    x_partial: _Vector | None = None,
    y_partial: _Vector | None = None,
    truncate: bool = True,
    dropna: bool = True,
    x_jitter: float | None = None,
    y_jitter: float | None = None,
    label: str | None = None,
    color: ColorType | None = None,
    marker: str = "o",
    scatter_kws: dict[str, Any] | None = None,
    line_kws: dict[str, Any] | None = None,
    ax: Axes | None = None,
) -> Axes: ...
@overload
def regplot(
    data: pd.DataFrame,
    *,
    x: str | _Vector | None = None,
    y: str | _Vector | None = None,
    x_estimator: Callable[[Incomplete], Incomplete] | None = None,
    x_bins: int | _Vector | None = None,
    x_ci: Literal["ci", "sd"] | int | None = "ci",
    scatter: bool = True,
    fit_reg: bool = True,
    ci: int | None = 95,
    n_boot: int = 1000,
    units: str | _Vector | None = None,
    seed: _Seed | None = None,
    order: int = 1,
    logistic: bool = False,
    lowess: bool = False,
    robust: bool = False,
    logx: bool = False,
    x_partial: str | _Vector | None = None,
    y_partial: str | _Vector | None = None,
    truncate: bool = True,
    dropna: bool = True,
    x_jitter: float | None = None,
    y_jitter: float | None = None,
    label: str | None = None,
    color: ColorType | None = None,
    marker: str = "o",
    scatter_kws: dict[str, Any] | None = None,
    line_kws: dict[str, Any] | None = None,
    ax: Axes | None = None,
) -> Axes: ...
@overload
def residplot(
    data: None = None,
    *,
    x: _Vector | None = None,
    y: _Vector | None = None,
    x_partial: _Vector | None = None,
    y_partial: _Vector | None = None,
    lowess: bool = False,
    order: int = 1,
    robust: bool = False,
    dropna: bool = True,
    label: str | None = None,
    color: ColorType | None = None,
    scatter_kws: dict[str, Any] | None = None,
    line_kws: dict[str, Any] | None = None,
    ax: Axes | None = None,
) -> Axes: ...
@overload
def residplot(
    data: pd.DataFrame,
    *,
    x: str | _Vector | None = None,
    y: str | _Vector | None = None,
    x_partial: str | _Vector | None = None,
    y_partial: str | _Vector | None = None,
    lowess: bool = False,
    order: int = 1,
    robust: bool = False,
    dropna: bool = True,
    label: str | None = None,
    color: ColorType | None = None,
    scatter_kws: dict[str, Any] | None = None,
    line_kws: dict[str, Any] | None = None,
    ax: Axes | None = None,
) -> Axes: ...
