from collections.abc import Callable, Sequence
from typing import Literal, NamedTuple, TypeAlias

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["binned_statistic", "binned_statistic_2d", "binned_statistic_dd"]

_Statistic: TypeAlias = Literal["mean", "std", "median", "count", "sum", "min", "max"]

class BinnedStatisticResult(NamedTuple):
    statistic: onp.Array1D[npc.inexact]
    bin_edges: onp.Array1D[np.float64]
    binnumber: onp.Array1D[np.intp]

def binned_statistic(
    x: onp.ToComplex1D,
    values: onp.ToComplex1D | Sequence[onp.ToComplex1D],
    statistic: _Statistic | Callable[[onp.Array1D[np.float64]], onp.ToFloat] = "mean",
    bins: onp.ToInt | onp.ToFloat1D = 10,
    range: tuple[float, float] | Sequence[tuple[float, float]] | None = None,
) -> BinnedStatisticResult: ...

class BinnedStatistic2dResult(NamedTuple):
    statistic: onp.Array2D[npc.inexact]
    x_edge: onp.Array1D[np.float64]
    y_edge: onp.Array1D[np.float64]
    binnumber: onp.Array1D[np.intp]

def binned_statistic_2d(
    x: onp.ToComplex1D,
    y: onp.ToComplex1D,
    values: onp.ToComplex1D | Sequence[onp.ToComplex1D],
    statistic: _Statistic | Callable[[onp.ArrayND[np.float64]], onp.ToFloat] = "mean",
    bins: onp.ToInt | onp.ToFloat1D | Sequence[onp.ToFloat1D] = 10,
    range: tuple[int, int] | None = None,
    expand_binnumbers: bool = False,
) -> BinnedStatistic2dResult: ...

class BinnedStatisticddResult(NamedTuple):
    statistic: onp.ArrayND[npc.inexact]
    bin_edges: list[onp.Array1D[np.float64]]
    binnumber: onp.Array1D[np.intp] | onp.Array2D[np.intp]

def binned_statistic_dd(
    sample: onp.ToComplex2D,
    values: onp.ToComplex1D | Sequence[onp.ToComplex1D],
    statistic: _Statistic | Callable[[onp.ArrayND[np.float64]], onp.ToFloat] = "mean",
    bins: onp.ToInt | onp.ToFloat1D = 10,
    range: tuple[int, int] | None = None,
    expand_binnumbers: bool = False,
    binned_statistic_result: BinnedStatisticddResult | None = None,
) -> BinnedStatisticddResult: ...
