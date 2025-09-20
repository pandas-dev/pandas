from collections.abc import Callable, Sequence
from typing import Concatenate, Literal, NamedTuple, TypeAlias, overload
from typing_extensions import TypeVar

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.stats.qmc import QMCEngine

__all__ = ["cumulative_simpson", "cumulative_trapezoid", "fixed_quad", "newton_cotes", "qmc_quad", "romb", "simpson", "trapezoid"]

_NDT_f = TypeVar("_NDT_f", bound=_QuadFuncOut)
_QuadFuncOut: TypeAlias = onp.ArrayND[npc.floating] | Sequence[float]

###

class QMCQuadResult(NamedTuple):
    integral: float
    standard_error: float

# sample-based integration
@overload
def trapezoid(
    y: onp.ToFloatND, x: onp.ToFloatND | None = None, dx: onp.ToFloat = 1.0, axis: op.CanIndex = -1
) -> npc.floating | onp.ArrayND[npc.floating]: ...
@overload
def trapezoid(
    y: onp.ToComplexND, x: onp.ToFloatND | None = None, dx: onp.ToFloat = 1.0, axis: op.CanIndex = -1
) -> npc.inexact | onp.ArrayND[npc.inexact]: ...

#
@overload
def simpson(
    y: onp.ToFloatND, x: onp.ToFloatND | None = None, *, dx: onp.ToFloat = 1.0, axis: op.CanIndex = -1
) -> npc.floating | onp.ArrayND[npc.floating]: ...
@overload
def simpson(
    y: onp.ToComplexND, x: onp.ToFloatND | None = None, *, dx: onp.ToFloat = 1.0, axis: op.CanIndex = -1
) -> npc.inexact | onp.ArrayND[npc.inexact]: ...

#
@overload
def romb(
    y: onp.ToFloatND, dx: onp.ToFloat = 1.0, axis: op.CanIndex = -1, show: bool = False
) -> npc.floating | onp.ArrayND[npc.floating]: ...
@overload
def romb(
    y: onp.ToComplexND, dx: onp.ToFloat = 1.0, axis: op.CanIndex = -1, show: bool = False
) -> npc.inexact | onp.ArrayND[npc.inexact]: ...

# sample-based cumulative integration
@overload
def cumulative_trapezoid(
    y: onp.ToFloatND,
    x: onp.ToFloatND | None = None,
    dx: onp.ToFloat = 1.0,
    axis: op.CanIndex = -1,
    initial: Literal[0] | None = None,
) -> onp.ArrayND[npc.floating]: ...
@overload
def cumulative_trapezoid(
    y: onp.ToComplexND,
    x: onp.ToFloatND | None = None,
    dx: onp.ToFloat = 1.0,
    axis: op.CanIndex = -1,
    initial: Literal[0] | None = None,
) -> onp.ArrayND[npc.inexact]: ...

#
@overload
def cumulative_simpson(
    y: onp.ToFloatND,
    *,
    x: onp.ToFloatND | None = None,
    dx: onp.ToFloat = 1.0,
    axis: op.CanIndex = -1,
    initial: onp.ToFloatND | None = None,
) -> onp.ArrayND[npc.floating]: ...
@overload
def cumulative_simpson(
    y: onp.ToComplexND,
    *,
    x: onp.ToFloatND | None = None,
    dx: onp.ToFloat = 1.0,
    axis: op.CanIndex = -1,
    initial: onp.ToComplexND | None = None,
) -> onp.ArrayND[npc.inexact]: ...

# function-based
@overload
def fixed_quad(
    func: Callable[[onp.Array1D[np.float64]], _NDT_f], a: onp.ToFloat, b: onp.ToFloat, args: tuple[()] = (), n: op.CanIndex = 5
) -> _NDT_f: ...
@overload
def fixed_quad(
    func: Callable[Concatenate[onp.Array1D[np.float64], ...], _NDT_f],
    a: onp.ToFloat,
    b: onp.ToFloat,
    args: tuple[object, ...],
    n: op.CanIndex = 5,
) -> _NDT_f: ...

#
def qmc_quad(
    func: Callable[[onp.Array2D[np.float64]], onp.ArrayND[npc.floating]],
    a: onp.ToFloat1D,
    b: onp.ToFloat1D,
    *,
    n_estimates: int = 8,
    n_points: int = 1024,
    qrng: QMCEngine | None = None,
    log: bool = False,
) -> QMCQuadResult: ...

# low-level
def newton_cotes(rn: int, equal: int = 0) -> tuple[onp.Array1D[np.float64], float]: ...
