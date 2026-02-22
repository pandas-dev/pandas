from collections.abc import Callable
from typing import Concatenate, Literal, TypeAlias, TypeVar

import numpy as np
import optype.numpy as onp

__all__ = [
    "LineSearchWarning",
    "line_search_armijo",
    "line_search_wolfe1",
    "line_search_wolfe2",
    "scalar_search_wolfe1",
    "scalar_search_wolfe2",
]

_Float: TypeAlias = float | np.float64
_Float1D: TypeAlias = onp.Array1D[np.float64]

_RT = TypeVar("_RT")
_Fun1D: TypeAlias = Callable[Concatenate[_Float1D, ...], _RT]

###

class LineSearchWarning(RuntimeWarning): ...

def line_search_wolfe1(
    f: _Fun1D[onp.ToFloat],
    fprime: _Fun1D[onp.ToFloat1D],
    xk: onp.ToFloat1D,
    pk: onp.ToFloat1D,
    gfk: onp.ToFloat1D | None = None,
    old_fval: onp.ToFloat | None = None,
    old_old_fval: onp.ToFloat | None = None,
    args: tuple[object, ...] = (),
    c1: onp.ToFloat = 1e-4,
    c2: onp.ToFloat = 0.9,
    amax: onp.ToJustInt = 50,
    amin: onp.ToFloat = 1e-08,
    xtol: onp.ToFloat = 1e-14,
) -> tuple[_Float | None, int, int, _Float | None, _Float, _Float | None]: ...

# NOTE: exported as `scipy.optimize.line_search`
def line_search_wolfe2(
    f: _Fun1D[onp.ToFloat],
    myfprime: _Fun1D[onp.ToFloat1D],
    xk: onp.ToFloat1D,
    pk: onp.ToFloat1D,
    gfk: onp.ToFloat1D | None = None,
    old_fval: onp.ToFloat | None = None,
    old_old_fval: onp.ToFloat | None = None,
    args: tuple[object, ...] = (),
    c1: onp.ToFloat = 1e-4,
    c2: onp.ToFloat = 0.9,
    amax: onp.ToFloat | None = None,
    extra_condition: Callable[[float, _Float1D, float, _Float1D], onp.ToBool] | None = None,
    maxiter: onp.ToJustInt = 10,
) -> tuple[_Float | None, int, int, _Float | None, _Float, _Float | None]: ...

#
def line_search_armijo(
    f: _Fun1D[onp.ToFloat],
    xk: onp.ToFloat1D,
    pk: onp.ToFloat1D,
    gfk: onp.ToFloat1D,
    old_fval: onp.ToFloat,
    args: tuple[object, ...] = (),
    c1: onp.ToFloat = 1e-4,
    alpha0: onp.ToFloat = 1,
) -> tuple[_Float | None, int, _Float]: ...

# undocumented
def line_search_BFGS(
    f: _Fun1D[onp.ToFloat],
    xk: onp.ToFloat1D,
    pk: onp.ToFloat1D,
    gfk: onp.ToFloat1D,
    old_fval: onp.ToFloat,
    args: tuple[object, ...] = (),
    c1: onp.ToFloat = 1e-4,
    alpha0: onp.ToFloat = 1,
) -> tuple[_Float | None, int, Literal[0], _Float]: ...

#
def scalar_search_wolfe1(
    phi: Callable[[float], onp.ToFloat],
    derphi: Callable[[float], onp.ToFloat],
    phi0: onp.ToFloat | None = None,
    old_phi0: onp.ToFloat | None = None,
    derphi0: onp.ToFloat | None = None,
    c1: onp.ToFloat = 1e-4,
    c2: onp.ToFloat = 0.9,
    amax: onp.ToJustInt = 50,
    amin: onp.ToFloat = 1e-08,
    xtol: onp.ToFloat = 1e-14,
) -> tuple[_Float | None, _Float, _Float]: ...

#
def scalar_search_wolfe2(
    phi: Callable[[float], onp.ToFloat],
    derphi: Callable[[float], onp.ToFloat],
    phi0: onp.ToFloat | None = None,
    old_phi0: onp.ToFloat | None = None,
    derphi0: onp.ToFloat | None = None,
    c1: onp.ToFloat = 1e-4,
    c2: onp.ToFloat = 0.9,
    amax: onp.ToFloat | None = None,
    extra_condition: Callable[[float, float], onp.ToBool] | None = None,
    maxiter: onp.ToJustInt = 10,
) -> tuple[_Float | None, _Float, _Float, _Float | None]: ...

# undocumented
def scalar_search_armijo(
    phi: Callable[[float], onp.ToFloat],
    phi0: onp.ToFloat,
    derphi0: onp.ToFloat,
    c1: onp.ToFloat = 1e-4,
    alpha0: onp.ToFloat = 1,
    amin: onp.ToJustInt = 0,
) -> tuple[_Float | None, _Float]: ...
