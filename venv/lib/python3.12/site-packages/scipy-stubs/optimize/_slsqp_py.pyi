from collections.abc import Callable, Sequence
from typing import Any, Concatenate, Final, Literal, TypeAlias, TypeVar, TypedDict, overload, type_check_only

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["approx_jacobian", "fmin_slsqp"]

###

_FT = TypeVar("_FT", bound=onp.ToFloat | onp.ToFloatND)
_Fun: TypeAlias = Callable[Concatenate[onp.Array1D[np.float64], ...], _FT]
_Fun0D: TypeAlias = _Fun[onp.ToFloat]
_Fun1D: TypeAlias = _Fun[onp.ToFloat1D]
_Fun2D: TypeAlias = _Fun[onp.ToFloat2D]

_Ignored: TypeAlias = object

_ExitMode: TypeAlias = Literal[-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
_ExitDesc: TypeAlias = Literal[
    "Gradient evaluation required (g & a)",  # -1
    "Optimization terminated successfully",  # 0
    "Function evaluation required (f & c)",  # 1
    "More equality constraints than independent variables",  # 2
    "More than 3*n iterations in LSQ subproblem",  # 3
    "Inequality constraints incompatible",  # 4
    "Singular matrix E in LSQ subproblem",  # 5
    "Singular matrix C in LSQ subproblem",  # 6
    "Rank-deficient equality constraint subproblem HFTI",  # 7
    "Positive directional derivative for linesearch",  # 8
    "Iteration limit reached",  # 9
]

@type_check_only
class _ConDict(TypedDict):
    fun: Callable[..., onp.ToFloat]
    jac: Callable[..., onp.ToFloat2D] | None
    args: tuple[Any, ...]

@type_check_only
class _ConsDict(TypedDict):
    eq: tuple[_ConDict, ...]
    ineq: tuple[_ConDict, ...]

###

__docformat__: Final = "restructuredtext en"  # undocumented

# private
def approx_jacobian(x: onp.ToFloat1D, func: _Fun0D, epsilon: onp.ToFloat, *args: object) -> onp.Array2D[np.float64]: ...

# public
@overload
def fmin_slsqp(
    func: _Fun0D,
    x0: onp.ToFloat1D,
    eqcons: Sequence[_Fun0D] = (),
    f_eqcons: _Fun1D | None = None,
    ieqcons: Sequence[_Fun0D] = (),
    f_ieqcons: _Fun1D | None = None,
    bounds: Sequence[tuple[onp.ToFloat, onp.ToFloat]] = (),
    fprime: _Fun1D | None = None,
    fprime_eqcons: _Fun2D | None = None,
    fprime_ieqcons: _Fun2D | None = None,
    args: Sequence[object] = (),
    iter: onp.ToJustInt = 100,
    acc: onp.ToFloat = 1e-06,
    iprint: onp.ToJustInt = 1,
    disp: onp.ToInt | None = None,
    full_output: onp.ToFalse = 0,
    epsilon: onp.ToFloat = ...,  # = np.sqrt(np.finfo(float).eps)
    callback: Callable[[onp.Array1D[np.float64]], _Ignored] | None = None,
) -> onp.Array1D[np.float64]: ...
@overload
def fmin_slsqp(
    func: _Fun0D,
    x0: onp.ToFloat1D,
    eqcons: Sequence[_Fun0D] = (),
    f_eqcons: _Fun1D | None = None,
    ieqcons: Sequence[_Fun0D] = (),
    f_ieqcons: _Fun1D | None = None,
    bounds: Sequence[tuple[onp.ToFloat, onp.ToFloat]] = (),
    fprime: _Fun1D | None = None,
    fprime_eqcons: _Fun2D | None = None,
    fprime_ieqcons: _Fun2D | None = None,
    args: Sequence[object] = (),
    iter: onp.ToJustInt = 100,
    acc: onp.ToFloat = 1e-06,
    iprint: onp.ToJustInt = 1,
    disp: onp.ToInt | None = None,
    *,
    full_output: onp.ToTrue,
    epsilon: onp.ToFloat = ...,  # = np.sqrt(np.finfo(float).eps)
    callback: Callable[[onp.Array1D[np.float64]], _Ignored] | None = None,
) -> tuple[onp.Array1D[np.float64], float | np.float64, int, _ExitMode, _ExitDesc]: ...

#
def _eval_constraint(d: onp.Array1D[np.float64], x: onp.Array1D[npc.floating], cons: _ConsDict, m: int, meq: int) -> None: ...
def _eval_con_normals(C: onp.Array2D[np.float64], x: onp.Array1D[npc.floating], cons: _ConsDict, m: int, meq: int) -> None: ...
