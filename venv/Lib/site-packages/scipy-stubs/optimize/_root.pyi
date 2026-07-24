from _typeshed import Unused
from collections.abc import Callable, Mapping
from typing import Any, Concatenate, Generic, Literal, TypedDict, overload, type_check_only
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._nonlin import InverseJacobian
from ._optimize import OptimizeResult as _OptimizeResult
from scipy.sparse.linalg import LinearOperator

__all__ = ["root"]

###

_ScalarT_co = TypeVar("_ScalarT_co", bound=npc.inexact, default=np.float64 | Any, covariant=True)
_ShapeT_co = TypeVar("_ShapeT_co", bound=tuple[int, ...], default=tuple[Any, ...], covariant=True)

@type_check_only
class _RootOptionsHybr(TypedDict, total=False):
    col_deriv: bool
    xtol: float
    maxffev: int
    band: tuple[int, int]
    eps: float
    factor: float
    diag: onp.ToFloat1D

@type_check_only
class _RootOptionsLM(TypedDict, total=False):
    col_deriv: bool
    ftol: float
    xtol: float
    gtol: float
    maxiter: int
    eps: float
    factor: float
    diag: onp.ToFloat1D

@type_check_only
class _RootOptionsNonlin[JacOptionsT: Mapping[str, object]](TypedDict, total=False):
    nit: int
    disp: bool
    maxiter: int
    ftol: float
    fatol: float
    xtol: float
    xatol: float
    tol_norm: Callable[[onp.Array1D[np.float64]], float]
    line_search: Literal["armijo", "wolfe"] | None
    jac_options: JacOptionsT

@type_check_only
class _JacOptionsBase(TypedDict, total=False):
    alpha: float

@type_check_only
class _JacOptionsBroyden(_JacOptionsBase, total=False):
    reduction_method: Literal["restart", "simple", "svd"] | tuple[Literal["svd"], int]
    max_rank: int

@type_check_only
class _JacOptionsAnderson(_JacOptionsBase, total=False):
    M: float
    w0: float

@type_check_only
class _JacOptionsExcitingMixing(_JacOptionsBase, total=False):
    alphamax: float

@type_check_only
class _JacOptionsKrylov(TypedDict, total=False):
    rdiff: float
    method: Literal["lgmres", "gmres", "bicgstab", "cgs", "minres", "tfqmr"] | Callable[..., onp.ToFloatND]
    inner_M: LinearOperator | InverseJacobian
    inner_rtol: float
    inner_atol: float
    inner_maxiter: int
    outer_k: int

type _JacOptionsNonlin = _JacOptionsBroyden | _JacOptionsAnderson | _JacOptionsExcitingMixing | _JacOptionsKrylov

@type_check_only
class _RootOptionsDFSane(TypedDict, total=False):
    ftol: float
    fatol: float
    fnorm: Callable[[onp.ArrayND[np.float64]], float]
    maxfev: int
    disp: bool
    eta_strategy: Callable[[int, onp.ArrayND[np.float64], onp.ArrayND[np.float64]], onp.ToFloat1D]
    sigma_eps: float
    sigma_0: float
    M: int
    line_search: Literal["cruz", "cheng"]

type _MethodNonlin = Literal["broyden1", "broyden2", "anderson", "linearmixing", "diagbroyden", "excitingmixing", "krylov"]
type _Method = Literal["hybr", "lm", "df-sane", _MethodNonlin]

type _CallbackFn[ScalarT: npc.inexact, ShapeT: tuple[int, ...]] = Callable[
    [onp.ArrayND[ScalarT, ShapeT], onp.ArrayND[ScalarT, ShapeT]], Unused
]

type _ToFloatOrND = onp.ToFloat | onp.ToFloatND
type _ToComplexOrND = onp.ToComplex | onp.ToComplexND

###

class OptimizeResult(_OptimizeResult, Generic[_ScalarT_co, _ShapeT_co]):
    message: str
    success: bool
    fun: onp.ArrayND[_ScalarT_co, _ShapeT_co]
    x: onp.ArrayND[_ScalarT_co, _ShapeT_co]
    nfev: int
    method: _Method

@overload  # hybr | lm
def root(
    fun: Callable[Concatenate[onp.Array1D[np.float64], ...], _ToFloatOrND],
    x0: _ToFloatOrND,
    args: tuple[object, ...] = (),
    method: Literal["hybr", "lm"] = "hybr",
    jac: Callable[Concatenate[onp.Array1D[np.float64], ...], onp.ToFloat2D] | Literal[False] | None = None,
    tol: float | None = None,
    callback: _CallbackFn[np.float64, tuple[int]] | None = None,
    options: _RootOptionsHybr | _RootOptionsLM | None = None,
) -> OptimizeResult[np.float64, tuple[int]]: ...
@overload  # hybr | lm, jac=True
def root(
    fun: Callable[Concatenate[onp.Array1D[np.float64], ...], tuple[_ToFloatOrND, _ToFloatOrND]],
    x0: _ToFloatOrND,
    args: tuple[object, ...] = (),
    method: Literal["hybr", "lm"] = "hybr",
    *,
    jac: Literal[True],
    tol: float | None = None,
    callback: _CallbackFn[np.float64, tuple[int]] | None = None,
    options: _RootOptionsHybr | _RootOptionsLM | None = None,
) -> OptimizeResult[np.float64, tuple[int]]: ...
@overload  # df-sane, complex
def root[ScalarT: npc.inexact, ShapeT: tuple[int, ...]](
    fun: Callable[Concatenate[onp.ArrayND[ScalarT, ShapeT], ...], _ToComplexOrND],
    x0: _ToComplexOrND,
    args: tuple[object, ...] = (),
    *,
    method: Literal["df-sane"],
    jac: None = None,
    tol: float | None = None,
    callback: _CallbackFn[ScalarT, ShapeT] | None = None,
    options: _RootOptionsDFSane | None = None,
) -> OptimizeResult[ScalarT, ShapeT]: ...
@overload  # broyden1 | broyden2 | anderson | linearmixing | diagbroyden | excitingmixing | krylov
def root[ScalarT: npc.inexact, ShapeT: tuple[int, ...]](
    fun: Callable[Concatenate[onp.ArrayND[ScalarT, ShapeT], ...], _ToComplexOrND],
    x0: _ToComplexOrND,
    args: tuple[object, ...] = (),
    *,
    method: _MethodNonlin,
    jac: Callable[Concatenate[onp.ArrayND[ScalarT, ShapeT], ...], _ToComplexOrND] | Literal[False] | None = None,
    tol: float | None = None,
    callback: _CallbackFn[ScalarT, ShapeT] | None = None,
    options: _RootOptionsNonlin[_JacOptionsNonlin] | None = None,
) -> OptimizeResult[ScalarT, ShapeT]: ...
@overload  # broyden1 | broyden2 | anderson | linearmixing | diagbroyden | excitingmixing | krylov, jac=True
def root[ScalarT: npc.inexact, ShapeT: tuple[int, ...]](
    fun: Callable[Concatenate[onp.ArrayND[ScalarT, ShapeT], ...], tuple[_ToComplexOrND, _ToComplexOrND]],
    x0: _ToComplexOrND,
    args: tuple[object, ...] = (),
    *,
    method: _MethodNonlin,
    jac: Literal[True],
    tol: float | None = None,
    callback: _CallbackFn[ScalarT, ShapeT] | None = None,
    options: _RootOptionsNonlin[_JacOptionsNonlin] | None = None,
) -> OptimizeResult[ScalarT, ShapeT]: ...
