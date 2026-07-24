from collections.abc import Callable, Iterable
from typing import Any, Concatenate, Final, Generic, Literal, Protocol, overload, type_check_only
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._linesearch import line_search_wolfe2 as line_search
from ._typing import Brack, MethodAll, Solver
from scipy._lib._util import _RichResult

__all__ = [
    "OptimizeResult",
    "OptimizeWarning",
    "approx_fprime",
    "bracket",
    "brent",
    "brute",
    "check_grad",
    "fmin",
    "fmin_bfgs",
    "fmin_cg",
    "fmin_ncg",
    "fmin_powell",
    "fminbound",
    "golden",
    "line_search",
    "rosen",
    "rosen_der",
    "rosen_hess",
    "rosen_hess_prod",
    "show_options",
]

###

type _Fn1[_XT, _YT] = Callable[Concatenate[_XT, ...], _YT]
type _Fn1_0d[_YT] = _Fn1[float, _YT] | _Fn1[np.float64, _YT]
type _Fn1_1d[_YT] = _Fn1[_Float1D, _YT]
type _Fn2[_XT, _PT, _YT] = Callable[Concatenate[_XT, _PT, ...], _YT]
type _Callback_1d = Callable[[_Float1D], None]

type _Int1D = onp.Array1D[np.intp]
type _Float = float | np.float64  # equivalent to `np.float64` in `numpy>=2.2`
type _Float1D = onp.Array1D[np.float64]
type _Float2D = onp.Array2D[np.float64]
type _ComplexCo1D = onp.Array1D[npc.number | np.bool]
type _FloatingND = onp.ArrayND[npc.floating]
type _FloatingCoND = onp.ArrayND[npc.floating | npc.integer | np.bool]
type _NumericND = onp.ArrayND[npc.number | np.bool | np.timedelta64 | np.object_]

type _Args = tuple[object, ...]
type _Brack = tuple[float, float] | tuple[float, float, float]
type _Disp = Literal[0, 1, 2, 3] | bool | np.bool
type _BracketInfo = tuple[
    _Float, _Float, _Float,  # xa, xb, xc
    _Float, _Float, _Float,  # fa, fb, fx
    int,  # funcalls
]  # fmt: skip
type _WarnFlag = Literal[0, 1, 2, 3, 4]
type _AllVecs = list[_Int1D | _Float1D]

_ResultValueT_co = TypeVar("_ResultValueT_co", default=Any, covariant=True)
_XT_contra = TypeVar("_XT_contra", bound=_ComplexCo1D, default=_Float1D, contravariant=True)
_ValueT_co = TypeVar("_ValueT_co", bound=float | npc.floating, default=_Float, covariant=True)
_JacT_co = TypeVar("_JacT_co", bound=onp.Array1D[npc.floating] | onp.Array2D[npc.floating], default=_Float1D, covariant=True)

@type_check_only
class _DoesFMin(Protocol):
    def __call__(self, func: _Fn1_1d[onp.ToFloat], x0: _Float1D, /, *, args: _Args) -> _FloatingND: ...

@type_check_only
class _DoesMap(Protocol):
    def __call__[VT, RT](self, func: Callable[[VT], RT], iterable: Iterable[VT], /) -> Iterable[RT]: ...

###

# NOTE: Unlike the docs suggest, `OptimizeResult` has no attributes by default:
#   For example, `RootResult` does not have any of the documented attributes,
#   even though it is a subclass of `OptimizeResult`
class OptimizeResult(_RichResult[_ResultValueT_co], Generic[_ResultValueT_co]): ...

#
class OptimizeWarning(UserWarning): ...
class BracketError(RuntimeError): ...  # undocumented

# undocumented
class MemoizeJac(Generic[_XT_contra, _ValueT_co, _JacT_co]):
    fun: _Fn1[_XT_contra, tuple[_ValueT_co, _JacT_co]]  # readonly
    x: _XT_contra  # readonly
    _value: _ValueT_co  # readonly
    jac: _JacT_co  # readonly

    def __init__(self, /, fun: _Fn1[_XT_contra, tuple[_ValueT_co, _JacT_co]]) -> None: ...
    def __call__(self, /, x: _XT_contra, *args: object) -> _ValueT_co: ...
    def derivative(self, /, x: _XT_contra, *args: object) -> _JacT_co: ...

# undocumented
class Brent(Generic[_ValueT_co]):
    _mintol: Final[float]  # 1e-11
    _cg: Final[float]  # 0.3819660

    func: _Fn1[np.float64, _ValueT_co]  # `float & np.float64`
    args: Final[_Args]
    tol: Final[float]
    maxiter: Final[int]

    xmin: _Float | None
    fval: _Float | None
    iter: int
    funcalls: int
    disp: _Disp
    brack: _Brack | None  # might be undefined; set by `set_bracket`

    def __init__(
        self,
        /,
        func: _Fn1_0d[onp.ToFloat],
        args: _Args = (),
        tol: onp.ToFloat = 1.48e-08,
        maxiter: int = 500,
        full_output: onp.ToBool = 0,  # ignored
        disp: _Disp = 0,
    ) -> None: ...
    def set_bracket(self, /, brack: _Brack | None = None) -> None: ...
    def get_bracket_info(self, /) -> _BracketInfo: ...
    def optimize(self, /) -> None: ...
    @overload
    def get_result(self, /, full_output: onp.ToFalse = False) -> _Float: ...
    @overload
    def get_result(self, /, full_output: onp.ToTrue) -> tuple[_Float, _Float, int, int]: ...  # xmin, fval, itere, funcalls

# undocumented
@overload
def is_finite_scalar(x: onp.ToScalar) -> np.bool: ...
@overload  # returns a `np.ndarray` of `size = 1`, but could have any `ndim`
def is_finite_scalar(x: _NumericND) -> Literal[False] | onp.Array[onp.AtLeast1D[Any], np.bool]: ...

# undocumented
@overload
def vecnorm(x: onp.ToComplex, ord: onp.ToFloat = 2) -> onp.ToFloat: ...
@overload
def vecnorm(x: _NumericND, ord: onp.ToInt = 2) -> _FloatingCoND: ...
@overload
def vecnorm(x: onp.ToFloatND, ord: onp.ToInt = 2) -> onp.ToFloat: ...
@overload
def vecnorm(x: onp.ToComplexND, ord: onp.ToFloat = 2) -> onp.ToFloat | _FloatingCoND: ...

# undocumented
def approx_fhess_p(
    x0: onp.ToFloat | onp.ToFloat1D,
    p: onp.ToFloat,
    fprime: _Fn1[_Float1D, _FloatingCoND],
    epsilon: onp.ToFloat | _FloatingCoND,  # scalar or 1d ndarray
    *args: object,
) -> _Float1D: ...

#
@overload  # full_output: False = ..., retall: False = ...
def fmin(
    func: _Fn1_1d[onp.ToFloat],
    x0: onp.ToFloat | onp.ToFloat1D,
    args: _Args = (),
    xtol: onp.ToFloat = 1e-4,
    ftol: onp.ToFloat = 1e-4,
    maxiter: int | None = None,
    maxfun: int | None = None,
    full_output: onp.ToFalse = 0,
    disp: _Disp = 1,
    retall: onp.ToFalse = 0,
    callback: _Callback_1d | None = None,
    initial_simplex: onp.ToFloat2D | None = None,
) -> _Float1D: ...
@overload  # full_output: False = ..., retall: True (keyword)
def fmin(
    func: _Fn1_1d[onp.ToFloat],
    x0: onp.ToFloat | onp.ToFloat1D,
    args: _Args = (),
    xtol: onp.ToFloat = 1e-4,
    ftol: onp.ToFloat = 1e-4,
    maxiter: int | None = None,
    maxfun: int | None = None,
    full_output: onp.ToFalse = 0,
    disp: _Disp = 1,
    *,
    retall: onp.ToTrue,
    callback: _Callback_1d | None = None,
    initial_simplex: onp.ToFloat2D | None = None,
) -> tuple[_Float1D, _AllVecs]: ...
@overload  # full_output: True  (keyword), retall: False = ...
def fmin(
    func: _Fn1_1d[onp.ToFloat],
    x0: onp.ToFloat | onp.ToFloat1D,
    args: _Args = (),
    xtol: onp.ToFloat = 1e-4,
    ftol: onp.ToFloat = 1e-4,
    maxiter: int | None = None,
    maxfun: int | None = None,
    *,
    full_output: onp.ToTrue,
    disp: _Disp = 1,
    retall: onp.ToFalse = 0,
    callback: _Callback_1d | None = None,
    initial_simplex: onp.ToFloat2D | None = None,
) -> tuple[_Float1D, onp.ToFloat, int, int, _WarnFlag]: ...  # x, fun, nit, nfev, status
@overload  # full_output: True  (keyword), retall: True
def fmin(
    func: _Fn1_1d[onp.ToFloat],
    x0: onp.ToFloat | onp.ToFloat1D,
    args: _Args = (),
    xtol: onp.ToFloat = 1e-4,
    ftol: onp.ToFloat = 1e-4,
    maxiter: int | None = None,
    maxfun: int | None = None,
    *,
    full_output: onp.ToTrue,
    disp: _Disp = 1,
    retall: onp.ToTrue,
    callback: _Callback_1d | None = None,
    initial_simplex: onp.ToFloat2D | None = None,
) -> tuple[_Float1D, onp.ToFloat, int, int, _WarnFlag, _AllVecs]: ...  # x, fun, nit, nfev, status, allvecs

#
@overload  # full_output: False = ..., retall: False = ...
def fmin_bfgs(
    f: _Fn1_1d[onp.ToFloat],
    x0: onp.ToFloat | onp.ToFloat1D,
    fprime: _Fn1_1d[_FloatingCoND] | None = None,
    args: _Args = (),
    gtol: onp.ToFloat = 1e-05,
    norm: onp.ToFloat = ...,  # inf
    epsilon: onp.ToFloat | onp.ToFloat1D = ...,
    maxiter: int | None = None,
    full_output: onp.ToFalse = 0,
    disp: _Disp = 1,
    retall: onp.ToFalse = 0,
    callback: _Callback_1d | None = None,
    xrtol: onp.ToFloat = 0,
    c1: onp.ToFloat = 1e-4,
    c2: onp.ToFloat = 0.9,
    hess_inv0: onp.ToFloat2D | None = None,
) -> _Float1D: ...
@overload  # full_output: False = ..., retall: True  (keyword)
def fmin_bfgs(
    f: _Fn1_1d[onp.ToFloat],
    x0: onp.ToFloat | onp.ToFloat1D,
    fprime: _Fn1_1d[_FloatingCoND] | None = None,
    args: _Args = (),
    gtol: onp.ToFloat = 1e-05,
    norm: onp.ToFloat = ...,  # inf
    epsilon: onp.ToFloat | onp.ToFloat1D = ...,
    maxiter: int | None = None,
    full_output: onp.ToFalse = 0,
    disp: _Disp = 1,
    *,
    retall: onp.ToTrue,
    callback: _Callback_1d | None = None,
    xrtol: onp.ToFloat = 0,
    c1: onp.ToFloat = 1e-4,
    c2: onp.ToFloat = 0.9,
    hess_inv0: onp.ToFloat2D | None = None,
) -> tuple[_Float1D, _AllVecs]: ...
@overload  # full_output: True (keyword), retall: False = ...
def fmin_bfgs(
    f: _Fn1_1d[onp.ToFloat],
    x0: onp.ToFloat | onp.ToFloat1D,
    fprime: _Fn1_1d[_FloatingCoND] | None = None,
    args: _Args = (),
    gtol: onp.ToFloat = 1e-05,
    norm: onp.ToFloat = ...,  # inf
    epsilon: onp.ToFloat | onp.ToFloat1D = ...,
    maxiter: int | None = None,
    *,
    full_output: onp.ToTrue,
    disp: _Disp = 1,
    retall: onp.ToFalse = 0,
    callback: _Callback_1d | None = None,
    xrtol: onp.ToFloat = 0,
    c1: onp.ToFloat = 1e-4,
    c2: onp.ToFloat = 0.9,
    hess_inv0: onp.ToFloat2D | None = None,
) -> tuple[_Float1D, _Float, _Float1D, _Float2D, int, int, _WarnFlag]: ...
@overload  # full_output: True (keyword), retall: True (keyword)
def fmin_bfgs(
    f: _Fn1_1d[onp.ToFloat],
    x0: onp.ToFloat | onp.ToFloat1D,
    fprime: _Fn1_1d[_FloatingCoND] | None = None,
    args: _Args = (),
    gtol: onp.ToFloat = 1e-05,
    norm: onp.ToFloat = ...,  # inf
    epsilon: onp.ToFloat | onp.ToFloat1D = ...,
    maxiter: int | None = None,
    *,
    full_output: onp.ToTrue,
    disp: _Disp = 1,
    retall: onp.ToTrue,
    callback: _Callback_1d | None = None,
    xrtol: onp.ToFloat = 0,
    c1: onp.ToFloat = 1e-4,
    c2: onp.ToFloat = 0.9,
    hess_inv0: onp.ToFloat2D | None = None,
) -> tuple[_Float1D, _Float, _Float1D, _Float2D, int, int, _WarnFlag, _AllVecs]: ...

#
@overload  # full_output: False = ..., retall: False = ...
def fmin_cg(
    f: _Fn1_1d[onp.ToFloat],
    x0: onp.ToFloat | onp.ToFloat1D,
    fprime: _Fn1_1d[_FloatingCoND] | None = None,
    args: _Args = (),
    gtol: onp.ToFloat = 1e-05,
    norm: onp.ToFloat = ...,  # inf
    epsilon: onp.ToFloat | _FloatingCoND = ...,
    maxiter: int | None = None,
    full_output: onp.ToFalse = 0,
    disp: _Disp = 1,
    retall: onp.ToFalse = 0,
    callback: _Callback_1d | None = None,
    c1: onp.ToFloat = 1e-4,
    c2: onp.ToFloat = 0.4,
) -> _Float1D: ...
@overload  # full_output: False = ..., retall: True  (keyword)
def fmin_cg(
    f: _Fn1_1d[onp.ToFloat],
    x0: onp.ToFloat | onp.ToFloat1D,
    fprime: _Fn1_1d[_FloatingCoND] | None = None,
    args: _Args = (),
    gtol: onp.ToFloat = 1e-05,
    norm: onp.ToFloat = ...,  # inf
    epsilon: onp.ToFloat | _FloatingCoND = ...,
    maxiter: int | None = None,
    full_output: onp.ToFalse = 0,
    disp: _Disp = 1,
    *,
    retall: onp.ToTrue,
    callback: _Callback_1d | None = None,
    c1: onp.ToFloat = 1e-4,
    c2: onp.ToFloat = 0.4,
) -> tuple[_Float1D, _AllVecs]: ...
@overload  # full_output: True (keyword), retall: False = ...
def fmin_cg(
    f: _Fn1_1d[onp.ToFloat],
    x0: onp.ToFloat | onp.ToFloat1D,
    fprime: _Fn1_1d[_FloatingCoND] | None = None,
    args: _Args = (),
    gtol: onp.ToFloat = 1e-05,
    norm: onp.ToFloat = ...,  # inf
    epsilon: onp.ToFloat | _FloatingCoND = ...,
    maxiter: int | None = None,
    *,
    full_output: onp.ToTrue,
    disp: _Disp = 1,
    retall: onp.ToFalse = 0,
    callback: _Callback_1d | None = None,
    c1: onp.ToFloat = 1e-4,
    c2: onp.ToFloat = 0.4,
) -> tuple[_Float1D, _Float, int, int, _WarnFlag]: ...
@overload  # full_output: True (keyword), retall: True (keyword)
def fmin_cg(
    f: _Fn1_1d[onp.ToFloat],
    x0: onp.ToFloat | onp.ToFloat1D,
    fprime: _Fn1_1d[_FloatingCoND] | None = None,
    args: _Args = (),
    gtol: onp.ToFloat = 1e-05,
    norm: onp.ToFloat = ...,  # inf
    epsilon: onp.ToFloat | _FloatingCoND = ...,
    maxiter: int | None = None,
    *,
    full_output: onp.ToTrue,
    disp: _Disp = 1,
    retall: onp.ToTrue,
    callback: _Callback_1d | None = None,
    c1: onp.ToFloat = 1e-4,
    c2: onp.ToFloat = 0.4,
) -> tuple[_Float1D, _Float, int, int, _WarnFlag, _AllVecs]: ...
@overload  # full_output: False = ..., retall: False = ...
def fmin_ncg(
    f: _Fn1_1d[onp.ToFloat],
    x0: onp.ToFloat | onp.ToFloat1D,
    fprime: _Fn1_1d[_FloatingCoND],
    fhess_p: _Fn2[_Float1D, _Float1D, onp.ToFloat] | None = None,
    fhess: _Fn1_1d[_FloatingCoND] | None = None,
    args: _Args = (),
    avextol: onp.ToFloat = 1e-5,
    epsilon: onp.ToFloat | _FloatingCoND = ...,
    maxiter: int | None = None,
    full_output: onp.ToFalse = 0,
    disp: _Disp = 1,
    retall: onp.ToFalse = 0,
    callback: _Callback_1d | None = None,
    c1: onp.ToFloat = 1e-4,
    c2: onp.ToFloat = 0.9,
) -> _Float1D: ...
@overload  # full_output: False = ..., retall: True  (keyword)
def fmin_ncg(
    f: _Fn1_1d[onp.ToFloat],
    x0: onp.ToFloat | onp.ToFloat1D,
    fprime: _Fn1_1d[_FloatingCoND],
    fhess_p: _Fn2[_Float1D, _Float1D, onp.ToFloat] | None = None,
    fhess: _Fn1_1d[_FloatingCoND] | None = None,
    args: _Args = (),
    avextol: onp.ToFloat = 1e-5,
    epsilon: onp.ToFloat | _FloatingCoND = ...,
    maxiter: int | None = None,
    full_output: onp.ToFalse = 0,
    disp: _Disp = 1,
    *,
    retall: onp.ToTrue,
    callback: _Callback_1d | None = None,
    c1: onp.ToFloat = 1e-4,
    c2: onp.ToFloat = 0.9,
) -> tuple[_Float1D, _AllVecs]: ...
@overload  # full_output: True (keyword), retall: False = ...
def fmin_ncg(
    f: _Fn1_1d[onp.ToFloat],
    x0: onp.ToFloat | onp.ToFloat1D,
    fprime: _Fn1_1d[_FloatingCoND],
    fhess_p: _Fn2[_Float1D, _Float1D, onp.ToFloat] | None = None,
    fhess: _Fn1_1d[_FloatingCoND] | None = None,
    args: _Args = (),
    avextol: onp.ToFloat = 1e-5,
    epsilon: onp.ToFloat | _FloatingCoND = ...,
    maxiter: int | None = None,
    *,
    full_output: onp.ToTrue,
    disp: _Disp = 1,
    retall: onp.ToFalse = 0,
    callback: _Callback_1d | None = None,
    c1: onp.ToFloat = 1e-4,
    c2: onp.ToFloat = 0.9,
) -> tuple[_Float1D, _Float, int, int, int, _WarnFlag]: ...
@overload  # full_output: True (keyword), retall: True (keyword)
def fmin_ncg(
    f: _Fn1_1d[onp.ToFloat],
    x0: onp.ToFloat | onp.ToFloat1D,
    fprime: _Fn1_1d[_FloatingCoND],
    fhess_p: _Fn2[_Float1D, _Float1D, onp.ToFloat] | None = None,
    fhess: _Fn1_1d[_FloatingCoND] | None = None,
    args: _Args = (),
    avextol: onp.ToFloat = 1e-5,
    epsilon: onp.ToFloat | _FloatingCoND = ...,
    maxiter: int | None = None,
    *,
    full_output: onp.ToTrue,
    disp: _Disp = 1,
    retall: onp.ToTrue,
    callback: _Callback_1d | None = None,
    c1: onp.ToFloat = 1e-4,
    c2: onp.ToFloat = 0.9,
) -> tuple[_Float1D, _Float, int, int, int, _WarnFlag, _AllVecs]: ...

#
@overload  # full_output: False = ..., retall: False = ...
def fmin_powell(
    func: _Fn1_1d[onp.ToFloat],
    x0: onp.ToFloat | onp.ToFloat1D,
    args: _Args = (),
    xtol: onp.ToFloat = 1e-4,
    ftol: onp.ToFloat = 1e-4,
    maxiter: int | None = None,
    maxfun: int | None = None,
    full_output: onp.ToFalse = 0,
    disp: _Disp = 1,
    retall: onp.ToFalse = 0,
    callback: _Callback_1d | None = None,
    direc: onp.ToFloat2D | None = None,
) -> _Float1D: ...
@overload  # full_output: False = ..., retall: True  (keyword)
def fmin_powell(
    func: _Fn1_1d[onp.ToFloat],
    x0: onp.ToFloat | onp.ToFloat1D,
    args: _Args = (),
    xtol: onp.ToFloat = 1e-4,
    ftol: onp.ToFloat = 1e-4,
    maxiter: int | None = None,
    maxfun: int | None = None,
    full_output: onp.ToFalse = 0,
    disp: _Disp = 1,
    *,
    retall: onp.ToTrue,
    callback: _Callback_1d | None = None,
    direc: onp.ToFloat2D | None = None,
) -> tuple[_Float1D, _AllVecs]: ...
@overload  # full_output: True (keyword), retall: False = ...
def fmin_powell(
    func: _Fn1_1d[onp.ToFloat],
    x0: onp.ToFloat | onp.ToFloat1D,
    args: _Args = (),
    xtol: onp.ToFloat = 1e-4,
    ftol: onp.ToFloat = 1e-4,
    maxiter: int | None = None,
    maxfun: int | None = None,
    *,
    full_output: onp.ToTrue,
    disp: _Disp = 1,
    retall: onp.ToFalse = 0,
    callback: _Callback_1d | None = None,
    direc: onp.ToFloat2D | None = None,
) -> tuple[_Float1D, _Float, _Float2D, int, int, _WarnFlag]: ...
@overload  # full_output: True (keyword), retall: True (keyword)
def fmin_powell(
    func: _Fn1_1d[onp.ToFloat],
    x0: onp.ToFloat | onp.ToFloat1D,
    args: _Args = (),
    xtol: onp.ToFloat = 1e-4,
    ftol: onp.ToFloat = 1e-4,
    maxiter: int | None = None,
    maxfun: int | None = None,
    *,
    full_output: onp.ToTrue,
    disp: _Disp = 1,
    retall: onp.ToTrue,
    callback: _Callback_1d | None = None,
    direc: onp.ToFloat2D | None = None,
) -> tuple[_Float1D, _Float, _Float2D, int, int, _WarnFlag, _AllVecs]: ...

#
@overload  # full_output: False = ...
def fminbound(
    func: _Fn1_0d[onp.ToFloat],
    x1: onp.ToFloat,
    x2: onp.ToFloat,
    args: _Args = (),
    xtol: onp.ToFloat = 1e-05,
    maxfun: int = 500,
    full_output: onp.ToFalse = 0,
    disp: _Disp = 1,
) -> _Float: ...
@overload  # full_output: True (keyword)
def fminbound(
    func: _Fn1_0d[onp.ToFloat],
    x1: onp.ToFloat,
    x2: onp.ToFloat,
    args: _Args = (),
    xtol: onp.ToFloat = 1e-05,
    maxfun: int = 500,
    *,
    full_output: onp.ToTrue,
    disp: _Disp = 1,
) -> tuple[_Float, _Float, _WarnFlag, int]: ...  # x, fun, status, nfev

#
@overload  # full_output: False = ...
def brute(
    func: _Fn1_1d[onp.ToFloat],
    ranges: tuple[tuple[onp.ToFloat, onp.ToFloat] | slice, ...],
    args: _Args = (),
    Ns: int = 20,
    full_output: onp.ToFalse = 0,
    finish: _DoesFMin | None = ...,  # default: `fmin`
    disp: bool = False,
    workers: int | _DoesMap = 1,
) -> _Float1D: ...
@overload  # full_output: True (keyword)
def brute(
    func: _Fn1_1d[onp.ToFloat],
    ranges: tuple[tuple[onp.ToFloat, onp.ToFloat] | slice, ...],
    args: _Args = (),
    Ns: int = 20,
    *,
    full_output: onp.ToTrue,
    finish: _DoesFMin | None = ...,  # default: `fmin`
    disp: bool = False,
    workers: int | _DoesMap = 1,
) -> tuple[_Float1D, np.float64, onp.Array3D[np.float64], onp.Array2D[npc.floating]]: ...

#
@overload  # full_output: False = ...
def brent(
    func: _Fn1_0d[onp.ToFloat],
    args: _Args = (),
    brack: Brack | None = None,
    tol: onp.ToFloat = 1.48e-08,
    full_output: onp.ToFalse = 0,
    maxiter: int = 500,
) -> _Float: ...
@overload  # full_output: True (positional)
def brent(
    func: _Fn1_0d[onp.ToFloat], args: _Args, brack: Brack | None, tol: onp.ToFloat, full_output: onp.ToTrue, maxiter: int = 500
) -> tuple[_Float, _Float, int, int]: ...
@overload  # full_output: True (keyword)
def brent(
    func: _Fn1_0d[onp.ToFloat],
    args: _Args = (),
    brack: Brack | None = None,
    tol: onp.ToFloat = 1.48e-08,
    *,
    full_output: onp.ToTrue,
    maxiter: int = 500,
) -> tuple[_Float, _Float, int, int]: ...

#
@overload  # full_output: False = ...
def golden(
    func: _Fn1_0d[onp.ToFloat],
    args: _Args = (),
    brack: Brack | None = None,
    tol: onp.ToFloat = ...,
    full_output: onp.ToFalse = 0,
    maxiter: int = 5_000,
) -> _Float: ...
@overload  # full_output: True (positional)
def golden(
    func: _Fn1_0d[onp.ToFloat], args: _Args, brack: Brack | None, tol: onp.ToFloat, full_output: onp.ToTrue, maxiter: int = 5_000
) -> tuple[_Float, _Float, int]: ...
@overload  # full_output: True (keyword)
def golden(
    func: _Fn1_0d[onp.ToFloat],
    args: _Args = (),
    brack: Brack | None = None,
    tol: onp.ToFloat = ...,
    *,
    full_output: onp.ToTrue,
    maxiter: int = 5_000,
) -> tuple[_Float, _Float, int]: ...

#
def bracket(
    func: _Fn1_0d[onp.ToFloat],
    xa: onp.ToFloat = 0.0,
    xb: onp.ToFloat = 1.0,
    args: _Args = (),
    grow_limit: onp.ToFloat = 110.0,
    maxiter: int = 1_000,
) -> _BracketInfo: ...

# rosenbrock

@overload
def rosen[ScalarT: npc.inexact](x: onp.Array1D[ScalarT]) -> ScalarT: ...
@overload
def rosen[ScalarT: npc.inexact](x: onp.Array2D[ScalarT]) -> onp.Array1D[ScalarT]: ...
@overload
def rosen[ScalarT: npc.inexact](x: onp.Array3D[ScalarT]) -> onp.Array2D[ScalarT]: ...
@overload
def rosen[ScalarT: npc.inexact](x: onp.ArrayND[ScalarT]) -> onp.ArrayND[ScalarT] | ScalarT: ...
@overload
def rosen(x: onp.ToArrayStrict1D[float, npc.integer | np.bool]) -> np.float64: ...
@overload
def rosen(x: onp.ToArrayStrict2D[float, npc.integer | np.bool]) -> onp.Array1D[np.float64]: ...
@overload
def rosen(x: onp.ToArrayStrict3D[float, npc.integer | np.bool]) -> onp.Array2D[np.float64]: ...
@overload
def rosen(x: onp.ToArrayND[float, npc.integer | np.bool]) -> onp.ArrayND[np.float64] | np.float64: ...

#
@overload
def rosen_der[ShapeT: onp.AtLeast1D, ScalarT: npc.inexact](x: onp.Array[ShapeT, ScalarT]) -> onp.Array[ShapeT, ScalarT]: ...
@overload
def rosen_der[ShapeT: onp.AtLeast1D](x: onp.Array[ShapeT, npc.integer | np.bool]) -> onp.ArrayND[np.float64, ShapeT]: ...
@overload  # the weird shape-type avoids a pyright overload overlap error on `numpy<2.1`
def rosen_der(x: onp.ToArrayND[float, npc.integer | np.bool]) -> onp.ArrayND[np.float64, tuple[Any, ...] | tuple[int]]: ...

#
@overload
def rosen_hess[ScalarT: npc.inexact](x: onp.Array1D[ScalarT]) -> onp.Array2D[ScalarT]: ...
@overload
def rosen_hess(x: onp.ToArray1D[float, npc.integer | np.bool]) -> onp.Array2D[np.float64]: ...

#
@overload
def rosen_hess_prod[ScalarT: npc.inexact](x: onp.Array1D[ScalarT], p: onp.ToComplex1D) -> onp.Array1D[ScalarT]: ...
@overload
def rosen_hess_prod(x: onp.ToArray1D[float, npc.integer | np.bool], p: onp.ToFloat1D) -> onp.Array1D[np.float64]: ...

#
@overload  # disp: True = ...
def show_options(solver: Solver | None = None, method: MethodAll | None = None, disp: onp.ToTrue = True) -> None: ...
@overload  # disp: False  (positional)
def show_options(solver: Solver | None, method: MethodAll | None, disp: onp.ToFalse) -> str: ...
@overload  # disp: False  (keyword)
def show_options(solver: Solver | None = None, method: MethodAll | None = None, *, disp: onp.ToFalse) -> str: ...

#
@overload
def approx_fprime(
    xk: onp.ToFloat1D, f: _Fn1_1d[_Float], epsilon: onp.ToFloat | _FloatingCoND = ..., *args: object
) -> _Float1D: ...
@overload
def approx_fprime(
    xk: onp.ToFloat1D, f: _Fn1_1d[_Float1D], epsilon: onp.ToFloat | _FloatingCoND = ..., *args: object
) -> _Float2D: ...

#
def check_grad(
    func: _Fn1_1d[onp.ToFloat],
    grad: _Fn1_1d[_FloatingCoND],
    x0: onp.ToFloat | onp.ToFloat1D,
    *args: object,
    epsilon: onp.ToFloat = ...,
    direction: Literal["all", "random"] = "all",
    rng: onp.random.ToRNG | None = None,
    seed: onp.random.ToRNG | None = None,
) -> _Float: ...
