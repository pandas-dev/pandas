from _typeshed import Unused
from collections.abc import Callable
from typing import Literal, overload
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse._base import _spbase
from scipy.sparse.linalg import LinearOperator

__all__ = ["bicg", "bicgstab", "cg", "cgs", "gmres", "qmr"]

###

type _Float = np.float32 | np.float64
type _Complex = np.complex64 | np.complex128

type _ToInt = npc.integer | np.bool
type _ToFloat = _Float | _ToInt
type _ToComplex = _Complex | _ToFloat

type _ToLinearOperator[_ScalarT: npc.number | np.bool] = onp.CanArrayND[_ScalarT] | _spbase[_ScalarT] | LinearOperator[_ScalarT]

type _Callback[_ScalarT: npc.number | np.bool] = Callable[[onp.Array1D[_ScalarT]], Unused]

_FloatT = TypeVar("_FloatT", bound=_Float, default=np.float64)

###

@overload  # real
def bicg(
    A: _ToLinearOperator[_FloatT | _ToInt],
    b: onp.ToFloat1D,
    x0: onp.ToFloat1D | None = None,
    *,
    rtol: onp.ToFloat = 1e-5,
    atol: onp.ToFloat = 0.0,
    maxiter: int | None = None,
    M: _ToLinearOperator[_FloatT | _ToInt] | None = None,
    callback: _Callback[_FloatT] | None = None,
) -> tuple[onp.Array1D[_FloatT], int]: ...
@overload  # complex
def bicg[ComplexT: _Complex](
    A: _ToLinearOperator[ComplexT],
    b: onp.ToComplex1D,
    x0: onp.ToComplex1D | None = None,
    *,
    rtol: onp.ToFloat = 1e-5,
    atol: onp.ToFloat = 0.0,
    maxiter: int | None = None,
    M: _ToLinearOperator[ComplexT] | None = None,
    callback: _Callback[ComplexT] | None = None,
) -> tuple[onp.Array1D[ComplexT], int]: ...

#
@overload  # real
def bicgstab(
    A: _ToLinearOperator[_FloatT | _ToInt],
    b: onp.ToFloat1D,
    x0: onp.ToFloat1D | None = None,
    *,
    rtol: onp.ToFloat = 1e-5,
    atol: onp.ToFloat = 0.0,
    maxiter: int | None = None,
    M: _ToLinearOperator[_FloatT | _ToInt] | None = None,
    callback: _Callback[_FloatT] | None = None,
) -> tuple[onp.Array1D[_FloatT], int]: ...
@overload  # complex
def bicgstab[ComplexT: _Complex](
    A: _ToLinearOperator[ComplexT],
    b: onp.ToComplex1D,
    x0: onp.ToComplex1D | None = None,
    *,
    rtol: onp.ToFloat = 1e-5,
    atol: onp.ToFloat = 0.0,
    maxiter: int | None = None,
    M: _ToLinearOperator[ComplexT] | None = None,
    callback: _Callback[ComplexT] | None = None,
) -> tuple[onp.Array1D[ComplexT], int]: ...

#
@overload  # real
def cg(
    A: _ToLinearOperator[_FloatT | _ToInt],
    b: onp.ToFloat1D,
    x0: onp.ToFloat1D | None = None,
    *,
    rtol: onp.ToFloat = 1e-5,
    atol: onp.ToFloat = 0.0,
    maxiter: int | None = None,
    M: _ToLinearOperator[_FloatT | _ToInt] | None = None,
    callback: _Callback[_FloatT] | None = None,
) -> tuple[onp.Array1D[_FloatT], int]: ...
@overload  # complex
def cg[ComplexT: _Complex](
    A: _ToLinearOperator[ComplexT],
    b: onp.ToComplex1D,
    x0: onp.ToComplex1D | None = None,
    *,
    rtol: onp.ToFloat = 1e-5,
    atol: onp.ToFloat = 0.0,
    maxiter: int | None = None,
    M: _ToLinearOperator[ComplexT] | None = None,
    callback: _Callback[ComplexT] | None = None,
) -> tuple[onp.Array1D[ComplexT], int]: ...

#
def cgs(
    A: _ToLinearOperator[_FloatT | _ToInt],
    b: onp.ToFloat1D,
    x0: onp.ToFloat1D | None = None,
    *,
    rtol: onp.ToFloat = 1e-5,
    atol: onp.ToFloat = 0.0,
    maxiter: int | None = None,
    M: _ToLinearOperator[_FloatT | _ToInt] | None = None,
    callback: _Callback[_FloatT] | None = None,
) -> tuple[onp.Array1D[_FloatT], int]: ...

#
@overload  # real, callback_type: {"pr_norm", "legacy"} | None = ...
def gmres(
    A: _ToLinearOperator[_FloatT | _ToInt],
    b: onp.ToFloat1D,
    x0: onp.ToFloat1D | None = None,
    *,
    rtol: onp.ToFloat = 1e-5,
    atol: onp.ToFloat = 0.0,
    restart: int | None = None,
    maxiter: int | None = None,
    M: _ToLinearOperator[_ToFloat] | None = None,
    callback: Callable[[float], Unused] | Callable[[np.float64], Unused] | None = None,
    callback_type: Literal["pr_norm", "legacy"] | None = None,
) -> tuple[onp.Array1D[_FloatT], int]: ...
@overload  # real, callback_type: {"x"}
def gmres(
    A: _ToLinearOperator[_FloatT | _ToInt],
    b: onp.ToFloat1D,
    x0: onp.ToFloat1D | None = None,
    *,
    rtol: onp.ToFloat = 1e-5,
    atol: onp.ToFloat = 0.0,
    restart: int | None = None,
    maxiter: int | None = None,
    M: _ToLinearOperator[_ToFloat] | None = None,
    callback: _Callback[_FloatT] | None = None,
    callback_type: Literal["x"],
) -> tuple[onp.Array1D[_FloatT], int]: ...
@overload  # complex, callback_type: {"pr_norm", "legacy"} | None = ...
def gmres[ComplexT: _Complex](
    A: _ToLinearOperator[ComplexT],
    b: onp.ToComplex1D,
    x0: onp.ToComplex1D | None = None,
    *,
    rtol: onp.ToFloat = 1e-5,
    atol: onp.ToFloat = 0.0,
    restart: int | None = None,
    maxiter: int | None = None,
    M: _ToLinearOperator[_ToFloat] | None = None,
    callback: Callable[[float], Unused] | Callable[[np.float64], Unused] | None = None,
    callback_type: Literal["pr_norm", "legacy"] | None = None,
) -> tuple[onp.Array1D[ComplexT], int]: ...
@overload  # complex, callback_type: {"x"}
def gmres[ComplexT: _Complex](
    A: _ToLinearOperator[ComplexT],
    b: onp.ToComplex1D,
    x0: onp.ToComplex1D | None = None,
    *,
    rtol: onp.ToFloat = 1e-5,
    atol: onp.ToFloat = 0.0,
    restart: int | None = None,
    maxiter: int | None = None,
    M: _ToLinearOperator[_ToComplex] | None = None,
    callback: _Callback[ComplexT] | None = None,
    callback_type: Literal["x"],
) -> tuple[onp.Array1D[ComplexT], int]: ...

#
@overload  # real
def qmr(
    A: _ToLinearOperator[_FloatT | _ToInt],
    b: onp.ToFloat1D,
    x0: onp.ToFloat1D | None = None,
    *,
    rtol: onp.ToFloat = 1e-5,
    atol: onp.ToFloat = 0.0,
    maxiter: int | None = None,
    M1: _ToLinearOperator[_ToFloat] | None = None,
    M2: _ToLinearOperator[_ToFloat] | None = None,
    callback: _Callback[_FloatT] | None = None,
) -> tuple[onp.Array1D[_FloatT], int]: ...
@overload  # complex
def qmr[ComplexT: _Complex](
    A: _ToLinearOperator[ComplexT],
    b: onp.ToComplex1D,
    x0: onp.ToComplex1D | None = None,
    *,
    rtol: onp.ToFloat = 1e-5,
    atol: onp.ToFloat = 0.0,
    maxiter: int | None = None,
    M1: _ToLinearOperator[_ToComplex] | None = None,
    M2: _ToLinearOperator[_ToComplex] | None = None,
    callback: _Callback[ComplexT] | None = None,
) -> tuple[onp.Array1D[ComplexT], int]: ...
