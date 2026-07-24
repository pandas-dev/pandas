from _typeshed import Unused
from collections.abc import Callable
from typing import Any, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse._base import _spbase
from scipy.sparse.linalg import LinearOperator

__all__ = ["minres"]

###

type _ToLinearOperator[ScalarT: npc.number | np.bool] = onp.CanArrayND[ScalarT] | _spbase[ScalarT] | LinearOperator[ScalarT]

###

# mypy false positive (I'm guessing related to the unions in the next overload)
# mypy: disable-error-code=overload-overlap

@overload  # f64
def minres(
    A: _ToLinearOperator[npc.floating64 | npc.integer],
    b: onp.ToFloat1D,
    x0: onp.ToFloat1D | None = None,
    *,
    rtol: float = 1e-5,
    shift: float = 0.0,
    maxiter: int | None = None,
    M: _ToLinearOperator[npc.floating | npc.integer] | None = None,
    callback: Callable[[onp.Array1D[np.float64]], Unused] | None = None,
    show: bool = False,
    check: bool = False,
) -> tuple[onp.Array1D[np.float64], int]: ...
@overload  # f32
def minres(
    A: _ToLinearOperator[np.float32],
    b: onp.ToJustFloat32_1D,
    x0: onp.ToFloat1D | None = None,
    *,
    rtol: float = 1e-5,
    shift: float = 0.0,
    maxiter: int | None = None,
    M: _ToLinearOperator[np.float32] | None = None,
    callback: Callable[[onp.Array1D[np.float32]], Unused] | None = None,
    show: bool = False,
    check: bool = False,
) -> tuple[onp.Array1D[np.float32], int]: ...
@overload  # c128
def minres(
    A: _ToLinearOperator[np.complex128],
    b: onp.ToComplex128_1D,
    x0: onp.ToComplex1D | None = None,
    *,
    rtol: float = 1e-5,
    shift: float = 0.0,
    maxiter: int | None = None,
    M: _ToLinearOperator[npc.number] | None = None,
    callback: Callable[[onp.Array1D[np.complex128]], Unused] | None = None,
    show: bool = False,
    check: bool = False,
) -> tuple[onp.Array1D[np.complex128], int]: ...
@overload  # c64
def minres(
    A: _ToLinearOperator[np.complex64],
    b: onp.ToJustFloat32_1D | onp.ToJustComplex64_1D,
    x0: onp.ToJustFloat32_1D | onp.ToJustComplex64_1D | None = None,
    *,
    rtol: float = 1e-5,
    shift: float = 0.0,
    maxiter: int | None = None,
    M: _ToLinearOperator[npc.inexact32] | None = None,
    callback: Callable[[onp.Array1D[np.complex64]], Unused] | None = None,
    show: bool = False,
    check: bool = False,
) -> tuple[onp.Array1D[np.complex64], int]: ...
@overload  # fallback
def minres(
    A: _ToLinearOperator[npc.number],
    b: onp.ToFloat1D,
    x0: onp.ToFloat1D | None = None,
    *,
    rtol: float = 1e-5,
    shift: float = 0.0,
    maxiter: int | None = None,
    M: _ToLinearOperator[npc.number] | None = None,
    callback: Callable[[onp.Array1D[Any]], Unused] | None = None,
    show: bool = False,
    check: bool = False,
) -> tuple[onp.Array1D[Any], int]: ...
