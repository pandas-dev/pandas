from _typeshed import Unused
from collections.abc import Callable
from typing import overload
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse._base import _spbase
from scipy.sparse.linalg import LinearOperator

__all__ = ["lgmres"]

###

type _Float = np.float32 | np.float64
type _Complex = np.complex64 | np.complex128

type _ToInt = npc.integer | np.bool
type _ToLinearOperator[_ScalarT: npc.number | np.bool] = onp.CanArrayND[_ScalarT] | _spbase[_ScalarT] | LinearOperator[_ScalarT]

type _Callback[_ScalarT: npc.number | np.bool] = Callable[[onp.Array1D[_ScalarT]], Unused]

_FloatT = TypeVar("_FloatT", bound=_Float, default=np.float64)
_ComplexT = TypeVar("_ComplexT", bound=_Complex, default=np.complex128)

###

@overload
def lgmres(
    A: _ToLinearOperator[_FloatT | _ToInt],
    b: onp.ToFloat1D,
    x0: onp.ToFloat1D | None = None,
    *,
    rtol: onp.ToFloat = 1e-5,
    atol: onp.ToFloat = 0.0,
    maxiter: int = 1_000,
    M: _ToLinearOperator[_FloatT | _ToInt] | None = None,
    callback: _Callback[_FloatT] | None = None,
    inner_m: int = 30,
    outer_k: int = 3,
    outer_v: list[tuple[onp.ArrayND[_Float], onp.ArrayND[_Float] | None]] | None = None,
    store_outer_Av: bool = True,
    prepend_outer_v: bool = False,
) -> tuple[onp.Array1D[_FloatT], int]: ...
@overload
def lgmres(
    A: _ToLinearOperator[_ComplexT],
    b: onp.ToComplex1D,
    x0: onp.ToComplex1D | None = None,
    *,
    rtol: onp.ToFloat = 1e-5,
    atol: onp.ToFloat = 0.0,
    maxiter: int = 1_000,
    M: _ToLinearOperator[_ComplexT] | None = None,
    callback: _Callback[_ComplexT] | None = None,
    inner_m: int = 30,
    outer_k: int = 3,
    outer_v: list[tuple[onp.ArrayND[_Float | _Complex], onp.ArrayND[_Complex] | None]] | None = None,
    store_outer_Av: bool = True,
    prepend_outer_v: bool = False,
) -> tuple[onp.Array1D[_ComplexT], int]: ...
