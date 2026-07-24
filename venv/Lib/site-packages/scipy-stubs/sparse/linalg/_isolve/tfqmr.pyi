from _typeshed import Unused
from collections.abc import Callable
from typing import overload
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse._base import _spbase
from scipy.sparse.linalg import LinearOperator

__all__ = ["tfqmr"]

###

type _ToLinearOperator[_ScalarT: npc.number | np.bool] = onp.CanArrayND[_ScalarT] | _spbase[_ScalarT] | LinearOperator[_ScalarT]
type _Callback[_ScalarT: npc.number | np.bool] = Callable[[onp.Array1D[_ScalarT]], Unused]

_FloatT = TypeVar("_FloatT", bound=np.float32 | np.float64, default=np.float64)
_ComplexT = TypeVar("_ComplexT", bound=np.complex64 | np.complex128, default=np.complex128)

###

@overload  # real
def tfqmr(
    A: _ToLinearOperator[_FloatT | npc.integer | np.bool],
    b: onp.ToFloat1D,
    x0: onp.ToFloat1D | None = None,
    *,
    rtol: float = 1e-5,
    atol: float = 0.0,
    maxiter: int | None = None,
    M: _ToLinearOperator[_FloatT] | None = None,
    callback: _Callback[_FloatT] | None = None,
    show: bool = False,
) -> tuple[onp.Array1D[_FloatT], int]: ...
@overload  # complex
def tfqmr(
    A: _ToLinearOperator[_ComplexT],
    b: onp.ToComplex1D,
    x0: onp.ToComplex1D | None = None,
    *,
    rtol: float = 1e-5,
    atol: float = 0.0,
    maxiter: int | None = None,
    M: _ToLinearOperator[_ComplexT] | None = None,
    callback: _Callback[_ComplexT] | None = None,
    show: bool = False,
) -> tuple[onp.Array1D[_ComplexT], int]: ...
