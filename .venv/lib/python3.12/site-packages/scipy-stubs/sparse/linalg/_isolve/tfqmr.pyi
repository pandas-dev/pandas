from collections.abc import Callable
from typing import TypeAlias, overload
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse._base import _spbase
from scipy.sparse.linalg import LinearOperator

__all__ = ["tfqmr"]

_FloatT = TypeVar("_FloatT", bound=np.float32 | np.float64, default=np.float64)
_ComplexT = TypeVar("_ComplexT", bound=np.complex64 | np.complex128, default=np.complex128)
_ScalarT = TypeVar("_ScalarT", bound=npc.number | np.bool_)

_ToLinearOperator: TypeAlias = onp.CanArrayND[_ScalarT] | _spbase[_ScalarT] | LinearOperator[_ScalarT]

_Ignored: TypeAlias = object
_Callback: TypeAlias = Callable[[onp.Array1D[_ScalarT]], _Ignored]

###

@overload  # real
def tfqmr(
    A: _ToLinearOperator[_FloatT | npc.integer | np.bool_],
    b: onp.ToFloat1D,
    x0: onp.ToFloat1D | None = None,
    *,
    rtol: onp.ToFloat = 1e-5,
    atol: onp.ToFloat = 0.0,
    maxiter: int | None = None,
    M: _ToLinearOperator[_FloatT] | None = None,
    callback: _Callback[_FloatT] | None = None,
    show: onp.ToBool = False,
) -> tuple[onp.Array1D[_FloatT], int]: ...
@overload  # complex
def tfqmr(
    A: _ToLinearOperator[_ComplexT],
    b: onp.ToComplex1D,
    x0: onp.ToComplex1D | None = None,
    *,
    rtol: onp.ToFloat = 1e-5,
    atol: onp.ToFloat = 0.0,
    maxiter: int | None = None,
    M: _ToLinearOperator[_ComplexT] | None = None,
    callback: _Callback[_ComplexT] | None = None,
    show: onp.ToBool = False,
) -> tuple[onp.Array1D[_ComplexT], int]: ...
