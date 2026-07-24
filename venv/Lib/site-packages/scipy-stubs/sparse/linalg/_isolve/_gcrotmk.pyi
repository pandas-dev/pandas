from _typeshed import Unused
from collections.abc import Callable, Sequence
from typing import Literal, overload
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse._base import _spbase
from scipy.sparse.linalg import LinearOperator

__all__ = ["gcrotmk"]

###

type _Float = np.float32 | np.float64
type _Complex = np.complex64 | np.complex128
type _Inexact = _Float | _Complex

type _ToInt = npc.integer | np.bool
type _ToLinearOperator[_NumericT: npc.number | np.bool] = (
    onp.CanArrayND[_NumericT] | _spbase[_NumericT] | LinearOperator[_NumericT]
)

type _Callback[_NumericT: npc.number | np.bool] = Callable[[onp.Array1D[_NumericT]], Unused]

type _Truncate = Literal["oldest", "newest"]

_FloatT = TypeVar("_FloatT", bound=_Float, default=np.float64)
_ComplexT = TypeVar("_ComplexT", bound=_Complex, default=np.complex128)

###

@overload
def gcrotmk(
    A: _ToLinearOperator[_FloatT | _ToInt],
    b: onp.ToFloat1D,
    x0: onp.ToFloat1D | None = None,
    *,
    rtol: onp.ToFloat = 1e-5,
    atol: onp.ToFloat = 0.0,
    maxiter: int = 1_000,
    M: _ToLinearOperator[_FloatT | _ToInt] | None = None,
    callback: _Callback[_FloatT] | None = None,
    m: int = 20,
    k: int | None = None,
    CU: Sequence[tuple[Sequence[onp.ArrayND[_Float]], Sequence[onp.ArrayND[_Float]] | None]] | None = None,
    discard_C: bool = False,
    truncate: _Truncate = "oldest",
) -> tuple[onp.Array1D[_FloatT], int]: ...
@overload
def gcrotmk(
    A: _ToLinearOperator[_ComplexT],
    b: onp.ToComplex1D,
    x0: onp.ToComplex1D | None = None,
    *,
    rtol: onp.ToFloat = 1e-5,
    atol: onp.ToFloat = 0.0,
    maxiter: int = 1_000,
    M: _ToLinearOperator[_ComplexT] | None = None,
    callback: _Callback[_ComplexT] | None = None,
    m: int = 20,
    k: int | None = None,
    CU: Sequence[tuple[Sequence[onp.ArrayND[_Inexact]], Sequence[onp.ArrayND[_Inexact]] | None]] | None = None,
    discard_C: bool = False,
    truncate: _Truncate = "oldest",
) -> tuple[onp.Array1D[_ComplexT], int]: ...
