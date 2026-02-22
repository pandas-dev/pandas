from collections.abc import Callable, Sequence
from typing import Literal, TypeAlias, overload
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse._base import _spbase
from scipy.sparse.linalg import LinearOperator

__all__ = ["gcrotmk"]

_Float: TypeAlias = np.float32 | np.float64
_Complex: TypeAlias = np.complex64 | np.complex128
_Inexact: TypeAlias = _Float | _Complex

_FloatT = TypeVar("_FloatT", bound=_Float, default=np.float64)
_ComplexT = TypeVar("_ComplexT", bound=_Complex, default=np.complex128)
_NumericT = TypeVar("_NumericT", bound=npc.number | np.bool_)

_ToInt: TypeAlias = npc.integer | np.bool_
_ToLinearOperator: TypeAlias = onp.CanArrayND[_NumericT] | _spbase[_NumericT] | LinearOperator[_NumericT]

_Ignored: TypeAlias = object
_Callback: TypeAlias = Callable[[onp.Array1D[_NumericT]], _Ignored]

_Truncate: TypeAlias = Literal["oldest", "newest"]

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
    discard_C: onp.ToBool = False,
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
    discard_C: onp.ToBool = False,
    truncate: _Truncate = "oldest",
) -> tuple[onp.Array1D[_ComplexT], int]: ...
