from collections.abc import Callable, Iterable, Mapping
from typing import Concatenate, Literal, TypeAlias, TypedDict, overload, type_check_only

import numpy as np
import optype.numpy as onp

from ._differentiable_functions import _DoesMap
from scipy.sparse import csr_array, sparray, spmatrix
from scipy.sparse.linalg import LinearOperator

_ToSparse: TypeAlias = spmatrix | sparray

_Float1D: TypeAlias = onp.Array1D[np.float64]

_ToFloatOr1D: TypeAlias = onp.ToFloat | onp.ToFloat1D
_ToFloatMatrix: TypeAlias = onp.ToFloat2D | _ToSparse

_ToIntMatrix: TypeAlias = onp.ToInt2D | _ToSparse
_Sparsity: TypeAlias = _ToIntMatrix | tuple[_ToIntMatrix, onp.ToFloat1D]

_Fun: TypeAlias = Callable[Concatenate[_Float1D, ...], _ToFloatOr1D]
_Jac: TypeAlias = Callable[Concatenate[_Float1D, ...], onp.ToFloat1D | _ToFloatMatrix]

@type_check_only
class _InfoDict(TypedDict):
    nfev: int

###

def group_columns(A: _ToFloatMatrix, order: onp.ToInt | Iterable[onp.ToInt] = 0) -> onp.Array1D[np.intp]: ...

#
@overload  # sparsity: None (default), as_linear_operator: False (default), full_output: False (default)
def approx_derivative(
    fun: _Fun,
    x0: _ToFloatOr1D,
    method: Literal["2-point", "3-point", "cs"] = "3-point",
    rel_step: _ToFloatOr1D | None = None,
    abs_step: _ToFloatOr1D | None = None,
    f0: _ToFloatOr1D | None = None,
    bounds: tuple[_ToFloatOr1D, _ToFloatOr1D] = ...,
    sparsity: None = None,
    as_linear_operator: onp.ToFalse = False,
    args: tuple[object, ...] = (),
    kwargs: Mapping[str, object] | None = None,
    full_output: onp.ToFalse = False,
    workers: int | _DoesMap | None = None,
) -> _Float1D | onp.Array2D[np.float64]: ...
@overload  # sparsity: <given>, as_linear_operator: False (default), full_output: False (default)
def approx_derivative(
    fun: _Fun,
    x0: _ToFloatOr1D,
    method: Literal["2-point", "3-point", "cs"] = "3-point",
    rel_step: _ToFloatOr1D | None = None,
    abs_step: _ToFloatOr1D | None = None,
    f0: _ToFloatOr1D | None = None,
    bounds: tuple[_ToFloatOr1D, _ToFloatOr1D] = ...,
    *,
    sparsity: _Sparsity,
    as_linear_operator: onp.ToFalse = False,
    args: tuple[object, ...] = (),
    kwargs: Mapping[str, object] | None = None,
    full_output: onp.ToFalse = False,
    workers: int | _DoesMap | None = None,
) -> csr_array: ...
@overload  # as_linear_operator: True, full_output: False (default)
def approx_derivative(
    fun: _Fun,
    x0: _ToFloatOr1D,
    method: Literal["2-point", "3-point", "cs"] = "3-point",
    rel_step: _ToFloatOr1D | None = None,
    abs_step: _ToFloatOr1D | None = None,
    f0: _ToFloatOr1D | None = None,
    bounds: tuple[_ToFloatOr1D, _ToFloatOr1D] = ...,
    sparsity: _Sparsity | None = None,
    *,
    as_linear_operator: onp.ToTrue,
    args: tuple[object, ...] = (),
    kwargs: Mapping[str, object] | None = None,
    full_output: onp.ToFalse = False,
    workers: int | _DoesMap | None = None,
) -> LinearOperator: ...
@overload  # sparsity: None (default), as_linear_operator: False (default), full_output: True (keyword)
def approx_derivative(
    fun: _Fun,
    x0: _ToFloatOr1D,
    method: Literal["2-point", "3-point", "cs"] = "3-point",
    rel_step: _ToFloatOr1D | None = None,
    abs_step: _ToFloatOr1D | None = None,
    f0: _ToFloatOr1D | None = None,
    bounds: tuple[_ToFloatOr1D, _ToFloatOr1D] = ...,
    sparsity: None = None,
    as_linear_operator: onp.ToFalse = False,
    args: tuple[object, ...] = (),
    kwargs: Mapping[str, object] | None = None,
    *,
    full_output: onp.ToTrue,
    workers: int | _DoesMap | None = None,
) -> tuple[_Float1D | onp.Array2D[np.float64], _InfoDict]: ...
@overload  # sparsity: <given>, as_linear_operator: False (default), full_output: True
def approx_derivative(
    fun: _Fun,
    x0: _ToFloatOr1D,
    method: Literal["2-point", "3-point", "cs"] = "3-point",
    rel_step: _ToFloatOr1D | None = None,
    abs_step: _ToFloatOr1D | None = None,
    f0: _ToFloatOr1D | None = None,
    bounds: tuple[_ToFloatOr1D, _ToFloatOr1D] = ...,
    *,
    sparsity: _Sparsity,
    as_linear_operator: onp.ToFalse = False,
    args: tuple[object, ...] = (),
    kwargs: Mapping[str, object] | None = None,
    full_output: onp.ToTrue,
    workers: int | _DoesMap | None = None,
) -> tuple[csr_array, _InfoDict]: ...
@overload  # as_linear_operator: True, full_output: True
def approx_derivative(
    fun: _Fun,
    x0: _ToFloatOr1D,
    method: Literal["2-point", "3-point", "cs"] = "3-point",
    rel_step: _ToFloatOr1D | None = None,
    abs_step: _ToFloatOr1D | None = None,
    f0: _ToFloatOr1D | None = None,
    bounds: tuple[_ToFloatOr1D, _ToFloatOr1D] = ...,
    sparsity: _Sparsity | None = None,
    *,
    as_linear_operator: onp.ToTrue,
    args: tuple[object, ...] = (),
    kwargs: Mapping[str, object] | None = None,
    full_output: onp.ToTrue,
    workers: int | _DoesMap | None = None,
) -> tuple[LinearOperator, _InfoDict]: ...

#
def check_derivative(
    fun: _Fun,
    jac: _Jac,
    x0: _ToFloatOr1D,
    bounds: tuple[_ToFloatOr1D, _ToFloatOr1D] = ...,
    args: tuple[object, ...] = (),
    kwargs: Mapping[str, object] | None = None,
) -> float | np.float64: ...
