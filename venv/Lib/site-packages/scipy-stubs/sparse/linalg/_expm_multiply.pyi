from _typeshed import Incomplete
from typing import Any, Never, SupportsIndex, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._interface import LinearOperator
from scipy.sparse._base import _spbase, sparray

__all__ = ["expm_multiply"]

###

type _ToLinearOperator[ScalarT: npc.number | np.bool] = (
    LinearOperator[ScalarT] | _spbase[ScalarT, tuple[int, int]] | onp.ArrayND[ScalarT]
)
type _SparseOrDense[ScalarT: npc.number | np.bool, ShapeT: tuple[Any, ...]] = (
    sparray[ScalarT, ShapeT] | onp.ArrayND[ScalarT, ShapeT]
)

type _AsFloat64 = np.float64 | npc.integer | np.bool
type _ToFloat64 = _AsFloat64 | np.float32 | np.float16

# workaround for mypy's and pyright's typing spec non-compliance regarding overloads
type _JustAnyShape = tuple[Never, Never, Never]

###

@overload
def expm_multiply(
    A: _ToLinearOperator[_AsFloat64],
    B: _SparseOrDense[_ToFloat64, _JustAnyShape],
    start: onp.ToFloat | None = None,
    stop: onp.ToFloat | None = None,
    num: SupportsIndex | None = None,
    endpoint: bool | None = None,
    traceA: onp.ToFloat | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload
def expm_multiply[InexactT: npc.inexact](
    A: _ToLinearOperator[InexactT],
    B: _SparseOrDense[npc.integer | np.bool, _JustAnyShape],
    start: onp.ToFloat | None = None,
    stop: onp.ToFloat | None = None,
    num: SupportsIndex | None = None,
    endpoint: bool | None = None,
    traceA: onp.ToComplex | None = None,
) -> onp.ArrayND[InexactT]: ...
@overload
def expm_multiply[InexactT: npc.inexact](
    A: _ToLinearOperator[InexactT],
    B: _SparseOrDense[InexactT, _JustAnyShape],
    start: onp.ToFloat | None = None,
    stop: onp.ToFloat | None = None,
    num: SupportsIndex | None = None,
    endpoint: bool | None = None,
    traceA: onp.ToComplex | None = None,
) -> onp.ArrayND[InexactT]: ...
@overload  # 1-d
def expm_multiply(
    A: _ToLinearOperator[_AsFloat64],
    B: _SparseOrDense[_ToFloat64, tuple[int]],
    start: onp.ToFloat | None = None,
    stop: onp.ToFloat | None = None,
    num: SupportsIndex | None = None,
    endpoint: bool | None = None,
    traceA: onp.ToFloat | None = None,
) -> onp.Array1D[np.float64]: ...
@overload
def expm_multiply[InexactT: npc.inexact](
    A: _ToLinearOperator[InexactT],
    B: _SparseOrDense[npc.integer | np.bool, tuple[int]],
    start: onp.ToFloat | None = None,
    stop: onp.ToFloat | None = None,
    num: SupportsIndex | None = None,
    endpoint: bool | None = None,
    traceA: onp.ToComplex | None = None,
) -> onp.Array1D[InexactT]: ...
@overload
def expm_multiply[InexactT: npc.inexact](
    A: _ToLinearOperator[InexactT],
    B: _SparseOrDense[InexactT, tuple[int]],
    start: onp.ToFloat | None = None,
    stop: onp.ToFloat | None = None,
    num: SupportsIndex | None = None,
    endpoint: bool | None = None,
    traceA: onp.ToComplex | None = None,
) -> onp.Array1D[InexactT]: ...
@overload  # 2-d
def expm_multiply(
    A: _ToLinearOperator[_AsFloat64],
    B: _SparseOrDense[_ToFloat64, tuple[int, int]],
    start: onp.ToFloat | None = None,
    stop: onp.ToFloat | None = None,
    num: SupportsIndex | None = None,
    endpoint: bool | None = None,
    traceA: onp.ToFloat | None = None,
) -> onp.Array2D[np.float64]: ...
@overload
def expm_multiply[InexactT: npc.inexact](
    A: _ToLinearOperator[InexactT],
    B: _SparseOrDense[npc.integer | np.bool, tuple[int, int]],
    start: onp.ToFloat | None = None,
    stop: onp.ToFloat | None = None,
    num: SupportsIndex | None = None,
    endpoint: bool | None = None,
    traceA: onp.ToComplex | None = None,
) -> onp.Array2D[InexactT]: ...
@overload
def expm_multiply[InexactT: npc.inexact](
    A: _ToLinearOperator[InexactT],
    B: _SparseOrDense[InexactT, tuple[int, int]],
    start: onp.ToFloat | None = None,
    stop: onp.ToFloat | None = None,
    num: SupportsIndex | None = None,
    endpoint: bool | None = None,
    traceA: onp.ToComplex | None = None,
) -> onp.Array2D[InexactT]: ...
@overload  # 1-d or 2-d
def expm_multiply(
    A: _ToLinearOperator[_AsFloat64],
    B: _SparseOrDense[_ToFloat64, tuple[Any, ...]],
    start: onp.ToFloat | None = None,
    stop: onp.ToFloat | None = None,
    num: SupportsIndex | None = None,
    endpoint: bool | None = None,
    traceA: onp.ToFloat | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload
def expm_multiply[InexactT: npc.inexact](
    A: _ToLinearOperator[InexactT],
    B: _SparseOrDense[npc.integer | np.bool, tuple[Any, ...]],
    start: onp.ToFloat | None = None,
    stop: onp.ToFloat | None = None,
    num: SupportsIndex | None = None,
    endpoint: bool | None = None,
    traceA: onp.ToComplex | None = None,
) -> onp.ArrayND[InexactT]: ...
@overload
def expm_multiply[InexactT: npc.inexact](
    A: _ToLinearOperator[InexactT],
    B: _SparseOrDense[InexactT, tuple[Any, ...]],
    start: onp.ToFloat | None = None,
    stop: onp.ToFloat | None = None,
    num: SupportsIndex | None = None,
    endpoint: bool | None = None,
    traceA: onp.ToComplex | None = None,
) -> onp.ArrayND[InexactT]: ...
@overload  # fallback
def expm_multiply(
    A: _ToLinearOperator[npc.number],
    B: _SparseOrDense[npc.number, tuple[Any, ...]],
    start: onp.ToFloat | None = None,
    stop: onp.ToFloat | None = None,
    num: SupportsIndex | None = None,
    endpoint: bool | None = None,
    traceA: onp.ToComplex | None = None,
) -> onp.ArrayND[Incomplete]: ...
