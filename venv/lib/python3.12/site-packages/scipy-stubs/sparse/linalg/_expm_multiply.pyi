from _typeshed import Incomplete
from typing import Any, Never, TypeAlias, TypeVar, overload
from typing_extensions import TypeAliasType

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from ._interface import LinearOperator
from scipy.sparse._base import _spbase, sparray

__all__ = ["expm_multiply"]

_ScalarT = TypeVar("_ScalarT", bound=npc.number | np.bool_)
_InexactT = TypeVar("_InexactT", bound=npc.inexact)
_ShapeT = TypeVar("_ShapeT", bound=tuple[Any, ...])

_ToLinearOperator: TypeAlias = LinearOperator[_ScalarT] | _spbase[_ScalarT, tuple[int, int]] | onp.ArrayND[_ScalarT]
_SparseOrDense = TypeAliasType(
    "_SparseOrDense", sparray[_ScalarT, _ShapeT] | onp.ArrayND[_ScalarT, _ShapeT], type_params=(_ScalarT, _ShapeT)
)

_AsFloat64: TypeAlias = np.float64 | npc.integer | np.bool_
_ToFloat64: TypeAlias = _AsFloat64 | np.float32 | np.float16

###

@overload  # workaround for mypy's and pyright's typing spec non-compliance regarding overloads
def expm_multiply(
    A: _ToLinearOperator[_AsFloat64],
    B: _SparseOrDense[_ToFloat64, tuple[Never] | tuple[Never, Never]],
    start: onp.ToFloat | None = None,
    stop: onp.ToFloat | None = None,
    num: op.CanIndex | None = None,
    endpoint: bool | None = None,
    traceA: onp.ToComplex | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload
def expm_multiply(
    A: _ToLinearOperator[_InexactT],
    B: _SparseOrDense[_InexactT | npc.integer | np.bool_, tuple[Never] | tuple[Never, Never]],
    start: onp.ToFloat | None = None,
    stop: onp.ToFloat | None = None,
    num: op.CanIndex | None = None,
    endpoint: bool | None = None,
    traceA: onp.ToComplex | None = None,
) -> onp.ArrayND[_InexactT]: ...
@overload  # 1-d
def expm_multiply(
    A: _ToLinearOperator[_AsFloat64],
    B: _SparseOrDense[_ToFloat64, tuple[int]],
    start: onp.ToFloat | None = None,
    stop: onp.ToFloat | None = None,
    num: op.CanIndex | None = None,
    endpoint: bool | None = None,
    traceA: onp.ToComplex | None = None,
) -> onp.Array1D[np.float64]: ...
@overload
def expm_multiply(
    A: _ToLinearOperator[_InexactT],
    B: _SparseOrDense[_InexactT | npc.integer | np.bool_, tuple[int]],
    start: onp.ToFloat | None = None,
    stop: onp.ToFloat | None = None,
    num: op.CanIndex | None = None,
    endpoint: bool | None = None,
    traceA: onp.ToComplex | None = None,
) -> onp.Array1D[_InexactT]: ...
@overload  # 2-d
def expm_multiply(
    A: _ToLinearOperator[_AsFloat64],
    B: _SparseOrDense[_ToFloat64, tuple[int, int]],
    start: onp.ToFloat | None = None,
    stop: onp.ToFloat | None = None,
    num: op.CanIndex | None = None,
    endpoint: bool | None = None,
    traceA: onp.ToComplex | None = None,
) -> onp.Array2D[np.float64]: ...
@overload
def expm_multiply(
    A: _ToLinearOperator[_InexactT],
    B: _SparseOrDense[_InexactT | npc.integer | np.bool_, tuple[int, int]],
    start: onp.ToFloat | None = None,
    stop: onp.ToFloat | None = None,
    num: op.CanIndex | None = None,
    endpoint: bool | None = None,
    traceA: onp.ToComplex | None = None,
) -> onp.Array2D[_InexactT]: ...
@overload  # 1-d or 2-d
def expm_multiply(
    A: _ToLinearOperator[_AsFloat64],
    B: _SparseOrDense[_ToFloat64, tuple[Any, ...]],
    start: onp.ToFloat | None = None,
    stop: onp.ToFloat | None = None,
    num: op.CanIndex | None = None,
    endpoint: bool | None = None,
    traceA: onp.ToComplex | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload
def expm_multiply(
    A: _ToLinearOperator[_InexactT],
    B: _SparseOrDense[_InexactT | npc.integer | np.bool_, tuple[Any, ...]],
    start: onp.ToFloat | None = None,
    stop: onp.ToFloat | None = None,
    num: op.CanIndex | None = None,
    endpoint: bool | None = None,
    traceA: onp.ToComplex | None = None,
) -> onp.ArrayND[_InexactT]: ...
@overload  # fallback
def expm_multiply(
    A: _ToLinearOperator[npc.number],
    B: _SparseOrDense[npc.number, tuple[Any, ...]],
    start: onp.ToFloat | None = None,
    stop: onp.ToFloat | None = None,
    num: op.CanIndex | None = None,
    endpoint: bool | None = None,
    traceA: onp.ToComplex | None = None,
) -> onp.ArrayND[Incomplete]: ...
