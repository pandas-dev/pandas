from collections.abc import Iterable, Sequence
from typing import Literal, TypeAlias, overload
from typing_extensions import TypeVar

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc
import optype.typing as opt
from numpy._typing import _ArrayLike

__all__ = [
    "cdf2rdf",
    "eig",
    "eig_banded",
    "eigh",
    "eigh_tridiagonal",
    "eigvals",
    "eigvals_banded",
    "eigvalsh",
    "eigvalsh_tridiagonal",
    "hessenberg",
]

_FloatT = TypeVar("_FloatT", bound=_Floating, default=_Float)
_FloatT2 = TypeVar("_FloatT2", bound=_Floating, default=_Float)

# scalar types
_Integer: TypeAlias = npc.integer
_Floating: TypeAlias = npc.floating
_Float: TypeAlias = np.float32 | np.float64
_Complex: TypeAlias = np.complex64 | np.complex128

# input types
# NOTE: only "a", "v" and "i" are documented for the `select` params, but internally 0, 1, and 2 are used, respectively.
_SelectA: TypeAlias = Literal["a", "all", 0]
_SelectV: TypeAlias = Literal["v", "value", 1]
_SelectI: TypeAlias = Literal["i", "index", 2]

# NOTE: `_check_select()` requires the `select_range` array-like to be of `int{16,32,64}` when `select: _SelectIndex`
# https://github.com/scipy/scipy-stubs/issues/154
# NOTE: This `select_range` parameter type must be of shape `(2,)` and in nondescending order
_SelectRange: TypeAlias = Sequence[float | _Integer | _Floating]
_SelectRangeI: TypeAlias = Sequence[int | np.int16 | np.int32 | np.int64]  # no bool, int8 or unsigned ints

_EigHType: TypeAlias = Literal[1, 2, 3]
_EigHSubsetByIndex: TypeAlias = Iterable[opt.AnyInt]
_EigHSubsetByValue: TypeAlias = Iterable[onp.ToFloat]

# LAPACK drivers
_DriverGV: TypeAlias = Literal["gv", "gvd", "gvx"]
_DriverEV: TypeAlias = Literal["ev", "evd", "evx", "evr"]
_DriverSTE: TypeAlias = Literal["stemr", "stebz", "sterf", "stev"]
_DriverAuto: TypeAlias = Literal["auto"]

# output types
_FloatND: TypeAlias = onp.ArrayND[_Float]
_ComplexND: TypeAlias = onp.ArrayND[_Complex]
_InexactND: TypeAlias = onp.ArrayND[_Float | _Complex]

###

# `eig` has `2 * (6 + 7) + 1 == 27` overloads...
@overload  # float, left: True (positional), right: True = ..., homogeneous_eigvals: False = ...
def eig(
    a: onp.ToFloatND,
    b: onp.ToFloatND | None,
    left: onp.ToTrue,
    right: onp.ToTrue = True,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    check_finite: op.CanBool = True,
    homogeneous_eigvals: onp.ToFalse = False,
) -> tuple[_ComplexND, _FloatND]: ...
@overload  # float, left: True (positional), right: True = ..., homogeneous_eigvals: True (keyword)
def eig(
    a: onp.ToFloatND,
    b: onp.ToFloatND | None,
    left: onp.ToTrue,
    right: onp.ToTrue = True,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    check_finite: op.CanBool = True,
    *,
    homogeneous_eigvals: onp.ToTrue,
) -> tuple[_ComplexND, _FloatND]: ...
@overload  # float, left: True (keyword), right: True = ..., homogeneous_eigvals: False = ...
def eig(
    a: onp.ToFloatND,
    b: onp.ToFloatND | None = None,
    *,
    left: onp.ToTrue,
    right: onp.ToTrue = True,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    check_finite: op.CanBool = True,
    homogeneous_eigvals: onp.ToFalse = False,
) -> tuple[_ComplexND, _FloatND]: ...
@overload  # float, left: True (keyword), right: True = ..., homogeneous_eigvals: True
def eig(
    a: onp.ToFloatND,
    b: onp.ToFloatND | None = None,
    *,
    left: onp.ToTrue,
    right: onp.ToTrue = True,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    check_finite: op.CanBool = True,
    homogeneous_eigvals: onp.ToTrue,
) -> tuple[_ComplexND, _FloatND]: ...
@overload  # float, left: False, right: False (positional), homogeneous_eigvals: False = ...
def eig(
    a: onp.ToFloatND,
    b: onp.ToFloatND | None,
    left: onp.ToFalse,
    right: onp.ToFalse,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    check_finite: op.CanBool = True,
    homogeneous_eigvals: onp.ToFalse = False,
) -> tuple[_ComplexND, _FloatND]: ...
@overload  # float, left: False, right: False (positional), homogeneous_eigvals: True (keyword)
def eig(
    a: onp.ToFloatND,
    b: onp.ToFloatND | None,
    left: onp.ToFalse,
    right: onp.ToFalse,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    check_finite: op.CanBool = True,
    *,
    homogeneous_eigvals: onp.ToTrue,
) -> tuple[_ComplexND, _FloatND]: ...
@overload  # float, left: False = ..., right: False (keyword), homogeneous_eigvals: False = ...
def eig(
    a: onp.ToFloatND,
    b: onp.ToFloatND | None = None,
    left: onp.ToFalse = False,
    *,
    right: onp.ToFalse,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    check_finite: op.CanBool = True,
    homogeneous_eigvals: onp.ToFalse = False,
) -> tuple[_ComplexND, _FloatND]: ...
@overload  # float, left: False = ..., right: False (keyword), homogeneous_eigvals: True
def eig(
    a: onp.ToFloatND,
    b: onp.ToFloatND | None = None,
    left: onp.ToFalse = False,
    *,
    right: onp.ToFalse,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    check_finite: op.CanBool = True,
    homogeneous_eigvals: onp.ToTrue,
) -> tuple[_ComplexND, _FloatND]: ...
@overload  # float, left: True (positional), right: False, homogeneous_eigvals: False = ...
def eig(
    a: onp.ToFloatND,
    b: onp.ToFloatND | None,
    left: onp.ToTrue,
    right: onp.ToFalse,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    check_finite: op.CanBool = True,
    homogeneous_eigvals: onp.ToFalse = False,
) -> tuple[_ComplexND, _FloatND, _FloatND]: ...
@overload  # float, left: True (positional), right: False, homogeneous_eigvals: True (keyword)
def eig(
    a: onp.ToFloatND,
    b: onp.ToFloatND | None,
    left: onp.ToTrue,
    right: onp.ToFalse,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    check_finite: op.CanBool = True,
    *,
    homogeneous_eigvals: onp.ToTrue,
) -> tuple[_ComplexND, _FloatND, _FloatND]: ...
@overload  # float, left: True (keyword), right: False, homogeneous_eigvals: False = ...
def eig(
    a: onp.ToFloatND,
    b: onp.ToFloatND | None = None,
    *,
    left: onp.ToTrue,
    right: onp.ToFalse,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    check_finite: op.CanBool = True,
    homogeneous_eigvals: onp.ToFalse = False,
) -> tuple[_ComplexND, _FloatND, _FloatND]: ...
@overload  # float, left: True (keyword), right: False, homogeneous_eigvals: True
def eig(
    a: onp.ToFloatND,
    b: onp.ToFloatND | None = None,
    *,
    left: onp.ToTrue,
    right: onp.ToFalse,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    check_finite: op.CanBool = True,
    homogeneous_eigvals: onp.ToTrue,
) -> tuple[_ComplexND, _FloatND, _FloatND]: ...
@overload  # complex, left: False = ..., right: True = ..., homogeneous_eigvals: False = ...
def eig(
    a: onp.ToComplexND,
    b: onp.ToComplexND | None = None,
    left: onp.ToFalse = False,
    right: onp.ToTrue = True,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    check_finite: op.CanBool = True,
    homogeneous_eigvals: onp.ToFalse = False,
) -> _ComplexND: ...
@overload  # complex, left: False = ..., right: True = ..., homogeneous_eigvals: True (keyword)
def eig(
    a: onp.ToComplexND,
    b: onp.ToComplexND | None = None,
    left: onp.ToFalse = False,
    right: onp.ToTrue = True,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    check_finite: op.CanBool = True,
    *,
    homogeneous_eigvals: onp.ToTrue,
) -> _ComplexND: ...
@overload  # complex, left: True (positional), right: True = ..., homogeneous_eigvals: False = ...
def eig(
    a: onp.ToComplexND,
    b: onp.ToComplexND | None,
    left: onp.ToTrue,
    right: onp.ToTrue = True,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    check_finite: op.CanBool = True,
    homogeneous_eigvals: onp.ToFalse = False,
) -> tuple[_ComplexND, _InexactND]: ...
@overload  # complex, left: True (positional), right: True = ..., homogeneous_eigvals: True (keyword)
def eig(
    a: onp.ToComplexND,
    b: onp.ToComplexND | None,
    left: onp.ToTrue,
    right: onp.ToTrue = True,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    check_finite: op.CanBool = True,
    *,
    homogeneous_eigvals: onp.ToTrue,
) -> tuple[_ComplexND, _InexactND]: ...
@overload  # complex, left: True (keyword), right: True = ..., homogeneous_eigvals: False = ...
def eig(
    a: onp.ToComplexND,
    b: onp.ToComplexND | None = None,
    *,
    left: onp.ToTrue,
    right: onp.ToTrue = True,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    check_finite: op.CanBool = True,
    homogeneous_eigvals: onp.ToFalse = False,
) -> tuple[_ComplexND, _InexactND]: ...
@overload  # complex, left: True (keyword), right: True = ..., homogeneous_eigvals: True
def eig(
    a: onp.ToComplexND,
    b: onp.ToComplexND | None = None,
    *,
    left: onp.ToTrue,
    right: onp.ToTrue = True,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    check_finite: op.CanBool = True,
    homogeneous_eigvals: onp.ToTrue,
) -> tuple[_ComplexND, _InexactND]: ...
@overload  # complex, left: False, right: False (positional), homogeneous_eigvals: False = ...
def eig(
    a: onp.ToComplexND,
    b: onp.ToComplexND | None,
    left: onp.ToFalse,
    right: onp.ToFalse,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    check_finite: op.CanBool = True,
    homogeneous_eigvals: onp.ToFalse = False,
) -> tuple[_ComplexND, _InexactND]: ...
@overload  # complex, left: False, right: False (positional), homogeneous_eigvals: True (keyword)
def eig(
    a: onp.ToComplexND,
    b: onp.ToComplexND | None,
    left: onp.ToFalse,
    right: onp.ToFalse,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    check_finite: op.CanBool = True,
    *,
    homogeneous_eigvals: onp.ToTrue,
) -> tuple[_ComplexND, _InexactND]: ...
@overload  # complex, left: False = ..., right: False (keyword), homogeneous_eigvals: False = ...
def eig(
    a: onp.ToComplexND,
    b: onp.ToComplexND | None = None,
    left: onp.ToFalse = False,
    *,
    right: onp.ToFalse,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    check_finite: op.CanBool = True,
    homogeneous_eigvals: onp.ToFalse = False,
) -> tuple[_ComplexND, _InexactND]: ...
@overload  # complex, left: False = ..., right: False (keyword), homogeneous_eigvals: True
def eig(
    a: onp.ToComplexND,
    b: onp.ToComplexND | None = None,
    left: onp.ToFalse = False,
    *,
    right: onp.ToFalse,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    check_finite: op.CanBool = True,
    homogeneous_eigvals: onp.ToTrue,
) -> tuple[_ComplexND, _InexactND]: ...
@overload  # complex, left: True (positional), right: False, homogeneous_eigvals: False = ...
def eig(
    a: onp.ToComplexND,
    b: onp.ToComplexND | None,
    left: onp.ToTrue,
    right: onp.ToFalse,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    check_finite: op.CanBool = True,
    homogeneous_eigvals: onp.ToFalse = False,
) -> tuple[_ComplexND, _InexactND, _InexactND]: ...
@overload  # complex, left: True (positional), right: False, homogeneous_eigvals: True (keyword)
def eig(
    a: onp.ToComplexND,
    b: onp.ToComplexND | None,
    left: onp.ToTrue,
    right: onp.ToFalse,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    check_finite: op.CanBool = True,
    *,
    homogeneous_eigvals: onp.ToTrue,
) -> tuple[_ComplexND, _InexactND, _InexactND]: ...
@overload  # complex, left: True (keyword), right: False (keyword), homogeneous_eigvals: False = ...
def eig(
    a: onp.ToComplexND,
    b: onp.ToComplexND | None = None,
    *,
    left: onp.ToTrue,
    right: onp.ToFalse,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    check_finite: op.CanBool = True,
    homogeneous_eigvals: onp.ToFalse = False,
) -> tuple[_ComplexND, _InexactND, _InexactND]: ...
@overload  # complex, left: True (keyword), right: False (keyword), homogeneous_eigvals: True
def eig(
    a: onp.ToComplexND,
    b: onp.ToComplexND | None = None,
    *,
    left: onp.ToTrue,
    right: onp.ToFalse,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    check_finite: op.CanBool = True,
    homogeneous_eigvals: onp.ToTrue,
) -> tuple[_ComplexND, _InexactND, _InexactND]: ...
@overload  # catch-all
def eig(
    a: onp.ToComplexND,
    b: onp.ToComplexND | None = None,
    left: onp.ToBool = False,
    right: onp.ToBool = True,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    check_finite: op.CanBool = True,
    homogeneous_eigvals: onp.ToBool = False,
) -> (
    _ComplexND
    | tuple[_ComplexND, _InexactND]
    | tuple[_ComplexND, _InexactND, _InexactND]
):  # fmt: skip
    ...

#
@overload  # +float64, eigvals_only: False = ...
def eigh(  #
    a: onp.ToArrayND[float, np.float64 | npc.floating80 | npc.integer64 | npc.integer32],
    b: onp.ToFloat64_ND | None = None,
    *,
    lower: op.CanBool = True,
    eigvals_only: onp.ToFalse = False,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    type: _EigHType = 1,
    check_finite: op.CanBool = True,
    subset_by_index: _EigHSubsetByIndex | None = None,
    subset_by_value: _EigHSubsetByValue | None = None,
    driver: _DriverEV | _DriverGV | None = None,
) -> tuple[onp.ArrayND[np.float64], onp.ArrayND[np.float64]]: ...
@overload  # +float, eigvals_only: False = ...
def eigh(
    a: onp.ToFloatND,
    b: onp.ToFloatND | None = None,
    *,
    lower: op.CanBool = True,
    eigvals_only: onp.ToFalse = False,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    type: _EigHType = 1,
    check_finite: op.CanBool = True,
    subset_by_index: _EigHSubsetByIndex | None = None,
    subset_by_value: _EigHSubsetByValue | None = None,
    driver: _DriverEV | _DriverGV | None = None,
) -> tuple[_FloatND, _FloatND]: ...
@overload  # ~complex, eigvals_only: False = ...
def eigh(
    a: onp.ToJustComplexND,
    b: onp.ToComplexND | None = None,
    *,
    lower: op.CanBool = True,
    eigvals_only: onp.ToFalse = False,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    type: _EigHType = 1,
    check_finite: op.CanBool = True,
    subset_by_index: _EigHSubsetByIndex | None = None,
    subset_by_value: _EigHSubsetByValue | None = None,
    driver: _DriverEV | _DriverGV | None = None,
) -> tuple[_FloatND, _ComplexND]: ...
@overload  # +complex, eigvals_only: False = ...
def eigh(
    a: onp.ToComplexND,
    b: onp.ToComplexND | None = None,
    *,
    lower: op.CanBool = True,
    eigvals_only: onp.ToFalse = False,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    type: _EigHType = 1,
    check_finite: op.CanBool = True,
    subset_by_index: _EigHSubsetByIndex | None = None,
    subset_by_value: _EigHSubsetByValue | None = None,
    driver: _DriverEV | _DriverGV | None = None,
) -> tuple[_FloatND, _InexactND]: ...
@overload  # +complex128, eigvals_only: True
def eigh(
    a: onp.ToArrayND[float, npc.inexact80 | npc.number64 | npc.integer32],
    b: onp.ToComplex128_ND | None = None,
    *,
    lower: op.CanBool = True,
    eigvals_only: onp.ToTrue,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    type: _EigHType = 1,
    check_finite: op.CanBool = True,
    subset_by_index: _EigHSubsetByIndex | None = None,
    subset_by_value: _EigHSubsetByValue | None = None,
    driver: _DriverEV | _EigHSubsetByValue | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload  # +complex, eigvals_only: True
def eigh(
    a: onp.ToComplexND,
    b: onp.ToComplexND | None = None,
    *,
    lower: op.CanBool = True,
    eigvals_only: onp.ToTrue,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    type: _EigHType = 1,
    check_finite: op.CanBool = True,
    subset_by_index: _EigHSubsetByIndex | None = None,
    subset_by_value: _EigHSubsetByValue | None = None,
    driver: _DriverEV | _EigHSubsetByValue | None = None,
) -> _FloatND: ...

#
@overload  # float, eigvals_only: False = ..., select: _SelectA = ...
def eig_banded(
    a_band: onp.ToFloatND,
    lower: op.CanBool = False,
    eigvals_only: onp.ToFalse = False,
    overwrite_a_band: op.CanBool = False,
    select: _SelectA = "a",
    select_range: _SelectRange | None = None,
    max_ev: onp.ToInt = 0,
    check_finite: op.CanBool = True,
) -> tuple[_FloatND, _FloatND]: ...
@overload  # float, eigvals_only: False = ..., select: _SelectV (keyword)
def eig_banded(
    a_band: onp.ToFloatND,
    lower: op.CanBool = False,
    eigvals_only: onp.ToFalse = False,
    overwrite_a_band: op.CanBool = False,
    *,
    select: _SelectV,
    select_range: _SelectRange,
    max_ev: onp.ToInt = 0,
    check_finite: op.CanBool = True,
) -> tuple[_FloatND, _FloatND]: ...
@overload  # float, eigvals_only: False = ..., select: _SelectI (keyword)
def eig_banded(
    a_band: onp.ToFloatND,
    lower: op.CanBool = False,
    eigvals_only: onp.ToFalse = False,
    overwrite_a_band: op.CanBool = False,
    *,
    select: _SelectI,
    select_range: _SelectRangeI,
    max_ev: onp.ToInt = 0,
    check_finite: op.CanBool = True,
) -> tuple[_FloatND, _FloatND]: ...
@overload  # complex, eigvals_only: False = ..., select: _SelectA = ...
def eig_banded(
    a_band: onp.ToComplexND,
    lower: op.CanBool = False,
    eigvals_only: onp.ToFalse = False,
    overwrite_a_band: op.CanBool = False,
    select: _SelectA = "a",
    select_range: _SelectRange | None = None,
    max_ev: onp.ToInt = 0,
    check_finite: op.CanBool = True,
) -> tuple[_FloatND, _InexactND]: ...
@overload  # complex, eigvals_only: False = ..., select: _SelectV (keyword)
def eig_banded(
    a_band: onp.ToComplexND,
    lower: op.CanBool = False,
    eigvals_only: onp.ToFalse = False,
    overwrite_a_band: op.CanBool = False,
    *,
    select: _SelectV,
    select_range: _SelectRange,
    max_ev: onp.ToInt = 0,
    check_finite: op.CanBool = True,
) -> tuple[_FloatND, _InexactND]: ...
@overload  # complex, eigvals_only: False = ..., select: _SelectI (keyword)
def eig_banded(
    a_band: onp.ToComplexND,
    lower: op.CanBool = False,
    eigvals_only: onp.ToFalse = False,
    overwrite_a_band: op.CanBool = False,
    *,
    select: _SelectI,
    select_range: _SelectRangeI,
    max_ev: onp.ToInt = 0,
    check_finite: op.CanBool = True,
) -> tuple[_FloatND, _InexactND]: ...
@overload  # eigvals_only: True  (positional), select: _SelectA = ...
def eig_banded(
    a_band: onp.ToComplexND,
    lower: op.CanBool,
    eigvals_only: onp.ToTrue,
    overwrite_a_band: op.CanBool = False,
    select: _SelectA = "a",
    select_range: _SelectRange | None = None,
    max_ev: onp.ToInt = 0,
    check_finite: op.CanBool = True,
) -> _FloatND: ...
@overload  # eigvals_only: True  (keyword), select: _SelectA = ... (keyword)
def eig_banded(
    a_band: onp.ToComplexND,
    lower: op.CanBool = False,
    *,
    eigvals_only: onp.ToTrue,
    overwrite_a_band: op.CanBool = False,
    select: _SelectA = "a",
    select_range: _SelectRange | None = None,
    max_ev: onp.ToInt = 0,
    check_finite: op.CanBool = True,
) -> _FloatND: ...
@overload  # eigvals_only: True  (positional), select: _SelectV (keyword)
def eig_banded(
    a_band: onp.ToComplexND,
    lower: op.CanBool,
    eigvals_only: onp.ToTrue,
    overwrite_a_band: op.CanBool = False,
    *,
    select: _SelectV,
    select_range: _SelectRange,
    max_ev: onp.ToInt = 0,
    check_finite: op.CanBool = True,
) -> _FloatND: ...
@overload  # eigvals_only: True  (keyword), select: _SelectV (keyword)
def eig_banded(
    a_band: onp.ToComplexND,
    lower: op.CanBool = False,
    *,
    eigvals_only: onp.ToTrue,
    overwrite_a_band: op.CanBool = False,
    select: _SelectV,
    select_range: _SelectRange,
    max_ev: onp.ToInt = 0,
    check_finite: op.CanBool = True,
) -> _FloatND: ...
@overload  # eigvals_only: True (positional), select: _SelectI (keyword)
def eig_banded(
    a_band: onp.ToComplexND,
    lower: op.CanBool,
    eigvals_only: onp.ToTrue,
    overwrite_a_band: op.CanBool = False,
    *,
    select: _SelectI,
    select_range: _SelectRangeI,
    max_ev: onp.ToInt = 0,
    check_finite: op.CanBool = True,
) -> _FloatND: ...
@overload  # eigvals_only: True (keyword), select: _SelectI (keyword)
def eig_banded(
    a_band: onp.ToComplexND,
    lower: op.CanBool = False,
    *,
    eigvals_only: onp.ToTrue,
    overwrite_a_band: op.CanBool = False,
    select: _SelectI,
    select_range: _SelectRangeI,
    max_ev: onp.ToInt = 0,
    check_finite: op.CanBool = True,
) -> _FloatND: ...

#
@overload  # homogeneous_eigvals: False = ...
def eigvals(
    a: onp.ToComplexND,
    b: onp.ToComplexND | None = None,
    overwrite_a: op.CanBool = False,
    check_finite: op.CanBool = True,
    homogeneous_eigvals: onp.ToFalse = False,
) -> _ComplexND: ...
@overload  # homogeneous_eigvals: True (positional)
def eigvals(
    a: onp.ToComplexND,
    b: onp.ToComplexND | None,
    overwrite_a: op.CanBool,
    check_finite: op.CanBool,
    homogeneous_eigvals: onp.ToTrue,
) -> _ComplexND: ...
@overload  # homogeneous_eigvals: True (keyword)
def eigvals(
    a: onp.ToComplexND,
    b: onp.ToComplexND | None = None,
    overwrite_a: op.CanBool = False,
    check_finite: op.CanBool = True,
    *,
    homogeneous_eigvals: onp.ToTrue,
) -> _ComplexND: ...
@overload  # catch-all
def eigvals(
    a: onp.ToComplexND,
    b: onp.ToComplexND | None = None,
    overwrite_a: op.CanBool = False,
    check_finite: op.CanBool = True,
    homogeneous_eigvals: op.CanBool = False,
) -> _ComplexND: ...

#
def eigvalsh(
    a: onp.ToComplexND,
    b: onp.ToComplexND | None = None,
    *,
    lower: op.CanBool = True,
    overwrite_a: op.CanBool = False,
    overwrite_b: op.CanBool = False,
    type: _EigHType = 1,
    check_finite: op.CanBool = True,
    subset_by_index: _EigHSubsetByIndex | None = None,
    subset_by_value: _EigHSubsetByValue | None = None,
    driver: _DriverEV | _EigHSubsetByValue | None = None,
) -> _FloatND: ...

#
@overload  # select: _SelectA = ...
def eigvals_banded(
    a_band: onp.ToComplexND,
    lower: op.CanBool = False,
    overwrite_a_band: op.CanBool = False,
    select: _SelectA = "a",
    select_range: _SelectRange | None = None,
    check_finite: op.CanBool = True,
) -> _FloatND: ...
@overload  # select: _SelectV (positional)
def eigvals_banded(
    a_band: onp.ToComplexND,
    lower: op.CanBool,
    overwrite_a_band: op.CanBool,
    select: _SelectV,
    select_range: _SelectRange,
    check_finite: op.CanBool = True,
) -> _FloatND: ...
@overload  # select: _SelectV (keyword)
def eigvals_banded(
    a_band: onp.ToComplexND,
    lower: op.CanBool = False,
    overwrite_a_band: op.CanBool = False,
    *,
    select: _SelectV,
    select_range: _SelectRange,
    check_finite: op.CanBool = True,
) -> _FloatND: ...
@overload  # select: _SelectI (positional)
def eigvals_banded(
    a_band: onp.ToComplexND,
    lower: op.CanBool,
    overwrite_a_band: op.CanBool,
    select: _SelectI,
    select_range: _SelectRangeI,
    check_finite: op.CanBool = True,
) -> _FloatND: ...
@overload  # select: _SelectI (keyword)
def eigvals_banded(
    a_band: onp.ToComplexND,
    lower: op.CanBool = False,
    overwrite_a_band: op.CanBool = False,
    *,
    select: _SelectI,
    select_range: _SelectRangeI,
    check_finite: op.CanBool = True,
) -> _FloatND: ...

#
@overload  # select: _SelectA = ...
def eigvalsh_tridiagonal(
    d: onp.ToComplexND,
    e: onp.ToComplexND,
    select: _SelectA = "a",
    select_range: _SelectRange | None = None,
    check_finite: op.CanBool = True,
    tol: onp.ToFloat = 0.0,
    lapack_driver: _DriverSTE | _DriverAuto = "auto",
) -> _FloatND: ...
@overload  # select: _SelectV
def eigvalsh_tridiagonal(
    d: onp.ToComplexND,
    e: onp.ToComplexND,
    select: _SelectV,
    select_range: _SelectRange,
    check_finite: op.CanBool = True,
    tol: onp.ToFloat = 0.0,
    lapack_driver: _DriverSTE | _DriverAuto = "auto",
) -> _FloatND: ...
@overload  # select: _SelectI
def eigvalsh_tridiagonal(
    d: onp.ToComplexND,
    e: onp.ToComplexND,
    select: _SelectI,
    select_range: _SelectRangeI,
    check_finite: op.CanBool = True,
    tol: onp.ToFloat = 0.0,
    lapack_driver: _DriverSTE | _DriverAuto = "auto",
) -> _FloatND: ...

#
@overload  # eigvals_only: False = ..., select: _SelectA = ...
def eigh_tridiagonal(
    d: onp.ToFloatND,
    e: onp.ToFloatND,
    eigvals_only: onp.ToFalse = False,
    select: _SelectA = "a",
    select_range: _SelectRange | None = None,
    check_finite: op.CanBool = True,
    tol: onp.ToFloat = 0.0,
    lapack_driver: _DriverSTE | _DriverAuto = "auto",
) -> tuple[_FloatND, _FloatND]: ...
@overload  # eigvals_only: False, select: _SelectV (positional)
def eigh_tridiagonal(
    d: onp.ToFloatND,
    e: onp.ToFloatND,
    eigvals_only: onp.ToFalse,
    select: _SelectV,
    select_range: _SelectRange,
    check_finite: op.CanBool = True,
    tol: onp.ToFloat = 0.0,
    lapack_driver: _DriverSTE | _DriverAuto = "auto",
) -> tuple[_FloatND, _FloatND]: ...
@overload  # eigvals_only: False = ..., select: _SelectV (keyword)
def eigh_tridiagonal(
    d: onp.ToFloatND,
    e: onp.ToFloatND,
    eigvals_only: onp.ToFalse = False,
    *,
    select: _SelectV,
    select_range: _SelectRange,
    check_finite: op.CanBool = True,
    tol: onp.ToFloat = 0.0,
    lapack_driver: _DriverSTE | _DriverAuto = "auto",
) -> tuple[_FloatND, _FloatND]: ...
@overload  # eigvals_only: False, select: _SelectI (positional)
def eigh_tridiagonal(
    d: onp.ToFloatND,
    e: onp.ToFloatND,
    eigvals_only: onp.ToFalse,
    select: _SelectI,
    select_range: _SelectRangeI,
    check_finite: op.CanBool = True,
    tol: onp.ToFloat = 0.0,
    lapack_driver: _DriverSTE | _DriverAuto = "auto",
) -> tuple[_FloatND, _FloatND]: ...
@overload  # eigvals_only: False = ..., select: _SelectI (keyword)
def eigh_tridiagonal(
    d: onp.ToFloatND,
    e: onp.ToFloatND,
    eigvals_only: onp.ToFalse = False,
    *,
    select: _SelectI,
    select_range: _SelectRangeI,
    check_finite: op.CanBool = True,
    tol: onp.ToFloat = 0.0,
    lapack_driver: _DriverSTE | _DriverAuto = "auto",
) -> tuple[_FloatND, _FloatND]: ...
@overload  # eigvals_only: True, select: _SelectA = ...
def eigh_tridiagonal(
    d: onp.ToFloatND,
    e: onp.ToFloatND,
    eigvals_only: onp.ToTrue,
    select: _SelectA = "a",
    select_range: _SelectRange | None = None,
    check_finite: op.CanBool = True,
    tol: onp.ToFloat = 0.0,
    lapack_driver: _DriverSTE | _DriverAuto = "auto",
) -> _FloatND: ...
@overload  # eigvals_only: True, select: _SelectV
def eigh_tridiagonal(
    d: onp.ToFloatND,
    e: onp.ToFloatND,
    eigvals_only: onp.ToTrue,
    select: _SelectV,
    select_range: _SelectRange,
    check_finite: op.CanBool = True,
    tol: onp.ToFloat = 0.0,
    lapack_driver: _DriverSTE | _DriverAuto = "auto",
) -> _FloatND: ...
@overload  # eigvals_only: True, select: _SelectI
def eigh_tridiagonal(
    d: onp.ToFloatND,
    e: onp.ToFloatND,
    eigvals_only: onp.ToTrue,
    select: _SelectI,
    select_range: _SelectRangeI,
    check_finite: op.CanBool = True,
    tol: onp.ToFloat = 0.0,
    lapack_driver: _DriverSTE | _DriverAuto = "auto",
) -> _FloatND: ...

#
@overload  # float, calc_q: False = ...
def hessenberg(
    a: onp.ToFloatND, calc_q: onp.ToFalse = False, overwrite_a: op.CanBool = False, check_finite: op.CanBool = True
) -> _FloatND: ...
@overload  # float, calc_q: True
def hessenberg(
    a: onp.ToFloatND, calc_q: onp.ToTrue, overwrite_a: op.CanBool = False, check_finite: op.CanBool = True
) -> tuple[_FloatND, _FloatND]: ...
@overload  # complex, calc_q: False = ...
def hessenberg(
    a: onp.ToComplexND, calc_q: onp.ToFalse = False, overwrite_a: op.CanBool = False, check_finite: op.CanBool = True
) -> _InexactND: ...
@overload  # complex, calc_q: True
def hessenberg(
    a: onp.ToComplexND, calc_q: onp.ToTrue, overwrite_a: op.CanBool = False, check_finite: op.CanBool = True
) -> tuple[_InexactND, _InexactND]: ...
@overload  # complex, calc_q: CanBool (catch-all)
def hessenberg(
    a: onp.ToComplexND, calc_q: op.CanBool, overwrite_a: op.CanBool = False, check_finite: op.CanBool = True
) -> _InexactND | tuple[_InexactND, _InexactND]: ...

#
@overload
def cdf2rdf(w: _ArrayLike[_FloatT], v: _ArrayLike[_FloatT2]) -> tuple[onp.ArrayND[_FloatT], onp.ArrayND[_FloatT2]]: ...
@overload
def cdf2rdf(w: _ArrayLike[_FloatT], v: onp.ToComplexND) -> tuple[onp.ArrayND[_FloatT], _FloatND]: ...
@overload
def cdf2rdf(w: onp.ToComplexND, v: _ArrayLike[_FloatT2]) -> tuple[_FloatND, onp.ArrayND[_FloatT2]]: ...
@overload
def cdf2rdf(w: onp.ToComplexND, v: onp.ToComplexND) -> tuple[_FloatND, _FloatND]: ...
