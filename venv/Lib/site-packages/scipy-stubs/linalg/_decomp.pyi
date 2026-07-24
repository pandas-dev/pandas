from collections.abc import Iterable, Sequence
from typing import Any, Literal, overload

import numpy as np
import numpy_typing_compat as nptc
import optype.numpy as onp
import optype.numpy.compat as npc
import optype.typing as opt

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

# input types
# NOTE: only "a", "v" and "i" are documented for the `select` params, but internally 0, 1, and 2 are used, respectively.
type _SelectA = Literal["a", "all", 0]
type _SelectV = Literal["v", "value", 1]
type _SelectI = Literal["i", "index", 2]

# NOTE: `_check_select()` requires the `select_range` array-like to be of `int{16,32,64}` when `select: _SelectIndex`
# https://github.com/scipy/scipy-stubs/issues/154
# NOTE: This `select_range` parameter type must be of shape `(2,)` and in nondescending order
type _SelectRange = Sequence[float | npc.integer | npc.floating]
type _SelectRangeI = Sequence[int | np.int16 | np.int32 | np.int64]  # no bool, int8 or unsigned ints

type _EigHType = Literal[1, 2, 3]
type _EigHSubsetByIndex = Iterable[opt.AnyInt]
type _EigHSubsetByValue = Iterable[onp.ToFloat]

# LAPACK drivers
type _DriverGV = Literal["gv", "gvd", "gvx"]
type _DriverEV = Literal["ev", "evd", "evx", "evr"]
type _DriverSTE = Literal["stemr", "stebz", "sterf", "stev"]
type _DriverAuto = Literal["auto"]

# output types
type _FloatND = onp.ArrayND[np.float64 | np.float32]
type _ComplexND = onp.ArrayND[np.complex128 | np.complex64]
type _InexactND = onp.ArrayND[np.complex128 | np.complex64 | np.float64 | np.float32]

###

# NOTE: The eigenvectors of real `a` can be either real or complex, depending on its values.
# TODO(@jorenham): f32/f64/c64/c128-specific overloads
@overload  # complex, left: False, right: False (positional)
def eig(
    a: onp.ToComplexND,
    b: onp.ToComplexND | None,
    left: Literal[False],
    right: Literal[False],
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    homogeneous_eigvals: bool = False,
) -> _ComplexND: ...
@overload  # complex, left: False = ..., right: False (keyword)
def eig(
    a: onp.ToComplexND,
    b: onp.ToComplexND | None = None,
    left: Literal[False] = False,
    *,
    right: Literal[False],
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    homogeneous_eigvals: bool = False,
) -> _ComplexND: ...
@overload  # float, left: False = ..., right: True = ...
def eig(
    a: onp.ToFloatND,
    b: onp.ToFloatND | None = None,
    left: Literal[False] = False,
    right: Literal[True] = True,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    homogeneous_eigvals: bool = False,
) -> tuple[_ComplexND, _InexactND]: ...
@overload  # complex, left: False = ..., right: True = ...
def eig(
    a: onp.ToJustComplexND,
    b: onp.ToComplexND | None = None,
    left: Literal[False] = False,
    right: Literal[True] = True,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    homogeneous_eigvals: bool = False,
) -> tuple[_ComplexND, _ComplexND]: ...
@overload  # float, left: True (positional), right: False
def eig(
    a: onp.ToFloatND,
    b: onp.ToFloatND | None,
    left: Literal[True],
    right: Literal[False],
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    homogeneous_eigvals: bool = False,
) -> tuple[_ComplexND, _InexactND]: ...
@overload  # complex, left: True (positional), right: False
def eig(
    a: onp.ToJustComplexND,
    b: onp.ToComplexND | None,
    left: Literal[True],
    right: Literal[False],
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    homogeneous_eigvals: bool = False,
) -> tuple[_ComplexND, _ComplexND]: ...
@overload  # float, left: True (keyword), right: False
def eig(
    a: onp.ToFloatND,
    b: onp.ToFloatND | None = None,
    *,
    left: Literal[True],
    right: Literal[False],
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    homogeneous_eigvals: bool = False,
) -> tuple[_ComplexND, _InexactND]: ...
@overload  # complex, left: True (keyword), right: False (keyword)
def eig(
    a: onp.ToJustComplexND,
    b: onp.ToComplexND | None = None,
    *,
    left: Literal[True],
    right: Literal[False],
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    homogeneous_eigvals: bool = False,
) -> tuple[_ComplexND, _ComplexND]: ...
@overload  # float, left: True (positional), right: True = ...
def eig(
    a: onp.ToFloatND,
    b: onp.ToFloatND | None,
    left: Literal[True],
    right: Literal[True] = True,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    homogeneous_eigvals: bool = False,
) -> tuple[_ComplexND, _InexactND, _InexactND]: ...
@overload  # complex, left: True (positional), right: True = ...
def eig(
    a: onp.ToJustComplexND,
    b: onp.ToComplexND | None,
    left: Literal[True],
    right: Literal[True] = True,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    homogeneous_eigvals: bool = False,
) -> tuple[_ComplexND, _ComplexND, _ComplexND]: ...
@overload  # float, left: True (keyword), right: True = ...
def eig(
    a: onp.ToFloatND,
    b: onp.ToFloatND | None = None,
    *,
    left: Literal[True],
    right: Literal[True] = True,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    homogeneous_eigvals: bool = False,
) -> tuple[_ComplexND, _InexactND, _InexactND]: ...
@overload  # complex, left: True (keyword), right: True = ...
def eig(
    a: onp.ToJustComplexND,
    b: onp.ToComplexND | None = None,
    *,
    left: Literal[True],
    right: Literal[True] = True,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    homogeneous_eigvals: bool = False,
) -> tuple[_ComplexND, _ComplexND, _ComplexND]: ...
@overload  # catch-all
def eig(
    a: onp.ToComplexND,
    b: onp.ToComplexND | None = None,
    left: bool = False,
    right: bool = True,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    homogeneous_eigvals: bool = False,
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
    lower: bool = True,
    eigvals_only: Literal[False] = False,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    type: _EigHType = 1,
    check_finite: bool = True,
    subset_by_index: _EigHSubsetByIndex | None = None,
    subset_by_value: _EigHSubsetByValue | None = None,
    driver: _DriverEV | _DriverGV | None = None,
) -> tuple[onp.ArrayND[np.float64], onp.ArrayND[np.float64]]: ...
@overload  # +float, eigvals_only: False = ...
def eigh(
    a: onp.ToFloatND,
    b: onp.ToFloatND | None = None,
    *,
    lower: bool = True,
    eigvals_only: Literal[False] = False,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    type: _EigHType = 1,
    check_finite: bool = True,
    subset_by_index: _EigHSubsetByIndex | None = None,
    subset_by_value: _EigHSubsetByValue | None = None,
    driver: _DriverEV | _DriverGV | None = None,
) -> tuple[_FloatND, _FloatND]: ...
@overload  # ~complex, eigvals_only: False = ...
def eigh(
    a: onp.ToJustComplexND,
    b: onp.ToComplexND | None = None,
    *,
    lower: bool = True,
    eigvals_only: Literal[False] = False,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    type: _EigHType = 1,
    check_finite: bool = True,
    subset_by_index: _EigHSubsetByIndex | None = None,
    subset_by_value: _EigHSubsetByValue | None = None,
    driver: _DriverEV | _DriverGV | None = None,
) -> tuple[_FloatND, _ComplexND]: ...
@overload  # +complex, eigvals_only: False = ...
def eigh(
    a: onp.ToComplexND,
    b: onp.ToComplexND | None = None,
    *,
    lower: bool = True,
    eigvals_only: Literal[False] = False,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    type: _EigHType = 1,
    check_finite: bool = True,
    subset_by_index: _EigHSubsetByIndex | None = None,
    subset_by_value: _EigHSubsetByValue | None = None,
    driver: _DriverEV | _DriverGV | None = None,
) -> tuple[_FloatND, _InexactND]: ...
@overload  # +complex128, eigvals_only: True
def eigh(
    a: onp.ToArrayND[float, npc.inexact80 | npc.number64 | npc.integer32],
    b: onp.ToComplex128_ND | None = None,
    *,
    lower: bool = True,
    eigvals_only: Literal[True],
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    type: _EigHType = 1,
    check_finite: bool = True,
    subset_by_index: _EigHSubsetByIndex | None = None,
    subset_by_value: _EigHSubsetByValue | None = None,
    driver: _DriverEV | _EigHSubsetByValue | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload  # +complex, eigvals_only: True
def eigh(
    a: onp.ToComplexND,
    b: onp.ToComplexND | None = None,
    *,
    lower: bool = True,
    eigvals_only: Literal[True],
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    type: _EigHType = 1,
    check_finite: bool = True,
    subset_by_index: _EigHSubsetByIndex | None = None,
    subset_by_value: _EigHSubsetByValue | None = None,
    driver: _DriverEV | _EigHSubsetByValue | None = None,
) -> _FloatND: ...

#
@overload  # float, eigvals_only: False = ..., select: _SelectA = ...
def eig_banded(
    a_band: onp.ToFloatND,
    lower: bool = False,
    eigvals_only: Literal[False] = False,
    overwrite_a_band: bool = False,
    select: _SelectA = "a",
    select_range: _SelectRange | None = None,
    max_ev: onp.ToInt = 0,
    check_finite: bool = True,
) -> tuple[_FloatND, _FloatND]: ...
@overload  # float, eigvals_only: False = ..., select: _SelectV (keyword)
def eig_banded(
    a_band: onp.ToFloatND,
    lower: bool = False,
    eigvals_only: Literal[False] = False,
    overwrite_a_band: bool = False,
    *,
    select: _SelectV,
    select_range: _SelectRange,
    max_ev: onp.ToInt = 0,
    check_finite: bool = True,
) -> tuple[_FloatND, _FloatND]: ...
@overload  # float, eigvals_only: False = ..., select: _SelectI (keyword)
def eig_banded(
    a_band: onp.ToFloatND,
    lower: bool = False,
    eigvals_only: Literal[False] = False,
    overwrite_a_band: bool = False,
    *,
    select: _SelectI,
    select_range: _SelectRangeI,
    max_ev: onp.ToInt = 0,
    check_finite: bool = True,
) -> tuple[_FloatND, _FloatND]: ...
@overload  # complex, eigvals_only: False = ..., select: _SelectA = ...
def eig_banded(
    a_band: onp.ToComplexND,
    lower: bool = False,
    eigvals_only: Literal[False] = False,
    overwrite_a_band: bool = False,
    select: _SelectA = "a",
    select_range: _SelectRange | None = None,
    max_ev: onp.ToInt = 0,
    check_finite: bool = True,
) -> tuple[_FloatND, _InexactND]: ...
@overload  # complex, eigvals_only: False = ..., select: _SelectV (keyword)
def eig_banded(
    a_band: onp.ToComplexND,
    lower: bool = False,
    eigvals_only: Literal[False] = False,
    overwrite_a_band: bool = False,
    *,
    select: _SelectV,
    select_range: _SelectRange,
    max_ev: onp.ToInt = 0,
    check_finite: bool = True,
) -> tuple[_FloatND, _InexactND]: ...
@overload  # complex, eigvals_only: False = ..., select: _SelectI (keyword)
def eig_banded(
    a_band: onp.ToComplexND,
    lower: bool = False,
    eigvals_only: Literal[False] = False,
    overwrite_a_band: bool = False,
    *,
    select: _SelectI,
    select_range: _SelectRangeI,
    max_ev: onp.ToInt = 0,
    check_finite: bool = True,
) -> tuple[_FloatND, _InexactND]: ...
@overload  # eigvals_only: True  (positional), select: _SelectA = ...
def eig_banded(
    a_band: onp.ToComplexND,
    lower: bool,
    eigvals_only: Literal[True],
    overwrite_a_band: bool = False,
    select: _SelectA = "a",
    select_range: _SelectRange | None = None,
    max_ev: onp.ToInt = 0,
    check_finite: bool = True,
) -> _FloatND: ...
@overload  # eigvals_only: True  (keyword), select: _SelectA = ... (keyword)
def eig_banded(
    a_band: onp.ToComplexND,
    lower: bool = False,
    *,
    eigvals_only: Literal[True],
    overwrite_a_band: bool = False,
    select: _SelectA = "a",
    select_range: _SelectRange | None = None,
    max_ev: onp.ToInt = 0,
    check_finite: bool = True,
) -> _FloatND: ...
@overload  # eigvals_only: True  (positional), select: _SelectV (keyword)
def eig_banded(
    a_band: onp.ToComplexND,
    lower: bool,
    eigvals_only: Literal[True],
    overwrite_a_band: bool = False,
    *,
    select: _SelectV,
    select_range: _SelectRange,
    max_ev: onp.ToInt = 0,
    check_finite: bool = True,
) -> _FloatND: ...
@overload  # eigvals_only: True  (keyword), select: _SelectV (keyword)
def eig_banded(
    a_band: onp.ToComplexND,
    lower: bool = False,
    *,
    eigvals_only: Literal[True],
    overwrite_a_band: bool = False,
    select: _SelectV,
    select_range: _SelectRange,
    max_ev: onp.ToInt = 0,
    check_finite: bool = True,
) -> _FloatND: ...
@overload  # eigvals_only: True (positional), select: _SelectI (keyword)
def eig_banded(
    a_band: onp.ToComplexND,
    lower: bool,
    eigvals_only: Literal[True],
    overwrite_a_band: bool = False,
    *,
    select: _SelectI,
    select_range: _SelectRangeI,
    max_ev: onp.ToInt = 0,
    check_finite: bool = True,
) -> _FloatND: ...
@overload  # eigvals_only: True (keyword), select: _SelectI (keyword)
def eig_banded(
    a_band: onp.ToComplexND,
    lower: bool = False,
    *,
    eigvals_only: Literal[True],
    overwrite_a_band: bool = False,
    select: _SelectI,
    select_range: _SelectRangeI,
    max_ev: onp.ToInt = 0,
    check_finite: bool = True,
) -> _FloatND: ...

#
def eigvals(
    a: onp.ToComplexND,
    b: onp.ToComplexND | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    homogeneous_eigvals: bool = False,
) -> _ComplexND: ...

#
def eigvalsh(
    a: onp.ToComplexND,
    b: onp.ToComplexND | None = None,
    *,
    lower: bool = True,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    type: _EigHType = 1,
    check_finite: bool = True,
    subset_by_index: _EigHSubsetByIndex | None = None,
    subset_by_value: _EigHSubsetByValue | None = None,
    driver: _DriverEV | _EigHSubsetByValue | None = None,
) -> _FloatND: ...

#
@overload  # select: _SelectA = ...
def eigvals_banded(
    a_band: onp.ToComplexND,
    lower: bool = False,
    overwrite_a_band: bool = False,
    select: _SelectA = "a",
    select_range: _SelectRange | None = None,
    check_finite: bool = True,
) -> _FloatND: ...
@overload  # select: _SelectV (positional)
def eigvals_banded(
    a_band: onp.ToComplexND,
    lower: bool,
    overwrite_a_band: bool,
    select: _SelectV,
    select_range: _SelectRange,
    check_finite: bool = True,
) -> _FloatND: ...
@overload  # select: _SelectV (keyword)
def eigvals_banded(
    a_band: onp.ToComplexND,
    lower: bool = False,
    overwrite_a_band: bool = False,
    *,
    select: _SelectV,
    select_range: _SelectRange,
    check_finite: bool = True,
) -> _FloatND: ...
@overload  # select: _SelectI (positional)
def eigvals_banded(
    a_band: onp.ToComplexND,
    lower: bool,
    overwrite_a_band: bool,
    select: _SelectI,
    select_range: _SelectRangeI,
    check_finite: bool = True,
) -> _FloatND: ...
@overload  # select: _SelectI (keyword)
def eigvals_banded(
    a_band: onp.ToComplexND,
    lower: bool = False,
    overwrite_a_band: bool = False,
    *,
    select: _SelectI,
    select_range: _SelectRangeI,
    check_finite: bool = True,
) -> _FloatND: ...

#
@overload  # select: _SelectA = ...
def eigvalsh_tridiagonal(
    d: onp.ToComplexND,
    e: onp.ToComplexND,
    select: _SelectA = "a",
    select_range: _SelectRange | None = None,
    check_finite: bool = True,
    tol: onp.ToFloat = 0.0,
    lapack_driver: _DriverSTE | _DriverAuto = "auto",
) -> _FloatND: ...
@overload  # select: _SelectV
def eigvalsh_tridiagonal(
    d: onp.ToComplexND,
    e: onp.ToComplexND,
    select: _SelectV,
    select_range: _SelectRange,
    check_finite: bool = True,
    tol: onp.ToFloat = 0.0,
    lapack_driver: _DriverSTE | _DriverAuto = "auto",
) -> _FloatND: ...
@overload  # select: _SelectI
def eigvalsh_tridiagonal(
    d: onp.ToComplexND,
    e: onp.ToComplexND,
    select: _SelectI,
    select_range: _SelectRangeI,
    check_finite: bool = True,
    tol: onp.ToFloat = 0.0,
    lapack_driver: _DriverSTE | _DriverAuto = "auto",
) -> _FloatND: ...

#
@overload  # eigvals_only: False = ..., select: _SelectA = ...
def eigh_tridiagonal(
    d: onp.ToFloatND,
    e: onp.ToFloatND,
    eigvals_only: Literal[False] = False,
    select: _SelectA = "a",
    select_range: _SelectRange | None = None,
    check_finite: bool = True,
    tol: onp.ToFloat = 0.0,
    lapack_driver: _DriverSTE | _DriverAuto = "auto",
) -> tuple[_FloatND, _FloatND]: ...
@overload  # eigvals_only: False, select: _SelectV (positional)
def eigh_tridiagonal(
    d: onp.ToFloatND,
    e: onp.ToFloatND,
    eigvals_only: Literal[False],
    select: _SelectV,
    select_range: _SelectRange,
    check_finite: bool = True,
    tol: onp.ToFloat = 0.0,
    lapack_driver: _DriverSTE | _DriverAuto = "auto",
) -> tuple[_FloatND, _FloatND]: ...
@overload  # eigvals_only: False = ..., select: _SelectV (keyword)
def eigh_tridiagonal(
    d: onp.ToFloatND,
    e: onp.ToFloatND,
    eigvals_only: Literal[False] = False,
    *,
    select: _SelectV,
    select_range: _SelectRange,
    check_finite: bool = True,
    tol: onp.ToFloat = 0.0,
    lapack_driver: _DriverSTE | _DriverAuto = "auto",
) -> tuple[_FloatND, _FloatND]: ...
@overload  # eigvals_only: False, select: _SelectI (positional)
def eigh_tridiagonal(
    d: onp.ToFloatND,
    e: onp.ToFloatND,
    eigvals_only: Literal[False],
    select: _SelectI,
    select_range: _SelectRangeI,
    check_finite: bool = True,
    tol: onp.ToFloat = 0.0,
    lapack_driver: _DriverSTE | _DriverAuto = "auto",
) -> tuple[_FloatND, _FloatND]: ...
@overload  # eigvals_only: False = ..., select: _SelectI (keyword)
def eigh_tridiagonal(
    d: onp.ToFloatND,
    e: onp.ToFloatND,
    eigvals_only: Literal[False] = False,
    *,
    select: _SelectI,
    select_range: _SelectRangeI,
    check_finite: bool = True,
    tol: onp.ToFloat = 0.0,
    lapack_driver: _DriverSTE | _DriverAuto = "auto",
) -> tuple[_FloatND, _FloatND]: ...
@overload  # eigvals_only: True, select: _SelectA = ...
def eigh_tridiagonal(
    d: onp.ToFloatND,
    e: onp.ToFloatND,
    eigvals_only: Literal[True],
    select: _SelectA = "a",
    select_range: _SelectRange | None = None,
    check_finite: bool = True,
    tol: onp.ToFloat = 0.0,
    lapack_driver: _DriverSTE | _DriverAuto = "auto",
) -> _FloatND: ...
@overload  # eigvals_only: True, select: _SelectV
def eigh_tridiagonal(
    d: onp.ToFloatND,
    e: onp.ToFloatND,
    eigvals_only: Literal[True],
    select: _SelectV,
    select_range: _SelectRange,
    check_finite: bool = True,
    tol: onp.ToFloat = 0.0,
    lapack_driver: _DriverSTE | _DriverAuto = "auto",
) -> _FloatND: ...
@overload  # eigvals_only: True, select: _SelectI
def eigh_tridiagonal(
    d: onp.ToFloatND,
    e: onp.ToFloatND,
    eigvals_only: Literal[True],
    select: _SelectI,
    select_range: _SelectRangeI,
    check_finite: bool = True,
    tol: onp.ToFloat = 0.0,
    lapack_driver: _DriverSTE | _DriverAuto = "auto",
) -> _FloatND: ...

#
@overload  # float, calc_q: False = ...
def hessenberg(
    a: onp.ToFloatND, calc_q: Literal[False] = False, overwrite_a: bool = False, check_finite: bool = True
) -> _FloatND: ...
@overload  # float, calc_q: True
def hessenberg(
    a: onp.ToFloatND, calc_q: Literal[True], overwrite_a: bool = False, check_finite: bool = True
) -> tuple[_FloatND, _FloatND]: ...
@overload  # complex, calc_q: False = ...
def hessenberg(
    a: onp.ToComplexND, calc_q: Literal[False] = False, overwrite_a: bool = False, check_finite: bool = True
) -> _InexactND: ...
@overload  # complex, calc_q: True
def hessenberg(
    a: onp.ToComplexND, calc_q: Literal[True], overwrite_a: bool = False, check_finite: bool = True
) -> tuple[_InexactND, _InexactND]: ...
@overload  # complex, calc_q: CanBool (catch-all)
def hessenberg(
    a: onp.ToComplexND, calc_q: bool, overwrite_a: bool = False, check_finite: bool = True
) -> _InexactND | tuple[_InexactND, _InexactND]: ...

#
@overload
def cdf2rdf[FloatVT: npc.floating, FloatWT: npc.floating](
    w: nptc.CanArray[Any, np.dtype[FloatVT]], v: nptc.CanArray[Any, np.dtype[FloatWT]]
) -> tuple[onp.ArrayND[FloatVT], onp.ArrayND[FloatWT]]: ...
@overload
def cdf2rdf[FloatT: npc.floating](
    w: nptc.CanArray[Any, np.dtype[FloatT]], v: onp.ToComplexND
) -> tuple[onp.ArrayND[FloatT], _FloatND]: ...
@overload
def cdf2rdf[FloatT: npc.floating](
    w: onp.ToComplexND, v: nptc.CanArray[Any, np.dtype[FloatT]]
) -> tuple[_FloatND, onp.ArrayND[FloatT]]: ...
@overload
def cdf2rdf(w: onp.ToComplexND, v: onp.ToComplexND) -> tuple[_FloatND, _FloatND]: ...
