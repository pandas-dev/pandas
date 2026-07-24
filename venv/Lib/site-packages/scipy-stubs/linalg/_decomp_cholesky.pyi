from collections.abc import MutableSequence, Sequence
from typing import Any, overload

import numpy as np
import numpy_typing_compat as nptc
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["cho_factor", "cho_solve", "cho_solve_banded", "cholesky", "cholesky_banded"]

###

type _as_f32 = np.float32 | np.float16 | npc.integer16 | npc.integer8 | np.bool  # noqa: PYI042
type _as_f64 = npc.floating64 | npc.floating80 | npc.integer64 | npc.integer32  # noqa: PYI042
type _as_c64 = np.complex64  # noqa: PYI042
type _as_c128 = npc.complexfloating160 | npc.complexfloating128  # noqa: PYI042

type _Sequence2D[T] = Sequence[Sequence[T]]

###

# NOTE: The ignored `overload-overlap` mypy errors are false positives

# keep in sync with `cholesky_banded` and `cho_factor`
@overload  # Nd +f64
def cholesky[Shape2T: tuple[int, int, *tuple[int, ...]]](  # type: ignore[overload-overlap]
    a: nptc.CanArray[Shape2T, np.dtype[_as_f64]], lower: bool = False, overwrite_a: bool = False, check_finite: bool = True
) -> onp.ArrayND[np.float64, Shape2T]: ...
@overload  # Nd +f32
def cholesky[Shape2T: tuple[int, int, *tuple[int, ...]]](
    a: nptc.CanArray[Shape2T, np.dtype[_as_f32]], lower: bool = False, overwrite_a: bool = False, check_finite: bool = True
) -> onp.ArrayND[np.float32, Shape2T]: ...
@overload  # Nd +c128
def cholesky[Shape2T: tuple[int, int, *tuple[int, ...]]](
    a: nptc.CanArray[Shape2T, np.dtype[_as_c128]], lower: bool = False, overwrite_a: bool = False, check_finite: bool = True
) -> onp.ArrayND[np.complex128, Shape2T]: ...
@overload  # Nd ~c64
def cholesky[Shape2T: tuple[int, int, *tuple[int, ...]]](
    a: nptc.CanArray[Shape2T, np.dtype[_as_c64]], lower: bool = False, overwrite_a: bool = False, check_finite: bool = True
) -> onp.ArrayND[np.complex64, Shape2T]: ...
@overload  # Nd ~number
def cholesky[Shape2T: tuple[int, int, *tuple[int, ...]]](
    a: nptc.CanArray[Shape2T, np.dtype[npc.number]], lower: bool = False, overwrite_a: bool = False, check_finite: bool = True
) -> onp.ArrayND[Any, Shape2T]: ...
@overload  # 2d +f64
def cholesky(
    a: _Sequence2D[float], lower: bool = False, overwrite_a: bool = False, check_finite: bool = True
) -> onp.Array2D[np.float64]: ...
@overload  # 2d ~c128
def cholesky(
    a: Sequence[MutableSequence[complex]], lower: bool = False, overwrite_a: bool = False, check_finite: bool = True
) -> onp.Array2D[np.complex128]: ...
@overload  # ?d +f64
def cholesky(
    a: onp.SequenceND[_Sequence2D[float]], lower: bool = False, overwrite_a: bool = False, check_finite: bool = True
) -> onp.ArrayND[np.float64]: ...
@overload  # ?d ~c128
def cholesky(
    a: onp.SequenceND[Sequence[MutableSequence[complex]]],
    lower: bool = False,
    overwrite_a: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.complex128]: ...

# keep in sync with `cholesky` (but swap `lower` and `overwrite_*`)
@overload  # Nd +f64
def cholesky_banded[Shape2T: tuple[int, int, *tuple[int, ...]]](  # type: ignore[overload-overlap]
    ab: nptc.CanArray[Shape2T, np.dtype[_as_f64]], overwrite_ab: bool = False, lower: bool = False, check_finite: bool = True
) -> onp.ArrayND[np.float64, Shape2T]: ...
@overload  # Nd +f32
def cholesky_banded[Shape2T: tuple[int, int, *tuple[int, ...]]](
    ab: nptc.CanArray[Shape2T, np.dtype[_as_f32]], overwrite_ab: bool = False, lower: bool = False, check_finite: bool = True
) -> onp.ArrayND[np.float32, Shape2T]: ...
@overload  # Nd +c128
def cholesky_banded[Shape2T: tuple[int, int, *tuple[int, ...]]](
    ab: nptc.CanArray[Shape2T, np.dtype[_as_c128]], overwrite_ab: bool = False, lower: bool = False, check_finite: bool = True
) -> onp.ArrayND[np.complex128, Shape2T]: ...
@overload  # Nd ~c64
def cholesky_banded[Shape2T: tuple[int, int, *tuple[int, ...]]](
    ab: nptc.CanArray[Shape2T, np.dtype[_as_c64]], overwrite_ab: bool = False, lower: bool = False, check_finite: bool = True
) -> onp.ArrayND[np.complex64, Shape2T]: ...
@overload  # Nd ~number
def cholesky_banded[Shape2T: tuple[int, int, *tuple[int, ...]]](
    ab: nptc.CanArray[Shape2T, np.dtype[npc.number]], overwrite_ab: bool = False, lower: bool = False, check_finite: bool = True
) -> onp.ArrayND[Any, Shape2T]: ...
@overload  # 2d +f64
def cholesky_banded(
    ab: _Sequence2D[float], overwrite_ab: bool = False, lower: bool = False, check_finite: bool = True
) -> onp.Array2D[np.float64]: ...
@overload  # 2d ~c128
def cholesky_banded(
    ab: Sequence[MutableSequence[complex]], overwrite_ab: bool = False, lower: bool = False, check_finite: bool = True
) -> onp.Array2D[np.complex128]: ...
@overload  # ?d +f64
def cholesky_banded(
    ab: onp.SequenceND[_Sequence2D[float]], overwrite_ab: bool = False, lower: bool = False, check_finite: bool = True
) -> onp.ArrayND[np.float64]: ...
@overload  # ?d ~c128
def cholesky_banded(
    ab: onp.SequenceND[Sequence[MutableSequence[complex]]],
    overwrite_ab: bool = False,
    lower: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.complex128]: ...

# keep in sync with `cholesky`
@overload  # Nd +f64
def cho_factor[Shape2T: tuple[int, int, *tuple[int, ...]]](  # type: ignore[overload-overlap]
    a: nptc.CanArray[Shape2T, np.dtype[_as_f64]], lower: bool = False, overwrite_a: bool = False, check_finite: bool = True
) -> tuple[onp.ArrayND[np.float64, Shape2T], onp.ArrayND[np.bool]]: ...
@overload  # Nd +f32
def cho_factor[Shape2T: tuple[int, int, *tuple[int, ...]]](
    a: nptc.CanArray[Shape2T, np.dtype[_as_f32]], lower: bool = False, overwrite_a: bool = False, check_finite: bool = True
) -> tuple[onp.ArrayND[np.float32, Shape2T], onp.ArrayND[np.bool]]: ...
@overload  # Nd +c128
def cho_factor[Shape2T: tuple[int, int, *tuple[int, ...]]](
    a: nptc.CanArray[Shape2T, np.dtype[_as_c128]], lower: bool = False, overwrite_a: bool = False, check_finite: bool = True
) -> tuple[onp.ArrayND[np.complex128, Shape2T], onp.ArrayND[np.bool]]: ...
@overload  # Nd ~c64
def cho_factor[Shape2T: tuple[int, int, *tuple[int, ...]]](
    a: nptc.CanArray[Shape2T, np.dtype[_as_c64]], lower: bool = False, overwrite_a: bool = False, check_finite: bool = True
) -> tuple[onp.ArrayND[np.complex64, Shape2T], onp.ArrayND[np.bool]]: ...
@overload  # Nd ~number
def cho_factor[Shape2T: tuple[int, int, *tuple[int, ...]]](
    a: nptc.CanArray[Shape2T, np.dtype[npc.number]], lower: bool = False, overwrite_a: bool = False, check_finite: bool = True
) -> tuple[onp.ArrayND[Any, Shape2T], onp.ArrayND[np.bool]]: ...
@overload  # 2d +f64
def cho_factor(
    a: _Sequence2D[float], lower: bool = False, overwrite_a: bool = False, check_finite: bool = True
) -> tuple[onp.Array2D[np.float64], onp.ArrayND[np.bool]]: ...
@overload  # 2d ~c128
def cho_factor(
    a: Sequence[MutableSequence[complex]], lower: bool = False, overwrite_a: bool = False, check_finite: bool = True
) -> tuple[onp.Array2D[np.complex128], onp.ArrayND[np.bool]]: ...
@overload  # ?d +f64
def cho_factor(
    a: onp.SequenceND[_Sequence2D[float]], lower: bool = False, overwrite_a: bool = False, check_finite: bool = True
) -> tuple[onp.ArrayND[np.float64], onp.ArrayND[np.bool]]: ...
@overload  # ?d ~c128
def cho_factor(
    a: onp.SequenceND[Sequence[MutableSequence[complex]]],
    lower: bool = False,
    overwrite_a: bool = False,
    check_finite: bool = True,
) -> tuple[onp.ArrayND[np.complex128], onp.ArrayND[np.bool]]: ...

# keep in sync with `cho_solve_banded` and `lu_solve` in `_decomp_lu`
@overload  # ?d +f64\+f32, ?d +f64
def cho_solve(  # type: ignore[overload-overlap]
    c_and_lower: tuple[onp.ToArrayND[float, _as_f64], onp.ToBool | onp.ToBoolND],
    b: onp.ToFloatND,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.float64]: ...
@overload  # ?d +f64, ?d +f64\+f32
def cho_solve(  # type: ignore[overload-overlap]
    c_and_lower: tuple[onp.ToFloatND, onp.ToBool | onp.ToBoolND],
    b: onp.ToArrayND[float, _as_f64],
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.float64]: ...
@overload  # ?d +f32, ?d +f32
def cho_solve(
    c_and_lower: tuple[onp.ToFloat32_ND, onp.ToBool | onp.ToBoolND],
    b: onp.ToFloat32_ND,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.float32]: ...
@overload  # ?d ~c128|c160, ?d +c128
def cho_solve(
    c_and_lower: tuple[onp.ToJustComplex128_ND | onp.ToJustCLongDoubleND, onp.ToBool | onp.ToBoolND],
    b: onp.ToComplexND,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.complex128]: ...
@overload  # ?d +c128, ?d ~c128|c160
def cho_solve(
    c_and_lower: tuple[onp.ToComplexND, onp.ToBool | onp.ToBoolND],
    b: onp.ToJustComplex128_ND | onp.ToJustCLongDoubleND,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.complex128]: ...
@overload  # ?d ~c64, ?d +c64
def cho_solve(
    c_and_lower: tuple[onp.ToJustComplex64_ND, onp.ToBool | onp.ToBoolND],
    b: onp.ToComplex64_ND,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.complex64]: ...
@overload  # ?d +c64, ?d ~c64
def cho_solve(
    c_and_lower: tuple[onp.ToComplex64_ND, onp.ToBool | onp.ToBoolND],
    b: onp.ToJustComplex64_ND,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.complex64]: ...
@overload  # ?d +cfloating, ?d ~cfloating (fallback)
def cho_solve(
    c_and_lower: tuple[onp.ToComplexND, onp.ToBool | onp.ToBoolND],
    b: onp.ToComplexND,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[Any]: ...

# keep in sync with `cho_solve`
@overload  # ?d +f64\+f32, ?d +f64
def cho_solve_banded(  # type: ignore[overload-overlap]
    cb_and_lower: tuple[onp.ToArrayND[float, _as_f64], onp.ToBool | onp.ToBoolND],
    b: onp.ToFloatND,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.float64]: ...
@overload  # ?d +f64, ?d +f64\+f32
def cho_solve_banded(  # type: ignore[overload-overlap]
    cb_and_lower: tuple[onp.ToFloatND, onp.ToBool | onp.ToBoolND],
    b: onp.ToArrayND[float, _as_f64],
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.float64]: ...
@overload  # ?d +f32, ?d +f32
def cho_solve_banded(
    cb_and_lower: tuple[onp.ToFloat32_ND, onp.ToBool | onp.ToBoolND],
    b: onp.ToFloat32_ND,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.float32]: ...
@overload  # ?d ~c128|c160, ?d +c128
def cho_solve_banded(
    cb_and_lower: tuple[onp.ToJustComplex128_ND | onp.ToJustCLongDoubleND, onp.ToBool | onp.ToBoolND],
    b: onp.ToComplexND,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.complex128]: ...
@overload  # ?d +c128, ?d ~c128|c160
def cho_solve_banded(
    cb_and_lower: tuple[onp.ToComplexND, onp.ToBool | onp.ToBoolND],
    b: onp.ToJustComplex128_ND | onp.ToJustCLongDoubleND,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.complex128]: ...
@overload  # ?d ~c64, ?d +c64
def cho_solve_banded(
    cb_and_lower: tuple[onp.ToJustComplex64_ND, onp.ToBool | onp.ToBoolND],
    b: onp.ToComplex64_ND,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.complex64]: ...
@overload  # ?d +c64, ?d ~c64
def cho_solve_banded(
    cb_and_lower: tuple[onp.ToComplex64_ND, onp.ToBool | onp.ToBoolND],
    b: onp.ToJustComplex64_ND,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.complex64]: ...
@overload  # ?d +cfloating, ?d ~cfloating (fallback)
def cho_solve_banded(
    cb_and_lower: tuple[onp.ToComplexND, onp.ToBool | onp.ToBoolND],
    b: onp.ToComplexND,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[Any]: ...
