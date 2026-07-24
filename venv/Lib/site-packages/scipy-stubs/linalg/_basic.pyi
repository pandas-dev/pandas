# mypy: disable-error-code=overload-overlap

from _typeshed import Incomplete
from collections.abc import Sequence
from typing import Any, Literal, SupportsIndex, overload

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from ._misc import LinAlgError as LinAlgError, LinAlgWarning as LinAlgWarning
from .lapack import get_lapack_funcs as get_lapack_funcs

__all__ = [
    "det",
    "inv",
    "lstsq",
    "matmul_toeplitz",
    "matrix_balance",
    "pinv",
    "pinvh",
    "solve",
    "solve_banded",
    "solve_circulant",
    "solve_toeplitz",
    "solve_triangular",
    "solveh_banded",
]

###

type _Tuple2[T] = tuple[T, T]
type _COrCR[T] = T | _Tuple2[T]

type _Float1D = onp.Array1D[npc.floating]
type _FloatND = onp.ArrayND[npc.floating]

type _Inexact1D = onp.Array1D[npc.inexact]
type _InexactND = onp.ArrayND[npc.inexact]

# TODO(@jorenham): better naming

type _InputFloat = onp.ToArrayND[float, np.float64 | npc.floating80 | npc.integer | np.bool]
type _InputFloatStrict1D = onp.ToArrayStrict1D[float, np.float64 | npc.floating80 | npc.integer | np.bool]
type _InputFloatStrict2D = onp.ToArrayStrict2D[float, np.float64 | npc.floating80 | npc.integer | np.bool]

type _InputF64 = onp.ToArrayND[float, np.float64 | npc.integer | np.bool]
type _InputF64Strict1D = onp.ToArrayStrict1D[float, np.float64 | npc.integer | np.bool]
type _InputF64Strict2D = onp.ToArrayStrict2D[float, np.float64 | npc.integer | np.bool]

type _InputComplex = onp.ToArrayND[op.JustComplex, np.complex128 | npc.complexfloating160]
type _InputComplexStrict1D = onp.ToArrayStrict1D[op.JustComplex, np.complex128 | npc.complexfloating160]
type _InputComplexStrict2D = onp.ToArrayStrict2D[op.JustComplex, np.complex128 | npc.complexfloating160]

type _AssumeA = Literal[
    "diagonal",
    "tridiagonal",
    "banded",
    "upper triangular",
    "lower triangular",
    "symmetric", "sym",
    "hermitian", "her",
    "positive definite", "pos",
    "general", "gen",
]  # fmt: skip
type _TransSystem = Literal[0, "N", 1, "T", 2, "C"]
type _Singular = Literal["lstsq", "raise"]
type _LapackDriver = Literal["gelsd", "gelsy", "gelss"]
type _LapackDriverDS = Literal["gelsd", "gelss"]
type _LapackDriverY = Literal["gelsy"]

# `lstsq` inputs, named after the scalar-type of the result they produce
type _CoFloat64 = np.float64 | npc.floating80 | npc.integer | np.bool

type _AsFloat64Strict1D = onp.ToArrayStrict1D[float, _CoFloat64]
type _AsFloat64Strict2D = onp.ToArrayStrict2D[float, _CoFloat64]
type _AsFloat64_2D = onp.ToArray2D[float, _CoFloat64]
type _AsFloat64_ND = onp.ToArrayND[float, _CoFloat64]

type _AsFloat32Strict1D = onp.ToJustFloat32Strict1D | onp.ToJustFloat16Strict1D
type _AsFloat32Strict2D = onp.ToJustFloat32Strict2D | onp.ToJustFloat16Strict2D
type _AsFloat32_2D = onp.ToJustFloat32_2D | onp.ToJustFloat16_2D
type _AsFloat32_ND = onp.ToJustFloat32_ND | onp.ToJustFloat16_ND

type _AsComplex128Strict1D = onp.ToJustComplex128Strict1D | onp.ToJustCLongDoubleStrict1D
type _AsComplex128Strict2D = onp.ToJustComplex128Strict2D | onp.ToJustCLongDoubleStrict2D
type _AsComplex128_2D = onp.ToJustComplex128_2D | onp.ToJustCLongDouble2D
type _AsComplex128_ND = onp.ToJustComplex128_ND | onp.ToJustCLongDoubleND

type _AsComplex64Strict1D = onp.ToJustComplex64Strict1D | _AsFloat32Strict1D
type _AsComplex64Strict2D = onp.ToJustComplex64Strict2D | _AsFloat32Strict2D
type _AsComplex64_ND = onp.ToJustComplex64_ND | _AsFloat32_ND

type _LstSqResult11[ScalarT: np.generic] = tuple[onp.Array1D[ScalarT], onp.Array1D[ScalarT], np.int64, None]  # gelsy
type _LstSqResult1N[ScalarT: np.generic, T] = tuple[onp.Array1D[ScalarT], onp.ArrayND[ScalarT], np.int64, T]
type _LstSqResult21[ScalarT: np.generic, T] = tuple[onp.Array2D[ScalarT], onp.Array1D[ScalarT], np.int64, T]
type _LstSqResultND[ScalarT: np.generic, T] = tuple[onp.ArrayND[ScalarT], onp.ArrayND[ScalarT], np.int64 | Any, T]

###

@overload  # 2d ~float64, 1d +float64
def solve(
    a: _InputFloatStrict2D,
    b: onp.ToFloatStrict1D,
    lower: bool = False,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    assume_a: _AssumeA | None = None,
    transposed: bool = False,
) -> onp.Array1D[np.float64]: ...
@overload  # 2d ~float64, 2d +float64
def solve(
    a: _InputFloatStrict2D,
    b: onp.ToFloatStrict2D,
    lower: bool = False,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    assume_a: _AssumeA | None = None,
    transposed: bool = False,
) -> onp.Array2D[np.float64]: ...
@overload  # Nd ~float64, Nd +float64
def solve(
    a: _InputFloat,
    b: onp.ToFloatND,
    lower: bool = False,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    assume_a: _AssumeA | None = None,
    transposed: bool = False,
) -> onp.ArrayND[np.float64]: ...
@overload  # 2d +float64, 1d ~float64
def solve(
    a: onp.ToFloatStrict2D,
    b: _InputFloatStrict1D,
    lower: bool = False,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    assume_a: _AssumeA | None = None,
    transposed: bool = False,
) -> onp.Array1D[np.float64]: ...
@overload  # 2d +float64, 2d ~float64
def solve(
    a: onp.ToFloatStrict2D,
    b: _InputFloatStrict2D,
    lower: bool = False,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    assume_a: _AssumeA | None = None,
    transposed: bool = False,
) -> onp.Array2D[np.float64]: ...
@overload  # Nd +float64, Nd ~float64
def solve(
    a: onp.ToFloatND,
    b: _InputFloat,
    lower: bool = False,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    assume_a: _AssumeA | None = None,
    transposed: bool = False,
) -> onp.ArrayND[np.float64]: ...
@overload  # 2d ~complex128, 1d +complex128
def solve(
    a: _InputComplexStrict2D,
    b: onp.ToComplexStrict1D,
    lower: bool = False,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    assume_a: _AssumeA | None = None,
    transposed: bool = False,
) -> onp.Array1D[np.complex128]: ...
@overload  # 2d ~complex128, 2d +complex128
def solve(
    a: _InputComplexStrict2D,
    b: onp.ToComplexStrict2D,
    lower: bool = False,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    assume_a: _AssumeA | None = None,
    transposed: bool = False,
) -> onp.Array2D[np.complex128]: ...
@overload  # Nd ~complex128, Nd +complex128
def solve(
    a: _InputComplex,
    b: onp.ToComplexND,
    lower: bool = False,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    assume_a: _AssumeA | None = None,
    transposed: bool = False,
) -> onp.ArrayND[np.complex128]: ...
@overload  # 2d +complex128, 1d ~complex128
def solve(
    a: onp.ToComplexStrict2D,
    b: _InputComplexStrict1D,
    lower: bool = False,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    assume_a: _AssumeA | None = None,
    transposed: bool = False,
) -> onp.Array1D[np.complex128]: ...
@overload  # 2d +complex128, 2d ~complex128
def solve(
    a: onp.ToComplexStrict2D,
    b: _InputComplexStrict2D,
    lower: bool = False,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    assume_a: _AssumeA | None = None,
    transposed: bool = False,
) -> onp.Array2D[np.complex128]: ...
@overload  # Nd +complex128, Nd ~complex128
def solve(
    a: onp.ToComplexND,
    b: _InputComplex,
    lower: bool = False,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    assume_a: _AssumeA | None = None,
    transposed: bool = False,
) -> onp.ArrayND[np.complex128]: ...
@overload  # 2d +floating, 1d +floating
def solve(
    a: onp.ToFloatStrict2D,
    b: onp.ToFloatStrict1D,
    lower: bool = False,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    assume_a: _AssumeA | None = None,
    transposed: bool = False,
) -> onp.Array1D[np.float32 | np.float64]: ...
@overload  # 2d +floating, 2d +floating
def solve(
    a: onp.ToFloatStrict2D,
    b: onp.ToFloatStrict2D,
    lower: bool = False,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    assume_a: _AssumeA | None = None,
    transposed: bool = False,
) -> onp.Array2D[np.float32 | np.float64]: ...
@overload  # Nd +floating, Nd +floating
def solve(
    a: onp.ToFloatND,
    b: onp.ToFloatND,
    lower: bool = False,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    assume_a: _AssumeA | None = None,
    transposed: bool = False,
) -> onp.ArrayND[np.float32 | np.float64]: ...
@overload  # 2d +complexfloating, 1d ~complexfloating
def solve(
    a: onp.ToComplexStrict2D,
    b: onp.ToJustComplexStrict1D,
    lower: bool = False,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    assume_a: _AssumeA | None = None,
    transposed: bool = False,
) -> onp.Array1D[np.complex64 | np.complex128]: ...
@overload  # 2d +complexfloating, 2d ~complexfloating
def solve(
    a: onp.ToComplexStrict2D,
    b: onp.ToJustComplexStrict2D,
    lower: bool = False,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    assume_a: _AssumeA | None = None,
    transposed: bool = False,
) -> onp.Array2D[np.complex64 | np.complex128]: ...
@overload  # Nd +complexfloating, Nd ~complexfloating
def solve(
    a: onp.ToComplexND,
    b: onp.ToJustComplexND,
    lower: bool = False,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    assume_a: _AssumeA | None = None,
    transposed: bool = False,
) -> onp.ArrayND[np.complex64 | np.complex128]: ...
@overload  # Nd +complexfloating, Nd +complexfloating
def solve(
    a: onp.ToComplexND,
    b: onp.ToComplexND,
    lower: bool = False,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    assume_a: _AssumeA | None = None,
    transposed: bool = False,
) -> onp.ArrayND[np.float32 | np.float64 | np.complex64 | np.complex128]: ...

#
@overload  # 1D ~float64, +float64
def solve_triangular(
    a: _InputFloatStrict2D,
    b: onp.ToFloatStrict1D,
    trans: _TransSystem = 0,
    lower: bool = False,
    unit_diagonal: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.Array1D[np.float64]: ...
@overload  # 2D ~float64, +float64
def solve_triangular(
    a: _InputFloatStrict2D,
    b: onp.ToFloatStrict2D,
    trans: _TransSystem = 0,
    lower: bool = False,
    unit_diagonal: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.Array2D[np.float64]: ...
@overload  # Nd ~float64, +float64
def solve_triangular(
    a: _InputFloat,
    b: onp.ToFloatND,
    trans: _TransSystem = 0,
    lower: bool = False,
    unit_diagonal: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.float64]: ...
@overload  # 1d +float64, ~float64
def solve_triangular(
    a: onp.ToFloatStrict2D,
    b: _InputFloatStrict1D,
    trans: _TransSystem = 0,
    lower: bool = False,
    unit_diagonal: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.Array1D[np.float64]: ...
@overload  # 2d +float64, ~float64
def solve_triangular(
    a: onp.ToFloatStrict2D,
    b: _InputFloatStrict2D,
    trans: _TransSystem = 0,
    lower: bool = False,
    unit_diagonal: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.Array2D[np.float64]: ...
@overload  # Nd +float64, ~float64
def solve_triangular(
    a: onp.ToFloatND,
    b: _InputFloat,
    trans: _TransSystem = 0,
    lower: bool = False,
    unit_diagonal: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.float64]: ...
@overload  # 1d ~complex128, +complex128
def solve_triangular(
    a: _InputComplexStrict2D,
    b: onp.ToComplexStrict1D,
    trans: _TransSystem = 0,
    lower: bool = False,
    unit_diagonal: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.Array1D[np.complex128]: ...
@overload  # 2d ~complex128, +complex128
def solve_triangular(
    a: _InputComplexStrict2D,
    b: onp.ToComplexStrict2D,
    trans: _TransSystem = 0,
    lower: bool = False,
    unit_diagonal: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.Array2D[np.complex128]: ...
@overload  # Nd ~complex128, +complex128
def solve_triangular(
    a: _InputComplex,
    b: onp.ToComplexND,
    trans: _TransSystem = 0,
    lower: bool = False,
    unit_diagonal: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.complex128]: ...
@overload  # 1d +complex128, ~complex128
def solve_triangular(
    a: onp.ToComplexStrict2D,
    b: _InputComplexStrict1D,
    trans: _TransSystem = 0,
    lower: bool = False,
    unit_diagonal: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.Array1D[np.complex128]: ...
@overload  # 2d +complex128, ~complex128
def solve_triangular(
    a: onp.ToComplexStrict2D,
    b: _InputComplexStrict2D,
    trans: _TransSystem = 0,
    lower: bool = False,
    unit_diagonal: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.Array2D[np.complex128]: ...
@overload  # Nd +complex128, ~complex128
def solve_triangular(
    a: onp.ToComplexND,
    b: _InputComplex,
    trans: _TransSystem = 0,
    lower: bool = False,
    unit_diagonal: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.complex128]: ...
@overload  # 1d +floating, +floating
def solve_triangular(
    a: onp.ToFloatStrict2D,
    b: onp.ToFloatStrict1D,
    trans: _TransSystem = 0,
    lower: bool = False,
    unit_diagonal: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.Array1D[np.float32 | np.float64]: ...
@overload  # 2d +floating, +floating
def solve_triangular(
    a: onp.ToFloatStrict2D,
    b: onp.ToFloatStrict2D,
    trans: _TransSystem = 0,
    lower: bool = False,
    unit_diagonal: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.Array2D[np.float32 | np.float64]: ...
@overload  # Nd +floating, +floating
def solve_triangular(
    a: onp.ToFloatND,
    b: onp.ToFloatND,
    trans: _TransSystem = 0,
    lower: bool = False,
    unit_diagonal: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.float32 | np.float64]: ...
@overload  # 1d +complexfloating, ~complexfloating
def solve_triangular(
    a: onp.ToComplexStrict2D,
    b: onp.ToJustComplexStrict1D,
    trans: _TransSystem = 0,
    lower: bool = False,
    unit_diagonal: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.Array1D[np.complex64 | np.complex128]: ...
@overload  # 2d +complexfloating, ~complexfloating
def solve_triangular(
    a: onp.ToComplexStrict2D,
    b: onp.ToJustComplexStrict2D,
    trans: _TransSystem = 0,
    lower: bool = False,
    unit_diagonal: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.Array2D[np.complex64 | np.complex128]: ...
@overload  # Nd +complexfloating, ~complexfloating
def solve_triangular(
    a: onp.ToComplexND,
    b: onp.ToJustComplexND,
    trans: _TransSystem = 0,
    lower: bool = False,
    unit_diagonal: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.complex64 | np.complex128]: ...
@overload  # Nd +complexfloating, +complexfloating
def solve_triangular(
    a: onp.ToComplexND,
    b: onp.ToComplexND,
    trans: _TransSystem = 0,
    lower: bool = False,
    unit_diagonal: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.float32 | np.float64 | np.complex64 | np.complex128]: ...

# NOTE: keep overload structure consistent with `solveh_banded` below
@overload  # 1D ~float64, +float64
def solve_banded(
    l_and_u: tuple[int, int],
    ab: _InputFloatStrict2D,
    b: onp.ToFloatStrict1D,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.Array1D[np.float64]: ...
@overload  # 2D ~float64, +float64
def solve_banded(
    l_and_u: tuple[int, int],
    ab: _InputFloatStrict2D,
    b: onp.ToFloatStrict2D,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.Array2D[np.float64]: ...
@overload  # Nd ~float64, +float64
def solve_banded(
    l_and_u: tuple[int, int],
    ab: _InputFloat,
    b: onp.ToFloatND,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.float64]: ...
@overload  # 1d +float64, ~float64
def solve_banded(
    l_and_u: tuple[int, int],
    ab: onp.ToFloatStrict2D,
    b: _InputFloatStrict1D,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.Array1D[np.float64]: ...
@overload  # 2d +float64, ~float64
def solve_banded(
    l_and_u: tuple[int, int],
    ab: onp.ToFloatStrict2D,
    b: _InputFloatStrict2D,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.Array2D[np.float64]: ...
@overload  # Nd +float64, ~float64
def solve_banded(
    l_and_u: tuple[int, int],
    ab: onp.ToFloatND,
    b: _InputFloat,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.float64]: ...
@overload  # 1d ~complex128, +complex128
def solve_banded(
    l_and_u: tuple[int, int],
    ab: _InputComplexStrict2D,
    b: onp.ToComplexStrict1D,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.Array1D[np.complex128]: ...
@overload  # 2d ~complex128, +complex128
def solve_banded(
    l_and_u: tuple[int, int],
    ab: _InputComplexStrict2D,
    b: onp.ToComplexStrict2D,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.Array2D[np.complex128]: ...
@overload  # Nd ~complex128, +complex128
def solve_banded(
    l_and_u: tuple[int, int],
    ab: _InputComplex,
    b: onp.ToComplexND,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.complex128]: ...
@overload  # 1d +complex128, ~complex128
def solve_banded(
    l_and_u: tuple[int, int],
    ab: onp.ToComplexStrict2D,
    b: _InputComplexStrict1D,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.Array1D[np.complex128]: ...
@overload  # 2d +complex128, ~complex128
def solve_banded(
    l_and_u: tuple[int, int],
    ab: onp.ToComplexStrict2D,
    b: _InputComplexStrict2D,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.Array2D[np.complex128]: ...
@overload  # Nd +complex128, ~complex128
def solve_banded(
    l_and_u: tuple[int, int],
    ab: onp.ToComplexND,
    b: _InputComplex,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.complex128]: ...
@overload  # 1d +floating, +floating
def solve_banded(
    l_and_u: tuple[int, int],
    ab: onp.ToFloatStrict2D,
    b: onp.ToFloatStrict1D,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.Array1D[np.float32 | np.float64]: ...
@overload  # 2d +floating, +floating
def solve_banded(
    l_and_u: tuple[int, int],
    ab: onp.ToFloatStrict2D,
    b: onp.ToFloatStrict2D,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.Array2D[np.float32 | np.float64]: ...
@overload  # Nd +floating, +floating
def solve_banded(
    l_and_u: tuple[int, int],
    ab: onp.ToFloatND,
    b: onp.ToFloatND,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.float32 | np.float64]: ...
@overload  # 1d +complexfloating, ~complexfloating
def solve_banded(
    l_and_u: tuple[int, int],
    ab: onp.ToComplexStrict2D,
    b: onp.ToJustComplexStrict1D,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.Array1D[np.complex64 | np.complex128]: ...
@overload  # 2d +complexfloating, ~complexfloating
def solve_banded(
    l_and_u: tuple[int, int],
    ab: onp.ToComplexStrict2D,
    b: onp.ToJustComplexStrict2D,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.Array2D[np.complex64 | np.complex128]: ...
@overload  # Nd +complexfloating, ~complexfloating
def solve_banded(
    l_and_u: tuple[int, int],
    ab: onp.ToComplexND,
    b: onp.ToJustComplexND,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.complex64 | np.complex128]: ...
@overload  # Nd +complexfloating, +complexfloating
def solve_banded(
    l_and_u: tuple[int, int],
    ab: onp.ToComplexND,
    b: onp.ToComplexND,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.float32 | np.float64 | np.complex64 | np.complex128]: ...

# NOTE: keep overload structure consistent with `solve_banded` above
@overload  # 1D ~float64, +float64
def solveh_banded(
    ab: _InputFloatStrict2D,
    b: onp.ToFloatStrict1D,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    lower: bool = False,
    check_finite: bool = True,
) -> onp.Array1D[np.float64]: ...
@overload  # 2D ~float64, +float64
def solveh_banded(
    ab: _InputFloatStrict2D,
    b: onp.ToFloatStrict2D,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    lower: bool = False,
    check_finite: bool = True,
) -> onp.Array2D[np.float64]: ...
@overload  # Nd ~float64, +float64
def solveh_banded(
    ab: _InputFloat,
    b: onp.ToFloatND,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    lower: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.float64]: ...
@overload  # 1d +float64, ~float64
def solveh_banded(
    ab: onp.ToFloatStrict2D,
    b: _InputFloatStrict1D,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    lower: bool = False,
    check_finite: bool = True,
) -> onp.Array1D[np.float64]: ...
@overload  # 2d +float64, ~float64
def solveh_banded(
    ab: onp.ToFloatStrict2D,
    b: _InputFloatStrict2D,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    lower: bool = False,
    check_finite: bool = True,
) -> onp.Array2D[np.float64]: ...
@overload  # Nd +float64, ~float64
def solveh_banded(
    ab: onp.ToFloatND,
    b: _InputFloat,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    lower: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.float64]: ...
@overload  # 1d ~complex128, +complex128
def solveh_banded(
    ab: _InputComplexStrict2D,
    b: onp.ToComplexStrict1D,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    lower: bool = False,
    check_finite: bool = True,
) -> onp.Array1D[np.complex128]: ...
@overload  # 2d ~complex128, +complex128
def solveh_banded(
    ab: _InputComplexStrict2D,
    b: onp.ToComplexStrict2D,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    lower: bool = False,
    check_finite: bool = True,
) -> onp.Array2D[np.complex128]: ...
@overload  # Nd ~complex128, +complex128
def solveh_banded(
    ab: _InputComplex,
    b: onp.ToComplexND,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    lower: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.complex128]: ...
@overload  # 1d +complex128, ~complex128
def solveh_banded(
    ab: onp.ToComplexStrict2D,
    b: _InputComplexStrict1D,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    lower: bool = False,
    check_finite: bool = True,
) -> onp.Array1D[np.complex128]: ...
@overload  # 2d +complex128, ~complex128
def solveh_banded(
    ab: onp.ToComplexStrict2D,
    b: _InputComplexStrict2D,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    lower: bool = False,
    check_finite: bool = True,
) -> onp.Array2D[np.complex128]: ...
@overload  # Nd +complex128, ~complex128
def solveh_banded(
    ab: onp.ToComplexND,
    b: _InputComplex,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    lower: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.complex128]: ...
@overload  # 1d +floating, +floating
def solveh_banded(
    ab: onp.ToFloatStrict2D,
    b: onp.ToFloatStrict1D,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    lower: bool = False,
    check_finite: bool = True,
) -> onp.Array1D[np.float32 | np.float64]: ...
@overload  # 2d +floating, +floating
def solveh_banded(
    ab: onp.ToFloatStrict2D,
    b: onp.ToFloatStrict2D,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    lower: bool = False,
    check_finite: bool = True,
) -> onp.Array2D[np.float32 | np.float64]: ...
@overload  # Nd +floating, +floating
def solveh_banded(
    ab: onp.ToFloatND,
    b: onp.ToFloatND,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    lower: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.float32 | np.float64]: ...
@overload  # 1d +complexfloating, ~complexfloating
def solveh_banded(
    ab: onp.ToComplexStrict2D,
    b: onp.ToJustComplexStrict1D,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    lower: bool = False,
    check_finite: bool = True,
) -> onp.Array1D[np.complex64 | np.complex128]: ...
@overload  # 2d +complexfloating, ~complexfloating
def solveh_banded(
    ab: onp.ToComplexStrict2D,
    b: onp.ToJustComplexStrict2D,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    lower: bool = False,
    check_finite: bool = True,
) -> onp.Array2D[np.complex64 | np.complex128]: ...
@overload  # Nd +complexfloating, ~complexfloating
def solveh_banded(
    ab: onp.ToComplexND,
    b: onp.ToJustComplexND,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    lower: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.complex64 | np.complex128]: ...
@overload  # Nd +complexfloating, +complexfloating
def solveh_banded(
    ab: onp.ToComplexND,
    b: onp.ToComplexND,
    overwrite_ab: bool = False,
    overwrite_b: bool = False,
    lower: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.float32 | np.float64 | np.complex64 | np.complex128]: ...

#
@overload  # 1d +float, +float
def solve_toeplitz(
    c_or_cr: _COrCR[onp.ToFloatStrict1D], b: onp.ToFloatStrict1D, check_finite: bool = True
) -> onp.Array1D[np.float64]: ...
@overload  # 2d +float, +float
def solve_toeplitz(
    c_or_cr: _COrCR[onp.ToFloatStrict1D], b: onp.ToFloatStrict2D, check_finite: bool = True
) -> onp.Array2D[np.float64]: ...
@overload  # Nd +float, +float
def solve_toeplitz(c_or_cr: _COrCR[onp.ToFloatND], b: onp.ToFloatND, check_finite: bool = True) -> onp.ArrayND[np.float64]: ...
@overload  # 1d ~complex, +complex
def solve_toeplitz(
    c_or_cr: _COrCR[onp.ToJustComplexStrict1D], b: onp.ToComplexStrict1D, check_finite: bool = True
) -> onp.Array1D[np.complex128]: ...
@overload  # 2d ~complex, +complex
def solve_toeplitz(
    c_or_cr: _COrCR[onp.ToJustComplexStrict1D], b: onp.ToComplexStrict2D, check_finite: bool = True
) -> onp.Array2D[np.complex128]: ...
@overload  # Nd ~complex, +complex
def solve_toeplitz(
    c_or_cr: _COrCR[onp.ToJustComplexND], b: onp.ToComplexND, check_finite: bool = True
) -> onp.ArrayND[np.complex128]: ...
@overload  # 1d +complex, ~complex
def solve_toeplitz(
    c_or_cr: _COrCR[onp.ToComplexStrict1D], b: onp.ToJustComplexStrict1D, check_finite: bool = True
) -> onp.Array1D[np.complex128]: ...
@overload  # 2d +complex, ~complex
def solve_toeplitz(
    c_or_cr: _COrCR[onp.ToComplexStrict1D], b: onp.ToJustComplexStrict2D, check_finite: bool = True
) -> onp.Array2D[np.complex128]: ...
@overload  # Nd +complex, ~complex
def solve_toeplitz(
    c_or_cr: _COrCR[onp.ToComplexND], b: onp.ToJustComplexND, check_finite: bool = True
) -> onp.ArrayND[np.complex128]: ...
@overload  # 1d +complex, +complex
def solve_toeplitz(
    c_or_cr: _COrCR[onp.ToComplexStrict1D], b: onp.ToComplexStrict1D, check_finite: bool = True
) -> onp.Array1D[np.float64 | np.complex128]: ...
@overload  # 2d +complex, +complex
def solve_toeplitz(
    c_or_cr: _COrCR[onp.ToComplexStrict1D], b: onp.ToComplexStrict2D, check_finite: bool = True
) -> onp.Array2D[np.float64 | np.complex128]: ...
@overload  # Nd +complex, +complex
def solve_toeplitz(
    c_or_cr: _COrCR[onp.ToComplexND], b: onp.ToComplexND, check_finite: bool = True
) -> onp.ArrayND[np.float64 | np.complex128]: ...

#
@overload  # 1D ~float64, +float64
def solve_circulant(
    c: _InputF64Strict1D,
    b: onp.ToFloat64Strict1D,
    singular: _Singular = "raise",
    tol: float | None = None,
    caxis: int = -1,
    baxis: int = 0,
    outaxis: int = 0,
) -> onp.Array1D[np.float64]: ...
@overload  # 2D ~float64, +float64
def solve_circulant(
    c: _InputF64Strict1D,
    b: onp.ToFloat64Strict2D,
    singular: _Singular = "raise",
    tol: float | None = None,
    caxis: int = -1,
    baxis: int = 0,
    outaxis: int = 0,
) -> onp.Array2D[np.float64]: ...
@overload  # Nd ~float64, +float64
def solve_circulant(
    c: _InputF64,
    b: onp.ToFloat64_ND,
    singular: _Singular = "raise",
    tol: float | None = None,
    caxis: int = -1,
    baxis: int = 0,
    outaxis: int = 0,
) -> onp.ArrayND[np.float64]: ...
@overload  # 1d +float64, ~float64
def solve_circulant(
    c: onp.ToFloat64Strict1D,
    b: _InputF64Strict1D,
    singular: _Singular = "raise",
    tol: float | None = None,
    caxis: int = -1,
    baxis: int = 0,
    outaxis: int = 0,
) -> onp.Array1D[np.float64]: ...
@overload  # 2d +float64, ~float64
def solve_circulant(
    c: onp.ToFloat64Strict1D,
    b: _InputF64Strict2D,
    singular: _Singular = "raise",
    tol: float | None = None,
    caxis: int = -1,
    baxis: int = 0,
    outaxis: int = 0,
) -> onp.Array2D[np.float64]: ...
@overload  # Nd +float64, ~float64
def solve_circulant(
    c: onp.ToFloatND,
    b: _InputF64,
    singular: _Singular = "raise",
    tol: float | None = None,
    caxis: int = -1,
    baxis: int = 0,
    outaxis: int = 0,
) -> onp.ArrayND[np.float64]: ...
@overload  # 1d ~complex128, +complex128
def solve_circulant(
    c: onp.ToJustComplex128Strict1D,
    b: onp.ToComplex128Strict1D,
    singular: _Singular = "raise",
    tol: float | None = None,
    caxis: int = -1,
    baxis: int = 0,
    outaxis: int = 0,
) -> onp.Array1D[np.complex128]: ...
@overload  # 2d ~complex128, +complex128
def solve_circulant(
    c: onp.ToJustComplex128Strict1D,
    b: onp.ToComplex128Strict2D,
    singular: _Singular = "raise",
    tol: float | None = None,
    caxis: int = -1,
    baxis: int = 0,
    outaxis: int = 0,
) -> onp.Array2D[np.complex128]: ...
@overload  # Nd ~complex128, +complex128
def solve_circulant(
    c: onp.ToJustComplex128_ND,
    b: onp.ToComplex128_ND,
    singular: _Singular = "raise",
    tol: float | None = None,
    caxis: int = -1,
    baxis: int = 0,
    outaxis: int = 0,
) -> onp.ArrayND[np.complex128]: ...
@overload  # 1d +complex128, ~complex128
def solve_circulant(
    c: onp.ToComplex128Strict1D,
    b: onp.ToJustComplex128Strict1D,
    singular: _Singular = "raise",
    tol: float | None = None,
    caxis: int = -1,
    baxis: int = 0,
    outaxis: int = 0,
) -> onp.Array1D[np.complex128]: ...
@overload  # 2d +complex128, ~complex128
def solve_circulant(
    c: onp.ToComplex128Strict1D,
    b: onp.ToJustComplex128Strict2D,
    singular: _Singular = "raise",
    tol: float | None = None,
    caxis: int = -1,
    baxis: int = 0,
    outaxis: int = 0,
) -> onp.Array2D[np.complex128]: ...
@overload  # Nd +complex128, ~complex128
def solve_circulant(
    c: onp.ToComplex128_ND,
    b: onp.ToJustComplex128_ND,
    singular: _Singular = "raise",
    tol: float | None = None,
    caxis: int = -1,
    baxis: int = 0,
    outaxis: int = 0,
) -> onp.ArrayND[np.complex128]: ...
@overload  # 1d +floating, +floating
def solve_circulant(
    c: onp.ToFloatStrict1D,
    b: onp.ToFloatStrict1D,
    singular: _Singular = "raise",
    tol: onp.ToFloat | None = None,
    caxis: SupportsIndex = -1,
    baxis: SupportsIndex = 0,
    outaxis: SupportsIndex = 0,
) -> onp.Array1D[npc.floating]: ...
@overload  # 2d +floating, +floating
def solve_circulant(
    c: onp.ToFloatStrict1D,
    b: onp.ToFloatStrict2D,
    singular: _Singular = "raise",
    tol: onp.ToFloat | None = None,
    caxis: SupportsIndex = -1,
    baxis: SupportsIndex = 0,
    outaxis: SupportsIndex = 0,
) -> onp.Array2D[npc.floating]: ...
@overload  # Nd +floating, +floating
def solve_circulant(
    c: onp.ToFloatND,
    b: onp.ToFloatND,
    singular: _Singular = "raise",
    tol: float | None = None,
    caxis: int = -1,
    baxis: int = 0,
    outaxis: int = 0,
) -> onp.ArrayND[npc.floating]: ...
@overload  # 1d +complexfloating, ~complexfloating
def solve_circulant(
    c: onp.ToComplexStrict1D,
    b: onp.ToJustComplexStrict1D,
    singular: _Singular = "raise",
    tol: onp.ToFloat | None = None,
    caxis: SupportsIndex = -1,
    baxis: SupportsIndex = 0,
    outaxis: SupportsIndex = 0,
) -> onp.Array1D[npc.complexfloating]: ...
@overload  # 2d +complexfloating, ~complexfloating
def solve_circulant(
    c: onp.ToComplexStrict1D,
    b: onp.ToJustComplexStrict2D,
    singular: _Singular = "raise",
    tol: onp.ToFloat | None = None,
    caxis: SupportsIndex = -1,
    baxis: SupportsIndex = 0,
    outaxis: SupportsIndex = 0,
) -> onp.Array2D[npc.complexfloating]: ...
@overload  # Nd +complexfloating, ~complexfloating
def solve_circulant(
    c: onp.ToComplexND,
    b: onp.ToJustComplexND,
    singular: _Singular = "raise",
    tol: float | None = None,
    caxis: int = -1,
    baxis: int = 0,
    outaxis: int = 0,
) -> onp.ArrayND[npc.complexfloating]: ...
@overload  # Nd +complexfloating, +complexfloating
def solve_circulant(
    c: onp.ToComplexND,
    b: onp.ToComplexND,
    singular: _Singular = "raise",
    tol: float | None = None,
    caxis: int = -1,
    baxis: int = 0,
    outaxis: int = 0,
) -> onp.ArrayND[npc.inexact]: ...

#
@overload  # 2d bool sequence
def inv(
    a: Sequence[Sequence[bool]],
    overwrite_a: bool = False,
    check_finite: bool = True,
    *,
    assume_a: _AssumeA | None = None,
    lower: bool = False,
) -> onp.Array2D[np.float32]: ...
@overload  # Nd bool sequence
def inv(
    a: Sequence[onp.SequenceND[bool]],
    overwrite_a: bool = False,
    check_finite: bool = True,
    *,
    assume_a: _AssumeA | None = None,
    lower: bool = False,
) -> onp.ArrayND[np.float32]: ...
@overload  # 2d float or int sequence
def inv(
    a: Sequence[Sequence[op.JustFloat | op.JustInt]],
    overwrite_a: bool = False,
    check_finite: bool = True,
    *,
    assume_a: _AssumeA | None = None,
    lower: bool = False,
) -> onp.Array2D[np.float64]: ...
@overload  # Nd float or int sequence
def inv(
    a: Sequence[onp.SequenceND[op.JustFloat | op.JustInt]],
    overwrite_a: bool = False,
    check_finite: bool = True,
    *,
    assume_a: _AssumeA | None = None,
    lower: bool = False,
) -> onp.ArrayND[np.float64]: ...
@overload  # 2d complex sequence
def inv(
    a: Sequence[Sequence[op.JustComplex]],
    overwrite_a: bool = False,
    check_finite: bool = True,
    *,
    assume_a: _AssumeA | None = None,
    lower: bool = False,
) -> onp.Array2D[np.complex128]: ...
@overload  # Nd complex sequence
def inv(
    a: Sequence[onp.SequenceND[op.JustComplex]],
    overwrite_a: bool = False,
    check_finite: bool = True,
    *,
    assume_a: _AssumeA | None = None,
    lower: bool = False,
) -> onp.ArrayND[np.complex128]: ...
@overload  # generic shape, as float32
def inv[ShapeT: tuple[int, ...]](
    a: onp.CanArrayND[np.float32 | npc.number16 | npc.integer8 | np.bool, ShapeT],
    overwrite_a: bool = False,
    check_finite: bool = True,
    *,
    assume_a: _AssumeA | None = None,
    lower: bool = False,
) -> onp.ArrayND[np.float32, ShapeT]: ...
@overload  # generic shape, as float64
def inv[ShapeT: tuple[int, ...]](
    a: onp.CanArrayND[np.float64 | npc.floating80 | npc.integer64 | npc.integer32, ShapeT],
    overwrite_a: bool = False,
    check_finite: bool = True,
    *,
    assume_a: _AssumeA | None = None,
    lower: bool = False,
) -> onp.ArrayND[np.float64, ShapeT]: ...
@overload  # generic shape, as complex64
def inv[ShapeT: tuple[int, ...]](
    a: onp.CanArrayND[np.complex64, ShapeT],
    overwrite_a: bool = False,
    check_finite: bool = True,
    *,
    assume_a: _AssumeA | None = None,
    lower: bool = False,
) -> onp.ArrayND[np.complex64, ShapeT]: ...
@overload  # generic shape, as complex128
def inv[ShapeT: tuple[int, ...]](
    a: onp.CanArrayND[np.complex128 | npc.complexfloating160, ShapeT],
    overwrite_a: bool = False,
    check_finite: bool = True,
    *,
    assume_a: _AssumeA | None = None,
    lower: bool = False,
) -> onp.ArrayND[np.complex128, ShapeT]: ...

# NOTE: The order of the overloads has been carefully chosen to avoid triggering a Pyright bug.
@overload  # +float64 2d
def det(a: onp.ToFloat64Strict2D, overwrite_a: bool = False, check_finite: bool = True) -> np.float64: ...
@overload  # complex128 | complex64 2d
def det(
    a: onp.ToArrayStrict2D[op.JustComplex, np.complex128 | np.complex64], overwrite_a: bool = False, check_finite: bool = True
) -> np.complex128: ...
@overload  # +float64 3d
def det(a: onp.ToFloat64Strict3D, overwrite_a: bool = False, check_finite: bool = True) -> onp.Array1D[np.float64]: ...
@overload  # complex128 | complex64 3d
def det(
    a: onp.ToArrayStrict3D[op.JustComplex, np.complex128 | np.complex64], overwrite_a: bool = False, check_finite: bool = True
) -> onp.Array1D[np.complex128]: ...
@overload  # +float64 ND
def det(a: onp.ToFloat64_ND, overwrite_a: bool = False, check_finite: bool = True) -> np.float64 | onp.ArrayND[np.float64]: ...
@overload  # complex128 | complex64 Nd
def det(
    a: onp.ToArrayND[op.JustComplex, np.complex128 | np.complex64], overwrite_a: bool = False, check_finite: bool = True
) -> np.complex128 | onp.ArrayND[np.complex128]: ...
@overload  # +complex128 2d
def det(a: onp.ToComplex128Strict2D, overwrite_a: bool = False, check_finite: bool = True) -> np.float64 | np.complex128: ...
@overload  # +complex128 3d
def det(
    a: onp.ToComplex128Strict3D, overwrite_a: bool = False, check_finite: bool = True
) -> onp.Array1D[np.float64 | np.complex128]: ...
@overload  # +complex128 Nd
def det(
    a: onp.ToComplex128_ND, overwrite_a: bool = False, check_finite: bool = True
) -> np.float64 | np.complex128 | onp.ArrayND[np.float64 | np.complex128]: ...

#
@overload  # ~f64, +f64 1d
def lstsq(
    a: _AsFloat64Strict2D,
    b: onp.ToFloatStrict1D,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriverDS | None = None,
) -> _LstSqResult1N[np.float64, onp.Array1D[np.float64]]: ...
@overload  # ~f64, +f64 1d, lapack_driver='gelsy' (keyword)
def lstsq(
    a: _AsFloat64Strict2D,
    b: onp.ToFloatStrict1D,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    *,
    lapack_driver: _LapackDriverY,
) -> _LstSqResult11[np.float64]: ...
@overload  # ~f64, +f64 2d
def lstsq(
    a: _AsFloat64Strict2D,
    b: onp.ToFloatStrict2D,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriverDS | None = None,
) -> _LstSqResult21[np.float64, onp.Array1D[np.float64]]: ...
@overload  # ~f64, +f64 2d, lapack_driver='gelsy' (keyword)
def lstsq(
    a: _AsFloat64Strict2D,
    b: onp.ToFloatStrict2D,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    *,
    lapack_driver: _LapackDriverY,
) -> _LstSqResult21[np.float64, None]: ...
@overload  # ~f64, +f64 Nd
def lstsq(
    a: _AsFloat64_2D,
    b: onp.ToFloatND,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriverDS | None = None,
) -> _LstSqResultND[np.float64, onp.ArrayND[np.float64]]: ...
@overload  # ~f64, +f64 Nd, lapack_driver='gelsy' (keyword)
def lstsq(
    a: _AsFloat64_2D,
    b: onp.ToFloatND,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    *,
    lapack_driver: _LapackDriverY,
) -> _LstSqResultND[np.float64, None]: ...
@overload  # +f64, ~f64 1d
def lstsq(
    a: onp.ToFloatStrict2D,
    b: _AsFloat64Strict1D,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriverDS | None = None,
) -> _LstSqResult1N[np.float64, onp.Array1D[np.float64]]: ...
@overload  # +f64, ~f64 1d, lapack_driver='gelsy' (keyword)
def lstsq(
    a: onp.ToFloatStrict2D,
    b: _AsFloat64Strict1D,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    *,
    lapack_driver: _LapackDriverY,
) -> _LstSqResult11[np.float64]: ...
@overload  # +f64, ~f64 2d
def lstsq(
    a: onp.ToFloatStrict2D,
    b: _AsFloat64Strict2D,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriverDS | None = None,
) -> _LstSqResult21[np.float64, onp.Array1D[np.float64]]: ...
@overload  # +f64, ~f64 2d, lapack_driver='gelsy' (keyword)
def lstsq(
    a: onp.ToFloatStrict2D,
    b: _AsFloat64Strict2D,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    *,
    lapack_driver: _LapackDriverY,
) -> _LstSqResult21[np.float64, None]: ...
@overload  # +f64, ~f64 Nd
def lstsq(
    a: onp.ToFloat2D,
    b: _AsFloat64_ND,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriverDS | None = None,
) -> _LstSqResultND[np.float64, onp.ArrayND[np.float64]]: ...
@overload  # +f64, ~f64 Nd, lapack_driver='gelsy' (keyword)
def lstsq(
    a: onp.ToFloat2D,
    b: _AsFloat64_ND,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    *,
    lapack_driver: _LapackDriverY,
) -> _LstSqResultND[np.float64, None]: ...
@overload  # ~f32 1d
def lstsq(
    a: _AsFloat32Strict2D,
    b: _AsFloat32Strict1D,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriverDS | None = None,
) -> _LstSqResult1N[np.float32, onp.Array1D[np.float32]]: ...
@overload  # ~f32 1d, lapack_driver='gelsy' (keyword)
def lstsq(
    a: _AsFloat32Strict2D,
    b: _AsFloat32Strict1D,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    *,
    lapack_driver: _LapackDriverY,
) -> _LstSqResult11[np.float32]: ...
@overload  # ~f32 2d
def lstsq(
    a: _AsFloat32Strict2D,
    b: _AsFloat32Strict2D,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriverDS | None = None,
) -> _LstSqResult21[np.float32, onp.Array1D[np.float32]]: ...
@overload  # ~f32 2d, lapack_driver='gelsy' (keyword)
def lstsq(
    a: _AsFloat32Strict2D,
    b: _AsFloat32Strict2D,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    *,
    lapack_driver: _LapackDriverY,
) -> _LstSqResult21[np.float32, None]: ...
@overload  # ~f32 Nd
def lstsq(
    a: _AsFloat32_2D,
    b: _AsFloat32_ND,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriverDS | None = None,
) -> _LstSqResultND[np.float32, onp.ArrayND[np.float32]]: ...
@overload  # ~f32 Nd, lapack_driver='gelsy' (keyword)
def lstsq(
    a: _AsFloat32_2D,
    b: _AsFloat32_ND,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    *,
    lapack_driver: _LapackDriverY,
) -> _LstSqResultND[np.float32, None]: ...
@overload  # ~c128, +c128 1d
def lstsq(
    a: _AsComplex128Strict2D,
    b: onp.ToComplexStrict1D,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriverDS | None = None,
) -> _LstSqResult1N[np.complex128, onp.Array1D[np.float64]]: ...
@overload  # ~c128, +c128 1d, lapack_driver='gelsy' (keyword)
def lstsq(
    a: _AsComplex128Strict2D,
    b: onp.ToComplexStrict1D,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    *,
    lapack_driver: _LapackDriverY,
) -> _LstSqResult11[np.complex128]: ...
@overload  # ~c128, +c128 2d
def lstsq(
    a: _AsComplex128Strict2D,
    b: onp.ToComplexStrict2D,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriverDS | None = None,
) -> _LstSqResult21[np.complex128, onp.Array1D[np.float64]]: ...
@overload  # ~c128, +c128 2d, lapack_driver='gelsy' (keyword)
def lstsq(
    a: _AsComplex128Strict2D,
    b: onp.ToComplexStrict2D,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    *,
    lapack_driver: _LapackDriverY,
) -> _LstSqResult21[np.complex128, None]: ...
@overload  # ~c128, +c128 Nd
def lstsq(
    a: _AsComplex128_2D,
    b: onp.ToComplexND,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriverDS | None = None,
) -> _LstSqResultND[np.complex128, onp.ArrayND[np.float64]]: ...
@overload  # ~c128, +c128 Nd, lapack_driver='gelsy' (keyword)
def lstsq(
    a: _AsComplex128_2D,
    b: onp.ToComplexND,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    *,
    lapack_driver: _LapackDriverY,
) -> _LstSqResultND[np.complex128, None]: ...
@overload  # +c128, ~c128 1d
def lstsq(
    a: onp.ToComplexStrict2D,
    b: _AsComplex128Strict1D,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriverDS | None = None,
) -> _LstSqResult1N[np.complex128, onp.Array1D[np.float64]]: ...
@overload  # +c128, ~c128 1d, lapack_driver='gelsy' (keyword)
def lstsq(
    a: onp.ToComplexStrict2D,
    b: _AsComplex128Strict1D,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    *,
    lapack_driver: _LapackDriverY,
) -> _LstSqResult11[np.complex128]: ...
@overload  # +c128, ~c128 2d
def lstsq(
    a: onp.ToComplexStrict2D,
    b: _AsComplex128Strict2D,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriverDS | None = None,
) -> _LstSqResult21[np.complex128, onp.Array1D[np.float64]]: ...
@overload  # +c128, ~c128 2d, lapack_driver='gelsy' (keyword)
def lstsq(
    a: onp.ToComplexStrict2D,
    b: _AsComplex128Strict2D,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    *,
    lapack_driver: _LapackDriverY,
) -> _LstSqResult21[np.complex128, None]: ...
@overload  # +c128, ~c128 Nd
def lstsq(
    a: onp.ToComplex2D,
    b: _AsComplex128_ND,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriverDS | None = None,
) -> _LstSqResultND[np.complex128, onp.ArrayND[np.float64]]: ...
@overload  # +c128, ~c128 Nd, lapack_driver='gelsy' (keyword)
def lstsq(
    a: onp.ToComplex2D,
    b: _AsComplex128_ND,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    *,
    lapack_driver: _LapackDriverY,
) -> _LstSqResultND[np.complex128, None]: ...
@overload  # ~c64 1d
def lstsq(
    a: onp.ToJustComplex64Strict2D,
    b: _AsComplex64Strict1D,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriverDS | None = None,
) -> _LstSqResult1N[np.complex64, onp.Array1D[np.float32]]: ...
@overload  # ~c64 1d, lapack_driver='gelsy' (keyword)
def lstsq(
    a: onp.ToJustComplex64Strict2D,
    b: _AsComplex64Strict1D,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    *,
    lapack_driver: _LapackDriverY,
) -> _LstSqResult11[np.complex64]: ...
@overload  # ~c64 2d
def lstsq(
    a: onp.ToJustComplex64Strict2D,
    b: _AsComplex64Strict2D,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriverDS | None = None,
) -> _LstSqResult21[np.complex64, onp.Array1D[np.float32]]: ...
@overload  # ~c64 2d, lapack_driver='gelsy' (keyword)
def lstsq(
    a: onp.ToJustComplex64Strict2D,
    b: _AsComplex64Strict2D,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    *,
    lapack_driver: _LapackDriverY,
) -> _LstSqResult21[np.complex64, None]: ...
@overload  # ~c64 Nd
def lstsq(
    a: onp.ToJustComplex64_2D,
    b: _AsComplex64_ND,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriverDS | None = None,
) -> _LstSqResultND[np.complex64, onp.ArrayND[np.float32]]: ...
@overload  # ~c64 Nd, lapack_driver='gelsy' (keyword)
def lstsq(
    a: onp.ToJustComplex64_2D,
    b: _AsComplex64_ND,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    *,
    lapack_driver: _LapackDriverY,
) -> _LstSqResultND[np.complex64, None]: ...
@overload  # +fallback
def lstsq(
    a: onp.ToComplex2D,
    b: onp.ToComplexND,
    cond: float | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriver | None = None,
) -> _LstSqResultND[Incomplete, onp.ArrayND[np.float64 | Any] | Any]: ...

# TODO(jorenham): improve this
@overload
def pinv(  # (float[:, :], return_rank=False) -> float[:, :]
    a: onp.ToFloatND,
    *,
    atol: onp.ToFloat | None = None,
    rtol: onp.ToFloat | None = None,
    return_rank: Literal[False] = False,
    check_finite: bool = True,
) -> _FloatND: ...
@overload  # (float[:, :], return_rank=True) -> (float[:, :], int)
def pinv(
    a: onp.ToFloatND,
    *,
    atol: onp.ToFloat | None = None,
    rtol: onp.ToFloat | None = None,
    return_rank: Literal[True],
    check_finite: bool = True,
) -> tuple[_FloatND, np.int_]: ...
@overload  # (complex[:, :], return_rank=False) -> complex[:, :]
def pinv(
    a: onp.ToComplexND,
    *,
    atol: onp.ToFloat | None = None,
    rtol: onp.ToFloat | None = None,
    return_rank: Literal[False] = False,
    check_finite: bool = True,
) -> _InexactND: ...
@overload  # (complex[:, :], return_rank=True) -> (complex[:, :], int)
def pinv(
    a: onp.ToComplexND,
    *,
    atol: onp.ToFloat | None = None,
    rtol: onp.ToFloat | None = None,
    return_rank: Literal[True],
    check_finite: bool = True,
) -> tuple[_InexactND, np.int_]: ...

# TODO(jorenham): improve this
@overload  # (float[:, :], return_rank=False) -> float[:, :]
def pinvh(
    a: onp.ToFloatND,
    atol: onp.ToFloat | None = None,
    rtol: onp.ToFloat | None = None,
    lower: bool = True,
    return_rank: Literal[False] = False,
    check_finite: bool = True,
) -> _FloatND: ...
@overload  # (float[:, :], return_rank=True, /) -> (float[:, :], int)
def pinvh(
    a: onp.ToFloatND,
    atol: onp.ToFloat | None,
    rtol: onp.ToFloat | None,
    lower: bool,
    return_rank: Literal[True],
    check_finite: bool = True,
) -> tuple[_FloatND, int]: ...
@overload  # (float[:, :], *, return_rank=True) -> (float[:, :], int)
def pinvh(
    a: onp.ToFloatND,
    atol: onp.ToFloat | None = None,
    rtol: onp.ToFloat | None = None,
    lower: bool = True,
    *,
    return_rank: Literal[True],
    check_finite: bool = True,
) -> tuple[_FloatND, int]: ...
@overload  # (complex[:, :], return_rank=False) -> complex[:, :]
def pinvh(
    a: onp.ToComplexND,
    atol: onp.ToFloat | None = None,
    rtol: onp.ToFloat | None = None,
    lower: bool = True,
    return_rank: Literal[False] = False,
    check_finite: bool = True,
) -> _InexactND: ...
@overload  # (complex[:, :], return_rank=True, /) -> (complex[:, :], int)
def pinvh(
    a: onp.ToComplexND,
    atol: onp.ToFloat | None,
    rtol: onp.ToFloat | None,
    lower: bool,
    return_rank: Literal[True],
    check_finite: bool = True,
) -> tuple[_InexactND, int]: ...
@overload  # (complex[:, :], *, return_rank=True) -> (complex[:, :], int)
def pinvh(
    a: onp.ToComplexND,
    atol: onp.ToFloat | None = None,
    rtol: onp.ToFloat | None = None,
    lower: bool = True,
    *,
    return_rank: Literal[True],
    check_finite: bool = True,
) -> tuple[_InexactND, int]: ...

# TODO(jorenham): improve this
@overload  # (float[:, :], separate=True) -> (float[:, :], float[:, :])
def matrix_balance(
    A: onp.ToFloatND, permute: bool = True, scale: bool = True, separate: Literal[False] = False, overwrite_a: bool = False
) -> _Tuple2[_FloatND]: ...
@overload  # (float[:, :], separate=False, /) -> (float[:, :], (float[:], float[:]))
def matrix_balance(
    A: onp.ToFloatND, permute: bool, scale: bool, separate: Literal[True], overwrite_a: bool = False
) -> tuple[_FloatND, _Tuple2[_FloatND]]: ...
@overload  # (float[:, :], *, separate=False) -> (float[:, :], (float[:], float[:]))
def matrix_balance(
    A: onp.ToFloatND, permute: bool = True, scale: bool = True, *, separate: Literal[True], overwrite_a: bool = False
) -> tuple[_FloatND, _Tuple2[_FloatND]]: ...
@overload  # (complex[:, :], separate=True) -> (complex[:, :], complex[:, :])
def matrix_balance(
    A: onp.ToComplexND, permute: bool = True, scale: bool = True, separate: Literal[False] = False, overwrite_a: bool = False
) -> _Tuple2[_InexactND]: ...
@overload  # (complex[:, :], separate=False, /) -> (complex[:, :], (complex[:], complex[:]))
def matrix_balance(
    A: onp.ToComplexND, permute: bool, scale: bool, separate: Literal[True], overwrite_a: bool = False
) -> tuple[_InexactND, _Tuple2[_InexactND]]: ...
@overload  # (complex[:, :], *, separate=False) -> (complex[:, :], (complex[:], complex[:]))
def matrix_balance(
    A: onp.ToComplexND, permute: bool = True, scale: bool = True, *, separate: Literal[True], overwrite_a: bool = False
) -> tuple[_InexactND, _Tuple2[_InexactND]]: ...

# TODO(jorenham): improve this
@overload  # floating 1d, 1d
def matmul_toeplitz(
    c_or_cr: onp.ToFloatStrict1D | _Tuple2[onp.ToFloatStrict1D],
    x: onp.ToFloatStrict1D,
    check_finite: bool = False,
    workers: int | None = None,
) -> _Float1D: ...
@overload  # floating
def matmul_toeplitz(
    c_or_cr: onp.ToFloatND | _Tuple2[onp.ToFloatND], x: onp.ToFloatND, check_finite: bool = False, workers: int | None = None
) -> _FloatND: ...
@overload  # complexfloating 1d, 1d
def matmul_toeplitz(
    c_or_cr: onp.ToComplexStrict1D | _Tuple2[onp.ToComplexStrict1D],
    x: onp.ToComplexStrict1D,
    check_finite: bool = False,
    workers: int | None = None,
) -> _Inexact1D: ...
@overload  # complexfloating
def matmul_toeplitz(
    c_or_cr: onp.ToComplexND | _Tuple2[onp.ToComplexND],
    x: onp.ToComplexND,
    check_finite: bool = False,
    workers: int | None = None,
) -> _InexactND: ...
