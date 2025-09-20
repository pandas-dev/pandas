# mypy: disable-error-code=overload-overlap

from collections.abc import Sequence
from typing import Final, Literal, TypeAlias, TypeVar, overload

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

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

_ShapeT = TypeVar("_ShapeT", bound=tuple[int, ...])
_T = TypeVar("_T")

_Tuple2: TypeAlias = tuple[_T, _T]
_COrCR: TypeAlias = _T | _Tuple2[_T]

_Float: TypeAlias = npc.floating
_Float0D: TypeAlias = onp.Array0D[_Float]
_Float1D: TypeAlias = onp.Array1D[_Float]
_FloatND: TypeAlias = onp.ArrayND[_Float]

_Inexact: TypeAlias = npc.inexact
_Inexact0D: TypeAlias = onp.Array0D[_Inexact]
_Inexact1D: TypeAlias = onp.Array1D[_Inexact]
_InexactND: TypeAlias = onp.ArrayND[_Inexact]

_InputFloat: TypeAlias = onp.ToArrayND[float, np.float64 | npc.floating80 | npc.integer | np.bool_]
_InputFloatStrict1D: TypeAlias = onp.ToArrayStrict1D[float, np.float64 | npc.floating80 | npc.integer | np.bool_]
_InputFloatStrict2D: TypeAlias = onp.ToArrayStrict2D[float, np.float64 | npc.floating80 | npc.integer | np.bool_]

_InputF64: TypeAlias = onp.ToArrayND[float, np.float64 | npc.integer | np.bool_]
_InputF64Strict1D: TypeAlias = onp.ToArrayStrict1D[float, np.float64 | npc.integer | np.bool_]
_InputF64Strict2D: TypeAlias = onp.ToArrayStrict2D[float, np.float64 | npc.integer | np.bool_]

_InputComplex: TypeAlias = onp.ToArrayND[op.JustComplex, np.complex128 | npc.complexfloating160]
_InputComplexStrict1D: TypeAlias = onp.ToArrayStrict1D[op.JustComplex, np.complex128 | npc.complexfloating160]
_InputComplexStrict2D: TypeAlias = onp.ToArrayStrict2D[op.JustComplex, np.complex128 | npc.complexfloating160]

_AssumeA: TypeAlias = Literal[
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
_TransSystem: TypeAlias = Literal[0, "N", 1, "T", 2, "C"]
_Singular: TypeAlias = Literal["lstsq", "raise"]
_LapackDriver: TypeAlias = Literal["gelsd", "gelsy", "gelss"]

###

lapack_cast_dict: Final[dict[str, str]] = ...

@overload  # 2D ~float64, +float64
def solve(
    a: _InputFloatStrict2D,
    b: onp.ToFloatStrict1D | onp.ToFloatStrict2D,
    lower: bool = False,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    assume_a: _AssumeA | None = None,
    transposed: bool = False,
) -> onp.Array2D[np.float64]: ...
@overload  # Nd ~float64, +float64
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
@overload  # 2d +float64, ~float64
def solve(
    a: onp.ToFloatStrict2D,
    b: _InputFloatStrict1D | _InputFloatStrict2D,
    lower: bool = False,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    assume_a: _AssumeA | None = None,
    transposed: bool = False,
) -> onp.Array2D[np.float64]: ...
@overload  # Nd +float64, ~float64
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
@overload  # 2d ~complex128, +complex128
def solve(
    a: _InputComplexStrict2D,
    b: onp.ToComplexStrict1D | onp.ToComplexStrict2D,
    lower: bool = False,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    assume_a: _AssumeA | None = None,
    transposed: bool = False,
) -> onp.Array2D[np.complex128]: ...
@overload  # Nd ~complex128, +complex128
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
@overload  # 2d +complex128, ~complex128
def solve(
    a: onp.ToComplexStrict2D,
    b: _InputComplexStrict1D | _InputComplexStrict2D,
    lower: bool = False,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    assume_a: _AssumeA | None = None,
    transposed: bool = False,
) -> onp.Array2D[np.complex128]: ...
@overload  # Nd +complex128, ~complex128
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
@overload  # 2d +floating, +floating
def solve(
    a: onp.ToFloatStrict2D,
    b: onp.ToFloatStrict1D | onp.ToFloatStrict2D,
    lower: bool = False,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    assume_a: _AssumeA | None = None,
    transposed: bool = False,
) -> onp.Array2D[np.float32 | np.float64]: ...
@overload  # Nd +floating, +floating
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
@overload  # 2d +complexfloating, ~complexfloating
def solve(
    a: onp.ToComplexStrict2D,
    b: onp.ToJustComplexStrict1D | onp.ToJustComplexStrict2D,
    lower: bool = False,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    assume_a: _AssumeA | None = None,
    transposed: bool = False,
) -> onp.Array2D[np.complex64 | np.complex128]: ...
@overload  # Nd +complexfloating, ~complexfloating
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
@overload  # Nd +complexfloating, +complexfloating
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
    caxis: op.CanIndex = -1,
    baxis: op.CanIndex = 0,
    outaxis: op.CanIndex = 0,
) -> onp.Array1D[npc.floating]: ...
@overload  # 2d +floating, +floating
def solve_circulant(
    c: onp.ToFloatStrict1D,
    b: onp.ToFloatStrict2D,
    singular: _Singular = "raise",
    tol: onp.ToFloat | None = None,
    caxis: op.CanIndex = -1,
    baxis: op.CanIndex = 0,
    outaxis: op.CanIndex = 0,
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
    caxis: op.CanIndex = -1,
    baxis: op.CanIndex = 0,
    outaxis: op.CanIndex = 0,
) -> onp.Array1D[npc.complexfloating]: ...
@overload  # 2d +complexfloating, ~complexfloating
def solve_circulant(
    c: onp.ToComplexStrict1D,
    b: onp.ToJustComplexStrict2D,
    singular: _Singular = "raise",
    tol: onp.ToFloat | None = None,
    caxis: op.CanIndex = -1,
    baxis: op.CanIndex = 0,
    outaxis: op.CanIndex = 0,
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
def inv(a: Sequence[Sequence[bool]], overwrite_a: bool = False, check_finite: bool = True) -> onp.Array2D[np.float32]: ...
@overload  # Nd bool sequence
def inv(a: Sequence[onp.SequenceND[bool]], overwrite_a: bool = False, check_finite: bool = True) -> onp.ArrayND[np.float32]: ...
@overload  # 2d float or int sequence
def inv(
    a: Sequence[Sequence[op.JustFloat | op.JustInt]], overwrite_a: bool = False, check_finite: bool = True
) -> onp.Array2D[np.float64]: ...
@overload  # Nd float or int sequence
def inv(
    a: Sequence[onp.SequenceND[op.JustFloat | op.JustInt]], overwrite_a: bool = False, check_finite: bool = True
) -> onp.ArrayND[np.float64]: ...
@overload  # 2d complex sequence
def inv(
    a: Sequence[Sequence[op.JustComplex]], overwrite_a: bool = False, check_finite: bool = True
) -> onp.Array2D[np.complex128]: ...
@overload  # Nd complex sequence
def inv(
    a: Sequence[onp.SequenceND[op.JustComplex]], overwrite_a: bool = False, check_finite: bool = True
) -> onp.ArrayND[np.complex128]: ...
@overload  # generic shape, as float32
def inv(
    a: onp.CanArrayND[np.float32 | npc.number16 | npc.integer8 | np.bool_, _ShapeT],
    overwrite_a: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.float32, _ShapeT]: ...
@overload  # generic shape, as float64
def inv(
    a: onp.CanArrayND[np.float64 | npc.floating80 | npc.integer64 | npc.integer32, _ShapeT],
    overwrite_a: bool = False,
    check_finite: bool = True,
) -> onp.ArrayND[np.float64, _ShapeT]: ...
@overload  # generic shape, as complex64
def inv(
    a: onp.CanArrayND[np.complex64, _ShapeT], overwrite_a: bool = False, check_finite: bool = True
) -> onp.ArrayND[np.complex64, _ShapeT]: ...
@overload  # generic shape, as complex128
def inv(
    a: onp.CanArrayND[np.complex128 | npc.complexfloating160, _ShapeT], overwrite_a: bool = False, check_finite: bool = True
) -> onp.ArrayND[np.complex128, _ShapeT]: ...

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

# TODO(jorenham): improve this
@overload  # (float[:, :], float[:]) -> (float[:], float[], ...)
def lstsq(
    a: onp.ToFloatStrict2D,
    b: onp.ToFloatStrict1D,
    cond: onp.ToFloat | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriver | None = None,
) -> tuple[_Float1D, _Float0D, int, _Float1D | None]: ...
@overload  # (float[:, :], float[:, :]) -> (float[:, :], float[:], ...)
def lstsq(
    a: onp.ToFloatND,
    b: onp.ToFloatStrict2D,
    cond: onp.ToFloat | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriver | None = None,
) -> tuple[_FloatND, _FloatND, int, _FloatND | None]: ...
@overload  # (float[:, :], float[:, :?]) -> (float[:, :?], float[:?], ...)
def lstsq(
    a: onp.ToFloatND,
    b: onp.ToFloatND,
    cond: onp.ToFloat | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriver | None = None,
) -> tuple[_FloatND, _Float0D | _FloatND, int, _FloatND | None]: ...
@overload  # (complex[:, :], complex[:, :?]) -> (complex[:, :?], complex[:?], ...)
def lstsq(
    a: onp.ToComplexND,
    b: onp.ToComplexND,
    cond: onp.ToFloat | None = None,
    overwrite_a: bool = False,
    overwrite_b: bool = False,
    check_finite: bool = True,
    lapack_driver: _LapackDriver | None = None,
) -> tuple[_InexactND, _Inexact0D | _InexactND, int, _InexactND | None]: ...

# TODO(jorenham): improve this
@overload
def pinv(  # (float[:, :], return_rank=False) -> float[:, :]
    a: onp.ToFloatND,
    *,
    atol: onp.ToFloat | None = None,
    rtol: onp.ToFloat | None = None,
    return_rank: onp.ToFalse = False,
    check_finite: bool = True,
) -> _FloatND: ...
@overload  # (float[:, :], return_rank=True) -> (float[:, :], int)
def pinv(
    a: onp.ToFloatND,
    *,
    atol: onp.ToFloat | None = None,
    rtol: onp.ToFloat | None = None,
    return_rank: onp.ToTrue,
    check_finite: bool = True,
) -> tuple[_FloatND, int]: ...
@overload  # (complex[:, :], return_rank=False) -> complex[:, :]
def pinv(
    a: onp.ToComplexND,
    *,
    atol: onp.ToFloat | None = None,
    rtol: onp.ToFloat | None = None,
    return_rank: onp.ToFalse = False,
    check_finite: bool = True,
) -> _InexactND: ...
@overload  # (complex[:, :], return_rank=True) -> (complex[:, :], int)
def pinv(
    a: onp.ToComplexND,
    *,
    atol: onp.ToFloat | None = None,
    rtol: onp.ToFloat | None = None,
    return_rank: onp.ToTrue,
    check_finite: bool = True,
) -> tuple[_InexactND, int]: ...

# TODO(jorenham): improve this
@overload  # (float[:, :], return_rank=False) -> float[:, :]
def pinvh(
    a: onp.ToFloatND,
    atol: onp.ToFloat | None = None,
    rtol: onp.ToFloat | None = None,
    lower: bool = True,
    return_rank: onp.ToFalse = False,
    check_finite: bool = True,
) -> _FloatND: ...
@overload  # (float[:, :], return_rank=True, /) -> (float[:, :], int)
def pinvh(
    a: onp.ToFloatND,
    atol: onp.ToFloat | None,
    rtol: onp.ToFloat | None,
    lower: bool,
    return_rank: onp.ToTrue,
    check_finite: bool = True,
) -> tuple[_FloatND, int]: ...
@overload  # (float[:, :], *, return_rank=True) -> (float[:, :], int)
def pinvh(
    a: onp.ToFloatND,
    atol: onp.ToFloat | None = None,
    rtol: onp.ToFloat | None = None,
    lower: bool = True,
    *,
    return_rank: onp.ToTrue,
    check_finite: bool = True,
) -> tuple[_FloatND, int]: ...
@overload  # (complex[:, :], return_rank=False) -> complex[:, :]
def pinvh(
    a: onp.ToComplexND,
    atol: onp.ToFloat | None = None,
    rtol: onp.ToFloat | None = None,
    lower: bool = True,
    return_rank: onp.ToFalse = False,
    check_finite: bool = True,
) -> _InexactND: ...
@overload  # (complex[:, :], return_rank=True, /) -> (complex[:, :], int)
def pinvh(
    a: onp.ToComplexND,
    atol: onp.ToFloat | None,
    rtol: onp.ToFloat | None,
    lower: bool,
    return_rank: onp.ToTrue,
    check_finite: bool = True,
) -> tuple[_InexactND, int]: ...
@overload  # (complex[:, :], *, return_rank=True) -> (complex[:, :], int)
def pinvh(
    a: onp.ToComplexND,
    atol: onp.ToFloat | None = None,
    rtol: onp.ToFloat | None = None,
    lower: bool = True,
    *,
    return_rank: onp.ToTrue,
    check_finite: bool = True,
) -> tuple[_InexactND, int]: ...

# TODO(jorenham): improve this
@overload  # (float[:, :], separate=True) -> (float[:, :], float[:, :])
def matrix_balance(
    A: onp.ToFloatND,
    permute: onp.ToBool = True,
    scale: onp.ToBool = True,
    separate: onp.ToFalse = False,
    overwrite_a: bool = False,
) -> _Tuple2[_FloatND]: ...
@overload  # (float[:, :], separate=False, /) -> (float[:, :], (float[:], float[:]))
def matrix_balance(
    A: onp.ToFloatND, permute: onp.ToBool, scale: onp.ToBool, separate: onp.ToTrue, overwrite_a: bool = False
) -> tuple[_FloatND, _Tuple2[_FloatND]]: ...
@overload  # (float[:, :], *, separate=False) -> (float[:, :], (float[:], float[:]))
def matrix_balance(
    A: onp.ToFloatND, permute: onp.ToBool = True, scale: onp.ToBool = True, *, separate: onp.ToTrue, overwrite_a: bool = False
) -> tuple[_FloatND, _Tuple2[_FloatND]]: ...
@overload  # (complex[:, :], separate=True) -> (complex[:, :], complex[:, :])
def matrix_balance(
    A: onp.ToComplexND,
    permute: onp.ToBool = True,
    scale: onp.ToBool = True,
    separate: onp.ToFalse = False,
    overwrite_a: bool = False,
) -> _Tuple2[_InexactND]: ...
@overload  # (complex[:, :], separate=False, /) -> (complex[:, :], (complex[:], complex[:]))
def matrix_balance(
    A: onp.ToComplexND, permute: onp.ToBool, scale: onp.ToBool, separate: onp.ToTrue, overwrite_a: bool = False
) -> tuple[_InexactND, _Tuple2[_InexactND]]: ...
@overload  # (complex[:, :], *, separate=False) -> (complex[:, :], (complex[:], complex[:]))
def matrix_balance(
    A: onp.ToComplexND, permute: onp.ToBool = True, scale: onp.ToBool = True, *, separate: onp.ToTrue, overwrite_a: bool = False
) -> tuple[_InexactND, _Tuple2[_InexactND]]: ...

# TODO(jorenham): improve this
@overload  # floating 1d, 1d
def matmul_toeplitz(
    c_or_cr: onp.ToFloatStrict1D | _Tuple2[onp.ToFloatStrict1D],
    x: onp.ToFloatStrict1D,
    check_finite: bool = False,
    workers: onp.ToJustInt | None = None,
) -> _Float1D: ...
@overload  # floating
def matmul_toeplitz(
    c_or_cr: onp.ToFloatND | _Tuple2[onp.ToFloatND],
    x: onp.ToFloatND,
    check_finite: bool = False,
    workers: onp.ToJustInt | None = None,
) -> _FloatND: ...
@overload  # complexfloating 1d, 1d
def matmul_toeplitz(
    c_or_cr: onp.ToComplexStrict1D | _Tuple2[onp.ToComplexStrict1D],
    x: onp.ToComplexStrict1D,
    check_finite: bool = False,
    workers: onp.ToJustInt | None = None,
) -> _Inexact1D: ...
@overload  # complexfloating
def matmul_toeplitz(
    c_or_cr: onp.ToComplexND | _Tuple2[onp.ToComplexND],
    x: onp.ToComplexND,
    check_finite: bool = False,
    workers: onp.ToJustInt | None = None,
) -> _InexactND: ...
