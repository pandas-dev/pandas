from collections.abc import Callable
from typing import Literal, overload

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["rsf2csf", "schur"]

###

type _Tuple2[T] = tuple[T, T]
type _Tuple2i[T] = tuple[T, T, int]

type _OutputReal = Literal["real", "r"]
type _OutputComplex = Literal["complex", "c"]
type _Output = Literal[_OutputReal, _OutputComplex]

type _Sort = Literal["lhp", "rhp", "iuc", "ouc"] | Callable[[float, float], bool]

type _as_f32 = np.float32 | np.float16 | np.bool  # noqa: PYI042
type _as_f64 = npc.floating64 | npc.floating80 | npc.integer  # noqa: PYI042
type _as_c128 = npc.complexfloating128 | npc.complexfloating160  # noqa: PYI042

###

# NOTE: On numpy<2.1, pyright reports 12 false positive incompatible overload errors here.
# pyright: reportOverlappingOverload=false

# NOTE: 1/10 of mypy's false positive `overload-overlap` errors in this module are only reported with numpy<2.5
# mypy: disable-error-code=overload-overlap

#
@overload  # Nd f64
def schur[Shape2T: tuple[int, int, *tuple[int, ...]]](
    a: onp.ArrayND[_as_f64, Shape2T],
    output: _OutputReal = "real",
    lwork: int | None = None,
    overwrite_a: bool = False,
    sort: None = None,
    check_finite: bool = True,
) -> _Tuple2[onp.ArrayND[np.float64, Shape2T]]: ...
@overload  # Nd f64, sort=<given>
def schur[Shape2T: tuple[int, int, *tuple[int, ...]]](
    a: onp.ArrayND[_as_f64, Shape2T],
    output: _OutputReal = "real",
    lwork: int | None = None,
    overwrite_a: bool = False,
    *,
    sort: _Sort,
    check_finite: bool = True,
) -> _Tuple2i[onp.ArrayND[np.float64, Shape2T]]: ...
@overload  # Nd f64, output="complex"
def schur[Shape2T: tuple[int, int, *tuple[int, ...]]](
    a: onp.ArrayND[_as_f64, Shape2T],
    output: _OutputComplex,
    lwork: int | None = None,
    overwrite_a: bool = False,
    sort: None = None,
    check_finite: bool = True,
) -> _Tuple2[onp.ArrayND[np.complex128, Shape2T]]: ...
@overload  # Nd f64, output="complex", sort=<given>
def schur[Shape2T: tuple[int, int, *tuple[int, ...]]](
    a: onp.ArrayND[_as_f64, Shape2T],
    output: _OutputComplex,
    lwork: int | None = None,
    overwrite_a: bool = False,
    *,
    sort: _Sort,
    check_finite: bool = True,
) -> _Tuple2i[onp.ArrayND[np.complex128, Shape2T]]: ...
@overload  # Nd f32
def schur[Shape2T: tuple[int, int, *tuple[int, ...]]](
    a: onp.ArrayND[_as_f32, Shape2T],
    output: _OutputReal = "real",
    lwork: int | None = None,
    overwrite_a: bool = False,
    sort: None = None,
    check_finite: bool = True,
) -> _Tuple2[onp.ArrayND[np.float32, Shape2T]]: ...
@overload  # Nd f32, sort=<given>
def schur[Shape2T: tuple[int, int, *tuple[int, ...]]](
    a: onp.ArrayND[_as_f32, Shape2T],
    output: _OutputReal = "real",
    lwork: int | None = None,
    overwrite_a: bool = False,
    *,
    sort: _Sort,
    check_finite: bool = True,
) -> _Tuple2i[onp.ArrayND[np.float32, Shape2T]]: ...
@overload  # Nd f32, output="complex"
def schur[Shape2T: tuple[int, int, *tuple[int, ...]]](
    a: onp.ArrayND[_as_f32, Shape2T],
    output: _OutputComplex,
    lwork: int | None = None,
    overwrite_a: bool = False,
    sort: None = None,
    check_finite: bool = True,
) -> _Tuple2[onp.ArrayND[np.complex64, Shape2T]]: ...
@overload  # Nd f32, output="complex", sort=<given>
def schur[Shape2T: tuple[int, int, *tuple[int, ...]]](
    a: onp.ArrayND[_as_f32, Shape2T],
    output: _OutputComplex,
    lwork: int | None = None,
    overwrite_a: bool = False,
    *,
    sort: _Sort,
    check_finite: bool = True,
) -> _Tuple2i[onp.ArrayND[np.complex64, Shape2T]]: ...
@overload  # Nd c128
def schur[Shape2T: tuple[int, int, *tuple[int, ...]]](
    a: onp.ArrayND[_as_c128, Shape2T],
    output: _Output = "real",
    lwork: int | None = None,
    overwrite_a: bool = False,
    sort: None = None,
    check_finite: bool = True,
) -> _Tuple2[onp.ArrayND[np.complex128, Shape2T]]: ...
@overload  # Nd c128, sort=<given>
def schur[Shape2T: tuple[int, int, *tuple[int, ...]]](
    a: onp.ArrayND[_as_c128, Shape2T],
    output: _Output = "real",
    lwork: int | None = None,
    overwrite_a: bool = False,
    *,
    sort: _Sort,
    check_finite: bool = True,
) -> _Tuple2i[onp.ArrayND[np.complex128, Shape2T]]: ...
@overload  # Nd c64
def schur[Shape2T: tuple[int, int, *tuple[int, ...]]](
    a: onp.ArrayND[np.complex64, Shape2T],
    output: _Output = "real",
    lwork: int | None = None,
    overwrite_a: bool = False,
    sort: None = None,
    check_finite: bool = True,
) -> _Tuple2[onp.ArrayND[np.complex64, Shape2T]]: ...
@overload  # Nd c64, sort=<given>
def schur[Shape2T: tuple[int, int, *tuple[int, ...]]](
    a: onp.ArrayND[np.complex64, Shape2T],
    output: _Output = "real",
    lwork: int | None = None,
    overwrite_a: bool = False,
    *,
    sort: _Sort,
    check_finite: bool = True,
) -> _Tuple2i[onp.ArrayND[np.complex64, Shape2T]]: ...
@overload  # ?d f64
def schur(
    a: onp.ToArrayND[float, _as_f64],
    output: _OutputReal = "real",
    lwork: int | None = None,
    overwrite_a: bool = False,
    sort: None = None,
    check_finite: bool = True,
) -> _Tuple2[onp.ArrayND[np.float64]]: ...
@overload  # ?d f64, sort=<given>
def schur(
    a: onp.ToArrayND[float, _as_f64],
    output: _OutputReal = "real",
    lwork: int | None = None,
    overwrite_a: bool = False,
    *,
    sort: _Sort,
    check_finite: bool = True,
) -> _Tuple2i[onp.ArrayND[np.float64]]: ...
@overload  # ?d f64, output="complex"
def schur(
    a: onp.ToArrayND[float, _as_f64],
    output: _OutputComplex,
    lwork: int | None = None,
    overwrite_a: bool = False,
    sort: None = None,
    check_finite: bool = True,
) -> _Tuple2[onp.ArrayND[np.complex128]]: ...
@overload  # ?d f64, output="complex", sort=<given>
def schur(
    a: onp.ToArrayND[float, _as_f64],
    output: _OutputComplex,
    lwork: int | None = None,
    overwrite_a: bool = False,
    *,
    sort: _Sort,
    check_finite: bool = True,
) -> _Tuple2i[onp.ArrayND[np.complex128]]: ...
@overload  # ?d f32
def schur(
    a: onp.ToArrayND[np.float32, _as_f32],
    output: _OutputReal = "real",
    lwork: int | None = None,
    overwrite_a: bool = False,
    sort: None = None,
    check_finite: bool = True,
) -> _Tuple2[onp.ArrayND[np.float32]]: ...
@overload  # ?d f32, sort=<given>
def schur(
    a: onp.ToArrayND[np.float32, _as_f32],
    output: _OutputReal = "real",
    lwork: int | None = None,
    overwrite_a: bool = False,
    *,
    sort: _Sort,
    check_finite: bool = True,
) -> _Tuple2i[onp.ArrayND[np.float32]]: ...
@overload  # ?d f32, output="complex"
def schur(
    a: onp.ToArrayND[np.float32, _as_f32],
    output: _OutputComplex,
    lwork: int | None = None,
    overwrite_a: bool = False,
    sort: None = None,
    check_finite: bool = True,
) -> _Tuple2[onp.ArrayND[np.complex64]]: ...
@overload  # ?d f32, output="complex", sort=<given>
def schur(
    a: onp.ToArrayND[np.float32, _as_f32],
    output: _OutputComplex,
    lwork: int | None = None,
    overwrite_a: bool = False,
    *,
    sort: _Sort,
    check_finite: bool = True,
) -> _Tuple2i[onp.ArrayND[np.complex64]]: ...
@overload  # ?d c128
def schur(
    a: onp.ToArrayND[op.JustComplex, _as_c128],
    output: _Output = "real",
    lwork: int | None = None,
    overwrite_a: bool = False,
    sort: None = None,
    check_finite: bool = True,
) -> _Tuple2[onp.ArrayND[np.complex128]]: ...
@overload  # ?d c128, sort=<given>
def schur(
    a: onp.ToArrayND[op.JustComplex, _as_c128],
    output: _Output = "real",
    lwork: int | None = None,
    overwrite_a: bool = False,
    *,
    sort: _Sort,
    check_finite: bool = True,
) -> _Tuple2i[onp.ArrayND[np.complex128]]: ...
@overload  # ?d c64
def schur(
    a: onp.ToJustComplex64_ND,
    output: _Output = "real",
    lwork: int | None = None,
    overwrite_a: bool = False,
    sort: None = None,
    check_finite: bool = True,
) -> _Tuple2[onp.ArrayND[np.complex64]]: ...
@overload  # ?d c64, sort=<given>
def schur(
    a: onp.ToJustComplex64_ND,
    output: _Output = "real",
    lwork: int | None = None,
    overwrite_a: bool = False,
    *,
    sort: _Sort,
    check_finite: bool = True,
) -> _Tuple2i[onp.ArrayND[np.complex64]]: ...

# will raise for dtypes that don't have character code in `ilfdFD`
@overload  # ?d c128|f64|i64|i32, ?d c128|f64|i64|i32
def rsf2csf(
    T: onp.ToArrayND[complex, npc.inexact64 | npc.integer64 | npc.integer32],
    Z: onp.ToArrayND[complex, npc.inexact64 | npc.integer64 | npc.integer32],
    check_finite: bool = True,
) -> _Tuple2[onp.ArrayND[np.complex128]]: ...
@overload  # ?d f32|c64, ?d f32|c64 -> c64
def rsf2csf(
    T: onp.ToArrayND[npc.inexact32, npc.inexact32], Z: onp.ToArrayND[npc.inexact32, npc.inexact32], check_finite: bool = True
) -> _Tuple2[onp.ArrayND[np.complex64]]: ...
