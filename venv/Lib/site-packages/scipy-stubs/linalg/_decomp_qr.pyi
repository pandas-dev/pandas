from typing import Literal, overload

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["qr", "qr_multiply", "rq"]

type _Tuple2[T] = tuple[T, T]

type _Int1D = onp.Array1D[np.int32 | np.int64]
type _IntND = onp.ArrayND[np.int32 | np.int64]
type _Float1D = onp.Array1D[npc.floating]
type _Float2D = onp.Array2D[npc.floating]
type _FloatND = onp.ArrayND[npc.floating]
type _Inexact1D = onp.Array1D[npc.inexact]
type _Inexact2D = onp.Array2D[npc.inexact]
type _InexactND = onp.ArrayND[npc.inexact]

type _Side = Literal["left", "right"]
type _ModeFullEcon = Literal["full", "economic"]
type _ModeR = Literal["r"]
type _ModeRaw = Literal["raw"]

type _NoValueType = op.JustObject

###

# 2 * (3 + 4 + 4) = 22 overloads (10/22 handle the positional cases of `mode`/`pivoting`)
@overload  # float, mode: {full, economic}, pivoting: {False}
def qr(
    a: onp.ToFloatND,
    overwrite_a: bool = False,
    lwork: _NoValueType = ...,
    mode: _ModeFullEcon = "full",
    pivoting: onp.ToFalse = False,
    check_finite: bool = True,
) -> _Tuple2[_FloatND]: ...
@overload  # float, mode: {full, economic}, *, pivoting: {True}
def qr(
    a: onp.ToFloatND,
    overwrite_a: bool = False,
    lwork: _NoValueType = ...,
    mode: _ModeFullEcon = "full",
    *,
    pivoting: onp.ToTrue,
    check_finite: bool = True,
) -> tuple[_FloatND, _FloatND, _IntND]: ...
@overload  # float, *, mode: {r}, pivoting: {False}
def qr(
    a: onp.ToFloatND,
    overwrite_a: bool = False,
    lwork: _NoValueType = ...,
    *,
    mode: _ModeR,
    pivoting: onp.ToFalse = False,
    check_finite: bool = True,
) -> tuple[_FloatND]: ...
@overload  # float, *, mode: {r}, pivoting: {True}
def qr(
    a: onp.ToFloatND,
    overwrite_a: bool = False,
    lwork: _NoValueType = ...,
    *,
    mode: _ModeR,
    pivoting: onp.ToTrue,
    check_finite: bool = True,
) -> tuple[_FloatND, _IntND]: ...
@overload  # float, *, mode: {raw}, pivoting: {False}
def qr(
    a: onp.ToFloatND,
    overwrite_a: bool = False,
    lwork: _NoValueType = ...,
    *,
    mode: _ModeRaw,
    pivoting: onp.ToFalse = False,
    check_finite: bool = True,
) -> tuple[_Tuple2[_FloatND], _FloatND]: ...
@overload  # float, *, mode: {raw}, pivoting: {True}
def qr(
    a: onp.ToFloatND,
    overwrite_a: bool = False,
    lwork: _NoValueType = ...,
    *,
    mode: _ModeRaw,
    pivoting: onp.ToTrue,
    check_finite: bool = True,
) -> tuple[_Tuple2[_FloatND], _FloatND, _IntND]: ...
@overload  # complex, mode: {full, economic}, pivoting: {False}
def qr(
    a: onp.ToComplexND,
    overwrite_a: bool = False,
    lwork: _NoValueType = ...,
    mode: _ModeFullEcon = "full",
    pivoting: onp.ToFalse = False,
    check_finite: bool = True,
) -> _Tuple2[_InexactND]: ...
@overload  # complex, mode: {full, economic}, *, pivoting: {True}
def qr(
    a: onp.ToComplexND,
    overwrite_a: bool = False,
    lwork: _NoValueType = ...,
    mode: _ModeFullEcon = "full",
    *,
    pivoting: onp.ToTrue,
    check_finite: bool = True,
) -> tuple[_InexactND, _InexactND, _IntND]: ...
@overload  # complex, *, mode: {r}, pivoting: {False}
def qr(
    a: onp.ToComplexND,
    overwrite_a: bool = False,
    lwork: _NoValueType = ...,
    *,
    mode: _ModeR,
    pivoting: onp.ToFalse = False,
    check_finite: bool = True,
) -> tuple[_InexactND]: ...
@overload  # complex, *, mode: {r}, pivoting: {True}
def qr(
    a: onp.ToComplexND,
    overwrite_a: bool = False,
    lwork: _NoValueType = ...,
    *,
    mode: _ModeR,
    pivoting: onp.ToTrue,
    check_finite: bool = True,
) -> tuple[_InexactND, _IntND]: ...
@overload  # complex, *, mode: {raw}, pivoting: {False}
def qr(
    a: onp.ToComplexND,
    overwrite_a: bool = False,
    lwork: _NoValueType = ...,
    *,
    mode: _ModeRaw,
    pivoting: onp.ToFalse = False,
    check_finite: bool = True,
) -> tuple[_Tuple2[_InexactND], _InexactND]: ...
@overload  # complex, *, mode: {raw}, pivoting: {True}
def qr(
    a: onp.ToComplexND,
    overwrite_a: bool = False,
    lwork: _NoValueType = ...,
    *,
    mode: _ModeRaw,
    pivoting: onp.ToTrue,
    check_finite: bool = True,
) -> tuple[_Tuple2[_InexactND], _InexactND, _IntND]: ...

#
@overload  # (float[:, :], float[:], pivoting=False) -> (float[:], float[:, :])
def qr_multiply(
    a: onp.ToFloatStrict2D,
    c: onp.ToFloatStrict1D,
    mode: _Side = "right",
    pivoting: onp.ToFalse = False,
    conjugate: bool = False,
    overwrite_a: bool = False,
    overwrite_c: bool = False,
) -> tuple[_Float1D, _Inexact2D]: ...
@overload  # (float[:, :], float[:, :], pivoting=False) -> (float[:, :], float[:, :])
def qr_multiply(
    a: onp.ToFloatStrict2D,
    c: onp.ToFloatStrict2D,
    mode: _Side = "right",
    pivoting: onp.ToFalse = False,
    conjugate: bool = False,
    overwrite_a: bool = False,
    overwrite_c: bool = False,
) -> tuple[_Float2D, _Inexact2D]: ...
@overload  # (float[:, :], float[:, :?], pivoting=False) -> (float[:, :?], float[:, :])
def qr_multiply(
    a: onp.ToFloatStrict2D,
    c: onp.ToFloatStrict1D | onp.ToFloatStrict2D,
    mode: _Side = "right",
    pivoting: onp.ToFalse = False,
    conjugate: bool = False,
    overwrite_a: bool = False,
    overwrite_c: bool = False,
) -> tuple[_Float1D | _Float2D, _Inexact2D]: ...
@overload  # (float[:, :, ...], float[:, ...], pivoting=False) -> (float[:, ...], float[:, :, ...])
def qr_multiply(
    a: onp.ToFloatND,
    c: onp.ToFloatND,
    mode: _Side = "right",
    pivoting: onp.ToFalse = False,
    conjugate: bool = False,
    overwrite_a: bool = False,
    overwrite_c: bool = False,
) -> tuple[_FloatND, _InexactND]: ...
@overload  # (float[:, :], float[:, :?], pivoting=True, /) -> (float[:, :?], float[:, :], int[:])
def qr_multiply(
    a: onp.ToFloatStrict2D,
    c: onp.ToFloatStrict1D | onp.ToFloatStrict2D,
    mode: _Side,
    pivoting: onp.ToTrue,
    conjugate: bool = False,
    overwrite_a: bool = False,
    overwrite_c: bool = False,
) -> tuple[_Float1D | _Float2D, _Float2D, _Int1D]: ...
@overload  # (float[:, :], float[:, :?], *, pivoting=True) -> (float[:, :?], float[:, :], int[:])
def qr_multiply(
    a: onp.ToFloatStrict2D,
    c: onp.ToFloatStrict1D | onp.ToFloatStrict2D,
    mode: _Side = "right",
    *,
    pivoting: onp.ToTrue,
    conjugate: bool = False,
    overwrite_a: bool = False,
    overwrite_c: bool = False,
) -> tuple[_Float1D | _Float2D, _Float2D, _Int1D]: ...
@overload  # (float[:, :, ...], float[:, ...], *, pivoting=True) -> (float[:, ...], float[:, :, ...], int[:, ...])
def qr_multiply(
    a: onp.ToFloatND,
    c: onp.ToFloatND,
    mode: _Side = "right",
    *,
    pivoting: onp.ToTrue,
    conjugate: bool = False,
    overwrite_a: bool = False,
    overwrite_c: bool = False,
) -> tuple[_FloatND, _FloatND, _IntND]: ...
@overload  # (complex[:, :], complex[:, :?], pivoting=False) -> (complex[:, :?], complex[:, :])
def qr_multiply(
    a: onp.ToComplexStrict2D,
    c: onp.ToComplexStrict1D | onp.ToComplexStrict2D,
    mode: _Side = "right",
    pivoting: onp.ToFalse = False,
    conjugate: bool = False,
    overwrite_a: bool = False,
    overwrite_c: bool = False,
) -> tuple[_Inexact1D | _Inexact2D, _Inexact2D]: ...
@overload  # (complex[:, :, ...], complex[:, ...], pivoting=False) -> (complex[:, ...], complex[:, :, ...])
def qr_multiply(
    a: onp.ToComplexND,
    c: onp.ToComplexND,
    mode: _Side = "right",
    pivoting: onp.ToFalse = False,
    conjugate: bool = False,
    overwrite_a: bool = False,
    overwrite_c: bool = False,
) -> tuple[_InexactND, _InexactND]: ...
@overload  # (complex[:, :], complex[:, :?], pivoting=True, /) -> (complex[:, :?], complex[:, :], int[:])
def qr_multiply(
    a: onp.ToComplexStrict2D,
    c: onp.ToComplexStrict1D | onp.ToComplexStrict2D,
    mode: _Side,
    pivoting: onp.ToTrue,
    conjugate: bool = False,
    overwrite_a: bool = False,
    overwrite_c: bool = False,
) -> tuple[_Inexact1D | _Inexact2D, _Inexact2D, _Int1D]: ...
@overload  # (complex[:, :], complex[:, :?], *, pivoting=True) -> (complex[:, :?], complex[:, :], int[:])
def qr_multiply(
    a: onp.ToComplexStrict2D,
    c: onp.ToComplexStrict1D | onp.ToComplexStrict2D,
    mode: _Side = "right",
    *,
    pivoting: onp.ToTrue,
    conjugate: bool = False,
    overwrite_a: bool = False,
    overwrite_c: bool = False,
) -> tuple[_Inexact1D | _Inexact2D, _Inexact2D, _Int1D]: ...
@overload  # (complex[:, :, ...], complex[:, ...], *, pivoting=True) -> (complex[:, ...], complex[:, :, ...], int[:, ...])
def qr_multiply(
    a: onp.ToComplexND,
    c: onp.ToComplexND,
    mode: _Side = "right",
    *,
    pivoting: onp.ToTrue,
    conjugate: bool = False,
    overwrite_a: bool = False,
    overwrite_c: bool = False,
) -> tuple[_InexactND, _InexactND, _IntND]: ...

#
@overload  # (float[:, :], mode: {"full", "economic"}) -> (float[:, :], float[:, :])
def rq(
    a: onp.ToFloatStrict2D,
    overwrite_a: bool = False,
    lwork: onp.ToJustInt | None = None,
    mode: _ModeFullEcon = "full",
    check_finite: bool = True,
) -> tuple[_Float2D, _Float2D]: ...
@overload  # (float[:, :, ...], mode: {"full", "economic"}) -> (float[:, :, ...], float[:, :, ...])
def rq(
    a: onp.ToFloatND,
    overwrite_a: bool = False,
    lwork: onp.ToJustInt | None = None,
    mode: _ModeFullEcon = "full",
    check_finite: bool = True,
) -> tuple[_FloatND, _FloatND]: ...
@overload  # (float[:, :], mode: {"r"}, /) -> float[:, :]
def rq(
    a: onp.ToFloatStrict2D, overwrite_a: bool, lwork: onp.ToJustInt | None, mode: _ModeR, check_finite: bool = True
) -> _Float2D: ...
@overload  # (float[:, :, ...], mode: {"r"}, /) -> float[:, :, ...]
def rq(a: onp.ToFloatND, overwrite_a: bool, lwork: onp.ToJustInt | None, mode: _ModeR, check_finite: bool = True) -> _FloatND: ...
@overload  # (float[:, :], *, mode: {"r"}) -> float[:, :]
def rq(
    a: onp.ToFloatStrict2D,
    overwrite_a: bool = False,
    lwork: onp.ToJustInt | None = None,
    *,
    mode: _ModeR,
    check_finite: bool = True,
) -> _Float2D: ...
@overload  # (float[:, :, ...], *, mode: {"r"}) -> float[:, : ...]
def rq(
    a: onp.ToFloatND, overwrite_a: bool = False, lwork: onp.ToJustInt | None = None, *, mode: _ModeR, check_finite: bool = True
) -> _FloatND: ...
@overload  # (complex[:, :], mode: {"full", "economic"}) -> (complex[:, :], complex[:, :])
def rq(
    a: onp.ToComplexStrict2D,
    overwrite_a: bool = False,
    lwork: onp.ToJustInt | None = None,
    mode: _ModeFullEcon = "full",
    check_finite: bool = True,
) -> _Tuple2[_Inexact2D]: ...
@overload  # (complex[:, :], mode: {"full", "economic"}) -> (complex[:, :], complex[:, :])
def rq(
    a: onp.ToComplexND,
    overwrite_a: bool = False,
    lwork: onp.ToJustInt | None = None,
    mode: _ModeFullEcon = "full",
    check_finite: bool = True,
) -> _Tuple2[_InexactND]: ...
@overload  # (complex[:, :], mode: {"r"}, /) -> complex[:, :]
def rq(
    a: onp.ToComplexStrict2D, overwrite_a: bool, lwork: onp.ToJustInt | None, mode: _ModeR, check_finite: bool = True
) -> _Inexact2D: ...
@overload  # (complex[:, :, ...], mode: {"r"}, /) -> complex[:, :, ...]
def rq(
    a: onp.ToComplexND, overwrite_a: bool, lwork: onp.ToJustInt | None, mode: _ModeR, check_finite: bool = True
) -> _InexactND: ...
@overload  # (complex[:, :], *, mode: {"r"}) -> complex[:, :]
def rq(
    a: onp.ToComplexStrict2D,
    overwrite_a: bool = False,
    lwork: onp.ToJustInt | None = None,
    *,
    mode: _ModeR,
    check_finite: bool = True,
) -> _Inexact2D: ...
@overload  # (complex[:, :, ...], *, mode: {"r"}) -> complex[:, :, ...]
def rq(
    a: onp.ToComplexND, overwrite_a: bool = False, lwork: onp.ToJustInt | None = None, *, mode: _ModeR, check_finite: bool = True
) -> _InexactND: ...
