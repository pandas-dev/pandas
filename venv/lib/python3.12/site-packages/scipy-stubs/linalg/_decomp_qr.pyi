from typing import Literal, TypeAlias, TypeVar, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["qr", "qr_multiply", "rq"]

_T = TypeVar("_T")
_Tuple2: TypeAlias = tuple[_T, _T]

_Int1D: TypeAlias = onp.Array1D[np.int32 | np.int64]
_IntND: TypeAlias = onp.ArrayND[np.int32 | np.int64]
_Float1D: TypeAlias = onp.Array1D[npc.floating]
_Float2D: TypeAlias = onp.Array2D[npc.floating]
_FloatND: TypeAlias = onp.ArrayND[npc.floating]
_Inexact1D: TypeAlias = onp.Array1D[npc.inexact]
_Inexact2D: TypeAlias = onp.Array2D[npc.inexact]
_InexactND: TypeAlias = onp.ArrayND[npc.inexact]

_Side: TypeAlias = Literal["left", "right"]
_ModeFullEcon: TypeAlias = Literal["full", "economic"]
_ModeR: TypeAlias = Literal["r"]
_ModeRaw: TypeAlias = Literal["raw"]

###

# 2 * (3 + 4 + 4) = 22 overloads (10/22 handle the positional cases of `mode`/`pivoting`)
@overload  # float, mode: {full, economic}, pivoting: {False}
def qr(
    a: onp.ToFloatND,
    overwrite_a: onp.ToBool = False,
    lwork: onp.ToJustInt | None = None,
    mode: _ModeFullEcon = "full",
    pivoting: onp.ToFalse = False,
    check_finite: onp.ToBool = True,
) -> _Tuple2[_FloatND]: ...
@overload  # float, mode: {full, economic}, pivoting: {True}
def qr(
    a: onp.ToFloatND,
    overwrite_a: onp.ToBool,
    lwork: onp.ToJustInt | None,
    mode: _ModeFullEcon,
    pivoting: onp.ToTrue,
    check_finite: onp.ToBool = True,
) -> tuple[_FloatND, _FloatND, _IntND]: ...
@overload  # float, mode: {full, economic}, *, pivoting: {True}
def qr(
    a: onp.ToFloatND,
    overwrite_a: onp.ToBool = False,
    lwork: onp.ToJustInt | None = None,
    mode: _ModeFullEcon = "full",
    *,
    pivoting: onp.ToTrue,
    check_finite: onp.ToBool = True,
) -> tuple[_FloatND, _FloatND, _IntND]: ...
@overload  # float, mode: {r}, pivoting: {False}
def qr(
    a: onp.ToFloatND,
    overwrite_a: onp.ToBool,
    lwork: onp.ToJustInt | None,
    mode: _ModeR,
    pivoting: onp.ToFalse = False,
    check_finite: onp.ToBool = True,
) -> tuple[_FloatND]: ...
@overload  # float, mode: {r}, pivoting: {True}
def qr(
    a: onp.ToFloatND,
    overwrite_a: onp.ToBool,
    lwork: onp.ToJustInt | None,
    mode: _ModeR,
    pivoting: onp.ToTrue,
    check_finite: onp.ToBool = True,
) -> tuple[_FloatND, _IntND]: ...
@overload  # float, *, mode: {r}, pivoting: {False}
def qr(
    a: onp.ToFloatND,
    overwrite_a: onp.ToBool = False,
    lwork: onp.ToJustInt | None = None,
    *,
    mode: _ModeR,
    pivoting: onp.ToFalse = False,
    check_finite: onp.ToBool = True,
) -> tuple[_FloatND]: ...
@overload  # float, *, mode: {r}, pivoting: {True}
def qr(
    a: onp.ToFloatND,
    overwrite_a: onp.ToBool = False,
    lwork: onp.ToJustInt | None = None,
    *,
    mode: _ModeR,
    pivoting: onp.ToTrue,
    check_finite: onp.ToBool = True,
) -> tuple[_FloatND, _IntND]: ...
@overload  # float, mode: {raw}, pivoting: {False}
def qr(
    a: onp.ToFloatND,
    overwrite_a: onp.ToBool,
    lwork: onp.ToJustInt | None,
    mode: _ModeRaw,
    pivoting: onp.ToFalse = False,
    check_finite: onp.ToBool = True,
) -> tuple[_Tuple2[_FloatND], _FloatND]: ...
@overload  # float, mode: {raw}, pivoting: {True}
def qr(
    a: onp.ToFloatND,
    overwrite_a: onp.ToBool,
    lwork: onp.ToJustInt | None,
    mode: _ModeRaw,
    pivoting: onp.ToTrue,
    check_finite: onp.ToBool = True,
) -> tuple[_Tuple2[_FloatND], _FloatND, _IntND]: ...
@overload  # float, *, mode: {raw}, pivoting: {False}
def qr(
    a: onp.ToFloatND,
    overwrite_a: onp.ToBool = False,
    lwork: onp.ToJustInt | None = None,
    *,
    mode: _ModeRaw,
    pivoting: onp.ToFalse = False,
    check_finite: onp.ToBool = True,
) -> tuple[_Tuple2[_FloatND], _FloatND]: ...
@overload  # float, *, mode: {raw}, pivoting: {True}
def qr(
    a: onp.ToFloatND,
    overwrite_a: onp.ToBool = False,
    lwork: onp.ToJustInt | None = None,
    *,
    mode: _ModeRaw,
    pivoting: onp.ToTrue,
    check_finite: onp.ToBool = True,
) -> tuple[_Tuple2[_FloatND], _FloatND, _IntND]: ...
@overload  # complex, mode: {full, economic}, pivoting: {False}
def qr(
    a: onp.ToComplexND,
    overwrite_a: onp.ToBool = False,
    lwork: onp.ToJustInt | None = None,
    mode: _ModeFullEcon = "full",
    pivoting: onp.ToFalse = False,
    check_finite: onp.ToBool = True,
) -> _Tuple2[_InexactND]: ...
@overload  # complex, mode: {full, economic}, pivoting: {True}
def qr(
    a: onp.ToComplexND,
    overwrite_a: onp.ToBool,
    lwork: onp.ToJustInt | None,
    mode: _ModeFullEcon,
    pivoting: onp.ToTrue,
    check_finite: onp.ToBool = True,
) -> tuple[_InexactND, _InexactND, _IntND]: ...
@overload  # complex, mode: {full, economic}, *, pivoting: {True}
def qr(
    a: onp.ToComplexND,
    overwrite_a: onp.ToBool = False,
    lwork: onp.ToJustInt | None = None,
    mode: _ModeFullEcon = "full",
    *,
    pivoting: onp.ToTrue,
    check_finite: onp.ToBool = True,
) -> tuple[_InexactND, _InexactND, _IntND]: ...
@overload  # complex, mode: {r}, pivoting: {False}
def qr(
    a: onp.ToComplexND,
    overwrite_a: onp.ToBool,
    lwork: onp.ToJustInt | None,
    mode: _ModeR,
    pivoting: onp.ToFalse = False,
    check_finite: onp.ToBool = True,
) -> tuple[_InexactND]: ...
@overload  # complex, mode: {r}, pivoting: {True}
def qr(
    a: onp.ToComplexND,
    overwrite_a: onp.ToBool,
    lwork: onp.ToJustInt | None,
    mode: _ModeR,
    pivoting: onp.ToTrue,
    check_finite: onp.ToBool = True,
) -> tuple[_InexactND, _IntND]: ...
@overload  # complex, *, mode: {r}, pivoting: {False}
def qr(
    a: onp.ToComplexND,
    overwrite_a: onp.ToBool = False,
    lwork: onp.ToJustInt | None = None,
    *,
    mode: _ModeR,
    pivoting: onp.ToFalse = False,
    check_finite: onp.ToBool = True,
) -> tuple[_InexactND]: ...
@overload  # complex, *, mode: {r}, pivoting: {True}
def qr(
    a: onp.ToComplexND,
    overwrite_a: onp.ToBool = False,
    lwork: onp.ToJustInt | None = None,
    *,
    mode: _ModeR,
    pivoting: onp.ToTrue,
    check_finite: onp.ToBool = True,
) -> tuple[_InexactND, _IntND]: ...
@overload  # complex, mode: {raw}, pivoting: {False}
def qr(
    a: onp.ToComplexND,
    overwrite_a: onp.ToBool,
    lwork: onp.ToJustInt | None,
    mode: _ModeRaw,
    pivoting: onp.ToFalse = False,
    check_finite: onp.ToBool = True,
) -> tuple[_Tuple2[_InexactND], _InexactND]: ...
@overload  # complex, mode: {raw}, pivoting: {True}
def qr(
    a: onp.ToComplexND,
    overwrite_a: onp.ToBool,
    lwork: onp.ToJustInt | None,
    mode: _ModeRaw,
    pivoting: onp.ToTrue,
    check_finite: onp.ToBool = True,
) -> tuple[_Tuple2[_InexactND], _InexactND, _IntND]: ...
@overload  # complex, *, mode: {raw}, pivoting: {False}
def qr(
    a: onp.ToComplexND,
    overwrite_a: onp.ToBool = False,
    lwork: onp.ToJustInt | None = None,
    *,
    mode: _ModeRaw,
    pivoting: onp.ToFalse = False,
    check_finite: onp.ToBool = True,
) -> tuple[_Tuple2[_InexactND], _InexactND]: ...
@overload  # complex, *, mode: {raw}, pivoting: {True}
def qr(
    a: onp.ToComplexND,
    overwrite_a: onp.ToBool = False,
    lwork: onp.ToJustInt | None = None,
    *,
    mode: _ModeRaw,
    pivoting: onp.ToTrue,
    check_finite: onp.ToBool = True,
) -> tuple[_Tuple2[_InexactND], _InexactND, _IntND]: ...

#
@overload  # (float[:, :], float[:], pivoting=False) -> (float[:], float[:, :])
def qr_multiply(
    a: onp.ToFloatStrict2D,
    c: onp.ToFloatStrict1D,
    mode: _Side = "right",
    pivoting: onp.ToFalse = False,
    conjugate: onp.ToBool = False,
    overwrite_a: onp.ToBool = False,
    overwrite_c: onp.ToBool = False,
) -> tuple[_Float1D, _Inexact2D]: ...
@overload  # (float[:, :], float[:, :], pivoting=False) -> (float[:, :], float[:, :])
def qr_multiply(
    a: onp.ToFloatStrict2D,
    c: onp.ToFloatStrict2D,
    mode: _Side = "right",
    pivoting: onp.ToFalse = False,
    conjugate: onp.ToBool = False,
    overwrite_a: onp.ToBool = False,
    overwrite_c: onp.ToBool = False,
) -> tuple[_Float2D, _Inexact2D]: ...
@overload  # (float[:, :], float[:, :?], pivoting=False) -> (float[:, :?], float[:, :])
def qr_multiply(
    a: onp.ToFloatStrict2D,
    c: onp.ToFloatStrict1D | onp.ToFloatStrict2D,
    mode: _Side = "right",
    pivoting: onp.ToFalse = False,
    conjugate: onp.ToBool = False,
    overwrite_a: onp.ToBool = False,
    overwrite_c: onp.ToBool = False,
) -> tuple[_Float1D | _Float2D, _Inexact2D]: ...
@overload  # (float[:, :, ...], float[:, ...], pivoting=False) -> (float[:, ...], float[:, :, ...])
def qr_multiply(
    a: onp.ToFloatND,
    c: onp.ToFloatND,
    mode: _Side = "right",
    pivoting: onp.ToFalse = False,
    conjugate: onp.ToBool = False,
    overwrite_a: onp.ToBool = False,
    overwrite_c: onp.ToBool = False,
) -> tuple[_FloatND, _InexactND]: ...
@overload  # (float[:, :], float[:, :?], pivoting=True, /) -> (float[:, :?], float[:, :], int[:])
def qr_multiply(
    a: onp.ToFloatStrict2D,
    c: onp.ToFloatStrict1D | onp.ToFloatStrict2D,
    mode: _Side,
    pivoting: onp.ToTrue,
    conjugate: onp.ToBool = False,
    overwrite_a: onp.ToBool = False,
    overwrite_c: onp.ToBool = False,
) -> tuple[_Float1D | _Float2D, _Float2D, _Int1D]: ...
@overload  # (float[:, :], float[:, :?], *, pivoting=True) -> (float[:, :?], float[:, :], int[:])
def qr_multiply(
    a: onp.ToFloatStrict2D,
    c: onp.ToFloatStrict1D | onp.ToFloatStrict2D,
    mode: _Side = "right",
    *,
    pivoting: onp.ToTrue,
    conjugate: onp.ToBool = False,
    overwrite_a: onp.ToBool = False,
    overwrite_c: onp.ToBool = False,
) -> tuple[_Float1D | _Float2D, _Float2D, _Int1D]: ...
@overload  # (float[:, :, ...], float[:, ...], *, pivoting=True) -> (float[:, ...], float[:, :, ...], int[:, ...])
def qr_multiply(
    a: onp.ToFloatND,
    c: onp.ToFloatND,
    mode: _Side = "right",
    *,
    pivoting: onp.ToTrue,
    conjugate: onp.ToBool = False,
    overwrite_a: onp.ToBool = False,
    overwrite_c: onp.ToBool = False,
) -> tuple[_FloatND, _FloatND, _IntND]: ...
@overload  # (complex[:, :], complex[:, :?], pivoting=False) -> (complex[:, :?], complex[:, :])
def qr_multiply(
    a: onp.ToComplexStrict2D,
    c: onp.ToComplexStrict1D | onp.ToComplexStrict2D,
    mode: _Side = "right",
    pivoting: onp.ToFalse = False,
    conjugate: onp.ToBool = False,
    overwrite_a: onp.ToBool = False,
    overwrite_c: onp.ToBool = False,
) -> tuple[_Inexact1D | _Inexact2D, _Inexact2D]: ...
@overload  # (complex[:, :, ...], complex[:, ...], pivoting=False) -> (complex[:, ...], complex[:, :, ...])
def qr_multiply(
    a: onp.ToComplexND,
    c: onp.ToComplexND,
    mode: _Side = "right",
    pivoting: onp.ToFalse = False,
    conjugate: onp.ToBool = False,
    overwrite_a: onp.ToBool = False,
    overwrite_c: onp.ToBool = False,
) -> tuple[_InexactND, _InexactND]: ...
@overload  # (complex[:, :], complex[:, :?], pivoting=True, /) -> (complex[:, :?], complex[:, :], int[:])
def qr_multiply(
    a: onp.ToComplexStrict2D,
    c: onp.ToComplexStrict1D | onp.ToComplexStrict2D,
    mode: _Side,
    pivoting: onp.ToTrue,
    conjugate: onp.ToBool = False,
    overwrite_a: onp.ToBool = False,
    overwrite_c: onp.ToBool = False,
) -> tuple[_Inexact1D | _Inexact2D, _Inexact2D, _Int1D]: ...
@overload  # (complex[:, :], complex[:, :?], *, pivoting=True) -> (complex[:, :?], complex[:, :], int[:])
def qr_multiply(
    a: onp.ToComplexStrict2D,
    c: onp.ToComplexStrict1D | onp.ToComplexStrict2D,
    mode: _Side = "right",
    *,
    pivoting: onp.ToTrue,
    conjugate: onp.ToBool = False,
    overwrite_a: onp.ToBool = False,
    overwrite_c: onp.ToBool = False,
) -> tuple[_Inexact1D | _Inexact2D, _Inexact2D, _Int1D]: ...
@overload  # (complex[:, :, ...], complex[:, ...], *, pivoting=True) -> (complex[:, ...], complex[:, :, ...], int[:, ...])
def qr_multiply(
    a: onp.ToComplexND,
    c: onp.ToComplexND,
    mode: _Side = "right",
    *,
    pivoting: onp.ToTrue,
    conjugate: onp.ToBool = False,
    overwrite_a: onp.ToBool = False,
    overwrite_c: onp.ToBool = False,
) -> tuple[_InexactND, _InexactND, _IntND]: ...

#
@overload  # (float[:, :], mode: {"full", "economic"}) -> (float[:, :], float[:, :])
def rq(
    a: onp.ToFloatStrict2D,
    overwrite_a: onp.ToBool = False,
    lwork: onp.ToJustInt | None = None,
    mode: _ModeFullEcon = "full",
    check_finite: onp.ToBool = True,
) -> tuple[_Float2D, _Float2D]: ...
@overload  # (float[:, :, ...], mode: {"full", "economic"}) -> (float[:, :, ...], float[:, :, ...])
def rq(
    a: onp.ToFloatND,
    overwrite_a: onp.ToBool = False,
    lwork: onp.ToJustInt | None = None,
    mode: _ModeFullEcon = "full",
    check_finite: onp.ToBool = True,
) -> tuple[_FloatND, _FloatND]: ...
@overload  # (float[:, :], mode: {"r"}, /) -> float[:, :]
def rq(
    a: onp.ToFloatStrict2D, overwrite_a: onp.ToBool, lwork: onp.ToJustInt | None, mode: _ModeR, check_finite: onp.ToBool = True
) -> _Float2D: ...
@overload  # (float[:, :, ...], mode: {"r"}, /) -> float[:, :, ...]
def rq(
    a: onp.ToFloatND, overwrite_a: onp.ToBool, lwork: onp.ToJustInt | None, mode: _ModeR, check_finite: onp.ToBool = True
) -> _FloatND: ...
@overload  # (float[:, :], *, mode: {"r"}) -> float[:, :]
def rq(
    a: onp.ToFloatStrict2D,
    overwrite_a: onp.ToBool = False,
    lwork: onp.ToJustInt | None = None,
    *,
    mode: _ModeR,
    check_finite: onp.ToBool = True,
) -> _Float2D: ...
@overload  # (float[:, :, ...], *, mode: {"r"}) -> float[:, : ...]
def rq(
    a: onp.ToFloatND,
    overwrite_a: onp.ToBool = False,
    lwork: onp.ToJustInt | None = None,
    *,
    mode: _ModeR,
    check_finite: onp.ToBool = True,
) -> _FloatND: ...
@overload  # (complex[:, :], mode: {"full", "economic"}) -> (complex[:, :], complex[:, :])
def rq(
    a: onp.ToComplexStrict2D,
    overwrite_a: onp.ToBool = False,
    lwork: onp.ToJustInt | None = None,
    mode: _ModeFullEcon = "full",
    check_finite: onp.ToBool = True,
) -> _Tuple2[_Inexact2D]: ...
@overload  # (complex[:, :], mode: {"full", "economic"}) -> (complex[:, :], complex[:, :])
def rq(
    a: onp.ToComplexND,
    overwrite_a: onp.ToBool = False,
    lwork: onp.ToJustInt | None = None,
    mode: _ModeFullEcon = "full",
    check_finite: onp.ToBool = True,
) -> _Tuple2[_InexactND]: ...
@overload  # (complex[:, :], mode: {"r"}, /) -> complex[:, :]
def rq(
    a: onp.ToComplexStrict2D, overwrite_a: onp.ToBool, lwork: onp.ToJustInt | None, mode: _ModeR, check_finite: onp.ToBool = True
) -> _Inexact2D: ...
@overload  # (complex[:, :, ...], mode: {"r"}, /) -> complex[:, :, ...]
def rq(
    a: onp.ToComplexND, overwrite_a: onp.ToBool, lwork: onp.ToJustInt | None, mode: _ModeR, check_finite: onp.ToBool = True
) -> _InexactND: ...
@overload  # (complex[:, :], *, mode: {"r"}) -> complex[:, :]
def rq(
    a: onp.ToComplexStrict2D,
    overwrite_a: onp.ToBool = False,
    lwork: onp.ToJustInt | None = None,
    *,
    mode: _ModeR,
    check_finite: onp.ToBool = True,
) -> _Inexact2D: ...
@overload  # (complex[:, :, ...], *, mode: {"r"}) -> complex[:, :, ...]
def rq(
    a: onp.ToComplexND,
    overwrite_a: onp.ToBool = False,
    lwork: onp.ToJustInt | None = None,
    *,
    mode: _ModeR,
    check_finite: onp.ToBool = True,
) -> _InexactND: ...
