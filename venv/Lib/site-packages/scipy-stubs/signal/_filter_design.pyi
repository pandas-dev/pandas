from _typeshed import Incomplete
from collections.abc import Callable, Sequence
from types import ModuleType
from typing import Any, Literal as L, SupportsIndex, TypeVar, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = [
    "BadCoefficients",
    "band_stop_obj",
    "bessel",
    "besselap",
    "bilinear",
    "bilinear_zpk",
    "buttap",
    "butter",
    "buttord",
    "cheb1ap",
    "cheb1ord",
    "cheb2ap",
    "cheb2ord",
    "cheby1",
    "cheby2",
    "ellip",
    "ellipap",
    "ellipord",
    "findfreqs",
    "freqs",
    "freqs_zpk",
    "freqz",
    "freqz_sos",
    "freqz_zpk",
    "gammatone",
    "group_delay",
    "iircomb",
    "iirdesign",
    "iirfilter",
    "iirnotch",
    "iirpeak",
    "lp2bp",
    "lp2bp_zpk",
    "lp2bs",
    "lp2bs_zpk",
    "lp2hp",
    "lp2hp_zpk",
    "lp2lp",
    "lp2lp_zpk",
    "normalize",
    "sos2tf",
    "sos2zpk",
    "sosfreqz",
    "tf2sos",
    "tf2zpk",
    "zpk2sos",
    "zpk2tf",
]

###

type _Float1D = onp.Array1D[np.float64]
type _Float2D = onp.Array2D[np.float64]
type _FloatND = onp.ArrayND[np.float64]
type _Complex1D = onp.Array1D[np.complex128]
type _ComplexND = onp.ArrayND[np.complex128]

type _Order = L[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24]

type _ZPK[ZT: np.generic, PT: np.generic, KT: np.generic | float] = tuple[onp.Array1D[ZT], onp.Array1D[PT], KT]

type _Ba1D[InexactT: npc.inexact] = tuple[onp.Array1D[InexactT], onp.Array1D[InexactT]]
type _BaND[InexactT: npc.inexact] = tuple[onp.ArrayND[InexactT], onp.Array1D[InexactT]]

# excludes `float16` and `longdouble`
type _ToFloat = float | np.float32 | np.float64 | npc.integer
type _ToFloat1D = Sequence[_ToFloat] | onp.CanArrayND[np.float32 | np.float64 | npc.integer]

type _BandStop = L["butter", "cheby", "ellip"]
type _FType0 = L["butter", "cheby1", "cheby2", "ellip"]
type _FType = _FType0 | L["bessel"]
type _BType0 = L["bandpass", "lowpass", "highpass", "bandstop"]
type _BType = _BType0 | L["band", "pass", "bp", "bands", "stop", "bs", "low", "lp", "l", "high", "hp", "h"]
type _BTypeSingle = L["l", "lp", "low", "lowpass", "h", "hp", "high", "highpass"]
type _BTypeDouble = L["bp", "band", "pass", "bandpass", "bs", "bands", "stop", "bandstop"]
type _Pairing = L["nearest", "keep_odd", "minimal"]
type _Norm = L["phase", "delay", "mag"]

type _WorNReal = int | onp.ToFloat1D | None

_AnyInexactT = TypeVar(
    "_AnyInexactT",
    np.float16,
    np.float32,
    np.float64,
    np.float96,
    np.float128,
    np.complex64,
    np.complex128,
    np.complex192,
    np.complex256,
)

###

class BadCoefficients(UserWarning): ...

#
def findfreqs(num: onp.ToComplex1D, den: onp.ToComplex1D, N: SupportsIndex, kind: L["ba", "zp"] = "ba") -> _Float1D: ...

#
@overload  # worN: real
def freqs(
    b: onp.ToComplex1D, a: onp.ToComplex1D, worN: _WorNReal = 200, plot: Callable[[_Float1D, _Complex1D], object] | None = None
) -> tuple[_Float1D, _Complex1D]: ...
@overload  # worN: complex
def freqs(
    b: onp.ToComplex1D,
    a: onp.ToComplex1D,
    worN: onp.ToJustComplex1D,
    plot: Callable[[_Float1D, _Complex1D], object] | Callable[[_Complex1D, _Complex1D], object] | None = None,
) -> tuple[_Complex1D, _Complex1D]: ...

#
@overload  # worN: real
def freqs_zpk(
    z: onp.ToComplex1D, p: onp.ToComplex1D, k: onp.ToComplex1D, worN: _WorNReal = 200
) -> tuple[_Float1D, _Complex1D]: ...
@overload  # worN: complex
def freqs_zpk(
    z: onp.ToComplex1D, p: onp.ToComplex1D, k: onp.ToComplex1D, worN: onp.ToJustComplex1D
) -> tuple[_Complex1D, _Complex1D]: ...

#
@overload  # worN: real
def freqz(
    b: onp.ToComplex | onp.ToComplexND,
    a: onp.ToComplex | onp.ToComplexND = 1,
    worN: _WorNReal = 512,
    whole: bool = False,
    plot: Callable[[_FloatND, _ComplexND], object] | None = None,
    fs: float = 6.283185307179586,
    include_nyquist: bool = False,
) -> tuple[_FloatND, _ComplexND]: ...
@overload  # worN: complex  (keyword)
def freqz(
    b: onp.ToComplex | onp.ToComplexND,
    a: onp.ToComplex | onp.ToComplexND = 1,
    *,
    worN: onp.ToJustComplex1D,
    whole: bool = False,
    plot: Callable[[_FloatND, _ComplexND], object] | None = None,
    fs: float = 6.283185307179586,
    include_nyquist: bool = False,
) -> tuple[_ComplexND, _ComplexND]: ...
@overload  # worN: complex  (positional)
def freqz(
    b: onp.ToComplex | onp.ToComplexND,
    a: onp.ToComplex | onp.ToComplexND,
    worN: onp.ToJustComplex1D,
    whole: bool = False,
    plot: Callable[[_FloatND, _ComplexND], object] | None = None,
    fs: float = 6.283185307179586,
    include_nyquist: bool = False,
) -> tuple[_ComplexND, _ComplexND]: ...

#
@overload  # worN: real
def freqz_zpk(
    z: onp.ToComplex1D,
    p: onp.ToComplex1D,
    k: onp.ToComplex1D,
    worN: _WorNReal = 512,
    whole: bool = False,
    fs: float = 6.283185307179586,
) -> tuple[_FloatND, _ComplexND]: ...
@overload  # worN: complex
def freqz_zpk(
    z: onp.ToComplex1D,
    p: onp.ToComplex1D,
    k: onp.ToComplex1D,
    worN: onp.ToJustComplex1D,
    whole: bool = False,
    fs: float = 6.283185307179586,
) -> tuple[_ComplexND, _ComplexND]: ...

#
@overload  # w: real
def group_delay(
    system: tuple[onp.ToComplex1D, onp.ToComplex1D], w: _WorNReal = 512, whole: bool = False, fs: float = 6.283185307179586
) -> _Ba1D[np.float64]: ...
@overload  # w: complex
def group_delay(
    system: tuple[onp.ToComplex1D, onp.ToComplex1D], w: onp.ToJustComplex1D, whole: bool = False, fs: float = 6.283185307179586
) -> tuple[_Complex1D, _Float1D]: ...

#
@overload  # worN: real
def freqz_sos(
    sos: onp.ToFloat2D, worN: _WorNReal = 512, whole: bool = False, fs: float = 6.283185307179586
) -> tuple[_Float1D, _Complex1D]: ...
@overload  # worN: real
def freqz_sos(
    sos: onp.ToFloat2D, worN: onp.ToJustComplex1D, whole: bool = False, fs: float = 6.283185307179586
) -> tuple[_Complex1D, _Complex1D]: ...

sosfreqz = freqz_sos

# 2/4 of mypy's false positive `overload-overlap` errors for `tf2zpk` are only reported with numpy<2.5
# mypy: disable-error-code=overload-overlap

# for real input the return dtype is either real or complex, depending on the values
@overload  # ~f64, +f64
def tf2zpk(
    b: onp.ToJustFloat64_1D | onp.ToInt1D, a: onp.ToFloat64_1D
) -> _ZPK[np.float64 | np.complex128, np.complex128, np.float64]: ...
@overload  # +f64, ~f64
def tf2zpk(
    b: onp.ToFloat64_1D, a: onp.ToJustFloat64_1D | onp.ToInt1D
) -> _ZPK[np.float64 | np.complex128, np.complex128, np.float64]: ...
@overload  # ~c128, +c128
def tf2zpk(b: onp.ToJustComplex128_1D, a: onp.ToComplex128_1D) -> _ZPK[np.complex128, np.complex128, np.float64]: ...
@overload  # +c128, ~c128
def tf2zpk(b: onp.ToComplex128_1D, a: onp.ToJustComplex128_1D) -> _ZPK[np.complex128, np.complex128, np.float64]: ...
@overload  # ~f32, +f32
def tf2zpk(b: onp.ToJustFloat32_1D, a: onp.ToFloat32_1D) -> _ZPK[np.float32 | np.complex64, np.complex64, np.float32]: ...
@overload  # +f32, ~f32
def tf2zpk(b: onp.ToFloat32_1D, a: onp.ToJustFloat32_1D) -> _ZPK[np.float32 | np.complex64, np.complex64, np.float32]: ...
@overload  # ~c64, +c64
def tf2zpk(b: onp.ToJustComplex64_1D, a: onp.ToComplex64_1D) -> _ZPK[np.complex64, np.complex64, np.float32]: ...
@overload  # +c64, ~c64
def tf2zpk(b: onp.ToComplex64_1D, a: onp.ToJustComplex64_1D) -> _ZPK[np.complex64, np.complex64, np.float32]: ...
@overload  # fallback
def tf2zpk(b: onp.ToComplex1D, a: onp.ToComplex1D) -> _ZPK[Any, Any, Any]: ...

#
def tf2sos(b: _ToFloat1D, a: _ToFloat1D, pairing: _Pairing | None = None, *, analog: bool = False) -> _Float2D: ...

#
def zpk2tf(z: onp.ToFloat1D, p: onp.ToFloat1D, k: float) -> _Ba1D[np.float64]: ...

#
def zpk2sos(
    z: onp.ToFloat1D, p: onp.ToFloat1D, k: float, pairing: _Pairing | None = None, *, analog: bool = False
) -> _Float2D: ...

#
@overload  # ~f64, +f64
def normalize(b: onp.ToIntND | onp.ToJustFloat64_ND, a: onp.ToFloat64_ND) -> _BaND[np.float64]: ...
@overload  # +f64, ~f64
def normalize(b: onp.ToFloat64_ND, a: onp.ToIntND | onp.ToJustFloat64_ND) -> _BaND[np.float64]: ...
@overload  # +c128, ~c128
def normalize(b: onp.ToComplex128_ND, a: onp.ToJustComplex128_ND) -> _BaND[np.complex128]: ...
@overload  # ~c128, +c128
def normalize(b: onp.ToJustComplex128_ND, a: onp.ToComplex128_ND) -> _BaND[np.complex128]: ...
@overload  # ~T, ~T
def normalize(b: onp.ArrayND[_AnyInexactT], a: onp.ArrayND[_AnyInexactT]) -> _BaND[_AnyInexactT]: ...  # noqa: UP047
@overload  # fallback
def normalize(b: onp.ToComplexND, a: onp.ToComplexND) -> _BaND[Any]: ...

#
@overload  # f64
def sos2tf(sos: onp.ToInt2D | onp.ToJustFloat64_2D) -> tuple[_Float1D, _Float1D]: ...
@overload  # c128
def sos2tf(sos: onp.ToJustComplex128_2D) -> tuple[_Complex1D, _Complex1D]: ...
@overload  # T: inexact
def sos2tf[InexactT: npc.inexact](sos: onp.Array2D[InexactT]) -> tuple[onp.Array1D[InexactT], onp.Array1D[InexactT]]: ...
@overload  # fallback
def sos2tf(sos: onp.ToComplex2D) -> tuple[onp.Array1D, onp.Array1D]: ...

# unlike `sos2tf`, `sos2zpk` does not accept `float16` or `[c]longdouble`
@overload  # f64
def sos2zpk(sos: onp.ToInt2D | onp.ToJustFloat64_2D) -> _ZPK[np.complex128, np.complex128, np.float64]: ...
@overload  # c128
def sos2zpk(sos: onp.ToJustComplex128_2D) -> _ZPK[np.complex128, np.complex128, np.complex128]: ...
@overload  # f32
def sos2zpk(sos: onp.ToJustFloat32_2D) -> _ZPK[np.complex128, np.complex128, np.float32]: ...
@overload  # c64
def sos2zpk(sos: onp.ToJustComplex64_2D) -> _ZPK[np.complex128, np.complex128, np.complex64]: ...
@overload  # fallback
def sos2zpk(sos: onp.ToJustComplex2D) -> _ZPK[np.complex128, np.complex128, Any]: ...

# upcasts to at least 64-bit inexact
@overload  # +f64, +f64
def lp2lp(b: onp.ToFloat64_ND, a: onp.ToFloat64_1D, wo: float = 1.0) -> tuple[_FloatND, _Float1D]: ...
@overload  # +c128, ~c128
def lp2lp(
    b: onp.ToComplex128_ND, a: onp.ToJustComplex64_1D | onp.ToJustComplex128_1D, wo: float = 1.0
) -> tuple[_ComplexND, _Complex1D]: ...
@overload  # ~c128, +c128
def lp2lp(
    b: onp.ToJustComplex64_ND | onp.ToJustComplex128_ND, a: onp.ToComplex128_1D, wo: float = 1.0
) -> tuple[_ComplexND, _Complex1D]: ...
@overload  # +f80, ~f80
def lp2lp(
    b: onp.ToFloatND, a: onp.ToJustLongDouble1D, wo: float = 1.0
) -> tuple[onp.ArrayND[np.longdouble], onp.Array1D[np.longdouble]]: ...
@overload  # ~f80, +f80
def lp2lp(
    b: onp.ToJustLongDoubleND, a: onp.ToFloat1D, wo: float = 1.0
) -> tuple[onp.ArrayND[np.longdouble], onp.Array1D[np.longdouble]]: ...
@overload  # +c160, ~c160
def lp2lp(
    b: onp.ToComplexND, a: onp.ToJustCLongDouble1D, wo: float = 1.0
) -> tuple[onp.ArrayND[np.clongdouble], onp.Array1D[np.clongdouble]]: ...
@overload  # ~c160, +c160
def lp2lp(
    b: onp.ToJustCLongDoubleND, a: onp.ToComplex1D, wo: float = 1.0
) -> tuple[onp.ArrayND[np.clongdouble], onp.Array1D[np.clongdouble]]: ...
@overload  # fallback
def lp2lp(b: onp.ToComplexND, a: onp.ToComplex1D, wo: float = 1.0) -> tuple[onp.ArrayND[Any], onp.Array1D[Any]]: ...

# lp2hp
@overload  # +f64, +f64
def lp2hp(b: onp.ToFloat64_1D, a: onp.ToFloat64_1D, wo: float = 1.0) -> tuple[_Float1D, _Float1D]: ...
@overload  # +c128, ~c128
def lp2hp(
    b: onp.ToComplex128_1D, a: onp.ToJustComplex64_1D | onp.ToJustComplex128_1D, wo: float = 1.0
) -> tuple[_Complex1D, _Complex1D]: ...
@overload  # ~c128, +c128
def lp2hp(
    b: onp.ToJustComplex64_1D | onp.ToJustComplex128_1D, a: onp.ToComplex128_1D, wo: float = 1.0
) -> tuple[_Complex1D, _Complex1D]: ...
@overload  # +f80, ~f80
def lp2hp(
    b: onp.ToFloat1D, a: onp.ToJustLongDouble1D, wo: float = 1.0
) -> tuple[onp.Array1D[np.longdouble], onp.Array1D[np.longdouble]]: ...
@overload  # ~f80, +f80
def lp2hp(
    b: onp.ToJustLongDouble1D, a: onp.ToFloat1D, wo: float = 1.0
) -> tuple[onp.Array1D[np.longdouble], onp.Array1D[np.longdouble]]: ...
@overload  # +c160, ~c160
def lp2hp(
    b: onp.ToComplex1D, a: onp.ToJustCLongDouble1D, wo: float = 1.0
) -> tuple[onp.Array1D[np.clongdouble], onp.Array1D[np.clongdouble]]: ...
@overload  # ~c160, +c160
def lp2hp(
    b: onp.ToJustCLongDouble1D, a: onp.ToComplex1D, wo: float = 1.0
) -> tuple[onp.Array1D[np.clongdouble], onp.Array1D[np.clongdouble]]: ...
@overload  # fallback
def lp2hp(b: onp.ToComplex1D, a: onp.ToComplex1D, wo: float = 1.0) -> tuple[onp.Array1D[Any], onp.Array1D[Any]]: ...

# keep in sync with `lp2hp` (above)
@overload  # +f64, +f64
def lp2bp(b: onp.ToFloat64_1D, a: onp.ToFloat64_1D, wo: float = 1.0, bw: float = 1.0) -> tuple[_Float1D, _Float1D]: ...
@overload  # +c128, ~c128
def lp2bp(
    b: onp.ToComplex128_1D, a: onp.ToJustComplex64_1D | onp.ToJustComplex128_1D, wo: float = 1.0, bw: float = 1.0
) -> tuple[_Complex1D, _Complex1D]: ...
@overload  # ~c128, +c128
def lp2bp(
    b: onp.ToJustComplex64_1D | onp.ToJustComplex128_1D, a: onp.ToComplex128_1D, wo: float = 1.0, bw: float = 1.0
) -> tuple[_Complex1D, _Complex1D]: ...
@overload  # +f80, ~f80
def lp2bp(
    b: onp.ToFloat1D, a: onp.ToJustLongDouble1D, wo: float = 1.0, bw: float = 1.0
) -> tuple[onp.Array1D[np.longdouble], onp.Array1D[np.longdouble]]: ...
@overload  # ~f80, +f80
def lp2bp(
    b: onp.ToJustLongDouble1D, a: onp.ToFloat1D, wo: float = 1.0, bw: float = 1.0
) -> tuple[onp.Array1D[np.longdouble], onp.Array1D[np.longdouble]]: ...
@overload  # +c160, ~c160
def lp2bp(
    b: onp.ToComplex1D, a: onp.ToJustCLongDouble1D, wo: float = 1.0, bw: float = 1.0
) -> tuple[onp.Array1D[np.clongdouble], onp.Array1D[np.clongdouble]]: ...
@overload  # ~c160, +c160
def lp2bp(
    b: onp.ToJustCLongDouble1D, a: onp.ToComplex1D, wo: float = 1.0, bw: float = 1.0
) -> tuple[onp.Array1D[np.clongdouble], onp.Array1D[np.clongdouble]]: ...
@overload  # fallback
def lp2bp(
    b: onp.ToComplex1D, a: onp.ToComplex1D, wo: float = 1.0, bw: float = 1.0
) -> tuple[onp.Array1D[Any], onp.Array1D[Any]]: ...

# keep in sync with `lp2bp` (above)
@overload  # +f64, +f64
def lp2bs(b: onp.ToFloat64_1D, a: onp.ToFloat64_1D, wo: float = 1.0, bw: float = 1.0) -> tuple[_Float1D, _Float1D]: ...
@overload  # +c128, ~c128
def lp2bs(
    b: onp.ToComplex128_1D, a: onp.ToJustComplex64_1D | onp.ToJustComplex128_1D, wo: float = 1.0, bw: float = 1.0
) -> tuple[_Complex1D, _Complex1D]: ...
@overload  # ~c128, +c128
def lp2bs(
    b: onp.ToJustComplex64_1D | onp.ToJustComplex128_1D, a: onp.ToComplex128_1D, wo: float = 1.0, bw: float = 1.0
) -> tuple[_Complex1D, _Complex1D]: ...
@overload  # +f80, ~f80
def lp2bs(
    b: onp.ToFloat1D, a: onp.ToJustLongDouble1D, wo: float = 1.0, bw: float = 1.0
) -> tuple[onp.Array1D[np.longdouble], onp.Array1D[np.longdouble]]: ...
@overload  # ~f80, +f80
def lp2bs(
    b: onp.ToJustLongDouble1D, a: onp.ToFloat1D, wo: float = 1.0, bw: float = 1.0
) -> tuple[onp.Array1D[np.longdouble], onp.Array1D[np.longdouble]]: ...
@overload  # +c160, ~c160
def lp2bs(
    b: onp.ToComplex1D, a: onp.ToJustCLongDouble1D, wo: float = 1.0, bw: float = 1.0
) -> tuple[onp.Array1D[np.clongdouble], onp.Array1D[np.clongdouble]]: ...
@overload  # ~c160, +c160
def lp2bs(
    b: onp.ToJustCLongDouble1D, a: onp.ToComplex1D, wo: float = 1.0, bw: float = 1.0
) -> tuple[onp.Array1D[np.clongdouble], onp.Array1D[np.clongdouble]]: ...
@overload  # fallback
def lp2bs(
    b: onp.ToComplex1D, a: onp.ToComplex1D, wo: float = 1.0, bw: float = 1.0
) -> tuple[onp.Array1D[Any], onp.Array1D[Any]]: ...

# lp2lp_zpk
@overload
def lp2lp_zpk(
    z: onp.ToInt1D | onp.ToJustFloat64_1D, p: onp.ToInt1D | onp.ToJustFloat64_1D, k: onp.ToFloat64, wo: float = 1.0
) -> _ZPK[np.float64, np.float64, float]: ...
@overload
def lp2lp_zpk(
    z: onp.ToJustComplex128_1D, p: onp.ToInt1D | onp.ToJustFloat64_1D, k: onp.ToFloat64, wo: float = 1.0
) -> _ZPK[np.complex128, np.float64, float]: ...
@overload
def lp2lp_zpk(
    z: onp.ToInt1D | onp.ToJustFloat64_1D, p: onp.ToJustComplex128_1D, k: onp.ToFloat64, wo: float = 1.0
) -> _ZPK[np.float64, np.complex128, float]: ...
@overload
def lp2lp_zpk(
    z: onp.ToJustComplex128_1D, p: onp.ToJustComplex128_1D, k: onp.ToFloat64, wo: float = 1.0
) -> _ZPK[np.complex128, np.complex128, float]: ...
@overload
def lp2lp_zpk(
    z: onp.ToJustFloat32_1D, p: onp.ToJustFloat32_1D, k: onp.ToFloat32, wo: float = 1.0
) -> _ZPK[np.float32, np.float32, float]: ...
@overload
def lp2lp_zpk(
    z: onp.ToJustComplex64_1D, p: onp.ToJustFloat32_1D, k: onp.ToFloat32, wo: float = 1.0
) -> _ZPK[np.complex64, np.float32, float]: ...
@overload
def lp2lp_zpk(
    z: onp.ToJustFloat32_1D, p: onp.ToJustComplex64_1D, k: onp.ToFloat32, wo: float = 1.0
) -> _ZPK[np.float32, np.complex64, float]: ...
@overload
def lp2lp_zpk(
    z: onp.ToJustComplex64_1D, p: onp.ToJustComplex64_1D, k: onp.ToFloat32, wo: float = 1.0
) -> _ZPK[np.complex64, np.complex64, float]: ...
@overload
def lp2lp_zpk(z: onp.ToFloat1D, p: onp.ToFloat1D, k: onp.ToFloat, wo: float = 1.0) -> _ZPK[npc.floating, npc.floating, float]: ...
@overload
def lp2lp_zpk(
    z: onp.ToComplex1D, p: onp.ToComplex1D, k: onp.ToFloat, wo: float = 1.0
) -> _ZPK[npc.inexact, npc.inexact, float]: ...

#
@overload
def lp2hp_zpk(
    z: onp.ToInt1D | onp.ToJustFloat64_1D, p: onp.ToInt1D | onp.ToJustFloat64_1D, k: onp.ToFloat64, wo: float = 1.0
) -> _ZPK[np.float64, np.float64, np.float64]: ...
@overload
def lp2hp_zpk(
    z: onp.ToJustComplex128_1D, p: onp.ToInt1D | onp.ToJustFloat64_1D, k: onp.ToFloat64, wo: float = 1.0
) -> _ZPK[np.complex128, np.float64, np.float64]: ...
@overload
def lp2hp_zpk(
    z: onp.ToInt1D | onp.ToJustFloat64_1D, p: onp.ToJustComplex128_1D, k: onp.ToFloat64, wo: float = 1.0
) -> _ZPK[np.float64, np.complex128, np.float64]: ...
@overload
def lp2hp_zpk(
    z: onp.ToJustComplex128_1D, p: onp.ToJustComplex128_1D, k: onp.ToFloat64, wo: float = 1.0
) -> _ZPK[np.complex128, np.complex128, np.float64]: ...
@overload
def lp2hp_zpk(
    z: onp.ToJustFloat32_1D, p: onp.ToJustFloat32_1D, k: onp.ToFloat32, wo: float = 1.0
) -> _ZPK[np.float64, np.float64, np.float64]: ...
@overload
def lp2hp_zpk(
    z: onp.ToJustComplex64_1D, p: onp.ToJustFloat32_1D, k: onp.ToFloat32, wo: float = 1.0
) -> _ZPK[np.complex128, np.float64, np.float64]: ...
@overload
def lp2hp_zpk(
    z: onp.ToJustFloat32_1D, p: onp.ToJustComplex64_1D, k: onp.ToFloat32, wo: float = 1.0
) -> _ZPK[np.float64, np.complex128, np.float64]: ...
@overload
def lp2hp_zpk(
    z: onp.ToJustComplex64_1D, p: onp.ToJustComplex64_1D, k: onp.ToFloat32, wo: float = 1.0
) -> _ZPK[np.complex128, np.complex128, np.float64]: ...
@overload
def lp2hp_zpk(
    z: onp.ToComplex1D, p: onp.ToComplex1D, k: onp.ToFloat, wo: float = 1.0
) -> _ZPK[npc.inexact, npc.inexact, npc.floating]: ...

# lp2bp_zpk
@overload
def lp2bp_zpk(
    z: onp.ToInt1D | onp.ToJustFloat64_1D | onp.ToJustComplex128_1D,
    p: onp.ToInt1D | onp.ToJustFloat64_1D | onp.ToJustComplex128_1D,
    k: onp.ToFloat64,
    wo: float = 1.0,
    bw: float = 1.0,
) -> _ZPK[np.complex128, np.complex128, float]: ...
@overload
def lp2bp_zpk(
    z: onp.ToJustFloat32_1D | onp.ToJustComplex64_1D,
    p: onp.ToJustFloat32_1D | onp.ToJustComplex64_1D,
    k: onp.ToFloat32,
    wo: float = 1.0,
    bw: float = 1.0,
) -> _ZPK[np.complex64, np.complex64, float]: ...
@overload
def lp2bp_zpk(
    z: onp.ToComplex1D, p: onp.ToComplex1D, k: onp.ToFloat, wo: float = 1.0, bw: float = 1.0
) -> _ZPK[npc.complexfloating, npc.complexfloating, float]: ...

# lp2bs_zpk
@overload
def lp2bs_zpk(
    z: onp.ToInt1D | onp.ToJustFloat64_1D | onp.ToJustComplex128_1D,
    p: onp.ToInt1D | onp.ToJustFloat64_1D | onp.ToJustComplex128_1D,
    k: onp.ToFloat64,
    wo: float = 1.0,
    bw: float = 1.0,
) -> _ZPK[np.complex128, np.complex128, np.float64]: ...
@overload
def lp2bs_zpk(
    z: onp.ToJustFloat32_1D | onp.ToJustComplex64_1D,
    p: onp.ToJustFloat32_1D | onp.ToJustComplex64_1D,
    k: onp.ToFloat32,
    wo: float = 1.0,
    bw: float = 1.0,
) -> _ZPK[np.complex64, np.complex64, np.float32]: ...
@overload
def lp2bs_zpk(
    z: onp.ToComplex1D, p: onp.ToComplex1D, k: onp.ToFloat, wo: float = 1.0, bw: float = 1.0
) -> _ZPK[npc.complexfloating, npc.complexfloating, npc.floating]: ...

#
def bilinear(b: onp.ToFloat1D, a: onp.ToFloat1D, fs: float = 1.0) -> _Ba1D[np.float64]: ...
@overload
def bilinear_zpk(
    z: onp.ToInt1D | onp.ToJustFloat64_1D, p: onp.ToInt1D | onp.ToJustFloat64_1D, k: onp.ToFloat64, fs: float
) -> _ZPK[np.float64, np.float64, np.float64]: ...
@overload
def bilinear_zpk(
    z: onp.ToJustComplex128_1D, p: onp.ToInt1D | onp.ToJustFloat64_1D, k: onp.ToFloat64, fs: float
) -> _ZPK[np.complex128, np.float64, np.float64]: ...
@overload
def bilinear_zpk(
    z: onp.ToInt1D | onp.ToJustFloat64_1D, p: onp.ToJustComplex128_1D, k: onp.ToFloat64, fs: float
) -> _ZPK[np.float64, np.complex128, np.float64]: ...
@overload
def bilinear_zpk(
    z: onp.ToJustComplex128_1D, p: onp.ToJustComplex128_1D, k: onp.ToFloat64, fs: float
) -> _ZPK[np.complex128, np.complex128, np.float64]: ...
@overload
def bilinear_zpk(
    z: onp.ToJustFloat32_1D, p: onp.ToJustFloat32_1D, k: onp.ToFloat32, fs: float
) -> _ZPK[np.float64, np.float64, np.float64]: ...
@overload
def bilinear_zpk(
    z: onp.ToJustComplex64_1D, p: onp.ToJustFloat32_1D, k: onp.ToFloat32, fs: float
) -> _ZPK[np.complex128, np.float64, np.float64]: ...
@overload
def bilinear_zpk(
    z: onp.ToJustFloat32_1D, p: onp.ToJustComplex64_1D, k: onp.ToFloat32, fs: float
) -> _ZPK[np.float64, np.complex128, np.float64]: ...
@overload
def bilinear_zpk(
    z: onp.ToJustComplex64_1D, p: onp.ToJustComplex64_1D, k: onp.ToFloat32, fs: float
) -> _ZPK[np.complex128, np.complex128, np.float64]: ...
@overload
def bilinear_zpk(
    z: onp.ToComplex1D, p: onp.ToComplex1D, k: onp.ToFloat, fs: float
) -> _ZPK[npc.inexact, npc.inexact, npc.floating]: ...

#
@overload  # output="ba" (default)
def iirdesign(
    wp: float | onp.ToFloat1D,
    ws: float | onp.ToFloat1D,
    gpass: float,
    gstop: float,
    analog: bool = False,
    ftype: _FType0 = "ellip",
    output: L["ba"] = "ba",
    fs: float | None = None,
) -> _Ba1D[np.float64]: ...
@overload  # ftype={"cheby2", "ellip"} (default), output="zpk"
def iirdesign(
    wp: float | onp.ToFloat1D,
    ws: float | onp.ToFloat1D,
    gpass: float,
    gstop: float,
    analog: bool = False,
    ftype: L["cheby2", "ellip"] = "ellip",
    *,
    output: L["zpk"],
    fs: float | None = None,
) -> _ZPK[np.complex128, np.complex128, np.float64]: ...
@overload  # ftype={"butter", "cheby1"}, lowpass/highpass (scalar wp/ws), output="zpk"
def iirdesign(
    wp: float,
    ws: float,
    gpass: float,
    gstop: float,
    analog: bool = False,
    *,
    ftype: L["butter", "cheby1"],
    output: L["zpk"],
    fs: float | None = None,
) -> _ZPK[np.float64, np.complex128, np.float64]: ...
@overload  # ftype={"butter", "cheby1"}, bandpass/bandstop (array wp/ws), output="zpk"
def iirdesign(
    wp: onp.ToFloat1D,
    ws: onp.ToFloat1D,
    gpass: float,
    gstop: float,
    analog: bool = False,
    *,
    ftype: L["butter", "cheby1"],
    output: L["zpk"],
    fs: float | None = None,
) -> _ZPK[np.complex128, np.complex128, np.float64]: ...
@overload  # output="sos"
def iirdesign(
    wp: float | onp.ToFloat1D,
    ws: float | onp.ToFloat1D,
    gpass: float,
    gstop: float,
    analog: bool = False,
    ftype: _FType0 = "ellip",
    *,
    output: L["sos"],
    fs: float | None = None,
) -> _Float2D: ...

#
@overload  # output="ba" (default)
def iirfilter(
    N: int,
    Wn: float | onp.ToFloat1D,
    rp: float | None = None,
    rs: float | None = None,
    btype: _BType = "band",
    analog: bool = False,
    ftype: _FType = "butter",
    output: L["ba"] = "ba",
    fs: float | None = None,
) -> _Ba1D[np.float64]: ...
@overload  # btype={"bandpass", "bandstop"} (default), output="zpk"
def iirfilter(
    N: int,
    Wn: onp.ToFloat1D,
    rp: float | None = None,
    rs: float | None = None,
    btype: _BTypeDouble = "band",
    analog: bool = False,
    ftype: _FType = "butter",
    *,
    output: L["zpk"],
    fs: float | None = None,
) -> _ZPK[np.complex128, np.complex128, np.float64]: ...
@overload  # btype={"lowpass", "highpass"}, ftype={"butter", "cheby1", "bessel"} (default), output="zpk"
def iirfilter(
    N: int,
    Wn: float,
    rp: float | None = None,
    rs: float | None = None,
    *,
    btype: _BTypeSingle,
    analog: bool = False,
    ftype: L["butter", "cheby1", "bessel"] = "butter",
    output: L["zpk"],
    fs: float | None = None,
) -> _ZPK[np.float64, np.complex128, np.float64]: ...
@overload  # btype={"lowpass", "highpass"}, ftype={"cheby2", "ellip"}, output="zpk"
def iirfilter(
    N: int,
    Wn: float,
    rp: float | None = None,
    rs: float | None = None,
    *,
    btype: _BTypeSingle,
    analog: bool = False,
    ftype: L["cheby2", "ellip"],
    output: L["zpk"],
    fs: float | None = None,
) -> _ZPK[np.complex128, np.complex128, np.float64]: ...
@overload  # output="sos"
def iirfilter(
    N: int,
    Wn: float | onp.ToFloat1D,
    rp: float | None = None,
    rs: float | None = None,
    btype: _BType = "band",
    analog: bool = False,
    ftype: _FType = "butter",
    *,
    output: L["sos"],
    fs: float | None = None,
) -> _Float2D: ...

#
@overload  # output="ba" (default)
def butter(
    N: int,
    Wn: float | onp.ToFloat1D,
    btype: _BType = "low",
    analog: bool = False,
    output: L["ba"] = "ba",
    fs: float | None = None,
) -> _Ba1D[np.float64]: ...
@overload  # btype={"lowpass", "highpass"} (default), output="zpk"
def butter(
    N: int, Wn: float, btype: _BTypeSingle = "low", analog: bool = False, *, output: L["zpk"], fs: float | None = None
) -> _ZPK[np.float64, np.complex128, np.float64]: ...
@overload  # btype={"bandpass", "bandstop"}, output="zpk"
def butter(
    N: int, Wn: onp.ToFloat1D, btype: _BTypeDouble, analog: bool = False, *, output: L["zpk"], fs: float | None = None
) -> _ZPK[np.complex128, np.complex128, np.float64]: ...
@overload  # output="sos"
def butter(
    N: int, Wn: float | onp.ToFloat1D, btype: _BType = "low", analog: bool = False, *, output: L["sos"], fs: float | None = None
) -> _Float2D: ...

#
@overload  # output="ba" (default)
def cheby1(
    N: int,
    rp: float,
    Wn: float | onp.ToFloat1D,
    btype: _BType = "low",
    analog: bool = False,
    output: L["ba"] = "ba",
    fs: float | None = None,
) -> _Ba1D[np.float64]: ...
@overload  # btype={"lowpass", "highpass"} (default), output="zpk"
def cheby1(
    N: int, rp: float, Wn: float, btype: _BTypeSingle = "low", analog: bool = False, *, output: L["zpk"], fs: float | None = None
) -> _ZPK[np.float64, np.complex128, np.float64]: ...
@overload  # btype={"bandpass", "bandstop"}, output="zpk"
def cheby1(
    N: int, rp: float, Wn: onp.ToFloat1D, btype: _BTypeDouble, analog: bool = False, *, output: L["zpk"], fs: float | None = None
) -> _ZPK[np.complex128, np.complex128, np.float64]: ...
@overload  # output="sos"
def cheby1(
    N: int,
    rp: float,
    Wn: float | onp.ToFloat1D,
    btype: _BType = "low",
    analog: bool = False,
    *,
    output: L["sos"],
    fs: float | None = None,
) -> _Float2D: ...

#
@overload  # output="ba" (default)
def cheby2(
    N: int,
    rs: float,
    Wn: float | onp.ToFloat1D,
    btype: _BType = "low",
    analog: bool = False,
    output: L["ba"] = "ba",
    fs: float | None = None,
) -> _Ba1D[np.float64]: ...
@overload  # output="zpk"
def cheby2(
    N: int,
    rs: float,
    Wn: float | onp.ToFloat1D,
    btype: _BType = "low",
    analog: bool = False,
    *,
    output: L["zpk"],
    fs: float | None = None,
) -> _ZPK[np.complex128, np.complex128, np.float64]: ...
@overload  # output="sos"
def cheby2(
    N: int,
    rs: float,
    Wn: float | onp.ToFloat1D,
    btype: _BType = "low",
    analog: bool = False,
    *,
    output: L["sos"],
    fs: float | None = None,
) -> _Float2D: ...

#
@overload  # output="ba" (default)
def ellip(
    N: int,
    rp: float,
    rs: float,
    Wn: float | onp.ToFloat1D,
    btype: _BType = "low",
    analog: bool = False,
    output: L["ba"] = "ba",
    fs: float | None = None,
) -> _Ba1D[np.float64]: ...
@overload  # output="zpk"
def ellip(
    N: int,
    rp: float,
    rs: float,
    Wn: float | onp.ToFloat1D,
    btype: _BType = "low",
    analog: bool = False,
    *,
    output: L["zpk"],
    fs: float | None = None,
) -> _ZPK[np.complex128, np.complex128, np.float64]: ...
@overload  # output="sos"
def ellip(
    N: int,
    rp: float,
    rs: float,
    Wn: float | onp.ToFloat1D,
    btype: _BType = "low",
    analog: bool = False,
    *,
    output: L["sos"],
    fs: float | None = None,
) -> _Float2D: ...

#
@overload  # output="ba" (default)
def bessel(
    N: int,
    Wn: float | onp.ToFloat1D,
    btype: _BType = "low",
    analog: bool = False,
    output: L["ba"] = "ba",
    norm: _Norm = "phase",
    fs: float | None = None,
) -> _Ba1D[np.float64]: ...
@overload  # btype={"lowpass", "highpass"} (default), output="zpk"
def bessel(
    N: int,
    Wn: float,
    btype: _BTypeSingle = "low",
    analog: bool = False,
    *,
    output: L["zpk"],
    norm: _Norm = "phase",
    fs: float | None = None,
) -> _ZPK[np.float64, np.complex128, np.float64]: ...
@overload  # btype={"bandpass", "bandstop"}, output="zpk"
def bessel(
    N: int,
    Wn: onp.ToFloat1D,
    btype: _BTypeDouble,
    analog: bool = False,
    *,
    output: L["zpk"],
    norm: _Norm = "phase",
    fs: float | None = None,
) -> _ZPK[np.complex128, np.complex128, np.float64]: ...
@overload  # output="sos"
def bessel(
    N: int,
    Wn: float | onp.ToFloat1D,
    btype: _BType = "low",
    analog: bool = False,
    *,
    output: L["sos"],
    norm: _Norm = "phase",
    fs: float | None = None,
) -> _Float2D: ...

#
@overload  # +float64, +float64
def band_stop_obj(
    wp: float,
    ind: L[0, 1] | npc.integer,
    passb: onp.ArrayND[np.float64 | np.float32 | np.float16 | npc.integer],
    stopb: onp.ArrayND[np.float64 | np.float32 | np.float16 | npc.integer],
    gpass: float,
    gstop: float,
    type: _BandStop,
) -> np.float64: ...
@overload  # +longdouble, +longdouble  (we can't have specific longdouble overloads due to numpy <2.2 compatibility)
def band_stop_obj(
    wp: float,
    ind: L[0, 1] | npc.integer,
    passb: onp.ArrayND[npc.floating],
    stopb: onp.ArrayND[npc.floating],
    gpass: float,
    gstop: float,
    type: _BandStop,
) -> np.longdouble | Any: ...

#
@overload
def buttord(
    wp: float, ws: float | onp.ToFloat64_ND, gpass: float, gstop: float, analog: bool = False, fs: float | None = None
) -> tuple[int, np.float64]: ...
@overload
def buttord(
    wp: float, ws: onp.ToJustLongDoubleND, gpass: float, gstop: float, analog: bool = False, fs: float | None = None
) -> tuple[int, np.longdouble]: ...
@overload
def buttord(  # N-d longdouble gets downcast to float64 for some reason
    wp: onp.ToFloatND, ws: float | onp.ToFloatND, gpass: float, gstop: float, analog: bool = False, fs: float | None = None
) -> tuple[int, onp.Array1D[np.float64]]: ...

#
@overload
def cheb1ord(
    wp: float, ws: float | onp.ToFloat64_1D, gpass: float, gstop: float, analog: bool = False, fs: float | None = None
) -> tuple[int, np.float64]: ...
@overload
def cheb1ord(
    wp: float, ws: onp.ToJustLongDouble1D, gpass: float, gstop: float, analog: bool = False, fs: float | None = None
) -> tuple[int, np.longdouble]: ...
@overload  # N-d longdouble gets downcast to float64 for some reason
def cheb1ord(
    wp: onp.ToFloatND, ws: float | onp.ToFloat1D, gpass: float, gstop: float, analog: bool = False, fs: float | None = None
) -> tuple[int, onp.Array1D[np.float64]]: ...

#
@overload
def cheb2ord(
    wp: float, ws: float | onp.ToFloat64_1D, gpass: float, gstop: float, analog: bool = False, fs: float | None = None
) -> tuple[int, np.float64]: ...
@overload
def cheb2ord(
    wp: float, ws: onp.ToJustLongDouble1D, gpass: float, gstop: float, analog: bool = False, fs: float | None = None
) -> tuple[int, np.longdouble]: ...
@overload  # N-d longdouble gets downcast to float64 for some reason
def cheb2ord(
    wp: onp.ToFloatND, ws: float | onp.ToFloat1D, gpass: float, gstop: float, analog: bool = False, fs: float | None = None
) -> tuple[int, onp.Array1D[np.float64]]: ...

# unlike the order `*ord` functions, `ellipord` does not support `longdouble` input
@overload
def ellipord(
    wp: float, ws: float | onp.ToFloat64_1D, gpass: float, gstop: float, analog: bool = False, fs: float | None = None
) -> tuple[int, np.float64]: ...
@overload
def ellipord(
    wp: onp.ToFloatND, ws: float | onp.ToFloat64_1D, gpass: float, gstop: float, analog: bool = False, fs: float | None = None
) -> tuple[int, onp.Array1D[np.float64]]: ...

#
@overload
def buttap(N: int, *, xp: None = None, device: None = None) -> tuple[_Float1D, _Complex1D, L[1]]: ...
@overload
def buttap(N: int, *, xp: ModuleType, device: object = None) -> tuple[Incomplete, Incomplete, L[1]]: ...

#
@overload
def cheb1ap(N: int, rp: float, *, xp: None = None, device: None = None) -> tuple[_Float1D, _Complex1D, float]: ...
@overload
def cheb1ap(N: int, rp: float, *, xp: ModuleType, device: object = None) -> tuple[Incomplete, Incomplete, float]: ...

#
@overload
def cheb2ap(N: int, rs: float, *, xp: None = None, device: None = None) -> tuple[_Complex1D, _Complex1D, float]: ...
@overload
def cheb2ap(N: int, rs: float, *, xp: ModuleType, device: object = None) -> tuple[Incomplete, Incomplete, float]: ...

#
@overload
def ellipap(N: int, rp: float, rs: float, *, xp: None = None, device: None = None) -> tuple[_Complex1D, _Complex1D, float]: ...
@overload
def ellipap(N: int, rp: float, rs: float, *, xp: ModuleType, device: object = None) -> tuple[Incomplete, Incomplete, float]: ...

#
@overload
def besselap(N: int, norm: _Norm = "phase", *, xp: None = None, device: None = None) -> tuple[_Float1D, _Complex1D, float]: ...
@overload
def besselap(N: int, norm: _Norm = "phase", *, xp: ModuleType, device: object = None) -> tuple[Incomplete, Incomplete, float]: ...

#
@overload
def iirnotch(w0: float, Q: float, fs: float = 2.0, *, xp: None = None, device: None = None) -> _Ba1D[np.float64]: ...
@overload
def iirnotch(w0: float, Q: float, fs: float = 2.0, *, xp: ModuleType, device: object = None) -> tuple[Incomplete, Incomplete]: ...

#
@overload
def iirpeak(w0: float, Q: float, fs: float = 2.0, *, xp: None = None, device: None = None) -> _Ba1D[np.float64]: ...
@overload
def iirpeak(w0: float, Q: float, fs: float = 2.0, *, xp: ModuleType, device: object = None) -> tuple[Incomplete, Incomplete]: ...

#
@overload
def iircomb(
    w0: float,
    Q: float,
    ftype: L["notch", "peak"] = "notch",
    fs: float = 2.0,
    *,
    pass_zero: bool = False,
    xp: None = None,
    device: None = None,
) -> _Ba1D[np.float64]: ...
@overload
def iircomb(
    w0: float,
    Q: float,
    ftype: L["notch", "peak"] = "notch",
    fs: float = 2.0,
    *,
    pass_zero: bool = False,
    xp: ModuleType,
    device: object = None,
) -> tuple[Incomplete, Incomplete]: ...

#
@overload
def gammatone(
    freq: float,
    ftype: L["fir", "iir"],
    order: _Order | None = None,
    numtaps: int | None = None,
    fs: float | None = None,
    *,
    xp: None = None,
    device: None = None,
) -> _Ba1D[np.float64]: ...
@overload
def gammatone(
    freq: float,
    ftype: L["fir", "iir"],
    order: _Order | None = None,
    numtaps: int | None = None,
    fs: float | None = None,
    *,
    xp: ModuleType,
    device: object = None,
) -> tuple[Incomplete, Incomplete]: ...

# ???
def maxflat() -> None: ...
def yulewalk() -> None: ...
