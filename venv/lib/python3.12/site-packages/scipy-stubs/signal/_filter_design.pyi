from collections.abc import Callable, Sequence
from typing import Literal as L, TypeAlias, overload
from typing_extensions import TypeVar

import numpy as np
import optype as op
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

_Floating: TypeAlias = npc.floating
_CFloating: TypeAlias = npc.complexfloating

_Floating1D: TypeAlias = onp.Array1D[npc.floating]
_FloatingND: TypeAlias = onp.ArrayND[npc.floating]
_Float1D: TypeAlias = onp.Array1D[np.float64]
_Float2D: TypeAlias = onp.Array2D[np.float64]
_FloatND: TypeAlias = onp.ArrayND[np.float64]
_Complex1D: TypeAlias = onp.Array1D[np.complex128]
_ComplexND: TypeAlias = onp.ArrayND[np.complex128]
_Inexact1D: TypeAlias = onp.Array1D[np.float64 | np.complex128]
_InexactND: TypeAlias = onp.ArrayND[np.float64 | np.complex128]

_SCT_z = TypeVar("_SCT_z", bound=np.generic)
_SCT_p = TypeVar("_SCT_p", bound=np.generic, default=np.complex128)
_SCT_k = TypeVar("_SCT_k", bound=np.generic | float, default=np.float64)
_ZPK: TypeAlias = tuple[onp.Array1D[_SCT_z], onp.Array1D[_SCT_p], _SCT_k]

_SCT_ba = TypeVar("_SCT_ba", bound=npc.floating, default=np.float64)
_Ba1D: TypeAlias = tuple[onp.Array1D[_SCT_ba], onp.Array1D[_SCT_ba]]
_Ba2D: TypeAlias = tuple[onp.Array2D[_SCT_ba], onp.Array1D[_SCT_ba]]

# excludes `float16` and `longdouble`
_ToFloat: TypeAlias = float | np.float32 | np.float64 | npc.integer
_ToFloat1D: TypeAlias = Sequence[_ToFloat] | onp.CanArrayND[np.float32 | np.float64 | npc.integer]
_ToFloat2D: TypeAlias = Sequence[_ToFloat1D] | onp.CanArrayND[np.float32 | np.float64 | npc.integer]

_FType0: TypeAlias = L["butter", "cheby1", "cheby2", "ellip"]
_FType: TypeAlias = _FType0 | L["bessel"]
_BType0: TypeAlias = L["bandpass", "lowpass", "highpass", "bandstop"]
_BType: TypeAlias = _BType0 | L["band", "pass", "bp", "bands", "stop", "bs", "low", "lp", "l", "high", "hp", "h"]
_Pairing: TypeAlias = L["nearest", "keep_odd", "minimal"]
_Normalization: TypeAlias = L["phase", "delay", "mag"]

###

class BadCoefficients(UserWarning): ...

#
def findfreqs(num: onp.ToComplex1D, den: onp.ToComplex1D, N: op.CanIndex, kind: L["ba", "zp"] = "ba") -> _Float1D: ...

#
@overload  # worN: real
def freqs(
    b: onp.ToComplex1D,
    a: onp.ToComplex1D,
    worN: op.CanIndex | onp.ToFloat1D = 200,
    plot: Callable[[_Float1D, _Complex1D], object] | None = None,
) -> tuple[_Float1D, _Complex1D]: ...
@overload  # worN: complex
def freqs(
    b: onp.ToComplex1D,
    a: onp.ToComplex1D,
    worN: onp.ToComplex1D,
    plot: Callable[[_Float1D, _Complex1D], object] | Callable[[_Complex1D, _Complex1D], object] | None = None,
) -> tuple[_Inexact1D, _Complex1D]: ...

#
@overload  # worN: real
def freqs_zpk(
    z: onp.ToComplex1D, p: onp.ToComplex1D, k: onp.ToComplex1D, worN: op.CanIndex | onp.ToFloat1D = 200
) -> tuple[_Float1D, _Complex1D]: ...
@overload  # worN: complex
def freqs_zpk(
    z: onp.ToComplex1D, p: onp.ToComplex1D, k: onp.ToComplex1D, worN: onp.ToComplex1D
) -> tuple[_Inexact1D, _Complex1D]: ...

#
@overload  # worN: real
def freqz(
    b: onp.ToComplex | onp.ToComplexND,
    a: onp.ToComplex | onp.ToComplexND = 1,
    worN: op.CanIndex | onp.ToFloat1D = 512,
    whole: op.CanBool = False,
    plot: Callable[[_FloatND, _ComplexND], object] | None = None,
    fs: onp.ToFloat = ...,  # 2 * pi
    include_nyquist: bool = False,
) -> tuple[_FloatND, _ComplexND]: ...
@overload  # worN: complex
def freqz(
    b: onp.ToComplex | onp.ToComplexND,
    a: onp.ToComplex | onp.ToComplexND = 1,
    worN: op.CanIndex | onp.ToComplex1D = 512,
    whole: op.CanBool = False,
    plot: Callable[[_FloatND, _ComplexND], object] | None = None,
    fs: onp.ToFloat = ...,  # 2 * pi
    include_nyquist: bool = False,
) -> tuple[_InexactND, _ComplexND]: ...

#
@overload  # worN: real
def freqz_zpk(
    z: onp.ToComplex1D,
    p: onp.ToComplex1D,
    k: onp.ToComplex1D,
    worN: op.CanIndex | onp.ToFloat1D = 512,
    whole: op.CanBool = False,
    fs: onp.ToFloat = ...,  # 2 * pi
) -> tuple[_FloatND, _ComplexND]: ...
@overload  # worN: complex
def freqz_zpk(
    z: onp.ToComplex1D,
    p: onp.ToComplex1D,
    k: onp.ToComplex1D,
    worN: onp.ToComplex1D,
    whole: op.CanBool = False,
    fs: onp.ToFloat = ...,  # 2 * pi
) -> tuple[_InexactND, _ComplexND]: ...

#
@overload  # w: real
def group_delay(
    system: tuple[onp.ToComplex1D, onp.ToComplex1D],
    w: op.CanIndex | onp.ToFloat1D | None = 512,
    whole: op.CanBool = False,
    fs: onp.ToFloat = ...,  # 2 * pi
) -> _Ba1D: ...
@overload  # w: complex
def group_delay(
    system: tuple[onp.ToComplex1D, onp.ToComplex1D],
    w: onp.ToComplex1D,
    whole: op.CanBool = False,
    fs: onp.ToFloat = ...,  # 2 * pi
) -> tuple[_Inexact1D, _Float1D]: ...

#
@overload  # worN: real
def freqz_sos(
    sos: onp.ToFloat2D,
    worN: op.CanIndex | onp.ToFloat1D = 512,
    whole: op.CanBool = False,
    fs: onp.ToFloat = ...,  # 2 * pi
) -> tuple[_Float1D, _Complex1D]: ...
@overload  # worN: real
def freqz_sos(
    sos: onp.ToFloat2D,
    worN: onp.ToComplex1D,
    whole: op.CanBool = False,
    fs: onp.ToFloat = ...,  # 2 * pi
) -> tuple[_Inexact1D, _Complex1D]: ...

sosfreqz = freqz_sos

#
def tf2zpk(b: _ToFloat1D, a: _ToFloat1D) -> _ZPK[np.float64 | np.complex128] | _ZPK[np.complex64, np.complex64, np.float32]: ...
def tf2sos(b: _ToFloat1D, a: _ToFloat1D, pairing: _Pairing | None = None, *, analog: onp.ToBool = False) -> _Float2D: ...

#
def zpk2tf(z: onp.ToFloat1D, p: onp.ToFloat1D, k: onp.ToFloat) -> _Ba1D: ...
def zpk2sos(
    z: onp.ToFloat1D, p: onp.ToFloat1D, k: onp.ToFloat, pairing: _Pairing | None = None, *, analog: onp.ToBool = False
) -> _Float2D: ...

#
@overload
def normalize(b: onp.ToFloatStrict1D, a: onp.ToFloat1D) -> _Ba1D[_Floating]: ...
@overload
def normalize(b: onp.ToFloatStrict2D, a: onp.ToFloat1D) -> _Ba2D[_Floating]: ...
@overload
def normalize(b: onp.ToFloat1D | onp.ToFloat2D, a: onp.ToFloat1D) -> _Ba1D[_Floating] | _Ba2D[_Floating]: ...

#
def sos2tf(sos: onp.ToFloat2D) -> tuple[_Floating1D, _Floating1D]: ...
def sos2zpk(sos: _ToFloat2D) -> _ZPK[np.complex128, np.complex128, np.float32 | np.float64]: ...

#
@overload
def lp2lp(b: onp.ToFloatStrict1D, a: onp.ToFloat1D, wo: onp.ToFloat = 1.0) -> _Ba1D | _Ba1D[np.longdouble]: ...
@overload
def lp2lp(b: onp.ToFloatStrict2D, a: onp.ToFloat1D, wo: onp.ToFloat = 1.0) -> _Ba2D | _Ba2D[np.longdouble]: ...
@overload
def lp2lp(
    b: onp.ToFloat1D | onp.ToFloat2D, a: onp.ToFloat1D, wo: onp.ToFloat = 1.0
) -> _Ba1D | _Ba1D[np.longdouble] | _Ba2D | _Ba2D[np.longdouble]: ...

#
def lp2hp(
    b: onp.ToFloat1D, a: onp.ToFloat1D, wo: onp.ToFloat = 1.0
) -> _Ba1D | _Ba1D[np.float16] | _Ba1D[np.float32] | _Ba1D[np.longdouble]: ...

#
def lp2bp(
    b: onp.ToFloat1D, a: onp.ToFloat1D, wo: onp.ToFloat = 1.0, bw: onp.ToFloat = 1.0
) -> _Ba1D | _Ba1D[np.float32] | _Ba1D[np.longdouble]: ...

#
def lp2bs(
    b: onp.ToFloat1D, a: onp.ToFloat1D, wo: onp.ToFloat = 1.0, bw: onp.ToFloat = 1.0
) -> _Ba1D | _Ba1D[np.float32] | _Ba1D[np.longdouble]: ...

#
def lp2lp_zpk(
    z: onp.ToComplex1D, p: onp.ToComplex1D, k: onp.ToFloat, wo: onp.ToFloat = 1.0
) -> _ZPK[npc.inexact, _CFloating, _Floating]: ...

#
def lp2hp_zpk(
    z: onp.ToComplex1D, p: onp.ToComplex1D, k: onp.ToFloat, wo: onp.ToFloat = 1.0
) -> _ZPK[npc.inexact, _CFloating, _Floating]: ...

#
def lp2bp_zpk(
    z: onp.ToComplex1D, p: onp.ToComplex1D, k: onp.ToFloat, wo: onp.ToFloat = 1.0, bw: onp.ToFloat = 1.0
) -> _ZPK[npc.inexact, _CFloating, _Floating]: ...

#
def lp2bs_zpk(
    z: onp.ToComplex1D, p: onp.ToComplex1D, k: onp.ToFloat, wo: onp.ToFloat = 1.0, bw: onp.ToFloat = 1.0
) -> _ZPK[npc.inexact, _CFloating, _Floating]: ...

#
def bilinear(b: onp.ToFloat1D, a: onp.ToFloat1D, fs: onp.ToFloat = 1.0) -> _Ba1D: ...

#
@overload
def bilinear_zpk(
    z: onp.ToFloat1D, p: onp.ToComplex1D, k: onp.ToFloat, fs: onp.ToFloat
) -> _ZPK[np.float64, _CFloating] | _ZPK[np.longdouble, _CFloating, np.longdouble]: ...
@overload
def bilinear_zpk(
    z: onp.ToComplex1D, p: onp.ToComplex1D, k: onp.ToFloat, fs: onp.ToFloat
) -> _ZPK[np.float64 | np.complex128, _CFloating] | _ZPK[np.longdouble | np.clongdouble, _CFloating, np.longdouble]: ...

#
@overload  # output="ba" (default)
def iirdesign(
    wp: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    ws: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    gpass: onp.ToFloat,
    gstop: onp.ToFloat,
    analog: onp.ToBool = False,
    ftype: _FType0 = "ellip",
    output: L["ba"] = "ba",
    fs: onp.ToFloat | None = None,
) -> _Ba1D: ...
@overload  # output="zpk" (positional)
def iirdesign(
    wp: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    ws: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    gpass: onp.ToFloat,
    gstop: onp.ToFloat,
    analog: onp.ToBool,
    ftype: _FType0,
    output: L["zpk"],
    fs: onp.ToFloat | None = None,
) -> _ZPK[np.complex128]: ...
@overload  # output="zpk" (keyword)
def iirdesign(
    wp: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    ws: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    gpass: onp.ToFloat,
    gstop: onp.ToFloat,
    analog: onp.ToBool = False,
    ftype: _FType0 = "ellip",
    *,
    output: L["zpk"],
    fs: onp.ToFloat | None = None,
) -> _ZPK[np.complex128]: ...
@overload  # output="sos" (positional)
def iirdesign(
    wp: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    ws: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    gpass: onp.ToFloat,
    gstop: onp.ToFloat,
    analog: onp.ToBool,
    ftype: _FType0,
    output: L["sos"],
    fs: onp.ToFloat | None = None,
) -> _Float2D: ...
@overload  # output="sos" (keyword)
def iirdesign(
    wp: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    ws: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    gpass: onp.ToFloat,
    gstop: onp.ToFloat,
    analog: onp.ToBool = False,
    ftype: _FType0 = "ellip",
    *,
    output: L["sos"],
    fs: onp.ToFloat | None = None,
) -> _Float2D: ...

#
@overload  # output="ba" (default)
def iirfilter(
    N: onp.ToJustInt,
    Wn: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    rp: onp.ToFloat | None = None,
    rs: onp.ToFloat | None = None,
    btype: _BType = "band",
    analog: onp.ToBool = False,
    ftype: _FType = "butter",
    output: L["ba"] = "ba",
    fs: onp.ToFloat | None = None,
) -> _Ba1D: ...
@overload  # output="zpk" (positional)
def iirfilter(
    N: onp.ToJustInt,
    Wn: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    rp: onp.ToFloat | None,
    rs: onp.ToFloat | None,
    btype: _BType,
    analog: onp.ToBool,
    ftype: _FType,
    output: L["zpk"],
    fs: onp.ToFloat | None = None,
) -> _ZPK[np.complex128]: ...
@overload  # output="zpk" (keyword)
def iirfilter(
    N: onp.ToJustInt,
    Wn: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    rp: onp.ToFloat | None = None,
    rs: onp.ToFloat | None = None,
    btype: _BType = "band",
    analog: onp.ToBool = False,
    ftype: _FType = "butter",
    *,
    output: L["zpk"],
    fs: onp.ToFloat | None = None,
) -> _ZPK[np.complex128]: ...
@overload  # output="sos" (positional)
def iirfilter(
    N: onp.ToJustInt,
    Wn: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    rp: onp.ToFloat | None,
    rs: onp.ToFloat | None,
    btype: _BType,
    analog: onp.ToBool,
    ftype: _FType,
    output: L["sos"],
    fs: onp.ToFloat | None = None,
) -> _Float2D: ...
@overload  # output="sos" (keyword)
def iirfilter(
    N: onp.ToJustInt,
    Wn: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    rp: onp.ToFloat | None = None,
    rs: onp.ToFloat | None = None,
    btype: _BType = "band",
    analog: onp.ToBool = False,
    ftype: _FType = "butter",
    *,
    output: L["sos"],
    fs: onp.ToFloat | None = None,
) -> _Float2D: ...

#
@overload  # output="ba" (default)
def butter(
    N: onp.ToJustInt,
    Wn: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    btype: _BType = "low",
    analog: onp.ToBool = False,
    output: L["ba"] = "ba",
    fs: onp.ToFloat | None = None,
) -> _Ba1D: ...
@overload  # output="zpk" (keyword)
def butter(
    N: onp.ToJustInt,
    Wn: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    btype: _BType = "low",
    analog: onp.ToBool = False,
    *,
    output: L["zpk"],
    fs: onp.ToFloat | None = None,
) -> _ZPK[np.float64, np.complex128, float]: ...
@overload  # output="sos" (keyword)
def butter(
    N: onp.ToJustInt,
    Wn: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    btype: _BType = "low",
    analog: onp.ToBool = False,
    *,
    output: L["sos"],
    fs: onp.ToFloat | None = None,
) -> _Float2D: ...

#
@overload  # output="ba" (default)
def cheby1(
    N: onp.ToJustInt,
    rp: onp.ToFloat,
    Wn: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    btype: _BType = "low",
    analog: onp.ToBool = False,
    output: L["ba"] = "ba",
    fs: onp.ToFloat | None = None,
) -> _Ba1D: ...
@overload  # output="zpk" (positional)
def cheby1(
    N: onp.ToJustInt,
    rp: onp.ToFloat,
    Wn: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    btype: _BType,
    analog: onp.ToBool,
    output: L["zpk"],
    fs: onp.ToFloat | None = None,
) -> _ZPK[np.complex128]: ...
@overload  # output="zpk" (keyword)
def cheby1(
    N: onp.ToJustInt,
    rp: onp.ToFloat,
    Wn: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    btype: _BType = "low",
    analog: onp.ToBool = False,
    *,
    output: L["zpk"],
    fs: onp.ToFloat | None = None,
) -> _ZPK[np.complex128]: ...
@overload  # output="sos" (positional)
def cheby1(
    N: onp.ToJustInt,
    rp: onp.ToFloat,
    Wn: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    btype: _BType,
    analog: onp.ToBool,
    output: L["sos"],
    fs: onp.ToFloat | None = None,
) -> _Float2D: ...
@overload  # output="sos" (keyword)
def cheby1(
    N: onp.ToJustInt,
    rp: onp.ToFloat,
    Wn: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    btype: _BType = "low",
    analog: onp.ToBool = False,
    *,
    output: L["sos"],
    fs: onp.ToFloat | None = None,
) -> _Float2D: ...

#
@overload  # output="ba" (default)
def cheby2(
    N: onp.ToJustInt,
    rs: onp.ToFloat,
    Wn: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    btype: _BType = "low",
    analog: onp.ToBool = False,
    output: L["ba"] = "ba",
    fs: onp.ToFloat | None = None,
) -> _Ba1D: ...
@overload  # output="zpk" (positional)
def cheby2(
    N: onp.ToJustInt,
    rs: onp.ToFloat,
    Wn: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    btype: _BType,
    analog: onp.ToBool,
    output: L["zpk"],
    fs: onp.ToFloat | None = None,
) -> _ZPK[np.complex128]: ...
@overload  # output="zpk" (keyword)
def cheby2(
    N: onp.ToJustInt,
    rs: onp.ToFloat,
    Wn: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    btype: _BType = "low",
    analog: onp.ToBool = False,
    *,
    output: L["zpk"],
    fs: onp.ToFloat | None = None,
) -> _ZPK[np.complex128]: ...
@overload  # output="sos" (positional)
def cheby2(
    N: onp.ToJustInt,
    rs: onp.ToFloat,
    Wn: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    btype: _BType,
    analog: onp.ToBool,
    output: L["sos"],
    fs: onp.ToFloat | None = None,
) -> _Float2D: ...
@overload  # output="sos" (keyword)
def cheby2(
    N: onp.ToJustInt,
    rs: onp.ToFloat,
    Wn: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    btype: _BType = "low",
    analog: onp.ToBool = False,
    *,
    output: L["sos"],
    fs: onp.ToFloat | None = None,
) -> _Float2D: ...

#
@overload  # output="ba" (default)
def ellip(
    N: onp.ToJustInt,
    rp: onp.ToFloat,
    rs: onp.ToFloat,
    Wn: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    btype: _BType = "low",
    analog: onp.ToBool = False,
    output: L["ba"] = "ba",
    fs: onp.ToFloat | None = None,
) -> _Ba1D: ...
@overload  # output="zpk" (postitional)
def ellip(
    N: onp.ToJustInt,
    rp: onp.ToFloat,
    rs: onp.ToFloat,
    Wn: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    btype: _BType,
    analog: onp.ToBool,
    output: L["zpk"],
    fs: onp.ToFloat | None = None,
) -> _ZPK[np.complex128]: ...
@overload  # output="zpk" (keyword)
def ellip(
    N: onp.ToJustInt,
    rp: onp.ToFloat,
    rs: onp.ToFloat,
    Wn: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    btype: _BType = "low",
    analog: onp.ToBool = False,
    *,
    output: L["zpk"],
    fs: onp.ToFloat | None = None,
) -> _ZPK[np.complex128]: ...
@overload  # output="sos" (postitional)
def ellip(
    N: onp.ToJustInt,
    rp: onp.ToFloat,
    rs: onp.ToFloat,
    Wn: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    btype: _BType,
    analog: onp.ToBool,
    output: L["sos"],
    fs: onp.ToFloat | None = None,
) -> _Float2D: ...
@overload  # output="sos" (keyword)
def ellip(
    N: onp.ToJustInt,
    rp: onp.ToFloat,
    rs: onp.ToFloat,
    Wn: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    btype: _BType = "low",
    analog: onp.ToBool = False,
    *,
    output: L["sos"],
    fs: onp.ToFloat | None = None,
) -> _Float2D: ...

#
@overload  # output="ba" (default)
def bessel(
    N: onp.ToJustInt,
    Wn: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    btype: _BType = "low",
    analog: onp.ToBool = False,
    output: L["ba"] = "ba",
    norm: _Normalization = "phase",
    fs: onp.ToFloat | None = None,
) -> _Ba1D: ...
@overload  # output="zpk" (postitional)
def bessel(
    N: onp.ToJustInt,
    Wn: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    btype: _BType,
    analog: onp.ToBool,
    output: L["zpk"],
    norm: _Normalization = "phase",
    fs: onp.ToFloat | None = None,
) -> _ZPK[np.float64]: ...
@overload  # output="zpk" (keyword)
def bessel(
    N: onp.ToJustInt,
    Wn: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    btype: _BType = "low",
    analog: onp.ToBool = False,
    *,
    output: L["zpk"],
    norm: _Normalization = "phase",
    fs: onp.ToFloat | None = None,
) -> _ZPK[np.float64]: ...
@overload  # output="sos" (postitional)
def bessel(
    N: onp.ToJustInt,
    Wn: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    btype: _BType,
    analog: onp.ToBool,
    output: L["sos"],
    norm: _Normalization = "phase",
    fs: onp.ToFloat | None = None,
) -> _Float2D: ...
@overload  # output="sos" (keyword)
def bessel(
    N: onp.ToJustInt,
    Wn: onp.ToFloat | onp.ToFloat1D,  # scalar or length-2
    btype: _BType = "low",
    analog: onp.ToBool = False,
    *,
    output: L["sos"],
    norm: _Normalization = "phase",
    fs: onp.ToFloat | None = None,
) -> _Float2D: ...

#
def band_stop_obj(
    wp: onp.ToFloat,
    ind: L[0, 1] | npc.integer,  # bool doesn't work
    passb: onp.ArrayND[_Floating | npc.integer],  # 1-d
    stopb: onp.ArrayND[_Floating | npc.integer],  # 1-d
    gpass: onp.ToFloat,
    gstop: onp.ToFloat,
    type: L["butter", "cheby", "ellip"],
) -> np.float64 | np.longdouble: ...

#
@overload
def buttord(
    wp: onp.ToFloat,
    ws: onp.ToFloat | onp.ToFloatND,
    gpass: onp.ToFloat,
    gstop: onp.ToFloat,
    analog: onp.ToBool = False,
    fs: onp.ToFloat | None = None,
) -> tuple[int, np.float64 | np.longdouble]: ...
@overload
def buttord(
    wp: onp.ToFloatND,
    ws: onp.ToFloat | onp.ToFloatND,
    gpass: onp.ToFloat,
    gstop: onp.ToFloat,
    analog: onp.ToBool = False,
    fs: onp.ToFloat | None = None,
) -> tuple[int, onp.Array1D[np.float64 | np.longdouble]]: ...

#
@overload
def cheb1ord(
    wp: onp.ToFloat,
    ws: onp.ToFloat | onp.ToFloatND,
    gpass: onp.ToFloat,
    gstop: onp.ToFloat,
    analog: onp.ToBool = False,
    fs: onp.ToFloat | None = None,
) -> tuple[int, _Floating]: ...
@overload
def cheb1ord(
    wp: onp.ToFloatND,
    ws: onp.ToFloat | onp.ToFloatND,
    gpass: onp.ToFloat,
    gstop: onp.ToFloat,
    analog: onp.ToBool = False,
    fs: onp.ToFloat | None = None,
) -> tuple[int, _FloatingND]: ...

#
@overload
def cheb2ord(
    wp: onp.ToFloat,
    ws: onp.ToFloat | onp.ToFloatND,
    gpass: onp.ToFloat,
    gstop: onp.ToFloat,
    analog: onp.ToBool = False,
    fs: onp.ToFloat | None = None,
) -> tuple[int, _Floating]: ...
@overload
def cheb2ord(
    wp: onp.ToFloatND,
    ws: onp.ToFloat | onp.ToFloatND,
    gpass: onp.ToFloat,
    gstop: onp.ToFloat,
    analog: onp.ToBool = False,
    fs: onp.ToFloat | None = None,
) -> tuple[int, _FloatND]: ...  # only nd-output is cast to float64

#
@overload
def ellipord(
    wp: onp.ToFloat,
    ws: onp.ToFloat | onp.ToFloatND,
    gpass: onp.ToFloat,
    gstop: onp.ToFloat,
    analog: onp.ToBool = False,
    fs: onp.ToFloat | None = None,
) -> tuple[int, _Floating]: ...
@overload
def ellipord(
    wp: onp.ToFloatND,
    ws: onp.ToFloat | onp.ToFloatND,
    gpass: onp.ToFloat,
    gstop: onp.ToFloat,
    analog: onp.ToBool = False,
    fs: onp.ToFloat | None = None,
) -> tuple[int, _FloatingND]: ...

#

def buttap(N: onp.ToJustInt) -> tuple[_Float1D, _Complex1D, L[1]]: ...
def cheb1ap(N: onp.ToJustInt, rp: onp.ToFloat) -> tuple[_Float1D, _Complex1D, np.float64]: ...
def cheb2ap(N: onp.ToJustInt, rs: onp.ToFloat) -> tuple[_Complex1D, _Complex1D, np.float64]: ...
def ellipap(N: onp.ToJustInt, rp: onp.ToFloat, rs: onp.ToFloat) -> tuple[_Complex1D, _Complex1D, np.float64]: ...
def besselap(N: onp.ToJustInt, norm: _Normalization = "phase") -> tuple[_Float1D, _Complex1D, float]: ...

#
def iirnotch(w0: onp.ToFloat, Q: onp.ToFloat, fs: onp.ToFloat = 2.0) -> _Ba1D: ...
def iirpeak(w0: onp.ToFloat, Q: onp.ToFloat, fs: onp.ToFloat = 2.0) -> _Ba1D: ...
def iircomb(
    w0: onp.ToFloat, Q: onp.ToFloat, ftype: L["notch", "peak"] = "notch", fs: onp.ToFloat = 2.0, *, pass_zero: onp.ToBool = False
) -> _Ba1D: ...

#
def gammatone(
    freq: onp.ToFloat,
    ftype: L["fir", "iir"],
    order: L[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24] | None = None,
    numtaps: int | None = None,
    fs: onp.ToFloat | None = None,
) -> _Ba1D: ...

# ???
def maxflat() -> None: ...
def yulewalk() -> None: ...
