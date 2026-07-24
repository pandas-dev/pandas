from typing import Any, Literal, TypeVar, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._ltisys import (
    StateSpaceContinuous,
    StateSpaceDiscrete,
    TransferFunctionContinuous,
    TransferFunctionDiscrete,
    ZerosPolesGainContinuous,
    ZerosPolesGainDiscrete,
    dlti,
    lti,
)

__all__ = ["abcd_normalize", "cont2discrete", "ss2tf", "ss2zpk", "tf2ss", "zpk2ss"]

###

_InexactAT = TypeVar("_InexactAT", bound=npc.inexact, default=np.float64)
_InexactBT = TypeVar("_InexactBT", bound=npc.inexact, default=np.float64)
_InexactCT = TypeVar("_InexactCT", bound=npc.inexact, default=np.float64)
_InexactDT = TypeVar("_InexactDT", bound=npc.inexact, default=np.float64)

type _Tuple4[T] = tuple[T, T, T, T]

type _SystemTF[InexactT: npc.inexact] = tuple[onp.Array2D[InexactT], onp.Array1D[InexactT]]
type _SystemSS[InexactT: npc.inexact] = tuple[
    onp.Array2D[InexactT],
    onp.Array2D[InexactT],
    onp.Array2D[InexactT],
    onp.Array2D[InexactT],
]  # fmt: skip

type _DiscretizeMethod = Literal["gbt", "bilinear", "euler", "backward_diff", "foh", "impulse", "zoh"]
type _DiscreteSS[SafeInexactT: np.float32 | np.float64 | np.complex64 | np.complex128, InexactT: npc.inexact] = tuple[
    onp.Array2D[SafeInexactT],
    onp.Array2D[SafeInexactT],
    onp.Array2D[InexactT],
    onp.Array2D[InexactT],
    float,
]  # fmt: skip

###

@overload  # ~f64, +f64
def tf2ss(
    num: onp.ToArray1D[float, npc.floating64 | npc.integer] | onp.ToArray2D[float, npc.floating64 | npc.integer],
    den: onp.ToFloat64_1D,
) -> _SystemSS[np.float64]: ...
@overload  # +floating, +floating
def tf2ss(
    num: onp.ToFloat1D | onp.ToFloat2D, den: onp.ToFloat1D
) -> tuple[
    onp.Array2D[np.float64 | Any], onp.Array2D[np.float64], onp.Array2D[np.float64 | Any], onp.Array2D[np.float64 | Any]
]: ...
@overload  # ~c128, +c128
def tf2ss(
    num: onp.ToJustComplex128_1D | onp.ToJustComplex128_2D, den: onp.ToComplex128_1D
) -> tuple[onp.Array2D[np.complex128], onp.Array2D[np.float64], onp.Array2D[np.complex128], onp.Array2D[np.complex128]]: ...
@overload  # ~complexfloating, ~complexfloating
def tf2ss(
    num: onp.ToJustComplex1D | onp.ToJustComplex2D, den: onp.ToJustComplex1D
) -> tuple[
    onp.Array2D[np.complex128 | Any], onp.Array2D[np.float64], onp.Array2D[np.complex128 | Any], onp.Array2D[np.complex128 | Any]
]: ...
@overload  # +complexfloating, +complexfloating
def tf2ss(
    num: onp.ToComplex1D | onp.ToComplex2D, den: onp.ToComplex1D
) -> tuple[onp.Array2D[Any], onp.Array2D[np.float64], onp.Array2D[Any], onp.Array2D[Any]]: ...

#
@overload
def abcd_normalize(
    A: onp.ArrayND[_InexactAT] | _InexactAT | onp.ToIntND | float | None = None,
    B: onp.ArrayND[_InexactBT] | _InexactAT | onp.ToIntND | float | None = None,
    C: onp.ArrayND[_InexactCT] | _InexactAT | onp.ToIntND | float | None = None,
    D: onp.ArrayND[_InexactDT] | _InexactAT | onp.ToIntND | float | None = None,
) -> tuple[onp.Array2D[_InexactAT], onp.Array2D[_InexactBT], onp.Array2D[_InexactCT], onp.Array2D[_InexactDT]]: ...
@overload
def abcd_normalize(
    A: onp.ToFloatND | onp.ToFloat | None = None,
    B: onp.ToFloatND | onp.ToFloat | None = None,
    C: onp.ToFloatND | onp.ToFloat | None = None,
    D: onp.ToFloatND | onp.ToFloat | None = None,
) -> _SystemSS[np.float64 | Any]: ...
@overload
def abcd_normalize(
    A: onp.ToComplexND | onp.ToComplex | None = None,
    B: onp.ToComplexND | onp.ToComplex | None = None,
    C: onp.ToComplexND | onp.ToComplex | None = None,
    D: onp.ToComplexND | onp.ToComplex | None = None,
) -> _SystemSS[Any]: ...

#
@overload  # ~f64, +f64, +f64, +f64
def ss2tf(
    A: onp.ToArray2D[float, npc.floating64 | npc.integer],
    B: onp.ToFloat64_2D,
    C: onp.ToFloat64_2D,
    D: onp.ToFloat64_2D,
    input: int = 0,
) -> _SystemTF[np.float64]: ...
@overload  # +f64, ~f64, +f64, +f64
def ss2tf(
    A: onp.ToFloat64_2D,
    B: onp.ToArray2D[float, npc.floating64 | npc.integer],
    C: onp.ToFloat64_2D,
    D: onp.ToFloat64_2D,
    input: int = 0,
) -> _SystemTF[np.float64]: ...
@overload  # +f64, +f64, ~f64, +f64
def ss2tf(
    A: onp.ToFloat64_2D,
    B: onp.ToFloat64_2D,
    C: onp.ToArray2D[float, npc.floating64 | npc.integer],
    D: onp.ToFloat64_2D,
    input: int = 0,
) -> _SystemTF[np.float64]: ...
@overload  # +f64, +f64, +f64, ~f64
def ss2tf(
    A: onp.ToFloat64_2D,
    B: onp.ToFloat64_2D,
    C: onp.ToFloat64_2D,
    D: onp.ToArray2D[float, npc.floating64 | npc.integer],
    input: int = 0,
) -> _SystemTF[np.float64]: ...
@overload  # +floating, +floating, +floating, +floating
def ss2tf(
    A: onp.ToFloat2D, B: onp.ToFloat2D, C: onp.ToFloat2D, D: onp.ToFloat2D, input: int = 0
) -> _SystemTF[np.float64 | Any]: ...
@overload  # ~c128, +c128, +c128, +c128
def ss2tf(
    A: onp.ToJustComplex128_2D, B: onp.ToComplex128_2D, C: onp.ToComplex128_2D, D: onp.ToComplex128_2D, input: int = 0
) -> _SystemTF[np.complex128]: ...
@overload  # +c128, ~c128, +c128, +c128
def ss2tf(
    A: onp.ToComplex128_2D, B: onp.ToJustComplex128_2D, C: onp.ToComplex128_2D, D: onp.ToComplex128_2D, input: int = 0
) -> _SystemTF[np.complex128]: ...
@overload  # +c128, +c128, ~c128, +c128
def ss2tf(
    A: onp.ToComplex128_2D, B: onp.ToComplex128_2D, C: onp.ToJustComplex128_2D, D: onp.ToComplex128_2D, input: int = 0
) -> _SystemTF[np.complex128]: ...
@overload  # +c128, +c128, +c128, ~c128
def ss2tf(
    A: onp.ToComplex128_2D, B: onp.ToComplex128_2D, C: onp.ToComplex128_2D, D: onp.ToJustComplex128_2D, input: int = 0
) -> _SystemTF[np.complex128]: ...
@overload  # +complexfloating, +complexfloating, +complexfloating, +complexfloating
def ss2tf(
    A: onp.ToComplex2D, B: onp.ToComplex2D, C: onp.ToComplex2D, D: onp.ToComplex2D, input: int = 0
) -> _SystemTF[np.complex128 | Any]: ...

#
@overload  # ~f64, +f64
def zpk2ss(z: onp.ToArray1D[float, npc.floating64 | npc.integer], p: onp.ToFloat64_1D, k: float) -> _SystemSS[np.float64]: ...
@overload  # +f64, ~f64
def zpk2ss(z: onp.ToFloat64_1D, p: onp.ToArray1D[float, npc.floating64 | npc.integer], k: float) -> _SystemSS[np.float64]: ...
@overload  # +floating, +floating
def zpk2ss(
    z: onp.ToFloat1D, p: onp.ToFloat1D, k: float
) -> tuple[
    onp.Array2D[np.float64 | Any], onp.Array2D[np.float64], onp.Array2D[np.float64 | Any], onp.Array2D[np.float64 | Any]
]: ...
@overload  # ~c128, +c128
def zpk2ss(
    z: onp.ToJustComplex128_1D, p: onp.ToComplex128_1D, k: float
) -> tuple[onp.Array2D[np.complex128], onp.Array2D[np.float64], onp.Array2D[np.complex128], onp.Array2D[np.complex128]]: ...
@overload  # +c128, ~c128
def zpk2ss(
    z: onp.ToComplex128_1D, p: onp.ToJustComplex128_1D, k: float
) -> tuple[onp.Array2D[np.complex128], onp.Array2D[np.float64], onp.Array2D[np.complex128], onp.Array2D[np.complex128]]: ...
@overload  # +complexfloating, +complexfloating
def zpk2ss(
    z: onp.ToComplex1D, p: onp.ToComplex1D, k: float
) -> tuple[onp.Array2D[Any], onp.Array2D[np.float64], onp.Array2D[Any], onp.Array2D[Any]]: ...

#
@overload  # ~f64, +f64, +f64, +f64
def ss2zpk(
    A: onp.ToArray2D[float, npc.floating64 | npc.integer],
    B: onp.ToFloat64_2D,
    C: onp.ToFloat64_2D,
    D: onp.ToFloat64_2D,
    input: int = 0,
) -> tuple[onp.Array1D[np.complex128], onp.Array1D[np.float64], np.float64]: ...
@overload  # +f64, ~f64, +f64, +f64
def ss2zpk(
    A: onp.ToFloat64_2D,
    B: onp.ToArray2D[float, npc.floating64 | npc.integer],
    C: onp.ToFloat64_2D,
    D: onp.ToFloat64_2D,
    input: int = 0,
) -> tuple[onp.Array1D[np.complex128], onp.Array1D[np.float64], np.float64]: ...
@overload  # +f64, +f64, ~f64, +f64
def ss2zpk(
    A: onp.ToFloat64_2D,
    B: onp.ToFloat64_2D,
    C: onp.ToArray2D[float, npc.floating64 | npc.integer],
    D: onp.ToFloat64_2D,
    input: int = 0,
) -> tuple[onp.Array1D[np.complex128], onp.Array1D[np.float64], np.float64]: ...
@overload  # +f64, +f64, +f64, ~f64
def ss2zpk(
    A: onp.ToFloat64_2D,
    B: onp.ToFloat64_2D,
    C: onp.ToFloat64_2D,
    D: onp.ToArray2D[float, npc.floating64 | npc.integer],
    input: int = 0,
) -> tuple[onp.Array1D[np.complex128], onp.Array1D[np.float64], np.float64]: ...
@overload  # +floating, +floating, +floating, +floating
def ss2zpk(
    A: onp.ToFloat2D, B: onp.ToFloat2D, C: onp.ToFloat2D, D: onp.ToFloat2D, input: int = 0
) -> tuple[onp.Array1D[np.complex128 | Any], onp.Array1D[np.float64 | Any], np.float64 | Any]: ...
@overload  # ~c128, +c128, +c128, +c128
def ss2zpk(
    A: onp.ToJustComplex128_2D, B: onp.ToComplex128_2D, C: onp.ToComplex128_2D, D: onp.ToComplex128_2D, input: int = 0
) -> tuple[onp.Array1D[np.complex128], onp.Array1D[np.complex128], np.complex128]: ...
@overload  # +c128, ~c128, +c128, +c128
def ss2zpk(
    A: onp.ToComplex128_2D, B: onp.ToJustComplex128_2D, C: onp.ToComplex128_2D, D: onp.ToComplex128_2D, input: int = 0
) -> tuple[onp.Array1D[np.complex128], onp.Array1D[np.complex128], np.complex128]: ...
@overload  # +c128, +c128, ~c128, +c128
def ss2zpk(
    A: onp.ToComplex128_2D, B: onp.ToComplex128_2D, C: onp.ToJustComplex128_2D, D: onp.ToComplex128_2D, input: int = 0
) -> tuple[onp.Array1D[np.complex128], onp.Array1D[np.complex128], np.complex128]: ...
@overload  # +c128, +c128, +c128, ~c128
def ss2zpk(
    A: onp.ToComplex128_2D, B: onp.ToComplex128_2D, C: onp.ToComplex128_2D, D: onp.ToJustComplex128_2D, input: int = 0
) -> tuple[onp.Array1D[np.complex128], onp.Array1D[np.complex128], np.complex128]: ...
@overload  # +complexfloating, +complexfloating, +complexfloating, +complexfloating
def ss2zpk(
    A: onp.ToComplex2D, B: onp.ToComplex2D, C: onp.ToComplex2D, D: onp.ToComplex2D, input: int = 0
) -> tuple[onp.Array1D[np.complex128 | Any], onp.Array1D[np.complex128 | Any], np.complex128 | Any]: ...

#
@overload  # TransferFunction
def cont2discrete[SafeFloatT: np.float32 | np.float64](
    system: TransferFunctionContinuous[SafeFloatT], dt: float, method: _DiscretizeMethod = "zoh", alpha: float | None = None
) -> TransferFunctionDiscrete[SafeFloatT]: ...
@overload  # ZerosPolesGain
def cont2discrete[SafeInexactT: np.float32 | np.float64 | np.complex64 | np.complex128, SafeFloatT: np.float32 | np.float64](
    system: ZerosPolesGainContinuous[SafeInexactT, SafeFloatT],
    dt: float,
    method: _DiscretizeMethod = "zoh",
    alpha: float | None = None,
) -> ZerosPolesGainDiscrete[SafeInexactT, SafeFloatT]: ...
@overload  # StateSpace
def cont2discrete[SafeInexactT: np.float32 | np.float64 | np.complex64 | np.complex128, SafeFloatT: np.float32 | np.float64](
    system: StateSpaceContinuous[SafeInexactT, SafeFloatT],
    dt: float,
    method: _DiscretizeMethod = "zoh",
    alpha: float | None = None,
) -> StateSpaceDiscrete[SafeInexactT, SafeFloatT]: ...
@overload  # lti
def cont2discrete[SafeInexactT: np.float32 | np.float64 | np.complex64 | np.complex128, SafeFloatT: np.float32 | np.float64](
    system: lti[SafeInexactT, SafeFloatT], dt: float, method: _DiscretizeMethod = "zoh", alpha: float | None = None
) -> dlti[SafeInexactT, SafeFloatT]: ...
@overload  # (+f64, +f64)
def cont2discrete(
    system: tuple[float | onp.ToFloat64_2D, float | onp.ToFloat64_1D],
    dt: float,
    method: _DiscretizeMethod = "zoh",
    alpha: float | None = None,
) -> tuple[onp.Array2D[np.float64], onp.Array1D[np.float64], float]: ...
@overload  # (+complex, +complex)
def cont2discrete(
    system: tuple[complex | onp.ToComplex2D, complex | onp.ToComplex1D],
    dt: float,
    method: _DiscretizeMethod = "zoh",
    alpha: float | None = None,
) -> tuple[onp.Array2D[np.complex128 | Any], onp.Array1D[np.float64], float]: ...
@overload  # (+f64, +f64, float)
def cont2discrete(
    system: tuple[float | onp.ToFloat64_1D, float | onp.ToFloat64_1D, float],
    dt: float,
    method: _DiscretizeMethod = "zoh",
    alpha: float | None = None,
) -> tuple[onp.Array1D[np.float64], onp.Array1D[np.float64], np.float64, float]: ...
@overload  # (+complex, +complex, complex)
def cont2discrete(
    system: tuple[complex | onp.ToComplex1D, complex | onp.ToComplex1D, complex],
    dt: float,
    method: _DiscretizeMethod = "zoh",
    alpha: float | None = None,
) -> tuple[onp.Array1D[np.complex128 | Any], onp.Array1D[np.complex128 | Any], np.complex128 | Any, float]: ...
@overload  # (~f64; 4)
def cont2discrete(
    system: _Tuple4[onp.ToInt2D | onp.ToJustFloat64_2D], dt: float, method: _DiscretizeMethod = "zoh", alpha: float | None = None
) -> _DiscreteSS[np.float64, np.float64]: ...
@overload  # (~f16; 4)
def cont2discrete(
    system: _Tuple4[onp.ToJustFloat16_2D], dt: float, method: _DiscretizeMethod = "zoh", alpha: float | None = None
) -> _DiscreteSS[np.float64, np.float16]: ...
@overload  # (~f32; 4)
def cont2discrete(
    system: _Tuple4[onp.ToJustFloat32_2D], dt: float, method: _DiscretizeMethod = "zoh", alpha: float | None = None
) -> _DiscreteSS[np.float64, np.float32]: ...
@overload  # (~c128; 4)
def cont2discrete(
    system: _Tuple4[onp.ToJustComplex128_2D], dt: float, method: _DiscretizeMethod = "zoh", alpha: float | None = None
) -> _DiscreteSS[np.complex128, np.complex128]: ...
@overload  # (~c64; 4)
def cont2discrete(
    system: _Tuple4[onp.ToJustComplex64_2D], dt: float, method: _DiscretizeMethod = "zoh", alpha: float | None = None
) -> _DiscreteSS[np.complex128, np.complex64]: ...
