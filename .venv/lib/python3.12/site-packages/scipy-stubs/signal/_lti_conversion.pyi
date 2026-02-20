from typing import Any, Literal, TypeAlias, TypeVar, overload
from typing_extensions import TypeAliasType

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

_SafeFloatT = TypeVar("_SafeFloatT", bound=np.float32 | np.float64)
_SafeInexactT = TypeVar("_SafeInexactT", bound=np.float32 | np.float64 | np.complex64 | np.complex128)
_InexactT = TypeVar("_InexactT", bound=npc.inexact, default=npc.inexact)

_T = TypeVar("_T")
_Tuple4: TypeAlias = tuple[_T, _T, _T, _T]

_SystemTF: TypeAlias = tuple[onp.Array2D[_InexactT], onp.Array1D[_InexactT]]
_SystemZPK: TypeAlias = tuple[onp.Array1D[_InexactT], onp.Array1D[_InexactT], float | np.float64]
_SystemSS: TypeAlias = tuple[onp.Array2D[_InexactT], onp.Array2D[_InexactT], onp.Array2D[_InexactT], onp.Array2D[_InexactT]]

_DiscretizeMethod: TypeAlias = Literal["gbt", "bilinear", "euler", "backward_diff", "foh", "impulse", "zoh"]
# ty (0.0.9) doens't seem to understand `TypeAlias` with multiple `TypeVar`s
_DiscreteSS = TypeAliasType(
    "_DiscreteSS",
    tuple[onp.Array2D[_SafeInexactT], onp.Array2D[_SafeInexactT], onp.Array2D[_InexactT], onp.Array2D[_InexactT], float],
    type_params=(_SafeInexactT, _InexactT),
)

###

# TODO(@jorenham): refine return dtypes
@overload
def abcd_normalize(
    A: onp.ToFloat2D | None = None, B: onp.ToFloat2D | None = None, C: onp.ToFloat2D | None = None, D: onp.ToFloat2D | None = None
) -> _SystemTF[npc.floating]: ...
@overload
def abcd_normalize(
    A: onp.ToComplex2D | None = None,
    B: onp.ToComplex2D | None = None,
    C: onp.ToComplex2D | None = None,
    D: onp.ToComplex2D | None = None,
) -> _SystemTF: ...

# TODO(@jorenham): refine return dtypes
@overload
def tf2ss(num: onp.ToFloat2D, den: onp.ToFloat1D) -> _SystemSS[npc.floating]: ...
@overload
def tf2ss(num: onp.ToComplex2D, den: onp.ToComplex1D) -> _SystemSS: ...

# TODO(@jorenham): refine return dtypes
@overload
def ss2tf(
    A: onp.ToFloat2D, B: onp.ToFloat2D, C: onp.ToFloat2D, D: onp.ToFloat2D, input: onp.ToInt = 0
) -> _SystemTF[npc.floating]: ...
@overload
def ss2tf(A: onp.ToComplex2D, B: onp.ToComplex2D, C: onp.ToComplex2D, D: onp.ToComplex2D, input: onp.ToInt = 0) -> _SystemTF: ...

# TODO(@jorenham): refine return dtypes
@overload
def zpk2ss(z: onp.ToFloat1D, p: onp.ToFloat1D, k: onp.ToFloat) -> _SystemSS[npc.floating]: ...
@overload
def zpk2ss(z: onp.ToComplex1D, p: onp.ToComplex1D, k: onp.ToFloat) -> _SystemSS: ...

# TODO(@jorenham): refine return dtypes
@overload
def ss2zpk(
    A: onp.ToFloat2D, B: onp.ToFloat2D, C: onp.ToFloat2D, D: onp.ToFloat2D, input: onp.ToInt = 0
) -> _SystemZPK[npc.floating]: ...
@overload
def ss2zpk(
    A: onp.ToComplex2D, B: onp.ToComplex2D, C: onp.ToComplex2D, D: onp.ToComplex2D, input: onp.ToInt = 0
) -> _SystemZPK: ...

#
@overload  # TransferFunction
def cont2discrete(
    system: TransferFunctionContinuous[_SafeFloatT], dt: float, method: _DiscretizeMethod = "zoh", alpha: float | None = None
) -> TransferFunctionDiscrete[_SafeFloatT]: ...
@overload  # ZerosPolesGain
def cont2discrete(
    system: ZerosPolesGainContinuous[_SafeInexactT, _SafeFloatT],
    dt: float,
    method: _DiscretizeMethod = "zoh",
    alpha: float | None = None,
) -> ZerosPolesGainDiscrete[_SafeInexactT, _SafeFloatT]: ...
@overload  # StateSpace
def cont2discrete(
    system: StateSpaceContinuous[_SafeInexactT, _SafeFloatT],
    dt: float,
    method: _DiscretizeMethod = "zoh",
    alpha: float | None = None,
) -> StateSpaceDiscrete[_SafeInexactT, _SafeFloatT]: ...
@overload  # lti
def cont2discrete(
    system: lti[_SafeInexactT, _SafeFloatT], dt: float, method: _DiscretizeMethod = "zoh", alpha: float | None = None
) -> dlti[_SafeInexactT, _SafeFloatT]: ...
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
