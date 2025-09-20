from typing import Literal, TypeAlias, TypeVar, overload

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

_ToSystemTF: TypeAlias = tuple[onp.ToComplex2D, onp.ToComplex1D]  # (num, den)
_ToSystemZPK: TypeAlias = tuple[onp.ToComplex1D, onp.ToComplex1D, onp.ToFloat]  # (z, p, k)
_ToSystemSS: TypeAlias = tuple[onp.ToComplex2D, onp.ToComplex2D, onp.ToComplex2D, onp.ToComplex2D]  # (A, B, C, D)

_InexactT = TypeVar("_InexactT", bound=npc.inexact, default=npc.inexact)
_Inexact1D: TypeAlias = onp.Array1D[_InexactT]
_Inexact2D: TypeAlias = onp.Array2D[_InexactT]

_SystemTF: TypeAlias = tuple[_Inexact2D[_InexactT], _Inexact1D[_InexactT]]
_SystemZPK: TypeAlias = tuple[_Inexact1D[_InexactT], _Inexact1D[_InexactT], float | np.float64]
_SystemSS: TypeAlias = tuple[_Inexact2D[_InexactT], _Inexact2D[_InexactT], _Inexact2D[_InexactT], _Inexact2D[_InexactT]]

_DiscretizeMethod: TypeAlias = Literal["gbt", "bilinear", "euler", "backward_diff", "foh", "impulse", "zoh"]

###
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

#
@overload
def tf2ss(num: onp.ToFloat2D, den: onp.ToFloat1D) -> _SystemSS[npc.floating]: ...
@overload
def tf2ss(num: onp.ToComplex2D, den: onp.ToComplex1D) -> _SystemSS: ...

#
@overload
def ss2tf(
    A: onp.ToFloat2D, B: onp.ToFloat2D, C: onp.ToFloat2D, D: onp.ToFloat2D, input: onp.ToInt = 0
) -> _SystemTF[npc.floating]: ...
@overload
def ss2tf(A: onp.ToComplex2D, B: onp.ToComplex2D, C: onp.ToComplex2D, D: onp.ToComplex2D, input: onp.ToInt = 0) -> _SystemTF: ...

#
@overload
def zpk2ss(z: onp.ToFloat1D, p: onp.ToFloat1D, k: onp.ToFloat) -> _SystemSS[npc.floating]: ...
@overload
def zpk2ss(z: onp.ToComplex1D, p: onp.ToComplex1D, k: onp.ToFloat) -> _SystemSS: ...

#
@overload
def ss2zpk(
    A: onp.ToFloat2D, B: onp.ToFloat2D, C: onp.ToFloat2D, D: onp.ToFloat2D, input: onp.ToInt = 0
) -> _SystemZPK[npc.floating]: ...
@overload
def ss2zpk(
    A: onp.ToComplex2D, B: onp.ToComplex2D, C: onp.ToComplex2D, D: onp.ToComplex2D, input: onp.ToInt = 0
) -> _SystemZPK: ...

#
@overload
def cont2discrete(
    system: TransferFunctionContinuous, dt: float, method: _DiscretizeMethod = "zoh", alpha: onp.ToJustFloat | None = None
) -> TransferFunctionDiscrete: ...
@overload
def cont2discrete(
    system: ZerosPolesGainContinuous, dt: float, method: _DiscretizeMethod = "zoh", alpha: onp.ToJustFloat | None = None
) -> ZerosPolesGainDiscrete: ...
@overload
def cont2discrete(
    system: StateSpaceContinuous, dt: float, method: _DiscretizeMethod = "zoh", alpha: onp.ToJustFloat | None = None
) -> StateSpaceDiscrete: ...
@overload
def cont2discrete(system: lti, dt: float, method: _DiscretizeMethod = "zoh", alpha: onp.ToJustFloat | None = None) -> dlti: ...
@overload
def cont2discrete(
    system: _ToSystemTF, dt: float, method: _DiscretizeMethod = "zoh", alpha: onp.ToJustFloat | None = None
) -> tuple[_Inexact2D[_InexactT], _Inexact1D[_InexactT], float]: ...
@overload
def cont2discrete(
    system: _ToSystemZPK, dt: float, method: _DiscretizeMethod = "zoh", alpha: onp.ToJustFloat | None = None
) -> tuple[_Inexact1D[_InexactT], _Inexact1D[_InexactT], float, float]: ...
@overload
def cont2discrete(
    system: _ToSystemSS, dt: float, method: _DiscretizeMethod = "zoh", alpha: onp.ToJustFloat | None = None
) -> tuple[_Inexact2D[_InexactT], _Inexact2D[_InexactT], _Inexact2D[_InexactT], _Inexact2D[_InexactT], float]: ...
