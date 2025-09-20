from _typeshed import Incomplete
from types import ModuleType
from typing import Literal, TypeAlias, overload
from typing_extensions import TypeAliasType

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = [
    "barthann",
    "bartlett",
    "blackman",
    "blackmanharris",
    "bohman",
    "boxcar",
    "chebwin",
    "cosine",
    "dpss",
    "exponential",
    "flattop",
    "gaussian",
    "general_cosine",
    "general_gaussian",
    "general_hamming",
    "get_window",
    "hamming",
    "hann",
    "kaiser",
    "kaiser_bessel_derived",
    "lanczos",
    "nuttall",
    "parzen",
    "taylor",
    "triang",
    "tukey",
]

###

_Xp = TypeAliasType("_Xp", ModuleType)
_Device = TypeAliasType("_Device", Incomplete)

_Float64_1D: TypeAlias = onp.Array1D[np.float64]
_Float64_2D: TypeAlias = onp.Array2D[np.float64]

_ToBool: TypeAlias = bool | np.bool_
_ToInt: TypeAlias = int | np.int16 | np.int32 | np.int64
_ToFloat: TypeAlias = float | npc.floating | npc.integer

_Norm: TypeAlias = Literal[2, "approximate", "subsample"]

_ToWindow = TypeAliasType(
    "_ToWindow",
    _ToFloat
    | str
    | tuple[str]
    | tuple[str, _ToFloat | onp.ToFloat1D]
    | tuple[str, _ToFloat, _ToFloat | op.CanIndex]
    | tuple[str, _ToInt, _ToInt, bool],
)

###

def get_window(
    window: _ToWindow, Nx: _ToInt, fftbins: _ToBool = True, *, xp: _Xp | None = None, device: _Device | None = None
) -> _Float64_1D: ...

#
def barthann(M: _ToInt, sym: _ToBool = True, *, xp: _Xp | None = None, device: _Device | None = None) -> _Float64_1D: ...
def bartlett(M: _ToInt, sym: _ToBool = True, *, xp: _Xp | None = None, device: _Device | None = None) -> _Float64_1D: ...
def blackman(M: _ToInt, sym: _ToBool = True, *, xp: _Xp | None = None, device: _Device | None = None) -> _Float64_1D: ...
def blackmanharris(M: _ToInt, sym: _ToBool = True, *, xp: _Xp | None = None, device: _Device | None = None) -> _Float64_1D: ...
def bohman(M: _ToInt, sym: _ToBool = True, *, xp: _Xp | None = None, device: _Device | None = None) -> _Float64_1D: ...
def boxcar(M: _ToInt, sym: _ToBool = True, *, xp: _Xp | None = None, device: _Device | None = None) -> _Float64_1D: ...
def cosine(M: _ToInt, sym: _ToBool = True, *, xp: _Xp | None = None, device: _Device | None = None) -> _Float64_1D: ...
def flattop(M: _ToInt, sym: _ToBool = True, *, xp: _Xp | None = None, device: _Device | None = None) -> _Float64_1D: ...
def hamming(M: _ToInt, sym: _ToBool = True, *, xp: _Xp | None = None, device: _Device | None = None) -> _Float64_1D: ...
def hann(M: _ToInt, sym: _ToBool = True, *, xp: _Xp | None = None, device: _Device | None = None) -> _Float64_1D: ...
def lanczos(M: _ToInt, *, sym: _ToBool = True, xp: _Xp | None = None, device: _Device | None = None) -> _Float64_1D: ...
def nuttall(M: _ToInt, sym: _ToBool = True, *, xp: _Xp | None = None, device: _Device | None = None) -> _Float64_1D: ...
def parzen(M: _ToInt, sym: _ToBool = True, *, xp: _Xp | None = None, device: _Device | None = None) -> _Float64_1D: ...
def triang(M: _ToInt, sym: _ToBool = True, *, xp: _Xp | None = None, device: _Device | None = None) -> _Float64_1D: ...

#
def chebwin(
    M: _ToInt, at: _ToFloat, sym: _ToBool = True, *, xp: _Xp | None = None, device: _Device | None = None
) -> _Float64_1D: ...
def gaussian(
    M: _ToInt, std: _ToFloat, sym: _ToBool = True, *, xp: _Xp | None = None, device: _Device | None = None
) -> _Float64_1D: ...
def general_hamming(
    M: _ToInt, alpha: _ToFloat, sym: _ToBool = True, *, xp: _Xp | None = None, device: _Device | None = None
) -> _Float64_1D: ...
def kaiser(
    M: _ToInt, beta: _ToFloat, sym: _ToBool = True, *, xp: _Xp | None = None, device: _Device | None = None
) -> _Float64_1D: ...
def kaiser_bessel_derived(
    M: _ToInt, beta: _ToFloat, *, sym: _ToBool = True, xp: _Xp | None = None, device: _Device | None = None
) -> _Float64_1D: ...
def tukey(
    M: _ToInt, alpha: _ToFloat = 0.5, sym: _ToBool = True, *, xp: _Xp | None = None, device: _Device | None = None
) -> _Float64_1D: ...

#
def general_cosine(M: _ToInt, a: onp.ToFloat1D, sym: _ToBool = True) -> _Float64_1D: ...

#
def exponential(
    M: _ToInt,
    center: _ToFloat | None = None,
    tau: _ToFloat = 1.0,
    sym: _ToBool = True,
    *,
    xp: _Xp | None = None,
    device: _Device | None = None,
) -> _Float64_1D: ...
def general_gaussian(
    M: _ToInt, p: _ToFloat, sig: _ToFloat, sym: _ToBool = True, *, xp: _Xp | None = None, device: _Device | None = None
) -> _Float64_1D: ...

#
def taylor(
    M: _ToInt,
    nbar: onp.ToInt = 4,
    sll: onp.ToInt = 30,
    norm: _ToBool = True,
    sym: _ToBool = True,
    *,
    xp: _Xp | None = None,
    device: _Device | None = None,
) -> _Float64_1D: ...

#
@overload  # Kmax=None, return_ratios=False
def dpss(
    M: _ToInt,
    NW: _ToFloat,
    Kmax: None = None,
    sym: _ToBool = True,
    norm: _Norm | None = None,
    return_ratios: onp.ToFalse = False,
    *,
    xp: _Xp | None = None,
    device: _Device | None = None,
) -> _Float64_1D: ...
@overload  # Kmax=None, return_ratios=True, /
def dpss(
    M: _ToInt,
    NW: _ToFloat,
    Kmax: None,
    sym: _ToBool,
    norm: _Norm | None,
    return_ratios: onp.ToTrue,
    *,
    xp: _Xp | None = None,
    device: _Device | None = None,
) -> tuple[_Float64_1D, np.float64]: ...
@overload  # Kmax=None, *, return_ratios=True
def dpss(
    M: _ToInt,
    NW: _ToFloat,
    Kmax: None = None,
    sym: _ToBool = True,
    norm: _Norm | None = None,
    *,
    return_ratios: onp.ToTrue,
    xp: _Xp | None = None,
    device: _Device | None = None,
) -> tuple[_Float64_1D, np.float64]: ...
@overload  # Kmax, return_ratios=False
def dpss(
    M: _ToInt,
    NW: _ToFloat,
    Kmax: op.CanIndex,
    sym: _ToBool = True,
    norm: _Norm | None = None,
    return_ratios: onp.ToFalse = False,
    *,
    xp: _Xp | None = None,
    device: _Device | None = None,
) -> _Float64_2D: ...
@overload  # Kmax, return_ratios=True, /
def dpss(
    M: _ToInt,
    NW: _ToFloat,
    Kmax: op.CanIndex,
    sym: _ToBool,
    norm: _Norm | None,
    return_ratios: onp.ToTrue,
    *,
    xp: _Xp | None = None,
    device: _Device | None = None,
) -> tuple[_Float64_2D, _Float64_1D]: ...
@overload  # Kmax, *, return_ratios=True
def dpss(
    M: _ToInt,
    NW: _ToFloat,
    Kmax: op.CanIndex,
    sym: _ToBool = True,
    norm: _Norm | None = None,
    *,
    return_ratios: onp.ToTrue,
    xp: _Xp | None = None,
    device: _Device | None = None,
) -> tuple[_Float64_2D, _Float64_1D]: ...
