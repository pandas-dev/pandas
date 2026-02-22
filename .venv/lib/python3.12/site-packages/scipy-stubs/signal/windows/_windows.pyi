from _typeshed import Incomplete
from types import ModuleType
from typing import Literal, TypeAlias, overload
from typing_extensions import TypeAliasType

import numpy as np
import optype.numpy as onp

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

_Float64_1D: TypeAlias = onp.Array1D[np.float64]
_Float64_2D: TypeAlias = onp.Array2D[np.float64]

_DeviceNP: TypeAlias = Literal["cpu"]
_Norm: TypeAlias = Literal[2, "approximate", "subsample"]

_ToWindow = TypeAliasType(
    "_ToWindow",
    float | str | tuple[str] | tuple[str, float | onp.ToFloat1D] | tuple[str, float, float | int] | tuple[str, int, int, bool],
)

###

@overload
def get_window(
    window: _ToWindow, Nx: int, fftbins: bool = True, *, xp: None = None, device: _DeviceNP | None = None
) -> _Float64_1D: ...
@overload
def get_window(
    window: _ToWindow, Nx: int, fftbins: bool = True, *, xp: ModuleType, device: Incomplete | None = None
) -> Incomplete: ...

#
@overload
def barthann(M: int, sym: bool = True, *, xp: None = None, device: _DeviceNP | None = None) -> _Float64_1D: ...
@overload
def barthann(M: int, sym: bool = True, *, xp: ModuleType, device: Incomplete | None = None) -> Incomplete: ...

#
@overload
def bartlett(M: int, sym: bool = True, *, xp: None = None, device: _DeviceNP | None = None) -> _Float64_1D: ...
@overload
def bartlett(M: int, sym: bool = True, *, xp: ModuleType, device: Incomplete | None = None) -> Incomplete: ...

#
@overload
def blackman(M: int, sym: bool = True, *, xp: None = None, device: _DeviceNP | None = None) -> _Float64_1D: ...
@overload
def blackman(M: int, sym: bool = True, *, xp: ModuleType, device: Incomplete | None = None) -> Incomplete: ...

#
@overload
def blackmanharris(M: int, sym: bool = True, *, xp: None = None, device: _DeviceNP | None = None) -> _Float64_1D: ...
@overload
def blackmanharris(M: int, sym: bool = True, *, xp: ModuleType, device: Incomplete | None = None) -> Incomplete: ...

#
@overload
def bohman(M: int, sym: bool = True, *, xp: None = None, device: _DeviceNP | None = None) -> _Float64_1D: ...
@overload
def bohman(M: int, sym: bool = True, *, xp: ModuleType, device: Incomplete | None = None) -> Incomplete: ...

#
@overload
def boxcar(M: int, sym: bool = True, *, xp: None = None, device: _DeviceNP | None = None) -> _Float64_1D: ...
@overload
def boxcar(M: int, sym: bool = True, *, xp: ModuleType, device: Incomplete | None = None) -> Incomplete: ...

#
@overload
def cosine(M: int, sym: bool = True, *, xp: None = None, device: _DeviceNP | None = None) -> _Float64_1D: ...
@overload
def cosine(M: int, sym: bool = True, *, xp: ModuleType, device: Incomplete | None = None) -> Incomplete: ...

#
@overload
def flattop(M: int, sym: bool = True, *, xp: None = None, device: _DeviceNP | None = None) -> _Float64_1D: ...
@overload
def flattop(M: int, sym: bool = True, *, xp: ModuleType, device: Incomplete | None = None) -> Incomplete: ...

#
@overload
def hamming(M: int, sym: bool = True, *, xp: None = None, device: _DeviceNP | None = None) -> _Float64_1D: ...
@overload
def hamming(M: int, sym: bool = True, *, xp: ModuleType, device: Incomplete | None = None) -> Incomplete: ...

#
@overload
def hann(M: int, sym: bool = True, *, xp: None = None, device: _DeviceNP | None = None) -> _Float64_1D: ...
@overload
def hann(M: int, sym: bool = True, *, xp: ModuleType, device: Incomplete | None = None) -> Incomplete: ...

#
@overload
def lanczos(M: int, *, sym: bool = True, xp: None = None, device: _DeviceNP | None = None) -> _Float64_1D: ...
@overload
def lanczos(M: int, *, sym: bool = True, xp: ModuleType, device: Incomplete | None = None) -> Incomplete: ...

#
@overload
def nuttall(M: int, sym: bool = True, *, xp: None = None, device: _DeviceNP | None = None) -> _Float64_1D: ...
@overload
def nuttall(M: int, sym: bool = True, *, xp: ModuleType, device: Incomplete | None = None) -> Incomplete: ...

#
@overload
def parzen(M: int, sym: bool = True, *, xp: None = None, device: _DeviceNP | None = None) -> _Float64_1D: ...
@overload
def parzen(M: int, sym: bool = True, *, xp: ModuleType, device: Incomplete | None = None) -> Incomplete: ...

#
@overload
def triang(M: int, sym: bool = True, *, xp: None = None, device: _DeviceNP | None = None) -> _Float64_1D: ...
@overload
def triang(M: int, sym: bool = True, *, xp: ModuleType, device: Incomplete | None = None) -> Incomplete: ...

#
@overload
def chebwin(M: int, at: float, sym: bool = True, *, xp: None = None, device: _DeviceNP | None = None) -> _Float64_1D: ...
@overload
def chebwin(M: int, at: float, sym: bool = True, *, xp: ModuleType, device: Incomplete | None = None) -> Incomplete: ...

#
@overload
def gaussian(M: int, std: float, sym: bool = True, *, xp: None = None, device: _DeviceNP | None = None) -> _Float64_1D: ...
@overload
def gaussian(M: int, std: float, sym: bool = True, *, xp: ModuleType, device: Incomplete | None = None) -> Incomplete: ...

#
@overload
def general_hamming(
    M: int, alpha: float, sym: bool = True, *, xp: None = None, device: _DeviceNP | None = None
) -> _Float64_1D: ...
@overload
def general_hamming(
    M: int, alpha: float, sym: bool = True, *, xp: ModuleType, device: Incomplete | None = None
) -> Incomplete: ...

#
@overload
def kaiser(M: int, beta: float, sym: bool = True, *, xp: None = None, device: _DeviceNP | None = None) -> _Float64_1D: ...
@overload
def kaiser(M: int, beta: float, sym: bool = True, *, xp: ModuleType, device: Incomplete | None = None) -> Incomplete: ...

#
@overload
def kaiser_bessel_derived(
    M: int, beta: float, *, sym: bool = True, xp: None = None, device: _DeviceNP | None = None
) -> _Float64_1D: ...
@overload
def kaiser_bessel_derived(
    M: int, beta: float, *, sym: bool = True, xp: ModuleType, device: Incomplete | None = None
) -> Incomplete: ...

#
@overload
def tukey(M: int, alpha: float = 0.5, sym: bool = True, *, xp: None = None, device: _DeviceNP | None = None) -> _Float64_1D: ...
@overload
def tukey(M: int, alpha: float = 0.5, sym: bool = True, *, xp: ModuleType, device: Incomplete | None = None) -> Incomplete: ...

#
def general_cosine(M: int, a: onp.ToFloat1D, sym: bool = True) -> _Float64_1D: ...

#
@overload
def exponential(
    M: int, center: float | None = None, tau: float = 1.0, sym: bool = True, *, xp: None = None, device: _DeviceNP | None = None
) -> _Float64_1D: ...
@overload
def exponential(
    M: int, center: float | None = None, tau: float = 1.0, sym: bool = True, *, xp: ModuleType, device: Incomplete | None = None
) -> Incomplete: ...

#
@overload
def general_gaussian(
    M: int, p: float, sig: float, sym: bool = True, *, xp: None = None, device: _DeviceNP | None = None
) -> _Float64_1D: ...
@overload
def general_gaussian(
    M: int, p: float, sig: float, sym: bool = True, *, xp: ModuleType, device: Incomplete | None = None
) -> Incomplete: ...

#
@overload
def taylor(
    M: int,
    nbar: onp.ToInt = 4,
    sll: onp.ToInt = 30,
    norm: bool = True,
    sym: bool = True,
    *,
    xp: None = None,
    device: _DeviceNP | None = None,
) -> _Float64_1D: ...
@overload
def taylor(
    M: int,
    nbar: onp.ToInt = 4,
    sll: onp.ToInt = 30,
    norm: bool = True,
    sym: bool = True,
    *,
    xp: ModuleType,
    device: Incomplete | None = None,
) -> Incomplete: ...

#
@overload  # Kmax=None, return_ratios=False
def dpss(
    M: int,
    NW: float,
    Kmax: None = None,
    sym: bool = True,
    norm: _Norm | None = None,
    return_ratios: Literal[False] = False,
    *,
    xp: None = None,
    device: _DeviceNP | None = None,
) -> _Float64_1D: ...
@overload  # Kmax=None, return_ratios=True, /
def dpss(
    M: int,
    NW: float,
    Kmax: None,
    sym: bool,
    norm: _Norm | None,
    return_ratios: Literal[True],
    *,
    xp: None = None,
    device: _DeviceNP | None = None,
) -> tuple[_Float64_1D, np.float64]: ...
@overload  # Kmax=None, *, return_ratios=True
def dpss(
    M: int,
    NW: float,
    Kmax: None = None,
    sym: bool = True,
    norm: _Norm | None = None,
    *,
    return_ratios: Literal[True],
    xp: None = None,
    device: _DeviceNP | None = None,
) -> tuple[_Float64_1D, np.float64]: ...
@overload  # Kmax, return_ratios=False
def dpss(
    M: int,
    NW: float,
    Kmax: int,
    sym: bool = True,
    norm: _Norm | None = None,
    return_ratios: Literal[False] = False,
    *,
    xp: None = None,
    device: _DeviceNP | None = None,
) -> _Float64_2D: ...
@overload  # Kmax, return_ratios=True, /
def dpss(
    M: int,
    NW: float,
    Kmax: int,
    sym: bool,
    norm: _Norm | None,
    return_ratios: Literal[True],
    *,
    xp: None = None,
    device: _DeviceNP | None = None,
) -> tuple[_Float64_2D, _Float64_1D]: ...
@overload  # Kmax, *, return_ratios=True
def dpss(
    M: int,
    NW: float,
    Kmax: int,
    sym: bool = True,
    norm: _Norm | None = None,
    *,
    return_ratios: Literal[True],
    xp: None = None,
    device: _DeviceNP | None = None,
) -> tuple[_Float64_2D, _Float64_1D]: ...
@overload  # xp: ModuleType
def dpss(
    M: int,
    NW: float,
    Kmax: int | None = None,
    sym: bool = True,
    norm: _Norm | None = None,
    return_ratios: bool = False,
    *,
    xp: ModuleType,
    device: Incomplete | None = None,
) -> tuple[Incomplete, Incomplete]: ...
