from typing import Literal, TypeAlias, overload
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from .windows._windows import _ToWindow

__all__ = ["firls", "firwin", "firwin2", "firwin_2d", "kaiser_atten", "kaiser_beta", "kaiserord", "minimum_phase", "remez"]

###

_InexactT = TypeVar("_InexactT", bound=float | npc.inexact)

_IIRFilterType: TypeAlias = Literal["bandpass", "lowpass", "highpass", "bandstop"]
_RemezFilterType: TypeAlias = Literal["bandpass", "differentiator", "hilbert"]
_LinearPhaseFIRMethod: TypeAlias = Literal["homomorphic", "hilbert"]

###

#
def kaiser_beta(a: float) -> float: ...
def kaiser_atten(numtaps: int, width: _InexactT) -> _InexactT: ...
def kaiserord(ripple: float, width: float) -> tuple[int, float]: ...

#
def firwin(
    numtaps: int,
    cutoff: float | onp.ToFloat1D,
    *,
    width: float | None = None,
    window: _ToWindow = "hamming",
    pass_zero: _IIRFilterType | bool = True,
    scale: bool = True,
    fs: float | None = None,
) -> onp.Array1D[np.float64]: ...

#
@overload  # `fc` required, `circular=True`
def firwin_2d(
    hsize: tuple[int, int],
    window: _ToWindow,
    *,
    fc: float | onp.ToFloat1D,
    fs: float = 2,
    circular: Literal[True],
    pass_zero: _IIRFilterType | bool = True,
    scale: bool = True,
) -> onp.Array2D[np.float64]: ...
@overload  # `fc` optional, `circular=False` (default)
def firwin_2d(
    hsize: tuple[int, int],
    window: _ToWindow,
    *,
    fc: float | onp.ToFloat1D | None = None,
    fs: float = 2,
    circular: Literal[False] = False,
    pass_zero: _IIRFilterType | bool = True,
    scale: bool = True,
) -> onp.Array2D[np.float64]: ...

#
def firwin2(
    numtaps: int,
    freq: onp.ToFloat1D,
    gain: onp.ToFloat1D,
    *,
    nfreqs: int | None = None,
    window: _ToWindow | None = "hamming",
    antisymmetric: bool = False,
    fs: float | None = None,
) -> onp.Array1D[np.float64]: ...

#
@overload
def firls(
    numtaps: int, bands: onp.ToFloat1D, desired: onp.ToFloat1D, *, weight: onp.ToFloat1D | None = None, fs: float | None = None
) -> onp.Array1D[np.float64]: ...
@overload
def firls(
    numtaps: int, bands: onp.ToFloat2D, desired: onp.ToFloat2D, *, weight: onp.ToFloat1D | None = None, fs: float | None = None
) -> onp.Array1D[np.float64]: ...

#
def remez(
    numtaps: int,
    bands: onp.ToFloat1D,
    desired: onp.ToFloat1D,
    *,
    weight: onp.ToFloat1D | None = None,
    type: _RemezFilterType = "bandpass",
    maxiter: int = 25,
    grid_density: int = 16,
    fs: float | None = None,
) -> onp.Array1D[np.float64]: ...

#
def minimum_phase(
    h: onp.ToFloat1D, method: _LinearPhaseFIRMethod = "homomorphic", n_fft: int | None = None, *, half: bool = True
) -> onp.Array1D[np.float64]: ...
