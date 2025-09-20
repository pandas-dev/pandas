from collections.abc import Callable
from typing import Literal, TypeAlias, overload

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from .windows._windows import _ToWindow

__all__ = ["check_COLA", "check_NOLA", "coherence", "csd", "istft", "lombscargle", "periodogram", "spectrogram", "stft", "welch"]

_Float1D: TypeAlias = onp.Array1D[np.float64]
_FloatND: TypeAlias = onp.ArrayND[np.float64]
_FloatingND: TypeAlias = onp.ArrayND[np.float32 | np.float64 | np.longdouble]
_CFloatingND: TypeAlias = onp.ArrayND[npc.complexfloating]

_Detrend: TypeAlias = Literal["literal", "constant", False] | Callable[[onp.ArrayND], onp.ArrayND]
_Scaling: TypeAlias = Literal["density", "spectrum"]
_LegacyScaling: TypeAlias = Literal["psd", "spectrum"]
_Average: TypeAlias = Literal["mean", "median"]
_Boundary: TypeAlias = Literal["even", "odd", "constant", "zeros"] | None

###

def lombscargle(
    x: onp.ToFloat1D,
    y: onp.ToFloat1D,
    freqs: onp.ToFloat1D,
    precenter: op.CanBool = False,
    normalize: op.CanBool = False,
    *,
    weights: onp.ToFloat1D | None = None,
    floating_mean: bool = False,
) -> _Float1D: ...

#
def periodogram(
    x: onp.ToComplexND,
    fs: onp.ToFloat = 1.0,
    window: _ToWindow | None = "boxcar",
    nfft: onp.ToInt | None = None,
    detrend: _Detrend = "constant",
    return_onesided: op.CanBool = True,
    scaling: _Scaling = "density",
    axis: op.CanIndex = -1,
) -> tuple[_FloatND, _FloatingND]: ...

#
def welch(
    x: onp.ToComplexND,
    fs: onp.ToFloat = 1.0,
    window: _ToWindow = "hann",
    nperseg: onp.ToInt | None = None,
    noverlap: onp.ToInt | None = None,
    nfft: onp.ToInt | None = None,
    detrend: _Detrend = "constant",
    return_onesided: op.CanBool = True,
    scaling: _Scaling = "density",
    axis: op.CanIndex = -1,
    average: _Average = "mean",
) -> tuple[_FloatND, _FloatingND]: ...

#
def csd(
    x: onp.ToComplexND,
    y: onp.ToComplexND,
    fs: onp.ToFloat = 1.0,
    window: _ToWindow = "hann",
    nperseg: onp.ToInt | None = None,
    noverlap: onp.ToInt | None = None,
    nfft: onp.ToInt | None = None,
    detrend: _Detrend = "constant",
    return_onesided: op.CanBool = True,
    scaling: _Scaling = "density",
    axis: op.CanIndex = -1,
    average: _Average = "mean",
) -> tuple[_FloatND, _CFloatingND]: ...

#
@overload  # non-complex mode (positional and keyword)
def spectrogram(
    x: onp.ToComplexND,
    fs: onp.ToFloat = 1.0,
    window: _ToWindow = ("tukey", 0.25),
    nperseg: onp.ToInt | None = None,
    noverlap: onp.ToInt | None = None,
    nfft: onp.ToInt | None = None,
    detrend: _Detrend = "constant",
    return_onesided: op.CanBool = True,
    scaling: _Scaling = "density",
    axis: op.CanIndex = -1,
    mode: Literal["psd", "magnitude", "angle", "phase"] = "psd",
) -> tuple[_FloatND, _FloatND, _FloatingND]: ...
@overload  # complex mode (positional)
def spectrogram(
    x: onp.ToComplexND,
    fs: onp.ToFloat,
    window: _ToWindow,
    nperseg: onp.ToInt | None,
    noverlap: onp.ToInt | None,
    nfft: onp.ToInt | None,
    detrend: _Detrend,
    return_onesided: op.CanBool,
    scaling: _Scaling,
    axis: op.CanIndex,
    mode: Literal["complex"],
) -> tuple[_FloatND, _FloatND, _CFloatingND]: ...
@overload  # complex mode (keyword)
def spectrogram(
    x: onp.ToComplexND,
    fs: onp.ToFloat = 1.0,
    window: _ToWindow = ("tukey", 0.25),
    nperseg: onp.ToInt | None = None,
    noverlap: onp.ToInt | None = None,
    nfft: onp.ToInt | None = None,
    detrend: _Detrend = "constant",
    return_onesided: op.CanBool = True,
    scaling: _Scaling = "density",
    axis: op.CanIndex = -1,
    *,
    mode: Literal["complex"],
) -> tuple[_FloatND, _FloatND, _CFloatingND]: ...

#
def check_COLA(window: _ToWindow, nperseg: onp.ToInt, noverlap: onp.ToInt, tol: onp.ToFloat = 1e-10) -> np.bool_: ...
def check_NOLA(window: _ToWindow, nperseg: onp.ToInt, noverlap: onp.ToInt, tol: onp.ToFloat = 1e-10) -> np.bool_: ...

#
def stft(
    x: onp.ToComplexND,
    fs: onp.ToFloat = 1.0,
    window: _ToWindow = "hann",
    nperseg: onp.ToInt = 256,
    noverlap: onp.ToInt | None = None,
    nfft: onp.ToInt | None = None,
    detrend: _Detrend = False,
    return_onesided: op.CanBool = True,
    boundary: _Boundary = "zeros",
    padded: op.CanBool = True,
    axis: op.CanIndex = -1,
    scaling: _LegacyScaling = "spectrum",
) -> tuple[_FloatND, _FloatND, _CFloatingND]: ...

#
@overload  # input_onesided is `True`
def istft(
    Zxx: onp.ToComplexND,
    fs: onp.ToFloat = 1.0,
    window: _ToWindow = "hann",
    nperseg: onp.ToInt | None = None,
    noverlap: onp.ToInt | None = None,
    nfft: onp.ToInt | None = None,
    input_onesided: onp.ToTrue = True,
    boundary: op.CanBool = True,
    time_axis: op.CanIndex = -1,
    freq_axis: op.CanIndex = -2,
    scaling: _LegacyScaling = "spectrum",
) -> tuple[_FloatND, _FloatingND]: ...
@overload  # input_onesided is `False` (positional)
def istft(
    Zxx: onp.ToComplexND,
    fs: onp.ToFloat,
    window: _ToWindow,
    nperseg: onp.ToInt | None,
    noverlap: onp.ToInt | None,
    nfft: onp.ToInt | None,
    input_onesided: onp.ToFalse,
    boundary: op.CanBool = True,
    time_axis: op.CanIndex = -1,
    freq_axis: op.CanIndex = -2,
    scaling: _LegacyScaling = "spectrum",
) -> tuple[_FloatND, _CFloatingND]: ...
@overload  # input_onesided is `False` (keyword)
def istft(
    Zxx: onp.ToComplexND,
    fs: onp.ToFloat = 1.0,
    window: _ToWindow = "hann",
    nperseg: onp.ToInt | None = None,
    noverlap: onp.ToInt | None = None,
    nfft: onp.ToInt | None = None,
    *,
    input_onesided: onp.ToFalse,
    boundary: op.CanBool = True,
    time_axis: op.CanIndex = -1,
    freq_axis: op.CanIndex = -2,
    scaling: _LegacyScaling = "spectrum",
) -> tuple[_FloatND, _CFloatingND]: ...

#
def coherence(
    x: onp.ToComplexND,
    y: onp.ToComplexND,
    fs: onp.ToFloat = 1.0,
    window: _ToWindow = "hann",
    nperseg: onp.ToInt | None = None,
    noverlap: onp.ToInt | None = None,
    nfft: onp.ToInt | None = None,
    detrend: _Detrend = "constant",
    axis: op.CanIndex = -1,
) -> tuple[_FloatND, _FloatingND]: ...
