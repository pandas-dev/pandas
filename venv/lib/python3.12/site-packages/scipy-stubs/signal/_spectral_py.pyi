from collections.abc import Callable
from typing import Literal, TypeAlias, overload
from typing_extensions import deprecated

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from .windows._windows import _ToWindow

__all__ = ["check_COLA", "check_NOLA", "coherence", "csd", "istft", "lombscargle", "periodogram", "spectrogram", "stft", "welch"]

###

_float64_1d: TypeAlias = onp.Array1D[np.float64]  # noqa: PYI042
_float32_nd: TypeAlias = onp.ArrayND[np.float32]  # noqa: PYI042
_float64_nd: TypeAlias = onp.ArrayND[np.float64]  # noqa: PYI042
_float80_nd: TypeAlias = onp.ArrayND[np.float96 | np.float128]  # noqa: PYI042
_complex64_nd: TypeAlias = onp.ArrayND[np.complex64]  # noqa: PYI042
_complex128_nd: TypeAlias = onp.ArrayND[np.complex128]  # noqa: PYI042
_complex160_nd: TypeAlias = onp.ArrayND[np.complex192 | np.complex256]  # noqa: PYI042

_ToInexact32ND: TypeAlias = onp.ToArrayND[npc.inexact32, npc.inexact32 | npc.floating16]
_ToInexact64ND: TypeAlias = onp.ToArrayND[complex, npc.inexact64 | npc.integer | np.bool_]
_ToInexact80ND: TypeAlias = onp.ToArrayND[npc.inexact80, npc.inexact80]

_Detrend: TypeAlias = Literal["literal", "constant", False] | Callable[[onp.ArrayND], onp.ArrayND]
_Scaling: TypeAlias = Literal["density", "spectrum"]
_LegacyScaling: TypeAlias = Literal["psd", "spectrum"]
_Average: TypeAlias = Literal["mean", "median"]
_Boundary: TypeAlias = Literal["even", "odd", "constant", "zeros"] | None
_Normalize: TypeAlias = Literal["power", "normalize", "amplitude"] | bool

###

@overload
def lombscargle(
    x: onp.ToFloat1D,
    y: onp.ToFloat1D,
    freqs: onp.ToFloat1D,
    *,
    precenter: op.JustObject = ...,
    normalize: _Normalize = False,
    weights: onp.ToFloat1D | None = None,
    floating_mean: bool = False,
) -> _float64_1d: ...
@overload
@deprecated(
    "The `precenter` argument is deprecated and will be removed in SciPy 1.19.0. "
    "The functionality can be substituted by passing `y - y.mean()` to `y`."
)
def lombscargle(
    x: onp.ToFloat1D,
    y: onp.ToFloat1D,
    freqs: onp.ToFloat1D,
    *,
    precenter: bool,
    normalize: _Normalize = False,
    weights: onp.ToFloat1D | None = None,
    floating_mean: bool = False,
) -> _float64_1d: ...

#
@overload  # f64
def periodogram(  # type: ignore[overload-overlap]
    x: _ToInexact64ND,
    fs: float = 1.0,
    window: _ToWindow | None = "boxcar",
    nfft: int | None = None,
    detrend: _Detrend = "constant",
    return_onesided: bool = True,
    scaling: _Scaling = "density",
    axis: int = -1,
) -> tuple[_float64_1d, _float64_nd]: ...
@overload  # f32
def periodogram(
    x: _ToInexact32ND,
    fs: float = 1.0,
    window: _ToWindow | None = "boxcar",
    nfft: int | None = None,
    detrend: _Detrend = "constant",
    return_onesided: bool = True,
    scaling: _Scaling = "density",
    axis: int = -1,
) -> tuple[_float64_1d, _float32_nd]: ...
@overload  # f80
def periodogram(
    x: _ToInexact80ND,
    fs: float = 1.0,
    window: _ToWindow | None = "boxcar",
    nfft: int | None = None,
    detrend: _Detrend = "constant",
    return_onesided: bool = True,
    scaling: _Scaling = "density",
    axis: int = -1,
) -> tuple[_float64_1d, _float80_nd]: ...

#
@overload  # f64
def welch(  # type: ignore[overload-overlap]
    x: _ToInexact64ND,
    fs: float = 1.0,
    window: _ToWindow = "hann_periodic",
    nperseg: int | None = None,
    noverlap: int | None = None,
    nfft: int | None = None,
    detrend: _Detrend = "constant",
    return_onesided: bool = True,
    scaling: _Scaling = "density",
    axis: int = -1,
    average: _Average = "mean",
) -> tuple[_float64_1d, _float64_nd]: ...
@overload  # f64
def welch(
    x: _ToInexact32ND,
    fs: float = 1.0,
    window: _ToWindow = "hann_periodic",
    nperseg: int | None = None,
    noverlap: int | None = None,
    nfft: int | None = None,
    detrend: _Detrend = "constant",
    return_onesided: bool = True,
    scaling: _Scaling = "density",
    axis: int = -1,
    average: _Average = "mean",
) -> tuple[_float64_1d, _float32_nd]: ...
@overload  # f64
def welch(
    x: _ToInexact80ND,
    fs: float = 1.0,
    window: _ToWindow = "hann_periodic",
    nperseg: int | None = None,
    noverlap: int | None = None,
    nfft: int | None = None,
    detrend: _Detrend = "constant",
    return_onesided: bool = True,
    scaling: _Scaling = "density",
    axis: int = -1,
    average: _Average = "mean",
) -> tuple[_float64_1d, _float80_nd]: ...

# NOTE: We assume that `x is not y` always holds here.
# See https://github.com/scipy/scipy/issues/24285 for details.
@overload  # c128
def csd(  # type: ignore[overload-overlap]
    x: _ToInexact64ND,
    y: _ToInexact64ND,
    fs: float = 1.0,
    window: _ToWindow = "hann_periodic",
    nperseg: int | None = None,
    noverlap: int | None = None,
    nfft: int | None = None,
    detrend: _Detrend = "constant",
    return_onesided: bool = True,
    scaling: _Scaling = "density",
    axis: int = -1,
    average: _Average = "mean",
) -> tuple[_float64_1d, _complex128_nd]: ...
@overload  # c64
def csd(
    x: _ToInexact32ND,
    y: _ToInexact32ND,
    fs: float = 1.0,
    window: _ToWindow = "hann_periodic",
    nperseg: int | None = None,
    noverlap: int | None = None,
    nfft: int | None = None,
    detrend: _Detrend = "constant",
    return_onesided: bool = True,
    scaling: _Scaling = "density",
    axis: int = -1,
    average: _Average = "mean",
) -> tuple[_float64_1d, _complex64_nd]: ...
@overload  # c160
def csd(
    x: _ToInexact80ND,
    y: _ToInexact80ND,
    fs: float = 1.0,
    window: _ToWindow = "hann_periodic",
    nperseg: int | None = None,
    noverlap: int | None = None,
    nfft: int | None = None,
    detrend: _Detrend = "constant",
    return_onesided: bool = True,
    scaling: _Scaling = "density",
    axis: int = -1,
    average: _Average = "mean",
) -> tuple[_float64_1d, _complex160_nd]: ...
@overload  # fallback
def csd(
    x: onp.ToComplexND,
    y: onp.ToComplexND,
    fs: float = 1.0,
    window: _ToWindow = "hann_periodic",
    nperseg: int | None = None,
    noverlap: int | None = None,
    nfft: int | None = None,
    detrend: _Detrend = "constant",
    return_onesided: bool = True,
    scaling: _Scaling = "density",
    axis: int = -1,
    average: _Average = "mean",
) -> tuple[_float64_1d, onp.Array]: ...

# NOTE: Even though it is theoretically possible to pass `mode` as positional argument, it's unlikely to be done in practice,
# and would significantly complicate the overloads. Thus, we only support passing `mode` as keyword argument here (if "complex").
@overload  # f64, mode != "complex"
def spectrogram(  # type: ignore[overload-overlap]
    x: _ToInexact64ND,
    fs: float = 1.0,
    window: _ToWindow = ("tukey_periodic", 0.25),
    nperseg: int | None = None,
    noverlap: int | None = None,
    nfft: int | None = None,
    detrend: _Detrend = "constant",
    return_onesided: bool = True,
    scaling: _Scaling = "density",
    axis: int = -1,
    mode: Literal["psd", "magnitude", "angle", "phase"] = "psd",
) -> tuple[_float64_1d, _float64_1d, _float64_nd]: ...
@overload  # c128, mode == "complex"
def spectrogram(  # type: ignore[overload-overlap]
    x: _ToInexact64ND,
    fs: float = 1.0,
    window: _ToWindow = ("tukey_periodic", 0.25),
    nperseg: int | None = None,
    noverlap: int | None = None,
    nfft: int | None = None,
    detrend: _Detrend = "constant",
    return_onesided: bool = True,
    scaling: _Scaling = "density",
    axis: int = -1,
    *,
    mode: Literal["complex"],
) -> tuple[_float64_1d, _float64_1d, _complex128_nd]: ...
@overload  # f32, mode != "complex"
def spectrogram(
    x: _ToInexact32ND,
    fs: float = 1.0,
    window: _ToWindow = ("tukey_periodic", 0.25),
    nperseg: int | None = None,
    noverlap: int | None = None,
    nfft: int | None = None,
    detrend: _Detrend = "constant",
    return_onesided: bool = True,
    scaling: _Scaling = "density",
    axis: int = -1,
    mode: Literal["psd", "magnitude", "angle", "phase"] = "psd",
) -> tuple[_float64_1d, _float64_1d, _float32_nd]: ...
@overload  # c64, mode == "complex"
def spectrogram(
    x: _ToInexact32ND,
    fs: float = 1.0,
    window: _ToWindow = ("tukey_periodic", 0.25),
    nperseg: int | None = None,
    noverlap: int | None = None,
    nfft: int | None = None,
    detrend: _Detrend = "constant",
    return_onesided: bool = True,
    scaling: _Scaling = "density",
    axis: int = -1,
    *,
    mode: Literal["complex"],
) -> tuple[_float64_1d, _float64_1d, _complex64_nd]: ...
@overload  # f80, mode != "complex"
def spectrogram(
    x: _ToInexact80ND,
    fs: float = 1.0,
    window: _ToWindow = ("tukey_periodic", 0.25),
    nperseg: int | None = None,
    noverlap: int | None = None,
    nfft: int | None = None,
    detrend: _Detrend = "constant",
    return_onesided: bool = True,
    scaling: _Scaling = "density",
    axis: int = -1,
    mode: Literal["psd", "magnitude", "angle", "phase"] = "psd",
) -> tuple[_float64_1d, _float64_1d, _float80_nd]: ...
@overload  # c80, mode == "complex"
def spectrogram(
    x: _ToInexact80ND,
    fs: float = 1.0,
    window: _ToWindow = ("tukey_periodic", 0.25),
    nperseg: int | None = None,
    noverlap: int | None = None,
    nfft: int | None = None,
    detrend: _Detrend = "constant",
    return_onesided: bool = True,
    scaling: _Scaling = "density",
    axis: int = -1,
    *,
    mode: Literal["complex"],
) -> tuple[_float64_1d, _float64_1d, _complex160_nd]: ...
@overload  # fallback
def spectrogram(
    x: onp.ToComplexND,
    fs: float = 1.0,
    window: _ToWindow = ("tukey_periodic", 0.25),
    nperseg: int | None = None,
    noverlap: int | None = None,
    nfft: int | None = None,
    detrend: _Detrend = "constant",
    return_onesided: bool = True,
    scaling: _Scaling = "density",
    axis: int = -1,
    mode: str = "psd",
) -> tuple[_float64_1d, _float64_1d, onp.Array]: ...

#
def check_COLA(window: _ToWindow, nperseg: int, noverlap: int, tol: float = 1e-10) -> np.bool_: ...
def check_NOLA(window: _ToWindow, nperseg: int, noverlap: int, tol: float = 1e-10) -> np.bool_: ...

#
@overload  # c128
def stft(  # type: ignore[overload-overlap]
    x: _ToInexact64ND,
    fs: float = 1.0,
    window: _ToWindow = "hann_periodic",
    nperseg: int = 256,
    noverlap: int | None = None,
    nfft: int | None = None,
    detrend: _Detrend = False,
    return_onesided: bool = True,
    boundary: _Boundary = "zeros",
    padded: bool = True,
    axis: int = -1,
    scaling: _LegacyScaling = "spectrum",
) -> tuple[_float64_1d, _float64_1d, _complex128_nd]: ...
@overload  # c64
def stft(
    x: _ToInexact32ND,
    fs: float = 1.0,
    window: _ToWindow = "hann_periodic",
    nperseg: int = 256,
    noverlap: int | None = None,
    nfft: int | None = None,
    detrend: _Detrend = False,
    return_onesided: bool = True,
    boundary: _Boundary = "zeros",
    padded: bool = True,
    axis: int = -1,
    scaling: _LegacyScaling = "spectrum",
) -> tuple[_float64_1d, _float64_1d, _complex64_nd]: ...
@overload  # c160
def stft(
    x: _ToInexact80ND,
    fs: float = 1.0,
    window: _ToWindow = "hann_periodic",
    nperseg: int = 256,
    noverlap: int | None = None,
    nfft: int | None = None,
    detrend: _Detrend = False,
    return_onesided: bool = True,
    boundary: _Boundary = "zeros",
    padded: bool = True,
    axis: int = -1,
    scaling: _LegacyScaling = "spectrum",
) -> tuple[_float64_1d, _float64_1d, _complex160_nd]: ...
@overload  # fallback
def stft(
    x: onp.ToComplexND,
    fs: float = 1.0,
    window: _ToWindow = "hann_periodic",
    nperseg: int = 256,
    noverlap: int | None = None,
    nfft: int | None = None,
    detrend: _Detrend = False,
    return_onesided: bool = True,
    boundary: _Boundary = "zeros",
    padded: bool = True,
    axis: int = -1,
    scaling: _LegacyScaling = "spectrum",
) -> tuple[_float64_1d, _float64_1d, onp.Array]: ...

# NOTE: Even though it is theoretically possible to pass `input_onesided` positionally, it's unlikely to be done in practice
# and would significantly complicate the overloads. Thus, we only support passing it as keyword argument here (if `False`).
@overload  # f64, input_onesided=True (default)
def istft(  # type: ignore[overload-overlap]
    Zxx: _ToInexact64ND,
    fs: float = 1.0,
    window: _ToWindow = "hann_periodic",
    nperseg: int | None = None,
    noverlap: int | None = None,
    nfft: int | None = None,
    input_onesided: Literal[True] = True,
    boundary: bool = True,
    time_axis: int = -1,
    freq_axis: int = -2,
    scaling: _LegacyScaling = "spectrum",
) -> tuple[_float64_1d, _float64_nd]: ...
@overload  # c128, input_onesided=False
def istft(  # type: ignore[overload-overlap]
    Zxx: _ToInexact64ND,
    fs: float = 1.0,
    window: _ToWindow = "hann_periodic",
    nperseg: int | None = None,
    noverlap: int | None = None,
    nfft: int | None = None,
    *,
    input_onesided: Literal[False],
    boundary: bool = True,
    time_axis: int = -1,
    freq_axis: int = -2,
    scaling: _LegacyScaling = "spectrum",
) -> tuple[_float64_1d, _complex128_nd]: ...
@overload  # f32, input_onesided=True (default)
def istft(
    Zxx: _ToInexact32ND,
    fs: float = 1.0,
    window: _ToWindow = "hann_periodic",
    nperseg: int | None = None,
    noverlap: int | None = None,
    nfft: int | None = None,
    input_onesided: Literal[True] = True,
    boundary: bool = True,
    time_axis: int = -1,
    freq_axis: int = -2,
    scaling: _LegacyScaling = "spectrum",
) -> tuple[_float64_1d, _float32_nd]: ...
@overload  # c64, input_onesided=False
def istft(
    Zxx: _ToInexact32ND,
    fs: float = 1.0,
    window: _ToWindow = "hann_periodic",
    nperseg: int | None = None,
    noverlap: int | None = None,
    nfft: int | None = None,
    *,
    input_onesided: Literal[False],
    boundary: bool = True,
    time_axis: int = -1,
    freq_axis: int = -2,
    scaling: _LegacyScaling = "spectrum",
) -> tuple[_float64_1d, _complex64_nd]: ...
@overload  # f80, input_onesided=True (default)
def istft(
    Zxx: _ToInexact80ND,
    fs: float = 1.0,
    window: _ToWindow = "hann_periodic",
    nperseg: int | None = None,
    noverlap: int | None = None,
    nfft: int | None = None,
    input_onesided: Literal[True] = True,
    boundary: bool = True,
    time_axis: int = -1,
    freq_axis: int = -2,
    scaling: _LegacyScaling = "spectrum",
) -> tuple[_float64_1d, _float80_nd]: ...
@overload  # c160, input_onesided=False
def istft(
    Zxx: _ToInexact80ND,
    fs: float = 1.0,
    window: _ToWindow = "hann_periodic",
    nperseg: int | None = None,
    noverlap: int | None = None,
    nfft: int | None = None,
    *,
    input_onesided: Literal[False],
    boundary: bool = True,
    time_axis: int = -1,
    freq_axis: int = -2,
    scaling: _LegacyScaling = "spectrum",
) -> tuple[_float64_1d, _complex160_nd]: ...
@overload  # fallback
def istft(
    Zxx: onp.ToComplexND,
    fs: float = 1.0,
    window: _ToWindow = "hann_periodic",
    nperseg: int | None = None,
    noverlap: int | None = None,
    nfft: int | None = None,
    *,
    input_onesided: bool,
    boundary: bool = True,
    time_axis: int = -1,
    freq_axis: int = -2,
    scaling: _LegacyScaling = "spectrum",
) -> tuple[_float64_1d, onp.Array]: ...

#
@overload  # f64
def coherence(  # type: ignore[overload-overlap]
    x: _ToInexact64ND,
    y: _ToInexact64ND,
    fs: float = 1.0,
    window: _ToWindow = "hann_periodic",
    nperseg: int | None = None,
    noverlap: int | None = None,
    nfft: int | None = None,
    detrend: _Detrend = "constant",
    axis: int = -1,
) -> tuple[_float64_1d, _float64_nd]: ...
@overload  # f32
def coherence(
    x: _ToInexact32ND,
    y: _ToInexact32ND,
    fs: float = 1.0,
    window: _ToWindow = "hann_periodic",
    nperseg: int | None = None,
    noverlap: int | None = None,
    nfft: int | None = None,
    detrend: _Detrend = "constant",
    axis: int = -1,
) -> tuple[_float64_1d, _float32_nd]: ...
@overload  # f80
def coherence(
    x: _ToInexact80ND,
    y: _ToInexact80ND,
    fs: float = 1.0,
    window: _ToWindow = "hann_periodic",
    nperseg: int | None = None,
    noverlap: int | None = None,
    nfft: int | None = None,
    detrend: _Detrend = "constant",
    axis: int = -1,
) -> tuple[_float64_1d, _float80_nd]: ...
@overload  # fallback
def coherence(
    x: onp.ToComplexND,
    y: onp.ToComplexND,
    fs: float = 1.0,
    window: _ToWindow = "hann_periodic",
    nperseg: int | None = None,
    noverlap: int | None = None,
    nfft: int | None = None,
    detrend: _Detrend = "constant",
    axis: int = -1,
) -> tuple[_float64_1d, onp.Array]: ...
