# This module is not meant for public use and will be removed in SciPy v2.0.0.
from typing_extensions import deprecated

from .windows.windows import get_window  # pyrefly: ignore[deprecated] # ty: ignore[deprecated]

__all__ = [
    "check_COLA",
    "check_NOLA",
    "coherence",
    "csd",
    "get_window",
    "istft",
    "lombscargle",
    "periodogram",
    "spectrogram",
    "stft",
    "welch",
]

@deprecated("will be removed in SciPy v2.0.0")
def check_COLA(window: object, nperseg: object, noverlap: object, tol: object = 1e-10) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def check_NOLA(window: object, nperseg: object, noverlap: object, tol: object = 1e-10) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def coherence(
    x: object,
    y: object,
    fs: object = 1.0,
    window: object = "hann_periodic",
    nperseg: object = None,
    noverlap: object = None,
    nfft: object = None,
    detrend: object = "constant",
    axis: object = -1,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def csd(
    x: object,
    y: object,
    fs: object = 1.0,
    window: object = "hann_periodic",
    nperseg: object = None,
    noverlap: object = None,
    nfft: object = None,
    detrend: object = "constant",
    return_onesided: object = True,
    scaling: object = "density",
    axis: object = -1,
    average: object = "mean",
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def istft(
    Zxx: object,
    fs: object = 1.0,
    window: object = "hann_periodic",
    nperseg: object = None,
    noverlap: object = None,
    nfft: object = None,
    input_onesided: object = True,
    boundary: object = True,
    time_axis: object = -1,
    freq_axis: object = -2,
    scaling: object = "spectrum",
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lombscargle(
    x: object,
    y: object,
    freqs: object,
    *,
    precenter: object = ...,
    normalize: object = False,
    weights: object | None = None,
    floating_mean: bool = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def periodogram(
    x: object,
    fs: object = 1.0,
    window: object = "boxcar",
    nfft: object = None,
    detrend: object = "constant",
    return_onesided: object = True,
    scaling: object = "density",
    axis: object = -1,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def welch(
    x: object,
    fs: object = 1.0,
    window: object = "hann_periodic",
    nperseg: object = None,
    noverlap: object = None,
    nfft: object = None,
    detrend: object = "constant",
    return_onesided: object = True,
    scaling: object = "density",
    axis: object = -1,
    average: object = "mean",
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def stft(
    x: object,
    fs: object = 1.0,
    window: object = "hann_periodic",
    nperseg: object = 256,
    noverlap: object = None,
    nfft: object = None,
    detrend: object = False,
    return_onesided: object = True,
    boundary: object = "zeros",
    padded: object = True,
    axis: object = -1,
    scaling: object = "spectrum",
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def spectrogram(
    x: object,
    fs: object = 1.0,
    window: object = ("tukey_periodic", 0.25),
    nperseg: object = None,
    noverlap: object = None,
    nfft: object = None,
    detrend: object = "constant",
    return_onesided: object = True,
    scaling: object = "density",
    axis: object = -1,
    mode: object = "psd",
) -> object: ...
