# This module is not meant for public use and will be removed in SciPy v2.0.0.
from typing_extensions import deprecated

from .windows.windows import get_window

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
def check_COLA(window: object, nperseg: object, noverlap: object, tol: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def check_NOLA(window: object, nperseg: object, noverlap: object, tol: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def coherence(
    x: object,
    y: object,
    fs: object = ...,
    window: object = ...,
    nperseg: object = ...,
    noverlap: object = ...,
    nfft: object = ...,
    detrend: object = ...,
    axis: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def csd(
    x: object,
    y: object,
    fs: object = ...,
    window: object = ...,
    nperseg: object = ...,
    noverlap: object = ...,
    nfft: object = ...,
    detrend: object = ...,
    return_onesided: object = ...,
    scaling: object = ...,
    axis: object = ...,
    average: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def istft(
    Zxx: object,
    fs: object = ...,
    window: object = ...,
    nperseg: object = ...,
    noverlap: object = ...,
    nfft: object = ...,
    input_onesided: object = ...,
    boundary: object = ...,
    time_axis: object = ...,
    freq_axis: object = ...,
    scaling: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lombscargle(
    x: object,
    y: object,
    freqs: object,
    precenter: object = ...,
    normalize: object = ...,
    *,
    weights: object | None = None,
    floating_mean: bool = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def periodogram(
    x: object,
    fs: object = ...,
    window: object = ...,
    nfft: object = ...,
    detrend: object = ...,
    return_onesided: object = ...,
    scaling: object = ...,
    axis: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def welch(
    x: object,
    fs: object = ...,
    window: object = ...,
    nperseg: object = ...,
    noverlap: object = ...,
    nfft: object = ...,
    detrend: object = ...,
    return_onesided: object = ...,
    scaling: object = ...,
    axis: object = ...,
    average: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def stft(
    x: object,
    fs: object = ...,
    window: object = ...,
    nperseg: object = ...,
    noverlap: object = ...,
    nfft: object = ...,
    detrend: object = ...,
    return_onesided: object = ...,
    boundary: object = ...,
    padded: object = ...,
    axis: object = ...,
    scaling: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def spectrogram(
    x: object,
    fs: object = ...,
    window: object = ...,
    nperseg: object = ...,
    noverlap: object = ...,
    nfft: object = ...,
    detrend: object = ...,
    return_onesided: object = ...,
    scaling: object = ...,
    axis: object = ...,
    mode: object = ...,
) -> object: ...
