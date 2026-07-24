# This module is not meant for public use and will be removed in SciPy v2.0.0.
from typing_extensions import deprecated

from .filter_design import cheby1  # pyrefly: ignore[deprecated] # ty: ignore[deprecated]
from .fir_filter_design import firwin  # pyrefly: ignore[deprecated] # ty: ignore[deprecated]
from .ltisys import dlti  # pyrefly: ignore[deprecated] # ty: ignore[deprecated]
from .windows.windows import get_window  # pyrefly: ignore[deprecated] # ty: ignore[deprecated]

__all__ = [
    "cheby1",
    "choose_conv_method",
    "convolve",
    "convolve2d",
    "correlate",
    "correlate2d",
    "correlation_lags",
    "decimate",
    "deconvolve",
    "detrend",
    "dlti",
    "fftconvolve",
    "filtfilt",
    "firwin",
    "get_window",
    "hilbert",
    "hilbert2",
    "invres",
    "invresz",
    "lfilter",
    "lfilter_zi",
    "lfiltic",
    "medfilt",
    "medfilt2d",
    "oaconvolve",
    "order_filter",
    "resample",
    "resample_poly",
    "residue",
    "residuez",
    "sosfilt",
    "sosfilt_zi",
    "sosfiltfilt",
    "unique_roots",
    "upfirdn",
    "vectorstrength",
    "wiener",
]

# ._upfirdn
def upfirdn(
    h: object, x: object, up: object = 1, down: object = 1, axis: object = -1, mode: object = "constant", cval: object = 0
) -> object: ...

# ._signaltools
@deprecated("will be removed in SciPy v2.0.0")
def correlate(in1: object, in2: object, mode: object = "full", method: object = "auto") -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def correlation_lags(in1_len: object, in2_len: object, mode: object = "full") -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def fftconvolve(in1: object, in2: object, mode: object = "full", axes: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def oaconvolve(in1: object, in2: object, mode: object = "full", axes: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def choose_conv_method(in1: object, in2: object, mode: object = "full", measure: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def convolve(in1: object, in2: object, mode: object = "full", method: object = "auto") -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def order_filter(a: object, domain: object, rank: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def medfilt(volume: object, kernel_size: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def wiener(im: object, mysize: object = None, noise: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def convolve2d(in1: object, in2: object, mode: object = "full", boundary: object = "fill", fillvalue: object = 0) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def correlate2d(in1: object, in2: object, mode: object = "full", boundary: object = "fill", fillvalue: object = 0) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def medfilt2d(input: object, kernel_size: object = 3) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lfilter(b: object, a: object, x: object, axis: object = -1, zi: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lfiltic(b: object, a: object, y: object, x: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def deconvolve(signal: object, divisor: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def hilbert(x: object, N: object = None, axis: object = -1) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def hilbert2(x: object, N: object = None, axes: object = (-2, -1)) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def unique_roots(p: object, tol: object = 1e-3, rtype: object = "min") -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def invres(r: object, p: object, k: object, tol: object = 1e-3, rtype: object = "avg") -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def residue(b: object, a: object, tol: object = 1e-3, rtype: object = "avg") -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def residuez(b: object, a: object, tol: object = 1e-3, rtype: object = "avg") -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def invresz(r: object, p: object, k: object, tol: object = 1e-3, rtype: object = "avg") -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def resample(
    x: object, num: object, t: object = None, axis: object = 0, window: object = None, domain: object = "time"
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def resample_poly(
    x: object,
    up: object,
    down: object,
    axis: object = 0,
    window: object = ("kaiser", 5.0),
    padtype: object = "constant",
    cval: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def vectorstrength(events: object, period: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def detrend(
    data: object, axis: object = -1, type: object = "linear", bp: object = 0, overwrite_data: object = False
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lfilter_zi(b: object, a: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def sosfilt_zi(sos: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def filtfilt(
    b: object,
    a: object,
    x: object,
    axis: object = -1,
    padtype: object = "odd",
    padlen: object = None,
    method: object = "pad",
    irlen: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def sosfilt(sos: object, x: object, axis: object = -1, zi: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def sosfiltfilt(sos: object, x: object, axis: object = -1, padtype: object = "odd", padlen: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def decimate(
    x: object, q: object, n: object = None, ftype: object = "iir", axis: object = -1, zero_phase: object = True
) -> object: ...
