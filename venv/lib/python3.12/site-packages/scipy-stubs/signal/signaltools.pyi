# This module is not meant for public use and will be removed in SciPy v2.0.0.
from typing_extensions import deprecated

from .filter_design import cheby1
from .fir_filter_design import firwin
from .ltisys import dlti
from .windows.windows import get_window

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
    h: object, x: object, up: object = ..., down: object = ..., axis: object = ..., mode: object = ..., cval: object = ...
) -> object: ...

# ._signaltools
@deprecated("will be removed in SciPy v2.0.0")
def correlate(in1: object, in2: object, mode: object = ..., method: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def correlation_lags(in1_len: object, in2_len: object, mode: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def fftconvolve(in1: object, in2: object, mode: object = ..., axes: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def oaconvolve(in1: object, in2: object, mode: object = ..., axes: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def choose_conv_method(in1: object, in2: object, mode: object = ..., measure: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def convolve(in1: object, in2: object, mode: object = ..., method: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def order_filter(a: object, domain: object, rank: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def medfilt(volume: object, kernel_size: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def wiener(im: object, mysize: object = ..., noise: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def convolve2d(in1: object, in2: object, mode: object = ..., boundary: object = ..., fillvalue: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def correlate2d(in1: object, in2: object, mode: object = ..., boundary: object = ..., fillvalue: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def medfilt2d(input: object, kernel_size: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lfilter(b: object, a: object, x: object, axis: object = ..., zi: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lfiltic(b: object, a: object, y: object, x: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def deconvolve(signal: object, divisor: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def hilbert(x: object, N: object = ..., axis: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def hilbert2(x: object, N: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def unique_roots(p: object, tol: object = ..., rtype: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def invres(r: object, p: object, k: object, tol: object = ..., rtype: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def residue(b: object, a: object, tol: object = ..., rtype: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def residuez(b: object, a: object, tol: object = ..., rtype: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def invresz(r: object, p: object, k: object, tol: object = ..., rtype: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def resample(
    x: object, num: object, t: object = ..., axis: object = ..., window: object = ..., domain: object = ...
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def resample_poly(
    x: object, up: object, down: object, axis: object = ..., window: object = ..., padtype: object = ..., cval: object = ...
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def vectorstrength(events: object, period: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def detrend(data: object, axis: object = ..., type: object = ..., bp: object = ..., overwrite_data: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lfilter_zi(b: object, a: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def sosfilt_zi(sos: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def filtfilt(
    b: object,
    a: object,
    x: object,
    axis: object = ...,
    padtype: object = ...,
    padlen: object = ...,
    method: object = ...,
    irlen: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def sosfilt(sos: object, x: object, axis: object = ..., zi: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def sosfiltfilt(sos: object, x: object, axis: object = ..., padtype: object = ..., padlen: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def decimate(
    x: object, q: object, n: object = ..., ftype: object = ..., axis: object = ..., zero_phase: object = ...
) -> object: ...
