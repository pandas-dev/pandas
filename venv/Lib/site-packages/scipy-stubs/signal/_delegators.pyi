from collections.abc import Callable
from types import ModuleType
from typing import Any, overload

import numpy as np

###

type _ = object  # ignored parameters

type _WindowSpec = str | tuple[object, ...] | Callable[..., object]

# these pretend that we can a module (`numpy` or `array_api_compat.np_compat`) in Python's type system
type _NumpyModule = ModuleType

# this pretends that we can express `array_namespace` in Python's type system
type _ArrayNamespace[ArrayT] = ModuleType

###

@overload
def _skip_if_lti[TupleT: tuple[object, ...]](arg: TupleT) -> TupleT: ...
@overload
def _skip_if_lti[T](arg: T) -> tuple[T] | Any: ...

#
@overload
def _skip_if_str_or_tuple(window: _WindowSpec) -> None: ...
@overload
def _skip_if_str_or_tuple[T](window: T) -> T: ...

#
@overload
def _skip_if_poly1d(arg: np.poly1d) -> None: ...
@overload
def _skip_if_poly1d[T](arg: T) -> T: ...

###

def abcd_normalize_signature[ArrayT](
    A: ArrayT | None = None, B: ArrayT | None = None, C: ArrayT | None = None, D: ArrayT | None = None
) -> _ArrayNamespace[ArrayT]: ...

#
def argrelextrema_signature[ArrayT](data: ArrayT, *args: _, **kwargs: _) -> _ArrayNamespace[ArrayT]: ...

argrelmax_signature = argrelextrema_signature
argrelmin_signature = argrelextrema_signature

def band_stop_obj_signature[ArrayT](
    wp: _, ind: _, passb: ArrayT, stopb: ArrayT, gpass: _, gstop: _, type: _
) -> _ArrayNamespace[ArrayT]: ...

#
def bessel_signature[ArrayT](N: _, Wn: ArrayT, *args: _, **kwds: _) -> _ArrayNamespace[ArrayT]: ...

butter_signature = bessel_signature

def cheby2_signature[ArrayT](N: _, rs: _, Wn: ArrayT, *args: _, **kwds: _) -> _ArrayNamespace[ArrayT]: ...
def cheby1_signature[ArrayT](N: _, rp: _, Wn: ArrayT, *args: _, **kwds: _) -> _ArrayNamespace[ArrayT]: ...
def ellip_signature[ArrayT](N: _, rp: _, rs: _, Wn: ArrayT, *args: _, **kwds: _) -> _ArrayNamespace[ArrayT]: ...

###

@overload
def besselap_signature(N: _, norm: _ = "phase", *, xp: None = None, device: None = None) -> _NumpyModule: ...
@overload
def besselap_signature[ModuleT: ModuleType](N: _, norm: _ = "phase", *, xp: ModuleT, device: _ = None) -> ModuleT: ...

#
@overload
def buttap_signature(N: _, *, xp: None = None, device: None = None) -> _NumpyModule: ...
@overload
def buttap_signature[ModuleT: ModuleType](N: _, *, xp: ModuleT, device: _ = None) -> ModuleT: ...

#
@overload
def cheb1ap_signature(N: _, rp: _, *, xp: None = None, device: None = None) -> _NumpyModule: ...
@overload
def cheb1ap_signature[ModuleT: ModuleType](N: _, rp: _, *, xp: ModuleT, device: _ = None) -> ModuleT: ...

#
@overload
def cheb2ap_signature(N: _, rs: _, *, xp: None = None, device: None = None) -> _NumpyModule: ...
@overload
def cheb2ap_signature[ModuleT: ModuleType](N: _, rs: _, *, xp: ModuleT, device: _ = None) -> ModuleT: ...

#
@overload
def ellipap_signature(N: _, rp: _, rs: _, *, xp: None = None, device: None = None) -> _NumpyModule: ...
@overload
def ellipap_signature[ModuleT: ModuleType](N: _, rp: _, rs: _, *, xp: ModuleT, device: _ = None) -> ModuleT: ...

#
def correlation_lags_signature(in1_len: _, in2_len: _, mode: _ = "full") -> _NumpyModule: ...

#
def czt_points_signature(m: _, w: _ = None, a: _ = 1 + 0j) -> _NumpyModule: ...

#
@overload
def gammatone_signature(
    freq: _, ftype: _, order: _ = None, numtaps: _ = None, fs: _ = None, *, xp: None = None, device: None = None
) -> _NumpyModule: ...
@overload
def gammatone_signature[ModuleT: ModuleType](
    freq: _, ftype: _, order: _ = None, numtaps: _ = None, fs: _ = None, *, xp: ModuleT, device: _ = None
) -> ModuleT: ...

#
@overload
def iircomb_signature(
    w0: _, Q: _, ftype: _ = "notch", fs: _ = 2.0, *, pass_zero: _ = False, xp: None = None, device: None = None
) -> _NumpyModule: ...
@overload
def iircomb_signature[ModuleT: ModuleType](
    w0: _, Q: _, ftype: _ = "notch", fs: _ = 2.0, *, pass_zero: _ = False, xp: ModuleT, device: _ = None
) -> ModuleT: ...

#
@overload
def iirnotch_signature(w0: _, Q: _, fs: _ = 2.0, *, xp: None = None, device: None = None) -> _NumpyModule: ...
@overload
def iirnotch_signature[ModuleT: ModuleType](w0: _, Q: _, fs: _ = 2.0, *, xp: ModuleT, device: _ = None) -> ModuleT: ...

iirpeak_signature = iirnotch_signature

@overload
def savgol_coeffs_signature(
    window_length: _,
    polyorder: _,
    deriv: _ = 0,
    delta: _ = 1.0,
    pos: _ = None,
    use: _ = "conv",
    *,
    xp: None = None,
    device: None = None,
) -> _NumpyModule: ...
@overload
def savgol_coeffs_signature[ModuleT: ModuleType](
    window_length: _, polyorder: _, deriv: _ = 0, delta: _ = 1.0, pos: _ = None, use: _ = "conv", *, xp: ModuleT, device: _ = None
) -> ModuleT: ...

#
def unit_impulse_signature(shape: _, idx: _ = None, dtype: _ = float) -> _NumpyModule: ...  # noqa: PYI011

#
def buttord_signature[ArrayT](
    wp: ArrayT, ws: ArrayT, gpass: _, gstop: _, analog: _ = False, fs: _ = None
) -> _ArrayNamespace[ArrayT]: ...

cheb1ord_signature = buttord_signature
cheb2ord_signature = buttord_signature
ellipord_signature = buttord_signature

#
def kaiser_atten_signature(numtaps: _, width: _) -> _NumpyModule: ...
def kaiser_beta_signature(a: _) -> _NumpyModule: ...
def kaiserord_signature(ripple: _, width: _) -> _NumpyModule: ...

#
@overload
def get_window_signature(window: _, Nx: _, fftbins: _ = True, *, xp: None = None, device: None = None) -> _NumpyModule: ...
@overload
def get_window_signature[ModuleT: ModuleType](
    window: _, Nx: _, fftbins: _ = True, *, xp: ModuleT, device: _ = None
) -> ModuleT: ...

#
def bode_signature[ArrayT](
    system: ArrayT | tuple[ArrayT, ...], w: ArrayT | None = None, n: _ = 100
) -> _ArrayNamespace[ArrayT]: ...

dbode_signature = bode_signature

def freqresp_signature[ArrayT](
    system: ArrayT | tuple[ArrayT, ...], w: ArrayT | None = None, n: _ = 10_000
) -> _ArrayNamespace[ArrayT]: ...

dfreqresp_signature = freqresp_signature

def impulse_signature[ArrayT](
    system: ArrayT | tuple[ArrayT, ...], X0: ArrayT | None = None, T: ArrayT | None = None, N: _ = None
) -> _ArrayNamespace[ArrayT]: ...

#
def dimpulse_signature[ArrayT](
    system: ArrayT | tuple[ArrayT, ...], x0: ArrayT | None = None, t: ArrayT | None = None, n: _ = None
) -> _ArrayNamespace[ArrayT]: ...

#
def lsim_signature[ArrayT](
    system: ArrayT | tuple[ArrayT, ...], U: ArrayT, T: ArrayT, X0: ArrayT | None = None, interp: _ = True
) -> _ArrayNamespace[ArrayT]: ...

#
def dlsim_signature[ArrayT](
    system: ArrayT | tuple[ArrayT, ...], u: ArrayT, t: ArrayT | None = None, x0: ArrayT | None = None
) -> _ArrayNamespace[ArrayT]: ...

#
def step_signature[ArrayT](
    system: ArrayT | tuple[ArrayT, ...], X0: ArrayT | None = None, T: ArrayT | None = None, N: _ = None
) -> _ArrayNamespace[ArrayT]: ...

#
def dstep_signature[ArrayT](
    system: ArrayT | tuple[ArrayT, ...], x0: ArrayT | None = None, t: ArrayT | None = None, n: _ = None
) -> _ArrayNamespace[ArrayT]: ...

#
def cont2discrete_signature[ArrayT](
    system: ArrayT | tuple[ArrayT, ...], dt: _, method: _ = "zoh", alpha: _ = None
) -> _ArrayNamespace[ArrayT]: ...

#
def bilinear_signature[ArrayT](b: ArrayT, a: ArrayT, fs: _ = 1.0) -> _ArrayNamespace[ArrayT]: ...
def bilinear_zpk_signature[ArrayT](z: ArrayT, p: ArrayT, k: _, fs: _) -> _ArrayNamespace[ArrayT]: ...

#
def chirp_signature[ArrayT](t: ArrayT, *args: _, **kwds: _) -> _ArrayNamespace[ArrayT]: ...
def choose_conv_method_signature[ArrayT](in1: ArrayT, in2: ArrayT, *args: _, **kwds: _) -> _ArrayNamespace[ArrayT]: ...
def convolve_signature[ArrayT](in1: ArrayT, in2: ArrayT, *args: _, **kwds: _) -> _ArrayNamespace[ArrayT]: ...

fftconvolve_signature = convolve_signature
oaconvolve_signature = convolve_signature
correlate_signature = convolve_signature
correlate_signature = convolve_signature
convolve2d_signature = convolve_signature
correlate2d_signature = convolve_signature

def coherence_signature[ArrayT](
    x: ArrayT, y: ArrayT, fs: _ = 1.0, window: ArrayT | _WindowSpec = "hann_periodic", *args: _, **kwds: _
) -> _ArrayNamespace[ArrayT]: ...

#
def csd_signature[ArrayT](
    x: ArrayT, y: ArrayT, fs: _ = 1.0, window: ArrayT | _WindowSpec = "hann_periodic", *args: _, **kwds: _
) -> _ArrayNamespace[ArrayT]: ...

#
def periodogram_signature[ArrayT](
    x: ArrayT, fs: _ = 1.0, window: ArrayT | _WindowSpec = "boxcar", *args: _, **kwds: _
) -> _ArrayNamespace[ArrayT]: ...

#
def welch_signature[ArrayT](
    x: ArrayT, fs: _ = 1.0, window: ArrayT | _WindowSpec = "hann_periodic", *args: _, **kwds: _
) -> _ArrayNamespace[ArrayT]: ...

#
def spectrogram_signature[ArrayT](
    x: ArrayT, fs: _ = 1.0, window: ArrayT | _WindowSpec = ("tukey_periodic", 0.25), *args: _, **kwds: _
) -> _ArrayNamespace[ArrayT]: ...

#
def stft_signature[ArrayT](
    x: ArrayT, fs: _ = 1.0, window: ArrayT | _WindowSpec = "hann_periodic", *args: _, **kwds: _
) -> _ArrayNamespace[ArrayT]: ...

#
def istft_signature[ArrayT](
    Zxx: ArrayT, fs: _ = 1.0, window: ArrayT | _WindowSpec = "hann_periodic", *args: _, **kwds: _
) -> _ArrayNamespace[ArrayT]: ...

#
def resample_signature[ArrayT](
    x: ArrayT, num: _, t: ArrayT | None = None, axis: _ = 0, window: ArrayT | _WindowSpec | None = None, domain: _ = "time"
) -> _ArrayNamespace[ArrayT]: ...

#
def resample_poly_signature[ArrayT](
    x: ArrayT, up: _, down: _, axis: _ = 0, window: ArrayT | _WindowSpec = ("kaiser", 5.0), *args: _, **kwds: _
) -> _ArrayNamespace[ArrayT]: ...

#
def check_COLA_signature[ArrayT](
    window: ArrayT | _WindowSpec, nperseg: _, noverlap: _, tol: _ = 1e-10
) -> _ArrayNamespace[ArrayT]: ...

#
def check_NOLA_signature[ArrayT](
    window: ArrayT | _WindowSpec, nperseg: _, noverlap: _, tol: _ = 1e-10
) -> _ArrayNamespace[ArrayT]: ...

#
def czt_signature[ArrayT](x: ArrayT, *args: _, **kwds: _) -> _ArrayNamespace[ArrayT]: ...

decimate_signature = czt_signature
gauss_spline_signature = czt_signature

def deconvolve_signature[ArrayT](signal: ArrayT, divisor: ArrayT) -> _ArrayNamespace[ArrayT]: ...

#
def detrend_signature[ArrayT](
    data: ArrayT, axis: _ = 1, type: _ = "linear", bp: ArrayT | _ = 0, *args: _, **kwds: _
) -> _ArrayNamespace[ArrayT]: ...

#
def filtfilt_signature[ArrayT](b: ArrayT, a: ArrayT, x: ArrayT, *args: _, **kwds: _) -> _ArrayNamespace[ArrayT]: ...

#
def lfilter_signature[ArrayT](
    b: ArrayT, a: ArrayT, x: ArrayT, axis: _ = -1, zi: ArrayT | None = None
) -> _ArrayNamespace[ArrayT]: ...

#
def envelope_signature[ArrayT](z: ArrayT, *args: _, **kwds: _) -> _ArrayNamespace[ArrayT]: ...

#
def find_peaks_signature(
    x: _,
    height: _ = None,
    threshold: _ = None,
    distance: _ = None,
    prominence: _ = None,
    width: _ = None,
    wlen: _ = None,
    rel_height: _ = 0.5,
    plateau_size: _ = None,
) -> _NumpyModule: ...  # np_compat

#
def find_peaks_cwt_signature[ArrayT](
    vector: ArrayT, widths: ArrayT, wavelet: _ = None, max_distances: ArrayT | None = None, *args: _, **kwds: _
) -> _ArrayNamespace[ArrayT]: ...

#
def findfreqs_signature[ArrayT](num: ArrayT, den: ArrayT, N: _, kind: _ = "ba") -> _ArrayNamespace[ArrayT]: ...

#
def firls_signature[ArrayT](
    numtaps: _, bands: ArrayT, desired: ArrayT, *, weight: ArrayT | None = None, fs: _ = None
) -> _ArrayNamespace[ArrayT]: ...

#
@overload
def firwin_signature(numtaps: _, cutoff: float, *args: _, **kwds: _) -> _NumpyModule: ...  # np_compat
@overload
def firwin_signature[ArrayT](numtaps: _, cutoff: ArrayT, *args: _, **kwds: _) -> _ArrayNamespace[ArrayT]: ...

#
def firwin2_signature[ArrayT](numtaps: _, freq: ArrayT, gain: ArrayT, *args: _, **kwds: _) -> _ArrayNamespace[ArrayT]: ...

#
def freqs_zpk_signature[ArrayT](
    z: ArrayT, p: ArrayT, k: _, worN: ArrayT | int = 200, *args: _, **kwds: _
) -> _ArrayNamespace[ArrayT]: ...

freqz_zpk_signature = freqs_zpk_signature

def freqs_signature[ArrayT](b: ArrayT, a: ArrayT, worN: ArrayT | int = 200, *args: _, **kwds: _) -> _ArrayNamespace[ArrayT]: ...
def freqz_signature[ArrayT](
    b: ArrayT, a: ArrayT | int = 1, worN: ArrayT | int = 512, *args: _, **kwds: _
) -> _ArrayNamespace[ArrayT]: ...
def freqz_sos_signature[ArrayT](sos: ArrayT, worN: ArrayT | int = 512, *args: _, **kwds: _) -> _ArrayNamespace[ArrayT]: ...

sosfreqz_signature = freqz_sos_signature

@overload
def gausspulse_signature(t: str, *args: _, **kwds: _) -> _NumpyModule: ...
@overload
def gausspulse_signature[ArrayT](t: ArrayT, *args: _, **kwds: _) -> _ArrayNamespace[ArrayT]: ...

#
def group_delay_signature[ArrayT](
    system: tuple[ArrayT, ArrayT], w: ArrayT | _ = 512, whole: _ = False, fs: _ = 6.283185307179586
) -> _ArrayNamespace[ArrayT]: ...

#
def hilbert_signature[ArrayT](x: ArrayT, *args: _, **kwds: _) -> _ArrayNamespace[ArrayT]: ...

hilbert2_signature = hilbert_signature

def iirdesign_signature[ArrayT](wp: ArrayT, ws: ArrayT, *args: _, **kwds: _) -> _ArrayNamespace[ArrayT]: ...
def iirfilter_signature[ArrayT](N: _, Wn: ArrayT, *args: _, **kwds: _) -> _ArrayNamespace[ArrayT]: ...
def invres_signature[ArrayT](r: ArrayT, p: ArrayT, k: ArrayT, tol: _ = 0.001, rtype: _ = "avg") -> _ArrayNamespace[ArrayT]: ...

invresz_signature = invres_signature

def lfilter_zi_signature[ArrayT](b: ArrayT, a: ArrayT) -> _ArrayNamespace[ArrayT]: ...
def sosfilt_zi_signature[ArrayT](sos: ArrayT) -> _ArrayNamespace[ArrayT]: ...

#
def remez_signature[ArrayT](
    numtaps: _, bands: ArrayT, desired: ArrayT, *, weight: ArrayT | None = None, **kwds: _
) -> _ArrayNamespace[ArrayT]: ...

#
def lfiltic_signature[ArrayT](b: ArrayT, a: ArrayT, y: ArrayT, x: ArrayT | None = None) -> _ArrayNamespace[ArrayT]: ...

#
def lombscargle_signature[ArrayT](
    x: ArrayT,
    y: ArrayT,
    freqs: ArrayT,
    precenter: _ = False,
    normalize: _ = False,
    *,
    weights: ArrayT | None = None,
    floating_mean: _ = False,
) -> _ArrayNamespace[ArrayT]: ...

#
def lp2bp_signature[ArrayT](b: ArrayT, a: ArrayT, *args: _, **kwds: _) -> _ArrayNamespace[ArrayT]: ...

lp2bs_signature = lp2bp_signature
lp2hp_signature = lp2bp_signature
lp2lp_signature = lp2bp_signature

tf2zpk_signature = lp2bp_signature
tf2sos_signature = lp2bp_signature

normalize_signature = lp2bp_signature
residue_signature = lp2bp_signature
residuez_signature = residue_signature

def lp2bp_zpk_signature[ArrayT](z: ArrayT, p: ArrayT, k: _, *args: _, **kwds: _) -> _ArrayNamespace[ArrayT]: ...

lp2bs_zpk_signature = lp2bp_zpk_signature
lp2hp_zpk_signature = lp2bs_zpk_signature
lp2lp_zpk_signature = lp2bs_zpk_signature

def zpk2sos_signature[ArrayT](z: ArrayT, p: ArrayT, k: _, *args: _, **kwds: _) -> _ArrayNamespace[ArrayT]: ...

zpk2ss_signature = zpk2sos_signature
zpk2tf_signature = zpk2sos_signature

def max_len_seq_signature[ArrayT](
    nbits: _, state: ArrayT | None = None, length: _ = None, taps: ArrayT | None = None
) -> _ArrayNamespace[ArrayT]: ...

#
def medfilt_signature[ArrayT](volume: ArrayT, kernel_size: _ = None) -> _ArrayNamespace[ArrayT]: ...
def medfilt2d_signature[ArrayT](input: ArrayT, kernel_size: _ = 3) -> _ArrayNamespace[ArrayT]: ...

#
def minimum_phase_signature[ArrayT](h: ArrayT, *args: _, **kwds: _) -> _ArrayNamespace[ArrayT]: ...

#
def order_filter_signature[ArrayT](a: ArrayT, domain: ArrayT, rank: _) -> _ArrayNamespace[ArrayT]: ...

#
def peak_prominences_signature[ArrayT](x: ArrayT, peaks: ArrayT, *args: _, **kwds: _) -> _ArrayNamespace[ArrayT]: ...

peak_widths_signature = peak_prominences_signature

def place_poles_signature[ArrayT](
    A: ArrayT, B: ArrayT, poles: ArrayT, method: _ = "YT", rtol: _ = 0.001, maxiter: _ = 30
) -> _ArrayNamespace[ArrayT]: ...

#
def savgol_filter_signature[ArrayT](x: ArrayT, *args: _, **kwds: _) -> _ArrayNamespace[ArrayT]: ...
def sawtooth_signature[ArrayT](t: ArrayT, width: _ = 1) -> _ArrayNamespace[ArrayT]: ...
def sepfir2d_signature[ArrayT](input: ArrayT, hrow: ArrayT, hcol: ArrayT) -> _ArrayNamespace[ArrayT]: ...
def sos2tf_signature[ArrayT](sos: ArrayT) -> _ArrayNamespace[ArrayT]: ...

sos2zpk_signature = sos2tf_signature

def sosfilt_signature[ArrayT](sos: ArrayT, x: ArrayT, axis: _ = -1, zi: ArrayT | None = None) -> _ArrayNamespace[ArrayT]: ...
def sosfiltfilt_signature[ArrayT](sos: ArrayT, x: ArrayT, *args: _, **kwds: _) -> _ArrayNamespace[ArrayT]: ...
def spline_filter_signature[ArrayT](Iin: ArrayT, lmbda: _ = 5.0) -> _ArrayNamespace[ArrayT]: ...
def square_signature[ArrayT](t: ArrayT, duty: _ = 0.5) -> _ArrayNamespace[ArrayT]: ...
def ss2tf_signature[ArrayT](A: ArrayT, B: ArrayT, C: ArrayT, D: ArrayT, input: _ = 0) -> _ArrayNamespace[ArrayT]: ...

ss2zpk_signature = ss2tf_signature

def sweep_poly_signature[ArrayT](t: ArrayT, poly: ArrayT | np.poly1d, phi: _ = 0) -> _ArrayNamespace[ArrayT]: ...

#
def symiirorder1_signature[ArrayT](signal: ArrayT, c0: _, z1: _, precision: _ = -1.0) -> _ArrayNamespace[ArrayT]: ...
def symiirorder2_signature[ArrayT](input: ArrayT, r: _, omega: _, precision: _ = -1.0) -> _ArrayNamespace[ArrayT]: ...
def cspline1d_signature[ArrayT](signal: ArrayT, *args: _, **kwds: _) -> _ArrayNamespace[ArrayT]: ...

qspline1d_signature = cspline1d_signature
cspline2d_signature = cspline1d_signature
qspline2d_signature = qspline1d_signature

def cspline1d_eval_signature[ArrayT](cj: ArrayT, newx: ArrayT, *args: _, **kwds: _) -> _ArrayNamespace[ArrayT]: ...

qspline1d_eval_signature = cspline1d_eval_signature

def tf2ss_signature[ArrayT](num: ArrayT, den: ArrayT) -> _ArrayNamespace[ArrayT]: ...
def unique_roots_signature[ArrayT](p: ArrayT, tol: _ = 0.001, rtype: _ = "min") -> _ArrayNamespace[ArrayT]: ...
def upfirdn_signature[ArrayT](
    h: ArrayT, x: ArrayT, up: _ = 1, down: _ = 1, axis: _ = -1, mode: _ = "constant", cval: _ = 0
) -> _ArrayNamespace[ArrayT]: ...
def vectorstrength_signature[ArrayT](events: ArrayT, period: ArrayT) -> _ArrayNamespace[ArrayT]: ...
def wiener_signature[ArrayT](im: ArrayT, mysize: _ = None, noise: _ = None) -> _ArrayNamespace[ArrayT]: ...
def zoom_fft_signature[ArrayT](
    x: ArrayT, fn: ArrayT, m: _ = None, *, fs: _ = 2, endpoint: _ = False, axis: _ = -1
) -> _ArrayNamespace[ArrayT]: ...
