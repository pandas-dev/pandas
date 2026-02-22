# This module is not meant for public use and will be removed in SciPy v2.0.0.
from typing_extensions import deprecated

from ._filter_design import BadCoefficients as _BadCoefficients

__all__ = [
    "BadCoefficients",
    "band_stop_obj",
    "bessel",
    "besselap",
    "bilinear",
    "bilinear_zpk",
    "buttap",
    "butter",
    "buttord",
    "cheb1ap",
    "cheb1ord",
    "cheb2ap",
    "cheb2ord",
    "cheby1",
    "cheby2",
    "ellip",
    "ellipap",
    "ellipord",
    "findfreqs",
    "freqs",
    "freqs_zpk",
    "freqz",
    "freqz_sos",
    "freqz_zpk",
    "gammatone",
    "group_delay",
    "iircomb",
    "iirdesign",
    "iirfilter",
    "iirnotch",
    "iirpeak",
    "lp2bp",
    "lp2bp_zpk",
    "lp2bs",
    "lp2bs_zpk",
    "lp2hp",
    "lp2hp_zpk",
    "lp2lp",
    "lp2lp_zpk",
    "normalize",
    "sos2tf",
    "sos2zpk",
    "sosfreqz",
    "tf2sos",
    "tf2zpk",
    "zpk2sos",
    "zpk2tf",
]

@deprecated("will be removed in SciPy v2.0.0")
class BadCoefficients(_BadCoefficients): ...

@deprecated("will be removed in SciPy v2.0.0")
def findfreqs(num: object, den: object, N: object, kind: object = "ba") -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def freqs(b: object, a: object, worN: object = 200, plot: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def freqs_zpk(z: object, p: object, k: object, worN: object = 200) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def freqz(
    b: object,
    a: object = 1,
    worN: object = 512,
    whole: object = False,
    plot: object = None,
    fs: object = ...,
    include_nyquist: object = False,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def freqz_zpk(z: object, p: object, k: object, worN: object = 512, whole: object = False, fs: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def group_delay(system: object, w: object = 512, whole: object = False, fs: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def freqz_sos(sos: object, worN: object = 512, whole: object = False, fs: object = ...) -> object: ...

sosfreqz = freqz_sos  # pyright: ignore[reportDeprecated] # pyrefly: ignore[deprecated] # ty: ignore[deprecated]

@deprecated("will be removed in SciPy v2.0.0")
def tf2zpk(b: object, a: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def zpk2tf(z: object, p: object, k: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def tf2sos(b: object, a: object, pairing: object = None, *, analog: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def sos2tf(sos: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def sos2zpk(sos: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def zpk2sos(z: object, p: object, k: object, pairing: object = None, *, analog: object = False) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def normalize(b: object, a: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lp2lp(b: object, a: object, wo: object = 1.0) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lp2hp(b: object, a: object, wo: object = 1.0) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lp2bp(b: object, a: object, wo: object = 1.0, bw: object = 1.0) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lp2bs(b: object, a: object, wo: object = 1.0, bw: object = 1.0) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def bilinear(b: object, a: object, fs: object = 1.0) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def iirdesign(
    wp: object,
    ws: object,
    gpass: object,
    gstop: object,
    analog: object = False,
    ftype: object = "ellip",
    output: object = "ba",
    fs: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def iirfilter(
    N: object,
    Wn: object,
    rp: object = None,
    rs: object = None,
    btype: object = "band",
    analog: object = False,
    ftype: object = "butter",
    output: object = "ba",
    fs: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def bilinear_zpk(z: object, p: object, k: object, fs: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lp2lp_zpk(z: object, p: object, k: object, wo: object = 1.0) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lp2hp_zpk(z: object, p: object, k: object, wo: object = 1.0) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lp2bp_zpk(z: object, p: object, k: object, wo: object = 1.0, bw: object = 1.0) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def lp2bs_zpk(z: object, p: object, k: object, wo: object = 1.0, bw: object = 1.0) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def butter(
    N: object, Wn: object, btype: object = "low", analog: object = False, output: object = "ba", fs: object = None
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def cheby1(
    N: object, rp: object, Wn: object, btype: object = "low", analog: object = False, output: object = "ba", fs: object = None
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def cheby2(
    N: object, rs: object, Wn: object, btype: object = "low", analog: object = False, output: object = "ba", fs: object = None
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ellip(
    N: object,
    rp: object,
    rs: object,
    Wn: object,
    btype: object = "low",
    analog: object = False,
    output: object = "ba",
    fs: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def bessel(
    N: object,
    Wn: object,
    btype: object = "low",
    analog: object = False,
    output: object = "ba",
    norm: object = "phase",
    fs: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def band_stop_obj(
    wp: object, ind: object, passb: object, stopb: object, gpass: object, gstop: object, type: object
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def buttord(wp: object, ws: object, gpass: object, gstop: object, analog: object = False, fs: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def cheb1ord(wp: object, ws: object, gpass: object, gstop: object, analog: object = False, fs: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def cheb2ord(wp: object, ws: object, gpass: object, gstop: object, analog: object = False, fs: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ellipord(wp: object, ws: object, gpass: object, gstop: object, analog: object = False, fs: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def buttap(N: object, *, xp: object = None, device: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def cheb1ap(N: object, rp: object, *, xp: object = None, device: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def cheb2ap(N: object, rs: object, *, xp: object = None, device: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def ellipap(N: object, rp: object, rs: object, *, xp: object = None, device: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def besselap(N: object, norm: object = "phase", *, xp: object = None, device: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def iirnotch(w0: object, Q: object, fs: object = 2.0, *, xp: object = None, device: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def iirpeak(w0: object, Q: object, fs: object = 2.0, *, xp: object = None, device: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def iircomb(
    w0: object,
    Q: object,
    ftype: object = "notch",
    fs: object = 2.0,
    *,
    pass_zero: object = False,
    xp: object = None,
    device: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def gammatone(
    freq: object,
    ftype: object,
    order: object = None,
    numtaps: object = None,
    fs: object = None,
    *,
    xp: object = None,
    device: object = None,
) -> object: ...
