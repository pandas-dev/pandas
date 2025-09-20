# This module is not meant for public use and will be removed in SciPy v2.0.0.
from typing_extensions import deprecated

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
    "nuttall",
    "parzen",
    "taylor",
    "triang",
    "tukey",
]

@deprecated("will be removed in SciPy v2.0.0")
def general_cosine(M: object, a: object, sym: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def boxcar(M: object, sym: object = ..., *, xp: object = None, device: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def triang(M: object, sym: object = ..., *, xp: object = None, device: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def parzen(M: object, sym: object = ..., *, xp: object = None, device: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def bohman(M: object, sym: object = ..., *, xp: object = None, device: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def blackman(M: object, sym: object = ..., *, xp: object = None, device: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def nuttall(M: object, sym: object = ..., *, xp: object = None, device: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def blackmanharris(M: object, sym: object = ..., *, xp: object = None, device: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def flattop(M: object, sym: object = ..., *, xp: object = None, device: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def bartlett(M: object, sym: object = ..., *, xp: object = None, device: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def hann(M: object, sym: object = ..., *, xp: object = None, device: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def tukey(M: object, alpha: object = ..., sym: object = ..., *, xp: object = None, device: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def barthann(M: object, sym: object = ..., *, xp: object = None, device: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def general_hamming(M: object, alpha: object, sym: object = ..., *, xp: object = None, device: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def hamming(M: object, sym: object = ..., *, xp: object = None, device: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kaiser(M: object, beta: object, sym: object = ..., *, xp: object = None, device: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def gaussian(M: object, std: object, sym: object = ..., *, xp: object = None, device: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def general_gaussian(
    M: object, p: object, sig: object, sym: object = ..., *, xp: object = None, device: object = None
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def chebwin(M: object, at: object, sym: object = ..., *, xp: object = None, device: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def cosine(M: object, sym: object = ..., *, xp: object = None, device: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def exponential(
    M: object, center: object = ..., tau: object = ..., sym: object = ..., *, xp: object = None, device: object = None
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def taylor(
    M: object,
    nbar: object = ...,
    sll: object = ...,
    norm: object = ...,
    sym: object = ...,
    *,
    xp: object = None,
    device: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def dpss(
    M: object,
    NW: object,
    Kmax: object = ...,
    sym: object = ...,
    norm: object = ...,
    return_ratios: object = ...,
    *,
    xp: object = None,
    device: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def get_window(window: object, Nx: object, fftbins: object = ..., *, xp: object = ..., device: object = ...) -> object: ...
