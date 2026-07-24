# This module is not meant for public use and will be removed in SciPy v2.0.0.
from typing_extensions import deprecated

__all__ = ["firls", "firwin", "firwin2", "firwin_2d", "kaiser_atten", "kaiser_beta", "kaiserord", "minimum_phase", "remez"]

@deprecated("will be removed in SciPy v2.0.0")
def firwin(
    numtaps: object,
    cutoff: object,
    *,
    width: object = None,
    window: object = "hamming",
    pass_zero: object = True,
    scale: object = True,
    fs: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def firwin_2d(
    hsize: object,
    window: object,
    *,
    fc: object | None = None,
    fs: object = 2,
    circular: object = False,
    pass_zero: object = True,
    scale: object = True,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def firwin2(
    numtaps: object,
    freq: object,
    gain: object,
    *,
    nfreqs: object = None,
    window: object = "hamming",
    antisymmetric: object = False,
    fs: object = None,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def firls(numtaps: object, bands: object, desired: object, *, weight: object = None, fs: object = None) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kaiser_beta(a: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kaiser_atten(numtaps: object, width: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kaiserord(ripple: object, width: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def minimum_phase(h: object, method: object = "homomorphic", n_fft: object = None, *, half: object = True) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def remez(
    numtaps: object,
    bands: object,
    desired: object,
    *,
    weight: object = None,
    type: object = "bandpass",
    maxiter: object = 25,
    grid_density: object = 16,
    fs: object = None,
) -> object: ...
