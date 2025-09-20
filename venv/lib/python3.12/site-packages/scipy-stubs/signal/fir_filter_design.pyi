# This module is not meant for public use and will be removed in SciPy v2.0.0.
from typing_extensions import deprecated

__all__ = ["firls", "firwin", "firwin2", "firwin_2d", "kaiser_atten", "kaiser_beta", "kaiserord", "minimum_phase", "remez"]

@deprecated("will be removed in SciPy v2.0.0")
def firwin(
    numtaps: object,
    cutoff: object,
    *,
    width: object = ...,
    window: object = ...,
    pass_zero: object = ...,
    scale: object = ...,
    fs: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def firwin_2d(
    hsize: object,
    window: object,
    *,
    fc: object | None = ...,
    fs: object = ...,
    circular: object = ...,
    pass_zero: object = ...,
    scale: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def firwin2(
    numtaps: object,
    freq: object,
    gain: object,
    *,
    nfreqs: object = ...,
    window: object = ...,
    antisymmetric: object = ...,
    fs: object = ...,
) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def firls(numtaps: object, bands: object, desired: object, *, weight: object = ..., fs: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kaiser_beta(a: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kaiser_atten(numtaps: object, width: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def kaiserord(ripple: object, width: object) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def minimum_phase(h: object, method: object = ..., n_fft: object = ..., *, half: object = ...) -> object: ...
@deprecated("will be removed in SciPy v2.0.0")
def remez(
    numtaps: object,
    bands: object,
    desired: object,
    *,
    weight: object = ...,
    type: object = ...,
    maxiter: object = ...,
    grid_density: object = ...,
    fs: object = ...,
) -> object: ...
