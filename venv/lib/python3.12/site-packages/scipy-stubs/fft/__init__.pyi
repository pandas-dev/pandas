from contextlib import _GeneratorContextManager

from ._backend import register_backend, set_backend, set_global_backend, skip_backend
from ._basic import (
    fft,
    fft2,
    fftn,
    hfft,
    hfft2,
    hfftn,
    ifft,
    ifft2,
    ifftn,
    ihfft,
    ihfft2,
    ihfftn,
    irfft,
    irfft2,
    irfftn,
    rfft,
    rfft2,
    rfftn,
)
from ._fftlog import fht, fhtoffset, ifht
from ._helper import fftfreq, fftshift, ifftshift, next_fast_len, prev_fast_len, rfftfreq
from ._realtransforms import dct, dctn, dst, dstn, idct, idctn, idst, idstn

__all__ = [
    "dct",
    "dctn",
    "dst",
    "dstn",
    "fft",
    "fft2",
    "fftfreq",
    "fftn",
    "fftshift",
    "fht",
    "fhtoffset",
    "get_workers",
    "hfft",
    "hfft2",
    "hfftn",
    "idct",
    "idctn",
    "idst",
    "idstn",
    "ifft",
    "ifft2",
    "ifftn",
    "ifftshift",
    "ifht",
    "ihfft",
    "ihfft2",
    "ihfftn",
    "irfft",
    "irfft2",
    "irfftn",
    "next_fast_len",
    "prev_fast_len",
    "register_backend",
    "rfft",
    "rfft2",
    "rfftfreq",
    "rfftn",
    "set_backend",
    "set_global_backend",
    "set_workers",
    "skip_backend",
]

# originally defined in `scipy.fft._pocketfft.helper`
def set_workers(workers: int) -> _GeneratorContextManager[None]: ...
def get_workers() -> int: ...
