from typing import Final

from ._basic import (
    fft as fft,
    fft2 as fft2,
    fftn as fftn,
    hfft as hfft,
    hfft2 as hfft2,
    hfftn as hfftn,
    ifft as ifft,
    ifft2 as ifft2,
    ifftn as ifftn,
    ihfft as ihfft,
    ihfft2 as ihfft2,
    ihfftn as ihfftn,
    irfft as irfft,
    irfft2 as irfft2,
    irfftn as irfftn,
    rfft as rfft,
    rfft2 as rfft2,
    rfftn as rfftn,
)

complex_funcs: Final = {"fft", "ifft", "fftn", "ifftn", "hfft", "irfft", "irfftn"}
