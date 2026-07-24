# This module is not meant for public use and will be removed in SciPy v2.0.0.
from typing import Final
from typing_extensions import deprecated

from numpy.fft import fftfreq

__all__ = ["fftfreq", "fftshift", "ifftshift", "next_fast_len", "rfftfreq"]

__MESSAGE: Final = "will be removed in SciPy v2.0.0"

@deprecated(__MESSAGE)
def fftshift(x: object, axes: object = None) -> object: ...
@deprecated(__MESSAGE)
def ifftshift(x: object, axes: object = None) -> object: ...
@deprecated(__MESSAGE)
def next_fast_len(target: object) -> object: ...
@deprecated(__MESSAGE)
def rfftfreq(n: object, d: object = 1.0) -> object: ...
