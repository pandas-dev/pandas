from typing import Protocol, SupportsIndex, type_check_only

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy._typing import AnyShape

__all__ = ["fft", "fft2", "fftn", "ifft", "ifft2", "ifftn", "irfft", "rfft"]

type _ArrayReal = onp.ArrayND[np.float32 | np.float64 | np.longdouble]  # no float16
type _ArrayComplex = onp.ArrayND[npc.complexfloating]

@type_check_only
class _OrderedIndex(SupportsIndex, Protocol):
    def __lt__(self, other: SupportsIndex, /) -> bool: ...
    def __le__(self, other: SupportsIndex, /) -> bool: ...

###

def fft(
    x: onp.ToComplexND, n: _OrderedIndex | None = None, axis: SupportsIndex = -1, overwrite_x: bool = False
) -> _ArrayComplex: ...
def ifft(
    x: onp.ToComplexND, n: _OrderedIndex | None = None, axis: SupportsIndex = -1, overwrite_x: bool = False
) -> _ArrayComplex: ...
def rfft(x: onp.ToFloatND, n: _OrderedIndex | None = None, axis: SupportsIndex = -1, overwrite_x: bool = False) -> _ArrayReal: ...
def irfft(
    x: onp.ToFloatND, n: _OrderedIndex | None = None, axis: SupportsIndex = -1, overwrite_x: bool = False
) -> _ArrayReal: ...
def fftn(
    x: onp.ToComplexND, shape: AnyShape | None = None, axes: AnyShape | None = None, overwrite_x: bool = False
) -> _ArrayComplex: ...
def ifftn(
    x: onp.ToComplexND, shape: AnyShape | None = None, axes: AnyShape | None = None, overwrite_x: bool = False
) -> _ArrayComplex: ...
def fft2(
    x: onp.ToComplexND,
    shape: AnyShape | None = None,
    axes: tuple[SupportsIndex, SupportsIndex] = (-2, -1),
    overwrite_x: bool = False,
) -> _ArrayComplex: ...
def ifft2(
    x: onp.ToComplexND,
    shape: AnyShape | None = None,
    axes: tuple[SupportsIndex, SupportsIndex] = (-2, -1),
    overwrite_x: bool = False,
) -> _ArrayComplex: ...
