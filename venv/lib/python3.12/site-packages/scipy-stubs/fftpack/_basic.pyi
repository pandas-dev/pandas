from typing import Protocol, TypeAlias, type_check_only

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy._typing import AnyShape

__all__ = ["fft", "fft2", "fftn", "ifft", "ifft2", "ifftn", "irfft", "rfft"]

_ArrayReal: TypeAlias = onp.ArrayND[np.float32 | np.float64 | np.longdouble]  # no float16
_ArrayComplex: TypeAlias = onp.ArrayND[npc.complexfloating]

@type_check_only
class _OrderedIndex(op.CanIndex, Protocol):
    def __lt__(self, other: op.CanIndex, /) -> bool: ...
    def __le__(self, other: op.CanIndex, /) -> bool: ...

###

def fft(
    x: onp.ToComplexND, n: _OrderedIndex | None = None, axis: op.CanIndex = -1, overwrite_x: onp.ToBool = False
) -> _ArrayComplex: ...
def ifft(
    x: onp.ToComplexND, n: _OrderedIndex | None = None, axis: op.CanIndex = -1, overwrite_x: onp.ToBool = False
) -> _ArrayComplex: ...
def rfft(
    x: onp.ToFloatND, n: _OrderedIndex | None = None, axis: op.CanIndex = -1, overwrite_x: onp.ToBool = False
) -> _ArrayReal: ...
def irfft(
    x: onp.ToFloatND, n: _OrderedIndex | None = None, axis: op.CanIndex = -1, overwrite_x: onp.ToBool = False
) -> _ArrayReal: ...
def fftn(
    x: onp.ToComplexND, shape: AnyShape | None = None, axes: AnyShape | None = None, overwrite_x: onp.ToBool = False
) -> _ArrayComplex: ...
def ifftn(
    x: onp.ToComplexND, shape: AnyShape | None = None, axes: AnyShape | None = None, overwrite_x: onp.ToBool = False
) -> _ArrayComplex: ...
def fft2(
    x: onp.ToComplexND,
    shape: AnyShape | None = None,
    axes: tuple[op.CanIndex, op.CanIndex] = (-2, -1),
    overwrite_x: onp.ToBool = False,
) -> _ArrayComplex: ...
def ifft2(
    x: onp.ToComplexND,
    shape: AnyShape | None = None,
    axes: tuple[op.CanIndex, op.CanIndex] = (-2, -1),
    overwrite_x: onp.ToBool = False,
) -> _ArrayComplex: ...
