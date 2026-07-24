from typing import Literal, SupportsIndex, overload

import numpy as np
import optype.numpy as onp

from scipy.sparse._base import _spbase

__all__ = ["norm"]

###

type _Ord = Literal["fro", 0, 1, 2, -1] | float
type _Real = np.int32 | np.int64 | np.float64

###

@overload  # no axis, two axes
def norm(x: _spbase, ord: _Ord | None = None, axis: tuple[SupportsIndex, SupportsIndex] | None = None) -> _Real: ...
@overload  # single axis (positional)
def norm(x: _spbase, ord: _Ord | None, axis: SupportsIndex) -> onp.Array1D[_Real]: ...
@overload  # single axis (keyword)
def norm(x: _spbase, ord: _Ord | None = None, *, axis: SupportsIndex) -> onp.Array1D[_Real]: ...
