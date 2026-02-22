from typing import Literal, TypeAlias, overload

import numpy as np
import optype as op
import optype.numpy as onp

from scipy.sparse._base import _spbase

__all__ = ["norm"]

_Ord: TypeAlias = Literal["fro", 0, 1, 2, -1] | float
_Real: TypeAlias = np.int32 | np.int64 | np.float64

###

@overload  # no axis, two axes
def norm(x: _spbase, ord: _Ord | None = None, axis: tuple[op.CanIndex, op.CanIndex] | None = None) -> _Real: ...
@overload  # single axis (positional)
def norm(x: _spbase, ord: _Ord | None, axis: op.CanIndex) -> onp.Array1D[_Real]: ...
@overload  # single axis (keyword)
def norm(x: _spbase, ord: _Ord | None = None, *, axis: op.CanIndex) -> onp.Array1D[_Real]: ...
