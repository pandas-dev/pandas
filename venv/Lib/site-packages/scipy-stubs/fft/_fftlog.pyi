from collections.abc import Sequence
from typing import Any, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["fht", "fhtoffset", "ifht"]

###

type _Float = np.float32 | np.float64 | npc.floating80

# workaround for a strange bug in pyright's overlapping overload detection with `numpy<2.1`
type _WorkaroundForPyright = tuple[int] | tuple[Any, ...]

###

@overload
def fht[FloatT: _Float, ShapeT: tuple[int, ...]](
    a: onp.ArrayND[FloatT, ShapeT], dln: onp.ToFloat, mu: onp.ToFloat, offset: onp.ToFloat = 0.0, bias: onp.ToFloat = 0.0
) -> onp.ArrayND[FloatT, ShapeT]: ...
@overload
def fht(
    a: Sequence[float], dln: onp.ToFloat, mu: onp.ToFloat, offset: onp.ToFloat = 0.0, bias: onp.ToFloat = 0.0
) -> onp.Array1D[np.float64]: ...
@overload
def fht(
    a: Sequence[Sequence[float]], dln: onp.ToFloat, mu: onp.ToFloat, offset: onp.ToFloat = 0.0, bias: onp.ToFloat = 0.0
) -> onp.Array2D[np.float64]: ...
@overload
def fht(
    a: Sequence[Sequence[Sequence[float]]], dln: onp.ToFloat, mu: onp.ToFloat, offset: onp.ToFloat = 0.0, bias: onp.ToFloat = 0.0
) -> onp.Array3D[np.float64]: ...
@overload
def fht(
    a: onp.ToFloatND, dln: onp.ToFloat, mu: onp.ToFloat, offset: onp.ToFloat = 0.0, bias: onp.ToFloat = 0.0
) -> onp.ArrayND[npc.floating, _WorkaroundForPyright]: ...

#
@overload
def ifht[FloatT: _Float, ShapeT: tuple[int, ...]](
    A: onp.ArrayND[FloatT, ShapeT], dln: onp.ToFloat, mu: onp.ToFloat, offset: onp.ToFloat = 0.0, bias: onp.ToFloat = 0.0
) -> onp.ArrayND[FloatT, ShapeT]: ...
@overload
def ifht(
    A: Sequence[float], dln: onp.ToFloat, mu: onp.ToFloat, offset: onp.ToFloat = 0.0, bias: onp.ToFloat = 0.0
) -> onp.Array1D[np.float64]: ...
@overload
def ifht(
    A: Sequence[Sequence[float]], dln: onp.ToFloat, mu: onp.ToFloat, offset: onp.ToFloat = 0.0, bias: onp.ToFloat = 0.0
) -> onp.Array2D[np.float64]: ...
@overload
def ifht(
    A: Sequence[Sequence[Sequence[float]]], dln: onp.ToFloat, mu: onp.ToFloat, offset: onp.ToFloat = 0.0, bias: onp.ToFloat = 0.0
) -> onp.Array3D[np.float64]: ...
@overload
def ifht(
    A: onp.ToFloatND, dln: onp.ToFloat, mu: onp.ToFloat, offset: onp.ToFloat = 0.0, bias: onp.ToFloat = 0.0
) -> onp.ArrayND[npc.floating, _WorkaroundForPyright]: ...

#
def fhtoffset(dln: onp.ToFloat, mu: onp.ToFloat, initial: onp.ToFloat = 0.0, bias: onp.ToFloat = 0.0) -> np.float64: ...
