from collections.abc import Sequence
from typing import Any, TypeAlias, TypeVar, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["fht", "fhtoffset", "ifht"]

_FloatT = TypeVar("_FloatT", bound=np.float32 | np.float64 | npc.floating80)
_ShapeT = TypeVar("_ShapeT", bound=tuple[int, ...])

# workaround for a strange bug in pyright's overlapping overload detection with `numpy<2.1`
_WorkaroundForPyright: TypeAlias = tuple[int] | tuple[Any, ...]

###

@overload
def fht(
    a: onp.ArrayND[_FloatT, _ShapeT], dln: onp.ToFloat, mu: onp.ToFloat, offset: onp.ToFloat = 0.0, bias: onp.ToFloat = 0.0
) -> onp.ArrayND[_FloatT, _ShapeT]: ...
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
def ifht(
    A: onp.ArrayND[_FloatT, _ShapeT], dln: onp.ToFloat, mu: onp.ToFloat, offset: onp.ToFloat = 0.0, bias: onp.ToFloat = 0.0
) -> onp.ArrayND[_FloatT, _ShapeT]: ...
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
