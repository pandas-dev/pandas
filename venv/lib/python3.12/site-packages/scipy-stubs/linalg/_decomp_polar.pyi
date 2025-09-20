from typing import Literal, TypeAlias, TypeVar, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = ["polar"]

_T = TypeVar("_T")
_Tuple2: TypeAlias = tuple[_T, _T]
_Side: TypeAlias = Literal["left", "right"]

###

@overload
def polar(a: onp.ToIntStrict2D | onp.ToJustFloat64Strict2D, side: _Side = "right") -> _Tuple2[onp.Array2D[np.float64]]: ...
@overload
def polar(a: onp.ToIntND | onp.ToJustFloat64_ND, side: _Side = "right") -> _Tuple2[onp.ArrayND[np.float64]]: ...
@overload
def polar(a: onp.ToFloatStrict2D, side: _Side = "right") -> _Tuple2[onp.Array2D[npc.floating]]: ...
@overload
def polar(a: onp.ToFloatND, side: _Side = "right") -> _Tuple2[onp.ArrayND[npc.floating]]: ...
@overload
def polar(a: onp.ToJustComplexStrict2D, side: _Side = "right") -> _Tuple2[onp.Array2D[npc.complexfloating]]: ...
@overload
def polar(a: onp.ToJustComplexND, side: _Side = "right") -> _Tuple2[onp.ArrayND[npc.complexfloating]]: ...
@overload
def polar(a: onp.ToComplexND, side: _Side = "right") -> _Tuple2[onp.ArrayND[npc.inexact]]: ...
