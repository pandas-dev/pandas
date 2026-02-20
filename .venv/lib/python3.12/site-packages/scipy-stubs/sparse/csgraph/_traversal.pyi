from typing import Final, Literal, TypeAlias, TypeVar, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse import csr_array
from scipy.sparse._base import _spbase

_T = TypeVar("_T")
_Pair: TypeAlias = tuple[_T, _T]

_Real: TypeAlias = npc.integer | npc.floating
_Int1D: TypeAlias = onp.Array1D[np.int32]

_ToGraph: TypeAlias = onp.ToFloat2D | _spbase[_Real, tuple[int, int]]

_RealT = TypeVar("_RealT", bound=_Real)
_Graph: TypeAlias = onp.CanArrayND[_RealT] | _spbase[_RealT, tuple[int, int]]

###

DTYPE: Final[type[np.float64]] = ...
ITYPE: Final[type[np.int32]] = ...

def connected_components(
    csgraph: _ToGraph, directed: bool = True, connection: Literal["weak", "strong"] = "weak", return_labels: bool = True
) -> tuple[int, _Int1D]: ...

#
@overload
def breadth_first_tree(csgraph: _Graph[_RealT], i_start: int, directed: bool = True) -> csr_array[_RealT, tuple[int, int]]: ...
@overload
def breadth_first_tree(csgraph: _ToGraph, i_start: int, directed: bool = True) -> csr_array[_Real, tuple[int, int]]: ...

#
@overload
def depth_first_tree(csgraph: _Graph[_RealT], i_start: int, directed: bool = True) -> csr_array[_RealT, tuple[int, int]]: ...
@overload
def depth_first_tree(csgraph: _ToGraph, i_start: int, directed: bool = True) -> csr_array[_Real, tuple[int, int]]: ...

#
@overload
def breadth_first_order(
    csgraph: _ToGraph, i_start: int, directed: bool = True, return_predecessors: onp.ToTrue = True
) -> _Pair[_Int1D]: ...
@overload
def breadth_first_order(csgraph: _ToGraph, i_start: int, directed: bool, return_predecessors: onp.ToFalse) -> _Int1D: ...
@overload
def breadth_first_order(
    csgraph: _ToGraph, i_start: int, directed: bool = True, *, return_predecessors: onp.ToFalse
) -> _Int1D: ...

#
@overload
def depth_first_order(
    csgraph: _ToGraph, i_start: int, directed: bool = True, return_predecessors: onp.ToTrue = True
) -> _Pair[_Int1D]: ...
@overload
def depth_first_order(csgraph: _ToGraph, i_start: int, directed: bool, return_predecessors: onp.ToFalse) -> _Int1D: ...
@overload
def depth_first_order(csgraph: _ToGraph, i_start: int, directed: bool = True, *, return_predecessors: onp.ToFalse) -> _Int1D: ...
