from typing import Final, Literal, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse import csr_array
from scipy.sparse._base import _spbase

###

type _Pair[T] = tuple[T, T]

type _Real = npc.integer | npc.floating
type _Int1D = onp.Array1D[np.int32]

type _ToGraph = onp.ToFloat2D | _spbase[_Real, tuple[int, int]]

###

DTYPE: Final[type[np.float64]] = ...
ITYPE: Final[type[np.int32]] = ...

@overload
def connected_components(
    csgraph: _ToGraph, directed: bool = True, connection: Literal["weak", "strong"] = "weak", return_labels: Literal[True] = True
) -> tuple[int, _Int1D]: ...
@overload
def connected_components(
    csgraph: _ToGraph, directed: bool = True, connection: Literal["weak", "strong"] = "weak", *, return_labels: Literal[False]
) -> int: ...

#
def breadth_first_tree(csgraph: _ToGraph, i_start: int, directed: bool = True) -> csr_array[np.float64, tuple[int, int]]: ...

#
def depth_first_tree(csgraph: _ToGraph, i_start: int, directed: bool = True) -> csr_array[np.float64, tuple[int, int]]: ...

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
