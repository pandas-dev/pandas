from collections.abc import Sequence
from typing import Any, Final, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse import csc_array, csc_matrix, csr_array, csr_matrix, lil_array, lil_matrix
from scipy.sparse._base import _spbase

###

type _Real = npc.integer | npc.floating

type _SparseGraph[RealT: _Real] = (
    csr_array[RealT] | csr_matrix[RealT]
    | csc_array[RealT] | csc_matrix[RealT]
    | lil_array[RealT] | lil_matrix[RealT]
)  # fmt: skip

type _ToGraph = onp.ToFloat2D | _spbase[_Real, tuple[int, int]]

###

DTYPE: Final[type[np.float64]] = ...
ITYPE: Final[type[np.int32]] = ...

#
def csgraph_from_masked(graph: onp.MArray2D[_Real]) -> csr_array[np.float64, tuple[int, int]]: ...

#
@overload
def csgraph_masked_from_dense[ScalarT: npc.number | np.bool](
    graph: onp.ArrayND[ScalarT], null_value: int | None = 0, nan_null: bool = True, infinity_null: bool = True, copy: bool = True
) -> onp.MArray2D[ScalarT]: ...
@overload
def csgraph_masked_from_dense(
    graph: Sequence[list[int]], null_value: int | None = 0, nan_null: bool = True, infinity_null: bool = True, copy: bool = True
) -> onp.MArray2D[np.int_]: ...
@overload
def csgraph_masked_from_dense(
    graph: Sequence[list[float]], null_value: int | None = 0, nan_null: bool = True, infinity_null: bool = True, copy: bool = True
) -> onp.MArray2D[np.float64]: ...
@overload
def csgraph_masked_from_dense(
    graph: onp.ToComplex2D, null_value: float | None = 0, nan_null: bool = True, infinity_null: bool = True, copy: bool = True
) -> onp.MArray2D[Any]: ...

#
def csgraph_from_dense(
    graph: onp.ToFloat2D, null_value: float | None = 0, nan_null: bool = True, infinity_null: bool = True
) -> csr_array[np.float64, tuple[int, int]]: ...

#
def csgraph_to_dense(csgraph: _SparseGraph[_Real], null_value: float | None = 0) -> onp.Array2D[np.float64]: ...

#
def csgraph_to_masked(csgraph: _SparseGraph[_Real]) -> onp.MArray2D[np.float64]: ...

#
def reconstruct_path(
    csgraph: _ToGraph, predecessors: onp.ToFloatND, directed: bool = True
) -> csr_array[np.float64, tuple[int, int]]: ...

#
def construct_dist_matrix(
    graph: _ToGraph, predecessors: onp.ToFloatND, directed: bool = True, null_value: float = ...
) -> onp.Array2D[np.float64]: ...
