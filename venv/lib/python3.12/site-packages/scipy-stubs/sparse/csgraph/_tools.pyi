from typing import Final, TypeAlias, overload
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse import csc_array, csc_matrix, csr_array, csr_matrix, lil_array, lil_matrix
from scipy.sparse._base import _spbase

_Real: TypeAlias = npc.integer | npc.floating
_RealT = TypeVar("_RealT", bound=_Real)

_SparseGraph: TypeAlias = (
    csr_array[_RealT] | csr_matrix[_RealT]
    | csc_array[_RealT] | csc_matrix[_RealT]
    | lil_array[_RealT] | lil_matrix[_RealT]
)  # fmt: skip

_ToGraph: TypeAlias = onp.ToFloat2D | _spbase[_Real, tuple[int, int]]
_Graph: TypeAlias = onp.CanArray2D[_RealT] | _spbase[_RealT, tuple[int, int]]

###

DTYPE: Final[type[np.float64]] = ...
ITYPE: Final[type[np.int32]] = ...

#
@overload
def csgraph_from_masked(graph: onp.MArray2D[npc.integer]) -> csr_array[np.int32, tuple[int, int]]: ...
@overload
def csgraph_from_masked(graph: onp.MArray2D[npc.floating]) -> csr_array[np.float64, tuple[int, int]]: ...
@overload
def csgraph_from_masked(graph: onp.MArray2D[_Real]) -> csr_array[np.float64 | np.int32, tuple[int, int]]: ...

#
@overload
def csgraph_masked_from_dense(
    graph: onp.ToInt2D, null_value: int | None = 0, nan_null: bool = True, infinity_null: bool = True, copy: bool = True
) -> onp.MArray2D[np.int32]: ...
@overload
def csgraph_masked_from_dense(
    graph: onp.ToJustFloat2D, null_value: float | None = 0, nan_null: bool = True, infinity_null: bool = True, copy: bool = True
) -> onp.MArray2D[np.float64]: ...
@overload
def csgraph_masked_from_dense(
    graph: onp.ToFloat2D, null_value: float | None = 0, nan_null: bool = True, infinity_null: bool = True, copy: bool = True
) -> onp.MArray2D[np.float64 | np.int32]: ...

#
@overload
def csgraph_from_dense(
    graph: onp.ToInt2D, null_value: int | None = 0, nan_null: bool = True, infinity_null: bool = True
) -> csr_array[np.int32, tuple[int, int]]: ...
@overload
def csgraph_from_dense(
    graph: onp.ToJustFloat2D, null_value: float | None = 0, nan_null: bool = True, infinity_null: bool = True
) -> csr_array[np.float64, tuple[int, int]]: ...
@overload
def csgraph_from_dense(
    graph: onp.ToFloat2D, null_value: float | None = 0, nan_null: bool = True, infinity_null: bool = True
) -> csr_array[np.float64 | np.int32, tuple[int, int]]: ...

#
@overload
def csgraph_to_dense(csgraph: _SparseGraph[npc.integer], null_value: int | None = 0) -> onp.Array2D[np.int32]: ...
@overload
def csgraph_to_dense(csgraph: _SparseGraph[npc.floating], null_value: float | None = 0) -> onp.Array2D[np.float64]: ...
@overload
def csgraph_to_dense(csgraph: _SparseGraph[_Real], null_value: float | None = 0) -> onp.Array2D[np.float64 | np.int32]: ...

#
@overload
def csgraph_to_masked(csgraph: _SparseGraph[npc.integer]) -> onp.MArray2D[np.int32]: ...
@overload
def csgraph_to_masked(csgraph: _SparseGraph[npc.floating]) -> onp.MArray2D[np.float64]: ...
@overload
def csgraph_to_masked(csgraph: _SparseGraph[_Real]) -> onp.MArray2D[np.float64 | np.int32]: ...

#
@overload
def reconstruct_path(
    csgraph: _Graph[npc.integer], predecessors: onp.ToIntND, directed: bool = True
) -> csr_array[np.int32, tuple[int, int]]: ...
@overload
def reconstruct_path(
    csgraph: _Graph[npc.floating], predecessors: onp.ToFloatND, directed: bool = True
) -> csr_array[np.float64, tuple[int, int]]: ...
@overload
def reconstruct_path(
    csgraph: _ToGraph, predecessors: onp.ToFloatND, directed: bool = True
) -> csr_array[np.float64 | np.int32, tuple[int, int]]: ...

#
@overload
def construct_dist_matrix(
    graph: _Graph[npc.integer], predecessors: onp.ToIntND, directed: bool = True, null_value: float = ...
) -> onp.Array2D[np.int32]: ...
@overload
def construct_dist_matrix(
    graph: _Graph[npc.floating], predecessors: onp.ToFloatND, directed: bool = True, null_value: float = ...
) -> onp.Array2D[np.float64]: ...
@overload
def construct_dist_matrix(
    graph: _ToGraph, predecessors: onp.ToFloatND, directed: bool = True, null_value: float = ...
) -> onp.Array2D[np.float64 | np.int32]: ...
