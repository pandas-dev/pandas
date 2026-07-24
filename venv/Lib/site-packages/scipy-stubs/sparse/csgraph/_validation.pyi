from typing import Final, overload

import numpy as np
import numpy.typing as npt
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse._base import _spbase

###

type _Real = npc.integer | npc.floating
type _ToGraph = onp.ToFloat2D | _spbase[_Real, tuple[int, int]]

type _Graph[RealT: _Real] = onp.CanArrayND[RealT] | _spbase[RealT, tuple[int, int]]

###

DTYPE: Final[type[np.float64]] = ...

@overload  # no dtype
def validate_graph[RealT: _Real](
    csgraph: _Graph[RealT],
    directed: bool,
    dtype: None,
    csr_output: bool = True,
    dense_output: bool = True,
    copy_if_dense: bool = False,
    copy_if_sparse: bool = False,
    null_value_in: float = 0,
    null_value_out: float = ...,  # inf
    infinity_null: bool = True,
    nan_null: bool = True,
) -> onp.Array2D[RealT]: ...
@overload  # default dtype  (float64)
def validate_graph(
    csgraph: _ToGraph,
    directed: bool,
    dtype: type[np.float64] | np.dtype[np.float64] = ...,
    csr_output: bool = True,
    dense_output: bool = True,
    copy_if_dense: bool = False,
    copy_if_sparse: bool = False,
    null_value_in: float = 0,
    null_value_out: float = ...,  # inf
    infinity_null: bool = True,
    nan_null: bool = True,
) -> onp.Array2D[np.float64]: ...
@overload  # unknown dtype
def validate_graph(
    csgraph: _ToGraph,
    directed: bool,
    dtype: npt.DTypeLike = ...,  # np.float64
    csr_output: bool = True,
    dense_output: bool = True,
    copy_if_dense: bool = False,
    copy_if_sparse: bool = False,
    null_value_in: float = 0,
    null_value_out: float = ...,  # inf
    infinity_null: bool = True,
    nan_null: bool = True,
) -> onp.Array2D[_Real]: ...
