from typing import Final, TypeAlias, overload
from typing_extensions import TypeVar

import numpy as np
import numpy.typing as npt
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse._base import _spbase

_Real: TypeAlias = npc.integer | npc.floating
_ToGraph: TypeAlias = onp.ToFloat2D | _spbase[_Real, tuple[int, int]]

_RealT = TypeVar("_RealT", bound=_Real)
_Graph: TypeAlias = onp.CanArrayND[_RealT] | _spbase[_RealT, tuple[int, int]]

###

DTYPE: Final[type[np.float64]] = ...

@overload  # no dtype
def validate_graph(
    csgraph: _Graph[_RealT],
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
) -> onp.Array2D[_RealT]: ...
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
