from typing import Literal, TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse._base import _spbase

_NumberT = TypeVar("_NumberT", bound=npc.number)
_NumberT2 = TypeVar("_NumberT2", bound=npc.number)

###

def _row_count(A: onp.Array2D[npc.number]) -> onp.Array1D[np.intp]: ...
def _get_densest(A: onp.Array2D[_NumberT], eligibleRows: onp.Array1D[np.bool_]) -> np.intp: ...
def _remove_zero_rows(
    A: onp.Array2D[_NumberT], b: onp.Array1D[_NumberT2]
) -> tuple[onp.Array2D[_NumberT], onp.Array1D[_NumberT2], Literal[0, 2], str]: ...

#
def bg_update_dense(
    plu: tuple[onp.Array1D[_NumberT], onp.Array1D[_NumberT2]], perm_r: onp.Array1D[np.intp], v: onp.Array1D[npc.number], j: int
) -> tuple[onp.Array1D[_NumberT], onp.Array1D[_NumberT2]]: ...

#
def _remove_redundancy_pivot_dense(
    A: onp.Array2D[_NumberT], rhs: onp.Array1D[_NumberT2], true_rank: int | None = None
) -> tuple[onp.Array2D[_NumberT], onp.Array1D[_NumberT2], int, str]: ...
def _remove_redundancy_pivot_sparse(
    A: _spbase[_NumberT, tuple[int, int]], rhs: onp.Array1D[_NumberT2]
) -> tuple[onp.Array2D[_NumberT], onp.Array1D[_NumberT2], int, str]: ...
def _remove_redundancy_svd(
    A: onp.Array2D[_NumberT], b: onp.Array1D[_NumberT2]
) -> tuple[onp.Array2D[_NumberT], onp.Array1D[_NumberT2], int, str]: ...
def _remove_redundancy_id(
    A: onp.Array2D[_NumberT], rhs: onp.Array1D[_NumberT2], rank: int | None = None, randomized: bool = True
) -> tuple[onp.Array2D[_NumberT], onp.Array1D[_NumberT2], int, str]: ...
