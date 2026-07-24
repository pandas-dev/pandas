from typing import Literal

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse._base import _spbase

def _row_count(A: onp.Array2D[npc.number]) -> onp.Array1D[np.intp]: ...
def _get_densest[NumberT: npc.number](A: onp.Array2D[NumberT], eligibleRows: onp.Array1D[np.bool]) -> np.intp: ...
def _remove_zero_rows[NumberT1: npc.number, NumberT2: npc.number](
    A: onp.Array2D[NumberT1], b: onp.Array1D[NumberT2]
) -> tuple[onp.Array2D[NumberT1], onp.Array1D[NumberT2], Literal[0, 2], str]: ...

#
def bg_update_dense[NumberT1: npc.number, NumberT2: npc.number](
    plu: tuple[onp.Array1D[NumberT1], onp.Array1D[NumberT2]], perm_r: onp.Array1D[np.intp], v: onp.Array1D[npc.number], j: int
) -> tuple[onp.Array1D[NumberT1], onp.Array1D[NumberT2]]: ...

#
def _remove_redundancy_pivot_dense[NumberT1: npc.number, NumberT2: npc.number](
    A: onp.Array2D[NumberT1], rhs: onp.Array1D[NumberT2], true_rank: int | None = None
) -> tuple[onp.Array2D[NumberT1], onp.Array1D[NumberT2], int, str]: ...
def _remove_redundancy_pivot_sparse[NumberT1: npc.number, NumberT2: npc.number](
    A: _spbase[NumberT1, tuple[int, int]], rhs: onp.Array1D[NumberT2]
) -> tuple[onp.Array2D[NumberT1], onp.Array1D[NumberT2], int, str]: ...
def _remove_redundancy_svd[NumberT1: npc.number, NumberT2: npc.number](
    A: onp.Array2D[NumberT1], b: onp.Array1D[NumberT2]
) -> tuple[onp.Array2D[NumberT1], onp.Array1D[NumberT2], int, str]: ...
def _remove_redundancy_id[NumberT: npc.number, NumberT2: npc.number](
    A: onp.Array2D[NumberT], rhs: onp.Array1D[NumberT2], rank: int | None = None, randomized: bool = True
) -> tuple[onp.Array2D[NumberT], onp.Array1D[NumberT2], int, str]: ...
