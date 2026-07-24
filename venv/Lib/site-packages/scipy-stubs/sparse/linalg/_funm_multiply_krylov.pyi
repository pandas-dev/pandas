from collections.abc import Callable
from typing import Literal

import optype.numpy as onp
import optype.numpy.compat as npc

from ._interface import LinearOperator
from scipy.sparse._base import _spbase

__all__ = ["funm_multiply_krylov"]

def funm_multiply_krylov[
    MatrixT: onp.Array2D[npc.inexact] | LinearOperator[npc.inexact] | _spbase[npc.inexact, tuple[int, int]],
    ScalarT: npc.inexact,
](
    f: Callable[[MatrixT], onp.ArrayND[ScalarT]],
    A: MatrixT,
    b: onp.Array1D[ScalarT],
    *,
    assume_a: Literal["general", "gen", "hermitian", "her"] = "general",
    t: float = 1.0,
    atol: float = 0.0,
    rtol: float = 1e-6,
    restart_every_m: int | None = None,
    max_restarts: int = 20,
) -> onp.Array1D[ScalarT]: ...
