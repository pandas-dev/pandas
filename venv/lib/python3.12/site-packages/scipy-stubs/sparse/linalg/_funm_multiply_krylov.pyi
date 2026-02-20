from collections.abc import Callable
from typing import Literal, TypeVar

import optype.numpy as onp
import optype.numpy.compat as npc

from ._interface import LinearOperator
from scipy.sparse._base import _spbase

__all__ = ["funm_multiply_krylov"]

_MatrixT = TypeVar(
    "_MatrixT", bound=onp.Array2D[npc.inexact] | LinearOperator[npc.inexact] | _spbase[npc.inexact, tuple[int, int]]
)
_ScalarT = TypeVar("_ScalarT", bound=npc.inexact)

def funm_multiply_krylov(
    f: Callable[[_MatrixT], onp.ArrayND[_ScalarT]],
    A: _MatrixT,
    b: onp.Array1D[_ScalarT],
    *,
    assume_a: Literal["general", "gen", "hermitian", "her"] = "general",
    t: float = 1.0,
    atol: float = 0.0,
    rtol: float = 1e-6,
    restart_every_m: int | None = None,
    max_restarts: int = 20,
) -> onp.Array1D[_ScalarT]: ...
