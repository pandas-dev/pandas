from _typeshed import Incomplete
from typing import Any, Generic, Literal
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._interface import LinearOperator
from scipy.sparse._base import _spbase

__all__ = ["_svdp"]

_NumberT = TypeVar("_NumberT", bound=npc.number)
_NumberT_co = TypeVar("_NumberT_co", bound=npc.number, default=np.float64 | Any, covariant=True)

class _AProd(Generic[_NumberT_co]):
    A: LinearOperator[_NumberT_co]

    def __init__(self, /, A: onp.ToArray2D[_NumberT_co] | _spbase[_NumberT_co, tuple[int, int]]) -> None: ...
    def __call__(
        self,
        /,
        transa: str,
        m: int,
        n: int,
        x: onp.Array1D[npc.number],
        y: onp.Array1D[npc.number],
        sparm: object,  # unused
        iparm: object,  # unused
    ) -> None: ...
    @property
    def shape(self) -> tuple[int, int]: ...
    @property
    def dtype(self) -> np.dtype[_NumberT_co]: ...

def _svdp(
    A: onp.ToArray2D[_NumberT] | _spbase[_NumberT, tuple[int, int]] | LinearOperator[_NumberT],
    k: int,
    which: Literal["LM", "SM"] = "LM",
    irl_mode: bool = True,
    kmax: int | None = None,
    compute_u: bool = True,
    compute_v: bool = True,
    v0: onp.ToFloatND | None = None,
    full_output: bool = False,
    tol: int = 0,
    delta: float | None = None,
    eta: float | None = None,
    anorm: int = 0,
    cgs: bool = False,
    elr: bool = True,
    min_relgap: float = 0.002,
    shifts: int | None = None,
    maxiter: int | None = None,
    rng: onp.random.ToRNG | None = None,
) -> Incomplete: ...  # complicated -- 8 different return types depending on parameters
