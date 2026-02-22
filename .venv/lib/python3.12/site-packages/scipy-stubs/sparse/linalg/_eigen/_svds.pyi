from collections.abc import Mapping
from typing import Literal, TypeAlias, TypeVar

import numpy as np
import optype.numpy as onp

from scipy.sparse._base import _spbase
from scipy.sparse.linalg import LinearOperator

__all__ = ["svds"]

_SCT = TypeVar("_SCT", bound=np.float32 | np.float64 | np.complex64 | np.complex128)
_ToMatrix: TypeAlias = onp.ArrayND[_SCT] | LinearOperator[_SCT] | _spbase

_Which: TypeAlias = Literal["LM", "SM"]
_ReturnSingularVectors: TypeAlias = Literal["u", "v"] | bool
_Solver: TypeAlias = Literal["arpack", "propack", "lobpcg"]

###

def svds(
    A: _ToMatrix[_SCT],
    k: int = 6,
    ncv: int | None = None,
    tol: float = 0,
    which: _Which = "LM",
    v0: onp.ArrayND[_SCT] | None = None,
    maxiter: int | None = None,
    return_singular_vectors: _ReturnSingularVectors = True,
    solver: _Solver = "arpack",
    rng: onp.random.ToRNG | None = None,
    options: Mapping[str, object] | None = None,
    *,
    random_state: onp.random.ToRNG | None = None,
) -> tuple[onp.Array2D[_SCT], onp.ArrayND[np.float32 | np.float64], onp.ArrayND[_SCT]]: ...
