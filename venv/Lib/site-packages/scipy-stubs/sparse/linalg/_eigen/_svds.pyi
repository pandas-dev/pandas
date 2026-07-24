from collections.abc import Mapping
from typing import Literal, overload

import numpy as np
import optype.numpy as onp

from scipy.sparse._base import _spbase
from scipy.sparse.linalg import LinearOperator

__all__ = ["svds"]

###

type _Inexact = np.float32 | np.float64 | np.complex64 | np.complex128
type _ToMatrix[ScalarT: _Inexact] = onp.ArrayND[ScalarT] | LinearOperator[ScalarT] | _spbase[ScalarT]

type _Which = Literal["LM", "SM"]
type _Solver = Literal["arpack", "propack", "lobpcg"]

###

# mypy reports 4 false positive `overload-overlap` errors on numpy>=2.2
# mypy: disable-error-code=overload-overlap

@overload  # return_singular_vectors=True  (default), f64|c128
def svds[ScalarT: np.float64 | np.complex128](
    A: _ToMatrix[ScalarT],
    k: int = 6,
    ncv: int | None = None,
    tol: float = 0,
    which: _Which = "LM",
    v0: onp.ArrayND[ScalarT] | None = None,
    maxiter: int | None = None,
    return_singular_vectors: Literal[True] = True,
    solver: _Solver = "arpack",
    rng: onp.random.ToRNG | None = None,
    options: Mapping[str, object] | None = None,
    *,
    random_state: onp.random.ToRNG | None = None,
) -> tuple[onp.Array2D[ScalarT], onp.Array1D[np.float64], onp.Array2D[ScalarT]]: ...
@overload  # return_singular_vectors=True  (default), f32|c64
def svds[ScalarT: np.float32 | np.complex64](
    A: _ToMatrix[ScalarT],
    k: int = 6,
    ncv: int | None = None,
    tol: float = 0,
    which: _Which = "LM",
    v0: onp.ArrayND[ScalarT] | None = None,
    maxiter: int | None = None,
    return_singular_vectors: Literal[True] = True,
    solver: _Solver = "arpack",
    rng: onp.random.ToRNG | None = None,
    options: Mapping[str, object] | None = None,
    *,
    random_state: onp.random.ToRNG | None = None,
) -> tuple[onp.Array2D[ScalarT], onp.Array1D[np.float32], onp.Array2D[ScalarT]]: ...
@overload  # return_singular_vectors=False, f64|c128
def svds(
    A: _ToMatrix[np.float64 | np.complex128],
    k: int = 6,
    ncv: int | None = None,
    tol: float = 0,
    which: _Which = "LM",
    v0: onp.ArrayND[np.float64 | np.complex128] | None = None,
    maxiter: int | None = None,
    *,
    return_singular_vectors: Literal[False],
    solver: _Solver = "arpack",
    rng: onp.random.ToRNG | None = None,
    options: Mapping[str, object] | None = None,
    random_state: onp.random.ToRNG | None = None,
) -> onp.Array1D[np.float64]: ...
@overload  # return_singular_vectors=False, f32|c64
def svds(
    A: _ToMatrix[np.float32 | np.complex64],
    k: int = 6,
    ncv: int | None = None,
    tol: float = 0,
    which: _Which = "LM",
    v0: onp.ArrayND[np.float32 | np.complex64] | None = None,
    maxiter: int | None = None,
    *,
    return_singular_vectors: Literal[False],
    solver: _Solver = "arpack",
    rng: onp.random.ToRNG | None = None,
    options: Mapping[str, object] | None = None,
    random_state: onp.random.ToRNG | None = None,
) -> onp.Array1D[np.float32]: ...
@overload  # return_singular_vectors="u", f64|c128
def svds[ScalarT: np.float64 | np.complex128](
    A: _ToMatrix[ScalarT],
    k: int = 6,
    ncv: int | None = None,
    tol: float = 0,
    which: _Which = "LM",
    v0: onp.ArrayND[ScalarT] | None = None,
    maxiter: int | None = None,
    *,
    return_singular_vectors: Literal["u"],
    solver: _Solver = "arpack",
    rng: onp.random.ToRNG | None = None,
    options: Mapping[str, object] | None = None,
    random_state: onp.random.ToRNG | None = None,
) -> tuple[onp.Array2D[ScalarT], onp.Array1D[np.float64], None]: ...
@overload  # return_singular_vectors="u", f32|c64
def svds[ScalarT: np.float32 | np.complex64](
    A: _ToMatrix[ScalarT],
    k: int = 6,
    ncv: int | None = None,
    tol: float = 0,
    which: _Which = "LM",
    v0: onp.ArrayND[ScalarT] | None = None,
    maxiter: int | None = None,
    *,
    return_singular_vectors: Literal["u"],
    solver: _Solver = "arpack",
    rng: onp.random.ToRNG | None = None,
    options: Mapping[str, object] | None = None,
    random_state: onp.random.ToRNG | None = None,
) -> tuple[onp.Array2D[ScalarT], onp.Array1D[np.float32], None]: ...
@overload  # return_singular_vectors="vh", f64|c128
def svds[ScalarT: np.float64 | np.complex128](
    A: _ToMatrix[ScalarT],
    k: int = 6,
    ncv: int | None = None,
    tol: float = 0,
    which: _Which = "LM",
    v0: onp.ArrayND[ScalarT] | None = None,
    maxiter: int | None = None,
    *,
    return_singular_vectors: Literal["vh"],
    solver: _Solver = "arpack",
    rng: onp.random.ToRNG | None = None,
    options: Mapping[str, object] | None = None,
    random_state: onp.random.ToRNG | None = None,
) -> tuple[None, onp.Array1D[np.float64], onp.Array2D[ScalarT]]: ...
@overload  # return_singular_vectors="vh", f32|c64
def svds[ScalarT: np.float32 | np.complex64](
    A: _ToMatrix[ScalarT],
    k: int = 6,
    ncv: int | None = None,
    tol: float = 0,
    which: _Which = "LM",
    v0: onp.ArrayND[ScalarT] | None = None,
    maxiter: int | None = None,
    *,
    return_singular_vectors: Literal["vh"],
    solver: _Solver = "arpack",
    rng: onp.random.ToRNG | None = None,
    options: Mapping[str, object] | None = None,
    random_state: onp.random.ToRNG | None = None,
) -> tuple[None, onp.Array1D[np.float32], onp.Array2D[ScalarT]]: ...
