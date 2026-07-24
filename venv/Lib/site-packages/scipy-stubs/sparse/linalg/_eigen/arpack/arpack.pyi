from collections.abc import Mapping
from typing import Final, Literal, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse import sparray, spmatrix
from scipy.sparse.linalg import LinearOperator

__all__ = ["ArpackError", "ArpackNoConvergence", "eigs", "eigsh"]

###

type _Numeric = npc.number | np.bool
type _ToFloat = npc.floating | npc.integer | np.bool

type _MatrixOperator[_SCT: _Numeric] = spmatrix[_SCT] | sparray[_SCT, tuple[int, int]] | LinearOperator[_SCT]

type _ToRealMatrix = onp.ToFloat2D | _MatrixOperator[_ToFloat]
type _ToJustComplexMatrix = onp.ToJustComplex2D | _MatrixOperator[npc.complexfloating]
type _ToComplexMatrix = onp.ToComplex2D | _MatrixOperator[_Numeric]

type _Which_eigs = Literal["LM", "SM", "LR", "SR", "LI", "SI"]
type _Which_eigsh = Literal["LM", "SM", "LA", "SA", "BE"]
type _OPpart = Literal["r", "i"]
type _Mode = Literal["normal", "buckling", "cayley"]

###

class ArpackError(RuntimeError):
    def __init__[T](self, /, info: T, infodict: Mapping[T, str] | None = None) -> None: ...

class ArpackNoConvergence(ArpackError):
    eigenvalues: Final[onp.Array1D[np.float64 | np.complex128]]
    eigenvectors: Final[onp.Array2D[np.float64 | np.complex128]]
    def __init__(
        self,
        /,
        msg: str,
        eigenvalues: onp.Array1D[np.float64 | np.complex128],
        eigenvectors: onp.Array2D[np.float64 | np.complex128],
    ) -> None: ...

#
@overload  # returns_eigenvectors: truthy (default)
def eigs(
    A: _ToComplexMatrix,
    k: int = 6,
    M: _ToComplexMatrix | None = None,
    sigma: onp.ToComplex | None = None,
    which: _Which_eigs = "LM",
    v0: onp.ToComplex1D | None = None,
    ncv: int | None = None,
    maxiter: int | None = None,
    tol: float = 0,
    return_eigenvectors: onp.ToTrue = True,
    Minv: _ToComplexMatrix | None = None,
    OPinv: _ToComplexMatrix | None = None,
    OPpart: _OPpart | None = None,
    rng: onp.random.ToRNG | None = None,
) -> tuple[onp.Array1D[np.complex128], onp.Array2D[np.complex128]]: ...
@overload  # returns_eigenvectors: falsy (positional)
def eigs(
    A: _ToComplexMatrix,
    k: int,
    M: _ToComplexMatrix | None,
    sigma: onp.ToComplex | None,
    which: _Which_eigs,
    v0: onp.ToComplex1D | None,
    ncv: int | None,
    maxiter: int | None,
    tol: float,
    return_eigenvectors: onp.ToFalse,
    Minv: _ToComplexMatrix | None = None,
    OPinv: _ToComplexMatrix | None = None,
    OPpart: _OPpart | None = None,
    rng: onp.random.ToRNG | None = None,
) -> onp.Array1D[np.complex128]: ...
@overload  # returns_eigenvectors: falsy (keyword)
def eigs(
    A: _ToComplexMatrix,
    k: int = 6,
    M: _ToComplexMatrix | None = None,
    sigma: onp.ToComplex | None = None,
    which: _Which_eigs = "LM",
    v0: onp.ToComplex1D | None = None,
    ncv: int | None = None,
    maxiter: int | None = None,
    tol: float = 0,
    *,
    return_eigenvectors: onp.ToFalse,
    Minv: _ToComplexMatrix | None = None,
    OPinv: _ToComplexMatrix | None = None,
    OPpart: _OPpart | None = None,
    rng: onp.random.ToRNG | None = None,
) -> onp.Array1D[np.complex128]: ...

#
@overload  # real, returns_eigenvectors: truthy (default)
def eigsh(
    A: _ToRealMatrix,
    k: int = 6,
    M: _ToRealMatrix | None = None,
    sigma: onp.ToFloat | None = None,
    which: _Which_eigsh = "LM",
    v0: onp.ToFloat1D | None = None,
    ncv: int | None = None,
    maxiter: int | None = None,
    tol: float = 0,
    return_eigenvectors: onp.ToTrue = True,
    Minv: _ToRealMatrix | None = None,
    OPinv: _ToRealMatrix | None = None,
    mode: _Mode = "normal",
    rng: onp.random.ToRNG | None = None,
) -> tuple[onp.Array1D[np.float64], onp.Array2D[np.float64]]: ...
@overload  # complex, returns_eigenvectors: truthy (default)
def eigsh(
    A: _ToJustComplexMatrix,
    k: int = 6,
    M: _ToComplexMatrix | None = None,
    sigma: onp.ToFloat | None = None,
    which: _Which_eigsh = "LM",
    v0: onp.ToComplex1D | None = None,
    ncv: int | None = None,
    maxiter: int | None = None,
    tol: float = 0,
    return_eigenvectors: onp.ToTrue = True,
    Minv: _ToComplexMatrix | None = None,
    OPinv: _ToComplexMatrix | None = None,
    mode: _Mode = "normal",
    rng: onp.random.ToRNG | None = None,
) -> tuple[onp.Array1D[np.float64], onp.Array2D[np.complex128]]: ...
@overload  # real or complex (catch-all), returns_eigenvectors: truthy (default)
def eigsh(
    A: _ToComplexMatrix,
    k: int = 6,
    M: _ToComplexMatrix | None = None,
    sigma: onp.ToFloat | None = None,
    which: _Which_eigsh = "LM",
    v0: onp.ToComplex1D | None = None,
    ncv: int | None = None,
    maxiter: int | None = None,
    tol: float = 0,
    return_eigenvectors: onp.ToTrue = True,
    Minv: _ToComplexMatrix | None = None,
    OPinv: _ToComplexMatrix | None = None,
    mode: _Mode = "normal",
    rng: onp.random.ToRNG | None = None,
) -> tuple[onp.Array1D[np.float64], onp.Array2D[np.float64 | np.complex128]]: ...
@overload  # real or complex, returns_eigenvectors: falsy (positional)
def eigsh(
    A: _ToComplexMatrix,
    k: int,
    M: _ToComplexMatrix | None,
    sigma: onp.ToFloat | None,
    which: _Which_eigsh,
    v0: onp.ToComplex1D | None,
    ncv: int | None,
    maxiter: int | None,
    tol: float,
    return_eigenvectors: onp.ToFalse,
    Minv: _ToComplexMatrix | None = None,
    OPinv: _ToComplexMatrix | None = None,
    mode: _Mode = "normal",
    rng: onp.random.ToRNG | None = None,
) -> onp.Array1D[np.float64]: ...
@overload  # real or complex, returns_eigenvectors: falsy (keyword)
def eigsh(
    A: _ToComplexMatrix,
    k: int = 6,
    M: _ToComplexMatrix | None = None,
    sigma: onp.ToFloat | None = None,
    which: _Which_eigsh = "LM",
    v0: onp.ToComplex1D | None = None,
    ncv: int | None = None,
    maxiter: int | None = None,
    tol: float = 0,
    *,
    return_eigenvectors: onp.ToFalse,
    Minv: _ToComplexMatrix | None = None,
    OPinv: _ToComplexMatrix | None = None,
    mode: _Mode = "normal",
    rng: onp.random.ToRNG | None = None,
) -> onp.Array1D[np.float64]: ...
