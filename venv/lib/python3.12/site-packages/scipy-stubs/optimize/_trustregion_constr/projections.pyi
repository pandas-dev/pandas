from collections.abc import Callable
from typing import Final, Literal, TypeAlias

import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse import sparray, spmatrix
from scipy.sparse.linalg import LinearOperator

__all__ = ["orthogonality", "projections"]

_ToMatrix: TypeAlias = onp.ToFloat2D | spmatrix | sparray

_Projections: TypeAlias = tuple[
    Callable[[onp.ToFloat1D], onp.ToFloat1D],
    Callable[[onp.ToFloat1D], onp.ToFloat1D],
    Callable[[onp.ToFloat1D], onp.ToFloat1D],
]  # fmt: skip

###

sksparse_available: Final[bool] = ...

# undocumented
def orthogonality(A: _ToMatrix, g: onp.ToFloat1D) -> npc.floating: ...

# undocumented
def normal_equation_projections(
    A: _ToMatrix, m: int, n: int, orth_tol: float | npc.floating, max_refin: int, tol: float
) -> _Projections: ...

# undocumented
def augmented_system_projections(
    A: _ToMatrix, m: int, n: int, orth_tol: float | npc.floating, max_refin: int, tol: float
) -> _Projections: ...

# undocumented
def qr_factorization_projections(
    A: _ToMatrix, m: int, n: int, orth_tol: float | npc.floating, max_refin: int, tol: float
) -> _Projections: ...

# undocumented
def svd_factorization_projections(
    A: _ToMatrix, m: int, n: int, orth_tol: float | npc.floating, max_refin: int, tol: float
) -> _Projections: ...

#
def projections(
    A: _ToMatrix,
    method: Literal["NormalEquation", "AugmentedSystem", "QRFactorization", "SVDFactorization"] | None = None,
    orth_tol: float | npc.floating = 1e-12,
    max_refin: onp.ToJustInt = 3,
    tol: float | npc.floating = 1e-15,
) -> tuple[LinearOperator, LinearOperator, LinearOperator]: ...
