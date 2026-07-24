from collections.abc import Callable, Iterable
from typing import Self, SupportsIndex

import numpy as np
import optype.numpy as onp

from scipy.optimize._constraints import PreparedConstraint
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import LinearOperator

###

type _Tuple2[T] = tuple[T, T]

type _FunConstr = Callable[[onp.Array1D[np.float64]], _Tuple2[onp.Array1D[np.float64]]]
type _FunJac = Callable[[onp.Array1D[np.float64]], _Tuple2[onp.Array2D[np.float64] | csr_matrix]]
type _FunHess = Callable[
    [onp.Array1D[np.float64], onp.Array1D[np.float64], onp.Array1D[np.float64]],
    _Tuple2[onp.Array2D[np.float64] | csr_matrix | LinearOperator],
]

type _PreparedConstraints = Iterable[CanonicalConstraint]

###

class CanonicalConstraint:
    n_eq: int
    n_ineq: int
    fun: _FunConstr
    jac: _FunJac
    hess: _FunHess

    keep_feasible: onp.Array1D[np.bool]

    def __init__(
        self, /, n_eq: int, n_ineq: int, fun: _FunConstr, jac: _FunJac, hess: _FunHess, keep_feasible: onp.Array1D[np.bool]
    ) -> None: ...
    @classmethod
    def from_PreparedConstraint(cls, constraint: PreparedConstraint) -> Self: ...
    @classmethod
    def empty(cls, n: SupportsIndex) -> Self: ...
    @classmethod
    def concatenate(cls, canonical_constraints: _PreparedConstraints, sparse_jacobian: bool) -> Self: ...

def initial_constraints_as_canonical(
    n: SupportsIndex, prepared_constraints: _PreparedConstraints, sparse_jacobian: bool
) -> tuple[
    onp.Array[onp.AtMost2D, np.float64],
    onp.Array[onp.AtMost2D, np.float64],
    onp.Array2D[np.float64] | csr_matrix[np.float64],
    onp.Array2D[np.float64] | csr_matrix[np.float64],
]: ...
