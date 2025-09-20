from collections.abc import Callable, Iterable
from typing import Self, TypeAlias, TypeVar

import numpy as np
import optype as op
import optype.numpy as onp

from scipy.optimize._constraints import PreparedConstraint
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import LinearOperator

_T = TypeVar("_T")
_Tuple2: TypeAlias = tuple[_T, _T]

_FunConstr: TypeAlias = Callable[[onp.Array1D[np.float64]], _Tuple2[onp.Array1D[np.float64]]]
_FunJac: TypeAlias = Callable[[onp.Array1D[np.float64]], _Tuple2[onp.Array2D[np.float64] | csr_matrix]]
_FunHess: TypeAlias = Callable[
    [onp.Array1D[np.float64], onp.Array1D[np.float64], onp.Array1D[np.float64]],
    _Tuple2[onp.Array2D[np.float64] | csr_matrix | LinearOperator],
]

_PreparedConstraints: TypeAlias = Iterable[CanonicalConstraint]

class CanonicalConstraint:
    n_eq: int
    n_ineq: int
    fun: _FunConstr
    jac: _FunJac
    hess: _FunHess

    keep_feasible: onp.Array1D[np.bool_]

    def __init__(
        self, /, n_eq: int, n_ineq: int, fun: _FunConstr, jac: _FunJac, hess: _FunHess, keep_feasible: onp.Array1D[np.bool_]
    ) -> None: ...
    @classmethod
    def from_PreparedConstraint(cls, constraint: PreparedConstraint) -> Self: ...
    @classmethod
    def empty(cls, n: op.CanIndex) -> Self: ...
    @classmethod
    def concatenate(cls, canonical_constraints: _PreparedConstraints, sparse_jacobian: bool | np.bool_) -> Self: ...

def initial_constraints_as_canonical(
    n: op.CanIndex, prepared_constraints: _PreparedConstraints, sparse_jacobian: bool | np.bool_
) -> tuple[
    onp.Array[onp.AtMost2D, np.float64],
    onp.Array[onp.AtMost2D, np.float64],
    onp.Array2D[np.float64] | csr_matrix,
    onp.Array2D[np.float64] | csr_matrix,
]: ...
