from collections.abc import Callable
from typing import TypeAlias, TypeVar, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse._base import _spbase
from scipy.sparse.linalg import LinearOperator

__all__ = ["lobpcg"]

_Float: TypeAlias = np.float32 | np.float64
_Complex: TypeAlias = np.complex64 | np.complex128
_FloatT = TypeVar("_FloatT", bound=_Float)

_ToRealMatrix: TypeAlias = (
    onp.ToFloat2D
    | LinearOperator[npc.integer | npc.floating]
    | _spbase
    | Callable[[onp.Array2D[_FloatT]], onp.ArrayND[_Float | _Complex]]
)
_ToComplexMatrix: TypeAlias = (
    onp.ToComplex2D
    | LinearOperator
    | _spbase
    | Callable[[onp.Array2D[_FloatT]], onp.ArrayND[_Float | _Complex]]
)  # fmt: skip

###

@overload  # retLambdaHistory: falsy = ..., retResidualNormsHistory: falsy = ...
def lobpcg(
    A: _ToComplexMatrix[_FloatT],
    X: onp.ArrayND[_FloatT],  # 2d
    B: _ToRealMatrix[_FloatT] | None = None,
    M: _ToRealMatrix[_FloatT] | None = None,
    Y: onp.ArrayND[_FloatT] | None = None,  # 2d
    tol: float | None = None,
    maxiter: int | None = None,
    largest: bool = True,
    verbosityLevel: int = 0,
    retLambdaHistory: onp.ToFalse = False,
    retResidualNormsHistory: onp.ToFalse = False,
    restartControl: int = 20,
) -> tuple[onp.Array1D[_FloatT], onp.Array2D[_FloatT | _Complex]]: ...
@overload  # retLambdaHistory: falsy = ..., retResidualNormsHistory: truthy  (positional)
def lobpcg(
    A: _ToComplexMatrix[_FloatT],
    X: onp.ArrayND[_FloatT],  # 2d
    B: _ToRealMatrix[_FloatT] | None,
    M: _ToRealMatrix[_FloatT] | None,
    Y: onp.ArrayND[_FloatT] | None,  # 2d
    tol: float | None,
    maxiter: int | None,
    largest: bool,
    verbosityLevel: int,
    retLambdaHistory: onp.ToFalse,
    retResidualNormsHistory: onp.ToTrue,
    restartControl: int = 20,
) -> tuple[onp.Array1D[_FloatT], onp.Array2D[_FloatT | _Complex], list[onp.Array0D[_FloatT]]]: ...
@overload  # retLambdaHistory: falsy = ..., retResidualNormsHistory: truthy  (keyword)
def lobpcg(
    A: _ToComplexMatrix[_FloatT],
    X: onp.ArrayND[_FloatT],  # 2d
    B: _ToRealMatrix[_FloatT] | None = None,
    M: _ToRealMatrix[_FloatT] | None = None,
    Y: onp.ArrayND[_FloatT] | None = None,  # 2d
    tol: float | None = None,
    maxiter: int | None = None,
    largest: bool = True,
    verbosityLevel: int = 0,
    retLambdaHistory: onp.ToFalse = False,
    *,
    retResidualNormsHistory: onp.ToTrue,
    restartControl: int = 20,
) -> tuple[onp.Array1D[_FloatT], onp.Array2D[_FloatT | _Complex], list[onp.Array0D[_FloatT]]]: ...
@overload  # retLambdaHistory: truthy  (positional), retResidualNormsHistory: falsy = ...
def lobpcg(
    A: _ToComplexMatrix[_FloatT],
    X: onp.ArrayND[_FloatT],  # 2d
    B: _ToRealMatrix[_FloatT] | None,
    M: _ToRealMatrix[_FloatT] | None,
    Y: onp.ArrayND[_FloatT] | None,  # 2d
    tol: float | None,
    maxiter: int | None,
    largest: bool,
    verbosityLevel: int,
    retLambdaHistory: onp.ToTrue,
    retResidualNormsHistory: onp.ToFalse = False,
    restartControl: int = 20,
) -> tuple[onp.Array1D[_FloatT], onp.Array2D[_FloatT | _Complex], list[onp.Array0D[_FloatT]]]: ...
@overload  # retLambdaHistory: truthy  (keyword), retResidualNormsHistory: falsy = ...
def lobpcg(
    A: _ToComplexMatrix[_FloatT],
    X: onp.ArrayND[_FloatT],  # 2d
    B: _ToRealMatrix[_FloatT] | None = None,
    M: _ToRealMatrix[_FloatT] | None = None,
    Y: onp.ArrayND[_FloatT] | None = None,  # 2d
    tol: float | None = None,
    maxiter: int | None = None,
    largest: bool = True,
    verbosityLevel: int = 0,
    *,
    retLambdaHistory: onp.ToTrue,
    retResidualNormsHistory: onp.ToFalse = False,
    restartControl: int = 20,
) -> tuple[onp.Array1D[_FloatT], onp.Array2D[_FloatT | _Complex], list[onp.Array0D[_FloatT]]]: ...
@overload  # retLambdaHistory: truthy  (positional), retResidualNormsHistory: truthy
def lobpcg(
    A: _ToComplexMatrix[_FloatT],
    X: onp.ArrayND[_FloatT],  # 2d
    B: _ToRealMatrix[_FloatT] | None,
    M: _ToRealMatrix[_FloatT] | None,
    Y: onp.ArrayND[_FloatT] | None,  # 2d
    tol: float | None,
    maxiter: int | None,
    largest: bool,
    verbosityLevel: int,
    retLambdaHistory: onp.ToTrue,
    retResidualNormsHistory: onp.ToTrue,
    restartControl: int = 20,
) -> tuple[onp.Array1D[_FloatT], onp.Array2D[_FloatT | _Complex], list[onp.Array0D[_FloatT]], list[onp.Array0D[_FloatT]]]: ...
@overload  # retLambdaHistory: truthy  (keyword), retResidualNormsHistory: truthy
def lobpcg(
    A: _ToComplexMatrix[_FloatT],
    X: onp.ArrayND[_FloatT],  # 2d
    B: _ToRealMatrix[_FloatT] | None = None,
    M: _ToRealMatrix[_FloatT] | None = None,
    Y: onp.ArrayND[_FloatT] | None = None,  # 2d
    tol: float | None = None,
    maxiter: int | None = None,
    largest: bool = True,
    verbosityLevel: int = 0,
    *,
    retLambdaHistory: onp.ToTrue,
    retResidualNormsHistory: onp.ToTrue,
    restartControl: int = 20,
) -> tuple[onp.Array1D[_FloatT], onp.Array2D[_FloatT | _Complex], list[onp.Array0D[_FloatT]], list[onp.Array0D[_FloatT]]]: ...
