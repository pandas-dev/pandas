from collections.abc import Callable
from typing import overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse._base import _spbase
from scipy.sparse.linalg import LinearOperator

__all__ = ["lobpcg"]

###

type _Float = np.float32 | np.float64
type _Complex = np.complex64 | np.complex128

type _ToRealMatrix[FloatT: _Float] = (
    onp.ToFloat2D
    | LinearOperator[npc.integer | npc.floating]
    | _spbase
    | Callable[[onp.Array2D[FloatT]], onp.ArrayND[_Float | _Complex]]
)
type _ToComplexMatrix[FloatT: _Float] = (
    onp.ToComplex2D
    | LinearOperator
    | _spbase
    | Callable[[onp.Array2D[FloatT]], onp.ArrayND[_Float | _Complex]]
)  # fmt: skip

###

@overload  # retLambdaHistory: falsy = ..., retResidualNormsHistory: falsy = ...
def lobpcg[FloatT: _Float](
    A: _ToComplexMatrix[FloatT],
    X: onp.ArrayND[FloatT],  # 2d
    B: _ToRealMatrix[FloatT] | None = None,
    M: _ToRealMatrix[FloatT] | None = None,
    Y: onp.ArrayND[FloatT] | None = None,  # 2d
    tol: float | None = None,
    maxiter: int | None = None,
    largest: bool = True,
    verbosityLevel: int = 0,
    retLambdaHistory: onp.ToFalse = False,
    retResidualNormsHistory: onp.ToFalse = False,
    restartControl: int = 20,
) -> tuple[onp.Array1D[FloatT], onp.Array2D[FloatT | _Complex]]: ...
@overload  # retLambdaHistory: falsy = ..., retResidualNormsHistory: truthy  (positional)
def lobpcg[FloatT: _Float](
    A: _ToComplexMatrix[FloatT],
    X: onp.ArrayND[FloatT],  # 2d
    B: _ToRealMatrix[FloatT] | None,
    M: _ToRealMatrix[FloatT] | None,
    Y: onp.ArrayND[FloatT] | None,  # 2d
    tol: float | None,
    maxiter: int | None,
    largest: bool,
    verbosityLevel: int,
    retLambdaHistory: onp.ToFalse,
    retResidualNormsHistory: onp.ToTrue,
    restartControl: int = 20,
) -> tuple[onp.Array1D[FloatT], onp.Array2D[FloatT | _Complex], list[onp.Array0D[FloatT]]]: ...
@overload  # retLambdaHistory: falsy = ..., retResidualNormsHistory: truthy  (keyword)
def lobpcg[FloatT: _Float](
    A: _ToComplexMatrix[FloatT],
    X: onp.ArrayND[FloatT],  # 2d
    B: _ToRealMatrix[FloatT] | None = None,
    M: _ToRealMatrix[FloatT] | None = None,
    Y: onp.ArrayND[FloatT] | None = None,  # 2d
    tol: float | None = None,
    maxiter: int | None = None,
    largest: bool = True,
    verbosityLevel: int = 0,
    retLambdaHistory: onp.ToFalse = False,
    *,
    retResidualNormsHistory: onp.ToTrue,
    restartControl: int = 20,
) -> tuple[onp.Array1D[FloatT], onp.Array2D[FloatT | _Complex], list[onp.Array0D[FloatT]]]: ...
@overload  # retLambdaHistory: truthy  (positional), retResidualNormsHistory: falsy = ...
def lobpcg[FloatT: _Float](
    A: _ToComplexMatrix[FloatT],
    X: onp.ArrayND[FloatT],  # 2d
    B: _ToRealMatrix[FloatT] | None,
    M: _ToRealMatrix[FloatT] | None,
    Y: onp.ArrayND[FloatT] | None,  # 2d
    tol: float | None,
    maxiter: int | None,
    largest: bool,
    verbosityLevel: int,
    retLambdaHistory: onp.ToTrue,
    retResidualNormsHistory: onp.ToFalse = False,
    restartControl: int = 20,
) -> tuple[onp.Array1D[FloatT], onp.Array2D[FloatT | _Complex], list[onp.Array0D[FloatT]]]: ...
@overload  # retLambdaHistory: truthy  (keyword), retResidualNormsHistory: falsy = ...
def lobpcg[FloatT: _Float](
    A: _ToComplexMatrix[FloatT],
    X: onp.ArrayND[FloatT],  # 2d
    B: _ToRealMatrix[FloatT] | None = None,
    M: _ToRealMatrix[FloatT] | None = None,
    Y: onp.ArrayND[FloatT] | None = None,  # 2d
    tol: float | None = None,
    maxiter: int | None = None,
    largest: bool = True,
    verbosityLevel: int = 0,
    *,
    retLambdaHistory: onp.ToTrue,
    retResidualNormsHistory: onp.ToFalse = False,
    restartControl: int = 20,
) -> tuple[onp.Array1D[FloatT], onp.Array2D[FloatT | _Complex], list[onp.Array0D[FloatT]]]: ...
@overload  # retLambdaHistory: truthy  (positional), retResidualNormsHistory: truthy
def lobpcg[FloatT: _Float](
    A: _ToComplexMatrix[FloatT],
    X: onp.ArrayND[FloatT],  # 2d
    B: _ToRealMatrix[FloatT] | None,
    M: _ToRealMatrix[FloatT] | None,
    Y: onp.ArrayND[FloatT] | None,  # 2d
    tol: float | None,
    maxiter: int | None,
    largest: bool,
    verbosityLevel: int,
    retLambdaHistory: onp.ToTrue,
    retResidualNormsHistory: onp.ToTrue,
    restartControl: int = 20,
) -> tuple[onp.Array1D[FloatT], onp.Array2D[FloatT | _Complex], list[onp.Array0D[FloatT]], list[onp.Array0D[FloatT]]]: ...
@overload  # retLambdaHistory: truthy  (keyword), retResidualNormsHistory: truthy
def lobpcg[FloatT: _Float](
    A: _ToComplexMatrix[FloatT],
    X: onp.ArrayND[FloatT],  # 2d
    B: _ToRealMatrix[FloatT] | None = None,
    M: _ToRealMatrix[FloatT] | None = None,
    Y: onp.ArrayND[FloatT] | None = None,  # 2d
    tol: float | None = None,
    maxiter: int | None = None,
    largest: bool = True,
    verbosityLevel: int = 0,
    *,
    retLambdaHistory: onp.ToTrue,
    retResidualNormsHistory: onp.ToTrue,
    restartControl: int = 20,
) -> tuple[onp.Array1D[FloatT], onp.Array2D[FloatT | _Complex], list[onp.Array0D[FloatT]], list[onp.Array0D[FloatT]]]: ...
