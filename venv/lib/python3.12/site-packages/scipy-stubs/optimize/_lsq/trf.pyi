from collections.abc import Callable, Mapping
from typing import Any, Literal, TypeAlias

import numpy as np
import optype.numpy as onp

from scipy.optimize import OptimizeResult as _OptimizeResult
from scipy.optimize._typing import TRSolver
from scipy.sparse import sparray, spmatrix
from scipy.sparse.linalg import LinearOperator

_Ignored: TypeAlias = object
_ValueFloat: TypeAlias = float | np.float64

_ArrayFloat: TypeAlias = onp.ArrayND[np.float64]
_MatrixFloat: TypeAlias = onp.ArrayND[np.float64] | sparray | spmatrix | LinearOperator

_FunObj: TypeAlias = Callable[[onp.Array1D[np.float64], _ValueFloat], _MatrixFloat]
_FunJac: TypeAlias = Callable[[onp.Array1D[np.float64], _ValueFloat], _MatrixFloat]
_FunLoss: TypeAlias = Callable[[_ValueFloat], _ValueFloat]
_FunCallback: TypeAlias = Callable[[OptimizeResult], _Ignored]

class OptimizeResult(_OptimizeResult):
    x: _ArrayFloat
    jac: _MatrixFloat
    cost: _ValueFloat
    f0: _ValueFloat
    J0: _MatrixFloat
    ftol: _ValueFloat
    xtol: _ValueFloat
    gtol: _ValueFloat
    max_nfev: int
    x_scale: _ValueFloat | _ArrayFloat
    loss_function: _ValueFloat | _ArrayFloat
    tr_solver: TRSolver | None
    tr_options: Mapping[str, Any]

# undocumented
def trf(
    fun: _FunObj,
    jac: _FunJac,
    x0: _ArrayFloat,
    f0: _ValueFloat,
    J0: _MatrixFloat,
    lb: _ArrayFloat,
    ub: _ArrayFloat,
    ftol: _ValueFloat,
    xtol: _ValueFloat,
    gtol: _ValueFloat,
    max_nfev: int,
    x_scale: Literal["jac"] | _ValueFloat | _ArrayFloat,
    loss_function: _FunLoss,
    tr_solver: TRSolver | None,
    tr_options: Mapping[str, object],
    verbose: bool,
    callback: _FunCallback | None = None,
) -> OptimizeResult: ...

# undocumented
def trf_bounds(
    fun: _FunObj,
    jac: _FunJac,
    x0: _ArrayFloat,
    f0: _ValueFloat,
    J0: _MatrixFloat,
    lb: _ArrayFloat,
    ub: _ArrayFloat,
    ftol: _ValueFloat,
    xtol: _ValueFloat,
    gtol: _ValueFloat,
    max_nfev: int,
    x_scale: Literal["jac"] | _ValueFloat | _ArrayFloat,
    loss_function: _FunLoss,
    tr_solver: TRSolver | None,
    tr_options: Mapping[str, object],
    verbose: bool,
    callback: _FunCallback | None = None,
) -> OptimizeResult: ...

# undocumented
def trf_no_bounds(
    fun: _FunObj,
    jac: _FunJac,
    x0: _ArrayFloat,
    f0: _ValueFloat,
    J0: _MatrixFloat,
    ftol: _ValueFloat,
    xtol: _ValueFloat,
    gtol: _ValueFloat,
    max_nfev: int,
    x_scale: Literal["jac"] | _ValueFloat | _ArrayFloat,
    loss_function: _FunLoss,
    tr_solver: TRSolver | None,
    tr_options: Mapping[str, object],
    verbose: bool,
    callback: _FunCallback | None = None,
) -> OptimizeResult: ...

# undocumented
def select_step(
    x: _ArrayFloat,
    J_h: _MatrixFloat,
    diag_h: _ArrayFloat,
    g_h: _ArrayFloat,
    p: _ArrayFloat,
    p_h: _ArrayFloat,
    d: _ArrayFloat,
    Delta: _ValueFloat,
    lb: _ArrayFloat,
    ub: _ArrayFloat,
    theta: _ValueFloat,
) -> tuple[_ArrayFloat, _ArrayFloat, _ValueFloat]: ...
