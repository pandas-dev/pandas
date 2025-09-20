from collections.abc import Callable, Iterable, Mapping
from typing import Concatenate, Final, Literal, Protocol, TypeAlias, type_check_only

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.optimize import Bounds, OptimizeResult as _OptimizeResult
from scipy.optimize._differentiable_functions import _Workers
from scipy.sparse import csr_array
from scipy.sparse._base import _spbase
from scipy.sparse.linalg import LinearOperator

_Ignored: TypeAlias = object

_Float1D: TypeAlias = onp.Array1D[np.float64]
_Float1ND: TypeAlias = onp.Array[onp.AtLeast1D, np.float64]

_LeastSquaresMethod: TypeAlias = Literal["trf", "dogbox", "lm"]

_JacMethod: TypeAlias = Literal["2-point", "3-point", "cs"]
_ToJac2D: TypeAlias = onp.ToFloat2D | _spbase
_JacFunction: TypeAlias = Callable[Concatenate[_Float1D, ...], _ToJac2D | LinearOperator]
_ToJac: TypeAlias = _JacFunction | _JacMethod

_XScaleMethod: TypeAlias = Literal["jac"]
_XScale: TypeAlias = onp.ToFloat | onp.ToFloatND | _XScaleMethod

_LossMethod: TypeAlias = Literal["linear", "soft_l1", "huber", "cauchy", "arctan"]
_Loss: TypeAlias = _UserLossFunction | _LossMethod

_ResidFunction: TypeAlias = Callable[Concatenate[_Float1D, ...], onp.ToFloat1D]

_ResultStatus: TypeAlias = Literal[-2, -1, 0, 1, 2, 3, 4]

@type_check_only
class _UserLossFunction(Protocol):
    def __call__(self, x: _Float1D, /, *, cost_only: bool | None = None) -> onp.ToFloat1D: ...

@type_check_only
class _ImplementedLossFunction(Protocol):
    def __call__(self, /, z: onp.Array1D[np.float64], rho: onp.Array2D[np.float64], cost_only: bool) -> None: ...

@type_check_only
class _BaseOptimizeResult(_OptimizeResult):
    x: _Float1D
    cost: float
    fun: _Float1D
    jac: onp.Array2D[np.float64] | csr_array[np.float64] | LinearOperator[np.float64]
    grad: _Float1D
    optimality: float
    active_mask: onp.Array1D[np.int_]
    nfev: int
    njev: int | None
    status: _ResultStatus

@type_check_only
class _Callback(Protocol):
    def __call__(self, /, *, intermediate_result: _OptimizeResult) -> _Ignored: ...

_ToCallback: TypeAlias = _Callback | Callable[[_Float1D], _Ignored]

###
# undocumented internal machinery

TERMINATION_MESSAGES: Final[dict[_ResultStatus, str]] = ...
FROM_MINPACK_TO_COMMON: Final[dict[Literal[0, 1, 2, 3, 4, 5], Literal[-1, 2, 3, 4, 1, 0]]] = ...

def call_minpack(
    fun: _ResidFunction,
    x0: onp.ToFloat1D,
    jac: _ToJac | None,
    ftol: onp.ToFloat,
    xtol: onp.ToFloat,
    gtol: onp.ToFloat,
    max_nfev: onp.ToInt | None,
    x_scale: onp.ToFloat | _Float1ND,
    jac_method: _ToJac | None = None,
) -> _BaseOptimizeResult: ...
def prepare_bounds(bounds: Iterable[onp.ToFloat | onp.ToFloat1D], n: op.CanIndex) -> tuple[_Float1D, _Float1D]: ...
def check_tolerance(
    ftol: onp.ToFloat | None, xtol: onp.ToFloat | None, gtol: onp.ToFloat | None, method: _LeastSquaresMethod
) -> tuple[onp.ToFloat, onp.ToFloat, onp.ToFloat]: ...
def check_x_scale(x_scale: _XScale, x0: onp.ArrayND[npc.floating], method: _LeastSquaresMethod) -> _Float1ND: ...
def check_jac_sparsity(jac_sparsity: _ToJac2D | None, m: onp.ToInt, n: onp.ToInt) -> _Float1D: ...

#
def huber(z: onp.Array1D[np.float64], rho: onp.Array2D[np.float64], cost_only: bool) -> None: ...
def soft_l1(z: onp.Array1D[np.float64], rho: onp.Array2D[np.float64], cost_only: bool) -> None: ...
def cauchy(z: onp.Array1D[np.float64], rho: onp.Array2D[np.float64], cost_only: bool) -> None: ...
def arctan(z: onp.Array1D[np.float64], rho: onp.Array2D[np.float64], cost_only: bool) -> None: ...

IMPLEMENTED_LOSSES: Final[dict[Literal["linear", "huber", "soft_l1", "cauchy", "arctan"], _ImplementedLossFunction]] = ...

def construct_loss_function(m: op.CanIndex, loss: _Loss, f_scale: onp.ToFloat) -> _UserLossFunction: ...

###
# public API

class OptimizeResult(_BaseOptimizeResult):
    message: str
    success: bool

def least_squares(
    fun: _ResidFunction,
    x0: onp.ToFloat | onp.ToFloat1D,
    jac: _ToJac = "2-point",
    bounds: tuple[onp.ToFloat | onp.ToFloat1D, onp.ToFloat | onp.ToFloat1D] | Bounds = ...,
    method: _LeastSquaresMethod = "trf",
    ftol: onp.ToFloat | None = 1e-8,
    xtol: onp.ToFloat | None = 1e-8,
    gtol: onp.ToFloat | None = 1e-8,
    x_scale: _XScale | None = None,
    loss: _Loss = "linear",
    f_scale: onp.ToFloat = 1.0,
    diff_step: onp.ToFloat1D | None = None,
    tr_solver: Literal["exact", "lsmr"] | None = None,
    tr_options: Mapping[str, object] | None = None,
    jac_sparsity: _ToJac2D | None = None,
    max_nfev: onp.ToInt | None = None,
    verbose: Literal[0, 1, 2] = 0,
    args: Iterable[object] = (),
    kwargs: Mapping[str, object] | None = None,
    callback: _ToCallback | None = None,
    workers: _Workers | None = None,
) -> OptimizeResult: ...
