from _typeshed import Unused
from collections.abc import Callable, Iterable, Mapping
from typing import Any, Concatenate, Final, Generic, Literal, Protocol, SupportsIndex, overload, type_check_only
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.optimize import Bounds, OptimizeResult as _OptimizeResult
from scipy.optimize._differentiable_functions import _Workers
from scipy.sparse import csr_array
from scipy.sparse._base import _spbase
from scipy.sparse.linalg import LinearOperator

###

type _Float1D = onp.Array1D[np.float64]
type _Float1ND = onp.Array[onp.AtLeast1D[Any], np.float64]

type _LeastSquaresMethod = Literal["trf", "dogbox", "lm"]

type _JacMethod = Literal["2-point", "3-point", "cs"]
type _ToJac2D = onp.ToFloat2D | _spbase
type _JacFunction = Callable[Concatenate[_Float1D, ...], _ToJac2D | LinearOperator]
type _ToJac = _JacFunction | _JacMethod

type _ToBounds = tuple[onp.ToFloat | onp.ToFloat1D, onp.ToFloat | onp.ToFloat1D] | Bounds

type _XScaleMethod = Literal["jac"]
type _XScale = onp.ToFloat | onp.ToFloatND | _XScaleMethod

type _LossMethod = Literal["linear", "soft_l1", "huber", "cauchy", "arctan"]
type _Loss = _UserLossFunction | _LossMethod

type _ResidFunction = Callable[Concatenate[_Float1D, ...], onp.ToFloat1D | onp.ToFloat]

type _ResultStatus = Literal[-2, -1, 0, 1, 2, 3, 4]

@type_check_only
class _UserLossFunction(Protocol):
    def __call__(self, x: _Float1D, /, *, cost_only: bool | None = None) -> onp.ToFloat1D: ...

@type_check_only
class _ImplementedLossFunction(Protocol):
    def __call__(self, /, z: onp.Array1D[np.float64], rho: onp.Array2D[np.float64], cost_only: bool) -> None: ...

@type_check_only
class _Callback(Protocol):
    def __call__(self, /, *, intermediate_result: _OptimizeResult) -> Unused: ...

type _ToCallback = _Callback | Callable[[_Float1D], Unused]

_MaskT_co = TypeVar("_MaskT_co", bound=npc.number, default=np.int_ | np.float64, covariant=True)

@type_check_only
class _BaseOptimizeResult(_OptimizeResult, Generic[_MaskT_co]):
    x: _Float1D
    cost: float
    fun: _Float1D
    jac: onp.Array2D[np.float64] | csr_array[np.float64] | LinearOperator[np.float64]
    grad: _Float1D
    optimality: float
    active_mask: onp.Array1D[_MaskT_co]
    nfev: int
    njev: int | None
    status: _ResultStatus

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
    max_nfev: int | None,
    x_scale: onp.ToFloat | _Float1ND,
    jac_method: _ToJac | None = None,
) -> _BaseOptimizeResult[np.int_]: ...
def prepare_bounds(bounds: Iterable[onp.ToFloat | onp.ToFloat1D], n: SupportsIndex) -> tuple[_Float1D, _Float1D]: ...
def check_tolerance(
    ftol: onp.ToFloat | None, xtol: onp.ToFloat | None, gtol: onp.ToFloat | None, method: _LeastSquaresMethod
) -> tuple[onp.ToFloat, onp.ToFloat, onp.ToFloat]: ...
def check_x_scale(x_scale: _XScale, x0: onp.ArrayND[npc.floating], method: _LeastSquaresMethod) -> _Float1ND: ...
def check_jac_sparsity(jac_sparsity: _ToJac2D | None, m: int, n: int) -> _Float1D: ...

#
def huber(z: onp.Array1D[np.float64], rho: onp.Array2D[np.float64], cost_only: bool) -> None: ...
def soft_l1(z: onp.Array1D[np.float64], rho: onp.Array2D[np.float64], cost_only: bool) -> None: ...
def cauchy(z: onp.Array1D[np.float64], rho: onp.Array2D[np.float64], cost_only: bool) -> None: ...
def arctan(z: onp.Array1D[np.float64], rho: onp.Array2D[np.float64], cost_only: bool) -> None: ...

IMPLEMENTED_LOSSES: Final[dict[Literal["linear", "huber", "soft_l1", "cauchy", "arctan"], _ImplementedLossFunction]] = ...

def construct_loss_function(m: SupportsIndex, loss: _Loss, f_scale: onp.ToFloat) -> _UserLossFunction: ...

###
# public API

class OptimizeResult(_BaseOptimizeResult[_MaskT_co], Generic[_MaskT_co]):
    message: str
    success: bool

@overload  # method: "trf"
def least_squares(
    fun: _ResidFunction,
    x0: onp.ToFloat | onp.ToFloat1D,
    jac: _ToJac = "2-point",
    bounds: _ToBounds = ...,
    method: Literal["trf"] = "trf",
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
    max_nfev: int | None = None,
    verbose: Literal[0, 1, 2] = 0,
    args: Iterable[object] = (),
    kwargs: Mapping[str, object] | None = None,
    callback: _ToCallback | None = None,
    workers: _Workers | None = None,
) -> OptimizeResult[np.float64]: ...
@overload  # method: {"dogbox", "lm"}  (keyword)
def least_squares(
    fun: _ResidFunction,
    x0: onp.ToFloat | onp.ToFloat1D,
    jac: _ToJac = "2-point",
    bounds: _ToBounds = ...,
    *,
    method: Literal["dogbox", "lm"],
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
    max_nfev: int | None = None,
    verbose: Literal[0, 1, 2] = 0,
    args: Iterable[object] = (),
    kwargs: Mapping[str, object] | None = None,
    callback: _ToCallback | None = None,
    workers: _Workers | None = None,
) -> OptimizeResult[np.int_]: ...
@overload  # method: {"dogbox", "lm"}  (positional)
def least_squares(
    fun: _ResidFunction,
    x0: onp.ToFloat | onp.ToFloat1D,
    jac: _ToJac,
    bounds: _ToBounds,
    method: Literal["dogbox", "lm"],
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
    max_nfev: int | None = None,
    verbose: Literal[0, 1, 2] = 0,
    args: Iterable[object] = (),
    kwargs: Mapping[str, object] | None = None,
    callback: _ToCallback | None = None,
    workers: _Workers | None = None,
) -> OptimizeResult[np.int_]: ...
