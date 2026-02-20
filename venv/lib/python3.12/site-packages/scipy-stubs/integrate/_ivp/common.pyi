from collections.abc import Callable
from typing import Any, Final, Generic, Literal, TypeAlias, overload
from typing_extensions import TypeVar

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from .base import DenseOutput
from scipy.sparse import csc_matrix

_FloatingT = TypeVar("_FloatingT", bound=npc.floating)
_ToFloatT = TypeVar("_ToFloatT", bound=onp.ToFloat)
_InterpT_co = TypeVar("_InterpT_co", bound=DenseOutput[npc.inexact], default=DenseOutput[Any], covariant=True)

_ToFloat64: TypeAlias = np.float16 | np.float32 | np.float64 | npc.integer | np.bool_

###

EPS: Final[float] = ...
NUM_JAC_DIFF_REJECT: Final[float] = ...
NUM_JAC_DIFF_SMALL: Final[float] = ...
NUM_JAC_DIFF_BIG: Final[float] = ...
NUM_JAC_MIN_FACTOR: Final[float] = ...
NUM_JAC_FACTOR_INCREASE: Final[float] = 10
NUM_JAC_FACTOR_DECREASE: Final[float] = 0.1

class OdeSolution(Generic[_InterpT_co]):
    interpolants: list[_InterpT_co]
    ts: onp.Array1D[np.float64]
    ts_sorted: onp.Array1D[np.float64]
    t_min: np.float64
    t_max: np.float64
    ascending: bool
    side: Literal["left", "right"]
    n_segments: int

    def __init__(self, /, ts: onp.ToFloat1D, interpolants: list[_InterpT_co], alt_segment: op.CanBool = False) -> None: ...

    #
    @overload
    def __call__(self, /, t: float | _ToFloat64) -> onp.Array1D[np.float64]: ...
    @overload
    def __call__(self, /, t: op.JustComplex | np.complex128 | np.complex64) -> onp.Array1D[np.complex128]: ...
    @overload
    def __call__(self, /, t: npc.floating80) -> onp.Array1D[np.longdouble]: ...
    @overload
    def __call__(self, /, t: npc.complexfloating160) -> onp.Array1D[np.clongdouble]: ...
    @overload
    def __call__(self, /, t: onp.ToArray1D[float, _ToFloat64]) -> onp.Array2D[np.float64]: ...
    @overload
    def __call__(self, /, t: onp.ToArray1D[op.JustComplex, np.complex128 | np.complex64]) -> onp.Array2D[np.complex128]: ...
    @overload
    def __call__(self, /, t: onp.CanArrayND[npc.complexfloating160]) -> onp.Array2D[np.clongdouble]: ...

def validate_first_step(first_step: _ToFloatT, t0: onp.ToFloat, t_bound: onp.ToFloat) -> _ToFloatT: ...  # undocumented
def validate_max_step(max_step: _ToFloatT) -> _ToFloatT: ...  # undocumented
def warn_extraneous(extraneous: dict[str, Any]) -> None: ...  # undocumented
def validate_tol(
    rtol: onp.ArrayND[_FloatingT], atol: onp.ArrayND[_FloatingT], n: int
) -> tuple[onp.Array1D[_FloatingT], onp.Array1D[_FloatingT]]: ...  # undocumented
def norm(x: onp.ToFloatND) -> npc.floating: ...  # undocumented
def select_initial_step(
    fun: Callable[[np.float64, onp.Array1D[np.float64]], onp.Array1D[np.float64]],
    t0: float | np.float64,
    y0: onp.ArrayND[np.float64],
    t_bound: float | np.float64,
    max_step: float | np.float64,
    f0: onp.ArrayND[np.float64],
    direction: float | np.float64,
    order: float | np.float64,
    rtol: float | np.float64,
    atol: float | np.float64,
) -> float: ...  # undocumented
def num_jac(
    fun: Callable[[np.float64, onp.Array1D[np.float64]], onp.Array1D[np.float64]],
    t: float | np.float64,
    y: onp.ArrayND[np.float64],
    f: onp.ArrayND[np.float64],
    threshold: float | np.float64,
    factor: onp.ArrayND[np.float64] | None,
    sparsity: tuple[csc_matrix, onp.ArrayND[np.intp]] | None = None,
) -> tuple[onp.Array2D[np.float64] | csc_matrix[np.float64], onp.Array1D[np.float64]]: ...  # undocumented
