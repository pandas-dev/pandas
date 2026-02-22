from collections.abc import Callable, Sequence
from typing import Any, Final, Generic, Literal, TypeAlias, overload, type_check_only
from typing_extensions import TypeAliasType, TypeVar, TypeVarTuple, TypedDict, Unpack

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from .base import DenseOutput, OdeSolver
from .common import OdeSolution
from scipy._lib._util import _RichResult
from scipy.sparse._typing import _Sparse2D

_Ts = TypeVarTuple("_Ts")
_Inexact64T = TypeVar("_Inexact64T", bound=np.float64 | np.complex128)
_Inexact64T_co = TypeVar("_Inexact64T_co", bound=np.float64 | np.complex128, default=np.float64 | np.complex128, covariant=True)

# numpy <2.2 workarounds
_Float: TypeAlias = np.float64 | float
_FloatT = TypeVar("_FloatT", bound=_Float)

_FuncSol: TypeAlias = Callable[[np.float64], onp.ArrayND[_Inexact64T]] | Callable[[float], onp.ArrayND[_Inexact64T]]
_FuncEvent = TypeAliasType(
    "_FuncEvent",
    Callable[[np.float64, onp.ArrayND[_Inexact64T], *_Ts], _Float] | Callable[[float, onp.ArrayND[_Inexact64T], *_Ts], _Float],
    type_params=(_Inexact64T, _Ts),
)
_Events: TypeAlias = Sequence[_FuncEvent[_Inexact64T, *_Ts]] | _FuncEvent[_Inexact64T, *_Ts]

_Int1D: TypeAlias = onp.Array1D[np.int_]
_Float1D: TypeAlias = onp.Array1D[np.float64]
_Float2D: TypeAlias = onp.Array2D[np.float64]
_Complex1D: TypeAlias = onp.Array1D[np.complex128]
_Complex2D: TypeAlias = onp.Array2D[np.complex128]

_ToFloatMax1D: TypeAlias = onp.ToFloat1D | onp.ToFloat
_ToComplexMax1D: TypeAlias = onp.ToComplex1D | onp.ToComplex

_ToJac: TypeAlias = onp.ToArray2D[complex, npc.inexact] | _Sparse2D[npc.inexact]

_IVPMethod: TypeAlias = Literal["RK23", "RK45", "DOP853", "Radau", "BDF", "LSODA"] | type[OdeSolver]

@type_check_only
class _SolverOptions(TypedDict, total=False):
    first_step: float | None
    max_step: float
    rtol: float | onp.ToFloat1D
    atol: float | onp.ToFloat1D
    jac: _ToJac | Callable[[np.float64, onp.Array1D], _ToJac] | None
    jac_sparsity: onp.ToFloat2D | _Sparse2D[npc.floating] | None
    lband: int | None
    uband: int | None
    min_step: float

###

METHODS: Final[dict[str, type]] = ...
MESSAGES: Final[dict[int, str]] = ...

class OdeResult(_RichResult[Any], Generic[_Inexact64T_co]):
    t: _Float1D
    y: onp.Array2D[_Inexact64T_co]
    sol: OdeSolution | None
    t_events: list[_Float1D] | None
    y_events: list[onp.ArrayND[_Inexact64T_co]] | None
    nfev: int
    njev: int
    nlu: int
    status: Literal[-1, 0, 1]
    message: str
    success: bool

def prepare_events(events: _Events[_Inexact64T]) -> tuple[_Events[_Inexact64T], _Float1D, _Float1D]: ...
def solve_event_equation(event: _FuncEvent[_Inexact64T], sol: _FuncSol[_Inexact64T], t_old: float, t: float) -> float: ...
def handle_events(
    sol: DenseOutput[Any],
    events: Sequence[_FuncEvent[_Inexact64T]],
    active_events: onp.ArrayND[np.intp],
    event_count: onp.ArrayND[np.intp | np.float64],
    max_events: onp.ArrayND[np.intp | np.float64],
    t_old: float,
    t: float,
) -> tuple[_Int1D, _Float1D, bool]: ...
def find_active_events(g: onp.ToFloat1D, g_new: onp.ToFloat1D, direction: onp.ArrayND[np.float64]) -> _Int1D: ...

# NOTE: The *free* `_FloatT` type variable works around `float64` not being a subtype of `float` on `numpy <2.2`.
@overload  # float, vectorized=False (default), args=None (default)
def solve_ivp(
    fun: Callable[[_FloatT, _Float1D], _ToFloatMax1D],
    t_span: Sequence[float],
    y0: onp.ToFloat1D,
    method: _IVPMethod = "RK45",
    t_eval: onp.ToFloat1D | None = None,
    dense_output: bool = False,
    events: _Events[np.float64] | None = None,
    vectorized: onp.ToFalse = False,
    args: None = None,
    **options: Unpack[_SolverOptions],
) -> OdeResult[np.float64]: ...
@overload  # float, vectorized=False (default), args=<given>
def solve_ivp(
    fun: Callable[[_FloatT, _Float1D, *_Ts], _ToFloatMax1D],
    t_span: Sequence[float],
    y0: onp.ToFloat1D,
    method: _IVPMethod = "RK45",
    t_eval: onp.ToFloat1D | None = None,
    dense_output: bool = False,
    events: _Events[np.float64] | None = None,
    vectorized: onp.ToFalse = False,
    *,
    args: tuple[*_Ts],
    **options: Unpack[_SolverOptions],
) -> OdeResult[np.float64]: ...
@overload  # float, vectorized=True, args=None (default)
def solve_ivp(
    fun: Callable[[_Float1D, _Float2D], onp.ToFloat2D],
    t_span: Sequence[float],
    y0: onp.ToFloat1D,
    method: _IVPMethod = "RK45",
    t_eval: onp.ToFloat1D | None = None,
    dense_output: bool = False,
    events: _Events[np.float64] | None = None,
    *,
    vectorized: onp.ToTrue,
    args: None = None,
    **options: Unpack[_SolverOptions],
) -> OdeResult[np.float64]: ...
@overload  # float, vectorized=True, args=<given>
def solve_ivp(
    fun: Callable[[_Float1D, _Float2D, *_Ts], onp.ToFloat2D],
    t_span: Sequence[float],
    y0: onp.ToFloat1D,
    method: _IVPMethod = "RK45",
    t_eval: onp.ToFloat1D | None = None,
    dense_output: bool = False,
    events: _Events[np.float64] | None = None,
    *,
    vectorized: onp.ToTrue,
    args: tuple[*_Ts],
    **options: Unpack[_SolverOptions],
) -> OdeResult[np.float64]: ...
@overload  # complex, vectorized=False (default), args=None (default)
def solve_ivp(
    fun: Callable[[_FloatT, _Complex1D], _ToComplexMax1D],
    t_span: Sequence[float],
    y0: onp.ToComplex1D,
    method: _IVPMethod = "RK45",
    t_eval: onp.ToFloat1D | None = None,
    dense_output: bool = False,
    events: _Events[np.complex128] | None = None,
    vectorized: onp.ToFalse = False,
    args: None = None,
    **options: Unpack[_SolverOptions],
) -> OdeResult[np.complex128]: ...
@overload  # complex, vectorized=False (default), args=<given>
def solve_ivp(
    fun: Callable[[_FloatT, _Complex1D, *_Ts], _ToComplexMax1D],
    t_span: Sequence[float],
    y0: onp.ToComplex1D,
    method: _IVPMethod = "RK45",
    t_eval: onp.ToFloat1D | None = None,
    dense_output: bool = False,
    events: _Events[np.complex128] | None = None,
    vectorized: onp.ToFalse = False,
    *,
    args: tuple[*_Ts],
    **options: Unpack[_SolverOptions],
) -> OdeResult[np.complex128]: ...
@overload  # complex, vectorized=True, args=None (default)
def solve_ivp(
    fun: Callable[[_Float1D, _Complex2D], onp.ToComplex2D],
    t_span: Sequence[float],
    y0: onp.ToComplex1D,
    method: _IVPMethod = "RK45",
    t_eval: onp.ToFloat1D | None = None,
    dense_output: bool = False,
    events: _Events[np.complex128] | None = None,
    *,
    vectorized: onp.ToTrue,
    args: None = None,
    **options: Unpack[_SolverOptions],
) -> OdeResult[np.complex128]: ...
@overload  # complex, vectorized=True, args=<given>
def solve_ivp(
    fun: Callable[[_Float1D, _Complex2D, *_Ts], onp.ToComplex2D],
    t_span: Sequence[float],
    y0: onp.ToComplex1D,
    method: _IVPMethod = "RK45",
    t_eval: onp.ToFloat1D | None = None,
    dense_output: bool = False,
    events: _Events[np.complex128] | None = None,
    *,
    vectorized: onp.ToTrue,
    args: tuple[*_Ts],
    **options: Unpack[_SolverOptions],
) -> OdeResult[np.complex128]: ...
