from collections.abc import Callable
from typing import Any, ClassVar, Final, Generic, Literal, TypeVar, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

_VT = TypeVar("_VT", bound=onp.ArrayND[npc.inexact], default=onp.ArrayND[Any])

class OdeSolver:
    TOO_SMALL_STEP: ClassVar[str] = ...

    t: float
    t_old: float
    t_bound: float
    vectorized: bool
    fun: Callable[[float, onp.ArrayND[np.float64]], onp.ArrayND[np.float64]]
    fun_single: Callable[[float, onp.ArrayND[np.float64]], onp.ArrayND[np.float64]]
    fun_vectorized: Callable[[float, onp.ArrayND[np.float64]], onp.ArrayND[np.float64]]
    direction: float
    n: int
    status: Literal["running", "finished", "failed"]
    nfev: int
    njev: int
    nlu: int

    @overload
    def __init__(
        self,
        /,
        fun: Callable[[float, onp.ArrayND[np.float64]], onp.ToFloatND],
        t0: onp.ToFloatND,
        y0: onp.ToFloatND,
        t_bound: onp.ToFloat,
        vectorized: bool,
        support_complex: onp.ToBool = False,
    ) -> None: ...
    @overload
    def __init__(
        self,
        /,
        fun: Callable[[float, onp.ArrayND[np.float64 | np.complex128]], onp.ToComplexND],
        t0: onp.ToFloat,
        y0: onp.ToComplexND,
        t_bound: onp.ToFloat,
        vectorized: bool,
        support_complex: onp.ToTrue,
    ) -> None: ...
    @property
    def step_size(self, /) -> float | None: ...
    def step(self, /) -> str | None: ...
    def dense_output(self, /) -> ConstantDenseOutput: ...

class DenseOutput:
    t_old: Final[float]
    t: Final[float]
    t_min: Final[float]
    t_max: Final[float]

    def __init__(self, /, t_old: onp.ToFloat, t: onp.ToFloat) -> None: ...
    @overload
    def __call__(self, /, t: onp.ToFloat) -> onp.Array1D[npc.inexact]: ...
    @overload
    def __call__(self, /, t: onp.ToFloatND) -> onp.ArrayND[npc.inexact]: ...

class ConstantDenseOutput(DenseOutput, Generic[_VT]):
    value: _VT
    def __init__(self, /, t_old: onp.ToFloat, t: onp.ToFloat, value: _VT) -> None: ...

def check_arguments(
    fun: Callable[[float, onp.ArrayND[np.float64]], onp.ToComplexND], y0: onp.ToComplexND, support_complex: bool
) -> (
    Callable[[float, onp.ArrayND[np.float64]], onp.ArrayND[np.float64]]
    | Callable[[float, onp.ArrayND[np.float64]], onp.ArrayND[np.complex128]]
): ...
