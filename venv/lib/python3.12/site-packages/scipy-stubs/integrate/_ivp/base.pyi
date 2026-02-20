import types
from collections.abc import Callable
from typing import Any, ClassVar, Final, Generic, Literal, TypeAlias, overload
from typing_extensions import TypeVar

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

_ScalarT = TypeVar("_ScalarT", bound=np.float64 | np.complex128, default=np.float64 | Any)
_ScalarT_co = TypeVar("_ScalarT_co", bound=npc.inexact, default=np.float64 | Any, covariant=True)

_ToFunReal: TypeAlias = Callable[[float, onp.ArrayND[np.float64]], onp.ToFloatND]
_ToFunComplex: TypeAlias = Callable[[float, onp.ArrayND[np.complex128]], onp.ToComplexND]

@overload
def check_arguments(
    fun: _ToFunReal, y0: onp.ToFloatND, support_complex: bool
) -> Callable[[float, onp.ArrayND[np.float64]], onp.ArrayND[np.float64]]: ...
@overload
def check_arguments(
    fun: _ToFunComplex, y0: onp.ToJustComplexND, support_complex: Literal[True]
) -> Callable[[float, onp.ArrayND[np.complex128]], onp.ArrayND[np.complex128]]: ...

class OdeSolver(Generic[_ScalarT]):
    TOO_SMALL_STEP: ClassVar[str] = ...

    t: float
    t_old: float
    t_bound: float
    y: onp.ArrayND[_ScalarT]
    vectorized: bool
    fun: Callable[[float, onp.ArrayND[_ScalarT]], onp.ArrayND[_ScalarT]]
    fun_single: Callable[[float, onp.Array1D[_ScalarT]], onp.Array1D[_ScalarT]]
    fun_vectorized: Callable[[float, onp.Array2D[_ScalarT]], onp.Array2D[_ScalarT]]
    direction: np.float64
    status: Literal["running", "finished", "failed"]
    n: int
    nfev: int
    njev: int
    nlu: int

    @classmethod
    def __class_getitem__(cls, arg: type | object, /) -> types.GenericAlias: ...

    #
    @overload
    def __init__(
        self: OdeSolver[np.float64],
        /,
        fun: _ToFunReal,
        t0: float,
        y0: onp.ToFloatND,
        t_bound: float,
        vectorized: bool,
        support_complex: onp.ToBool = False,
    ) -> None: ...
    @overload
    def __init__(
        self: OdeSolver[np.complex128],
        /,
        fun: _ToFunComplex,
        t0: float,
        y0: onp.ToJustComplexND,
        t_bound: onp.ToFloat,
        vectorized: bool,
        support_complex: onp.ToTrue,
    ) -> None: ...
    @property
    def step_size(self, /) -> np.float64 | None: ...
    def step(self, /) -> str | None: ...
    def dense_output(self, /) -> ConstantDenseOutput[_ScalarT]: ...

class DenseOutput(Generic[_ScalarT_co]):
    t_old: Final[float]
    t: Final[float]
    t_min: Final[float]
    t_max: Final[float]

    @classmethod
    def __class_getitem__(cls, arg: type | object, /) -> types.GenericAlias: ...

    #
    def __init__(self, /, t_old: float, t: float) -> None: ...

    #
    @overload
    def __call__(self, /, t: onp.ToFloat) -> onp.Array1D[_ScalarT_co]: ...
    @overload
    def __call__(self, /, t: onp.ToFloatND) -> onp.ArrayND[_ScalarT_co]: ...

class ConstantDenseOutput(DenseOutput[_ScalarT_co], Generic[_ScalarT_co]):
    value: onp.ArrayND[_ScalarT_co]

    def __init__(self, /, t_old: float, t: float, value: onp.ArrayND[_ScalarT_co]) -> None: ...
