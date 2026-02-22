import collections
from collections.abc import Callable
from typing import Any, Concatenate, Final, Generic, Literal, Never, NoReturn, Protocol, TypeAlias, overload, type_check_only
from typing_extensions import TypeVar, override

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

_S = TypeVar("_S")
_T = TypeVar("_T")
_VT = TypeVar("_VT", default=Any)
_NDT_co = TypeVar("_NDT_co", bound=_FloatingND, default=_FloatingND, covariant=True)
_InexactT = TypeVar("_InexactT", bound=npc.inexact)
_InexactT_co = TypeVar("_InexactT_co", bound=npc.inexact, default=Any, covariant=True)
_ShapeT_co = TypeVar("_ShapeT_co", bound=tuple[int, ...], default=tuple[Any, ...], covariant=True)

_Floating: TypeAlias = float | npc.floating
_FloatingND: TypeAlias = onp.ArrayND[npc.floating] | _Floating

_Fun: TypeAlias = Callable[Concatenate[float, ...], _T]

_Norm: TypeAlias = Literal["max", "2"]
_Quadrature: TypeAlias = Literal["gk21", "gk15", "trapezoid"]

@type_check_only
class _DoesMap(Protocol):
    def __call__(self, func: Callable[[_S], _T], iterable: op.CanIter[op.CanNext[_S]], /) -> op.CanIter[op.CanIterSelf[_T]]: ...

@type_check_only
class _InfiniteFunc(Protocol[_NDT_co]):
    def get_t(self, /, x: float) -> float: ...
    def __call__(self, /, t: float) -> _NDT_co: ...

###

# undocumented
class LRUDict(collections.OrderedDict[tuple[float, float], _VT], Generic[_VT]):
    def __init__(self, /, max_size: int) -> None: ...
    @override
    def update(self, other: Never) -> NoReturn: ...  # type: ignore[override]  # pyright: ignore[reportIncompatibleMethodOverride]  # ty: ignore[invalid-method-override]

# undocumented
class SemiInfiniteFunc(_InfiniteFunc[_NDT_co], Generic[_NDT_co]):
    def __init__(self, /, func: Callable[[float], _NDT_co], start: float, infty: bool) -> None: ...

# undocumented
class DoubleInfiniteFunc(_InfiniteFunc, Generic[_NDT_co]):
    def __init__(self, /, func: Callable[[float], _NDT_co]) -> None: ...

# NOTE: This is only used as "info dict" for `quad_vec(..., full_output=True)`,
# even though, confusingly, it is not even even a mapping.
# NOTE: Because this "bunch" is only used as "info dict" (and nowhere else),
# its the ~keys~ attributes have been annotated right here.
class _Bunch(Generic[_InexactT_co, _ShapeT_co]):  # undocumented
    def __init__(
        self,
        /,
        *,
        success: bool,
        status: Literal[0, 1, 2],
        neval: int,
        message: str,
        intervals: onp.Array2D[np.float64],
        errors: onp.Array1D[np.float64],
        integrals: onp.ArrayND[_InexactT_co],
    ) -> None: ...
    success: Final[bool]
    status: Final[Literal[0, 1, 2]]
    neval: Final[int]
    message: Final[str]
    intervals: Final[onp.Array2D[np.float64]]
    errors: Final[onp.Array1D[np.float64]]
    integrals: onp.ArrayND[_InexactT_co, _ShapeT_co]

#
@overload  # 0d +integer
def quad_vec(
    f: _Fun[npc.integer | np.bool_],
    a: float,
    b: float,
    epsabs: float = 1e-200,
    epsrel: float = 1e-08,
    norm: _Norm = "2",
    cache_size: float = 100_000_000,
    limit: float = 10_000,
    workers: int | _DoesMap = 1,
    points: onp.ToFloat1D | None = None,
    quadrature: _Quadrature | None = None,
    full_output: Literal[False] = False,
    *,
    args: tuple[object, ...] = (),
) -> tuple[np.float64, float]: ...
@overload  # 0d +integer, full_output=True
def quad_vec(
    f: _Fun[npc.integer | np.bool_],
    a: float,
    b: float,
    epsabs: float = 1e-200,
    epsrel: float = 1e-08,
    norm: _Norm = "2",
    cache_size: float = 100_000_000,
    limit: float = 10_000,
    workers: int | _DoesMap = 1,
    points: onp.ToFloat1D | None = None,
    quadrature: _Quadrature | None = None,
    *,
    full_output: Literal[True],
    args: tuple[object, ...] = (),
) -> tuple[np.float64, float, _Bunch[np.float64, tuple[int]]]: ...
@overload  # Nd +integer
def quad_vec(
    f: _Fun[onp.ArrayND[npc.integer | np.bool_]],
    a: float,
    b: float,
    epsabs: float = 1e-200,
    epsrel: float = 1e-08,
    norm: _Norm = "2",
    cache_size: float = 100_000_000,
    limit: float = 10_000,
    workers: int | _DoesMap = 1,
    points: onp.ToFloat1D | None = None,
    quadrature: _Quadrature | None = None,
    full_output: Literal[False] = False,
    *,
    args: tuple[object, ...] = (),
) -> tuple[onp.Array1D[np.float64], float]: ...
@overload  # 1d +integer, full_output=True
def quad_vec(
    f: _Fun[onp.ArrayND[npc.integer | np.bool_]],
    a: float,
    b: float,
    epsabs: float = 1e-200,
    epsrel: float = 1e-08,
    norm: _Norm = "2",
    cache_size: float = 100_000_000,
    limit: float = 10_000,
    workers: int | _DoesMap = 1,
    points: onp.ToFloat1D | None = None,
    quadrature: _Quadrature | None = None,
    *,
    full_output: Literal[True],
    args: tuple[object, ...] = (),
) -> tuple[onp.Array1D[np.float64], float, _Bunch[np.float64, tuple[int, int]]]: ...
@overload  # 0d T:inexact
def quad_vec(
    f: _Fun[_InexactT],
    a: float,
    b: float,
    epsabs: float = 1e-200,
    epsrel: float = 1e-08,
    norm: _Norm = "2",
    cache_size: float = 100_000_000,
    limit: float = 10_000,
    workers: int | _DoesMap = 1,
    points: onp.ToFloat1D | None = None,
    quadrature: _Quadrature | None = None,
    full_output: Literal[False] = False,
    *,
    args: tuple[object, ...] = (),
) -> tuple[_InexactT, float]: ...
@overload  # 0d T:inexact, full_output=True
def quad_vec(
    f: _Fun[_InexactT],
    a: float,
    b: float,
    epsabs: float = 1e-200,
    epsrel: float = 1e-08,
    norm: _Norm = "2",
    cache_size: float = 100_000_000,
    limit: float = 10_000,
    workers: int | _DoesMap = 1,
    points: onp.ToFloat1D | None = None,
    quadrature: _Quadrature | None = None,
    *,
    full_output: Literal[True],
    args: tuple[object, ...] = (),
) -> tuple[_InexactT, float, _Bunch[_InexactT, tuple[int]]]: ...
@overload  # Nd T:inexact
def quad_vec(
    f: _Fun[onp.ArrayND[_InexactT]],
    a: float,
    b: float,
    epsabs: float = 1e-200,
    epsrel: float = 1e-08,
    norm: _Norm = "2",
    cache_size: float = 100_000_000,
    limit: float = 10_000,
    workers: int | _DoesMap = 1,
    points: onp.ToFloat1D | None = None,
    quadrature: _Quadrature | None = None,
    full_output: Literal[False] = False,
    *,
    args: tuple[object, ...] = (),
) -> tuple[onp.Array1D[_InexactT], float]: ...
@overload  # Nd T:inexact, full_output=True
def quad_vec(
    f: _Fun[onp.ArrayND[_InexactT]],
    a: float,
    b: float,
    epsabs: float = 1e-200,
    epsrel: float = 1e-08,
    norm: _Norm = "2",
    cache_size: float = 100_000_000,
    limit: float = 10_000,
    workers: int | _DoesMap = 1,
    points: onp.ToFloat1D | None = None,
    quadrature: _Quadrature | None = None,
    *,
    full_output: Literal[True],
    args: tuple[object, ...] = (),
) -> tuple[onp.Array1D[_InexactT], float, _Bunch[_InexactT, tuple[int, int]]]: ...
@overload  # 0d float
def quad_vec(
    f: _Fun[float],
    a: float,
    b: float,
    epsabs: float = 1e-200,
    epsrel: float = 1e-08,
    norm: _Norm = "2",
    cache_size: float = 100_000_000,
    limit: float = 10_000,
    workers: int | _DoesMap = 1,
    points: onp.ToFloat1D | None = None,
    quadrature: _Quadrature | None = None,
    full_output: Literal[False] = False,
    *,
    args: tuple[object, ...] = (),
) -> tuple[float, float]: ...
@overload  # 0d float, full_output=True
def quad_vec(
    f: _Fun[float],
    a: float,
    b: float,
    epsabs: float = 1e-200,
    epsrel: float = 1e-08,
    norm: _Norm = "2",
    cache_size: float = 100_000_000,
    limit: float = 10_000,
    workers: int | _DoesMap = 1,
    points: onp.ToFloat1D | None = None,
    quadrature: _Quadrature | None = None,
    *,
    full_output: Literal[True],
    args: tuple[object, ...] = (),
) -> tuple[float, float, _Bunch[np.float64, tuple[int]]]: ...
@overload  # 0d complex
def quad_vec(
    f: _Fun[complex],
    a: float,
    b: float,
    epsabs: float = 1e-200,
    epsrel: float = 1e-08,
    norm: _Norm = "2",
    cache_size: float = 100_000_000,
    limit: float = 10_000,
    workers: int | _DoesMap = 1,
    points: onp.ToFloat1D | None = None,
    quadrature: _Quadrature | None = None,
    full_output: Literal[False] = False,
    *,
    args: tuple[object, ...] = (),
) -> tuple[complex, float]: ...
@overload  # 0d complex, full_output=True
def quad_vec(
    f: _Fun[op.JustComplex],
    a: float,
    b: float,
    epsabs: float = 1e-200,
    epsrel: float = 1e-08,
    norm: _Norm = "2",
    cache_size: float = 100_000_000,
    limit: float = 10_000,
    workers: int | _DoesMap = 1,
    points: onp.ToFloat1D | None = None,
    quadrature: _Quadrature | None = None,
    *,
    full_output: Literal[True],
    args: tuple[object, ...] = (),
) -> tuple[complex, float, _Bunch[np.complex128, tuple[int]]]: ...
