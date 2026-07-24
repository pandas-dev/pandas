import collections
from collections.abc import Callable
from typing import Any, Concatenate, Final, Generic, Literal, Never, NoReturn, Protocol, overload, override, type_check_only
from typing_extensions import TypeVar

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

###

type _Floating = float | npc.floating
type _FloatingND = onp.ArrayND[npc.floating] | _Floating

type _Fun[T] = Callable[Concatenate[float, ...], T]

type _Norm = Literal["max", "2"]
type _Quadrature = Literal["gk21", "gk15", "trapezoid"]

# workaround for mypy & pyright's failure to conform to the overload typing specification
type _JustAnyShape = tuple[Never, Never, Never]

_VT = TypeVar("_VT", default=Any)
_NDT_co = TypeVar("_NDT_co", bound=_FloatingND, default=_FloatingND, covariant=True)
_InexactT_co = TypeVar("_InexactT_co", bound=npc.inexact, default=Any, covariant=True)
_ShapeT_co = TypeVar("_ShapeT_co", bound=tuple[int, ...], default=tuple[Any, ...], covariant=True)

@type_check_only
class _DoesMap(Protocol):
    def __call__[S, T](self, func: Callable[[S], T], iterable: op.CanIter[op.CanNext[S]], /) -> op.CanIter[op.CanIterSelf[T]]: ...

@type_check_only
class _InfiniteFunc(Protocol[_NDT_co]):
    def get_t(self, /, x: float) -> float: ...
    def __call__(self, /, t: float) -> _NDT_co: ...

###

# undocumented
class LRUDict(collections.OrderedDict[tuple[float, float], _VT], Generic[_VT]):
    def __init__(self, /, max_size: int) -> None: ...
    @override
    # pyrefly: ignore [bad-override]
    def update(self, other: Never, /, **kwargs: Never) -> NoReturn: ...  # type: ignore[override]  # pyright: ignore[reportIncompatibleMethodOverride]  # ty: ignore[invalid-method-override]

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

# NOTE: Unlike `fixed_quad`, `quad_vec` evaluates `f` at scalar points, so the output
# shape equals the shape of the returned array of `f` (no axis is reduced), and dtypes
# below 64 bits precision are not upcast (the quadrature weights are python floats).
# The `integrals` array of the "info dict" has one additional leading axis.

#
@overload  # 0d +integer
def quad_vec(
    f: _Fun[npc.integer | np.bool],
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
    f: _Fun[npc.integer | np.bool],
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
@overload  # ?d +integer  (mypy & pyright workaround)
def quad_vec(
    f: _Fun[onp.ArrayND[npc.integer | np.bool, _JustAnyShape]],
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
) -> tuple[onp.ArrayND[np.float64], float]: ...
@overload  # ?d +integer, full_output=True  (mypy & pyright workaround)
def quad_vec(
    f: _Fun[onp.ArrayND[npc.integer | np.bool, _JustAnyShape]],
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
) -> tuple[onp.ArrayND[np.float64], float, _Bunch[np.float64]]: ...
@overload  # 1d +integer
def quad_vec(
    f: _Fun[onp.Array1D[npc.integer | np.bool]],
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
    f: _Fun[onp.Array1D[npc.integer | np.bool]],
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
@overload  # 2d +integer
def quad_vec(
    f: _Fun[onp.Array2D[npc.integer | np.bool]],
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
) -> tuple[onp.Array2D[np.float64], float]: ...
@overload  # 2d +integer, full_output=True
def quad_vec(
    f: _Fun[onp.Array2D[npc.integer | np.bool]],
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
) -> tuple[onp.Array2D[np.float64], float, _Bunch[np.float64, tuple[int, int, int]]]: ...
@overload  # Nd +integer
def quad_vec(
    f: _Fun[onp.ArrayND[npc.integer | np.bool]],
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
) -> tuple[onp.ArrayND[np.float64], float]: ...
@overload  # Nd +integer, full_output=True
def quad_vec(
    f: _Fun[onp.ArrayND[npc.integer | np.bool]],
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
) -> tuple[onp.ArrayND[np.float64], float, _Bunch[np.float64]]: ...
@overload  # 0d T:inexact
def quad_vec[InexactT: npc.inexact](
    f: _Fun[InexactT],
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
) -> tuple[InexactT, float]: ...
@overload  # 0d T:inexact, full_output=True
def quad_vec[InexactT: npc.inexact](
    f: _Fun[InexactT],
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
) -> tuple[InexactT, float, _Bunch[InexactT, tuple[int]]]: ...
@overload  # ?d T:inexact  (mypy & pyright workaround)
def quad_vec[InexactT: npc.inexact](
    f: _Fun[onp.ArrayND[InexactT, _JustAnyShape]],
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
) -> tuple[onp.ArrayND[InexactT], float]: ...
@overload  # ?d T:inexact, full_output=True  (mypy & pyright workaround)
def quad_vec[InexactT: npc.inexact](
    f: _Fun[onp.ArrayND[InexactT, _JustAnyShape]],
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
) -> tuple[onp.ArrayND[InexactT], float, _Bunch[InexactT]]: ...
@overload  # 1d T:inexact
def quad_vec[InexactT: npc.inexact](
    f: _Fun[onp.Array1D[InexactT]],
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
) -> tuple[onp.Array1D[InexactT], float]: ...
@overload  # 1d T:inexact, full_output=True
def quad_vec[InexactT: npc.inexact](
    f: _Fun[onp.Array1D[InexactT]],
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
) -> tuple[onp.Array1D[InexactT], float, _Bunch[InexactT, tuple[int, int]]]: ...
@overload  # 2d T:inexact
def quad_vec[InexactT: npc.inexact](
    f: _Fun[onp.Array2D[InexactT]],
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
) -> tuple[onp.Array2D[InexactT], float]: ...
@overload  # 2d T:inexact, full_output=True
def quad_vec[InexactT: npc.inexact](
    f: _Fun[onp.Array2D[InexactT]],
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
) -> tuple[onp.Array2D[InexactT], float, _Bunch[InexactT, tuple[int, int, int]]]: ...
@overload  # Nd T:inexact
def quad_vec[InexactT: npc.inexact](
    f: _Fun[onp.ArrayND[InexactT]],
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
) -> tuple[onp.ArrayND[InexactT], float]: ...
@overload  # Nd T:inexact, full_output=True
def quad_vec[InexactT: npc.inexact](
    f: _Fun[onp.ArrayND[InexactT]],
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
) -> tuple[onp.ArrayND[InexactT], float, _Bunch[InexactT]]: ...
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
