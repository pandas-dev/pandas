from collections.abc import Callable
from typing import Any, Concatenate, Literal, NamedTuple, Never, TypeAlias, overload
from typing_extensions import TypeVar

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.stats.qmc import QMCEngine

__all__ = ["cumulative_simpson", "cumulative_trapezoid", "fixed_quad", "newton_cotes", "qmc_quad", "romb", "simpson", "trapezoid"]

_T = TypeVar("_T")
_InexactT = TypeVar("_InexactT", bound=npc.inexact)
_Inexact80T = TypeVar("_Inexact80T", bound=npc.inexact80)
_ShapeT = TypeVar("_ShapeT", bound=tuple[Any, ...])

_FixedQuadFunc: TypeAlias = Callable[Concatenate[onp.Array1D[np.float64], ...], _T]

# workaround for mypy & pyright's failure to conform to the overload typing specification
_JustAnyShape: TypeAlias = tuple[Never, Never, Never]

###

class QMCQuadResult(NamedTuple):
    integral: np.float64
    standard_error: np.float64

# sample-based integration

#
@overload  # ?d +complex  (mypy & pyright workaround)
def trapezoid(
    y: onp.Array[_JustAnyShape, npc.number | np.bool_], x: onp.ToFloatND | None = None, dx: float = 1.0, axis: int = -1
) -> Any: ...
@overload  # 1d T:inexact
def trapezoid(y: onp.Array1D[_InexactT], x: onp.ToFloat1D | None = None, dx: float = 1.0, axis: int = -1) -> _InexactT: ...
@overload  # 1d +int
def trapezoid(
    y: onp.ToArrayStrict1D[float, npc.integer | np.bool_], x: onp.ToFloat1D | None = None, dx: float = 1.0, axis: int = -1
) -> np.float64: ...
@overload  # 1d ~complex
def trapezoid(
    y: onp.ToJustComplex128Strict1D, x: onp.ToFloat1D | None = None, dx: float = 1.0, axis: int = -1
) -> np.complex128: ...
@overload  # 2d T:inexact
def trapezoid(
    y: onp.Array2D[_InexactT], x: onp.ToFloatND | None = None, dx: float = 1.0, axis: int = -1
) -> onp.Array1D[_InexactT]: ...
@overload  # 2d +int
def trapezoid(
    y: onp.ToArrayStrict2D[float, npc.integer | np.bool_], x: onp.ToFloatND | None = None, dx: float = 1.0, axis: int = -1
) -> onp.Array1D[np.float64]: ...
@overload  # 2d ~complex
def trapezoid(
    y: onp.ToJustComplex128Strict2D, x: onp.ToFloatND | None = None, dx: float = 1.0, axis: int = -1
) -> onp.Array1D[np.complex128]: ...
@overload  # Nd T:inexact
def trapezoid(
    y: onp.ArrayND[_InexactT], x: onp.ToFloatND | None = None, dx: float = 1.0, axis: int = -1
) -> onp.ArrayND[_InexactT] | Any: ...
@overload  # Nd +int
def trapezoid(
    y: onp.ToArrayND[float, npc.integer | np.bool_], x: onp.ToFloatND | None = None, dx: float = 1.0, axis: int = -1
) -> onp.ArrayND[np.float64] | Any: ...
@overload  # Nd ~complex
def trapezoid(
    y: onp.ToJustComplex128_ND, x: onp.ToFloatND | None = None, dx: float = 1.0, axis: int = -1
) -> onp.ArrayND[np.complex128] | Any: ...
@overload  # +float (fallback)
def trapezoid(
    y: onp.ToFloatND, x: onp.ToFloatND | None = None, dx: float = 1.0, axis: int = -1
) -> onp.ArrayND[np.float64 | Any] | Any: ...
@overload  # +complex (fallback)
def trapezoid(
    y: onp.ToComplexND, x: onp.ToFloatND | None = None, dx: float = 1.0, axis: int = -1
) -> onp.ArrayND[np.complex128 | Any] | Any: ...

# NOTE: unlike `trapezoid`, this will upcast scalars below 64-bits precision in case of scalar output
@overload  # ?d +complex  (mypy & pyright workaround)
def simpson(
    y: onp.Array[_JustAnyShape, npc.number | np.bool_], x: onp.ToFloatND | None = None, *, dx: float = 1.0, axis: int = -1
) -> Any: ...
@overload  # 1d +f64
def simpson(y: onp.ToFloat64Strict1D, x: onp.ToFloatND | None = None, *, dx: float = 1.0, axis: int = -1) -> np.float64: ...
@overload  # 1d ~complex
def simpson(
    y: onp.ToJustComplex128Strict1D | onp.ToJustComplex64Strict1D,
    x: onp.ToFloatND | None = None,
    *,
    dx: float = 1.0,
    axis: int = -1,
) -> np.complex128: ...
@overload  # 1d T:inexact
def simpson(y: onp.Array1D[_Inexact80T], x: onp.ToFloatND | None = None, *, dx: float = 1.0, axis: int = -1) -> _Inexact80T: ...
@overload  # 2d T:inexact
def simpson(
    y: onp.Array2D[_InexactT], x: onp.ToFloatND | None = None, *, dx: float = 1.0, axis: int = -1
) -> onp.Array1D[_InexactT]: ...
@overload  # 2d +int
def simpson(
    y: onp.ToArrayStrict2D[float, npc.integer | np.bool_], x: onp.ToFloatND | None = None, *, dx: float = 1.0, axis: int = -1
) -> onp.Array1D[np.float64]: ...
@overload  # 2d ~complex
def simpson(
    y: onp.ToJustComplex128Strict2D, x: onp.ToFloatND | None = None, *, dx: float = 1.0, axis: int = -1
) -> onp.Array1D[np.complex128]: ...
@overload  # Nd T:inexact
def simpson(
    y: onp.ArrayND[_InexactT], x: onp.ToFloatND | None = None, *, dx: float = 1.0, axis: int = -1
) -> onp.ArrayND[_InexactT] | Any: ...
@overload  # Nd +int
def simpson(
    y: onp.ToArrayND[float, npc.integer | np.bool_], x: onp.ToFloatND | None = None, *, dx: float = 1.0, axis: int = -1
) -> onp.ArrayND[np.float64] | Any: ...
@overload  # Nd ~complex
def simpson(
    y: onp.ToJustComplex128_ND, x: onp.ToFloatND | None = None, *, dx: float = 1.0, axis: int = -1
) -> onp.ArrayND[np.complex128] | Any: ...
@overload  # +float (fallback)
def simpson(
    y: onp.ToFloatND, x: onp.ToFloatND | None = None, *, dx: float = 1.0, axis: int = -1
) -> onp.ArrayND[np.float64 | Any] | Any: ...
@overload  # +complex (fallback)
def simpson(
    y: onp.ToComplexND, x: onp.ToFloatND | None = None, *, dx: float = 1.0, axis: int = -1
) -> onp.ArrayND[np.complex128 | Any] | Any: ...

# NOTE: like `simpson`, but this will also upcast scalars below 64-bits precision in case of array output
@overload  # ?d +complex  (mypy & pyright workaround)
def romb(y: onp.Array[_JustAnyShape, npc.number | np.bool_], dx: float = 1.0, axis: int = -1, show: bool = False) -> Any: ...
@overload  # 1d +f64
def romb(y: onp.ToFloat64Strict1D, dx: float = 1.0, axis: int = -1) -> np.float64: ...
@overload  # 1d ~complex
def romb(
    y: onp.ToJustComplex128Strict1D | onp.ToJustComplex64Strict1D, dx: float = 1.0, axis: int = -1, show: bool = False
) -> np.complex128: ...
@overload  # 1d T:inexact
def romb(y: onp.Array1D[_Inexact80T], dx: float = 1.0, axis: int = -1, show: bool = False) -> _Inexact80T: ...
@overload  # 2d +int
def romb(y: onp.ToFloat64Strict2D, dx: float = 1.0, axis: int = -1, show: bool = False) -> onp.Array1D[np.float64]: ...
@overload  # 2d ~complex
def romb(
    y: onp.ToJustComplex128Strict2D | onp.ToJustComplex64Strict2D, dx: float = 1.0, axis: int = -1, show: bool = False
) -> onp.Array1D[np.complex128]: ...
@overload  # 2d T:inexact
def romb(y: onp.Array2D[_Inexact80T], dx: float = 1.0, axis: int = -1, show: bool = False) -> onp.Array1D[_Inexact80T]: ...
@overload  # Nd +int
def romb(y: onp.ToFloat64_ND, dx: float = 1.0, axis: int = -1, show: bool = False) -> onp.ArrayND[np.float64] | Any: ...
@overload  # Nd ~complex
def romb(
    y: onp.ToJustComplex128_ND | onp.ToJustComplex64_ND, dx: float = 1.0, axis: int = -1, show: bool = False
) -> onp.ArrayND[np.complex128] | Any: ...
@overload  # Nd T:inexact
def romb(y: onp.ArrayND[_Inexact80T], dx: float = 1.0, axis: int = -1, show: bool = False) -> onp.ArrayND[_Inexact80T] | Any: ...
@overload  # +float (fallback)
def romb(y: onp.ToFloatND, dx: float = 1.0, axis: int = -1, show: bool = False) -> onp.ArrayND[np.float64 | Any] | Any: ...
@overload  # +complex (fallback)
def romb(y: onp.ToComplexND, dx: float = 1.0, axis: int = -1, show: bool = False) -> onp.ArrayND[np.complex128 | Any] | Any: ...

# sample-based cumulative integration

# keep in sync with `cumulative_simpson`
@overload  # +int, shape known
def cumulative_trapezoid(
    y: onp.ArrayND[npc.integer | np.bool_, _ShapeT],
    x: onp.ToFloatND | None = None,
    dx: float = 1.0,
    axis: int = -1,
    initial: Literal[0] | None = None,
) -> onp.ArrayND[np.float64, _ShapeT]: ...
@overload  # +float, shape 1d
def cumulative_trapezoid(
    y: onp.ToArrayStrict1D[float, npc.integer | np.bool_],
    x: onp.ToFloatND | None = None,
    dx: float = 1.0,
    axis: int = -1,
    initial: Literal[0] | None = None,
) -> onp.Array1D[np.float64]: ...
@overload  # +float, shape 2d
def cumulative_trapezoid(
    y: onp.ToArrayStrict2D[float, npc.integer | np.bool_],
    x: onp.ToFloatND | None = None,
    dx: float = 1.0,
    axis: int = -1,
    initial: Literal[0] | None = None,
) -> onp.Array2D[np.float64]: ...
@overload  # +float, shape unknown
def cumulative_trapezoid(
    y: onp.SequenceND[float], x: onp.ToFloatND | None = None, dx: float = 1.0, axis: int = -1, initial: Literal[0] | None = None
) -> onp.ArrayND[np.float64]: ...
@overload  # T:inexact, shape known
def cumulative_trapezoid(
    y: onp.ArrayND[_InexactT, _ShapeT],
    x: onp.ToFloatND | None = None,
    dx: float = 1.0,
    axis: int = -1,
    initial: Literal[0] | None = None,
) -> onp.ArrayND[_InexactT, _ShapeT]: ...
@overload  # ~complex, shape 1d
def cumulative_trapezoid(
    y: onp.ToJustComplex128Strict1D,
    x: onp.ToFloatND | None = None,
    dx: float = 1.0,
    axis: int = -1,
    initial: Literal[0] | None = None,
) -> onp.Array1D[np.complex128]: ...
@overload  # ~complex, shape 2d
def cumulative_trapezoid(
    y: onp.ToJustComplex128Strict2D,
    x: onp.ToFloatND | None = None,
    dx: float = 1.0,
    axis: int = -1,
    initial: Literal[0] | None = None,
) -> onp.Array2D[np.complex128]: ...
@overload  # ~complex, shape unknown
def cumulative_trapezoid(
    y: onp.ToJustComplex128_ND, x: onp.ToFloatND | None = None, dx: float = 1.0, axis: int = -1, initial: Literal[0] | None = None
) -> onp.ArrayND[np.complex128]: ...
@overload  # fallback
def cumulative_trapezoid(
    y: onp.ToJustFloatND | onp.ToJustComplexND,  # `ToComplexND` would overlap with the first overload
    x: onp.ToFloatND | None = None,
    dx: float = 1.0,
    axis: int = -1,
    initial: Literal[0] | None = None,
) -> onp.Array: ...

# keep in sync with `cumulative_trapezoid`
@overload  # +int, shape known
def cumulative_simpson(
    y: onp.ArrayND[npc.integer | np.bool_, _ShapeT],
    *,
    x: onp.ToFloatND | None = None,
    dx: float = 1.0,
    axis: int = -1,
    initial: float | onp.ToFloatND | None = None,
) -> onp.ArrayND[np.float64, _ShapeT]: ...
@overload  # +float, shape 1d
def cumulative_simpson(
    y: onp.ToArrayStrict1D[float, npc.integer | np.bool_],
    *,
    x: onp.ToFloatND | None = None,
    dx: float = 1.0,
    axis: int = -1,
    initial: float | onp.ToFloatND | None = None,
) -> onp.Array1D[np.float64]: ...
@overload  # +float, shape 2d
def cumulative_simpson(
    y: onp.ToArrayStrict2D[float, npc.integer | np.bool_],
    *,
    x: onp.ToFloatND | None = None,
    dx: float = 1.0,
    axis: int = -1,
    initial: float | onp.ToFloatND | None = None,
) -> onp.Array2D[np.float64]: ...
@overload  # +float, shape unknown
def cumulative_simpson(
    y: onp.SequenceND[float],
    *,
    x: onp.ToFloatND | None = None,
    dx: float = 1.0,
    axis: int = -1,
    initial: float | onp.ToFloatND | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload  # T:inexact, shape known
def cumulative_simpson(
    y: onp.ArrayND[_InexactT, _ShapeT],
    *,
    x: onp.ToFloatND | None = None,
    dx: float = 1.0,
    axis: int = -1,
    initial: float | onp.ToFloatND | None = None,
) -> onp.ArrayND[_InexactT, _ShapeT]: ...
@overload  # ~complex, shape 1d
def cumulative_simpson(
    y: onp.ToJustComplex128Strict1D,
    *,
    x: onp.ToFloatND | None = None,
    dx: float = 1.0,
    axis: int = -1,
    initial: float | onp.ToFloatND | None = None,
) -> onp.Array1D[np.complex128]: ...
@overload  # ~complex, shape 2d
def cumulative_simpson(
    y: onp.ToJustComplex128Strict2D,
    *,
    x: onp.ToFloatND | None = None,
    dx: float = 1.0,
    axis: int = -1,
    initial: float | onp.ToFloatND | None = None,
) -> onp.Array2D[np.complex128]: ...
@overload  # ~complex, shape unknown
def cumulative_simpson(
    y: onp.ToJustComplex128_ND,
    *,
    x: onp.ToFloatND | None = None,
    dx: float = 1.0,
    axis: int = -1,
    initial: float | onp.ToFloatND | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload  # fallback
def cumulative_simpson(
    y: onp.ToJustFloatND | onp.ToJustComplexND,  # `ToComplexND` would overlap with the first overload
    *,
    x: onp.ToFloatND | None = None,
    dx: float = 1.0,
    axis: int = -1,
    initial: float | onp.ToFloatND | None = None,
) -> onp.Array: ...

# function-based

#
@overload  # (1d f64) -> +f64
def fixed_quad(
    func: _FixedQuadFunc[onp.ToFloat64_ND], a: float, b: float, args: tuple[object, ...] = (), n: int = 5
) -> tuple[np.float64, None]: ...
@overload  # (1d f64) -> ~c64 | ~c128
def fixed_quad(
    func: _FixedQuadFunc[onp.ToArray1D[op.JustComplex, np.complex64 | np.complex128]],
    a: float,
    b: float,
    args: tuple[object, ...] = (),
    n: int = 5,
) -> tuple[np.complex128, None]: ...
@overload  # (1d f64) -> ~f80 | ~c160
def fixed_quad(
    func: _FixedQuadFunc[onp.ArrayND[_Inexact80T]], a: float, b: float, args: tuple[object, ...] = (), n: int = 5
) -> tuple[_Inexact80T, None]: ...
@overload  # (1d f64) -> ~m64
def fixed_quad(
    func: _FixedQuadFunc[onp.ArrayND[np.timedelta64]], a: float, b: float, args: tuple[object, ...] = (), n: int = 5
) -> tuple[np.timedelta64, None]: ...
@overload  # (1d f64) -> ~object_
def fixed_quad(
    func: _FixedQuadFunc[onp.ArrayND[np.object_]], a: float, b: float, args: tuple[object, ...] = (), n: int = 5
) -> tuple[Any, None]: ...

#
def qmc_quad(
    func: Callable[[onp.ArrayND[np.float64]], onp.ArrayND[npc.number]],
    a: float | onp.ToFloat1D,
    b: float | onp.ToFloat1D,
    *,
    n_estimates: int = 8,
    n_points: int = 1_024,
    qrng: QMCEngine | None = None,
    log: bool = False,
) -> QMCQuadResult: ...

# low-level
def newton_cotes(rn: int, equal: int = 0) -> tuple[onp.Array1D[np.float64], float]: ...
