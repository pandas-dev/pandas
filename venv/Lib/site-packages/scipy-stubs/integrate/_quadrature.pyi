from collections.abc import Callable
from typing import Any, Concatenate, Literal, NamedTuple, Never, overload

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.stats.qmc import QMCEngine

__all__ = ["cumulative_simpson", "cumulative_trapezoid", "fixed_quad", "newton_cotes", "qmc_quad", "romb", "simpson", "trapezoid"]

###

type _FixedQuadFunc[T] = Callable[Concatenate[onp.Array1D[np.float64], ...], T]

# workaround for mypy & pyright's failure to conform to the overload typing specification
type _JustAnyShape = tuple[Never, Never, Never]

###

# mypy reports false positive `overload-overlap` errors for `cumulative_simpson` only on `numpy<2.1`
# mypy: disable-error-code=overload-overlap

class QMCQuadResult(NamedTuple):
    integral: np.float64
    standard_error: np.float64

# sample-based integration

#
@overload  # ?d +complex  (mypy & pyright workaround)
def trapezoid(
    y: onp.Array[_JustAnyShape, npc.number | np.bool], x: onp.ToFloatND | None = None, dx: float = 1.0, axis: int = -1
) -> Any: ...
@overload  # 1d T:inexact
def trapezoid[InexactT: npc.inexact](
    y: onp.Array1D[InexactT], x: onp.ToFloat1D | None = None, dx: float = 1.0, axis: int = -1
) -> InexactT: ...
@overload  # 1d +int
def trapezoid(
    y: onp.ToArrayStrict1D[float, npc.integer | np.bool], x: onp.ToFloat1D | None = None, dx: float = 1.0, axis: int = -1
) -> np.float64: ...
@overload  # 1d ~complex
def trapezoid(
    y: onp.ToJustComplex128Strict1D, x: onp.ToFloat1D | None = None, dx: float = 1.0, axis: int = -1
) -> np.complex128: ...
@overload  # 2d T:inexact
def trapezoid[InexactT: npc.inexact](
    y: onp.Array2D[InexactT], x: onp.ToFloatND | None = None, dx: float = 1.0, axis: int = -1
) -> onp.Array1D[InexactT]: ...
@overload  # 2d +int
def trapezoid(
    y: onp.ToArrayStrict2D[float, npc.integer | np.bool], x: onp.ToFloatND | None = None, dx: float = 1.0, axis: int = -1
) -> onp.Array1D[np.float64]: ...
@overload  # 2d ~complex
def trapezoid(
    y: onp.ToJustComplex128Strict2D, x: onp.ToFloatND | None = None, dx: float = 1.0, axis: int = -1
) -> onp.Array1D[np.complex128]: ...
@overload  # Nd T:inexact
def trapezoid[InexactT: npc.inexact](
    y: onp.ArrayND[InexactT], x: onp.ToFloatND | None = None, dx: float = 1.0, axis: int = -1
) -> onp.ArrayND[InexactT] | Any: ...
@overload  # Nd +int
def trapezoid(
    y: onp.ToArrayND[float, npc.integer | np.bool], x: onp.ToFloatND | None = None, dx: float = 1.0, axis: int = -1
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

#
@overload  # ?d +complex  (mypy & pyright workaround)
def simpson(
    y: onp.Array[_JustAnyShape, npc.number | np.bool], x: onp.ToFloatND | None = None, *, dx: float = 1.0, axis: int = -1
) -> Any: ...
@overload  # 1d T:inexact
def simpson[InexactT: npc.inexact](
    y: onp.Array1D[InexactT], x: onp.ToFloatND | None = None, *, dx: float = 1.0, axis: int = -1
) -> InexactT: ...
@overload  # 1d +int
def simpson(
    y: onp.ToArrayStrict1D[float, npc.integer | np.bool], x: onp.ToFloatND | None = None, *, dx: float = 1.0, axis: int = -1
) -> np.float64: ...
@overload  # 1d ~complex
def simpson(
    y: onp.ToJustComplex128Strict1D, x: onp.ToFloatND | None = None, *, dx: float = 1.0, axis: int = -1
) -> np.complex128: ...
@overload  # 2d T:inexact
def simpson[InexactT: npc.inexact](
    y: onp.Array2D[InexactT], x: onp.ToFloatND | None = None, *, dx: float = 1.0, axis: int = -1
) -> onp.Array1D[InexactT]: ...
@overload  # 2d +int
def simpson(
    y: onp.ToArrayStrict2D[float, npc.integer | np.bool], x: onp.ToFloatND | None = None, *, dx: float = 1.0, axis: int = -1
) -> onp.Array1D[np.float64]: ...
@overload  # 2d ~complex
def simpson(
    y: onp.ToJustComplex128Strict2D, x: onp.ToFloatND | None = None, *, dx: float = 1.0, axis: int = -1
) -> onp.Array1D[np.complex128]: ...
@overload  # Nd T:inexact
def simpson[InexactT: npc.inexact](
    y: onp.ArrayND[InexactT], x: onp.ToFloatND | None = None, *, dx: float = 1.0, axis: int = -1
) -> onp.ArrayND[InexactT] | Any: ...
@overload  # Nd +int
def simpson(
    y: onp.ToArrayND[float, npc.integer | np.bool], x: onp.ToFloatND | None = None, *, dx: float = 1.0, axis: int = -1
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

# NOTE: unlike `simpson`, this upcasts sub-64-bits sctypes
@overload  # ?d +complex  (mypy & pyright workaround)
def romb(y: onp.Array[_JustAnyShape, npc.number | np.bool], dx: float = 1.0, axis: int = -1, show: bool = False) -> Any: ...
@overload  # 1d +f64
def romb(y: onp.ToFloat64Strict1D, dx: float = 1.0, axis: int = -1) -> np.float64: ...
@overload  # 1d ~complex
def romb(
    y: onp.ToJustComplex128Strict1D | onp.ToJustComplex64Strict1D, dx: float = 1.0, axis: int = -1, show: bool = False
) -> np.complex128: ...
@overload  # 1d T:inexact
def romb[Inexact80T: npc.inexact80](
    y: onp.Array1D[Inexact80T], dx: float = 1.0, axis: int = -1, show: bool = False
) -> Inexact80T: ...
@overload  # 2d +int
def romb(y: onp.ToFloat64Strict2D, dx: float = 1.0, axis: int = -1, show: bool = False) -> onp.Array1D[np.float64]: ...
@overload  # 2d ~complex
def romb(
    y: onp.ToJustComplex128Strict2D | onp.ToJustComplex64Strict2D, dx: float = 1.0, axis: int = -1, show: bool = False
) -> onp.Array1D[np.complex128]: ...
@overload  # 2d T:inexact
def romb[Inexact80T: npc.inexact80](
    y: onp.Array2D[Inexact80T], dx: float = 1.0, axis: int = -1, show: bool = False
) -> onp.Array1D[Inexact80T]: ...
@overload  # Nd +int
def romb(y: onp.ToFloat64_ND, dx: float = 1.0, axis: int = -1, show: bool = False) -> onp.ArrayND[np.float64] | Any: ...
@overload  # Nd ~complex
def romb(
    y: onp.ToJustComplex128_ND | onp.ToJustComplex64_ND, dx: float = 1.0, axis: int = -1, show: bool = False
) -> onp.ArrayND[np.complex128] | Any: ...
@overload  # Nd T:inexact
def romb[Inexact80T: npc.inexact80](
    y: onp.ArrayND[Inexact80T], dx: float = 1.0, axis: int = -1, show: bool = False
) -> onp.ArrayND[Inexact80T] | Any: ...
@overload  # +float (fallback)
def romb(y: onp.ToFloatND, dx: float = 1.0, axis: int = -1, show: bool = False) -> onp.ArrayND[np.float64 | Any] | Any: ...
@overload  # +complex (fallback)
def romb(y: onp.ToComplexND, dx: float = 1.0, axis: int = -1, show: bool = False) -> onp.ArrayND[np.complex128 | Any] | Any: ...

# sample-based cumulative integration

# NOTE: unlike `cumulative_simpson`, this is dtype-preserving for dtypes below 64-bits precision
@overload  # +int, shape known
def cumulative_trapezoid[ShapeT: tuple[int, ...]](
    y: onp.ArrayND[npc.integer | np.bool, ShapeT],
    x: onp.ToFloatND | None = None,
    dx: float = 1.0,
    axis: int = -1,
    initial: Literal[0] | None = None,
) -> onp.ArrayND[np.float64, ShapeT]: ...
@overload  # +float, shape 1d
def cumulative_trapezoid(
    y: onp.ToArrayStrict1D[float, npc.integer | np.bool],
    x: onp.ToFloatND | None = None,
    dx: float = 1.0,
    axis: int = -1,
    initial: Literal[0] | None = None,
) -> onp.Array1D[np.float64]: ...
@overload  # +float, shape 2d
def cumulative_trapezoid(
    y: onp.ToArrayStrict2D[float, npc.integer | np.bool],
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
def cumulative_trapezoid[InexactT: npc.inexact, ShapeT: tuple[int, ...]](
    y: onp.ArrayND[InexactT, ShapeT],
    x: onp.ToFloatND | None = None,
    dx: float = 1.0,
    axis: int = -1,
    initial: Literal[0] | None = None,
) -> onp.ArrayND[InexactT, ShapeT]: ...
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

# NOTE: unlike `cumulative_trapezoid`, propagates 64bit sctypes if it matches `x`
@overload  # +int, shape known
def cumulative_simpson[ShapeT: tuple[int, ...]](
    y: onp.ArrayND[npc.integer | np.bool, ShapeT],
    *,
    x: onp.ToFloatND | None = None,
    dx: float = 1.0,
    axis: int = -1,
    initial: float | onp.ToFloatND | None = None,
) -> onp.ArrayND[np.float64, ShapeT]: ...
@overload  # +float, shape 1d
def cumulative_simpson(
    y: onp.ToArrayStrict1D[float, npc.integer | np.bool],
    *,
    x: onp.ToFloatND | None = None,
    dx: float = 1.0,
    axis: int = -1,
    initial: float | onp.ToFloatND | None = None,
) -> onp.Array1D[np.float64]: ...
@overload  # +float, shape 2d
def cumulative_simpson(
    y: onp.ToArrayStrict2D[float, npc.integer | np.bool],
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
@overload  # T:inexact64, shape known
def cumulative_simpson[InexactT: npc.inexact64 | npc.inexact80, ShapeT: tuple[int, ...]](
    y: onp.ArrayND[InexactT, ShapeT],
    *,
    x: onp.ToFloatND | None = None,
    dx: float = 1.0,
    axis: int = -1,
    initial: float | onp.ToFloatND | None = None,
) -> onp.ArrayND[InexactT, ShapeT]: ...
@overload  # T:float16/32, shape known, matching x
def cumulative_simpson[FloatT: (np.float16, np.float32), ShapeT: tuple[int, ...]](
    y: onp.ArrayND[FloatT, ShapeT],
    *,
    x: onp.ArrayND[FloatT],
    dx: float = 1.0,
    axis: int = -1,
    initial: FloatT | onp.ArrayND[FloatT] | None = None,
) -> onp.ArrayND[FloatT, ShapeT]: ...
@overload  # ~c64, shape known, narrow x
def cumulative_simpson[ShapeT: tuple[int, ...]](
    y: onp.ArrayND[np.complex64, ShapeT],
    *,
    x: onp.ArrayND[np.float32 | np.float16],
    dx: float = 1.0,
    axis: int = -1,
    initial: np.complex64 | onp.ArrayND[np.complex64] | None = None,
) -> onp.ArrayND[np.complex64, ShapeT]: ...
@overload  # ~f16/f32, shape known, promoting x
def cumulative_simpson[ShapeT: tuple[int, ...]](
    y: onp.ArrayND[np.float16 | np.float32, ShapeT],
    *,
    x: onp.ToJustFloat64_ND | None = None,
    dx: float = 1.0,
    axis: int = -1,
    initial: float | onp.ToFloatND | None = None,
) -> onp.ArrayND[np.float64, ShapeT]: ...
@overload  # ~c64, shape known, promoting x
def cumulative_simpson[ShapeT: tuple[int, ...]](
    y: onp.ArrayND[np.complex64, ShapeT],
    *,
    x: onp.ToJustFloat64_ND | None = None,
    dx: float = 1.0,
    axis: int = -1,
    initial: float | onp.ToFloatND | None = None,
) -> onp.ArrayND[np.complex128, ShapeT]: ...
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
    y: onp.ToComplexND,
    *,
    x: onp.ToFloatND | None = None,
    dx: float = 1.0,
    axis: int = -1,
    initial: float | onp.ToFloatND | None = None,
) -> onp.Array[tuple[Any, ...] | tuple[int]]: ...  # workaround for pyright on numpy<2.1

# function-based

#
@overload  # (?d f64) -> ?d f64  (mypy & pyright workaround)
def fixed_quad(
    func: _FixedQuadFunc[onp.ArrayND[np.float64 | np.float32 | np.float16 | npc.integer | np.bool, _JustAnyShape]],
    a: float,
    b: float,
    args: tuple[object, ...] = (),
    n: int = 5,
) -> tuple[np.float64 | onp.ArrayND[np.float64], None]: ...
@overload  # (1d f64) -> 1d +f64
def fixed_quad(
    func: _FixedQuadFunc[onp.ToFloat64Strict1D], a: float, b: float, args: tuple[object, ...] = (), n: int = 5
) -> tuple[np.float64, None]: ...
@overload  # (1d f64) -> 2d +f64
def fixed_quad(
    func: _FixedQuadFunc[onp.ToFloat64Strict2D], a: float, b: float, args: tuple[object, ...] = (), n: int = 5
) -> tuple[onp.Array1D[np.float64], None]: ...
@overload  # (1d f64) -> Nd f64
def fixed_quad(
    func: _FixedQuadFunc[onp.ToFloat64_ND], a: float, b: float, args: tuple[object, ...] = (), n: int = 5
) -> tuple[onp.ArrayND[np.float64] | Any, None]: ...
@overload  # (1d f64) -> ?d ~c64 | ~c128
def fixed_quad(
    func: _FixedQuadFunc[onp.ArrayND[np.complex64 | np.complex128, _JustAnyShape]],
    a: float,
    b: float,
    args: tuple[object, ...] = (),
    n: int = 5,
) -> tuple[np.complex128 | onp.ArrayND[np.complex128], None]: ...
@overload  # (1d f64) -> 1d ~c64 | ~c128
def fixed_quad(
    func: _FixedQuadFunc[onp.ToArrayStrict1D[op.JustComplex, np.complex64 | np.complex128]],
    a: float,
    b: float,
    args: tuple[object, ...] = (),
    n: int = 5,
) -> tuple[np.complex128, None]: ...
@overload  # (1d f64) -> 2d ~c64 | ~c128
def fixed_quad(
    func: _FixedQuadFunc[onp.ToArrayStrict2D[op.JustComplex, np.complex64 | np.complex128]],
    a: float,
    b: float,
    args: tuple[object, ...] = (),
    n: int = 5,
) -> tuple[onp.Array1D[np.complex128], None]: ...
@overload  # (1d f64) -> Nd ~c64 | ~c128
def fixed_quad(
    func: _FixedQuadFunc[onp.ToArrayND[op.JustComplex, np.complex64 | np.complex128]],
    a: float,
    b: float,
    args: tuple[object, ...] = (),
    n: int = 5,
) -> tuple[onp.ArrayND[np.complex128] | Any, None]: ...
@overload  # (1d f64) -> ?d ~f80 | ~c160 | timedelta64
def fixed_quad[ScalarT: npc.inexact80 | np.timedelta64](
    func: _FixedQuadFunc[onp.ArrayND[ScalarT, _JustAnyShape]], a: float, b: float, args: tuple[object, ...] = (), n: int = 5
) -> tuple[ScalarT | onp.ArrayND[ScalarT], None]: ...
@overload  # (1d f64) -> 1d ~f80 | ~c160 | timedelta64
def fixed_quad[ScalarT: npc.inexact80 | np.timedelta64](
    func: _FixedQuadFunc[onp.Array1D[ScalarT]], a: float, b: float, args: tuple[object, ...] = (), n: int = 5
) -> tuple[ScalarT, None]: ...
@overload  # (1d f64) -> 2d ~f80 | ~c160 | timedelta64
def fixed_quad[ScalarT: npc.inexact80 | np.timedelta64](
    func: _FixedQuadFunc[onp.Array2D[ScalarT]], a: float, b: float, args: tuple[object, ...] = (), n: int = 5
) -> tuple[onp.Array1D[ScalarT], None]: ...
@overload  # (1d f64) -> Nd ~f80 | ~c160 | timedelta64
def fixed_quad[ScalarT: npc.inexact80 | np.timedelta64](
    func: _FixedQuadFunc[onp.ArrayND[ScalarT]], a: float, b: float, args: tuple[object, ...] = (), n: int = 5
) -> tuple[onp.ArrayND[ScalarT] | Any, None]: ...
@overload  # (1d f64) -> Nd ~object_
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
