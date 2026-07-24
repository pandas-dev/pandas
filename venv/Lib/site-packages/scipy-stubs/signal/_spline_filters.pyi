from typing import Any, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

__all__ = [
    "cspline1d",
    "cspline1d_eval",
    "cspline2d",
    "gauss_spline",
    "qspline1d",
    "qspline1d_eval",
    "qspline2d",
    "spline_filter",
    "symiirorder1",
    "symiirorder2",
]

###

type _FloatQ = np.float64 | np.longdouble
type _ComplexQ = np.complex128 | np.clongdouble

###

#
def spline_filter[FloatDT: np.float32 | np.float64](
    Iin: onp.ArrayND[FloatDT], lmbda: onp.ToFloat = 5.0
) -> onp.Array2D[FloatDT]: ...

# NOTE: Mypy reports a false positive `overload-overlap` error with `numpy<2.1`.
# mypy: disable-error-code=overload-overlap

#
@overload
def gauss_spline[ShapeT: tuple[int, ...]](
    x: onp.ArrayND[npc.integer | np.bool, ShapeT], n: onp.ToFloat
) -> onp.ArrayND[np.float64, ShapeT]: ...
@overload
def gauss_spline[InexactQT: npc.inexact, ShapeT: tuple[int, ...]](
    x: onp.ArrayND[InexactQT, ShapeT], n: onp.ToFloat
) -> onp.ArrayND[InexactQT, ShapeT]: ...
@overload
def gauss_spline(x: onp.ToFloatStrict1D, n: onp.ToFloat) -> onp.Array1D[_FloatQ]: ...
@overload
def gauss_spline(x: onp.ToFloatStrict2D, n: onp.ToFloat) -> onp.Array2D[_FloatQ]: ...
@overload
def gauss_spline(x: onp.ToFloatStrict3D, n: onp.ToFloat) -> onp.Array3D[_FloatQ]: ...
@overload  # the weird shape-type is a workaround for a bug in pyright's overlapping overload detection
def gauss_spline(x: onp.ToFloatND, n: onp.ToFloat) -> onp.ArrayND[_FloatQ, tuple[int] | tuple[Any, ...]]: ...
@overload
def gauss_spline(x: onp.ToJustComplexStrict1D, n: onp.ToFloat) -> onp.Array1D[_ComplexQ]: ...
@overload
def gauss_spline(x: onp.ToJustComplexStrict2D, n: onp.ToFloat) -> onp.Array2D[_ComplexQ]: ...
@overload
def gauss_spline(x: onp.ToJustComplexStrict3D, n: onp.ToFloat) -> onp.Array3D[_ComplexQ]: ...
@overload
def gauss_spline(x: onp.ToJustComplexND, n: onp.ToFloat) -> onp.ArrayND[_ComplexQ]: ...

#
@overload
def cspline1d[InexactQT: npc.inexact64 | npc.inexact80](
    signal: onp.ArrayND[InexactQT], lamb: onp.ToFloat = 0.0
) -> onp.Array1D[InexactQT]: ...
@overload
def cspline1d(signal: onp.ToFloatND, lamb: onp.ToFloat = 0.0) -> onp.Array1D[_FloatQ]: ...
@overload
def cspline1d(signal: onp.ToJustComplexND, lamb: onp.ToFloat = 0.0) -> onp.Array1D[_ComplexQ]: ...

#
@overload
def qspline1d[InexactQT: npc.inexact64 | npc.inexact80](
    signal: onp.ArrayND[InexactQT], lamb: onp.ToFloat = 0.0
) -> onp.Array1D[InexactQT]: ...
@overload
def qspline1d(signal: onp.ToFloatND, lamb: onp.ToFloat = 0.0) -> onp.Array1D[_FloatQ]: ...
@overload
def qspline1d(signal: onp.ToJustComplexND, lamb: onp.ToFloat = 0.0) -> onp.Array1D[_ComplexQ]: ...

#
@overload
def cspline2d(signal: onp.ToFloatND, lamb: onp.ToFloat = 0.0, precision: onp.ToFloat = -1.0) -> onp.Array2D[np.float64]: ...
@overload
def cspline2d(
    signal: onp.ToJustComplexND, lamb: onp.ToFloat = 0.0, precision: onp.ToFloat = -1.0
) -> onp.Array2D[np.complex128]: ...

#
@overload
def qspline2d(signal: onp.ToFloatND, lamb: onp.ToFloat = 0.0, precision: onp.ToFloat = -1.0) -> onp.Array2D[np.float64]: ...
@overload
def qspline2d(
    signal: onp.ToJustComplexND, lamb: onp.ToFloat = 0.0, precision: onp.ToFloat = -1.0
) -> onp.Array2D[np.complex128]: ...

#
def cspline1d_eval[InexactT: npc.inexact](
    cj: onp.Array1D[InexactT], newx: onp.ToFloatND, dx: onp.ToFloat = 1.0, x0: onp.ToFloat = 0
) -> onp.Array1D[InexactT]: ...

#
def qspline1d_eval[InexactQT: npc.inexact64 | npc.inexact80](
    cj: onp.Array1D[InexactQT], newx: onp.ToFloatND, dx: onp.ToFloat = 1.0, x0: onp.ToFloat = 0
) -> onp.Array1D[InexactQT]: ...

#
def symiirorder1[InexactDT: npc.inexact32 | npc.inexact64, ShapeT: tuple[int, ...]](
    signal: onp.ArrayND[InexactDT, ShapeT], c0: onp.ToComplex, z1: onp.ToComplex, precision: onp.ToFloat = -1.0
) -> onp.ArrayND[InexactDT, ShapeT]: ...

#
def symiirorder2[FloatDT: np.float32 | np.float64, ShapeT: tuple[int, ...]](
    input: onp.ArrayND[FloatDT, ShapeT], r: onp.ToFloat, omega: onp.ToFloat, precision: onp.ToFloat = -1.0
) -> onp.ArrayND[FloatDT, ShapeT]: ...
