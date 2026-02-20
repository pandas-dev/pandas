from typing import Any, TypeAlias, overload
from typing_extensions import TypeVar

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

_SubFloat64: TypeAlias = np.bool_ | npc.integer | np.float16 | np.float32

_FloatQ: TypeAlias = np.float64 | np.longdouble
_ComplexQ: TypeAlias = np.complex128 | np.clongdouble

_FloatDT = TypeVar("_FloatDT", bound=np.float32 | np.float64)
_InexactDT = TypeVar("_InexactDT", bound=npc.inexact32 | npc.inexact64)
_InexactQT = TypeVar("_InexactQT", bound=npc.inexact64 | npc.inexact80)
_InexactT = TypeVar("_InexactT", bound=npc.inexact)
_ShapeT = TypeVar("_ShapeT", bound=tuple[int, ...])

###

#
def spline_filter(Iin: onp.ArrayND[_FloatDT], lmbda: onp.ToFloat = 5.0) -> onp.Array2D[_FloatDT]: ...

# NOTE: Mypy reports a false positive `overload-overlap` error with `numpy<2.1`.
# mypy: disable-error-code=overload-overlap

#
@overload
def gauss_spline(x: onp.ArrayND[_SubFloat64, _ShapeT], n: onp.ToFloat) -> onp.ArrayND[np.float64, _ShapeT]: ...
@overload
def gauss_spline(x: onp.ArrayND[_InexactQT, _ShapeT], n: onp.ToFloat) -> onp.ArrayND[_InexactQT, _ShapeT]: ...
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
def cspline1d(signal: onp.ArrayND[_InexactQT], lamb: onp.ToFloat = 0.0) -> onp.Array1D[_InexactQT]: ...
@overload
def cspline1d(signal: onp.ToFloatND, lamb: onp.ToFloat = 0.0) -> onp.Array1D[_FloatQ]: ...
@overload
def cspline1d(signal: onp.ToJustComplexND, lamb: onp.ToFloat = 0.0) -> onp.Array1D[_ComplexQ]: ...

#
@overload
def qspline1d(signal: onp.ArrayND[_InexactQT], lamb: onp.ToFloat = 0.0) -> onp.Array1D[_InexactQT]: ...
@overload
def qspline1d(signal: onp.ToFloatND, lamb: onp.ToFloat = 0.0) -> onp.Array1D[_FloatQ]: ...
@overload
def qspline1d(signal: onp.ToJustComplexND, lamb: onp.ToFloat = 0.0) -> onp.Array1D[_ComplexQ]: ...

#
@overload
def cspline2d(signal: onp.ToFloatND, lamb: onp.ToFloat = 0.0, precision: onp.ToFloat = -1.0) -> onp.Array1D[np.float64]: ...
@overload
def cspline2d(
    signal: onp.ToJustComplexND, lamb: onp.ToFloat = 0.0, precision: onp.ToFloat = -1.0
) -> onp.Array1D[np.complex128]: ...

#
@overload
def qspline2d(signal: onp.ToFloatND, lamb: onp.ToFloat = 0.0, precision: onp.ToFloat = -1.0) -> onp.Array1D[np.float64]: ...
@overload
def qspline2d(
    signal: onp.ToJustComplexND, lamb: onp.ToFloat = 0.0, precision: onp.ToFloat = -1.0
) -> onp.Array1D[np.complex128]: ...

#
def cspline1d_eval(
    cj: onp.Array1D[_InexactT], newx: onp.ToFloatND, dx: onp.ToFloat = 1.0, x0: onp.ToFloat = 0
) -> onp.Array1D[_InexactT]: ...

#
def qspline1d_eval(
    cj: onp.Array1D[_InexactQT], newx: onp.ToFloatND, dx: onp.ToFloat = 1.0, x0: onp.ToFloat = 0
) -> onp.Array1D[_InexactQT]: ...

#
def symiirorder1(
    signal: onp.ArrayND[_InexactDT, _ShapeT], c0: onp.ToComplex, z1: onp.ToComplex, precision: onp.ToFloat = -1.0
) -> onp.ArrayND[_InexactDT, _ShapeT]: ...

#
def symiirorder2(
    input: onp.ArrayND[_FloatDT, _ShapeT], r: onp.ToFloat, omega: onp.ToFloat, precision: onp.ToFloat = -1.0
) -> onp.ArrayND[_FloatDT, _ShapeT]: ...
