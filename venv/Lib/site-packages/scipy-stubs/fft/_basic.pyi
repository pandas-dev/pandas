from _typeshed import Unused
from collections.abc import Sequence
from typing import Any, Literal, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy._typing import AnyShape

###

type _Norm = Literal["backward", "ortho", "forward"]

type _CoInteger = npc.integer | np.bool

type _AsFloat32[ShapeT: tuple[int, ...]] = onp.CanArray[ShapeT, np.dtype[npc.floating32]]
type _AsFloat64[ShapeT: tuple[int, ...]] = onp.CanArray[ShapeT, np.dtype[npc.floating64 | _CoInteger]]
type _AsFloat80[ShapeT: tuple[int, ...]] = onp.CanArray[ShapeT, np.dtype[npc.floating80]]
type _AsComplex64[ShapeT: tuple[int, ...]] = onp.CanArray[ShapeT, np.dtype[npc.inexact32]]
type _AsComplex128[ShapeT: tuple[int, ...]] = onp.CanArray[ShapeT, np.dtype[npc.inexact64 | _CoInteger]]
type _AsComplex160[ShapeT: tuple[int, ...]] = onp.CanArray[ShapeT, np.dtype[npc.inexact80]]

type _ToFloat64_ND = onp.ToArrayND[float, npc.floating64 | _CoInteger]
type _ToComplex128_ND = onp.ToArrayND[complex, npc.inexact64 | _CoInteger]

# NOTE: The order of overloads has been carefully chosen to avoid triggering a pyright bug.

###
# 1-D

# keep in sync with `ifft`
@overload
def fft[ShapeT: tuple[int, ...]](
    x: _AsComplex128[ShapeT],
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128, ShapeT]: ...
@overload
def fft[ShapeT: tuple[int, ...]](
    x: _AsComplex64[ShapeT],
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex64, ShapeT]: ...
@overload
def fft[ShapeT: tuple[int, ...]](
    x: _AsComplex160[ShapeT],
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.clongdouble, ShapeT]: ...
@overload
def fft(
    x: Sequence[complex],
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.Array1D[np.complex128]: ...
@overload
def fft(
    x: _ToComplex128_ND,
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload
def fft(
    x: onp.ToComplexND,
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128 | Any]: ...

# keep in sync with `fft`
@overload
def ifft[ShapeT: tuple[int, ...]](
    x: _AsComplex128[ShapeT],
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128, ShapeT]: ...
@overload
def ifft[ShapeT: tuple[int, ...]](
    x: _AsComplex64[ShapeT],
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex64, ShapeT]: ...
@overload
def ifft[ShapeT: tuple[int, ...]](
    x: _AsComplex160[ShapeT],
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.clongdouble, ShapeT]: ...
@overload
def ifft(
    x: Sequence[complex],
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.Array1D[np.complex128]: ...
@overload
def ifft(
    x: _ToComplex128_ND,
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload
def ifft(
    x: onp.ToComplexND,
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128 | Any]: ...

# keep in sync with `ihfft`
@overload
def rfft[ShapeT: tuple[int, ...]](
    x: _AsFloat64[ShapeT],
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128, ShapeT]: ...
@overload
def rfft[ShapeT: tuple[int, ...]](
    x: _AsFloat32[ShapeT],
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex64, ShapeT]: ...
@overload
def rfft[ShapeT: tuple[int, ...]](
    x: _AsFloat80[ShapeT],
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.clongdouble, ShapeT]: ...
@overload
def rfft(
    x: Sequence[float],
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.Array1D[np.complex128]: ...
@overload
def rfft(
    x: _ToFloat64_ND,
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload
def rfft(
    x: onp.ToFloatND,
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128 | Any]: ...

# keep in sync with `hfft`
@overload
def irfft[ShapeT: tuple[int, ...]](
    x: _AsComplex128[ShapeT],
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.float64, ShapeT]: ...
@overload
def irfft[ShapeT: tuple[int, ...]](
    x: _AsComplex64[ShapeT],
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.float32, ShapeT]: ...
@overload
def irfft[ShapeT: tuple[int, ...]](
    x: _AsComplex160[ShapeT],
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.longdouble, ShapeT]: ...
@overload
def irfft(
    x: Sequence[complex],
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.Array1D[np.float64]: ...
@overload
def irfft(
    x: _ToComplex128_ND,
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload
def irfft(
    x: onp.ToComplexND,
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.float64 | Any]: ...

# keep in sync with `irfft`
@overload
def hfft[ShapeT: tuple[int, ...]](
    x: _AsComplex128[ShapeT],
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.float64, ShapeT]: ...
@overload
def hfft[ShapeT: tuple[int, ...]](
    x: _AsComplex64[ShapeT],
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.float32, ShapeT]: ...
@overload
def hfft[ShapeT: tuple[int, ...]](
    x: _AsComplex160[ShapeT],
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.longdouble, ShapeT]: ...
@overload
def hfft(
    x: Sequence[complex],
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.Array1D[np.float64]: ...
@overload
def hfft(
    x: _ToComplex128_ND,
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload
def hfft(
    x: onp.ToComplexND,
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.float64 | Any]: ...

#
# keep in sync with `rfft`
@overload
def ihfft[ShapeT: tuple[int, ...]](
    x: _AsFloat64[ShapeT],
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128, ShapeT]: ...
@overload
def ihfft[ShapeT: tuple[int, ...]](
    x: _AsFloat32[ShapeT],
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex64, ShapeT]: ...
@overload
def ihfft[ShapeT: tuple[int, ...]](
    x: _AsFloat80[ShapeT],
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.clongdouble, ShapeT]: ...
@overload
def ihfft(
    x: Sequence[float],
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.Array1D[np.complex128]: ...
@overload
def ihfft(
    x: _ToFloat64_ND,
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload
def ihfft(
    x: onp.ToFloatND,
    n: int | None = None,
    axis: int = -1,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128 | Any]: ...

###
# 2-D

# keep in sync with `ifft2`
@overload
def fft2[ShapeT: tuple[int, ...]](
    x: _AsComplex128[ShapeT],
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128, ShapeT]: ...
@overload
def fft2[ShapeT: tuple[int, ...]](
    x: _AsComplex64[ShapeT],
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex64, ShapeT]: ...
@overload
def fft2[ShapeT: tuple[int, ...]](
    x: _AsComplex160[ShapeT],
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.clongdouble, ShapeT]: ...
@overload
def fft2(
    x: Sequence[Sequence[complex]],
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.Array2D[np.complex128]: ...
@overload
def fft2(
    x: _ToComplex128_ND,
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload
def fft2(
    x: onp.ToComplexND,
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128 | Any]: ...

# keep in sync with `fft2`
@overload
def ifft2[ShapeT: tuple[int, ...]](
    x: _AsComplex128[ShapeT],
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128, ShapeT]: ...
@overload
def ifft2[ShapeT: tuple[int, ...]](
    x: _AsComplex64[ShapeT],
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex64, ShapeT]: ...
@overload
def ifft2[ShapeT: tuple[int, ...]](
    x: _AsComplex160[ShapeT],
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.clongdouble, ShapeT]: ...
@overload
def ifft2(
    x: Sequence[Sequence[complex]],
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.Array2D[np.complex128]: ...
@overload
def ifft2(
    x: _ToComplex128_ND,
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload
def ifft2(
    x: onp.ToComplexND,
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128 | Any]: ...

# keep in sync with `ihfft2`
@overload
def rfft2[ShapeT: tuple[int, ...]](
    x: _AsFloat64[ShapeT],
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128, ShapeT]: ...
@overload
def rfft2[ShapeT: tuple[int, ...]](
    x: _AsFloat32[ShapeT],
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex64, ShapeT]: ...
@overload
def rfft2[ShapeT: tuple[int, ...]](
    x: _AsFloat80[ShapeT],
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.clongdouble, ShapeT]: ...
@overload
def rfft2(
    x: Sequence[Sequence[float]],
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.Array2D[np.complex128]: ...
@overload
def rfft2(
    x: _ToFloat64_ND,
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload
def rfft2(
    x: onp.ToFloatND,
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128 | Any]: ...

# keep in sync with `hfft2`
@overload
def irfft2[ShapeT: tuple[int, ...]](
    x: _AsComplex128[ShapeT],
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.float64, ShapeT]: ...
@overload
def irfft2[ShapeT: tuple[int, ...]](
    x: _AsComplex64[ShapeT],
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.float32, ShapeT]: ...
@overload
def irfft2[ShapeT: tuple[int, ...]](
    x: _AsComplex160[ShapeT],
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.longdouble, ShapeT]: ...
@overload
def irfft2(
    x: Sequence[Sequence[complex]],
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.Array2D[np.float64]: ...
@overload
def irfft2(
    x: _ToComplex128_ND,
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload
def irfft2(
    x: onp.ToComplexND,
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.float64 | Any]: ...

# keep in sync with `irfft2`
@overload
def hfft2[ShapeT: tuple[int, ...]](
    x: _AsComplex128[ShapeT],
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.float64, ShapeT]: ...
@overload
def hfft2[ShapeT: tuple[int, ...]](
    x: _AsComplex64[ShapeT],
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.float32, ShapeT]: ...
@overload
def hfft2[ShapeT: tuple[int, ...]](
    x: _AsComplex160[ShapeT],
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.longdouble, ShapeT]: ...
@overload
def hfft2(
    x: Sequence[Sequence[complex]],
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.Array2D[np.float64]: ...
@overload
def hfft2(
    x: _ToComplex128_ND,
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload
def hfft2(
    x: onp.ToComplexND,
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.float64 | Any]: ...

# keep in sync with `rfft2`
@overload
def ihfft2[ShapeT: tuple[int, ...]](
    x: _AsFloat64[ShapeT],
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128, ShapeT]: ...
@overload
def ihfft2[ShapeT: tuple[int, ...]](
    x: _AsFloat32[ShapeT],
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex64, ShapeT]: ...
@overload
def ihfft2[ShapeT: tuple[int, ...]](
    x: _AsFloat80[ShapeT],
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.clongdouble, ShapeT]: ...
@overload
def ihfft2(
    x: Sequence[Sequence[float]],
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.Array2D[np.complex128]: ...
@overload
def ihfft2(
    x: _ToFloat64_ND,
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload
def ihfft2(
    x: onp.ToFloatND,
    s: onp.ToJustInt1D | None = None,
    axes: AnyShape = (-2, -1),
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128 | Any]: ...

###
# N-D

# keep in sync with `ifftn`
@overload
def fftn[ShapeT: tuple[int, ...]](
    x: _AsComplex128[ShapeT],
    s: onp.ToJustIntND | None = None,
    axes: AnyShape | None = None,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128, ShapeT]: ...
@overload
def fftn[ShapeT: tuple[int, ...]](
    x: _AsComplex64[ShapeT],
    s: onp.ToJustIntND | None = None,
    axes: AnyShape | None = None,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex64, ShapeT]: ...
@overload
def fftn[ShapeT: tuple[int, ...]](
    x: _AsComplex160[ShapeT],
    s: onp.ToJustIntND | None = None,
    axes: AnyShape | None = None,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.clongdouble, ShapeT]: ...
@overload
def fftn(
    x: _ToComplex128_ND,
    s: onp.ToJustIntND | None = None,
    axes: AnyShape | None = None,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload
def fftn(
    x: onp.ToComplexND,
    s: onp.ToJustIntND | None = None,
    axes: AnyShape | None = None,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128 | Any]: ...

# keep in sync with `fftn`
@overload
def ifftn[ShapeT: tuple[int, ...]](
    x: _AsComplex128[ShapeT],
    s: onp.ToJustIntND | None = None,
    axes: AnyShape | None = None,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128, ShapeT]: ...
@overload
def ifftn[ShapeT: tuple[int, ...]](
    x: _AsComplex64[ShapeT],
    s: onp.ToJustIntND | None = None,
    axes: AnyShape | None = None,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex64, ShapeT]: ...
@overload
def ifftn[ShapeT: tuple[int, ...]](
    x: _AsComplex160[ShapeT],
    s: onp.ToJustIntND | None = None,
    axes: AnyShape | None = None,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.clongdouble, ShapeT]: ...
@overload
def ifftn(
    x: _ToComplex128_ND,
    s: onp.ToJustIntND | None = None,
    axes: AnyShape | None = None,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload
def ifftn(
    x: onp.ToComplexND,
    s: onp.ToJustIntND | None = None,
    axes: AnyShape | None = None,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128 | Any]: ...

# keep in sync with `ihfftn`
@overload
def rfftn[ShapeT: tuple[int, ...]](
    x: _AsFloat64[ShapeT],
    s: onp.ToJustIntND | None = None,
    axes: AnyShape | None = None,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128, ShapeT]: ...
@overload
def rfftn[ShapeT: tuple[int, ...]](
    x: _AsFloat32[ShapeT],
    s: onp.ToJustIntND | None = None,
    axes: AnyShape | None = None,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex64, ShapeT]: ...
@overload
def rfftn[ShapeT: tuple[int, ...]](
    x: _AsFloat80[ShapeT],
    s: onp.ToJustIntND | None = None,
    axes: AnyShape | None = None,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.clongdouble, ShapeT]: ...
@overload
def rfftn(
    x: _ToFloat64_ND,
    s: onp.ToJustIntND | None = None,
    axes: AnyShape | None = None,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload
def rfftn(
    x: onp.ToFloatND,
    s: onp.ToJustIntND | None = None,
    axes: AnyShape | None = None,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128 | Any]: ...

# keep in sync with `hfftn`
@overload
def irfftn[ShapeT: tuple[int, ...]](
    x: _AsComplex128[ShapeT],
    s: onp.ToJustIntND | None = None,
    axes: AnyShape | None = None,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.float64, ShapeT]: ...
@overload
def irfftn[ShapeT: tuple[int, ...]](
    x: _AsComplex64[ShapeT],
    s: onp.ToJustIntND | None = None,
    axes: AnyShape | None = None,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.float32, ShapeT]: ...
@overload
def irfftn[ShapeT: tuple[int, ...]](
    x: _AsComplex160[ShapeT],
    s: onp.ToJustIntND | None = None,
    axes: AnyShape | None = None,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.longdouble, ShapeT]: ...
@overload
def irfftn(
    x: _ToComplex128_ND,
    s: onp.ToJustIntND | None = None,
    axes: AnyShape | None = None,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload
def irfftn(
    x: onp.ToComplexND,
    s: onp.ToJustIntND | None = None,
    axes: AnyShape | None = None,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.float64 | Any]: ...

# keep in sync with `irfftn`
@overload
def hfftn[ShapeT: tuple[int, ...]](
    x: _AsComplex128[ShapeT],
    s: onp.ToJustIntND | None = None,
    axes: AnyShape | None = None,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.float64, ShapeT]: ...
@overload
def hfftn[ShapeT: tuple[int, ...]](
    x: _AsComplex64[ShapeT],
    s: onp.ToJustIntND | None = None,
    axes: AnyShape | None = None,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.float32, ShapeT]: ...
@overload
def hfftn[ShapeT: tuple[int, ...]](
    x: _AsComplex160[ShapeT],
    s: onp.ToJustIntND | None = None,
    axes: AnyShape | None = None,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.longdouble, ShapeT]: ...
@overload
def hfftn(
    x: _ToComplex128_ND,
    s: onp.ToJustIntND | None = None,
    axes: AnyShape | None = None,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload
def hfftn(
    x: onp.ToComplexND,
    s: onp.ToJustIntND | None = None,
    axes: AnyShape | None = None,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.float64 | Any]: ...

# keep in sync with `rfftn`
@overload
def ihfftn[ShapeT: tuple[int, ...]](
    x: _AsFloat64[ShapeT],
    s: onp.ToJustIntND | None = None,
    axes: AnyShape | None = None,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128, ShapeT]: ...
@overload
def ihfftn[ShapeT: tuple[int, ...]](
    x: _AsFloat32[ShapeT],
    s: onp.ToJustIntND | None = None,
    axes: AnyShape | None = None,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex64, ShapeT]: ...
@overload
def ihfftn[ShapeT: tuple[int, ...]](
    x: _AsFloat80[ShapeT],
    s: onp.ToJustIntND | None = None,
    axes: AnyShape | None = None,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.clongdouble, ShapeT]: ...
@overload
def ihfftn(
    x: _ToFloat64_ND,
    s: onp.ToJustIntND | None = None,
    axes: AnyShape | None = None,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128]: ...
@overload
def ihfftn(
    x: onp.ToFloatND,
    s: onp.ToJustIntND | None = None,
    axes: AnyShape | None = None,
    norm: _Norm | None = None,
    overwrite_x: bool = False,
    workers: int | None = None,
    *,
    plan: Unused | None = None,
) -> onp.ArrayND[np.complex128 | Any]: ...
