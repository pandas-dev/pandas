from typing import Any, Literal, SupportsIndex, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy._typing import AnyShape
from scipy.fft._typing import DCTType

__all__ = ["dct", "dctn", "dst", "dstn", "idct", "idctn", "idst", "idstn"]

###

type _Inexact = np.float32 | np.float64 | npc.floating80 | npc.complexfloating

type _NormKind = Literal["ortho"] | None

# workaround for a strange bug in pyright's overlapping overload detection with `numpy<2.1`
type _WorkaroundForPyright = tuple[int] | tuple[Any, ...]
type _FloatND = onp.ArrayND[np.float32 | np.float64 | np.longdouble, _WorkaroundForPyright]

###

@overload
def dctn[ShapeT: tuple[int, ...]](
    x: onp.CanArrayND[npc.integer, ShapeT],
    type: DCTType = 2,
    shape: AnyShape | None = None,
    axes: AnyShape | None = None,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> onp.Array[ShapeT, np.float64]: ...
@overload
def dctn[ShapeT: tuple[int, ...]](
    x: onp.CanArrayND[np.float16, ShapeT],
    type: DCTType = 2,
    shape: AnyShape | None = None,
    axes: AnyShape | None = None,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> onp.Array[ShapeT, np.float32]: ...
@overload
def dctn[ShapeT: tuple[int, ...], DTypeT: np.dtype[_Inexact]](
    x: onp.CanArray[ShapeT, DTypeT],
    type: DCTType = 2,
    shape: AnyShape | None = None,
    axes: AnyShape | None = None,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> np.ndarray[ShapeT, DTypeT]: ...
@overload
def dctn(
    x: onp.SequenceND[float],
    type: DCTType = 2,
    shape: AnyShape | None = None,
    axes: AnyShape | None = None,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> onp.ArrayND[np.float64]: ...
@overload
def dctn(
    x: onp.SequenceND[list[complex]] | list[complex],
    type: DCTType = 2,
    shape: AnyShape | None = None,
    axes: AnyShape | None = None,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> onp.ArrayND[np.complex128]: ...
@overload
def dctn(
    x: onp.ToFloatND,
    type: DCTType = 2,
    shape: AnyShape | None = None,
    axes: AnyShape | None = None,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> _FloatND: ...

#
@overload
def idctn[ShapeT: tuple[int, ...]](
    x: onp.CanArrayND[npc.integer, ShapeT],
    type: DCTType = 2,
    shape: AnyShape | None = None,
    axes: AnyShape | None = None,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> onp.Array[ShapeT, np.float64]: ...
@overload
def idctn[ShapeT: tuple[int, ...]](
    x: onp.CanArrayND[np.float16, ShapeT],
    type: DCTType = 2,
    shape: AnyShape | None = None,
    axes: AnyShape | None = None,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> onp.Array[ShapeT, np.float32]: ...
@overload
def idctn[ShapeT: tuple[int, ...], DTypeT: np.dtype[_Inexact]](
    x: onp.CanArray[ShapeT, DTypeT],
    type: DCTType = 2,
    shape: AnyShape | None = None,
    axes: AnyShape | None = None,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> np.ndarray[ShapeT, DTypeT]: ...
@overload
def idctn(
    x: onp.SequenceND[float],
    type: DCTType = 2,
    shape: AnyShape | None = None,
    axes: AnyShape | None = None,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> onp.ArrayND[np.float64]: ...
@overload
def idctn(
    x: onp.SequenceND[list[complex]] | list[complex],
    type: DCTType = 2,
    shape: AnyShape | None = None,
    axes: AnyShape | None = None,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> onp.ArrayND[np.complex128]: ...
@overload
def idctn(
    x: onp.ToFloatND,
    type: DCTType = 2,
    shape: AnyShape | None = None,
    axes: AnyShape | None = None,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> _FloatND: ...

#
@overload
def dstn[ShapeT: tuple[int, ...]](
    x: onp.CanArrayND[npc.integer, ShapeT],
    type: DCTType = 2,
    shape: AnyShape | None = None,
    axes: AnyShape | None = None,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> onp.Array[ShapeT, np.float64]: ...
@overload
def dstn[ShapeT: tuple[int, ...]](
    x: onp.CanArrayND[np.float16, ShapeT],
    type: DCTType = 2,
    shape: AnyShape | None = None,
    axes: AnyShape | None = None,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> onp.Array[ShapeT, np.float32]: ...
@overload
def dstn[ShapeT: tuple[int, ...], DTypeT: np.dtype[_Inexact]](
    x: onp.CanArray[ShapeT, DTypeT],
    type: DCTType = 2,
    shape: AnyShape | None = None,
    axes: AnyShape | None = None,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> np.ndarray[ShapeT, DTypeT]: ...
@overload
def dstn(
    x: onp.SequenceND[float],
    type: DCTType = 2,
    shape: AnyShape | None = None,
    axes: AnyShape | None = None,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> onp.ArrayND[np.float64]: ...
@overload
def dstn(
    x: onp.SequenceND[list[complex]] | list[complex],
    type: DCTType = 2,
    shape: AnyShape | None = None,
    axes: AnyShape | None = None,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> onp.ArrayND[np.complex128]: ...
@overload
def dstn(
    x: onp.ToFloatND,
    type: DCTType = 2,
    shape: AnyShape | None = None,
    axes: AnyShape | None = None,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> _FloatND: ...

#
@overload
def idstn[ShapeT: tuple[int, ...]](
    x: onp.CanArrayND[npc.integer, ShapeT],
    type: DCTType = 2,
    shape: AnyShape | None = None,
    axes: AnyShape | None = None,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> onp.Array[ShapeT, np.float64]: ...
@overload
def idstn[ShapeT: tuple[int, ...]](
    x: onp.CanArrayND[np.float16, ShapeT],
    type: DCTType = 2,
    shape: AnyShape | None = None,
    axes: AnyShape | None = None,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> onp.Array[ShapeT, np.float32]: ...
@overload
def idstn[ShapeT: tuple[int, ...], DTypeT: np.dtype[_Inexact]](
    x: onp.CanArray[ShapeT, DTypeT],
    type: DCTType = 2,
    shape: AnyShape | None = None,
    axes: AnyShape | None = None,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> np.ndarray[ShapeT, DTypeT]: ...
@overload
def idstn(
    x: onp.SequenceND[float],
    type: DCTType = 2,
    shape: AnyShape | None = None,
    axes: AnyShape | None = None,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> onp.ArrayND[np.float64]: ...
@overload
def idstn(
    x: onp.SequenceND[list[complex]] | list[complex],
    type: DCTType = 2,
    shape: AnyShape | None = None,
    axes: AnyShape | None = None,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> onp.ArrayND[np.complex128]: ...
@overload
def idstn(
    x: onp.ToFloatND,
    type: DCTType = 2,
    shape: AnyShape | None = None,
    axes: AnyShape | None = None,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> _FloatND: ...

#
@overload
def dct[ShapeT: tuple[int, ...]](
    x: onp.CanArrayND[npc.integer, ShapeT],
    type: DCTType = 2,
    n: onp.ToInt | None = None,
    axis: SupportsIndex = -1,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> onp.Array[ShapeT, np.float64]: ...
@overload
def dct[ShapeT: tuple[int, ...]](
    x: onp.CanArrayND[np.float16, ShapeT],
    type: DCTType = 2,
    n: onp.ToInt | None = None,
    axis: SupportsIndex = -1,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> onp.Array[ShapeT, np.float32]: ...
@overload
def dct[ShapeT: tuple[int, ...], DTypeT: np.dtype[_Inexact]](
    x: onp.CanArray[ShapeT, DTypeT],
    type: DCTType = 2,
    n: onp.ToInt | None = None,
    axis: SupportsIndex = -1,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> np.ndarray[ShapeT, DTypeT]: ...
@overload
def dct(
    x: onp.SequenceND[float],
    type: DCTType = 2,
    n: onp.ToInt | None = None,
    axis: SupportsIndex = -1,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> onp.ArrayND[np.float64]: ...
@overload
def dct(
    x: onp.SequenceND[list[complex]] | list[complex],
    type: DCTType = 2,
    n: onp.ToInt | None = None,
    axis: SupportsIndex = -1,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> onp.ArrayND[np.complex128]: ...
@overload
def dct(
    x: onp.ToFloatND,
    type: DCTType = 2,
    n: onp.ToInt | None = None,
    axis: SupportsIndex = -1,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> _FloatND: ...

#
@overload
def idct[ShapeT: tuple[int, ...]](
    x: onp.CanArrayND[npc.integer, ShapeT],
    type: DCTType = 2,
    n: onp.ToInt | None = None,
    axis: SupportsIndex = -1,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> onp.Array[ShapeT, np.float64]: ...
@overload
def idct[ShapeT: tuple[int, ...]](
    x: onp.CanArrayND[np.float16, ShapeT],
    type: DCTType = 2,
    n: onp.ToInt | None = None,
    axis: SupportsIndex = -1,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> onp.Array[ShapeT, np.float32]: ...
@overload
def idct[ShapeT: tuple[int, ...], DTypeT: np.dtype[_Inexact]](
    x: onp.CanArray[ShapeT, DTypeT],
    type: DCTType = 2,
    n: onp.ToInt | None = None,
    axis: SupportsIndex = -1,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> np.ndarray[ShapeT, DTypeT]: ...
@overload
def idct(
    x: onp.SequenceND[float],
    type: DCTType = 2,
    n: onp.ToInt | None = None,
    axis: SupportsIndex = -1,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> onp.ArrayND[np.float64]: ...
@overload
def idct(
    x: onp.SequenceND[list[complex]] | list[complex],
    type: DCTType = 2,
    n: onp.ToInt | None = None,
    axis: SupportsIndex = -1,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> onp.ArrayND[np.complex128]: ...
@overload
def idct(
    x: onp.ToFloatND,
    type: DCTType = 2,
    n: onp.ToInt | None = None,
    axis: SupportsIndex = -1,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> _FloatND: ...

#
@overload
def dst[ShapeT: tuple[int, ...]](
    x: onp.CanArrayND[npc.integer, ShapeT],
    type: DCTType = 2,
    n: onp.ToInt | None = None,
    axis: SupportsIndex = -1,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> onp.Array[ShapeT, np.float64]: ...
@overload
def dst[ShapeT: tuple[int, ...]](
    x: onp.CanArrayND[np.float16, ShapeT],
    type: DCTType = 2,
    n: onp.ToInt | None = None,
    axis: SupportsIndex = -1,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> onp.Array[ShapeT, np.float32]: ...
@overload
def dst[ShapeT: tuple[int, ...], DTypeT: np.dtype[_Inexact]](
    x: onp.CanArray[ShapeT, DTypeT],
    type: DCTType = 2,
    n: onp.ToInt | None = None,
    axis: SupportsIndex = -1,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> np.ndarray[ShapeT, DTypeT]: ...
@overload
def dst(
    x: onp.SequenceND[float],
    type: DCTType = 2,
    n: onp.ToInt | None = None,
    axis: SupportsIndex = -1,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> onp.ArrayND[np.float64]: ...
@overload
def dst(
    x: onp.SequenceND[list[complex]] | list[complex],
    type: DCTType = 2,
    n: onp.ToInt | None = None,
    axis: SupportsIndex = -1,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> onp.ArrayND[np.complex128]: ...
@overload
def dst(
    x: onp.ToFloatND,
    type: DCTType = 2,
    n: onp.ToInt | None = None,
    axis: SupportsIndex = -1,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> _FloatND: ...

#
@overload
def idst[ShapeT: tuple[int, ...]](
    x: onp.CanArrayND[npc.integer, ShapeT],
    type: DCTType = 2,
    n: onp.ToInt | None = None,
    axis: SupportsIndex = -1,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> onp.Array[ShapeT, np.float64]: ...
@overload
def idst[ShapeT: tuple[int, ...]](
    x: onp.CanArrayND[np.float16, ShapeT],
    type: DCTType = 2,
    n: onp.ToInt | None = None,
    axis: SupportsIndex = -1,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> onp.Array[ShapeT, np.float32]: ...
@overload
def idst[ShapeT: tuple[int, ...], DTypeT: np.dtype[_Inexact]](
    x: onp.CanArray[ShapeT, DTypeT],
    type: DCTType = 2,
    n: onp.ToInt | None = None,
    axis: SupportsIndex = -1,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> np.ndarray[ShapeT, DTypeT]: ...
@overload
def idst(
    x: onp.SequenceND[float],
    type: DCTType = 2,
    n: onp.ToInt | None = None,
    axis: SupportsIndex = -1,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> onp.ArrayND[np.float64]: ...
@overload
def idst(
    x: onp.SequenceND[list[complex]] | list[complex],
    type: DCTType = 2,
    n: onp.ToInt | None = None,
    axis: SupportsIndex = -1,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> onp.ArrayND[np.complex128]: ...
@overload
def idst(
    x: onp.ToFloatND,
    type: DCTType = 2,
    n: onp.ToInt | None = None,
    axis: SupportsIndex = -1,
    norm: _NormKind = None,
    overwrite_x: bool = False,
) -> _FloatND: ...
