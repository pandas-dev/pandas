from typing import Any, overload

import numpy as np
import optype.numpy as onp
import optype.numpy.compat as npc

from ._realtransforms import dct, dctn, dst, dstn, idct, idst
from ._typing import DCTType, NormalizationMode
from scipy._typing import AnyShape

__all__ = ["dct", "dctn", "dst", "dstn", "idct", "idctn", "idst", "idstn"]

###

type _Inexact = np.float32 | np.float64 | npc.floating80 | npc.complexfloating

# workaround for a strange bug in pyright's overlapping overload detection with `numpy<2.1`
type _WorkaroundForPyright = tuple[int] | tuple[Any, ...]

###

# NOTE: Unlike the ones in `scipy.fft._realtransforms`, `orthogonalize` is keyword-only here.

#
@overload
def idctn[ShapeT: tuple[int, ...]](
    x: onp.CanArrayND[npc.integer, ShapeT],
    type: DCTType = 2,
    s: onp.ToInt | onp.ToIntND | None = None,
    axes: AnyShape | None = None,
    norm: NormalizationMode | None = None,
    overwrite_x: bool = False,
    workers: onp.ToInt | None = None,
    *,
    orthogonalize: bool | None = None,
) -> onp.Array[ShapeT, np.float64]: ...
@overload
def idctn[ShapeT: tuple[int, ...]](
    x: onp.CanArrayND[np.float16, ShapeT],
    type: DCTType = 2,
    s: onp.ToInt | onp.ToIntND | None = None,
    axes: AnyShape | None = None,
    norm: NormalizationMode | None = None,
    overwrite_x: bool = False,
    workers: onp.ToInt | None = None,
    *,
    orthogonalize: bool | None = None,
) -> onp.Array[ShapeT, np.float32]: ...
@overload
def idctn[ShapeT: tuple[int, ...], DTypeT: np.dtype[_Inexact]](
    x: onp.CanArray[ShapeT, DTypeT],
    type: DCTType = 2,
    s: onp.ToInt | onp.ToIntND | None = None,
    axes: AnyShape | None = None,
    norm: NormalizationMode | None = None,
    overwrite_x: bool = False,
    workers: onp.ToInt | None = None,
    *,
    orthogonalize: bool | None = None,
) -> np.ndarray[ShapeT, DTypeT]: ...
@overload
def idctn(
    x: onp.ToJustFloat64_ND,
    type: DCTType = 2,
    s: onp.ToInt | onp.ToIntND | None = None,
    axes: AnyShape | None = None,
    norm: NormalizationMode | None = None,
    overwrite_x: bool = False,
    workers: onp.ToInt | None = None,
    *,
    orthogonalize: bool | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload
def idctn(
    x: onp.ToFloatND,
    type: DCTType = 2,
    s: onp.ToInt | onp.ToIntND | None = None,
    axes: AnyShape | None = None,
    norm: NormalizationMode | None = None,
    overwrite_x: bool = False,
    workers: onp.ToInt | None = None,
    *,
    orthogonalize: bool | None = None,
) -> onp.Array[_WorkaroundForPyright, npc.floating]: ...

#
@overload
def idstn[ShapeT: tuple[int, ...]](
    x: onp.CanArrayND[npc.integer, ShapeT],
    type: DCTType = 2,
    s: onp.ToInt | onp.ToIntND | None = None,
    axes: AnyShape | None = None,
    norm: NormalizationMode | None = None,
    overwrite_x: bool = False,
    workers: onp.ToInt | None = None,
    *,
    orthogonalize: bool | None = None,
) -> onp.Array[ShapeT, np.float64]: ...
@overload
def idstn[ShapeT: tuple[int, ...]](
    x: onp.CanArrayND[np.float16, ShapeT],
    type: DCTType = 2,
    s: onp.ToInt | onp.ToIntND | None = None,
    axes: AnyShape | None = None,
    norm: NormalizationMode | None = None,
    overwrite_x: bool = False,
    workers: onp.ToInt | None = None,
    *,
    orthogonalize: bool | None = None,
) -> onp.Array[ShapeT, np.float32]: ...
@overload
def idstn[ShapeT: tuple[int, ...], DTypeT: np.dtype[_Inexact]](
    x: onp.CanArray[ShapeT, DTypeT],
    type: DCTType = 2,
    s: onp.ToInt | onp.ToIntND | None = None,
    axes: AnyShape | None = None,
    norm: NormalizationMode | None = None,
    overwrite_x: bool = False,
    workers: onp.ToInt | None = None,
    *,
    orthogonalize: bool | None = None,
) -> np.ndarray[ShapeT, DTypeT]: ...
@overload
def idstn(
    x: onp.ToJustFloat64_ND,
    type: DCTType = 2,
    s: onp.ToInt | onp.ToIntND | None = None,
    axes: AnyShape | None = None,
    norm: NormalizationMode | None = None,
    overwrite_x: bool = False,
    workers: onp.ToInt | None = None,
    *,
    orthogonalize: bool | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload
def idstn(
    x: onp.ToFloatND,
    type: DCTType = 2,
    s: onp.ToInt | onp.ToIntND | None = None,
    axes: AnyShape | None = None,
    norm: NormalizationMode | None = None,
    overwrite_x: bool = False,
    workers: onp.ToInt | None = None,
    *,
    orthogonalize: bool | None = None,
) -> onp.Array[_WorkaroundForPyright, npc.floating]: ...
