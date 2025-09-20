from typing import Any, TypeAlias, TypeVar, overload

import numpy as np
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from ._realtransforms import dct, dctn, dst, dstn, idct, idst
from ._typing import DCTType, NormalizationMode
from scipy._typing import AnyShape

__all__ = ["dct", "dctn", "dst", "dstn", "idct", "idctn", "idst", "idstn"]

_ShapeT = TypeVar("_ShapeT", bound=tuple[int, ...])
_DTypeT = TypeVar("_DTypeT", bound=np.dtype[np.float32 | np.float64 | npc.floating80 | npc.complexfloating])

# workaround for a strange bug in pyright's overlapping overload detection with `numpy<2.1`
_WorkaroundForPyright: TypeAlias = tuple[int] | tuple[Any, ...]

# NOTE: Unlike the ones in `scipy.fft._realtransforms`, `orthogonalize` is keyword-only here.

#
@overload
def idctn(
    x: onp.CanArrayND[npc.integer, _ShapeT],
    type: DCTType = 2,
    s: onp.ToInt | onp.ToIntND | None = None,
    axes: AnyShape | None = None,
    norm: NormalizationMode | None = None,
    overwrite_x: op.CanBool = False,
    workers: onp.ToInt | None = None,
    *,
    orthogonalize: op.CanBool | None = None,
) -> onp.Array[_ShapeT, np.float64]: ...
@overload
def idctn(
    x: onp.CanArrayND[np.float16, _ShapeT],
    type: DCTType = 2,
    s: onp.ToInt | onp.ToIntND | None = None,
    axes: AnyShape | None = None,
    norm: NormalizationMode | None = None,
    overwrite_x: op.CanBool = False,
    workers: onp.ToInt | None = None,
    *,
    orthogonalize: op.CanBool | None = None,
) -> onp.Array[_ShapeT, np.float32]: ...
@overload
def idctn(
    x: onp.CanArray[_ShapeT, _DTypeT],
    type: DCTType = 2,
    s: onp.ToInt | onp.ToIntND | None = None,
    axes: AnyShape | None = None,
    norm: NormalizationMode | None = None,
    overwrite_x: op.CanBool = False,
    workers: onp.ToInt | None = None,
    *,
    orthogonalize: op.CanBool | None = None,
) -> np.ndarray[_ShapeT, _DTypeT]: ...
@overload
def idctn(
    x: onp.ToJustFloat64_ND,
    type: DCTType = 2,
    s: onp.ToInt | onp.ToIntND | None = None,
    axes: AnyShape | None = None,
    norm: NormalizationMode | None = None,
    overwrite_x: op.CanBool = False,
    workers: onp.ToInt | None = None,
    *,
    orthogonalize: op.CanBool | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload
def idctn(
    x: onp.ToFloatND,
    type: DCTType = 2,
    s: onp.ToInt | onp.ToIntND | None = None,
    axes: AnyShape | None = None,
    norm: NormalizationMode | None = None,
    overwrite_x: op.CanBool = False,
    workers: onp.ToInt | None = None,
    *,
    orthogonalize: op.CanBool | None = None,
) -> onp.Array[_WorkaroundForPyright, npc.floating]: ...

#
@overload
def idstn(
    x: onp.CanArrayND[npc.integer, _ShapeT],
    type: DCTType = 2,
    s: onp.ToInt | onp.ToIntND | None = None,
    axes: AnyShape | None = None,
    norm: NormalizationMode | None = None,
    overwrite_x: op.CanBool = False,
    workers: onp.ToInt | None = None,
    *,
    orthogonalize: op.CanBool | None = None,
) -> onp.Array[_ShapeT, np.float64]: ...
@overload
def idstn(
    x: onp.CanArrayND[np.float16, _ShapeT],
    type: DCTType = 2,
    s: onp.ToInt | onp.ToIntND | None = None,
    axes: AnyShape | None = None,
    norm: NormalizationMode | None = None,
    overwrite_x: op.CanBool = False,
    workers: onp.ToInt | None = None,
    *,
    orthogonalize: op.CanBool | None = None,
) -> onp.Array[_ShapeT, np.float32]: ...
@overload
def idstn(
    x: onp.CanArray[_ShapeT, _DTypeT],
    type: DCTType = 2,
    s: onp.ToInt | onp.ToIntND | None = None,
    axes: AnyShape | None = None,
    norm: NormalizationMode | None = None,
    overwrite_x: op.CanBool = False,
    workers: onp.ToInt | None = None,
    *,
    orthogonalize: op.CanBool | None = None,
) -> np.ndarray[_ShapeT, _DTypeT]: ...
@overload
def idstn(
    x: onp.ToJustFloat64_ND,
    type: DCTType = 2,
    s: onp.ToInt | onp.ToIntND | None = None,
    axes: AnyShape | None = None,
    norm: NormalizationMode | None = None,
    overwrite_x: op.CanBool = False,
    workers: onp.ToInt | None = None,
    *,
    orthogonalize: op.CanBool | None = None,
) -> onp.ArrayND[np.float64]: ...
@overload
def idstn(
    x: onp.ToFloatND,
    type: DCTType = 2,
    s: onp.ToInt | onp.ToIntND | None = None,
    axes: AnyShape | None = None,
    norm: NormalizationMode | None = None,
    overwrite_x: op.CanBool = False,
    workers: onp.ToInt | None = None,
    *,
    orthogonalize: op.CanBool | None = None,
) -> onp.Array[_WorkaroundForPyright, npc.floating]: ...
