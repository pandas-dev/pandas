from collections.abc import Iterable
from typing import Any, Final, Literal as L, Protocol, SupportsIndex, TypeVar, overload, type_check_only
from typing_extensions import TypeIs

import numpy as np
import numpy.typing as npt
import optype as op
import optype.numpy as onp
import optype.numpy.compat as npc

from scipy.sparse import (
    bsr_array,
    bsr_matrix,
    coo_array,
    coo_matrix,
    csc_array,
    csc_matrix,
    csr_array,
    csr_matrix,
    dia_array,
    dia_matrix,
)

__all__ = [
    "broadcast_shapes",
    "get_sum_dtype",
    "getdata",
    "getdtype",
    "isdense",
    "isintlike",
    "ismatrix",
    "isscalarlike",
    "issequence",
    "isshape",
    "upcast",
]

###

type _Axis = L[-2, -1, 0, 1] | npc.integer
type _ShapeLike = Iterable[SupportsIndex]
type _ScalarLike = complex | bytes | str | np.generic | onp.Array0D
type _SequenceLike = tuple[_ScalarLike, ...] | list[_ScalarLike] | onp.Array1D
type _MatrixLike = tuple[_SequenceLike, ...] | list[_SequenceLike] | onp.Array2D

type _IntP = np.int32 | np.int64
type _UIntP = np.uint32 | np.uint64

_ShapeT = TypeVar("_ShapeT", bound=tuple[int, ...], default=Any)
_ScalarT = TypeVar("_ScalarT", bound=np.generic, default=Any)

@type_check_only
class _SizedIndexIterable(Protocol):
    def __len__(self, /) -> int: ...
    def __iter__(self, /) -> op.CanNext[SupportsIndex]: ...

###

supported_dtypes: Final[list[type[npc.number | np.bool]]] = ...

#
# NOTE: Technically any `numpy.generic` could be returned, but we only care about the supported scalar types in `scipy.sparse`.
def upcast(*args: npt.DTypeLike) -> npc.number | np.bool: ...
def upcast_char(*args: npt.DTypeLike) -> npc.number | np.bool: ...
@overload
def upcast_scalar(dtype: onp.ToDType[_ScalarT], scalar: onp.ToScalar) -> np.dtype[_ScalarT]: ...
@overload
def upcast_scalar(dtype: npt.DTypeLike, scalar: onp.ToScalar) -> np.dtype[Any]: ...

#
def downcast_intp_index(
    arr: onp.Array[_ShapeT, np.bool | npc.integer | npc.floating | np.timedelta64 | np.object_],
) -> onp.Array[_ShapeT, _IntP]: ...

#
@overload
def to_native(A: _ScalarT) -> onp.Array0D[_ScalarT]: ...
@overload
def to_native(A: onp.Array[_ShapeT, _ScalarT]) -> onp.Array[_ShapeT, _ScalarT]: ...
@overload
def to_native[DTypeT: np.dtype[Any]](A: onp.HasDType[DTypeT]) -> np.ndarray[Any, DTypeT]: ...

#
def getdtype(
    dtype: onp.ToDType[_ScalarT] | None,
    a: onp.HasDType[np.dtype[_ScalarT]] | None = None,
    default: onp.ToDType[_ScalarT] | None = None,
) -> np.dtype[_ScalarT]: ...

#
@overload
def getdata(obj: _ScalarT, dtype: onp.ToDType[_ScalarT] | None = None, copy: bool = False) -> onp.Array0D[_ScalarT]: ...
@overload
def getdata(obj: onp.ToComplex, dtype: onp.ToDType[_ScalarT], copy: bool = False) -> onp.Array0D[_ScalarT]: ...
@overload
def getdata(obj: onp.ToComplexStrict1D, dtype: onp.ToDType[_ScalarT], copy: bool = False) -> onp.Array1D[_ScalarT]: ...
@overload
def getdata(obj: onp.ToComplexStrict2D, dtype: onp.ToDType[_ScalarT], copy: bool = False) -> onp.Array2D[_ScalarT]: ...
@overload
def getdata(obj: onp.ToComplexStrict3D, dtype: onp.ToDType[_ScalarT], copy: bool = False) -> onp.Array3D[_ScalarT]: ...
@overload
def getdata(
    obj: onp.ToArrayND[_ScalarT, _ScalarT], dtype: onp.ToDType[_ScalarT] | None = None, copy: bool = False
) -> onp.ArrayND[_ScalarT]: ...

type _CoInt32 = np.bool | np.int8 | np.uint8 | np.int16 | np.uint16 | np.int32
type _ContraInt32 = np.uint32 | np.int64 | np.uint64

#
@overload
def get_index_dtype(
    arrays: tuple[()] = (), maxval: onp.ToFloat | None = None, check_contents: bool = False
) -> type[np.int32]: ...
@overload
def get_index_dtype(
    arrays: tuple[onp.CanArrayND[_CoInt32], *tuple[onp.CanArrayND[_CoInt32], ...]],
    maxval: onp.ToFloat | None = None,
    check_contents: bool = False,
) -> type[np.int32]: ...
@overload
def get_index_dtype(
    arrays: tuple[onp.CanArrayND[_ContraInt32], *tuple[onp.CanArrayND[_ContraInt32], ...]],
    maxval: onp.ToFloat | None = None,
    check_contents: bool = False,
) -> type[np.int64]: ...
@overload
def get_index_dtype(
    arrays: tuple[onp.ToInt | onp.ToIntND, ...], maxval: onp.ToFloat | None = None, check_contents: bool = False
) -> type[_IntP]: ...

# NOTE: The inline annotations (`(np.dtype) -> np.dtype`) are incorrect.
@overload
def get_sum_dtype(dtype: np.dtype[npc.unsignedinteger]) -> type[_UIntP]: ...
@overload
def get_sum_dtype(dtype: np.dtype[np.bool | npc.signedinteger]) -> type[_IntP]: ...
@overload
def get_sum_dtype[DTypeT: np.dtype[npc.inexact | np.flexible | np.datetime64 | np.timedelta64 | np.object_]](
    dtype: DTypeT,
) -> DTypeT: ...

#
# NOTE: all arrays implement `__index__` but if it raises this returns `False`, so `TypeIs` can't be used here
def isintlike(x: object) -> TypeIs[SupportsIndex]: ...
def isscalarlike(x: object) -> TypeIs[_ScalarLike]: ...
def isshape(x: _SizedIndexIterable, nonneg: bool = False, *, allow_nd: tuple[int, ...] = (2,), check_nd: bool = True) -> bool: ...
def issequence(t: object) -> TypeIs[_SequenceLike]: ...  # undocumented
def ismatrix(t: object) -> TypeIs[_MatrixLike]: ...  # undocumented
def isdense(x: object) -> TypeIs[onp.Array]: ...  # undocumented

#
@overload
def validateaxis(axis: None, *, ndim: int = 2) -> None: ...  # undocumented
@overload
def validateaxis(
    axis: _Axis | tuple[_Axis, *tuple[_Axis, ...]] | None, *, ndim: int = 2
) -> tuple[int, ...] | None: ...  # undocumented

#
def check_shape(
    args: _ShapeLike | tuple[_ShapeLike, ...], current_shape: tuple[int, ...] | None = None, *, allow_nd: tuple[int, ...] = (2,)
) -> tuple[int, ...]: ...

#
def matrix(
    object: onp.ToArray2D[_ScalarT],
    dtype: onp.ToDType[_ScalarT] | type | str | None = None,
    *,
    copy: L[0, 1, 2] | bool | None = True,
    order: L["K", "A", "C", "F"] = "K",
    subok: bool = False,
    ndmin: L[0, 1, 2] = 0,
    like: onp.CanArrayFunction | None = None,
) -> onp.Matrix[_ScalarT]: ...

#
@overload
def asmatrix(data: onp.ToArray2D[Any], dtype: onp.ToDType[_ScalarT]) -> onp.Matrix[_ScalarT]: ...
@overload
def asmatrix(data: onp.ToArray2D[_ScalarT], dtype: onp.ToDType[_ScalarT] | None = None) -> onp.Matrix[_ScalarT]: ...
@overload
def asmatrix(data: onp.ToArray2D[Any], dtype: npt.DTypeLike) -> onp.Matrix[Any]: ...

#
@overload  # BSR/CSC/CSR, dtype: <default>
def safely_cast_index_arrays(
    A: bsr_array | bsr_matrix | csc_array | csc_matrix | csr_array | csr_matrix,
    idx_dtype: onp.ToDType[np.int32] = ...,  # = np.int32
    msg: str = "",
) -> tuple[onp.Array1D[np.int32], onp.Array1D[np.int32]]: ...
@overload  # BSR/CSC/CSR, dtype: <known>
def safely_cast_index_arrays[IntT: npc.integer](
    A: bsr_array | bsr_matrix | csc_array | csc_matrix | csr_array | csr_matrix, idx_dtype: onp.ToDType[IntT], msg: str = ""
) -> tuple[onp.Array1D[IntT], onp.Array1D[IntT]]: ...
@overload  # 2d COO, dtype: <default>
def safely_cast_index_arrays(
    A: coo_array[Any, tuple[int, int]] | coo_matrix,
    idx_dtype: onp.ToDType[np.int32] = ...,  # = np.int32
    msg: str = "",
) -> tuple[onp.Array1D[np.int32], onp.Array1D[np.int32]]: ...
@overload  # 2d COO, dtype: <known>
def safely_cast_index_arrays[IntT: npc.integer](
    A: coo_array[Any, tuple[int, int]] | coo_matrix, idx_dtype: onp.ToDType[IntT], msg: str = ""
) -> tuple[onp.Array1D[IntT], onp.Array1D[IntT]]: ...
@overload  # nd COO, dtype: <default>
def safely_cast_index_arrays(
    A: coo_array,
    idx_dtype: onp.ToDType[np.int32] = ...,  # = np.int32
    msg: str = "",
) -> tuple[onp.Array1D[np.int32], ...]: ...
@overload  # nd COO, dtype: <known>
def safely_cast_index_arrays[IntT: npc.integer](
    A: coo_array, idx_dtype: onp.ToDType[IntT], msg: str = ""
) -> tuple[onp.Array1D[IntT], ...]: ...
@overload  # DIA, dtype: <default>
def safely_cast_index_arrays(
    A: dia_array | dia_matrix,
    idx_dtype: onp.ToDType[np.int32] = ...,  # = np.int32
    msg: str = "",
) -> onp.Array1D[np.int32]: ...
@overload  # DIA, dtype: <known>
def safely_cast_index_arrays[IntT: npc.integer](
    A: dia_array | dia_matrix, idx_dtype: onp.ToDType[IntT], msg: str = ""
) -> onp.Array1D[IntT]: ...

#
@overload
def broadcast_shapes() -> tuple[()]: ...
@overload
def broadcast_shapes(shape0: tuple[()], /, *shapes: tuple[()]) -> tuple[()]: ...
@overload
def broadcast_shapes(shape0: tuple[int], /, *shapes: onp.AtMost1D) -> tuple[int]: ...
@overload
def broadcast_shapes(shape0: tuple[int, int], /, *shapes: onp.AtMost2D) -> tuple[int, int]: ...
@overload
def broadcast_shapes(shape0: tuple[int, int, int], /, *shapes: onp.AtMost3D) -> tuple[int, int, int]: ...
@overload
def broadcast_shapes(shape0: _ShapeT, /, *shapes: tuple[()] | _ShapeT) -> _ShapeT: ...
