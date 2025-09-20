from collections.abc import Iterable
from typing import Any, Final, Literal as L, Protocol, TypeAlias, TypeVar, TypedDict, overload, type_check_only
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

_ShapeT = TypeVar("_ShapeT", bound=tuple[int, ...], default=Any)
_DTypeT = TypeVar("_DTypeT", bound=np.dtype[Any])
_ScalarT = TypeVar("_ScalarT", bound=np.generic, default=Any)
_IntT = TypeVar("_IntT", bound=npc.integer)
_NonIntDTypeT = TypeVar("_NonIntDTypeT", bound=np.dtype[npc.inexact | np.flexible | np.datetime64 | np.timedelta64 | np.object_])

_Axis: TypeAlias = L[-2, -1, 0, 1] | bool | np.bool_ | npc.integer
_ShapeLike: TypeAlias = Iterable[op.CanIndex]
_ScalarLike: TypeAlias = complex | bytes | str | np.generic | onp.Array0D
_SequenceLike: TypeAlias = tuple[_ScalarLike, ...] | list[_ScalarLike] | onp.Array1D
_MatrixLike: TypeAlias = tuple[_SequenceLike, ...] | list[_SequenceLike] | onp.Array2D

_IntP: TypeAlias = np.int32 | np.int64
_UIntP: TypeAlias = np.uint32 | np.uint64

@type_check_only
class _ReshapeKwargs(TypedDict, total=False):
    order: L["C", "F"]
    copy: bool

@type_check_only
class _SizedIndexIterable(Protocol):
    def __len__(self, /) -> int: ...
    def __iter__(self, /) -> op.CanNext[op.CanIndex]: ...

###

supported_dtypes: Final[list[type[npc.number | np.bool_]]] = ...

#
# NOTE: Technically any `numpy.generic` could be returned, but we only care about the supported scalar types in `scipy.sparse`.
def upcast(*args: npt.DTypeLike) -> npc.number | np.bool_: ...
def upcast_char(*args: npt.DTypeLike) -> npc.number | np.bool_: ...
@overload
def upcast_scalar(dtype: onp.ToDType[_ScalarT], scalar: onp.ToScalar) -> np.dtype[_ScalarT]: ...
@overload
def upcast_scalar(dtype: npt.DTypeLike, scalar: onp.ToScalar) -> np.dtype[Any]: ...

#
def downcast_intp_index(
    arr: onp.Array[_ShapeT, np.bool_ | npc.integer | npc.floating | np.timedelta64 | np.object_],
) -> onp.Array[_ShapeT, _IntP]: ...

#
@overload
def to_native(A: _ScalarT) -> onp.Array0D[_ScalarT]: ...
@overload
def to_native(A: onp.Array[_ShapeT, _ScalarT]) -> onp.Array[_ShapeT, _ScalarT]: ...
@overload
def to_native(A: onp.HasDType[_DTypeT]) -> np.ndarray[Any, _DTypeT]: ...

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

_CoInt32: TypeAlias = np.bool_ | np.int8 | np.uint8 | np.int16 | np.uint16 | np.int32
_ContraInt32: TypeAlias = np.uint32 | np.int64 | np.uint64

#
@overload
def get_index_dtype(
    arrays: tuple[()] = (), maxval: onp.ToFloat | None = None, check_contents: op.CanBool = False
) -> type[np.int32]: ...
@overload
def get_index_dtype(
    arrays: tuple[onp.CanArrayND[_CoInt32], *tuple[onp.CanArrayND[_CoInt32], ...]],
    maxval: onp.ToFloat | None = None,
    check_contents: op.CanBool = False,
) -> type[np.int32]: ...
@overload
def get_index_dtype(
    arrays: tuple[onp.CanArrayND[_ContraInt32], *tuple[onp.CanArrayND[_ContraInt32], ...]],
    maxval: onp.ToFloat | None = None,
    check_contents: op.CanBool = False,
) -> type[np.int64]: ...
@overload
def get_index_dtype(
    arrays: tuple[onp.ToInt | onp.ToIntND, ...], maxval: onp.ToFloat | None = None, check_contents: op.CanBool = False
) -> type[_IntP]: ...

# NOTE: The inline annotations (`(np.dtype) -> np.dtype`) are incorrect.
@overload
def get_sum_dtype(dtype: np.dtype[npc.unsignedinteger]) -> type[_UIntP]: ...
@overload
def get_sum_dtype(dtype: np.dtype[np.bool_ | npc.signedinteger]) -> type[_IntP]: ...
@overload
def get_sum_dtype(dtype: _NonIntDTypeT) -> _NonIntDTypeT: ...

#
# NOTE: all arrays implement `__index__` but if it raises this returns `False`, so `TypeIs` can't be used here
def isintlike(x: object) -> TypeIs[op.CanIndex]: ...
def isscalarlike(x: object) -> TypeIs[_ScalarLike]: ...
def isshape(x: _SizedIndexIterable, nonneg: bool = False, *, allow_nd: tuple[int, ...] = (2,)) -> bool: ...
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
def check_reshape_kwargs(kwargs: _ReshapeKwargs) -> L["C", "F"] | bool: ...

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
def safely_cast_index_arrays(
    A: bsr_array | bsr_matrix | csc_array | csc_matrix | csr_array | csr_matrix, idx_dtype: onp.ToDType[_IntT], msg: str = ""
) -> tuple[onp.Array1D[_IntT], onp.Array1D[_IntT]]: ...
@overload  # 2d COO, dtype: <default>
def safely_cast_index_arrays(
    A: coo_array[Any, tuple[int, int]] | coo_matrix,
    idx_dtype: onp.ToDType[np.int32] = ...,  # = np.int32
    msg: str = "",
) -> tuple[onp.Array1D[np.int32], onp.Array1D[np.int32]]: ...
@overload  # 2d COO, dtype: <known>
def safely_cast_index_arrays(
    A: coo_array[Any, tuple[int, int]] | coo_matrix, idx_dtype: onp.ToDType[_IntT], msg: str = ""
) -> tuple[onp.Array1D[_IntT], onp.Array1D[_IntT]]: ...
@overload  # nd COO, dtype: <default>
def safely_cast_index_arrays(
    A: coo_array,
    idx_dtype: onp.ToDType[np.int32] = ...,  # = np.int32
    msg: str = "",
) -> tuple[onp.Array1D[np.int32], ...]: ...
@overload  # nd COO, dtype: <known>
def safely_cast_index_arrays(A: coo_array, idx_dtype: onp.ToDType[_IntT], msg: str = "") -> tuple[onp.Array1D[_IntT], ...]: ...
@overload  # DIA, dtype: <default>
def safely_cast_index_arrays(
    A: dia_array | dia_matrix,
    idx_dtype: onp.ToDType[np.int32] = ...,  # = np.int32
    msg: str = "",
) -> onp.Array1D[np.int32]: ...
@overload  # DIA, dtype: <known>
def safely_cast_index_arrays(A: dia_array | dia_matrix, idx_dtype: onp.ToDType[_IntT], msg: str = "") -> onp.Array1D[_IntT]: ...

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
