from __future__ import annotations

from collections.abc import Callable, Hashable, Iterable, Mapping, Sequence
from enum import Enum
from types import EllipsisType, ModuleType
from typing import (
    TYPE_CHECKING,
    Any,
    Final,
    Literal,
    Protocol,
    SupportsIndex,
    TypeVar,
    Union,
    overload,
    runtime_checkable,
)

import numpy as np

try:
    from typing import TypeAlias
except ImportError:
    if TYPE_CHECKING:
        raise
    else:
        Self: Any = None


# Singleton type, as per https://github.com/python/typing/pull/240
class Default(Enum):
    token: Final = 0


_default = Default.token

# https://stackoverflow.com/questions/74633074/how-to-type-hint-a-generic-numpy-array
_T_co = TypeVar("_T_co", covariant=True)

_dtype = np.dtype
_DType = TypeVar("_DType", bound=np.dtype[Any])
_DType_co = TypeVar("_DType_co", covariant=True, bound=np.dtype[Any])
# A subset of `npt.DTypeLike` that can be parametrized w.r.t. `np.generic`

_ScalarType = TypeVar("_ScalarType", bound=np.generic)
_ScalarType_co = TypeVar("_ScalarType_co", bound=np.generic, covariant=True)


# A protocol for anything with the dtype attribute
@runtime_checkable
class _SupportsDType(Protocol[_DType_co]):
    @property
    def dtype(self) -> _DType_co: ...


_DTypeLike = Union[
    np.dtype[_ScalarType],
    type[_ScalarType],
    _SupportsDType[np.dtype[_ScalarType]],
]

# For unknown shapes Dask uses np.nan, array_api uses None:
_IntOrUnknown = int
_Shape = tuple[_IntOrUnknown, ...]
_ShapeLike = Union[SupportsIndex, Sequence[SupportsIndex]]
_ShapeType = TypeVar("_ShapeType", bound=Any)
_ShapeType_co = TypeVar("_ShapeType_co", bound=Any, covariant=True)

_Axis = int
_Axes = tuple[_Axis, ...]
_AxisLike = Union[_Axis, _Axes]

_Chunks = tuple[_Shape, ...]
_NormalizedChunks = tuple[tuple[int, ...], ...]
# FYI in some cases we don't allow `None`, which this doesn't take account of.
# # FYI the `str` is for a size string, e.g. "16MB", supported by dask.
T_ChunkDim: TypeAlias = str | int | Literal["auto"] | tuple[int, ...] | None  # noqa: PYI051
# We allow the tuple form of this (though arguably we could transition to named dims only)
T_Chunks: TypeAlias = T_ChunkDim | Mapping[Any, T_ChunkDim]

_Dim = Hashable
_Dims = tuple[_Dim, ...]

_DimsLike = Union[str, Iterable[_Dim]]

# https://data-apis.org/array-api/latest/API_specification/indexing.html
# TODO: np.array_api was bugged and didn't allow (None,), but should!
# https://github.com/numpy/numpy/pull/25022
# https://github.com/data-apis/array-api/pull/674
_IndexKey = Union[int, slice, EllipsisType]
_IndexKeys = tuple[_IndexKey, ...]  #  tuple[Union[_IndexKey, None], ...]
_IndexKeyLike = Union[_IndexKey, _IndexKeys]

_AttrsLike = Union[Mapping[Any, Any], None]


class _SupportsReal(Protocol[_T_co]):
    @property
    def real(self) -> _T_co: ...


class _SupportsImag(Protocol[_T_co]):
    @property
    def imag(self) -> _T_co: ...


@runtime_checkable
class _array(Protocol[_ShapeType_co, _DType_co]):
    """
    Minimal duck array named array uses.

    Corresponds to np.ndarray.
    """

    @property
    def shape(self) -> _Shape: ...

    @property
    def dtype(self) -> _DType_co: ...


@runtime_checkable
class _arrayfunction(
    _array[_ShapeType_co, _DType_co], Protocol[_ShapeType_co, _DType_co]
):
    """
    Duck array supporting NEP 18.

    Corresponds to np.ndarray.
    """

    @overload
    def __getitem__(
        self, key: _arrayfunction[Any, Any] | tuple[_arrayfunction[Any, Any], ...], /
    ) -> _arrayfunction[Any, _DType_co]: ...

    @overload
    def __getitem__(self, key: _IndexKeyLike, /) -> Any: ...

    def __getitem__(
        self,
        key: (
            _IndexKeyLike
            | _arrayfunction[Any, Any]
            | tuple[_arrayfunction[Any, Any], ...]
        ),
        /,
    ) -> _arrayfunction[Any, _DType_co] | Any: ...

    @overload
    def __array__(
        self, dtype: None = ..., /, *, copy: bool | None = ...
    ) -> np.ndarray[Any, _DType_co]: ...
    @overload
    def __array__(
        self, dtype: _DType, /, *, copy: bool | None = ...
    ) -> np.ndarray[Any, _DType]: ...

    def __array__(
        self, dtype: _DType | None = ..., /, *, copy: bool | None = ...
    ) -> np.ndarray[Any, _DType] | np.ndarray[Any, _DType_co]: ...

    # TODO: Should return the same subclass but with a new dtype generic.
    # https://github.com/python/typing/issues/548
    def __array_ufunc__(
        self,
        ufunc: Any,
        method: Any,
        *inputs: Any,
        **kwargs: Any,
    ) -> Any: ...

    # TODO: Should return the same subclass but with a new dtype generic.
    # https://github.com/python/typing/issues/548
    def __array_function__(
        self,
        func: Callable[..., Any],
        types: Iterable[type],
        args: Iterable[Any],
        kwargs: Mapping[str, Any],
    ) -> Any: ...

    @property
    def imag(self) -> _arrayfunction[_ShapeType_co, Any]: ...

    @property
    def real(self) -> _arrayfunction[_ShapeType_co, Any]: ...


@runtime_checkable
class _arrayapi(_array[_ShapeType_co, _DType_co], Protocol[_ShapeType_co, _DType_co]):
    """
    Duck array supporting NEP 47.

    Corresponds to np.ndarray.
    """

    def __getitem__(
        self,
        key: (
            _IndexKeyLike | Any
        ),  # TODO: Any should be _arrayapi[Any, _dtype[np.integer]]
        /,
    ) -> _arrayapi[Any, Any]: ...

    def __array_namespace__(self) -> ModuleType: ...


# NamedArray can most likely use both __array_function__ and __array_namespace__:
_arrayfunction_or_api = (_arrayfunction, _arrayapi)

duckarray = Union[
    _arrayfunction[_ShapeType_co, _DType_co], _arrayapi[_ShapeType_co, _DType_co]
]

# Corresponds to np.typing.NDArray:
DuckArray = _arrayfunction[Any, np.dtype[_ScalarType_co]]


@runtime_checkable
class _chunkedarray(
    _array[_ShapeType_co, _DType_co], Protocol[_ShapeType_co, _DType_co]
):
    """
    Minimal chunked duck array.

    Corresponds to np.ndarray.
    """

    @property
    def chunks(self) -> _Chunks: ...


@runtime_checkable
class _chunkedarrayfunction(
    _arrayfunction[_ShapeType_co, _DType_co], Protocol[_ShapeType_co, _DType_co]
):
    """
    Chunked duck array supporting NEP 18.

    Corresponds to np.ndarray.
    """

    @property
    def chunks(self) -> _Chunks: ...


@runtime_checkable
class _chunkedarrayapi(
    _arrayapi[_ShapeType_co, _DType_co], Protocol[_ShapeType_co, _DType_co]
):
    """
    Chunked duck array supporting NEP 47.

    Corresponds to np.ndarray.
    """

    @property
    def chunks(self) -> _Chunks: ...


# NamedArray can most likely use both __array_function__ and __array_namespace__:
_chunkedarrayfunction_or_api = (_chunkedarrayfunction, _chunkedarrayapi)
chunkedduckarray = Union[
    _chunkedarrayfunction[_ShapeType_co, _DType_co],
    _chunkedarrayapi[_ShapeType_co, _DType_co],
]


@runtime_checkable
class _sparsearray(
    _array[_ShapeType_co, _DType_co], Protocol[_ShapeType_co, _DType_co]
):
    """
    Minimal sparse duck array.

    Corresponds to np.ndarray.
    """

    def todense(self) -> np.ndarray[Any, _DType_co]: ...


@runtime_checkable
class _sparsearrayfunction(
    _arrayfunction[_ShapeType_co, _DType_co], Protocol[_ShapeType_co, _DType_co]
):
    """
    Sparse duck array supporting NEP 18.

    Corresponds to np.ndarray.
    """

    def todense(self) -> np.ndarray[Any, _DType_co]: ...


@runtime_checkable
class _sparsearrayapi(
    _arrayapi[_ShapeType_co, _DType_co], Protocol[_ShapeType_co, _DType_co]
):
    """
    Sparse duck array supporting NEP 47.

    Corresponds to np.ndarray.
    """

    def todense(self) -> np.ndarray[Any, _DType_co]: ...


# NamedArray can most likely use both __array_function__ and __array_namespace__:
_sparsearrayfunction_or_api = (_sparsearrayfunction, _sparsearrayapi)
sparseduckarray = Union[
    _sparsearrayfunction[_ShapeType_co, _DType_co],
    _sparsearrayapi[_ShapeType_co, _DType_co],
]

ErrorOptions = Literal["raise", "ignore"]
ErrorOptionsWithWarn = Literal["raise", "warn", "ignore"]
