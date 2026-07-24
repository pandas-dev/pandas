"""
This module contains various lazy array classes which can be wrapped and manipulated by xarray objects but will raise on data access.
"""

from collections.abc import Callable, Iterable
from typing import Any, Self

import numpy as np

from xarray.core import utils
from xarray.core.indexing import ExplicitlyIndexed


class UnexpectedDataAccess(Exception):
    pass


class InaccessibleArray(utils.NDArrayMixin, ExplicitlyIndexed):
    """Disallows any loading."""

    def __init__(self, array):
        self.array = array

    def get_duck_array(self):
        raise UnexpectedDataAccess("Tried accessing data")

    def __array__(
        self, dtype: np.typing.DTypeLike | None = None, /, *, copy: bool | None = None
    ) -> np.ndarray:
        raise UnexpectedDataAccess("Tried accessing data")

    def __getitem__(self, key):
        raise UnexpectedDataAccess("Tried accessing data.")


class FirstElementAccessibleArray(InaccessibleArray):
    def __getitem__(self, key):
        tuple_idxr = key.tuple
        if len(tuple_idxr) > 1:
            raise UnexpectedDataAccess("Tried accessing more than one element.")
        return self.array[tuple_idxr]


class IndexableArray(InaccessibleArray):
    """An InaccessibleArray subclass that supports indexing."""

    def __getitem__(self, key):
        return type(self)(self.array[key])

    def transpose(self, axes):
        return type(self)(self.array.transpose(axes))


class DuckArrayWrapper(utils.NDArrayMixin):
    """Array-like that prevents casting to array.
    Modeled after cupy."""

    def __init__(self, array: np.ndarray):
        self.array = array

    def __getitem__(self, key):
        return type(self)(self.array[key])

    def to_numpy(self) -> np.ndarray:
        """Allow explicit conversions to numpy in `to_numpy`, but disallow np.asarray etc."""
        return self.array

    def __array__(
        self, dtype: np.typing.DTypeLike | None = None, /, *, copy: bool | None = None
    ) -> np.ndarray:
        raise UnexpectedDataAccess("Tried accessing data")

    def __array_namespace__(self):
        """Present to satisfy is_duck_array test."""
        from xarray.tests import namespace

        return namespace


CONCATENATABLEARRAY_HANDLED_ARRAY_FUNCTIONS: dict[str, Callable] = {}


def implements(numpy_function):
    """Register an __array_function__ implementation for ConcatenatableArray objects."""

    def decorator(func):
        CONCATENATABLEARRAY_HANDLED_ARRAY_FUNCTIONS[numpy_function] = func
        return func

    return decorator


@implements(np.concatenate)
def concatenate(
    arrays: Iterable["ConcatenatableArray"], /, *, axis=0
) -> "ConcatenatableArray":
    if any(not isinstance(arr, ConcatenatableArray) for arr in arrays):
        raise TypeError

    result = np.concatenate([arr._array for arr in arrays], axis=axis)
    return ConcatenatableArray(result)


@implements(np.stack)
def stack(
    arrays: Iterable["ConcatenatableArray"], /, *, axis=0
) -> "ConcatenatableArray":
    if any(not isinstance(arr, ConcatenatableArray) for arr in arrays):
        raise TypeError

    result = np.stack([arr._array for arr in arrays], axis=axis)
    return ConcatenatableArray(result)


@implements(np.result_type)
def result_type(*arrays_and_dtypes) -> np.dtype:
    """Called by xarray to ensure all arguments to concat have the same dtype."""
    first_dtype, *other_dtypes = (np.dtype(obj) for obj in arrays_and_dtypes)
    for other_dtype in other_dtypes:
        if other_dtype != first_dtype:
            raise ValueError("dtypes not all consistent")
    return first_dtype


@implements(np.broadcast_to)
def broadcast_to(
    x: "ConcatenatableArray", /, shape: tuple[int, ...]
) -> "ConcatenatableArray":
    """
    Broadcasts an array to a specified shape, by either manipulating chunk keys or copying chunk manifest entries.
    """
    if not isinstance(x, ConcatenatableArray):
        raise TypeError

    result = np.broadcast_to(x._array, shape=shape)
    return ConcatenatableArray(result)


@implements(np.full_like)
def full_like(
    x: "ConcatenatableArray", /, fill_value, **kwargs
) -> "ConcatenatableArray":
    """
    Broadcasts an array to a specified shape, by either manipulating chunk keys or copying chunk manifest entries.
    """
    if not isinstance(x, ConcatenatableArray):
        raise TypeError
    return ConcatenatableArray(np.full(x.shape, fill_value=fill_value, **kwargs))


@implements(np.all)
def numpy_all(x: "ConcatenatableArray", **kwargs) -> "ConcatenatableArray":
    return type(x)(np.all(x._array, **kwargs))


class ConcatenatableArray:
    """Disallows loading or coercing to an index but does support concatenation / stacking."""

    def __init__(self, array):
        # use ._array instead of .array because we don't want this to be accessible even to xarray's internals (e.g. create_default_index_implicit)
        self._array = array

    @property
    def dtype(self: Any) -> np.dtype:
        return self._array.dtype

    @property
    def shape(self: Any) -> tuple[int, ...]:
        return self._array.shape

    @property
    def ndim(self: Any) -> int:
        return self._array.ndim

    def __repr__(self: Any) -> str:
        return f"{type(self).__name__}(array={self._array!r})"

    def get_duck_array(self):
        raise UnexpectedDataAccess("Tried accessing data")

    def __array__(
        self, dtype: np.typing.DTypeLike | None = None, /, *, copy: bool | None = None
    ) -> np.ndarray:
        raise UnexpectedDataAccess("Tried accessing data")

    def __getitem__(self, key) -> Self:
        """Some cases of concat require supporting expanding dims by dimensions of size 1"""
        # see https://data-apis.org/array-api/2022.12/API_specification/indexing.html#multi-axis-indexing
        arr = self._array
        for axis, indexer_1d in enumerate(key):
            if indexer_1d is None:
                arr = np.expand_dims(arr, axis)
            elif indexer_1d is Ellipsis:
                pass
            else:
                raise UnexpectedDataAccess("Tried accessing data.")
        return type(self)(arr)

    def __eq__(self, other: Self) -> Self:  # type: ignore[override]
        return type(self)(self._array == other._array)

    def __array_function__(self, func, types, args, kwargs) -> Any:
        if func not in CONCATENATABLEARRAY_HANDLED_ARRAY_FUNCTIONS:
            return NotImplemented

        # Note: this allows subclasses that don't override
        # __array_function__ to handle ManifestArray objects
        if not all(issubclass(t, ConcatenatableArray) for t in types):
            return NotImplemented

        return CONCATENATABLEARRAY_HANDLED_ARRAY_FUNCTIONS[func](*args, **kwargs)

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs) -> Any:
        """We have to define this in order to convince xarray that this class is a duckarray, even though we will never support ufuncs."""
        return NotImplemented

    def astype(self, dtype: np.dtype, /, *, copy: bool = True) -> Self:
        """Needed because xarray will call this even when it's a no-op"""
        if dtype != self.dtype:
            raise NotImplementedError()
        else:
            return self

    def __and__(self, other: Self) -> Self:
        return type(self)(self._array & other._array)

    def __or__(self, other: Self) -> Self:
        return type(self)(self._array | other._array)
