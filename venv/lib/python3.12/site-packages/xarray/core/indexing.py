from __future__ import annotations

import enum
import functools
import math
import operator
from collections import Counter, defaultdict
from collections.abc import Callable, Hashable, Iterable, Mapping
from contextlib import suppress
from dataclasses import dataclass, field
from datetime import timedelta
from typing import TYPE_CHECKING, Any, cast, overload

import numpy as np
import pandas as pd
from numpy.typing import DTypeLike
from packaging.version import Version

from xarray.core import duck_array_ops
from xarray.core.coordinate_transform import CoordinateTransform
from xarray.core.nputils import NumpyVIndexAdapter
from xarray.core.types import T_Xarray
from xarray.core.utils import (
    NDArrayMixin,
    either_dict_or_kwargs,
    get_valid_numpy_dtype,
    is_allowed_extension_array,
    is_allowed_extension_array_dtype,
    is_duck_array,
    is_duck_dask_array,
    is_full_slice,
    is_scalar,
    is_valid_numpy_dtype,
    to_0d_array,
)
from xarray.namedarray.parallelcompat import get_chunked_array_type
from xarray.namedarray.pycompat import array_type, integer_types, is_chunked_array

if TYPE_CHECKING:
    from xarray.core.extension_array import PandasExtensionArray
    from xarray.core.indexes import Index
    from xarray.core.types import Self
    from xarray.core.variable import Variable
    from xarray.namedarray._typing import _Shape, duckarray
    from xarray.namedarray.parallelcompat import ChunkManagerEntrypoint

BasicIndexerType = int | np.integer | slice
OuterIndexerType = BasicIndexerType | np.ndarray[Any, np.dtype[np.integer]]


@dataclass
class IndexSelResult:
    """Index query results.

    Attributes
    ----------
    dim_indexers: dict
        A dictionary where keys are array dimensions and values are
        location-based indexers.
    indexes: dict, optional
        New indexes to replace in the resulting DataArray or Dataset.
    variables : dict, optional
        New variables to replace in the resulting DataArray or Dataset.
    drop_coords : list, optional
        Coordinate(s) to drop in the resulting DataArray or Dataset.
    drop_indexes : list, optional
        Index(es) to drop in the resulting DataArray or Dataset.
    rename_dims : dict, optional
        A dictionary in the form ``{old_dim: new_dim}`` for dimension(s) to
        rename in the resulting DataArray or Dataset.

    """

    dim_indexers: dict[Any, Any]
    indexes: dict[Any, Index] = field(default_factory=dict)
    variables: dict[Any, Variable] = field(default_factory=dict)
    drop_coords: list[Hashable] = field(default_factory=list)
    drop_indexes: list[Hashable] = field(default_factory=list)
    rename_dims: dict[Any, Hashable] = field(default_factory=dict)

    def as_tuple(self):
        """Unlike ``dataclasses.astuple``, return a shallow copy.

        See https://stackoverflow.com/a/51802661

        """
        return (
            self.dim_indexers,
            self.indexes,
            self.variables,
            self.drop_coords,
            self.drop_indexes,
            self.rename_dims,
        )


def merge_sel_results(results: list[IndexSelResult]) -> IndexSelResult:
    all_dims_count = Counter([dim for res in results for dim in res.dim_indexers])
    duplicate_dims = {k: v for k, v in all_dims_count.items() if v > 1}

    if duplicate_dims:
        # TODO: this message is not right when combining indexe(s) queries with
        # location-based indexing on a dimension with no dimension-coordinate (failback)
        fmt_dims = [
            f"{dim!r}: {count} indexes involved"
            for dim, count in duplicate_dims.items()
        ]
        raise ValueError(
            "Xarray does not support label-based selection with more than one index "
            "over the following dimension(s):\n"
            + "\n".join(fmt_dims)
            + "\nSuggestion: use a multi-index for each of those dimension(s)."
        )

    dim_indexers = {}
    indexes = {}
    variables = {}
    drop_coords = []
    drop_indexes = []
    rename_dims = {}

    for res in results:
        dim_indexers.update(res.dim_indexers)
        indexes.update(res.indexes)
        variables.update(res.variables)
        drop_coords += res.drop_coords
        drop_indexes += res.drop_indexes
        rename_dims.update(res.rename_dims)

    return IndexSelResult(
        dim_indexers, indexes, variables, drop_coords, drop_indexes, rename_dims
    )


def group_indexers_by_index(
    obj: T_Xarray,
    indexers: Mapping[Any, Any],
    options: Mapping[str, Any],
) -> list[tuple[Index, dict[Any, Any]]]:
    """Returns a list of unique indexes and their corresponding indexers."""
    unique_indexes = {}
    grouped_indexers: Mapping[int | None, dict] = defaultdict(dict)

    for key, label in indexers.items():
        index: Index = obj.xindexes.get(key, None)

        if index is not None:
            index_id = id(index)
            unique_indexes[index_id] = index
            grouped_indexers[index_id][key] = label
        elif key in obj.coords:
            raise KeyError(f"no index found for coordinate {key!r}")
        elif key not in obj.dims:
            raise KeyError(
                f"{key!r} is not a valid dimension or coordinate for "
                f"{obj.__class__.__name__} with dimensions {obj.dims!r}"
            )
        elif len(options):
            raise ValueError(
                f"cannot supply selection options {options!r} for dimension {key!r}"
                "that has no associated coordinate or index"
            )
        else:
            # key is a dimension without a "dimension-coordinate"
            # failback to location-based selection
            # TODO: depreciate this implicit behavior and suggest using isel instead?
            unique_indexes[None] = None
            grouped_indexers[None][key] = label

    return [(unique_indexes[k], grouped_indexers[k]) for k in unique_indexes]


def map_index_queries(
    obj: T_Xarray,
    indexers: Mapping[Any, Any],
    method=None,
    tolerance: int | float | Iterable[int | float] | None = None,
    **indexers_kwargs: Any,
) -> IndexSelResult:
    """Execute index queries from a DataArray / Dataset and label-based indexers
    and return the (merged) query results.

    """
    from xarray.core.dataarray import DataArray

    # TODO benbovy - flexible indexes: remove when custom index options are available
    if method is None and tolerance is None:
        options = {}
    else:
        options = {"method": method, "tolerance": tolerance}

    indexers = either_dict_or_kwargs(indexers, indexers_kwargs, "map_index_queries")
    grouped_indexers = group_indexers_by_index(obj, indexers, options)

    results = []
    for index, labels in grouped_indexers:
        if index is None:
            # forward dimension indexers with no index/coordinate
            results.append(IndexSelResult(labels))
        else:
            results.append(index.sel(labels, **options))

    merged = merge_sel_results(results)

    # drop dimension coordinates found in dimension indexers
    # (also drop multi-index if any)
    # (.sel() already ensures alignment)
    for k, v in merged.dim_indexers.items():
        if isinstance(v, DataArray):
            if k in v._indexes:
                v = v.reset_index(k)
            drop_coords = [name for name in v._coords if name in merged.dim_indexers]
            merged.dim_indexers[k] = v.drop_vars(drop_coords)

    return merged


def expanded_indexer(key, ndim):
    """Given a key for indexing an ndarray, return an equivalent key which is a
    tuple with length equal to the number of dimensions.

    The expansion is done by replacing all `Ellipsis` items with the right
    number of full slices and then padding the key with full slices so that it
    reaches the appropriate dimensionality.
    """
    if not isinstance(key, tuple):
        # numpy treats non-tuple keys equivalent to tuples of length 1
        key = (key,)
    new_key = []
    # handling Ellipsis right is a little tricky, see:
    # https://numpy.org/doc/stable/reference/arrays.indexing.html#advanced-indexing
    found_ellipsis = False
    for k in key:
        if k is Ellipsis:
            if not found_ellipsis:
                new_key.extend((ndim + 1 - len(key)) * [slice(None)])
                found_ellipsis = True
            else:
                new_key.append(slice(None))
        else:
            new_key.append(k)
    if len(new_key) > ndim:
        raise IndexError("too many indices")
    new_key.extend((ndim - len(new_key)) * [slice(None)])
    return tuple(new_key)


def normalize_slice(sl: slice, size: int) -> slice:
    """
    Ensure that given slice only contains positive start and stop values
    (stop can be -1 for full-size slices with negative steps, e.g. [-10::-1])

    Examples
    --------
    >>> normalize_slice(slice(0, 9), 10)
    slice(0, 9, 1)
    >>> normalize_slice(slice(0, -1), 10)
    slice(0, 9, 1)
    """
    return slice(*sl.indices(size))


def _expand_slice(slice_: slice, size: int) -> np.ndarray[Any, np.dtype[np.integer]]:
    """
    Expand slice to an array containing only positive integers.

    Examples
    --------
    >>> _expand_slice(slice(0, 9), 10)
    array([0, 1, 2, 3, 4, 5, 6, 7, 8])
    >>> _expand_slice(slice(0, -1), 10)
    array([0, 1, 2, 3, 4, 5, 6, 7, 8])
    """
    sl = normalize_slice(slice_, size)
    return np.arange(sl.start, sl.stop, sl.step)


def slice_slice(old_slice: slice, applied_slice: slice, size: int) -> slice:
    """Given a slice and the size of the dimension to which it will be applied,
    index it with another slice to return a new slice equivalent to applying
    the slices sequentially
    """
    old_slice = normalize_slice(old_slice, size)

    size_after_old_slice = len(range(old_slice.start, old_slice.stop, old_slice.step))
    if size_after_old_slice == 0:
        # nothing left after applying first slice
        return slice(0)

    applied_slice = normalize_slice(applied_slice, size_after_old_slice)

    start = old_slice.start + applied_slice.start * old_slice.step
    if start < 0:
        # nothing left after applying second slice
        # (can only happen for old_slice.step < 0, e.g. [10::-1], [20:])
        return slice(0)

    stop = old_slice.start + applied_slice.stop * old_slice.step
    if stop < 0:
        stop = None

    step = old_slice.step * applied_slice.step

    return slice(start, stop, step)


def normalize_array(
    array: np.ndarray[Any, np.dtype[np.integer]], size: int
) -> np.ndarray[Any, np.dtype[np.integer]]:
    """
    Ensure that the given array only contains positive values.

    Examples
    --------
    >>> normalize_array(np.array([-1, -2, -3, -4]), 10)
    array([9, 8, 7, 6])
    >>> normalize_array(np.array([-5, 3, 5, -1, 8]), 12)
    array([ 7,  3,  5, 11,  8])
    """
    if np.issubdtype(array.dtype, np.unsignedinteger):
        return array

    return np.where(array >= 0, array, array + size)


def slice_slice_by_array(
    old_slice: slice,
    array: np.ndarray[Any, np.dtype[np.integer]],
    size: int,
) -> np.ndarray[Any, np.dtype[np.integer]]:
    """Given a slice and the size of the dimension to which it will be applied,
    index it with an array to return a new array equivalent to applying
    the slices sequentially

    Examples
    --------
    >>> slice_slice_by_array(slice(2, 10), np.array([1, 3, 5]), 12)
    array([3, 5, 7])
    >>> slice_slice_by_array(slice(1, None, 2), np.array([1, 3, 7, 8]), 20)
    array([ 3,  7, 15, 17])
    >>> slice_slice_by_array(slice(None, None, -1), np.array([2, 4, 7]), 20)
    array([17, 15, 12])
    """
    # to get a concrete slice, limited to the size of the array
    normalized_slice = normalize_slice(old_slice, size)

    size_after_slice = len(range(*normalized_slice.indices(size)))
    normalized_array = normalize_array(array, size_after_slice)

    new_indexer = normalized_array * normalized_slice.step + normalized_slice.start

    if np.any(new_indexer >= size):
        raise IndexError("indices out of bounds")  # TODO: more helpful error message

    return new_indexer


def _index_indexer_1d(
    old_indexer: OuterIndexerType,
    applied_indexer: OuterIndexerType,
    size: int,
) -> OuterIndexerType:
    if is_full_slice(applied_indexer):
        # shortcut for the usual case
        return old_indexer
    if is_full_slice(old_indexer):
        # shortcut for full slices
        return applied_indexer

    indexer: OuterIndexerType
    if isinstance(old_indexer, slice):
        if isinstance(applied_indexer, slice):
            indexer = slice_slice(old_indexer, applied_indexer, size)
        elif isinstance(applied_indexer, integer_types):
            indexer = range(*old_indexer.indices(size))[applied_indexer]
        else:
            indexer = slice_slice_by_array(old_indexer, applied_indexer, size)
    elif isinstance(old_indexer, np.ndarray):
        indexer = old_indexer[applied_indexer]
    else:
        # should be unreachable
        raise ValueError("cannot index integers. Please open an issuec-")

    return indexer


class ExplicitIndexer:
    """Base class for explicit indexer objects.

    ExplicitIndexer objects wrap a tuple of values given by their ``tuple``
    property. These tuples should always have length equal to the number of
    dimensions on the indexed array.

    Do not instantiate BaseIndexer objects directly: instead, use one of the
    sub-classes BasicIndexer, OuterIndexer or VectorizedIndexer.
    """

    __slots__ = ("_key",)

    def __init__(self, key: tuple[Any, ...]):
        if type(self) is ExplicitIndexer:
            raise TypeError("cannot instantiate base ExplicitIndexer objects")
        self._key = tuple(key)

    @property
    def tuple(self) -> tuple[Any, ...]:
        return self._key

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.tuple})"


@overload
def as_integer_or_none(value: int) -> int: ...
@overload
def as_integer_or_none(value: None) -> None: ...
def as_integer_or_none(value: int | None) -> int | None:
    return None if value is None else operator.index(value)


def as_integer_slice(value: slice) -> slice:
    start = as_integer_or_none(value.start)
    stop = as_integer_or_none(value.stop)
    step = as_integer_or_none(value.step)
    return slice(start, stop, step)


class IndexCallable:
    """Provide getitem and setitem syntax for callable objects."""

    __slots__ = ("getter", "setter")

    def __init__(
        self, getter: Callable[..., Any], setter: Callable[..., Any] | None = None
    ):
        self.getter = getter
        self.setter = setter

    def __getitem__(self, key: Any) -> Any:
        return self.getter(key)

    def __setitem__(self, key: Any, value: Any) -> None:
        if self.setter is None:
            raise NotImplementedError(
                "Setting values is not supported for this indexer."
            )
        self.setter(key, value)


class BasicIndexer(ExplicitIndexer):
    """Tuple for basic indexing.

    All elements should be int or slice objects. Indexing follows NumPy's
    rules for basic indexing: each axis is independently sliced and axes
    indexed with an integer are dropped from the result.
    """

    __slots__ = ()

    def __init__(self, key: tuple[BasicIndexerType, ...]):
        if not isinstance(key, tuple):
            raise TypeError(f"key must be a tuple: {key!r}")

        new_key = []
        for k in key:
            if isinstance(k, integer_types):
                k = int(k)
            elif isinstance(k, slice):
                k = as_integer_slice(k)
            else:
                raise TypeError(
                    f"unexpected indexer type for {type(self).__name__}: {k!r}"
                )
            new_key.append(k)

        super().__init__(tuple(new_key))


class OuterIndexer(ExplicitIndexer):
    """Tuple for outer/orthogonal indexing.

    All elements should be int, slice or 1-dimensional np.ndarray objects with
    an integer dtype. Indexing is applied independently along each axis, and
    axes indexed with an integer are dropped from the result. This type of
    indexing works like MATLAB/Fortran.
    """

    __slots__ = ()

    def __init__(
        self,
        key: tuple[BasicIndexerType | np.ndarray[Any, np.dtype[np.generic]], ...],
    ):
        if not isinstance(key, tuple):
            raise TypeError(f"key must be a tuple: {key!r}")

        new_key = []
        for k in key:
            if isinstance(k, integer_types) and not isinstance(k, bool):
                k = int(k)
            elif isinstance(k, slice):
                k = as_integer_slice(k)
            elif is_duck_array(k):
                if not np.issubdtype(k.dtype, np.integer):
                    raise TypeError(
                        f"invalid indexer array, does not have integer dtype: {k!r}"
                    )
                if k.ndim > 1:  # type: ignore[union-attr]
                    raise TypeError(
                        f"invalid indexer array for {type(self).__name__}; must be scalar "
                        f"or have 1 dimension: {k!r}"
                    )
                k = duck_array_ops.astype(k, np.int64, copy=False)
            else:
                raise TypeError(
                    f"unexpected indexer type for {type(self).__name__}: {k!r}, {type(k)}"
                )
            new_key.append(k)

        super().__init__(tuple(new_key))


class VectorizedIndexer(ExplicitIndexer):
    """Tuple for vectorized indexing.

    All elements should be slice or N-dimensional np.ndarray objects with an
    integer dtype and the same number of dimensions. Indexing follows proposed
    rules for np.ndarray.vindex, which matches NumPy's advanced indexing rules
    (including broadcasting) except sliced axes are always moved to the end:
    https://github.com/numpy/numpy/pull/6256
    """

    __slots__ = ()

    def __init__(self, key: tuple[slice | np.ndarray[Any, np.dtype[np.generic]], ...]):
        if not isinstance(key, tuple):
            raise TypeError(f"key must be a tuple: {key!r}")

        new_key = []
        ndim = None
        for k in key:
            if isinstance(k, slice):
                k = as_integer_slice(k)
            elif is_duck_array(k):
                if not np.issubdtype(k.dtype, np.integer):
                    raise TypeError(
                        f"invalid indexer array, does not have integer dtype: {k!r}"
                    )
                if ndim is None:
                    ndim = k.ndim  # type: ignore[union-attr]
                elif ndim != k.ndim:  # type: ignore[union-attr]
                    ndims = [k.ndim for k in key if isinstance(k, np.ndarray)]
                    raise ValueError(
                        "invalid indexer key: ndarray arguments "
                        f"have different numbers of dimensions: {ndims}"
                    )
                k = duck_array_ops.astype(k, np.int64, copy=False)
            else:
                raise TypeError(
                    f"unexpected indexer type for {type(self).__name__}: {k!r}"
                )
            new_key.append(k)

        super().__init__(tuple(new_key))


class ExplicitlyIndexed:
    """Mixin to mark support for Indexer subclasses in indexing."""

    __slots__ = ()

    def __array__(
        self, dtype: np.typing.DTypeLike = None, /, *, copy: bool | None = None
    ) -> np.ndarray:
        # Leave casting to an array up to the underlying array type.
        if Version(np.__version__) >= Version("2.0.0"):
            return np.asarray(self.get_duck_array(), dtype=dtype, copy=copy)
        else:
            return np.asarray(self.get_duck_array(), dtype=dtype)

    def get_duck_array(self):
        return self.array


class ExplicitlyIndexedNDArrayMixin(NDArrayMixin, ExplicitlyIndexed):
    __slots__ = ()

    def get_duck_array(self):
        raise NotImplementedError

    async def async_get_duck_array(self):
        raise NotImplementedError

    def _oindex_get(self, indexer: OuterIndexer):
        raise NotImplementedError(
            f"{self.__class__.__name__}._oindex_get method should be overridden"
        )

    def _vindex_get(self, indexer: VectorizedIndexer):
        raise NotImplementedError(
            f"{self.__class__.__name__}._vindex_get method should be overridden"
        )

    def _oindex_set(self, indexer: OuterIndexer, value: Any) -> None:
        raise NotImplementedError(
            f"{self.__class__.__name__}._oindex_set method should be overridden"
        )

    def _vindex_set(self, indexer: VectorizedIndexer, value: Any) -> None:
        raise NotImplementedError(
            f"{self.__class__.__name__}._vindex_set method should be overridden"
        )

    def _check_and_raise_if_non_basic_indexer(self, indexer: ExplicitIndexer) -> None:
        if isinstance(indexer, VectorizedIndexer | OuterIndexer):
            raise TypeError(
                "Vectorized indexing with vectorized or outer indexers is not supported. "
                "Please use .vindex and .oindex properties to index the array."
            )

    @property
    def oindex(self) -> IndexCallable:
        return IndexCallable(self._oindex_get, self._oindex_set)

    @property
    def vindex(self) -> IndexCallable:
        return IndexCallable(self._vindex_get, self._vindex_set)


class IndexingAdapter(ExplicitlyIndexedNDArrayMixin):
    """Marker class for indexing adapters.

    These classes translate between Xarray's indexing semantics and the underlying array's
    indexing semantics.
    """

    def get_duck_array(self):
        key = BasicIndexer((slice(None),) * self.ndim)
        return self[key]

    async def async_get_duck_array(self):
        """These classes are applied to in-memory arrays, so specific async support isn't needed."""
        return self.get_duck_array()


class ImplicitToExplicitIndexingAdapter(NDArrayMixin):
    """Wrap an array, converting tuples into the indicated explicit indexer."""

    __slots__ = ("array", "indexer_cls")

    def __init__(self, array, indexer_cls: type[ExplicitIndexer] = BasicIndexer):
        self.array = as_indexable(array)
        self.indexer_cls = indexer_cls

    def __array__(
        self, dtype: np.typing.DTypeLike = None, /, *, copy: bool | None = None
    ) -> np.ndarray:
        if Version(np.__version__) >= Version("2.0.0"):
            return np.asarray(self.get_duck_array(), dtype=dtype, copy=copy)
        else:
            return np.asarray(self.get_duck_array(), dtype=dtype)

    def get_duck_array(self):
        return self.array.get_duck_array()

    def __getitem__(self, key: Any):
        key = expanded_indexer(key, self.ndim)
        indexer = self.indexer_cls(key)

        result = apply_indexer(self.array, indexer)

        if isinstance(result, ExplicitlyIndexed):
            return type(self)(result, self.indexer_cls)
        else:
            # Sometimes explicitly indexed arrays return NumPy arrays or
            # scalars.
            return result


class LazilyIndexedArray(ExplicitlyIndexedNDArrayMixin):
    """Wrap an array to make basic and outer indexing lazy."""

    __slots__ = ("_shape", "array", "key")

    def __init__(self, array: Any, key: ExplicitIndexer | None = None):
        """
        Parameters
        ----------
        array : array_like
            Array like object to index.
        key : ExplicitIndexer, optional
            Array indexer. If provided, it is assumed to already be in
            canonical expanded form.
        """
        if isinstance(array, type(self)) and key is None:
            # unwrap
            key = array.key  # type: ignore[has-type]
            array = array.array  # type: ignore[has-type]

        if key is None:
            key = BasicIndexer((slice(None),) * array.ndim)

        self.array = as_indexable(array)
        self.key = key

        shape: _Shape = ()
        for size, k in zip(self.array.shape, self.key.tuple, strict=True):
            if isinstance(k, slice):
                shape += (len(range(*k.indices(size))),)
            elif isinstance(k, np.ndarray):
                shape += (k.size,)
        self._shape = shape

    def _updated_key(self, new_key: ExplicitIndexer) -> BasicIndexer | OuterIndexer:
        iter_new_key = iter(expanded_indexer(new_key.tuple, self.ndim))

        full_key: list[OuterIndexerType] = []
        for size, k in zip(self.array.shape, self.key.tuple, strict=True):
            if isinstance(k, integer_types):
                full_key.append(k)
            else:
                full_key.append(_index_indexer_1d(k, next(iter_new_key), size))
        full_key_tuple = tuple(full_key)

        if all(isinstance(k, integer_types + (slice,)) for k in full_key_tuple):
            return BasicIndexer(cast(tuple[BasicIndexerType, ...], full_key_tuple))
        return OuterIndexer(full_key_tuple)

    @property
    def shape(self) -> _Shape:
        return self._shape

    def get_duck_array(self):
        from xarray.backends.common import BackendArray

        if isinstance(self.array, BackendArray):
            array = self.array[self.key]
        else:
            array = apply_indexer(self.array, self.key)
            if isinstance(array, ExplicitlyIndexed):
                array = array.get_duck_array()
        return _wrap_numpy_scalars(array)

    async def async_get_duck_array(self):
        from xarray.backends.common import BackendArray

        if isinstance(self.array, BackendArray):
            array = await self.array.async_getitem(self.key)
        else:
            array = apply_indexer(self.array, self.key)
            if isinstance(array, ExplicitlyIndexed):
                array = await array.async_get_duck_array()
        return _wrap_numpy_scalars(array)

    def transpose(self, order):
        return LazilyVectorizedIndexedArray(self.array, self.key).transpose(order)

    def _oindex_get(self, indexer: OuterIndexer):
        return type(self)(self.array, self._updated_key(indexer))

    def _vindex_get(self, indexer: VectorizedIndexer):
        array = LazilyVectorizedIndexedArray(self.array, self.key)
        return array.vindex[indexer]

    def __getitem__(self, indexer: ExplicitIndexer):
        self._check_and_raise_if_non_basic_indexer(indexer)
        return type(self)(self.array, self._updated_key(indexer))

    def _vindex_set(self, key: VectorizedIndexer, value: Any) -> None:
        raise NotImplementedError(
            "Lazy item assignment with the vectorized indexer is not yet "
            "implemented. Load your data first by .load() or compute()."
        )

    def _oindex_set(self, key: OuterIndexer, value: Any) -> None:
        full_key = self._updated_key(key)
        self.array.oindex[full_key] = value

    def __setitem__(self, key: BasicIndexer, value: Any) -> None:
        self._check_and_raise_if_non_basic_indexer(key)
        full_key = self._updated_key(key)
        self.array[full_key] = value

    def __repr__(self) -> str:
        return f"{type(self).__name__}(array={self.array!r}, key={self.key!r})"


# keep an alias to the old name for external backends pydata/xarray#5111
LazilyOuterIndexedArray = LazilyIndexedArray


class LazilyVectorizedIndexedArray(ExplicitlyIndexedNDArrayMixin):
    """Wrap an array to make vectorized indexing lazy."""

    __slots__ = ("array", "key")

    def __init__(self, array: duckarray[Any, Any], key: ExplicitIndexer):
        """
        Parameters
        ----------
        array : array_like
            Array like object to index.
        key : VectorizedIndexer
        """
        if isinstance(key, BasicIndexer | OuterIndexer):
            self.key = _outer_to_vectorized_indexer(key, array.shape)
        elif isinstance(key, VectorizedIndexer):
            self.key = _arrayize_vectorized_indexer(key, array.shape)
        self.array = as_indexable(array)

    @property
    def shape(self) -> _Shape:
        return np.broadcast(*self.key.tuple).shape

    def get_duck_array(self):
        from xarray.backends.common import BackendArray

        if isinstance(self.array, BackendArray):
            array = self.array[self.key]
        else:
            array = apply_indexer(self.array, self.key)
            if isinstance(array, ExplicitlyIndexed):
                array = array.get_duck_array()
        return _wrap_numpy_scalars(array)

    async def async_get_duck_array(self):
        from xarray.backends.common import BackendArray

        if isinstance(self.array, BackendArray):
            array = await self.array.async_getitem(self.key)
        else:
            array = apply_indexer(self.array, self.key)
            if isinstance(array, ExplicitlyIndexed):
                array = await array.async_get_duck_array()
        return _wrap_numpy_scalars(array)

    def _updated_key(self, new_key: ExplicitIndexer):
        return _combine_indexers(self.key, self.shape, new_key)

    def _oindex_get(self, indexer: OuterIndexer):
        return type(self)(self.array, self._updated_key(indexer))

    def _vindex_get(self, indexer: VectorizedIndexer):
        return type(self)(self.array, self._updated_key(indexer))

    def __getitem__(self, indexer: ExplicitIndexer):
        self._check_and_raise_if_non_basic_indexer(indexer)
        # If the indexed array becomes a scalar, return LazilyIndexedArray
        if all(isinstance(ind, integer_types) for ind in indexer.tuple):
            key = BasicIndexer(tuple(k[indexer.tuple] for k in self.key.tuple))
            return LazilyIndexedArray(self.array, key)
        return type(self)(self.array, self._updated_key(indexer))

    def transpose(self, order):
        key = VectorizedIndexer(tuple(k.transpose(order) for k in self.key.tuple))
        return type(self)(self.array, key)

    def __setitem__(self, indexer: ExplicitIndexer, value: Any) -> None:
        raise NotImplementedError(
            "Lazy item assignment with the vectorized indexer is not yet "
            "implemented. Load your data first by .load() or compute()."
        )

    def __repr__(self) -> str:
        return f"{type(self).__name__}(array={self.array!r}, key={self.key!r})"


def _wrap_numpy_scalars(array):
    """Wrap NumPy scalars in 0d arrays."""
    ndim = duck_array_ops.ndim(array)
    if ndim == 0 and (
        isinstance(array, np.generic)
        or not (is_duck_array(array) or isinstance(array, NDArrayMixin))
    ):
        return np.array(array)
    elif hasattr(array, "dtype"):
        return array
    elif ndim == 0:
        return np.array(array)
    else:
        return array


class CopyOnWriteArray(ExplicitlyIndexedNDArrayMixin):
    __slots__ = ("_copied", "array")

    def __init__(self, array: duckarray[Any, Any]):
        self.array = as_indexable(array)
        self._copied = False

    def _ensure_copied(self):
        if not self._copied:
            self.array = as_indexable(np.array(self.array))
            self._copied = True

    def get_duck_array(self):
        return self.array.get_duck_array()

    async def async_get_duck_array(self):
        return await self.array.async_get_duck_array()

    def _oindex_get(self, indexer: OuterIndexer):
        return type(self)(_wrap_numpy_scalars(self.array.oindex[indexer]))

    def _vindex_get(self, indexer: VectorizedIndexer):
        return type(self)(_wrap_numpy_scalars(self.array.vindex[indexer]))

    def __getitem__(self, indexer: ExplicitIndexer):
        self._check_and_raise_if_non_basic_indexer(indexer)
        return type(self)(_wrap_numpy_scalars(self.array[indexer]))

    def transpose(self, order):
        return self.array.transpose(order)

    def _vindex_set(self, indexer: VectorizedIndexer, value: Any) -> None:
        self._ensure_copied()
        self.array.vindex[indexer] = value

    def _oindex_set(self, indexer: OuterIndexer, value: Any) -> None:
        self._ensure_copied()
        self.array.oindex[indexer] = value

    def __setitem__(self, indexer: ExplicitIndexer, value: Any) -> None:
        self._check_and_raise_if_non_basic_indexer(indexer)
        self._ensure_copied()

        self.array[indexer] = value

    def __deepcopy__(self, memo):
        # CopyOnWriteArray is used to wrap backend array objects, which might
        # point to files on disk, so we can't rely on the default deepcopy
        # implementation.
        return type(self)(self.array)


class MemoryCachedArray(ExplicitlyIndexedNDArrayMixin):
    __slots__ = ("array",)

    def __init__(self, array):
        self.array = _wrap_numpy_scalars(as_indexable(array))

    def get_duck_array(self):
        duck_array = self.array.get_duck_array()
        # ensure the array object is cached in-memory
        self.array = as_indexable(duck_array)
        return duck_array

    async def async_get_duck_array(self):
        duck_array = await self.array.async_get_duck_array()
        # ensure the array object is cached in-memory
        self.array = as_indexable(duck_array)
        return duck_array

    def _oindex_get(self, indexer: OuterIndexer):
        return type(self)(_wrap_numpy_scalars(self.array.oindex[indexer]))

    def _vindex_get(self, indexer: VectorizedIndexer):
        return type(self)(_wrap_numpy_scalars(self.array.vindex[indexer]))

    def __getitem__(self, indexer: ExplicitIndexer):
        self._check_and_raise_if_non_basic_indexer(indexer)
        return type(self)(_wrap_numpy_scalars(self.array[indexer]))

    def transpose(self, order):
        return self.array.transpose(order)

    def _vindex_set(self, indexer: VectorizedIndexer, value: Any) -> None:
        self.array.vindex[indexer] = value

    def _oindex_set(self, indexer: OuterIndexer, value: Any) -> None:
        self.array.oindex[indexer] = value

    def __setitem__(self, indexer: ExplicitIndexer, value: Any) -> None:
        self._check_and_raise_if_non_basic_indexer(indexer)
        self.array[indexer] = value


def as_indexable(array):
    """
    This function always returns a ExplicitlyIndexed subclass,
    so that the vectorized indexing is always possible with the returned
    object.
    """
    if isinstance(array, ExplicitlyIndexed):
        return array
    if isinstance(array, np.ndarray):
        return NumpyIndexingAdapter(array)
    if isinstance(array, pd.Index):
        return PandasIndexingAdapter(array)
    if is_duck_dask_array(array):
        return DaskIndexingAdapter(array)
    if hasattr(array, "__array_namespace__"):
        return ArrayApiIndexingAdapter(array)
    if hasattr(array, "__array_function__"):
        return NdArrayLikeIndexingAdapter(array)

    raise TypeError(f"Invalid array type: {type(array)}")


def _outer_to_vectorized_indexer(
    indexer: BasicIndexer | OuterIndexer, shape: _Shape
) -> VectorizedIndexer:
    """Convert an OuterIndexer into an vectorized indexer.

    Parameters
    ----------
    indexer : Outer/Basic Indexer
        An indexer to convert.
    shape : tuple
        Shape of the array subject to the indexing.

    Returns
    -------
    VectorizedIndexer
        Tuple suitable for use to index a NumPy array with vectorized indexing.
        Each element is an array: broadcasting them together gives the shape
        of the result.
    """
    key = indexer.tuple

    n_dim = len([k for k in key if not isinstance(k, integer_types)])
    i_dim = 0
    new_key = []
    for k, size in zip(key, shape, strict=True):
        if isinstance(k, integer_types):
            new_key.append(np.array(k).reshape((1,) * n_dim))
        else:  # np.ndarray or slice
            if isinstance(k, slice):
                k = np.arange(*k.indices(size))
            assert k.dtype.kind in {"i", "u"}
            new_shape = [(1,) * i_dim + (k.size,) + (1,) * (n_dim - i_dim - 1)]
            new_key.append(k.reshape(*new_shape))
            i_dim += 1
    return VectorizedIndexer(tuple(new_key))


def _outer_to_numpy_indexer(indexer: BasicIndexer | OuterIndexer, shape: _Shape):
    """Convert an OuterIndexer into an indexer for NumPy.

    Parameters
    ----------
    indexer : Basic/OuterIndexer
        An indexer to convert.
    shape : tuple
        Shape of the array subject to the indexing.

    Returns
    -------
    tuple
        Tuple suitable for use to index a NumPy array.
    """
    if len([k for k in indexer.tuple if not isinstance(k, slice)]) <= 1:
        # If there is only one vector and all others are slice,
        # it can be safely used in mixed basic/advanced indexing.
        # Boolean index should already be converted to integer array.
        return indexer.tuple
    else:
        return _outer_to_vectorized_indexer(indexer, shape).tuple


def _combine_indexers(old_key, shape: _Shape, new_key) -> VectorizedIndexer:
    """Combine two indexers.

    Parameters
    ----------
    old_key : ExplicitIndexer
        The first indexer for the original array
    shape : tuple of ints
        Shape of the original array to be indexed by old_key
    new_key
        The second indexer for indexing original[old_key]
    """
    if not isinstance(old_key, VectorizedIndexer):
        old_key = _outer_to_vectorized_indexer(old_key, shape)
    if len(old_key.tuple) == 0:
        return new_key

    new_shape = np.broadcast(*old_key.tuple).shape
    if isinstance(new_key, VectorizedIndexer):
        new_key = _arrayize_vectorized_indexer(new_key, new_shape)
    else:
        new_key = _outer_to_vectorized_indexer(new_key, new_shape)

    return VectorizedIndexer(
        tuple(o[new_key.tuple] for o in np.broadcast_arrays(*old_key.tuple))
    )


@enum.unique
class IndexingSupport(enum.Enum):
    # for backends that support only basic indexer
    BASIC = 0
    # for backends that support basic / outer indexer
    OUTER = 1
    # for backends that support outer indexer including at most 1 vector.
    OUTER_1VECTOR = 2
    # for backends that support full vectorized indexer.
    VECTORIZED = 3


def explicit_indexing_adapter(
    key: ExplicitIndexer,
    shape: _Shape,
    indexing_support: IndexingSupport,
    raw_indexing_method: Callable[..., Any],
) -> Any:
    """Support explicit indexing by delegating to a raw indexing method.

    Outer and/or vectorized indexers are supported by indexing a second time
    with a NumPy array.

    Parameters
    ----------
    key : ExplicitIndexer
        Explicit indexing object.
    shape : Tuple[int, ...]
        Shape of the indexed array.
    indexing_support : IndexingSupport enum
        Form of indexing supported by raw_indexing_method.
    raw_indexing_method : callable
        Function (like ndarray.__getitem__) that when called with indexing key
        in the form of a tuple returns an indexed array.

    Returns
    -------
    Indexing result, in the form of a duck numpy-array.
    """
    raw_key, numpy_indices = decompose_indexer(key, shape, indexing_support)
    result = raw_indexing_method(raw_key.tuple)
    if numpy_indices.tuple:
        # index the loaded duck array
        indexable = as_indexable(result)
        result = apply_indexer(indexable, numpy_indices)
    return result


async def async_explicit_indexing_adapter(
    key: ExplicitIndexer,
    shape: _Shape,
    indexing_support: IndexingSupport,
    raw_indexing_method: Callable[..., Any],
) -> Any:
    raw_key, numpy_indices = decompose_indexer(key, shape, indexing_support)
    result = await raw_indexing_method(raw_key.tuple)
    if numpy_indices.tuple:
        # index the loaded duck array
        indexable = as_indexable(result)
        result = apply_indexer(indexable, numpy_indices)
    return result


def apply_indexer(indexable, indexer: ExplicitIndexer):
    """Apply an indexer to an indexable object."""
    if isinstance(indexer, VectorizedIndexer):
        return indexable.vindex[indexer]
    elif isinstance(indexer, OuterIndexer):
        return indexable.oindex[indexer]
    else:
        return indexable[indexer]


def set_with_indexer(indexable, indexer: ExplicitIndexer, value: Any) -> None:
    """Set values in an indexable object using an indexer."""
    if isinstance(indexer, VectorizedIndexer):
        indexable.vindex[indexer] = value
    elif isinstance(indexer, OuterIndexer):
        indexable.oindex[indexer] = value
    else:
        indexable[indexer] = value


def decompose_indexer(
    indexer: ExplicitIndexer, shape: _Shape, indexing_support: IndexingSupport
) -> tuple[ExplicitIndexer, ExplicitIndexer]:
    if isinstance(indexer, VectorizedIndexer):
        return _decompose_vectorized_indexer(indexer, shape, indexing_support)
    if isinstance(indexer, BasicIndexer | OuterIndexer):
        return _decompose_outer_indexer(indexer, shape, indexing_support)
    raise TypeError(f"unexpected key type: {indexer}")


def _decompose_slice(key: slice, size: int) -> tuple[slice, slice]:
    """convert a slice to successive two slices. The first slice always has
    a positive step.

    >>> _decompose_slice(slice(2, 98, 2), 99)
    (slice(2, 98, 2), slice(None, None, None))

    >>> _decompose_slice(slice(98, 2, -2), 99)
    (slice(4, 99, 2), slice(None, None, -1))

    >>> _decompose_slice(slice(98, 2, -2), 98)
    (slice(3, 98, 2), slice(None, None, -1))

    >>> _decompose_slice(slice(360, None, -10), 361)
    (slice(0, 361, 10), slice(None, None, -1))
    """
    start, stop, step = key.indices(size)
    if step > 0:
        # If key already has a positive step, use it as is in the backend
        return key, slice(None)
    else:
        # determine stop precisely for step > 1 case
        # Use the range object to do the calculation
        # e.g. [98:2:-2] -> [98:3:-2]
        exact_stop = range(start, stop, step)[-1]
        return slice(exact_stop, start + 1, -step), slice(None, None, -1)


def _decompose_vectorized_indexer(
    indexer: VectorizedIndexer,
    shape: _Shape,
    indexing_support: IndexingSupport,
) -> tuple[ExplicitIndexer, ExplicitIndexer]:
    """
    Decompose vectorized indexer to the successive two indexers, where the
    first indexer will be used to index backend arrays, while the second one
    is used to index loaded on-memory np.ndarray.

    Parameters
    ----------
    indexer : VectorizedIndexer
    indexing_support : one of IndexerSupport entries

    Returns
    -------
    backend_indexer: OuterIndexer or BasicIndexer
    np_indexers: an ExplicitIndexer (VectorizedIndexer / BasicIndexer)

    Notes
    -----
    This function is used to realize the vectorized indexing for the backend
    arrays that only support basic or outer indexing.

    As an example, let us consider to index a few elements from a backend array
    with a vectorized indexer ([0, 3, 1], [2, 3, 2]).
    Even if the backend array only supports outer indexing, it is more
    efficient to load a subslice of the array than loading the entire array,

    >>> array = np.arange(36).reshape(6, 6)
    >>> backend_indexer = OuterIndexer((np.array([0, 1, 3]), np.array([2, 3])))
    >>> # load subslice of the array
    ... array = NumpyIndexingAdapter(array).oindex[backend_indexer]
    >>> np_indexer = VectorizedIndexer((np.array([0, 2, 1]), np.array([0, 1, 0])))
    >>> # vectorized indexing for on-memory np.ndarray.
    ... NumpyIndexingAdapter(array).vindex[np_indexer]
    array([ 2, 21,  8])
    """
    assert isinstance(indexer, VectorizedIndexer)

    if indexing_support is IndexingSupport.VECTORIZED:
        return indexer, BasicIndexer(())

    backend_indexer_elems = []
    np_indexer_elems = []
    # convert negative indices
    indexer_elems = [
        np.where(k < 0, k + s, k) if isinstance(k, np.ndarray) else k
        for k, s in zip(indexer.tuple, shape, strict=True)
    ]

    for k, s in zip(indexer_elems, shape, strict=True):
        if isinstance(k, slice):
            # If it is a slice, then we will slice it as-is
            # (but make its step positive) in the backend,
            # and then use all of it (slice(None)) for the in-memory portion.
            bk_slice, np_slice = _decompose_slice(k, s)
            backend_indexer_elems.append(bk_slice)
            np_indexer_elems.append(np_slice)
        else:
            # If it is a (multidimensional) np.ndarray, just pickup the used
            # keys without duplication and store them as a 1d-np.ndarray.
            oind, vind = np.unique(k, return_inverse=True)
            backend_indexer_elems.append(oind)
            np_indexer_elems.append(vind.reshape(*k.shape))

    backend_indexer = OuterIndexer(tuple(backend_indexer_elems))
    np_indexer = VectorizedIndexer(tuple(np_indexer_elems))

    if indexing_support is IndexingSupport.OUTER:
        return backend_indexer, np_indexer

    # If the backend does not support outer indexing,
    # backend_indexer (OuterIndexer) is also decomposed.
    backend_indexer1, np_indexer1 = _decompose_outer_indexer(
        backend_indexer, shape, indexing_support
    )
    np_indexer = _combine_indexers(np_indexer1, shape, np_indexer)
    return backend_indexer1, np_indexer


def _decompose_outer_indexer(
    indexer: BasicIndexer | OuterIndexer,
    shape: _Shape,
    indexing_support: IndexingSupport,
) -> tuple[ExplicitIndexer, ExplicitIndexer]:
    """
    Decompose outer indexer to the successive two indexers, where the
    first indexer will be used to index backend arrays, while the second one
    is used to index the loaded on-memory np.ndarray.

    Parameters
    ----------
    indexer : OuterIndexer or BasicIndexer
    indexing_support : One of the entries of IndexingSupport

    Returns
    -------
    backend_indexer: OuterIndexer or BasicIndexer
    np_indexers: an ExplicitIndexer (OuterIndexer / BasicIndexer)

    Notes
    -----
    This function is used to realize the vectorized indexing for the backend
    arrays that only support basic or outer indexing.

    As an example, let us consider to index a few elements from a backend array
    with a orthogonal indexer ([0, 3, 1], [2, 3, 2]).
    Even if the backend array only supports basic indexing, it is more
    efficient to load a subslice of the array than loading the entire array,

    >>> array = np.arange(36).reshape(6, 6)
    >>> backend_indexer = BasicIndexer((slice(0, 3), slice(2, 4)))
    >>> # load subslice of the array
    ... array = NumpyIndexingAdapter(array)[backend_indexer]
    >>> np_indexer = OuterIndexer((np.array([0, 2, 1]), np.array([0, 1, 0])))
    >>> # outer indexing for on-memory np.ndarray.
    ... NumpyIndexingAdapter(array).oindex[np_indexer]
    array([[ 2,  3,  2],
           [14, 15, 14],
           [ 8,  9,  8]])
    """
    backend_indexer: list[Any] = []
    np_indexer: list[Any] = []

    assert isinstance(indexer, OuterIndexer | BasicIndexer)

    if indexing_support == IndexingSupport.VECTORIZED:
        for k, s in zip(indexer.tuple, shape, strict=False):
            if isinstance(k, slice):
                # If it is a slice, then we will slice it as-is
                # (but make its step positive) in the backend,
                bk_slice, np_slice = _decompose_slice(k, s)
                backend_indexer.append(bk_slice)
                np_indexer.append(np_slice)
            else:
                backend_indexer.append(k)
                if not is_scalar(k):
                    np_indexer.append(slice(None))
        return type(indexer)(tuple(backend_indexer)), BasicIndexer(tuple(np_indexer))

    # make indexer positive
    pos_indexer: list[np.ndarray | int | np.number] = []
    for k, s in zip(indexer.tuple, shape, strict=False):
        if isinstance(k, np.ndarray):
            pos_indexer.append(np.where(k < 0, k + s, k))
        elif isinstance(k, integer_types) and k < 0:
            pos_indexer.append(k + s)
        else:
            pos_indexer.append(k)
    indexer_elems = pos_indexer

    if indexing_support is IndexingSupport.OUTER_1VECTOR:
        # some backends such as h5py supports only 1 vector in indexers
        # We choose the most efficient axis
        gains = [
            (
                (np.max(k) - np.min(k) + 1.0) / len(np.unique(k))
                if isinstance(k, np.ndarray)
                else 0
            )
            for k in indexer_elems
        ]
        array_index = np.argmax(np.array(gains)) if len(gains) > 0 else None

        for i, (k, s) in enumerate(zip(indexer_elems, shape, strict=False)):
            if isinstance(k, np.ndarray) and i != array_index:
                # np.ndarray key is converted to slice that covers the entire
                # entries of this key.
                backend_indexer.append(slice(np.min(k), np.max(k) + 1))
                np_indexer.append(k - np.min(k))
            elif isinstance(k, np.ndarray):
                # Remove duplicates and sort them in the increasing order
                pkey, ekey = np.unique(k, return_inverse=True)
                backend_indexer.append(pkey)
                np_indexer.append(ekey)
            elif isinstance(k, integer_types):
                backend_indexer.append(k)
            else:  # slice:  convert positive step slice for backend
                bk_slice, np_slice = _decompose_slice(cast(slice, k), s)
                backend_indexer.append(bk_slice)
                np_indexer.append(np_slice)

        return (OuterIndexer(tuple(backend_indexer)), OuterIndexer(tuple(np_indexer)))

    if indexing_support == IndexingSupport.OUTER:
        for k, s in zip(indexer_elems, shape, strict=False):
            if isinstance(k, slice):
                # slice:  convert positive step slice for backend
                bk_slice, np_slice = _decompose_slice(k, s)
                backend_indexer.append(bk_slice)
                np_indexer.append(np_slice)
            elif isinstance(k, integer_types):
                backend_indexer.append(k)
            elif isinstance(k, np.ndarray) and (np.diff(k) >= 0).all():
                backend_indexer.append(k)
                np_indexer.append(slice(None))
            else:
                # Remove duplicates and sort them in the increasing order
                oind, vind = np.unique(k, return_inverse=True)
                backend_indexer.append(oind)
                np_indexer.append(vind.reshape(*k.shape))

        return (OuterIndexer(tuple(backend_indexer)), OuterIndexer(tuple(np_indexer)))

    # basic indexer
    assert indexing_support == IndexingSupport.BASIC

    for k, s in zip(indexer_elems, shape, strict=False):
        if isinstance(k, np.ndarray):
            # np.ndarray key is converted to slice that covers the entire
            # entries of this key.
            backend_indexer.append(slice(np.min(k), np.max(k) + 1))
            np_indexer.append(k - np.min(k))
        elif isinstance(k, integer_types):
            backend_indexer.append(k)
        else:  # slice:  convert positive step slice for backend
            bk_slice, np_slice = _decompose_slice(cast(slice, k), s)
            backend_indexer.append(bk_slice)
            np_indexer.append(np_slice)

    return (BasicIndexer(tuple(backend_indexer)), OuterIndexer(tuple(np_indexer)))


def _posify_indices(indices: Any, size: int) -> np.ndarray:
    """Convert negative indices by their equivalent positive indices.

    Note: the resulting indices may still be out of bounds (< 0 or >= size).

    """
    return np.where(indices < 0, size + indices, indices)


def _check_bounds(indices: Any, size: int):
    """Check if the given indices are all within the array boundaries."""
    if np.any((indices < 0) | (indices >= size)):
        raise IndexError("out of bounds index")


def _arrayize_outer_indexer(indexer: OuterIndexer, shape) -> OuterIndexer:
    """Return a similar oindex with after replacing slices by arrays and
    negative indices by their corresponding positive indices.

    Also check if array indices are within bounds.

    """
    new_key = []

    for axis, value in enumerate(indexer.tuple):
        size = shape[axis]
        if isinstance(value, slice):
            value = _expand_slice(value, size)
        else:
            value = _posify_indices(value, size)
            _check_bounds(value, size)
        new_key.append(value)

    return OuterIndexer(tuple(new_key))


def _arrayize_vectorized_indexer(
    indexer: VectorizedIndexer, shape: _Shape
) -> VectorizedIndexer:
    """Return an identical vindex but slices are replaced by arrays"""
    slices = [v for v in indexer.tuple if isinstance(v, slice)]
    if len(slices) == 0:
        return indexer

    arrays = [v for v in indexer.tuple if isinstance(v, np.ndarray)]
    n_dim = arrays[0].ndim if len(arrays) > 0 else 0
    i_dim = 0
    new_key = []
    for v, size in zip(indexer.tuple, shape, strict=True):
        if isinstance(v, np.ndarray):
            new_key.append(np.reshape(v, v.shape + (1,) * len(slices)))
        else:  # slice
            shape = (1,) * (n_dim + i_dim) + (-1,) + (1,) * (len(slices) - i_dim - 1)
            new_key.append(np.arange(*v.indices(size)).reshape(shape))
            i_dim += 1
    return VectorizedIndexer(tuple(new_key))


def _chunked_array_with_chunks_hint(
    array, chunks, chunkmanager: ChunkManagerEntrypoint[Any]
):
    """Create a chunked array using the chunks hint for dimensions of size > 1."""

    if len(chunks) < array.ndim:
        raise ValueError("not enough chunks in hint")
    new_chunks = []
    for chunk, size in zip(chunks, array.shape, strict=False):
        new_chunks.append(chunk if size > 1 else (1,))
    return chunkmanager.from_array(array, new_chunks)  # type: ignore[arg-type]


def _logical_any(args):
    return functools.reduce(operator.or_, args)


def _masked_result_drop_slice(key, data: duckarray[Any, Any] | None = None):
    key = (k for k in key if not isinstance(k, slice))
    chunks_hint = getattr(data, "chunks", None)

    new_keys = []
    for k in key:
        if isinstance(k, np.ndarray):
            if is_chunked_array(data):  # type: ignore[arg-type]
                chunkmanager = get_chunked_array_type(data)
                new_keys.append(
                    _chunked_array_with_chunks_hint(k, chunks_hint, chunkmanager)
                )
            elif isinstance(data, array_type("sparse")):
                import sparse

                new_keys.append(sparse.COO.from_numpy(k))
            else:
                new_keys.append(k)
        else:
            new_keys.append(k)

    mask = _logical_any(k == -1 for k in new_keys)
    return mask


def create_mask(
    indexer: ExplicitIndexer, shape: _Shape, data: duckarray[Any, Any] | None = None
):
    """Create a mask for indexing with a fill-value.

    Parameters
    ----------
    indexer : ExplicitIndexer
        Indexer with -1 in integer or ndarray value to indicate locations in
        the result that should be masked.
    shape : tuple
        Shape of the array being indexed.
    data : optional
        Data for which mask is being created. If data is a dask arrays, its chunks
        are used as a hint for chunks on the resulting mask. If data is a sparse
        array, the returned mask is also a sparse array.

    Returns
    -------
    mask : bool, np.ndarray, SparseArray or dask.array.Array with dtype=bool
        Same type as data. Has the same shape as the indexing result.
    """
    if isinstance(indexer, OuterIndexer):
        key = _outer_to_vectorized_indexer(indexer, shape).tuple
        assert not any(isinstance(k, slice) for k in key)
        mask = _masked_result_drop_slice(key, data)

    elif isinstance(indexer, VectorizedIndexer):
        key = indexer.tuple
        base_mask = _masked_result_drop_slice(key, data)
        slice_shape = tuple(
            np.arange(*k.indices(size)).size
            for k, size in zip(key, shape, strict=False)
            if isinstance(k, slice)
        )
        expanded_mask = base_mask[(Ellipsis,) + (np.newaxis,) * len(slice_shape)]
        mask = duck_array_ops.broadcast_to(expanded_mask, base_mask.shape + slice_shape)

    elif isinstance(indexer, BasicIndexer):
        mask = any(k == -1 for k in indexer.tuple)

    else:
        raise TypeError(f"unexpected key type: {type(indexer)}")

    return mask


def _posify_mask_subindexer(
    index: np.ndarray[Any, np.dtype[np.generic]],
) -> np.ndarray[Any, np.dtype[np.generic]]:
    """Convert masked indices in a flat array to the nearest unmasked index.

    Parameters
    ----------
    index : np.ndarray
        One dimensional ndarray with dtype=int.

    Returns
    -------
    np.ndarray
        One dimensional ndarray with all values equal to -1 replaced by an
        adjacent non-masked element.
    """
    masked = index == -1
    unmasked_locs = np.flatnonzero(~masked)
    if not unmasked_locs.size:
        # indexing unmasked_locs is invalid
        return np.zeros_like(index)
    masked_locs = np.flatnonzero(masked)
    prev_value = np.maximum(0, np.searchsorted(unmasked_locs, masked_locs) - 1)
    new_index = index.copy()
    new_index[masked_locs] = index[unmasked_locs[prev_value]]
    return new_index


def posify_mask_indexer(indexer: ExplicitIndexer) -> ExplicitIndexer:
    """Convert masked values (-1) in an indexer to nearest unmasked values.

    This routine is useful for dask, where it can be much faster to index
    adjacent points than arbitrary points from the end of an array.

    Parameters
    ----------
    indexer : ExplicitIndexer
        Input indexer.

    Returns
    -------
    ExplicitIndexer
        Same type of input, with all values in ndarray keys equal to -1
        replaced by an adjacent non-masked element.
    """
    key = tuple(
        (
            _posify_mask_subindexer(k.ravel()).reshape(k.shape)
            if isinstance(k, np.ndarray)
            else k
        )
        for k in indexer.tuple
    )
    return type(indexer)(key)


def is_fancy_indexer(indexer: Any) -> bool:
    """Return False if indexer is a int, slice, a 1-dimensional list, or a 0 or
    1-dimensional ndarray; in all other cases return True
    """
    if isinstance(indexer, int | slice) and not isinstance(indexer, bool):
        return False
    if isinstance(indexer, np.ndarray):
        return indexer.ndim > 1
    if isinstance(indexer, list):
        return bool(indexer) and not isinstance(indexer[0], int)
    return True


class NumpyIndexingAdapter(IndexingAdapter):
    """Wrap a NumPy array to use explicit indexing."""

    __slots__ = ("array",)

    def __init__(self, array):
        # In NumpyIndexingAdapter we only allow to store bare np.ndarray
        if not isinstance(array, np.ndarray):
            raise TypeError(
                "NumpyIndexingAdapter only wraps np.ndarray. "
                f"Trying to wrap {type(array)}"
            )
        self.array = array

    def transpose(self, order):
        return self.array.transpose(order)

    def _oindex_get(self, indexer: OuterIndexer):
        key = _outer_to_numpy_indexer(indexer, self.array.shape)
        return self.array[key]

    def _vindex_get(self, indexer: VectorizedIndexer):
        _assert_not_chunked_indexer(indexer.tuple)
        array = NumpyVIndexAdapter(self.array)
        return array[indexer.tuple]

    def __getitem__(self, indexer: ExplicitIndexer):
        self._check_and_raise_if_non_basic_indexer(indexer)

        array = self.array
        # We want 0d slices rather than scalars. This is achieved by
        # appending an ellipsis (see
        # https://numpy.org/doc/stable/reference/arrays.indexing.html#detailed-notes).
        key = indexer.tuple + (Ellipsis,)
        return array[key]

    def _safe_setitem(self, array, key: tuple[Any, ...], value: Any) -> None:
        try:
            array[key] = value
        except ValueError as exc:
            # More informative exception if read-only view
            if not array.flags.writeable and not array.flags.owndata:
                raise ValueError(
                    "Assignment destination is a view.  "
                    "Do you want to .copy() array first?"
                ) from exc
            else:
                raise exc

    def _oindex_set(self, indexer: OuterIndexer, value: Any) -> None:
        key = _outer_to_numpy_indexer(indexer, self.array.shape)
        self._safe_setitem(self.array, key, value)

    def _vindex_set(self, indexer: VectorizedIndexer, value: Any) -> None:
        array = NumpyVIndexAdapter(self.array)
        self._safe_setitem(array, indexer.tuple, value)

    def __setitem__(self, indexer: ExplicitIndexer, value: Any) -> None:
        self._check_and_raise_if_non_basic_indexer(indexer)
        array = self.array
        # We want 0d slices rather than scalars. This is achieved by
        # appending an ellipsis (see
        # https://numpy.org/doc/stable/reference/arrays.indexing.html#detailed-notes).
        key = indexer.tuple + (Ellipsis,)
        self._safe_setitem(array, key, value)


class NdArrayLikeIndexingAdapter(NumpyIndexingAdapter):
    __slots__ = ("array",)

    def __init__(self, array):
        if not hasattr(array, "__array_function__"):
            raise TypeError(
                "NdArrayLikeIndexingAdapter must wrap an object that "
                "implements the __array_function__ protocol"
            )
        self.array = array


class ArrayApiIndexingAdapter(IndexingAdapter):
    """Wrap an array API array to use explicit indexing."""

    __slots__ = ("array",)

    def __init__(self, array):
        if not hasattr(array, "__array_namespace__"):
            raise TypeError(
                "ArrayApiIndexingAdapter must wrap an object that "
                "implements the __array_namespace__ protocol"
            )
        self.array = array

    def _oindex_get(self, indexer: OuterIndexer):
        # manual orthogonal indexing (implemented like DaskIndexingAdapter)
        key = indexer.tuple
        value = self.array
        for axis, subkey in reversed(list(enumerate(key))):
            value = value[(slice(None),) * axis + (subkey, Ellipsis)]
        return value

    def _vindex_get(self, indexer: VectorizedIndexer):
        raise TypeError("Vectorized indexing is not supported")

    def __getitem__(self, indexer: ExplicitIndexer):
        self._check_and_raise_if_non_basic_indexer(indexer)
        return self.array[indexer.tuple]

    def _oindex_set(self, indexer: OuterIndexer, value: Any) -> None:
        self.array[indexer.tuple] = value

    def _vindex_set(self, indexer: VectorizedIndexer, value: Any) -> None:
        raise TypeError("Vectorized indexing is not supported")

    def __setitem__(self, indexer: ExplicitIndexer, value: Any) -> None:
        self._check_and_raise_if_non_basic_indexer(indexer)
        self.array[indexer.tuple] = value

    def transpose(self, order):
        xp = self.array.__array_namespace__()
        return xp.permute_dims(self.array, order)


def _apply_vectorized_indexer_dask_wrapper(indices, coord):
    from xarray.core.indexing import VectorizedIndexer, apply_indexer, as_indexable

    return apply_indexer(
        as_indexable(coord), VectorizedIndexer((indices.squeeze(axis=-1),))
    )


def _assert_not_chunked_indexer(idxr: tuple[Any, ...]) -> None:
    if any(is_chunked_array(i) for i in idxr):
        raise ValueError(
            "Cannot index with a chunked array indexer. "
            "Please chunk the array you are indexing first, "
            "and drop any indexed dimension coordinate variables. "
            "Alternatively, call `.compute()` on any chunked arrays in the indexer."
        )


class DaskIndexingAdapter(IndexingAdapter):
    """Wrap a dask array to support explicit indexing."""

    __slots__ = ("array",)

    def __init__(self, array):
        """This adapter is created in Variable.__getitem__ in
        Variable._broadcast_indexes.
        """
        self.array = array

    def _oindex_get(self, indexer: OuterIndexer):
        key = indexer.tuple
        try:
            return self.array[key]
        except NotImplementedError:
            # manual orthogonal indexing
            value = self.array
            for axis, subkey in reversed(list(enumerate(key))):
                value = value[(slice(None),) * axis + (subkey,)]
            return value

    def _vindex_get(self, indexer: VectorizedIndexer):
        try:
            return self.array.vindex[indexer.tuple]
        except IndexError as e:
            # TODO: upstream to dask
            has_dask = any(is_duck_dask_array(i) for i in indexer.tuple)
            # this only works for "small" 1d coordinate arrays with one chunk
            # it is intended for idxmin, idxmax, and allows indexing with
            # the nD array output of argmin, argmax
            if (
                not has_dask
                or len(indexer.tuple) > 1
                or math.prod(self.array.numblocks) > 1
                or self.array.ndim > 1
            ):
                raise e
            (idxr,) = indexer.tuple
            if idxr.ndim == 0:
                return self.array[idxr.data]
            else:
                import dask.array

                return dask.array.map_blocks(
                    _apply_vectorized_indexer_dask_wrapper,
                    idxr[..., np.newaxis],
                    self.array,
                    chunks=idxr.chunks,
                    drop_axis=-1,
                    dtype=self.array.dtype,
                )

    def __getitem__(self, indexer: ExplicitIndexer):
        self._check_and_raise_if_non_basic_indexer(indexer)
        return self.array[indexer.tuple]

    def _oindex_set(self, indexer: OuterIndexer, value: Any) -> None:
        num_non_slices = sum(0 if isinstance(k, slice) else 1 for k in indexer.tuple)
        if num_non_slices > 1:
            raise NotImplementedError(
                "xarray can't set arrays with multiple array indices to dask yet."
            )
        self.array[indexer.tuple] = value

    def _vindex_set(self, indexer: VectorizedIndexer, value: Any) -> None:
        self.array.vindex[indexer.tuple] = value

    def __setitem__(self, indexer: ExplicitIndexer, value: Any) -> None:
        self._check_and_raise_if_non_basic_indexer(indexer)
        self.array[indexer.tuple] = value

    def transpose(self, order):
        return self.array.transpose(order)


class PandasIndexingAdapter(IndexingAdapter):
    """Wrap a pandas.Index to preserve dtypes and handle explicit indexing."""

    __slots__ = ("_dtype", "array")

    array: pd.Index
    _dtype: np.dtype | pd.api.extensions.ExtensionDtype

    def __init__(
        self,
        array: pd.Index,
        dtype: DTypeLike | pd.api.extensions.ExtensionDtype | None = None,
    ):
        from xarray.core.indexes import safe_cast_to_index

        self.array = safe_cast_to_index(array)

        if dtype is None:
            if is_allowed_extension_array(array):
                cast(pd.api.extensions.ExtensionDtype, array.dtype)
                self._dtype = array.dtype
            else:
                self._dtype = get_valid_numpy_dtype(array)
        elif is_allowed_extension_array_dtype(dtype):
            self._dtype = cast(pd.api.extensions.ExtensionDtype, dtype)
        else:
            self._dtype = np.dtype(cast(DTypeLike, dtype))

    @property
    def _in_memory(self) -> bool:
        # prevent costly conversion of a memory-saving pd.RangeIndex into a
        # large numpy array.
        return not isinstance(self.array, pd.RangeIndex)

    @property
    def dtype(self) -> np.dtype | pd.api.extensions.ExtensionDtype:  # type: ignore[override]
        return self._dtype

    def _get_numpy_dtype(self, dtype: np.typing.DTypeLike | None = None) -> np.dtype:
        if dtype is None:
            if is_valid_numpy_dtype(self.dtype):
                return cast(np.dtype, self.dtype)
            else:
                return get_valid_numpy_dtype(self.array)
        else:
            return np.dtype(dtype)

    def __array__(
        self,
        dtype: np.typing.DTypeLike | None = None,
        /,
        *,
        copy: bool | None = None,
    ) -> np.ndarray:
        dtype = self._get_numpy_dtype(dtype)
        array = self.array

        if isinstance(array, pd.PeriodIndex):
            with suppress(AttributeError):
                # this might not be public API
                array = array.astype("object")

        if Version(np.__version__) >= Version("2.0.0"):
            return np.asarray(array.values, dtype=dtype, copy=copy)
        else:
            return np.asarray(array.values, dtype=dtype)

    def get_duck_array(self) -> np.ndarray | PandasExtensionArray:
        # We return an PandasExtensionArray wrapper type that satisfies
        # duck array protocols.
        # `NumpyExtensionArray` is excluded
        if is_allowed_extension_array(self.array):
            from xarray.core.extension_array import PandasExtensionArray

            return PandasExtensionArray(self.array.array)
        return np.asarray(self)

    @property
    def shape(self) -> _Shape:
        return (len(self.array),)

    def _convert_scalar(self, item) -> np.ndarray:
        if item is pd.NaT:
            # work around the impossibility of casting NaT with asarray
            # note: it probably would be better in general to return
            # pd.Timestamp rather np.than datetime64 but this is easier
            # (for now)
            item = np.datetime64("NaT", "ns")
        elif isinstance(item, pd.Timedelta):
            item = item.to_numpy()
        elif isinstance(item, timedelta):
            item = np.timedelta64(item)
        elif isinstance(item, pd.Timestamp):
            # Work around for GH: pydata/xarray#1932 and numpy/numpy#10668
            # numpy fails to convert pd.Timestamp to np.datetime64[ns]
            item = np.asarray(item.to_datetime64())
        elif self.dtype != object:
            dtype = self._get_numpy_dtype()
            item = np.asarray(item, dtype=dtype)

        # as for numpy.ndarray indexing, we always want the result to be
        # a NumPy array.
        return to_0d_array(item)

    def _index_get(
        self, indexer: ExplicitIndexer, func_name: str
    ) -> PandasIndexingAdapter | np.ndarray:
        key = indexer.tuple

        if len(key) == 1:
            # unpack key so it can index a pandas.Index object (pandas.Index
            # objects don't like tuples)
            (key,) = key

        # if multidimensional key, convert the index to numpy array and index the latter
        if getattr(key, "ndim", 0) > 1:
            indexable = NumpyIndexingAdapter(np.asarray(self))
            return getattr(indexable, func_name)(indexer)

        # otherwise index the pandas index then re-wrap or convert the result
        result = self.array[key]

        if isinstance(result, pd.Index):
            return type(self)(result, dtype=self.dtype)
        else:
            return self._convert_scalar(result)

    def _oindex_get(self, indexer: OuterIndexer) -> PandasIndexingAdapter | np.ndarray:
        return self._index_get(indexer, "_oindex_get")

    def _vindex_get(
        self, indexer: VectorizedIndexer
    ) -> PandasIndexingAdapter | np.ndarray:
        _assert_not_chunked_indexer(indexer.tuple)
        return self._index_get(indexer, "_vindex_get")

    def __getitem__(
        self, indexer: ExplicitIndexer
    ) -> PandasIndexingAdapter | np.ndarray:
        return self._index_get(indexer, "__getitem__")

    def transpose(self, order) -> pd.Index:
        return self.array  # self.array should be always one-dimensional

    def _repr_inline_(self, max_width: int) -> str:
        # we want to display values in the inline repr for lazy coordinates too
        # (pd.RangeIndex and pd.MultiIndex). `format_array_flat` prevents loading
        # the whole array in memory.
        from xarray.core.formatting import format_array_flat

        return format_array_flat(self, max_width)

    def __repr__(self) -> str:
        return f"{type(self).__name__}(array={self.array!r}, dtype={self.dtype!r})"

    def copy(self, deep: bool = True) -> Self:
        # Not the same as just writing `self.array.copy(deep=deep)`, as
        # shallow copies of the underlying numpy.ndarrays become deep ones
        # upon pickling
        # >>> len(pickle.dumps((self.array, self.array)))
        # 4000281
        # >>> len(pickle.dumps((self.array, self.array.copy(deep=False))))
        # 8000341
        array = self.array.copy(deep=True) if deep else self.array
        return type(self)(array, self._dtype)

    @property
    def nbytes(self) -> int:
        if is_allowed_extension_array(self.array):
            return self.array.nbytes

        dtype = self._get_numpy_dtype()
        return dtype.itemsize * len(self.array)


class PandasMultiIndexingAdapter(PandasIndexingAdapter):
    """Handles explicit indexing for a pandas.MultiIndex.

    This allows creating one instance for each multi-index level while
    preserving indexing efficiency (memoized + might reuse another instance with
    the same multi-index).
    """

    __slots__ = ("_dtype", "adapter", "array", "level")

    array: pd.MultiIndex
    _dtype: np.dtype | pd.api.extensions.ExtensionDtype
    level: str | None

    def __init__(
        self,
        array: pd.MultiIndex,
        dtype: DTypeLike | pd.api.extensions.ExtensionDtype | None = None,
        level: str | None = None,
    ):
        super().__init__(array, dtype)
        self.level = level

    def __array__(
        self,
        dtype: DTypeLike | None = None,
        /,
        *,
        copy: bool | None = None,
    ) -> np.ndarray:
        dtype = self._get_numpy_dtype(dtype)

        if self.level is not None:
            return np.asarray(
                self.array.get_level_values(self.level).values, dtype=dtype
            )
        else:
            return super().__array__(dtype, copy=copy)

    @property
    def _in_memory(self) -> bool:
        # The pd.MultiIndex's data is fully in memory, but it has a different
        # layout than the level and dimension coordinate arrays. Marking this
        # adapter class as a "lazy" array will prevent costly conversion when,
        # e.g., formatting the Xarray reprs.
        return False

    def _convert_scalar(self, item: Any):
        if isinstance(item, tuple) and self.level is not None:
            idx = tuple(self.array.names).index(self.level)
            item = item[idx]
        return super()._convert_scalar(item)

    def _index_get(
        self, indexer: ExplicitIndexer, func_name: str
    ) -> PandasIndexingAdapter | np.ndarray:
        result = super()._index_get(indexer, func_name)
        if isinstance(result, type(self)):
            result.level = self.level
        return result

    def __repr__(self) -> str:
        if self.level is None:
            return super().__repr__()
        else:
            props = (
                f"(array={self.array!r}, level={self.level!r}, dtype={self.dtype!r})"
            )
            return f"{type(self).__name__}{props}"

    def _repr_inline_(self, max_width: int) -> str:
        if self.level is None:
            return "MultiIndex"
        else:
            return super()._repr_inline_(max_width=max_width)

    def copy(self, deep: bool = True) -> Self:
        # see PandasIndexingAdapter.copy
        array = self.array.copy(deep=True) if deep else self.array
        return type(self)(array, self._dtype, self.level)


class CoordinateTransformIndexingAdapter(IndexingAdapter):
    """Wrap a CoordinateTransform as a lazy coordinate array.

    Supports explicit indexing (both outer and vectorized).

    """

    _transform: CoordinateTransform
    _coord_name: Hashable
    _dims: tuple[str, ...]

    def __init__(
        self,
        transform: CoordinateTransform,
        coord_name: Hashable,
        dims: tuple[str, ...] | None = None,
    ):
        self._transform = transform
        self._coord_name = coord_name
        self._dims = dims or transform.dims

    @property
    def dtype(self) -> np.dtype:
        return self._transform.dtype

    @property
    def shape(self) -> tuple[int, ...]:
        return tuple(self._transform.dim_size.values())

    @property
    def _in_memory(self) -> bool:
        return False

    def get_duck_array(self) -> np.ndarray:
        all_coords = self._transform.generate_coords(dims=self._dims)
        return np.asarray(all_coords[self._coord_name])

    def _oindex_get(self, indexer: OuterIndexer):
        expanded_indexer_ = OuterIndexer(expanded_indexer(indexer.tuple, self.ndim))
        array_indexer = _arrayize_outer_indexer(expanded_indexer_, self.shape)

        positions = np.meshgrid(*array_indexer.tuple, indexing="ij")
        dim_positions = dict(zip(self._dims, positions, strict=False))

        result = self._transform.forward(dim_positions)
        return np.asarray(result[self._coord_name]).squeeze()

    def _oindex_set(self, indexer: OuterIndexer, value: Any) -> None:
        raise TypeError(
            "setting values is not supported on coordinate transform arrays."
        )

    def _vindex_get(self, indexer: VectorizedIndexer):
        expanded_indexer_ = VectorizedIndexer(
            expanded_indexer(indexer.tuple, self.ndim)
        )
        array_indexer = _arrayize_vectorized_indexer(expanded_indexer_, self.shape)

        dim_positions = {}
        for i, (dim, pos) in enumerate(
            zip(self._dims, array_indexer.tuple, strict=False)
        ):
            pos = _posify_indices(pos, self.shape[i])
            _check_bounds(pos, self.shape[i])
            dim_positions[dim] = pos

        result = self._transform.forward(dim_positions)
        return np.asarray(result[self._coord_name])

    def _vindex_set(self, indexer: VectorizedIndexer, value: Any) -> None:
        raise TypeError(
            "setting values is not supported on coordinate transform arrays."
        )

    def __getitem__(self, indexer: ExplicitIndexer):
        # TODO: make it lazy (i.e., re-calculate and re-wrap the transform) when possible?
        self._check_and_raise_if_non_basic_indexer(indexer)

        # also works with basic indexing
        return self._oindex_get(OuterIndexer(indexer.tuple))

    def __setitem__(self, indexer: ExplicitIndexer, value: Any) -> None:
        raise TypeError(
            "setting values is not supported on coordinate transform arrays."
        )

    def transpose(self, order: Iterable[int]) -> Self:
        new_dims = tuple(self._dims[i] for i in order)
        return type(self)(self._transform, self._coord_name, new_dims)

    def __repr__(self: Any) -> str:
        return f"{type(self).__name__}(transform={self._transform!r})"

    def _repr_inline_(self, max_width: int) -> str:
        # we want to display values in the inline repr for this lazy coordinate
        # `format_array_flat` prevents loading the whole array in memory.
        from xarray.core.formatting import format_array_flat

        return format_array_flat(self, max_width)
