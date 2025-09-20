"""Internal utilities; not for external use"""

# Some functions in this module are derived from functions in pandas. For
# reference, here is a copy of the pandas copyright notice:

# BSD 3-Clause License

# Copyright (c) 2008-2011, AQR Capital Management, LLC, Lambda Foundry, Inc. and PyData Development Team
# All rights reserved.

# Copyright (c) 2011-2022, Open source contributors.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.

# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.

# * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
from __future__ import annotations

import contextlib
import difflib
import functools
import importlib
import inspect
import io
import itertools
import math
import os
import re
import sys
import warnings
from collections.abc import (
    Callable,
    Collection,
    Container,
    Hashable,
    ItemsView,
    Iterable,
    Iterator,
    KeysView,
    Mapping,
    MutableMapping,
    MutableSet,
    Sequence,
    ValuesView,
)
from collections.abc import (
    Set as AbstractSet,
)
from enum import Enum
from pathlib import Path
from types import EllipsisType, ModuleType
from typing import (
    TYPE_CHECKING,
    Any,
    Generic,
    Literal,
    TypeGuard,
    TypeVar,
    cast,
    overload,
)

import numpy as np
import pandas as pd

from xarray.namedarray.utils import (  # noqa: F401
    ReprObject,
    drop_missing_dims,
    either_dict_or_kwargs,
    infix_dims,
    is_dask_collection,
    is_dict_like,
    is_duck_array,
    is_duck_dask_array,
    module_available,
    to_0d_object_array,
)

if TYPE_CHECKING:
    from xarray.core.types import Dims, ErrorOptionsWithWarn

K = TypeVar("K")
V = TypeVar("V")
T = TypeVar("T")


def is_allowed_extension_array_dtype(dtype: Any):
    return pd.api.types.is_extension_array_dtype(dtype) and not isinstance(  # noqa: TID251
        dtype, pd.StringDtype
    )


def is_allowed_extension_array(array: Any) -> bool:
    return (
        hasattr(array, "dtype")
        and is_allowed_extension_array_dtype(array.dtype)
        and not isinstance(array, pd.arrays.NumpyExtensionArray)  # type: ignore[attr-defined]
    )


def alias_message(old_name: str, new_name: str) -> str:
    return f"{old_name} has been deprecated. Use {new_name} instead."


def alias_warning(old_name: str, new_name: str, stacklevel: int = 3) -> None:
    warnings.warn(
        alias_message(old_name, new_name), FutureWarning, stacklevel=stacklevel
    )


def alias(obj: Callable[..., T], old_name: str) -> Callable[..., T]:
    assert isinstance(old_name, str)

    @functools.wraps(obj)
    def wrapper(*args, **kwargs):
        alias_warning(old_name, obj.__name__)
        return obj(*args, **kwargs)

    wrapper.__doc__ = alias_message(old_name, obj.__name__)
    return wrapper


def did_you_mean(
    word: Hashable, possibilities: Iterable[Hashable], *, n: int = 10
) -> str:
    """
    Suggest a few correct words based on a list of possibilities

    Parameters
    ----------
    word : Hashable
        Word to compare to a list of possibilities.
    possibilities : Iterable of Hashable
        The iterable of Hashable that contains the correct values.
    n : int, default: 10
        Maximum number of suggestions to show.

    Examples
    --------
    >>> did_you_mean("bluch", ("blech", "gray_r", 1, None, (2, 56)))
    "Did you mean one of ('blech',)?"
    >>> did_you_mean("none", ("blech", "gray_r", 1, None, (2, 56)))
    'Did you mean one of (None,)?'

    See also
    --------
    https://en.wikipedia.org/wiki/String_metric
    """
    # Convert all values to string, get_close_matches doesn't handle all hashables:
    possibilities_str: dict[str, Hashable] = {str(k): k for k in possibilities}

    msg = ""
    if len(
        best_str := difflib.get_close_matches(
            str(word), list(possibilities_str.keys()), n=n
        )
    ):
        best = tuple(possibilities_str[k] for k in best_str)
        msg = f"Did you mean one of {best}?"

    return msg


def get_valid_numpy_dtype(array: np.ndarray | pd.Index) -> np.dtype:
    """Return a numpy compatible dtype from either
    a numpy array or a pandas.Index.

    Used for wrapping a pandas.Index as an xarray.Variable.

    """
    if isinstance(array, pd.PeriodIndex):
        return np.dtype("O")

    if hasattr(array, "categories"):
        # category isn't a real numpy dtype
        dtype = array.categories.dtype
        if not is_valid_numpy_dtype(dtype):
            dtype = np.dtype("O")
        return dtype

    if not is_valid_numpy_dtype(array.dtype):
        return np.dtype("O")

    return array.dtype  # type: ignore[return-value]


def maybe_coerce_to_str(index, original_coords):
    """maybe coerce a pandas Index back to a nunpy array of type str

    pd.Index uses object-dtype to store str - try to avoid this for coords
    """
    from xarray.core import dtypes

    try:
        result_type = dtypes.result_type(*original_coords)
    except TypeError:
        pass
    else:
        if result_type.kind in "SU":
            index = np.asarray(index, dtype=result_type.type)

    return index


def maybe_wrap_array(original, new_array):
    """Wrap a transformed array with __array_wrap__ if it can be done safely.

    This lets us treat arbitrary functions that take and return ndarray objects
    like ufuncs, as long as they return an array with the same shape.
    """
    # in case func lost array's metadata
    if isinstance(new_array, np.ndarray) and new_array.shape == original.shape:
        return original.__array_wrap__(new_array)
    else:
        return new_array


def equivalent(first: T, second: T) -> bool:
    """Compare two objects for equivalence (identity or equality), using
    array_equiv if either object is an ndarray. If both objects are lists,
    equivalent is sequentially called on all the elements.
    """
    # TODO: refactor to avoid circular import
    from xarray.core import duck_array_ops

    if first is second:
        return True
    if isinstance(first, np.ndarray) or isinstance(second, np.ndarray):
        return duck_array_ops.array_equiv(first, second)
    if isinstance(first, list) or isinstance(second, list):
        return list_equiv(first, second)  # type: ignore[arg-type]
    return (first == second) or (pd.isnull(first) and pd.isnull(second))  # type: ignore[call-overload]


def list_equiv(first: Sequence[T], second: Sequence[T]) -> bool:
    if len(first) != len(second):
        return False
    return all(itertools.starmap(equivalent, zip(first, second, strict=True)))


def peek_at(iterable: Iterable[T]) -> tuple[T, Iterator[T]]:
    """Returns the first value from iterable, as well as a new iterator with
    the same content as the original iterable
    """
    gen = iter(iterable)
    peek = next(gen)
    return peek, itertools.chain([peek], gen)


def update_safety_check(
    first_dict: Mapping[K, V],
    second_dict: Mapping[K, V],
    compat: Callable[[V, V], bool] = equivalent,
) -> None:
    """Check the safety of updating one dictionary with another.

    Raises ValueError if dictionaries have non-compatible values for any key,
    where compatibility is determined by identity (they are the same item) or
    the `compat` function.

    Parameters
    ----------
    first_dict, second_dict : dict-like
        All items in the second dictionary are checked against for conflicts
        against items in the first dictionary.
    compat : function, optional
        Binary operator to determine if two values are compatible. By default,
        checks for equivalence.
    """
    for k, v in second_dict.items():
        if k in first_dict and not compat(v, first_dict[k]):
            raise ValueError(
                "unsafe to merge dictionaries without "
                f"overriding values; conflicting key {k!r}"
            )


def remove_incompatible_items(
    first_dict: MutableMapping[K, V],
    second_dict: Mapping[K, V],
    compat: Callable[[V, V], bool] = equivalent,
) -> None:
    """Remove incompatible items from the first dictionary in-place.

    Items are retained if their keys are found in both dictionaries and the
    values are compatible.

    Parameters
    ----------
    first_dict, second_dict : dict-like
        Mappings to merge.
    compat : function, optional
        Binary operator to determine if two values are compatible. By default,
        checks for equivalence.
    """
    for k in list(first_dict):
        if k not in second_dict or not compat(first_dict[k], second_dict[k]):
            del first_dict[k]


def is_full_slice(value: Any) -> bool:
    return isinstance(value, slice) and value == slice(None)


def is_list_like(value: Any) -> TypeGuard[list | tuple]:
    return isinstance(value, list | tuple)


def _is_scalar(value, include_0d):
    from xarray.core.variable import NON_NUMPY_SUPPORTED_ARRAY_TYPES

    if include_0d:
        include_0d = getattr(value, "ndim", None) == 0
    return (
        include_0d
        or isinstance(value, str | bytes)
        or not (
            isinstance(value, (Iterable,) + NON_NUMPY_SUPPORTED_ARRAY_TYPES)
            or hasattr(value, "__array_function__")
            or hasattr(value, "__array_namespace__")
        )
    )


def is_scalar(value: Any, include_0d: bool = True) -> TypeGuard[Hashable]:
    """Whether to treat a value as a scalar.

    Any non-iterable, string, or 0-D array
    """
    return _is_scalar(value, include_0d)


def is_valid_numpy_dtype(dtype: Any) -> bool:
    try:
        np.dtype(dtype)
    except (TypeError, ValueError):
        return False
    else:
        return True


def to_0d_array(value: Any) -> np.ndarray:
    """Given a value, wrap it in a 0-D numpy.ndarray."""
    if np.isscalar(value) or (isinstance(value, np.ndarray) and value.ndim == 0):
        return np.array(value)
    else:
        return to_0d_object_array(value)


def dict_equiv(
    first: Mapping[K, V],
    second: Mapping[K, V],
    compat: Callable[[V, V], bool] = equivalent,
) -> bool:
    """Test equivalence of two dict-like objects. If any of the values are
    numpy arrays, compare them correctly.

    Parameters
    ----------
    first, second : dict-like
        Dictionaries to compare for equality
    compat : function, optional
        Binary operator to determine if two values are compatible. By default,
        checks for equivalence.

    Returns
    -------
    equals : bool
        True if the dictionaries are equal
    """
    for k in first:
        if k not in second or not compat(first[k], second[k]):
            return False
    return all(k in first for k in second)


def compat_dict_intersection(
    first_dict: Mapping[K, V],
    second_dict: Mapping[K, V],
    compat: Callable[[V, V], bool] = equivalent,
) -> MutableMapping[K, V]:
    """Return the intersection of two dictionaries as a new dictionary.

    Items are retained if their keys are found in both dictionaries and the
    values are compatible.

    Parameters
    ----------
    first_dict, second_dict : dict-like
        Mappings to merge.
    compat : function, optional
        Binary operator to determine if two values are compatible. By default,
        checks for equivalence.

    Returns
    -------
    intersection : dict
        Intersection of the contents.
    """
    new_dict = dict(first_dict)
    remove_incompatible_items(new_dict, second_dict, compat)
    return new_dict


def compat_dict_union(
    first_dict: Mapping[K, V],
    second_dict: Mapping[K, V],
    compat: Callable[[V, V], bool] = equivalent,
) -> MutableMapping[K, V]:
    """Return the union of two dictionaries as a new dictionary.

    An exception is raised if any keys are found in both dictionaries and the
    values are not compatible.

    Parameters
    ----------
    first_dict, second_dict : dict-like
        Mappings to merge.
    compat : function, optional
        Binary operator to determine if two values are compatible. By default,
        checks for equivalence.

    Returns
    -------
    union : dict
        union of the contents.
    """
    new_dict = dict(first_dict)
    update_safety_check(first_dict, second_dict, compat)
    new_dict.update(second_dict)
    return new_dict


class Frozen(Mapping[K, V]):
    """Wrapper around an object implementing the mapping interface to make it
    immutable. If you really want to modify the mapping, the mutable version is
    saved under the `mapping` attribute.
    """

    __slots__ = ("mapping",)

    def __init__(self, mapping: Mapping[K, V]):
        self.mapping = mapping

    def __getitem__(self, key: K) -> V:
        return self.mapping[key]

    def __iter__(self) -> Iterator[K]:
        return iter(self.mapping)

    def __len__(self) -> int:
        return len(self.mapping)

    def __contains__(self, key: object) -> bool:
        return key in self.mapping

    def __repr__(self) -> str:
        return f"{type(self).__name__}({self.mapping!r})"


def FrozenDict(*args, **kwargs) -> Frozen:
    return Frozen(dict(*args, **kwargs))


class FrozenMappingWarningOnValuesAccess(Frozen[K, V]):
    """
    Class which behaves like a Mapping but warns if the values are accessed.

    Temporary object to aid in deprecation cycle of `Dataset.dims` (see GH issue #8496).
    `Dataset.dims` is being changed from returning a mapping of dimension names to lengths to just
    returning a frozen set of dimension names (to increase consistency with `DataArray.dims`).
    This class retains backwards compatibility but raises a warning only if the return value
    of ds.dims is used like a dictionary (i.e. it doesn't raise a warning if used in a way that
    would also be valid for a FrozenSet, e.g. iteration).
    """

    __slots__ = ("mapping",)

    def _warn(self) -> None:
        emit_user_level_warning(
            "The return type of `Dataset.dims` will be changed to return a set of dimension names in future, "
            "in order to be more consistent with `DataArray.dims`. To access a mapping from dimension names to lengths, "
            "please use `Dataset.sizes`.",
            FutureWarning,
        )

    def __getitem__(self, key: K) -> V:
        self._warn()
        return super().__getitem__(key)

    @overload
    def get(self, key: K, /) -> V | None: ...

    @overload
    def get(self, key: K, /, default: V | T) -> V | T: ...

    def get(self, key: K, default: T | None = None) -> V | T | None:
        self._warn()
        return super().get(key, default)

    def keys(self) -> KeysView[K]:
        self._warn()
        return super().keys()

    def items(self) -> ItemsView[K, V]:
        self._warn()
        return super().items()

    def values(self) -> ValuesView[V]:
        self._warn()
        return super().values()


class FilteredMapping(Mapping[K, V]):
    """Implements the Mapping interface. Uses the wrapped mapping for item lookup
    and a separate wrapped keys collection for iteration.

    Can be used to construct a mapping object from another dict-like object without
    eagerly accessing its items or when a mapping object is expected but only
    iteration over keys is actually used.

    Note: keys should be a subset of mapping, but FilteredMapping does not
    validate consistency of the provided `keys` and `mapping`. It is the
    caller's responsibility to ensure that they are suitable for the task at
    hand.
    """

    __slots__ = ("keys_", "mapping")

    def __init__(self, keys: Collection[K], mapping: Mapping[K, V]):
        self.keys_ = keys  # .keys is already a property on Mapping
        self.mapping = mapping

    def __getitem__(self, key: K) -> V:
        if key not in self.keys_:
            raise KeyError(key)
        return self.mapping[key]

    def __iter__(self) -> Iterator[K]:
        return iter(self.keys_)

    def __len__(self) -> int:
        return len(self.keys_)

    def __repr__(self) -> str:
        return f"{type(self).__name__}(keys={self.keys_!r}, mapping={self.mapping!r})"


class OrderedSet(MutableSet[T]):
    """A simple ordered set.

    The API matches the builtin set, but it preserves insertion order of elements, like
    a dict. Note that, unlike in an OrderedDict, equality tests are not order-sensitive.
    """

    _d: dict[T, None]

    __slots__ = ("_d",)

    def __init__(self, values: Iterable[T] | None = None):
        self._d = {}
        if values is not None:
            self.update(values)

    # Required methods for MutableSet

    def __contains__(self, value: Hashable) -> bool:
        return value in self._d

    def __iter__(self) -> Iterator[T]:
        return iter(self._d)

    def __len__(self) -> int:
        return len(self._d)

    def add(self, value: T) -> None:
        self._d[value] = None

    def discard(self, value: T) -> None:
        del self._d[value]

    # Additional methods

    def update(self, values: Iterable[T]) -> None:
        self._d.update(dict.fromkeys(values))

    def __repr__(self) -> str:
        return f"{type(self).__name__}({list(self)!r})"


class NdimSizeLenMixin:
    """Mixin class that extends a class that defines a ``shape`` property to
    one that also defines ``ndim``, ``size`` and ``__len__``.
    """

    __slots__ = ()

    @property
    def ndim(self: Any) -> int:
        """
        Number of array dimensions.

        See Also
        --------
        numpy.ndarray.ndim
        """
        return len(self.shape)

    @property
    def size(self: Any) -> int:
        """
        Number of elements in the array.

        Equal to ``np.prod(a.shape)``, i.e., the product of the array’s dimensions.

        See Also
        --------
        numpy.ndarray.size
        """
        return math.prod(self.shape)

    def __len__(self: Any) -> int:
        try:
            return self.shape[0]
        except IndexError as err:
            raise TypeError("len() of unsized object") from err


class NDArrayMixin(NdimSizeLenMixin):
    """Mixin class for making wrappers of N-dimensional arrays that conform to
    the ndarray interface required for the data argument to Variable objects.

    A subclass should set the `array` property and override one or more of
    `dtype`, `shape` and `__getitem__`.
    """

    __slots__ = ()

    @property
    def dtype(self: Any) -> np.dtype:
        return self.array.dtype

    @property
    def shape(self: Any) -> tuple[int, ...]:
        return self.array.shape

    def __getitem__(self: Any, key):
        return self.array[key]

    def __repr__(self: Any) -> str:
        return f"{type(self).__name__}(array={self.array!r})"


@contextlib.contextmanager
def close_on_error(f):
    """Context manager to ensure that a file opened by xarray is closed if an
    exception is raised before the user sees the file object.
    """
    try:
        yield
    except Exception:
        f.close()
        raise


def is_remote_uri(path: str) -> bool:
    """Finds URLs of the form protocol:// or protocol::

    This also matches for http[s]://, which were the only remote URLs
    supported in <=v0.16.2.
    """
    return bool(re.search(r"^[a-z][a-z0-9]*(\://|\:\:)", path))


def read_magic_number_from_file(filename_or_obj, count=8) -> bytes:
    # check byte header to determine file type
    if not isinstance(filename_or_obj, io.IOBase):
        raise TypeError(f"cannot read the magic number from {type(filename_or_obj)}")
    if filename_or_obj.tell() != 0:
        filename_or_obj.seek(0)
    magic_number = filename_or_obj.read(count)
    filename_or_obj.seek(0)
    return magic_number


def try_read_magic_number_from_path(pathlike, count=8) -> bytes | None:
    if isinstance(pathlike, str) or hasattr(pathlike, "__fspath__"):
        path = os.fspath(pathlike)
        try:
            with open(path, "rb") as f:
                return read_magic_number_from_file(f, count)
        except (FileNotFoundError, IsADirectoryError, TypeError):
            pass
    return None


def try_read_magic_number_from_file_or_path(filename_or_obj, count=8) -> bytes | None:
    magic_number = try_read_magic_number_from_path(filename_or_obj, count)
    if magic_number is None:
        with contextlib.suppress(TypeError):
            magic_number = read_magic_number_from_file(filename_or_obj, count)
    return magic_number


def is_uniform_spaced(arr, **kwargs) -> bool:
    """Return True if values of an array are uniformly spaced and sorted.

    >>> is_uniform_spaced(range(5))
    True
    >>> is_uniform_spaced([-4, 0, 100])
    False

    kwargs are additional arguments to ``np.isclose``
    """
    arr = np.array(arr, dtype=float)
    diffs = np.diff(arr)
    return bool(np.isclose(diffs.min(), diffs.max(), **kwargs))


def hashable(v: Any) -> TypeGuard[Hashable]:
    """Determine whether `v` can be hashed."""
    try:
        hash(v)
    except TypeError:
        return False
    return True


def iterable(v: Any) -> TypeGuard[Iterable[Any]]:
    """Determine whether `v` is iterable."""
    try:
        iter(v)
    except TypeError:
        return False
    return True


def iterable_of_hashable(v: Any) -> TypeGuard[Iterable[Hashable]]:
    """Determine whether `v` is an Iterable of Hashables."""
    try:
        it = iter(v)
    except TypeError:
        return False
    return all(hashable(elm) for elm in it)


def decode_numpy_dict_values(attrs: Mapping[K, V]) -> dict[K, V]:
    """Convert attribute values from numpy objects to native Python objects,
    for use in to_dict
    """
    attrs = dict(attrs)
    for k, v in attrs.items():
        if isinstance(v, np.ndarray):
            attrs[k] = cast(V, v.tolist())
        elif isinstance(v, np.generic):
            attrs[k] = v.item()
    return attrs


def ensure_us_time_resolution(val):
    """Convert val out of numpy time, for use in to_dict.
    Needed because of numpy bug GH#7619"""
    if np.issubdtype(val.dtype, np.datetime64):
        val = val.astype("datetime64[us]")
    elif np.issubdtype(val.dtype, np.timedelta64):
        val = val.astype("timedelta64[us]")
    return val


class HiddenKeyDict(MutableMapping[K, V]):
    """Acts like a normal dictionary, but hides certain keys."""

    __slots__ = ("_data", "_hidden_keys")

    # ``__init__`` method required to create instance from class.

    def __init__(self, data: MutableMapping[K, V], hidden_keys: Iterable[K]):
        self._data = data
        self._hidden_keys = frozenset(hidden_keys)

    def _raise_if_hidden(self, key: K) -> None:
        if key in self._hidden_keys:
            raise KeyError(f"Key `{key!r}` is hidden.")

    # The next five methods are requirements of the ABC.
    def __setitem__(self, key: K, value: V) -> None:
        self._raise_if_hidden(key)
        self._data[key] = value

    def __getitem__(self, key: K) -> V:
        self._raise_if_hidden(key)
        return self._data[key]

    def __delitem__(self, key: K) -> None:
        self._raise_if_hidden(key)
        del self._data[key]

    def __iter__(self) -> Iterator[K]:
        for k in self._data:
            if k not in self._hidden_keys:
                yield k

    def __len__(self) -> int:
        num_hidden = len(self._hidden_keys & self._data.keys())
        return len(self._data) - num_hidden


def get_temp_dimname(dims: Container[Hashable], new_dim: Hashable) -> Hashable:
    """Get an new dimension name based on new_dim, that is not used in dims.
    If the same name exists, we add an underscore(s) in the head.

    Example1:
        dims: ['a', 'b', 'c']
        new_dim: ['_rolling']
        -> ['_rolling']
    Example2:
        dims: ['a', 'b', 'c', '_rolling']
        new_dim: ['_rolling']
        -> ['__rolling']
    """
    while new_dim in dims:
        new_dim = "_" + str(new_dim)
    return new_dim


def drop_dims_from_indexers(
    indexers: Mapping[Any, Any],
    dims: Iterable[Hashable] | Mapping[Any, int],
    missing_dims: ErrorOptionsWithWarn,
) -> Mapping[Hashable, Any]:
    """Depending on the setting of missing_dims, drop any dimensions from indexers that
    are not present in dims.

    Parameters
    ----------
    indexers : dict
    dims : sequence
    missing_dims : {"raise", "warn", "ignore"}
    """

    if missing_dims == "raise":
        invalid = indexers.keys() - set(dims)
        if invalid:
            raise ValueError(
                f"Dimensions {invalid} do not exist. Expected one or more of {dims}"
            )

        return indexers

    elif missing_dims == "warn":
        # don't modify input
        indexers = dict(indexers)

        invalid = indexers.keys() - set(dims)
        if invalid:
            warnings.warn(
                f"Dimensions {invalid} do not exist. Expected one or more of {dims}",
                stacklevel=2,
            )
        for key in invalid:
            indexers.pop(key)

        return indexers

    elif missing_dims == "ignore":
        return {key: val for key, val in indexers.items() if key in dims}

    else:
        raise ValueError(
            f"Unrecognised option {missing_dims} for missing_dims argument"
        )


@overload
def parse_dims_as_tuple(
    dim: Dims,
    all_dims: tuple[Hashable, ...],
    *,
    check_exists: bool = True,
    replace_none: Literal[True] = True,
) -> tuple[Hashable, ...]: ...


@overload
def parse_dims_as_tuple(
    dim: Dims,
    all_dims: tuple[Hashable, ...],
    *,
    check_exists: bool = True,
    replace_none: Literal[False],
) -> tuple[Hashable, ...] | EllipsisType | None: ...


def parse_dims_as_tuple(
    dim: Dims,
    all_dims: tuple[Hashable, ...],
    *,
    check_exists: bool = True,
    replace_none: bool = True,
) -> tuple[Hashable, ...] | EllipsisType | None:
    """Parse one or more dimensions.

    A single dimension must be always a str, multiple dimensions
    can be Hashables. This supports e.g. using a tuple as a dimension.
    If you supply e.g. a set of dimensions the order cannot be
    conserved, but for sequences it will be.

    Parameters
    ----------
    dim : str, Iterable of Hashable, "..." or None
        Dimension(s) to parse.
    all_dims : tuple of Hashable
        All possible dimensions.
    check_exists: bool, default: True
        if True, check if dim is a subset of all_dims.
    replace_none : bool, default: True
        If True, return all_dims if dim is None or "...".

    Returns
    -------
    parsed_dims : tuple of Hashable
        Input dimensions as a tuple.
    """
    if dim is None or dim is ...:
        if replace_none:
            return all_dims
        return dim
    if isinstance(dim, str):
        dim = (dim,)
    if check_exists:
        _check_dims(set(dim), set(all_dims))
    return tuple(dim)


@overload
def parse_dims_as_set(
    dim: Dims,
    all_dims: set[Hashable],
    *,
    check_exists: bool = True,
    replace_none: Literal[True] = True,
) -> set[Hashable]: ...


@overload
def parse_dims_as_set(
    dim: Dims,
    all_dims: set[Hashable],
    *,
    check_exists: bool = True,
    replace_none: Literal[False],
) -> set[Hashable] | EllipsisType | None: ...


def parse_dims_as_set(
    dim: Dims,
    all_dims: set[Hashable],
    *,
    check_exists: bool = True,
    replace_none: bool = True,
) -> set[Hashable] | EllipsisType | None:
    """Like parse_dims_as_tuple, but returning a set instead of a tuple."""
    # TODO: Consider removing parse_dims_as_tuple?
    if dim is None or dim is ...:
        if replace_none:
            return all_dims
        return dim
    if isinstance(dim, str):
        dim = {dim}
    dim = set(dim)
    if check_exists:
        _check_dims(dim, all_dims)
    return dim


@overload
def parse_ordered_dims(
    dim: Dims,
    all_dims: tuple[Hashable, ...],
    *,
    check_exists: bool = True,
    replace_none: Literal[True] = True,
) -> tuple[Hashable, ...]: ...


@overload
def parse_ordered_dims(
    dim: Dims,
    all_dims: tuple[Hashable, ...],
    *,
    check_exists: bool = True,
    replace_none: Literal[False],
) -> tuple[Hashable, ...] | EllipsisType | None: ...


def parse_ordered_dims(
    dim: Dims,
    all_dims: tuple[Hashable, ...],
    *,
    check_exists: bool = True,
    replace_none: bool = True,
) -> tuple[Hashable, ...] | EllipsisType | None:
    """Parse one or more dimensions.

    A single dimension must be always a str, multiple dimensions
    can be Hashables. This supports e.g. using a tuple as a dimension.
    An ellipsis ("...") in a sequence of dimensions will be
    replaced with all remaining dimensions. This only makes sense when
    the input is a sequence and not e.g. a set.

    Parameters
    ----------
    dim : str, Sequence of Hashable or "...", "..." or None
        Dimension(s) to parse. If "..." appears in a Sequence
        it always gets replaced with all remaining dims
    all_dims : tuple of Hashable
        All possible dimensions.
    check_exists: bool, default: True
        if True, check if dim is a subset of all_dims.
    replace_none : bool, default: True
        If True, return all_dims if dim is None.

    Returns
    -------
    parsed_dims : tuple of Hashable
        Input dimensions as a tuple.
    """
    if dim is not None and dim is not ... and not isinstance(dim, str) and ... in dim:
        dims_set: set[Hashable | EllipsisType] = set(dim)
        all_dims_set = set(all_dims)
        if check_exists:
            _check_dims(dims_set, all_dims_set)
        if len(all_dims_set) != len(all_dims):
            raise ValueError("Cannot use ellipsis with repeated dims")
        dims = tuple(dim)
        if dims.count(...) > 1:
            raise ValueError("More than one ellipsis supplied")
        other_dims = tuple(d for d in all_dims if d not in dims_set)
        idx = dims.index(...)
        return dims[:idx] + other_dims + dims[idx + 1 :]
    else:
        # mypy cannot resolve that the sequence cannot contain "..."
        return parse_dims_as_tuple(  # type: ignore[call-overload]
            dim=dim,
            all_dims=all_dims,
            check_exists=check_exists,
            replace_none=replace_none,
        )


def _check_dims(dim: AbstractSet[Hashable], all_dims: AbstractSet[Hashable]) -> None:
    wrong_dims = (dim - all_dims) - {...}
    if wrong_dims:
        wrong_dims_str = ", ".join(f"'{d}'" for d in wrong_dims)
        raise ValueError(
            f"Dimension(s) {wrong_dims_str} do not exist. Expected one or more of {all_dims}"
        )


_Accessor = TypeVar("_Accessor")


class UncachedAccessor(Generic[_Accessor]):
    """Acts like a property, but on both classes and class instances

    This class is necessary because some tools (e.g. pydoc and sphinx)
    inspect classes for which property returns itself and not the
    accessor.
    """

    def __init__(self, accessor: type[_Accessor]) -> None:
        self._accessor = accessor

    @overload
    def __get__(self, obj: None, cls) -> type[_Accessor]: ...

    @overload
    def __get__(self, obj: object, cls) -> _Accessor: ...

    def __get__(self, obj: object | None, cls) -> type[_Accessor] | _Accessor:
        if obj is None:
            return self._accessor

        return self._accessor(obj)  # type: ignore[call-arg]  # assume it is a valid accessor!


# Singleton type, as per https://github.com/python/typing/pull/240
class Default(Enum):
    token = 0


_default = Default.token


def iterate_nested(nested_list):
    for item in nested_list:
        if isinstance(item, list):
            yield from iterate_nested(item)
        else:
            yield item


def contains_only_chunked_or_numpy(obj) -> bool:
    """Returns True if xarray object contains only numpy arrays or chunked arrays (i.e. pure dask or cubed).

    Expects obj to be Dataset or DataArray"""
    from xarray.core.dataarray import DataArray
    from xarray.core.indexing import ExplicitlyIndexed
    from xarray.namedarray.pycompat import is_chunked_array

    if isinstance(obj, DataArray):
        obj = obj._to_temp_dataset()

    return all(
        isinstance(var._data, ExplicitlyIndexed | np.ndarray)
        or is_chunked_array(var._data)
        for var in obj._variables.values()
    )


def find_stack_level(test_mode=False) -> int:
    """Find the first place in the stack that is not inside xarray or the Python standard library.

    This is unless the code emanates from a test, in which case we would prefer
    to see the xarray source.

    This function is taken from pandas and modified to exclude standard library paths.

    Parameters
    ----------
    test_mode : bool
        Flag used for testing purposes to switch off the detection of test
        directories in the stack trace.

    Returns
    -------
    stacklevel : int
        First level in the stack that is not part of xarray or the Python standard library.
    """
    import xarray as xr

    pkg_dir = Path(xr.__file__).parent
    test_dir = pkg_dir / "tests"

    std_lib_init = sys.modules["os"].__file__
    # Mostly to appease mypy; I don't think this can happen...
    if std_lib_init is None:
        return 0

    std_lib_dir = Path(std_lib_init).parent

    frame = inspect.currentframe()
    n = 0
    while frame:
        fname = inspect.getfile(frame)
        if (
            fname.startswith(str(pkg_dir))
            and (not fname.startswith(str(test_dir)) or test_mode)
        ) or (
            fname.startswith(str(std_lib_dir))
            and "site-packages" not in fname
            and "dist-packages" not in fname
        ):
            frame = frame.f_back
            n += 1
        else:
            break
    return n


def emit_user_level_warning(message, category=None) -> None:
    """Emit a warning at the user level by inspecting the stack trace."""
    stacklevel = find_stack_level()
    return warnings.warn(message, category=category, stacklevel=stacklevel)


def consolidate_dask_from_array_kwargs(
    from_array_kwargs: dict[Any, Any],
    name: str | None = None,
    lock: bool | None = None,
    inline_array: bool | None = None,
) -> dict[Any, Any]:
    """
    Merge dask-specific kwargs with arbitrary from_array_kwargs dict.

    Temporary function, to be deleted once explicitly passing dask-specific kwargs to .chunk() is deprecated.
    """

    from_array_kwargs = _resolve_doubly_passed_kwarg(
        from_array_kwargs,
        kwarg_name="name",
        passed_kwarg_value=name,
        default=None,
        err_msg_dict_name="from_array_kwargs",
    )
    from_array_kwargs = _resolve_doubly_passed_kwarg(
        from_array_kwargs,
        kwarg_name="lock",
        passed_kwarg_value=lock,
        default=False,
        err_msg_dict_name="from_array_kwargs",
    )
    from_array_kwargs = _resolve_doubly_passed_kwarg(
        from_array_kwargs,
        kwarg_name="inline_array",
        passed_kwarg_value=inline_array,
        default=False,
        err_msg_dict_name="from_array_kwargs",
    )

    return from_array_kwargs


def _resolve_doubly_passed_kwarg(
    kwargs_dict: dict[Any, Any],
    kwarg_name: str,
    passed_kwarg_value: str | bool | None,
    default: bool | None,
    err_msg_dict_name: str,
) -> dict[Any, Any]:
    # if in kwargs_dict but not passed explicitly then just pass kwargs_dict through unaltered
    if kwarg_name in kwargs_dict and passed_kwarg_value is None:
        pass
    # if passed explicitly but not in kwargs_dict then use that
    elif kwarg_name not in kwargs_dict and passed_kwarg_value is not None:
        kwargs_dict[kwarg_name] = passed_kwarg_value
    # if in neither then use default
    elif kwarg_name not in kwargs_dict and passed_kwarg_value is None:
        kwargs_dict[kwarg_name] = default
    # if in both then raise
    else:
        raise ValueError(
            f"argument {kwarg_name} cannot be passed both as a keyword argument and within "
            f"the {err_msg_dict_name} dictionary"
        )

    return kwargs_dict


def attempt_import(module: str) -> ModuleType:
    """Import an optional dependency, and raise an informative error on failure.

    Parameters
    ----------
    module : str
        Module to import. For example, ``'zarr'`` or ``'matplotlib.pyplot'``.

    Returns
    -------
    module : ModuleType
        The Imported module.

    Raises
    ------
    ImportError
        If the module could not be imported.

    Notes
    -----
    Static type checkers will not be able to infer the type of the returned module,
    so it is recommended to precede this function with a direct import of the module,
    guarded by an ``if TYPE_CHECKING`` block, to preserve type checker functionality.
    See the examples section below for a demonstration.

    Examples
    --------
    >>> from xarray.core.utils import attempt_import
    >>> if TYPE_CHECKING:
    ...     import zarr
    ... else:
    ...     zarr = attempt_import("zarr")
    ...
    """
    install_mapping = dict(nc_time_axis="nc-time-axis")
    package_purpose = dict(
        zarr="for working with Zarr stores",
        cftime="for working with non-standard calendars",
        matplotlib="for plotting",
        hypothesis="for the `xarray.testing.strategies` submodule",
    )
    package_name = module.split(".", maxsplit=1)[0]  # e.g. "zarr" from "zarr.storage"
    install_name = install_mapping.get(package_name, package_name)
    reason = package_purpose.get(package_name, "")
    try:
        return importlib.import_module(module)
    except ImportError as e:
        raise ImportError(
            f"The {install_name} package is required {reason}"
            " but could not be imported."
            " Please install it with your package manager (e.g. conda or pip)."
        ) from e


_DEFAULT_NAME = ReprObject("<default-name>")


def result_name(objects: Iterable[Any]) -> Any:
    # use the same naming heuristics as pandas:
    # https://github.com/blaze/blaze/issues/458#issuecomment-51936356
    names = {getattr(obj, "name", _DEFAULT_NAME) for obj in objects}
    names.discard(_DEFAULT_NAME)
    if len(names) == 1:
        (name,) = names
    else:
        name = None
    return name


def _get_func_args(func, param_names):
    """Use `inspect.signature` to try accessing `func` args. Otherwise, ensure
    they are provided by user.
    """
    try:
        func_args = inspect.signature(func).parameters
    except ValueError as err:
        func_args = {}
        if not param_names:
            raise ValueError(
                "Unable to inspect `func` signature, and `param_names` was not provided."
            ) from err
    if param_names:
        params = param_names
    else:
        params = list(func_args)[1:]
        if any(
            (p.kind in [p.VAR_POSITIONAL, p.VAR_KEYWORD]) for p in func_args.values()
        ):
            raise ValueError(
                "`param_names` must be provided because `func` takes variable length arguments."
            )
    return params, func_args
