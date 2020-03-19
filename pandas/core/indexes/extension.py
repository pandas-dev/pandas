"""
Shared methods for Index subclasses backed by ExtensionArray.
"""
from typing import List

import numpy as np

from pandas.compat.numpy import function as nv
from pandas.errors import AbstractMethodError
from pandas.util._decorators import Appender, cache_readonly

from pandas.core.dtypes.common import (
    ensure_platform_int,
    is_dtype_equal,
    is_object_dtype,
)
from pandas.core.dtypes.generic import ABCSeries

from pandas.core.arrays import ExtensionArray
from pandas.core.indexers import deprecate_ndim_indexing
from pandas.core.indexes.base import Index
from pandas.core.ops import get_op_result_name


def inherit_from_data(name: str, delegate, cache: bool = False, wrap: bool = False):
    """
    Make an alias for a method of the underlying ExtensionArray.

    Parameters
    ----------
    name : str
        Name of an attribute the class should inherit from its EA parent.
    delegate : class
    cache : bool, default False
        Whether to convert wrapped properties into cache_readonly
    wrap : bool, default False
        Whether to wrap the inherited result in an Index.

    Returns
    -------
    attribute, method, property, or cache_readonly
    """
    attr = getattr(delegate, name)

    if isinstance(attr, property):
        if cache:

            def cached(self):
                return getattr(self._data, name)

            cached.__name__ = name
            cached.__doc__ = attr.__doc__
            method = cache_readonly(cached)

        else:

            def fget(self):
                result = getattr(self._data, name)
                if wrap:
                    if isinstance(result, type(self._data)):
                        return type(self)._simple_new(result, name=self.name)
                    return Index(result, name=self.name)
                return result

            def fset(self, value):
                setattr(self._data, name, value)

            fget.__name__ = name
            fget.__doc__ = attr.__doc__

            method = property(fget, fset)

    elif not callable(attr):
        # just a normal attribute, no wrapping
        method = attr

    else:

        def method(self, *args, **kwargs):
            result = attr(self._data, *args, **kwargs)
            if wrap:
                if isinstance(result, type(self._data)):
                    return type(self)._simple_new(result, name=self.name)
                return Index(result, name=self.name)
            return result

        method.__name__ = name
        method.__doc__ = attr.__doc__
    return method


def inherit_names(names: List[str], delegate, cache: bool = False, wrap: bool = False):
    """
    Class decorator to pin attributes from an ExtensionArray to a Index subclass.

    Parameters
    ----------
    names : List[str]
    delegate : class
    cache : bool, default False
    wrap : bool, default False
        Whether to wrap the inherited result in an Index.
    """

    def wrapper(cls):
        for name in names:
            meth = inherit_from_data(name, delegate, cache=cache, wrap=wrap)
            setattr(cls, name, meth)

        return cls

    return wrapper


def _make_wrapped_comparison_op(opname: str):
    """
    Create a comparison method that dispatches to ``._data``.
    """

    def wrapper(self, other):
        if isinstance(other, ABCSeries):
            # the arrays defer to Series for comparison ops but the indexes
            #  don't, so we have to unwrap here.
            other = other._values

        other = _maybe_unwrap_index(other)

        op = getattr(self._data, opname)
        return op(other)

    wrapper.__name__ = opname
    return wrapper


def make_wrapped_arith_op(opname: str):
    def method(self, other):
        if (
            isinstance(other, Index)
            and is_object_dtype(other.dtype)
            and type(other) is not Index
        ):
            # We return NotImplemented for object-dtype index *subclasses* so they have
            # a chance to implement ops before we unwrap them.
            # See https://github.com/pandas-dev/pandas/issues/31109
            return NotImplemented
        meth = getattr(self._data, opname)
        result = meth(_maybe_unwrap_index(other))
        return _wrap_arithmetic_op(self, other, result)

    method.__name__ = opname
    return method


def _wrap_arithmetic_op(self, other, result):
    if result is NotImplemented:
        return NotImplemented

    if isinstance(result, tuple):
        # divmod, rdivmod
        assert len(result) == 2
        return (
            _wrap_arithmetic_op(self, other, result[0]),
            _wrap_arithmetic_op(self, other, result[1]),
        )

    if not isinstance(result, Index):
        # Index.__new__ will choose appropriate subclass for dtype
        result = Index(result)

    res_name = get_op_result_name(self, other)
    result.name = res_name
    return result


def _maybe_unwrap_index(obj):
    """
    If operating against another Index object, we need to unwrap the underlying
    data before deferring to the DatetimeArray/TimedeltaArray/PeriodArray
    implementation, otherwise we will incorrectly return NotImplemented.

    Parameters
    ----------
    obj : object

    Returns
    -------
    unwrapped object
    """
    if isinstance(obj, Index):
        return obj._data
    return obj


class ExtensionIndex(Index):
    """
    Index subclass for indexes backed by ExtensionArray.
    """

    # The base class already passes through to _data:
    #  size, __len__, dtype

    _data: ExtensionArray

    __eq__ = _make_wrapped_comparison_op("__eq__")
    __ne__ = _make_wrapped_comparison_op("__ne__")
    __lt__ = _make_wrapped_comparison_op("__lt__")
    __gt__ = _make_wrapped_comparison_op("__gt__")
    __le__ = _make_wrapped_comparison_op("__le__")
    __ge__ = _make_wrapped_comparison_op("__ge__")

    # ---------------------------------------------------------------------
    # NDarray-Like Methods

    def __getitem__(self, key):
        result = self._data[key]
        if isinstance(result, type(self._data)):
            return type(self)(result, name=self.name)

        # Includes cases where we get a 2D ndarray back for MPL compat
        deprecate_ndim_indexing(result)
        return result

    def __iter__(self):
        return self._data.__iter__()

    # ---------------------------------------------------------------------

    def __array__(self, dtype=None) -> np.ndarray:
        return np.asarray(self._data, dtype=dtype)

    def _get_engine_target(self) -> np.ndarray:
        return self._data._values_for_argsort()

    @Appender(Index.dropna.__doc__)
    def dropna(self, how="any"):
        if how not in ("any", "all"):
            raise ValueError(f"invalid how option: {how}")

        if self.hasnans:
            return self._shallow_copy(self._data[~self._isnan])
        return self._shallow_copy()

    def repeat(self, repeats, axis=None):
        nv.validate_repeat(tuple(), dict(axis=axis))
        result = self._data.repeat(repeats, axis=axis)
        return self._shallow_copy(result)

    def insert(self, loc: int, item):
        # ExtensionIndex subclasses must override Index.insert
        raise AbstractMethodError(self)

    def _concat_same_dtype(self, to_concat, name):
        arr = type(self._data)._concat_same_type(to_concat)
        return type(self)._simple_new(arr, name=name)

    @Appender(Index.take.__doc__)
    def take(self, indices, axis=0, allow_fill=True, fill_value=None, **kwargs):
        nv.validate_take(tuple(), kwargs)
        indices = ensure_platform_int(indices)

        taken = self._assert_take_fillable(
            self._data,
            indices,
            allow_fill=allow_fill,
            fill_value=fill_value,
            na_value=self._na_value,
        )
        return type(self)(taken, name=self.name)

    def unique(self, level=None):
        if level is not None:
            self._validate_index_level(level)

        result = self._data.unique()
        return self._shallow_copy(result)

    def _get_unique_index(self, dropna=False):
        if self.is_unique and not dropna:
            return self

        result = self._data.unique()
        if dropna and self.hasnans:
            result = result[~result.isna()]
        return self._shallow_copy(result)

    @Appender(Index.map.__doc__)
    def map(self, mapper, na_action=None):
        # Try to run function on index first, and then on elements of index
        # Especially important for group-by functionality
        try:
            result = mapper(self)

            # Try to use this result if we can
            if isinstance(result, np.ndarray):
                result = Index(result)

            if not isinstance(result, Index):
                raise TypeError("The map function must return an Index object")
            return result
        except Exception:
            return self.astype(object).map(mapper)

    @Appender(Index.astype.__doc__)
    def astype(self, dtype, copy=True):
        if is_dtype_equal(self.dtype, dtype) and copy is False:
            # Ensure that self.astype(self.dtype) is self
            return self

        new_values = self._data.astype(dtype, copy=copy)

        # pass copy=False because any copying will be done in the
        #  _data.astype call above
        return Index(new_values, dtype=new_values.dtype, name=self.name, copy=False)
