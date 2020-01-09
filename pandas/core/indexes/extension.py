"""
Shared methods for Index subclasses backed by ExtensionArray.
"""
from typing import List

from pandas.compat.numpy import function as nv
from pandas.util._decorators import cache_readonly

from pandas.core.dtypes.common import ensure_platform_int, is_dtype_equal
from pandas.core.dtypes.generic import ABCSeries

from pandas.core.arrays import ExtensionArray
from pandas.core.indexes.base import Index
from pandas.core.ops import get_op_result_name


def inherit_from_data(name: str, delegate, cache: bool = False):
    """
    Make an alias for a method of the underlying ExtensionArray.

    Parameters
    ----------
    name : str
        Name of an attribute the class should inherit from its EA parent.
    delegate : class
    cache : bool, default False
        Whether to convert wrapped properties into cache_readonly

    Returns
    -------
    attribute, method, property, or cache_readonly
    """

    attr = getattr(delegate, name)

    if isinstance(attr, property):
        if cache:
            method = cache_readonly(attr.fget)

        else:

            def fget(self):
                return getattr(self._data, name)

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
            return result

        method.__name__ = name
        method.__doc__ = attr.__doc__
    return method


def inherit_names(names: List[str], delegate, cache: bool = False):
    """
    Class decorator to pin attributes from an ExtensionArray to a Index subclass.

    Parameters
    ----------
    names : List[str]
    delegate : class
    cache : bool, default False
    """

    def wrapper(cls):
        for name in names:
            meth = inherit_from_data(name, delegate, cache=cache)
            setattr(cls, name, meth)

        return cls

    return wrapper


def _make_wrapped_comparison_op(opname):
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


def make_wrapped_arith_op(opname):
    def method(self, other):
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

    _data: ExtensionArray

    __eq__ = _make_wrapped_comparison_op("__eq__")
    __ne__ = _make_wrapped_comparison_op("__ne__")
    __lt__ = _make_wrapped_comparison_op("__lt__")
    __gt__ = _make_wrapped_comparison_op("__gt__")
    __le__ = _make_wrapped_comparison_op("__le__")
    __ge__ = _make_wrapped_comparison_op("__ge__")

    def repeat(self, repeats, axis=None):
        nv.validate_repeat(tuple(), dict(axis=axis))
        result = self._data.repeat(repeats, axis=axis)
        return self._shallow_copy(result)

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

    def astype(self, dtype, copy=True):
        if is_dtype_equal(self.dtype, dtype) and copy is False:
            # Ensure that self.astype(self.dtype) is self
            return self

        new_values = self._data.astype(dtype, copy=copy)

        # pass copy=False because any copying will be done in the
        #  _data.astype call above
        return Index(new_values, dtype=new_values.dtype, name=self.name, copy=False)
