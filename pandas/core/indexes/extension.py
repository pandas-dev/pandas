"""
Shared methods for Index subclasses backed by ExtensionArray.
"""
from typing import List

import numpy as np

from pandas.util._decorators import Appender, cache_readonly

from pandas.core.dtypes.common import is_dtype_equal

from pandas.core.arrays import ExtensionArray

from .base import Index, _index_shared_docs


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


class ExtensionIndex(Index):
    _data: ExtensionArray

    def __getitem__(self, key):
        result = self._data[key]
        if isinstance(result, type(self._data)):
            return type(self)(result, name=self.name)

        # Includes cases where we get a 2D ndarray back for MPL compat
        return result

    def __iter__(self):
        return self._data.__iter__()

    @property
    def _ndarray_values(self) -> np.ndarray:
        return self._data._ndarray_values

    @Appender(_index_shared_docs["astype"])
    def astype(self, dtype, copy=True):
        if is_dtype_equal(self.dtype, dtype) and copy is False:
            # Ensure that self.astype(self.dtype) is self
            return self

        new_values = self._data.astype(dtype, copy=copy)

        # pass copy=False because any copying will be done in the
        #  _data.astype call above
        return Index(new_values, dtype=new_values.dtype, name=self.name, copy=False)

    def dropna(self, how="any"):
        if how not in ("any", "all"):
            raise ValueError(f"invalid how option: {how}")

        if self.hasnans:
            return self._shallow_copy(self._data[~self._isnan])
        return self._shallow_copy()

    def _get_unique_index(self, dropna=False):
        if self.is_unique and not dropna:
            return self

        result = self._data.unique()
        if dropna and self.hasnans:
            result = result[~result.isna()]
        return self._shallow_copy(result)

    def unique(self, level=None):
        if level is not None:
            self._validate_index_level(level)

        result = self._data.unique()
        return self._shallow_copy(result)
