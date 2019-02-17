"""A boolean mask interface.

This module provides an interface to a numpy / pyarrow boolean mask.
This is limited as not all of the implementations can hold NA, so
for consistency this is an internal.
"""

import copy

import numpy as np

from pandas.api.extensions import ExtensionDtype
from pandas.api.types import is_scalar
from pandas.core.arrays.base import ExtensionArray
from pandas.core.missing import isna


class MaskDtype(ExtensionDtype):

    type = np.bool_
    kind = 'b'
    name = 'bool'

    @classmethod
    def construct_from_string(cls, string):
        if string == cls.name:
            return cls()
        else:
            raise TypeError("Cannot construct a '{}' from "
                            "'{}'".format(cls, string))

    def _is_boolean(self):
        return True

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        # compare == to np.dtype('bool')
        if isinstance(other, str):
            return other == self.name
        elif isinstance(other, type(self)):
            return True
        elif isinstance(other, np.dtype):
            return other == 'bool'
        else:
            return hash(self) == hash(other)


class MaskArray(ExtensionArray):
    """Common baseclass for both pyarrow and numpy masked arrays"""
    _typ = "maskarray"

    @classmethod
    def _from_sequence(cls, scalars, dtype=None, copy=False):
        return cls.from_scalars(scalars)

    @property
    def size(self):
        return len(self)

    def __eq__(self, other):
        return np.array(self, copy=False) == np.array(other, copy=False)

    def __len__(self):
        return len(self._data)

    def isna(self):
        nas = isna(np.array(self._data, copy=False))
        return type(self).from_scalars(nas)

    def __invert__(self):
        return type(self).from_scalars(
            ~np.array(self._data, copy=False)
        )

    def __or__(self, other):
        return type(self).from_scalars(np.array(
            self, copy=False).__or__(np.array(other, copy=False)))

    def __ior__(self, other):
        return type(self).from_scalars(
            np.array(self, copy=False) | np.array(other, copy=False))

    def __and__(self, other):
        return type(self).from_scalars(
            np.array(self, copy=False).__and__(np.array(other, copy=False)))

    def __iand__(self, other):
        return type(self).from_scalars(
            np.array(self, copy=False) & (np.array(other, copy=False)))

    def __getitem__(self, item):
        arr = np.array(self, copy=False)
        if is_scalar(item):
            return arr[item]
        else:
            arr = arr[item]
            return type(self).from_scalars(arr)

    def view(self, dtype=None):
        arr = np.array(self._data, copy=False)
        if dtype is not None:
            arr = arr.view(dtype=dtype)
        return arr

    def sum(self, axis=None, min_count=None):
        return np.array(self, copy=False).sum()

    def copy(self, deep=False):
        if deep:
            return type(self)(copy.deepcopy(self._data))
        else:
            return type(self)(copy.copy(self._data))

    def any(self, axis=0, out=None):
        return np.array(self._data, copy=False).any()

    def all(self, axis=0, out=None):
        return np.array(self._data, copy=False).all()

    def min(self, axis=0, out=None):
        return np.array(self._data, copy=False).min()

    def max(self, axis=0, out=None):
        return np.array(self._data, copy=False).max()

    def _reduce(self, method, skipna=True, **kwargs):
        if skipna:
            arr = self[~self.isna()]
        else:
            arr = self
        # we only allow explicity defined methods
        # ndarrays actually support: mean, var, prod, min, max
        try:
            op = getattr(arr, method)
            return op()
        except AttributeError:
            pass
        raise TypeError
