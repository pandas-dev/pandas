"""Test extension array for storing nested data in a pandas container.

The JSONArray stores lists of dictionaries. The storage mechanism is a list,
not an ndarray.

Note:

We currently store lists of UserDicts. Pandas has a few places
internally that specifically check for dicts, and does non-scalar things
in that case. We *want* the dictionaries to be treated as scalars, so we
hack around pandas by using UserDicts.
"""
from collections import UserDict, abc
import itertools
import numbers
import random
import string
import sys
from typing import Any, Mapping, Type

import numpy as np

import pandas as pd
from pandas.api.extensions import ExtensionArray, ExtensionDtype


class JSONDtype(ExtensionDtype):
    type = abc.Mapping
    name = "json"
    na_value: Mapping[str, Any] = UserDict()

    @classmethod
    def construct_array_type(cls) -> Type["JSONArray"]:
        """
        Return the array type associated with this dtype.

        Returns
        -------
        type
        """
        return JSONArray


class JSONArray(ExtensionArray):
    dtype = JSONDtype()
    __array_priority__ = 1000

    def __init__(self, values, dtype=None, copy=False):
        for val in values:
            if not isinstance(val, self.dtype.type):
                raise TypeError("All values must be of type " + str(self.dtype.type))
        self.data = values

        # Some aliases for common attribute names to ensure pandas supports
        # these
        self._items = self._data = self.data
        # those aliases are currently not working due to assumptions
        # in internal code (GH-20735)
        # self._values = self.values = self.data

    @classmethod
    def _from_sequence(cls, scalars, dtype=None, copy=False):
        return cls(scalars)

    @classmethod
    def _from_factorized(cls, values, original):
        return cls([UserDict(x) for x in values if x != ()])

    def __getitem__(self, item):
        if isinstance(item, numbers.Integral):
            return self.data[item]
        elif isinstance(item, slice) and item == slice(None):
            # Make sure we get a view
            return type(self)(self.data)
        elif isinstance(item, slice):
            # slice
            return type(self)(self.data[item])
        else:
            item = pd.api.indexers.check_array_indexer(self, item)
            if pd.api.types.is_bool_dtype(item.dtype):
                return self._from_sequence([x for x, m in zip(self, item) if m])
            # integer
            return type(self)([self.data[i] for i in item])

    def __setitem__(self, key, value):
        if isinstance(key, numbers.Integral):
            self.data[key] = value
        else:
            if not isinstance(value, (type(self), abc.Sequence)):
                # broadcast value
                value = itertools.cycle([value])

            if isinstance(key, np.ndarray) and key.dtype == "bool":
                # masking
                for i, (k, v) in enumerate(zip(key, value)):
                    if k:
                        assert isinstance(v, self.dtype.type)
                        self.data[i] = v
            else:
                for k, v in zip(key, value):
                    assert isinstance(v, self.dtype.type)
                    self.data[k] = v

    def __len__(self) -> int:
        return len(self.data)

    def __array__(self, dtype=None):
        if dtype is None:
            dtype = object
        return np.asarray(self.data, dtype=dtype)

    @property
    def nbytes(self) -> int:
        return sys.getsizeof(self.data)

    def isna(self):
        return np.array([x == self.dtype.na_value for x in self.data], dtype=bool)

    def take(self, indexer, allow_fill=False, fill_value=None):
        # re-implement here, since NumPy has trouble setting
        # sized objects like UserDicts into scalar slots of
        # an ndarary.
        indexer = np.asarray(indexer)
        msg = (
            "Index is out of bounds or cannot do a "
            "non-empty take from an empty array."
        )

        if allow_fill:
            if fill_value is None:
                fill_value = self.dtype.na_value
            # bounds check
            if (indexer < -1).any():
                raise ValueError
            try:
                output = [
                    self.data[loc] if loc != -1 else fill_value for loc in indexer
                ]
            except IndexError:
                raise IndexError(msg)
        else:
            try:
                output = [self.data[loc] for loc in indexer]
            except IndexError:
                raise IndexError(msg)

        return self._from_sequence(output)

    def copy(self):
        return type(self)(self.data[:])

    def astype(self, dtype, copy=True):
        # NumPy has issues when all the dicts are the same length.
        # np.array([UserDict(...), UserDict(...)]) fails,
        # but np.array([{...}, {...}]) works, so cast.

        # needed to add this check for the Series constructor
        if isinstance(dtype, type(self.dtype)) and dtype == self.dtype:
            if copy:
                return self.copy()
            return self
        return np.array([dict(x) for x in self], dtype=dtype, copy=copy)

    def unique(self):
        # Parent method doesn't work since np.array will try to infer
        # a 2-dim object.
        return type(self)(
            [dict(x) for x in list({tuple(d.items()) for d in self.data})]
        )

    @classmethod
    def _concat_same_type(cls, to_concat):
        data = list(itertools.chain.from_iterable([x.data for x in to_concat]))
        return cls(data)

    def _values_for_factorize(self):
        frozen = self._values_for_argsort()
        if len(frozen) == 0:
            # _factorize_array expects 1-d array, this is a len-0 2-d array.
            frozen = frozen.ravel()
        return frozen, ()

    def _values_for_argsort(self):
        # Disable NumPy's shape inference by including an empty tuple...
        # If all the elements of self are the same size P, NumPy will
        # cast them to an (N, P) array, instead of an (N,) array of tuples.
        frozen = [()] + [tuple(x.items()) for x in self]
        return np.array(frozen, dtype=object)[1:]


def make_data():
    # TODO: Use a regular dict. See _NDFrameIndexer._setitem_with_indexer
    return [
        UserDict(
            [
                (random.choice(string.ascii_letters), random.randint(0, 100))
                for _ in range(random.randint(0, 10))
            ]
        )
        for _ in range(100)
    ]
