import collections
import itertools
import numbers
import operator
import random
import string
import sys

import numpy as np
import pytest


from pandas.core.dtypes.base import ExtensionDtype
from pandas.core.arrays import ExtensionArray

from .base import BaseArrayTests, BaseDtypeTests


class JSONDtype(ExtensionDtype):
    type = collections.Mapping
    name = 'json'

    @classmethod
    def construct_from_string(cls, string):
        if string == cls.name:
            return cls()
        else:
            raise TypeError("Cannot construct a '{}' from "
                            "'{}'".format(cls, string))


class JSONArray(ExtensionArray):
    dtype = JSONDtype()

    def __init__(self, values):
        for val in values:
            if not isinstance(val, self.dtype.type):
                raise TypeError
        self.data = values

    def __getitem__(self, item):
        if isinstance(item, numbers.Integral):
            return self.data[item]
        elif isinstance(item, np.ndarray) and item.dtype == 'bool':
            return type(self)([x for x, m in zip(self, item) if m])
        else:
            return type(self)(self.data[item])

    def __setitem__(self, key, value):
        if isinstance(key, numbers.Integral):
            self.data[key] = value
        else:
            if not isinstance(value, (type(self),
                                      collections.Sequence)):
                # broadcast value
                value = itertools.cycle([value])

            if isinstance(key, np.ndarray) and key.dtype == 'bool':
                # masking
                for i, (k, v) in enumerate(zip(key, value)):
                    if k:
                        assert isinstance(v, self.dtype.type)
                        self.data[i] = v
            else:
                for k, v in zip(key, value):
                    assert isinstance(v, self.dtype.type)
                    self.data[k] = v

    def __len__(self):
        return len(self.data)

    def __repr__(self):
        return 'JSONArary({!r})'.format(self.data)

    @property
    def nbytes(self):
        return sys.getsizeof(self.data)

    def isna(self):
        return np.array([x == self._fill_value for x in self.data])

    def take(self, indexer, allow_fill=True, fill_value=None):
        output = [self.data[loc] if loc != -1 else self._fill_value
                  for loc in indexer]
        return type(self)(output)

    def copy(self, deep=False):
        return type(self)(self.data[:])

    @property
    def _fill_value(self):
        return {}

    @classmethod
    def _concat_same_type(cls, to_concat):
        data = list(itertools.chain.from_iterable([x.data for x in to_concat]))
        return cls(data)


def make_data():
    # TODO: Use a regular dict. See _NDFrameIndexer._setitem_with_indexer
    return [collections.UserDict([
        (random.choice(string.ascii_letters), random.randint(0, 100))
        for _ in range(random.randint(0, 10))]) for _ in range(100)]


class TestJSONDtype(BaseDtypeTests):
    @pytest.fixture
    def dtype(self):
        return JSONDtype()


class TestJSON(BaseArrayTests):

    @pytest.fixture
    def data(self):
        """Length-100 PeriodArray for semantics test."""
        return JSONArray(make_data())

    @pytest.fixture
    def data_missing(self):
        """Length 2 array with [NA, Valid]"""
        return JSONArray([{}, {'a': 10}])

    @pytest.fixture
    def na_cmp(self):
        return operator.eq

    @pytest.mark.skip(reason="Unhashable")
    def test_value_counts(self, all_data, dropna):
        pass

    # @pytest.mark.xfail(reason="Difficulty setting sized objects.")
    # def test_set_scalar(self):
    #     pass
    #

    @pytest.mark.xfail(reason="Difficulty setting sized objects.")
    def test_set_loc_scalar_mixed(self):
        # This fails on an np.ndarary(dict) call in _setitem_with_indexer
        pass

    # @pytest.mark.xfail(reason="Difficulty setting sized objects.")
    # def test_set_loc_scalar_single(self):
    #     pass
    #

    @pytest.mark.xfail(reason="Difficulty setting sized objects.")
    def test_set_loc_scalar_multiple_homogoneous(self):
        # This fails in _setitem_with_indexer with a
        # ValueError: Must have equal len keys and value when setting with
        # and iterable
        pass

    @pytest.mark.xfail(reason="Difficulty setting sized objects.")
    def test_set_iloc_scalar_mixed(self):
        # This fails in _setitem_with_indexer with a
        # ValueError: Must have equal len keys and value when setting with an
        # iterable
        pass

    # @pytest.mark.xfail(reason="Difficulty setting sized objects.")
    # def test_set_iloc_scalar_single(self):
    #     pass
    #
    @pytest.mark.xfail(reason="Difficulty setting sized objects.")
    def test_set_iloc_scalar_multiple_homogoneous(self):
        # this fails in _setitem_with_indexer with a
        # ValueError: Must have equal len keys and value when setting with an
        # iterable
        pass
