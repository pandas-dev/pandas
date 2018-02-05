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
            if not isinstance(val, collections.Mapping):
                raise TypeError
        self.data = values

    def __getitem__(self, item):
        # TDOO: fancy indexing
        if isinstance(item, numbers.Integral):
            return self.data[item]
        elif isinstance(item, np.ndarray) and item.dtype == 'bool':
            return type(self)([x for x, m in zip(self, item) if m])
        elif isinstance(item, collections.Sequence):
            return type(self)([self.data[i] for i in item])
        else:
            return type(self)(self.data[item])

    def __len__(self):
        return len(self.data)

    def __iter__(self):
        return iter(self.data)

    def __repr__(self):
        return 'JSONArary({!r})'.format(self.data)

    @property
    def nbytes(self):
        return sys.getsizeof(self.data)

    def isna(self):
        return np.array([x == {} for x in self.data])

    def take(self, indexer, allow_fill=True, fill_value=None):
        output = [self.data[loc] if loc != -1 else {}
                  for loc in indexer]
        return type(self)(output)

    def copy(self, deep=False):
        return type(self)(self.data.copy(deep=deep))

    @property
    def _fill_value(self):
        return {}

    @classmethod
    def _concat_same_type(cls, to_concat):
        data = list(itertools.chain.from_iterable([x.data for x in to_concat]))
        return cls(data)



def make_data():
    return [{random.choice(string.ascii_letters): random.randint(0, 100)
             for _ in range(random.randint(0, 10))}
             for _ in range(100)]


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

    @pytest.mark.skip(reason="Unorderable")
    def test_reduction_orderable(self, data, method):
        pass

    

# def test_concat_mixed_closed_raises():
#     one = IntervalArray.from_breaks([0, 1, 2], closed='left')
#     two = IntervalArray.from_breaks([1, 2, 3], closed='right')
#
#     with tm.assert_raises_regex(ValueError, "Intervals must all be closed"):
#         IntervalArray._concat_same_type([one, two])
#
#
# def test_series_constructor_intervalindex():
#     result = pd.Series(pd.IntervalIndex.from_breaks([0, 1, 2]))
#     assert result.dtype == 'interval[int64]'
#
#
# def dataframe_constructor_intervalindex():
#     result = pd.DataFrame({"A": pd.IntervalIndex.from_breaks([0, 1, 2])})
#     assert result.dtypes['A'] == 'interval[int64]'
#
#
# def dataframe_set_intervalarray():
#     df = pd.DataFrame({"A": [1, 2]})
#     arr = IntervalArray.from_breaks([0, 1, 2])
#     df['B'] = arr
#
#     expected = pd.DataFrame({"A": [1, 2], "B": arr})
#     tm.assert_frame_equal(df, expected)
#
#
# def dataframe_set_intervalindex():
#     df = pd.DataFrame({"A": [1, 2]})
#     arr = pd.IntervalIndex.from_breaks([0, 1, 2])
#     df['B'] = arr
#
#     expected = pd.DataFrame({"A": [1, 2], "B": arr})
#     tm.assert_frame_equal(df, expected)
