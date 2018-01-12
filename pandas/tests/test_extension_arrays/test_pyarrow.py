import collections

import pyarrow as pa
import pytest

import numpy as np
import pandas as pd
from pandas.core.extensions import ExtensionArray, ExtensionDtype
from .base import BaseArrayTests


class MyDtypeType(type):
    pass


class ArrowDtype(ExtensionDtype):
    _can_hold_na = True
    type = MyDtypeType
    base = None
    name = 'pa64'
    arrow_type = pa.int64()


class ArrowArray(ExtensionArray):
    dtype = ArrowDtype()
    ndim = 1
    can_hold_na = True

    def __init__(self, values):
        if not isinstance(values, pa.Array):
            values = pa.array(values)
        assert values.type == self.dtype.arrow_type
        self.data = values

    def __iter__(self):
        return iter(self.data)

    def __len__(self):
        return len(self.data)

    def __getitem__(self, item):
        result = self.data[item]
        if isinstance(item, (slice, collections.Sequence)):
            return type(self)(result)
        else:
            return result

    @property
    def nbytes(self):
        return 64 * len(self)

    @property
    def shape(self):
        return (len(self),)

    def take(self, indexer, allow_fill=True, fill_value=None):
        return type(self)(self.data.to_pandas().take(indexer))

    take_nd = take

    def copy(self):
        # TODO: Jira for pa.array(pyarrow.array)
        return pa.array(self.data.to_pandas())

    def isna(self):
        # https://github.com/apache/arrow/pull/1378
        return pd.isna(self.data.to_pandas())

    @classmethod
    def concat_same_type(cls, to_concat):
        return cls(np.concatenate([arr.data.to_pandas() for arr in to_concat]))

    def get_values(self):
        return self.data

    def to_dense(self):
        return self.data

    def formatting_values(self):
        return self.data.to_pandas()

    def slice(self, indexer):
        return self[indexer]


@pytest.fixture
def test_data():
    """Length-100 int64 arrow array for semantics test."""
    return ArrowArray(np.arange(100))


class TestArrow(BaseArrayTests):
    def test_iloc(self, test_data):
        ser = pd.Series(test_data)
        result = ser.iloc[:4]
        expected = test_data[:4]
        assert isinstance(result, pd.Series)
        assert result.values.data.equals(expected.data)

    def test_loc(self, test_data):
        ser = pd.Series(test_data)
        result = ser.loc[[0, 1, 2, 3]]
        expected = test_data[:4]
        assert isinstance(result, pd.Series)
        assert result.values.data.equals(expected.data)
