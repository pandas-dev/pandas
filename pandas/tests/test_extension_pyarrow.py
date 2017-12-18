import pyarrow as pa

import numpy as np
import pandas as pd
import pandas.util.testing as tm
from pandas.core.extensions import ExtensionArray, ExtensionDtype


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
        data = pa.array(values)
        assert data.type == self.dtype.arrow_type
        self.data = data

    def __iter__(self):
        return iter(self.data)

    def __len__(self):
        return len(self.data)

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


def test_series_ctor():
    arr = ArrowArray([1, 2, 3])
    result = pd.Series(arr)
    assert result.dtype == arr.dtype
    assert len(arr) == 3


def test_concat_works():
    arr1 = ArrowArray([1, 2, 3])
    arr2 = ArrowArray([4, 5, 6])
    result = pd.concat([pd.Series(arr1),
                        pd.Series(arr2)], ignore_index=True)
    expected = pa.array([1, 2, 3, 4, 5, 6])
    assert result.dtype == arr1.dtype
    assert isinstance(result.values, ArrowArray)
    assert result.values.data.equals(expected)


def test_slice_works():
    ser = pd.Series(ArrowArray([1, 2, 3]))
    result = ser.loc[[0, 1]]
    assert isinstance(result.values, ArrowArray)
