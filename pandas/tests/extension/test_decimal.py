import decimal
import numbers
import random
import sys

import numpy as np
import pandas as pd
import pandas.util.testing as tm
import pytest

from pandas.core.arrays import ExtensionArray
from pandas.core.dtypes.base import ExtensionDtype

from .base import BaseDtypeTests, BaseArrayTests


class DecimalDtype(ExtensionDtype):
    type = decimal.Decimal
    name = 'decimal'

    @classmethod
    def construct_from_string(cls, string):
        if string == cls.name:
            return cls()
        else:
            raise TypeError("Cannot construct a '{}' from "
                            "'{}'".format(cls, string))


class DecimalArray(ExtensionArray):
    dtype = DecimalDtype()

    def __init__(self, values):
        values = np.asarray(values, dtype=object)

        self.values = values

    def __getitem__(self, item):
        if isinstance(item, numbers.Integral):
            return self.values[item]
        elif isinstance(item, np.ndarray) and item.dtype == 'bool':
            return type(self)([x for x, m in zip(self, item) if m])
        else:
            return type(self)(self.values[item])

    def copy(self, deep=False):
        if deep:
            return type(self)(self.values.copy())
        return type(self)(self)

    def __setitem__(self, key, value):
        if pd.api.types.is_list_like(value):
            value = [decimal.Decimal(v) for v in value]
        else:
            value = decimal.Decimal(value)
        self.values[key] = value

    def __len__(self):
        return len(self.values)

    def __repr__(self):
        return repr(self.values)

    @property
    def nbytes(self):
        n = len(self)
        if n:
            return n * sys.getsizeof(self[0])
        return 0

    def isna(self):
        return np.array([x.is_nan() for x in self.values])

    def take(self, indexer, allow_fill=True, fill_value=None):
        mask = indexer == -1

        out = self.values.take(indexer)
        out[mask] = self._fill_value

        return type(self)(out)

    @property
    def _fill_value(self):
        return decimal.Decimal('NaN')

    @classmethod
    def _concat_same_type(cls, to_concat):
        return cls(np.concatenate([x.values for x in to_concat]))


def make_data():
    return [decimal.Decimal(random.random()) for _ in range(100)]


class TestDecimalDtype(BaseDtypeTests):

    @pytest.fixture
    def dtype(self):
        return DecimalDtype()


class TestDecimalArray(BaseArrayTests):

    @pytest.fixture
    def data(self):
        return DecimalArray(make_data())

    @pytest.fixture
    def data_missing(self):
        return DecimalArray([decimal.Decimal('NaN'), decimal.Decimal(1)])

    @pytest.fixture
    def na_cmp(self):
        return lambda x, y: x.is_nan() and y.is_nan()

    def test_align(self, data):
        a = data[:3]
        b = data[2:5]
        r1, r2 = pd.Series(a).align(pd.Series(b, index=[1, 2, 3]))

        # NaN handling
        e1 = pd.Series(type(data)(list(a) + [data._fill_value]))
        e2 = pd.Series(type(data)([data._fill_value] + list(b)))
        tm.assert_series_equal(r1.iloc[:3], e1.iloc[:3])
        assert r1[3].is_nan()
        assert e1[3].is_nan()

        tm.assert_series_equal(r2.iloc[1:], e2.iloc[1:])
        assert r2[0].is_nan()
        assert e2[0].is_nan()

    @pytest.mark.skip(reason="NaN Sorting")
    def test_value_counts(self, all_data, dropna):
        all_data = all_data[:10]
        if dropna:
            other = np.array(all_data[~all_data.isna()])
        else:
            other = all_data

        result = pd.Series(all_data).value_counts(dropna=dropna).sort_index()
        expected = pd.Series(other).value_counts(dropna=dropna).sort_index()

        tm.assert_series_equal(result, expected)


def test_series_constructor_with_dtype_coercion_raises():
    xpr = ("Cannot coerce data to extension dtype 'decimal'. Pass the "
           "extension array for 'decimal' directly instead.")
    with tm.assert_raises_regex(ValueError, xpr):
        pd.Series([0, 1, 2], dtype=DecimalDtype())


def test_series_constructor_with_same_dtype_ok():
    arr = DecimalArray([decimal.Decimal('10.0')])
    result = pd.Series(arr, dtype=DecimalDtype())
    expected = pd.Series(arr)
    tm.assert_series_equal(result, expected)


def test_series_constructor_with_different_dtype_raises():
    arr = DecimalArray([decimal.Decimal('10.0')])
    xpr = "Cannot specify a dtype 'int64' .* \('decimal'\)."

    with tm.assert_raises_regex(ValueError, xpr):
        pd.Series(arr, dtype='int64')


def test_dataframe_constructor_with_same_dtype_ok():
    arr = DecimalArray([decimal.Decimal('10.0')])

    result = pd.DataFrame({"A": arr}, dtype=DecimalDtype())
    expected = pd.DataFrame({"A": arr})
    tm.assert_frame_equal(result, expected)


def test_dataframe_constructor_with_different_dtype_raises():
    arr = DecimalArray([decimal.Decimal('10.0')])

    xpr = "Cannot coerce extension array to dtype 'int64'. "
    with tm.assert_raises_regex(ValueError, xpr):
        pd.DataFrame({"A": arr}, dtype='int64')
