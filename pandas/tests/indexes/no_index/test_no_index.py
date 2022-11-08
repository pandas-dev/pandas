# gosh, so much to update here...
# prob, inherit from rangeindex tests?

import numpy as np
import pytest

from pandas._config import (
    get_option,
    option_context,
    set_option,
)

from pandas.core.dtypes.common import ensure_platform_int

import pandas as pd
from pandas import (
    CategoricalIndex,
    DatetimeIndex,
    Index,
    IntervalIndex,
    MultiIndex,
    PeriodIndex,
    RangeIndex,
    Series,
    TimedeltaIndex,
    isna,
)
import pandas._testing as tm
from pandas.core.api import (  # noqa:F401
    Float64Index,
    Int64Index,
    NumericIndex,
    UInt64Index,
)
from pandas.core.arrays import BaseMaskedArray
from pandas.core.indexes.api import (
    Float64Index,
    Index,
    Int64Index,
    NoIndex,
    RangeIndex,
)
from pandas.tests.indexes.common import NumericBase
from pandas.tests.indexes.ranges.test_range import TestRangeIndex

# aliases to make some tests easier to read
NI = NoIndex
RI = RangeIndex
I64 = Int64Index
F64 = Float64Index
OI = Index


# class TestNoIndex(TestRangeIndex):
#     _index_cls = NoIndex


@pytest.fixture
def ser1():
    with option_context("mode.no_default_index", True):
        res = pd.Series([1, 2, 3])
    return res


@pytest.fixture
def ser2():
    with option_context("mode.no_default_index", True):
        res = pd.Series([4, 5, 6])
    return res


@pytest.fixture
def df1():
    with option_context("mode.no_default_index", True):
        res = pd.DataFrame({"a": [1, 2, 3], "b": [4, 5, 6]})
    return res


@pytest.fixture
def df2():
    with option_context("mode.no_default_index", True):
        res = pd.DataFrame({"a": [6, 5, 4], "b": [1, 4, 2]})
    return res


def test_boolean_mask(ser1):
    mask = ser1 > 1
    result = ser1[mask]
    expected = pd.Series([2, 3], index=NoIndex(2))
    tm.assert_series_equal(result, expected)

    ser1[mask] = 4
    expected = pd.Series([1, 4, 4], index=NoIndex(3))
    tm.assert_series_equal(ser1, expected)


def test_join(ser1, ser2, df1, df2):
    result = df1.join(df2, lsuffix="_df1")
    expected = pd.DataFrame(
        {
            "a_df1": [1, 2, 3],
            "b_df1": [1, 2, 3],
            "a": [6, 5, 4],
            "b": [1, 4, 2],
        },
        index=NoIndex(3),
    )


# def test_loc(ser1):

# ser1[mask] = 2
# print(ser1)

# df = pd.DataFrame({'a': [1,2,3]})
# try:
#     df.loc[0, 'a']  # raise!
# except TypeError as err:
#     print('error!')
#     print(repr(err))
#     pass
# else:
#     print('fail')
# try:
#     df.loc[0]  # raise!
# except TypeError as err:
#     print('error!')
#     print(repr(err))
#     pass
# else:
#     print('fail')
# print(df.loc[:, 'a'])  # work
# print()
# mask = df['a'] > 2
# print(df.loc[mask])
# print()
# print(df.loc[df['a']>2, 'a'])  # work!

# print(df.iloc[1:]) # work!

# try:
#     _sum =df + df.iloc[1:]
# except TypeError as err:
#     print(repr(err))
# else:
#     print('fail')

# df = pd.read_csv(io.StringIO('data\n1\n'))
# print(df.index)
# df = pd.DataFrame({'a': [1,2,3]*50})
# repr(df)
# print(df)

# df = DataFrame({'a': [1,2,3]})
# concatted = pd.concat([df, df])
# print(concatted.index)
# print(concatted)

# print(df.merge(df, on='a'))

# pd.DataFrame([[1, 2, 3]])
