# gosh, so much to update here...
# prob, inherit from rangeindex tests?

import numpy as np
import pytest

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


class TestNoIndex(TestRangeIndex):
    _index_cls = NoIndex


@pytest.fixture
def ser1():
    yield pd.Series([1, 2, 3])


def test_boolean_mask(ser1):
    # we need to register an option here...
    mask = ser1 > 1
    result = ser1[mask]
    expected = pd.Series([2, 3], index=NoIndex(2))

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
