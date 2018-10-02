import pytest

import numpy as np

from pandas import compat
import pandas.util.testing as tm
from pandas import DataFrame, date_range, NaT


@pytest.fixture
def float_frame():
    """
    Fixture for DataFrame of floats with index of unique strings

    Columns are ['A', 'B', 'C', 'D'].
    """
    return DataFrame(tm.getSeriesData())


@pytest.fixture
def float_frame2():
    """
    Fixture for DataFrame of floats with index of unique strings

    Columns are ['D', 'C', 'B', 'A']
    """
    return DataFrame(tm.getSeriesData(), columns=['D', 'C', 'B', 'A'])


@pytest.fixture
def int_frame():
    """
    Fixture for DataFrame of ints with index of unique strings

    Columns are ['A', 'B', 'C', 'D']
    """
    df = DataFrame({k: v.astype(int)
                   for k, v in compat.iteritems(tm.getSeriesData())})
    # force these all to int64 to avoid platform testing issues
    return DataFrame({c: s for c, s in compat.iteritems(df)}, dtype=np.int64)


@pytest.fixture
def datetime_frame():
    """
    Fixture for DataFrame of floats with DatetimeIndex

    Columns are ['A', 'B', 'C', 'D']
    """
    return DataFrame(tm.getTimeSeriesData())


@pytest.fixture
def float_string_frame():
    """
    Fixture for DataFrame of floats and strings with index of unique strings

    Columns are ['A', 'B', 'C', 'D', 'foo'].
    """
    df = DataFrame(tm.getSeriesData())
    df['foo'] = 'bar'
    return df


@pytest.fixture
def mixed_float_frame():
    """
    Fixture for DataFrame of different float types with index of unique strings

    Columns are ['A', 'B', 'C', 'D'].
    """
    df = DataFrame(tm.getSeriesData())
    df.A = df.A.astype('float32')
    df.B = df.B.astype('float32')
    df.C = df.C.astype('float16')
    df.D = df.D.astype('float64')
    return df


@pytest.fixture
def mixed_float_frame2():
    """
    Fixture for DataFrame of different float types with index of unique strings

    Columns are ['A', 'B', 'C', 'D'].
    """
    df = DataFrame(tm.getSeriesData())
    df.D = df.D.astype('float32')
    df.C = df.C.astype('float32')
    df.B = df.B.astype('float16')
    df.D = df.D.astype('float64')
    return df


@pytest.fixture
def mixed_int_frame():
    """
    Fixture for DataFrame of different int types with index of unique strings

    Columns are ['A', 'B', 'C', 'D'].
    """
    df = DataFrame({k: v.astype(int)
                   for k, v in compat.iteritems(tm.getSeriesData())})
    df.A = df.A.astype('int32')
    df.B = np.ones(len(df.B), dtype='uint64')
    df.C = df.C.astype('uint8')
    df.D = df.C.astype('int64')
    return df


@pytest.fixture
def mixed_type_frame():
    """
    Fixture for DataFrame of float/int/string columns with RangeIndex

    Columns are ['a', 'b', 'c', 'float32', 'int32'].
    """
    return DataFrame({'a': 1., 'b': 2, 'c': 'foo',
                      'float32': np.array([1.] * 10, dtype='float32'),
                      'int32': np.array([1] * 10, dtype='int32')},
                     index=np.arange(10))


@pytest.fixture
def timezone_frame():
    """
    Fixture for DataFrame of date_range Series with different time zones

    Columns are ['A', 'B', 'C']; some entries are missing
    """
    df = DataFrame({'A': date_range('20130101', periods=3),
                    'B': date_range('20130101', periods=3,
                                    tz='US/Eastern'),
                    'C': date_range('20130101', periods=3,
                                    tz='CET')})
    df.iloc[1, 1] = NaT
    df.iloc[1, 2] = NaT
    return df


@pytest.fixture
def empty_frame():
    """
    Fixture for empty DataFrame
    """
    return DataFrame({})


@pytest.fixture
def datetime_series():
    """
    Fixture for Series of floats with DatetimeIndex
    """
    return tm.makeTimeSeries(nper=30)


@pytest.fixture
def datetime_series_short():
    """
    Fixture for Series of floats with DatetimeIndex
    """
    return tm.makeTimeSeries(nper=30)[5:]


@pytest.fixture
def simple_frame():
    """
    Fixture for simple 3x3 DataFrame

    Columns are ['one', 'two', 'three'], index is ['a', 'b', 'c'].
    """
    arr = np.array([[1., 2., 3.],
                    [4., 5., 6.],
                    [7., 8., 9.]])

    return DataFrame(arr, columns=['one', 'two', 'three'],
                     index=['a', 'b', 'c'])


@pytest.fixture
def frame_of_index_cols():
    """
    Fixture for DataFrame of columns that can be used for indexing

    Columns are ['A', 'B', 'C', 'D', 'E', ('tuple', 'as', 'label')];
    'A' & 'B' contain duplicates (but are jointly unique), the rest are unique.
    """
    df = DataFrame({'A': ['foo', 'foo', 'foo', 'bar', 'bar'],
                    'B': ['one', 'two', 'three', 'one', 'two'],
                    'C': ['a', 'b', 'c', 'd', 'e'],
                    'D': np.random.randn(5),
                    'E': np.random.randn(5),
                    ('tuple', 'as', 'label'): np.random.randn(5)})
    return df
