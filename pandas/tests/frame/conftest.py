import pytest

import numpy as np

from pandas import compat
import pandas.util.testing as tm
from pandas import DataFrame, date_range, NaT


# This module is the start of transitioning from attributes of
# pandas/tests/frame/common.TestData towards fixtures (GH22471).
# Until all modules have been transitioned, it is advised not to change
# the (admittedly suboptimal) names of these fixtures.

@pytest.fixture
def frame():
    """
    Fixture for DataFrame of floats with index of unique strings

    After completing the fixturization of the frame tests (GH 22471), this
    fixture will be renamed to: float_frame_string_index

    Columns are ['A', 'B', 'C', 'D'], see pandas.util.testing.getSeriesData
    """
    return DataFrame(tm.getSeriesData())


@pytest.fixture
def frame2():
    """
    Fixture for DataFrame of floats with index of unique strings

    After completing the fixturization of the frame tests (GH 22471), this
    fixture will be renamed to: float_frame_string_index2

    Columns are ['D', 'C', 'B', 'A']. See pandas.util.testing.getSeriesData
    """
    return DataFrame(tm.getSeriesData(), columns=['D', 'C', 'B', 'A'])


@pytest.fixture
def intframe():
    """
    Fixture for DataFrame of ints with index of unique strings

    After completing the fixturization of the frame tests (GH 22471), this
    fixture will be renamed to: int_frame

    Columns are ['A', 'B', 'C', 'D']. See pandas.util.testing.getSeriesData
    """
    df = DataFrame({k: v.astype(int)
                   for k, v in compat.iteritems(tm.getSeriesData())})
    # force these all to int64 to avoid platform testing issues
    return DataFrame({c: s for c, s in compat.iteritems(df)}, dtype=np.int64)


@pytest.fixture
def tsframe():
    """
    Fixture for DataFrame of floats with DatetimeIndex

    After completing the fixturization of the frame tests (GH 22471), this
    fixture will be renamed to: float_frame_datetime_index

    Columns are ['A', 'B', 'C', 'D']. See pandas.util.testing.getTimeSeriesData
    """
    return DataFrame(tm.getTimeSeriesData())


@pytest.fixture
def mixed_frame():
    """
    Fixture for DataFrame of floats and strings with string index

    After completing the fixturization of the frame tests (GH 22471), this
    fixture will be renamed to: float_string_frame

    Columns are ['A', 'B', 'C', 'D', 'foo'].
    See pandas.util.testing.getSeriesData
    """
    df = DataFrame(tm.getSeriesData())
    df['foo'] = 'bar'
    return df


@pytest.fixture
def mixed_float():
    """
    Fixture for DataFrame of different float types with string index

    After completing the fixturization of the frame tests (GH 22471), this
    fixture will be renamed to: mixed_float_frame

    Columns are ['A', 'B', 'C', 'D']. See pandas.util.testing.getSeriesData
    """
    df = DataFrame(tm.getSeriesData())
    df.A = df.A.astype('float16')
    df.B = df.B.astype('float32')
    df.C = df.C.astype('float64')
    return df


@pytest.fixture
def mixed_float2():
    """
    Fixture for DataFrame of different float types with string index

    After completing the fixturization of the frame tests (GH 22471), this
    fixture will be renamed to: mixed_float_frame2

    Columns are ['A', 'B', 'C', 'D']. See pandas.util.testing.getSeriesData
    """
    df = DataFrame(tm.getSeriesData())
    df.D = df.D.astype('float16')
    df.C = df.C.astype('float32')
    df.B = df.B.astype('float64')
    return df


@pytest.fixture
def mixed_int():
    """
    Fixture for DataFrame of different int types with string index

    After completing the fixturization of the frame tests (GH 22471), this
    fixture will be renamed to: mixed_int_frame

    Columns are ['A', 'B', 'C', 'D']. See pandas.util.testing.getSeriesData
    """
    df = DataFrame({k: v.astype(int)
                   for k, v in compat.iteritems(tm.getSeriesData())})
    df.A = df.A.astype('uint8')
    df.B = df.B.astype('int32')
    df.C = df.C.astype('int64')
    df.D = np.ones(len(df.D), dtype='uint64')
    return df


@pytest.fixture
def all_mixed():
    """
    Fixture for DataFrame of float/int/string columns

    After completing the fixturization of the frame tests (GH 22471), this
    fixture will be renamed to: float_int_string_frame

    Columns are ['a', 'b', 'c', 'float32', 'int32'].
    """
    return DataFrame({'a': 1., 'b': 2, 'c': 'foo',
                      'float32': np.array([1.] * 10, dtype='float32'),
                      'int32': np.array([1] * 10, dtype='int32')},
                     index=np.arange(10))


@pytest.fixture
def tzframe():
    """
    Fixture for DataFrame of date_range Series with different timezones

    After completing the fixturization of the frame tests (GH 22471), this
    fixture will be renamed to: timezone_frame

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
def empty():
    """
    Fixture for empty DataFrame

    After completing the fixturization of the frame tests (GH 22471), this
    fixture will be renamed to: empty_frame
    """
    return DataFrame({})


@pytest.fixture
def ts1():
    """
    Fixture for Series of floats with DatetimeIndex

    After completing the fixturization of the frame tests (GH 22471), this
    fixture will be renamed to: datetime_series

    See pandas.util.testing.makeTimeSeries
    """
    return tm.makeTimeSeries(nper=30)


@pytest.fixture
def ts2():
    """
    Fixture for Series of floats with DatetimeIndex

    After completing the fixturization of the frame tests (GH 22471), this
    fixture will be renamed to: datetime_series_short

    See pandas.util.testing.makeTimeSeries
    """
    return tm.makeTimeSeries(nper=30)[5:]


@pytest.fixture
def simple():
    """
    Fixture for simple 3x3 DataFrame

    After completing the fixturization of the frame tests (GH 22471), this
    fixture will be renamed to: simple_frame

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

    Columns are ['A', 'B', 'C', 'D', 'E']; 'A' & 'B' contain duplicates,
    the rest are unique.
    """
    df = DataFrame({'A': ['foo', 'foo', 'foo', 'bar', 'bar'],
                    'B': ['one', 'two', 'three', 'one', 'two'],
                    'C': ['a', 'b', 'c', 'd', 'e'],
                    'D': np.random.randn(5),
                    'E': np.random.randn(5)})
    return df
