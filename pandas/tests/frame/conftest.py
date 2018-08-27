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
    return DataFrame(tm.getSeriesData())


@pytest.fixture
def frame2():
    return DataFrame(tm.getSeriesData(), columns=['D', 'C', 'B', 'A'])


@pytest.fixture
def intframe():
    df = DataFrame({k: v.astype(int)
                   for k, v in compat.iteritems(tm.getSeriesData())})
    # force these all to int64 to avoid platform testing issues
    return DataFrame({c: s for c, s in compat.iteritems(df)}, dtype=np.int64)


@pytest.fixture
def tsframe():
    return DataFrame(tm.getTimeSeriesData())


@pytest.fixture
def mixed_frame():
    df = DataFrame(tm.getSeriesData())
    df['foo'] = 'bar'
    return df


@pytest.fixture
def mixed_float():
    df = DataFrame(tm.getSeriesData())
    df.A = df.A.astype('float16')
    df.B = df.B.astype('float32')
    df.C = df.C.astype('float64')
    return df


@pytest.fixture
def mixed_float2():
    df = DataFrame(tm.getSeriesData())
    df.D = df.D.astype('float16')
    df.C = df.C.astype('float32')
    df.B = df.B.astype('float64')
    return df


@pytest.fixture
def mixed_int():
    df = DataFrame({k: v.astype(int)
                   for k, v in compat.iteritems(tm.getSeriesData())})
    df.A = df.A.astype('uint8')
    df.B = df.B.astype('int32')
    df.C = df.C.astype('int64')
    df.D = np.ones(len(df.D), dtype='uint64')
    return df


@pytest.fixture
def all_mixed():
    return DataFrame({'a': 1., 'b': 2, 'c': 'foo',
                      'float32': np.array([1.] * 10, dtype='float32'),
                      'int32': np.array([1] * 10, dtype='int32')},
                     index=np.arange(10))


@pytest.fixture
def tzframe():
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
    return DataFrame({})


@pytest.fixture
def ts1():
    return tm.makeTimeSeries(nper=30)


@pytest.fixture
def ts2():
    return tm.makeTimeSeries(nper=30)[5:]


@pytest.fixture
def simple():
    arr = np.array([[1., 2., 3.],
                    [4., 5., 6.],
                    [7., 8., 9.]])

    return DataFrame(arr, columns=['one', 'two', 'three'],
                     index=['a', 'b', 'c'])


@pytest.fixture
def frame_of_index_cols():
    df = DataFrame({'A': ['foo', 'foo', 'foo', 'bar', 'bar'],
                    'B': ['one', 'two', 'three', 'one', 'two'],
                    'C': ['a', 'b', 'c', 'd', 'e'],
                    'D': np.random.randn(5),
                    'E': np.random.randn(5)})
    return df
