import numpy as np
import pytest

from pandas import DataFrame, NaT, compat, date_range
import pandas.util.testing as tm


@pytest.fixture
def float_frame():
    """
    Fixture for DataFrame of floats with index of unique strings

    Columns are ['A', 'B', 'C', 'D'].
    """
    return DataFrame(tm.getSeriesData())


@pytest.fixture
def float_frame_with_na():
    """
    Fixture for DataFrame of floats with index of unique strings

    Columns are ['A', 'B', 'C', 'D']; some entries are missing
    """
    df = DataFrame(tm.getSeriesData())
    # set some NAs
    df.loc[5:10] = np.nan
    df.loc[15:20, -2:] = np.nan
    return df


@pytest.fixture
def bool_frame_with_na():
    """
    Fixture for DataFrame of booleans with index of unique strings

    Columns are ['A', 'B', 'C', 'D']; some entries are missing
    """
    df = DataFrame(tm.getSeriesData()) > 0
    df = df.astype(object)
    # set some NAs
    df.loc[5:10] = np.nan
    df.loc[15:20, -2:] = np.nan
    return df


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
def empty_frame():
    """
    Fixture for empty DataFrame
    """
    return DataFrame({})
