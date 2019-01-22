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
def empty_frame():
    """
    Fixture for empty DataFrame
    """
    return DataFrame({})
