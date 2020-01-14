import numpy as np
import pytest

from pandas import DataFrame, Index, MultiIndex
import pandas._testing as tm


@pytest.fixture
def multiindex_dataframe_random_data():
    """DataFrame with 2 level MultiIndex with random data"""
    index = MultiIndex(
        levels=[["foo", "bar", "baz", "qux"], ["one", "two", "three"]],
        codes=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3], [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
        names=["first", "second"],
    )
    return DataFrame(
        np.random.randn(10, 3), index=index, columns=Index(["A", "B", "C"], name="exp")
    )


@pytest.fixture
def multiindex_year_month_day_dataframe_random_data():
    """DataFrame with 3 level MultiIndex (year, month, day) covering
    first 100 business days from 2000-01-01 with random data"""
    tdf = tm.makeTimeDataFrame(100)
    ymd = tdf.groupby([lambda x: x.year, lambda x: x.month, lambda x: x.day]).sum()
    # use Int64Index, to make sure things work
    ymd.index.set_levels([lev.astype("i8") for lev in ymd.index.levels], inplace=True)
    ymd.index.set_names(["year", "month", "day"], inplace=True)
    return ymd
