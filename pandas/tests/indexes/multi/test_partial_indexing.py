import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, MultiIndex, date_range
import pandas._testing as tm


def test_partial_string_timestamp_multiindex():
    # GH10331
    dr = pd.date_range("2016-01-01", "2016-01-03", freq="12H")
    abc = ["a", "b", "c"]
    ix = pd.MultiIndex.from_product([dr, abc])
    df = pd.DataFrame({"c1": range(0, 15)}, index=ix)
    idx = pd.IndexSlice

    #                        c1
    # 2016-01-01 00:00:00 a   0
    #                     b   1
    #                     c   2
    # 2016-01-01 12:00:00 a   3
    #                     b   4
    #                     c   5
    # 2016-01-02 00:00:00 a   6
    #                     b   7
    #                     c   8
    # 2016-01-02 12:00:00 a   9
    #                     b  10
    #                     c  11
    # 2016-01-03 00:00:00 a  12
    #                     b  13
    #                     c  14

    # partial string matching on a single index
    for df_swap in (df.swaplevel(), df.swaplevel(0), df.swaplevel(0, 1)):
        df_swap = df_swap.sort_index()
        just_a = df_swap.loc["a"]
        result = just_a.loc["2016-01-01"]
        expected = df.loc[idx[:, "a"], :].iloc[0:2]
        expected.index = expected.index.droplevel(1)
        tm.assert_frame_equal(result, expected)

    # indexing with IndexSlice
    result = df.loc[idx["2016-01-01":"2016-02-01", :], :]
    expected = df
    tm.assert_frame_equal(result, expected)

    # match on secondary index
    result = df_swap.loc[idx[:, "2016-01-01":"2016-01-01"], :]
    expected = df_swap.iloc[[0, 1, 5, 6, 10, 11]]
    tm.assert_frame_equal(result, expected)

    # Even though this syntax works on a single index, this is somewhat
    # ambiguous and we don't want to extend this behavior forward to work
    # in multi-indexes. This would amount to selecting a scalar from a
    # column.
    with pytest.raises(KeyError, match="'2016-01-01'"):
        df["2016-01-01"]

    # partial string match on year only
    result = df.loc["2016"]
    expected = df
    tm.assert_frame_equal(result, expected)

    # partial string match on date
    result = df.loc["2016-01-01"]
    expected = df.iloc[0:6]
    tm.assert_frame_equal(result, expected)

    # partial string match on date and hour, from middle
    result = df.loc["2016-01-02 12"]
    expected = df.iloc[9:12]
    tm.assert_frame_equal(result, expected)

    # partial string match on secondary index
    result = df_swap.loc[idx[:, "2016-01-02"], :]
    expected = df_swap.iloc[[2, 3, 7, 8, 12, 13]]
    tm.assert_frame_equal(result, expected)

    # tuple selector with partial string match on date
    result = df.loc[("2016-01-01", "a"), :]
    expected = df.iloc[[0, 3]]
    tm.assert_frame_equal(result, expected)

    # Slicing date on first level should break (of course)
    with pytest.raises(KeyError, match="'2016-01-01'"):
        df_swap.loc["2016-01-01"]

    # GH12685 (partial string with daily resolution or below)
    dr = date_range("2013-01-01", periods=100, freq="D")
    ix = MultiIndex.from_product([dr, ["a", "b"]])
    df = DataFrame(np.random.randn(200, 1), columns=["A"], index=ix)

    result = df.loc[idx["2013-03":"2013-03", :], :]
    expected = df.iloc[118:180]
    tm.assert_frame_equal(result, expected)
