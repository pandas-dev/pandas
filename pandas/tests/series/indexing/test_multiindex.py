""" test get/set & misc """

import pytest

import pandas as pd
from pandas import MultiIndex, Series
import pandas._testing as tm


def test_access_none_value_in_multiindex():
    # GH34318: test that you can access a None value using .loc through a Multiindex

    s = Series([None], pd.MultiIndex.from_arrays([["Level1"], ["Level2"]]))
    result = s.loc[("Level1", "Level2")]
    assert result is None

    midx = MultiIndex.from_product([["Level1"], ["Level2_a", "Level2_b"]])
    s = Series([None] * len(midx), dtype=object, index=midx)
    result = s.loc[("Level1", "Level2_a")]
    assert result is None

    s = Series([1] * len(midx), dtype=object, index=midx)
    result = s.loc[("Level1", "Level2_a")]
    assert result == 1


@pytest.mark.parametrize(
    "ix_data, exp_data",
    [
        (
            [(pd.NaT, 1), (pd.NaT, 2)],
            {"a": [pd.NaT, pd.NaT], "b": [1, 2], "x": [11, 12]},
        ),
        (
            [(pd.NaT, 1), (pd.Timestamp("2020-01-01"), 2)],
            {"a": [pd.NaT, pd.Timestamp("2020-01-01")], "b": [1, 2], "x": [11, 12]},
        ),
        (
            [(pd.NaT, 1), (pd.Timedelta(123, "d"), 2)],
            {"a": [pd.NaT, pd.Timedelta(123, "d")], "b": [1, 2], "x": [11, 12]},
        ),
    ],
)
def test_nat_multi_index(ix_data, exp_data):
    # GH36541: that reset_index() does not raise ValueError
    ix = pd.MultiIndex.from_tuples(ix_data, names=["a", "b"])
    result = pd.DataFrame({"x": [11, 12]}, index=ix)
    result = result.reset_index()

    expected = pd.DataFrame(exp_data)
    tm.assert_frame_equal(result, expected)


def test_loc_getitem_multiindex_nonunique_len_zero():
    # GH#13691
    mi = pd.MultiIndex.from_product([[0], [1, 1]])
    ser = Series(0, index=mi)

    res = ser.loc[[]]

    expected = ser[:0]
    tm.assert_series_equal(res, expected)

    res2 = ser.loc[ser.iloc[0:0]]
    tm.assert_series_equal(res2, expected)
