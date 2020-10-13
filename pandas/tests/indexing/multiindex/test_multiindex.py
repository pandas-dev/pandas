import numpy as np
import pytest

import pandas._libs.index as _index
from pandas.errors import PerformanceWarning

import pandas as pd
from pandas import DataFrame, Index, MultiIndex, Series
import pandas._testing as tm


class TestMultiIndexBasic:
    def test_multiindex_perf_warn(self):

        df = DataFrame(
            {
                "jim": [0, 0, 1, 1],
                "joe": ["x", "x", "z", "y"],
                "jolie": np.random.rand(4),
            }
        ).set_index(["jim", "joe"])

        with tm.assert_produces_warning(PerformanceWarning):
            df.loc[(1, "z")]

        df = df.iloc[[2, 1, 3, 0]]
        with tm.assert_produces_warning(PerformanceWarning):
            df.loc[(0,)]

    def test_indexing_over_hashtable_size_cutoff(self):
        n = 10000

        old_cutoff = _index._SIZE_CUTOFF
        _index._SIZE_CUTOFF = 20000

        s = Series(np.arange(n), MultiIndex.from_arrays((["a"] * n, np.arange(n))))

        # hai it works!
        assert s[("a", 5)] == 5
        assert s[("a", 6)] == 6
        assert s[("a", 7)] == 7

        _index._SIZE_CUTOFF = old_cutoff

    def test_multi_nan_indexing(self):

        # GH 3588
        df = DataFrame(
            {
                "a": ["R1", "R2", np.nan, "R4"],
                "b": ["C1", "C2", "C3", "C4"],
                "c": [10, 15, np.nan, 20],
            }
        )
        result = df.set_index(["a", "b"], drop=False)
        expected = DataFrame(
            {
                "a": ["R1", "R2", np.nan, "R4"],
                "b": ["C1", "C2", "C3", "C4"],
                "c": [10, 15, np.nan, 20],
            },
            index=[
                Index(["R1", "R2", np.nan, "R4"], name="a"),
                Index(["C1", "C2", "C3", "C4"], name="b"),
            ],
        )
        tm.assert_frame_equal(result, expected)

    def test_nested_tuples_duplicates(self):
        # GH#30892

        dti = pd.to_datetime(["20190101", "20190101", "20190102"])
        idx = pd.Index(["a", "a", "c"])
        mi = pd.MultiIndex.from_arrays([dti, idx], names=["index1", "index2"])

        df = pd.DataFrame({"c1": [1, 2, 3], "c2": [np.nan, np.nan, np.nan]}, index=mi)

        expected = pd.DataFrame({"c1": df["c1"], "c2": [1.0, 1.0, np.nan]}, index=mi)

        df2 = df.copy(deep=True)
        df2.loc[(dti[0], "a"), "c2"] = 1.0
        tm.assert_frame_equal(df2, expected)

        df3 = df.copy(deep=True)
        df3.loc[[(dti[0], "a")], "c2"] = 1.0
        tm.assert_frame_equal(df3, expected)

    def test_multiindex_get_loc_list_raises(self):
        # https://github.com/pandas-dev/pandas/issues/35878
        idx = pd.MultiIndex.from_tuples([("a", 1), ("b", 2)])
        msg = "unhashable type"
        with pytest.raises(TypeError, match=msg):
            idx.get_loc([])


def test_combine_first_with_nan_index():
    mi1 = pd.MultiIndex.from_arrays(
        [["b", "b", "c", "a", "b", np.nan], [1, 2, 3, 4, 5, 6]], names=["a", "b"]
    )
    df = pd.DataFrame({"c": [1, 1, 1, 1, 1, 1]}, index=mi1)
    mi2 = pd.MultiIndex.from_arrays(
        [["a", "b", "c", "a", "b", "d"], [1, 1, 1, 1, 1, 1]], names=["a", "b"]
    )
    s = pd.Series([1, 2, 3, 4, 5, 6], index=mi2)
    df_combined = df.combine_first(pd.DataFrame({"col": s}))
    mi_expected = pd.MultiIndex.from_arrays(
        [
            ["a", "a", "a", "b", "b", "b", "b", "c", "c", "d", np.nan],
            [1, 1, 4, 1, 1, 2, 5, 1, 3, 1, 6],
        ],
        names=["a", "b"],
    )
    assert (df_combined.index == mi_expected).all()
    exp_col = np.asarray(
        [1.0, 4.0, np.nan, 2.0, 5.0, np.nan, np.nan, 3.0, np.nan, 6.0, np.nan]
    )
    act_col = df_combined["col"].values
    assert np.allclose(act_col, exp_col, rtol=0, atol=0, equal_nan=True)
