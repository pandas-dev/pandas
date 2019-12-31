import numpy as np
import pytest

from pandas.errors import PerformanceWarning

from pandas import DataFrame, MultiIndex
import pandas.util.testing as tm


class TestMultiIndex:
    def test_frame_setitem_loc(self, multiindex_dataframe_random_data):
        frame = multiindex_dataframe_random_data
        frame.loc[("bar", "two"), "B"] = 5
        assert frame.loc[("bar", "two"), "B"] == 5

        # with integer labels
        df = frame.copy()
        df.columns = list(range(3))
        df.loc[("bar", "two"), 1] = 7
        assert df.loc[("bar", "two"), 1] == 7

    def test_loc_general(self):

        # GH 2817
        data = {
            "amount": {0: 700, 1: 600, 2: 222, 3: 333, 4: 444},
            "col": {0: 3.5, 1: 3.5, 2: 4.0, 3: 4.0, 4: 4.0},
            "year": {0: 2012, 1: 2011, 2: 2012, 3: 2012, 4: 2012},
        }
        df = DataFrame(data).set_index(keys=["col", "year"])
        key = 4.0, 2012

        # emits a PerformanceWarning, ok
        with tm.assert_produces_warning(PerformanceWarning):
            tm.assert_frame_equal(df.loc[key], df.iloc[2:])

        # this is ok
        df.sort_index(inplace=True)
        res = df.loc[key]

        # col has float dtype, result should be Float64Index
        index = MultiIndex.from_arrays([[4.0] * 3, [2012] * 3], names=["col", "year"])
        expected = DataFrame({"amount": [222, 333, 444]}, index=index)
        tm.assert_frame_equal(res, expected)

    def test_loc_multiindex_missing_label_raises(self):
        # GH 21593
        df = DataFrame(
            np.random.randn(3, 3),
            columns=[[2, 2, 4], [6, 8, 10]],
            index=[[4, 4, 8], [8, 10, 12]],
        )

        with pytest.raises(KeyError, match=r"^2$"):
            df.loc[2]

    def test_series_loc_getitem_fancy(
        self, multiindex_year_month_day_dataframe_random_data
    ):
        s = multiindex_year_month_day_dataframe_random_data["A"]
        expected = s.reindex(s.index[49:51])
        result = s.loc[[(2000, 3, 10), (2000, 3, 13)]]
        tm.assert_series_equal(result, expected)
