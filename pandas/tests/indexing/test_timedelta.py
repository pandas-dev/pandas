import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


class TestTimedeltaIndexing:
    def test_loc_setitem_bool_mask(self):
        # GH 14946
        df = pd.DataFrame({"x": range(10)})
        df.index = pd.to_timedelta(range(10), unit="s")
        conditions = [df["x"] > 3, df["x"] == 3, df["x"] < 3]
        expected_data = [
            [0, 1, 2, 3, 10, 10, 10, 10, 10, 10],
            [0, 1, 2, 10, 4, 5, 6, 7, 8, 9],
            [10, 10, 10, 3, 4, 5, 6, 7, 8, 9],
        ]
        for cond, data in zip(conditions, expected_data):
            result = df.copy()
            result.loc[cond, "x"] = 10

            expected = pd.DataFrame(
                data,
                index=pd.to_timedelta(range(10), unit="s"),
                columns=["x"],
                dtype="int64",
            )
            tm.assert_frame_equal(expected, result)

    @pytest.mark.parametrize(
        "indexer, expected",
        [
            (0, [20, 1, 2, 3, 4, 5, 6, 7, 8, 9]),
            (slice(4, 8), [0, 1, 2, 3, 20, 20, 20, 20, 8, 9]),
            ([3, 5], [0, 1, 2, 20, 4, 20, 6, 7, 8, 9]),
        ],
    )
    def test_list_like_indexing(self, indexer, expected):
        # GH 16637
        df = pd.DataFrame({"x": range(10)}, dtype="int64")
        df.index = pd.to_timedelta(range(10), unit="s")

        df.loc[df.index[indexer], "x"] = 20

        expected = pd.DataFrame(
            expected,
            index=pd.to_timedelta(range(10), unit="s"),
            columns=["x"],
            dtype="int64",
        )

        tm.assert_frame_equal(expected, df)

    def test_string_indexing(self):
        # GH 16896
        df = pd.DataFrame({"x": range(3)}, index=pd.to_timedelta(range(3), unit="days"))
        expected = df.iloc[0]
        sliced = df.loc["0 days"]
        tm.assert_series_equal(sliced, expected)

    @pytest.mark.parametrize("value", [None, pd.NaT, np.nan])
    def test_setitem_mask_na_value_td64(self, value):
        # issue (#18586)
        series = pd.Series([0, 1, 2], dtype="timedelta64[ns]")
        series[series == series[0]] = value
        expected = pd.Series([pd.NaT, 1, 2], dtype="timedelta64[ns]")
        tm.assert_series_equal(series, expected)

    @pytest.mark.parametrize("value", [None, pd.NaT, np.nan])
    def test_listlike_setitem(self, value):
        # issue (#18586)
        series = pd.Series([0, 1, 2], dtype="timedelta64[ns]")
        series.iloc[0] = value
        expected = pd.Series([pd.NaT, 1, 2], dtype="timedelta64[ns]")
        tm.assert_series_equal(series, expected)

    @pytest.mark.parametrize(
        "start,stop, expected_slice",
        [
            [np.timedelta64(0, "ns"), None, slice(0, 11)],
            [np.timedelta64(1, "D"), np.timedelta64(6, "D"), slice(1, 7)],
            [None, np.timedelta64(4, "D"), slice(0, 5)],
        ],
    )
    def test_numpy_timedelta_scalar_indexing(self, start, stop, expected_slice):
        # GH 20393
        s = pd.Series(range(11), pd.timedelta_range("0 days", "10 days"))
        result = s.loc[slice(start, stop)]
        expected = s.iloc[expected_slice]
        tm.assert_series_equal(result, expected)

    def test_roundtrip_thru_setitem(self):
        # PR 23462
        dt1 = pd.Timedelta(0)
        dt2 = pd.Timedelta(28767471428571405)
        df = pd.DataFrame({"dt": pd.Series([dt1, dt2])})
        df_copy = df.copy()
        s = pd.Series([dt1])

        expected = df["dt"].iloc[1].value
        df.loc[[True, False]] = s
        result = df["dt"].iloc[1].value

        assert expected == result
        tm.assert_frame_equal(df, df_copy)

    def test_loc_str_slicing(self):
        ix = pd.timedelta_range(start="1 day", end="2 days", freq="1H")
        ser = ix.to_series()
        result = ser.loc[:"1 days"]
        expected = ser.iloc[:-1]

        tm.assert_series_equal(result, expected)

    def test_loc_slicing(self):
        ix = pd.timedelta_range(start="1 day", end="2 days", freq="1H")
        ser = ix.to_series()
        result = ser.loc[: ix[-2]]
        expected = ser.iloc[:-1]

        tm.assert_series_equal(result, expected)
