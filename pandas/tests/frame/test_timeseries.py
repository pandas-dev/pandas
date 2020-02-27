import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, Series, date_range, to_datetime
import pandas._testing as tm


class TestDataFrameTimeSeriesMethods:
    def test_frame_ctor_datetime64_column(self):
        rng = date_range("1/1/2000 00:00:00", "1/1/2000 1:59:50", freq="10s")
        dates = np.asarray(rng)

        df = DataFrame({"A": np.random.randn(len(rng)), "B": dates})
        assert np.issubdtype(df["B"].dtype, np.dtype("M8[ns]"))

    def test_frame_append_datetime64_column(self):
        rng = date_range("1/1/2000 00:00:00", "1/1/2000 1:59:50", freq="10s")
        df = DataFrame(index=np.arange(len(rng)))

        df["A"] = rng
        assert np.issubdtype(df["A"].dtype, np.dtype("M8[ns]"))

    def test_frame_datetime64_pre1900_repr(self):
        df = DataFrame({"year": date_range("1/1/1700", periods=50, freq="A-DEC")})
        # it works!
        repr(df)

    def test_frame_append_datetime64_col_other_units(self):
        n = 100

        units = ["h", "m", "s", "ms", "D", "M", "Y"]

        ns_dtype = np.dtype("M8[ns]")

        for unit in units:
            dtype = np.dtype(f"M8[{unit}]")
            vals = np.arange(n, dtype=np.int64).view(dtype)

            df = DataFrame({"ints": np.arange(n)}, index=np.arange(n))
            df[unit] = vals

            ex_vals = to_datetime(vals.astype("O")).values

            assert df[unit].dtype == ns_dtype
            assert (df[unit].values == ex_vals).all()

        # Test insertion into existing datetime64 column
        df = DataFrame({"ints": np.arange(n)}, index=np.arange(n))
        df["dates"] = np.arange(n, dtype=np.int64).view(ns_dtype)

        for unit in units:
            dtype = np.dtype(f"M8[{unit}]")
            vals = np.arange(n, dtype=np.int64).view(dtype)

            tmp = df.copy()

            tmp["dates"] = vals
            ex_vals = to_datetime(vals.astype("O")).values

            assert (tmp["dates"].values == ex_vals).all()

    @pytest.mark.parametrize(
        "data,idx,expected_first,expected_last",
        [
            ({"A": [1, 2, 3]}, [1, 1, 2], 1, 2),
            ({"A": [1, 2, 3]}, [1, 2, 2], 1, 2),
            ({"A": [1, 2, 3, 4]}, ["d", "d", "d", "d"], "d", "d"),
            ({"A": [1, np.nan, 3]}, [1, 1, 2], 1, 2),
            ({"A": [np.nan, np.nan, 3]}, [1, 1, 2], 2, 2),
            ({"A": [1, np.nan, 3]}, [1, 2, 2], 1, 2),
        ],
    )
    def test_first_last_valid(
        self, float_frame, data, idx, expected_first, expected_last
    ):
        N = len(float_frame.index)
        mat = np.random.randn(N)
        mat[:5] = np.nan
        mat[-5:] = np.nan

        frame = DataFrame({"foo": mat}, index=float_frame.index)
        index = frame.first_valid_index()

        assert index == frame.index[5]

        index = frame.last_valid_index()
        assert index == frame.index[-6]

        # GH12800
        empty = DataFrame()
        assert empty.last_valid_index() is None
        assert empty.first_valid_index() is None

        # GH17400: no valid entries
        frame[:] = np.nan
        assert frame.last_valid_index() is None
        assert frame.first_valid_index() is None

        # GH20499: its preserves freq with holes
        frame.index = date_range("20110101", periods=N, freq="B")
        frame.iloc[1] = 1
        frame.iloc[-2] = 1
        assert frame.first_valid_index() == frame.index[1]
        assert frame.last_valid_index() == frame.index[-2]
        assert frame.first_valid_index().freq == frame.index.freq
        assert frame.last_valid_index().freq == frame.index.freq

        # GH 21441
        df = DataFrame(data, index=idx)
        assert expected_first == df.first_valid_index()
        assert expected_last == df.last_valid_index()

    @pytest.mark.parametrize("klass", [Series, DataFrame])
    def test_first_valid_index_all_nan(self, klass):
        # GH#9752 Series/DataFrame should both return None, not raise
        obj = klass([np.nan])

        assert obj.first_valid_index() is None
        assert obj.iloc[:0].first_valid_index() is None

    def test_operation_on_NaT(self):
        # Both NaT and Timestamp are in DataFrame.
        df = pd.DataFrame({"foo": [pd.NaT, pd.NaT, pd.Timestamp("2012-05-01")]})

        res = df.min()
        exp = pd.Series([pd.Timestamp("2012-05-01")], index=["foo"])
        tm.assert_series_equal(res, exp)

        res = df.max()
        exp = pd.Series([pd.Timestamp("2012-05-01")], index=["foo"])
        tm.assert_series_equal(res, exp)

        # GH12941, only NaTs are in DataFrame.
        df = pd.DataFrame({"foo": [pd.NaT, pd.NaT]})

        res = df.min()
        exp = pd.Series([pd.NaT], index=["foo"])
        tm.assert_series_equal(res, exp)

        res = df.max()
        exp = pd.Series([pd.NaT], index=["foo"])
        tm.assert_series_equal(res, exp)

    def test_datetime_assignment_with_NaT_and_diff_time_units(self):
        # GH 7492
        data_ns = np.array([1, "nat"], dtype="datetime64[ns]")
        result = pd.Series(data_ns).to_frame()
        result["new"] = data_ns
        expected = pd.DataFrame(
            {0: [1, None], "new": [1, None]}, dtype="datetime64[ns]"
        )
        tm.assert_frame_equal(result, expected)
        # OutOfBoundsDatetime error shouldn't occur
        data_s = np.array([1, "nat"], dtype="datetime64[s]")
        result["new"] = data_s
        expected = pd.DataFrame(
            {0: [1, None], "new": [1e9, None]}, dtype="datetime64[ns]"
        )
        tm.assert_frame_equal(result, expected)
