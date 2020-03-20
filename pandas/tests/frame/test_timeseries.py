import numpy as np

import pandas as pd
from pandas import DataFrame, date_range, to_datetime
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
