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

    def test_slice_irregular_datetime_index_with_nan(self):
        # GH36953
        index = pd.to_datetime(["2012-01-01", "2012-01-02", "2012-01-03", None])
        df = pd.DataFrame(range(len(index)), index=index)
        expected = pd.DataFrame(range(len(index[:3])), index=index[:3])
        tm.assert_frame_equal(df["2012-01-01":"2012-01-04"], expected)

    def test_slice_datetime_index(self):
        # GH35509
        df = pd.DataFrame(
            {"col1": ["a", "b", "c"], "col2": [1, 2, 3]},
            index=pd.to_datetime(["2020-08-01", "2020-07-02", "2020-08-05"]),
        )
        expected = pd.DataFrame(
            {"col1": ["a", "c"], "col2": [1, 3]},
            index=pd.to_datetime(["2020-08-01", "2020-08-05"]),
        )
        tm.assert_frame_equal(df["2020-08"], expected)
