import numpy as np
import pytest

import pandas as pd
from pandas import DataFrame, to_datetime
import pandas._testing as tm


class TestDataFrameTimeSeriesMethods:
    @pytest.mark.parametrize("unit", ["h", "m", "s", "ms", "D", "M", "Y"])
    def test_frame_append_datetime64_col_other_units(self, unit):
        n = 100

        ns_dtype = np.dtype("M8[ns]")

        dtype = np.dtype(f"M8[{unit}]")
        vals = np.arange(n, dtype=np.int64).view(dtype)

        df = DataFrame({"ints": np.arange(n)}, index=np.arange(n))
        df[unit] = vals

        ex_vals = to_datetime(vals.astype("O")).values

        assert df[unit].dtype == ns_dtype
        assert (df[unit].values == ex_vals).all()

    @pytest.mark.parametrize("unit", ["h", "m", "s", "ms", "D", "M", "Y"])
    def test_frame_setitem_existing_datetime64_col_other_units(self, unit):
        # Test insertion into existing datetime64 column
        n = 100
        ns_dtype = np.dtype("M8[ns]")

        df = DataFrame({"ints": np.arange(n)}, index=np.arange(n))
        df["dates"] = np.arange(n, dtype=np.int64).view(ns_dtype)

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
        expected = DataFrame({0: [1, None], "new": [1, None]}, dtype="datetime64[ns]")
        tm.assert_frame_equal(result, expected)
        # OutOfBoundsDatetime error shouldn't occur
        data_s = np.array([1, "nat"], dtype="datetime64[s]")
        result["new"] = data_s
        expected = DataFrame({0: [1, None], "new": [1e9, None]}, dtype="datetime64[ns]")
        tm.assert_frame_equal(result, expected)
