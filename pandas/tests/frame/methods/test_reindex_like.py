import numpy as np
import pytest

from pandas.errors import Pandas4Warning

import pandas as pd
import pandas._testing as tm


class TestDataFrameReindexLike:
    def test_reindex_like(self, float_frame):
        other = float_frame.reindex(index=float_frame.index[:10], columns=["C", "B"])

        tm.assert_frame_equal(other, float_frame.reindex_like(other))

    @pytest.mark.parametrize(
        "method,expected_values",
        [
            ("nearest", [0, 1, 1, 2]),
            ("pad", [np.nan, 0, 1, 1]),
            ("backfill", [0, 1, 2, 2]),
        ],
    )
    def test_reindex_like_methods(self, method, expected_values):
        df = pd.DataFrame({"x": list(range(5))})

        msg = "the 'method' keyword is deprecated"
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            result = df.reindex_like(df, method=method, tolerance=0)
        tm.assert_frame_equal(df, result)
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            result = df.reindex_like(df, method=method, tolerance=[0, 0, 0, 0])
        tm.assert_frame_equal(df, result)

    @pytest.mark.parametrize("method", ["backfill", "bfill", "pad", "ffill", "nearest"])
    def test_reindex_like_method_not_applied_to_columns(self, method):
        # GH 31002
        df = DataFrame([4, 5, 6], index=[0.5, 1.5, 2.5], columns=["b"])
        other = DataFrame([1, 2, 3], columns=["a"])

        with tm.assert_produces_warning(Pandas4Warning):
            result = df.reindex_like(other, method=method)

        expected = DataFrame({"a": [np.nan, np.nan, np.nan]}, index=other.index)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "method, kwargs",
        [("bfill", {"limit": 1}), ("nearest", {"tolerance": 1})],
    )
    def test_reindex_like_method_partial_column_overlap(self, method, kwargs):
        # GH 31002
        df = DataFrame({"a": [4, 5, 6], "c": [7, 8, 9]}, index=[0.5, 1.5, 2.5])
        other = DataFrame({"a": [1, 2, 3], "b": [1, 2, 3]}, index=[0, 1, 2])

        with tm.assert_produces_warning(Pandas4Warning):
            result = df.reindex_like(other, method=method, **kwargs)

        expected = DataFrame(
            {"a": [4, 5, 6], "b": [np.nan, np.nan, np.nan]}, index=[0, 1, 2]
        )
        tm.assert_frame_equal(result, expected)

    def test_reindex_like_subclass(self):
        # https://github.com/pandas-dev/pandas/issues/31925
        class MyDataFrame(pd.DataFrame):
            pass

        expected = pd.DataFrame()
        df = MyDataFrame()
        result = df.reindex_like(expected)

        tm.assert_frame_equal(result, expected)
