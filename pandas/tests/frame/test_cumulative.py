"""
Tests for DataFrame cumulative operations

See also
--------
tests.series.test_cumulative
"""

import numpy as np
import pytest

from pandas import (
    DataFrame,
    Series,
    Timestamp,
)
import pandas._testing as tm


class TestDataFrameCumulativeOps:
    # ---------------------------------------------------------------------
    # Cumulative Operations - cumsum, cummax, ...

    def test_cumulative_ops_smoke(self):
        # it works
        df = DataFrame({"A": np.arange(20)}, index=np.arange(20))
        df.cummax()
        df.cummin()
        df.cumsum()

        dm = DataFrame(np.arange(20).reshape(4, 5), index=range(4), columns=range(5))
        # TODO(wesm): do something with this?
        dm.cumsum()

    def test_cumprod_smoke(self, datetime_frame):
        datetime_frame.iloc[5:10, 0] = np.nan
        datetime_frame.iloc[10:15, 1] = np.nan
        datetime_frame.iloc[15:, 2] = np.nan

        # ints
        df = datetime_frame.fillna(0).astype(int)
        df.cumprod(0)
        df.cumprod(1)

        # ints32
        df = datetime_frame.fillna(0).astype(np.int32)
        df.cumprod(0)
        df.cumprod(1)

    def test_cumulative_ops_match_series_apply(
        self, datetime_frame, all_numeric_accumulations
    ):
        datetime_frame.iloc[5:10, 0] = np.nan
        datetime_frame.iloc[10:15, 1] = np.nan
        datetime_frame.iloc[15:, 2] = np.nan

        # axis = 0
        result = getattr(datetime_frame, all_numeric_accumulations)()
        expected = datetime_frame.apply(getattr(Series, all_numeric_accumulations))
        tm.assert_frame_equal(result, expected)

        # axis = 1
        result = getattr(datetime_frame, all_numeric_accumulations)(axis=1)
        expected = datetime_frame.apply(
            getattr(Series, all_numeric_accumulations), axis=1
        )
        tm.assert_frame_equal(result, expected)

        # fix issue TODO: GH ref?
        assert np.shape(result) == np.shape(datetime_frame)

    def test_cumsum_preserve_dtypes(self):
        # GH#19296 dont incorrectly upcast to object
        df = DataFrame({"A": [1, 2, 3], "B": [1, 2, 3.0], "C": [True, False, False]})

        result = df.cumsum()

        expected = DataFrame(
            {
                "A": Series([1, 3, 6], dtype=np.int64),
                "B": Series([1, 3, 6], dtype=np.float64),
                "C": df["C"].cumsum(),
            }
        )
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize("method", ["cumsum", "cumprod", "cummin", "cummax"])
    @pytest.mark.parametrize("axis", [0, 1])
    def test_numeric_only_flag(self, method, axis):
        df = DataFrame(
            {
                "int": [1, 2, 3],
                "bool": [True, False, False],
                "string": ["a", "b", "c"],
                "float": [1.0, 3.5, 4.0],
                "datetime": [
                    Timestamp(2018, 1, 1),
                    Timestamp(2019, 1, 1),
                    Timestamp(2020, 1, 1),
                ],
            }
        )
        df_numeric_only = df.drop(["string", "datetime"], axis=1)

        result = getattr(df, method)(axis=axis, numeric_only=True)
        expected = getattr(df_numeric_only, method)(axis)
        tm.assert_frame_equal(result, expected)
