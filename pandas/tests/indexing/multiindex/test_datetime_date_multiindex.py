"""
Tests for GH#55969 - MultiIndex with datetime.date level
returns incorrect results when accessed with np.datetime64
"""

import datetime as dt

import numpy as np
import pytest

from pandas import (
    DataFrame,
    MultiIndex,
)
import pandas._testing as tm


@pytest.fixture
def mi_df_datetime_date():
    """DataFrame with 3-level MultiIndex where level 0 is datetime.date."""
    dates = [
        dt.date(2023, 11, 1),
        dt.date(2023, 11, 1),
        dt.date(2023, 11, 2),
    ]
    t1 = ["A", "B", "C"]
    t2 = ["C", "D", "E"]
    vals = [1.0, 2.0, 3.0]

    index = MultiIndex.from_arrays([dates, t1, t2], names=["dates", "t1", "t2"])
    return DataFrame({"vals": vals}, index=index)


class TestDatetimeDateMultiIndexLoc:
    """Tests for .loc with datetime.date MultiIndex and np.datetime64 key."""

    def test_loc_np_datetime64_partial_key_a(self, mi_df_datetime_date):
        # GH#55969
        df = mi_df_datetime_date
        result = df.loc[(np.datetime64("2023-11-01"), "A")]
        expected = DataFrame(
            {"vals": [1.0]},
            index=MultiIndex.from_arrays([["C"]], names=["t2"]),
        )
        tm.assert_frame_equal(result, expected)

    def test_loc_np_datetime64_partial_key_b(self, mi_df_datetime_date):
        # GH#55969
        df = mi_df_datetime_date
        result = df.loc[(np.datetime64("2023-11-01"), "B")]
        expected = DataFrame(
            {"vals": [2.0]},
            index=MultiIndex.from_arrays([["D"]], names=["t2"]),
        )
        tm.assert_frame_equal(result, expected)

    def test_loc_np_datetime64_nonexistent_second_level_raises(
        self, mi_df_datetime_date
    ):
        # GH#55969
        df = mi_df_datetime_date
        with pytest.raises(KeyError, match="C"):
            df.loc[(np.datetime64("2023-11-01"), "C")]

    def test_loc_np_datetime64_matches_native_date(self, mi_df_datetime_date):
        # GH#55969
        df = mi_df_datetime_date
        result_np = df.loc[(np.datetime64("2023-11-01"), "A")]
        result_native = df.loc[(dt.date(2023, 11, 1), "A")]
        tm.assert_frame_equal(result_np, result_native)

    def test_loc_np_datetime64_with_slice_none(self, mi_df_datetime_date):
        # GH#55969 - the known workaround should still work
        df = mi_df_datetime_date
        result = df.loc[(np.datetime64("2023-11-01"), "A", slice(None))]
        assert len(result) == 1
        assert result["vals"].iloc[0] == 1.0

    def test_loc_np_datetime64_full_key(self, mi_df_datetime_date):
        # GH#55969
        df = mi_df_datetime_date
        result = df.loc[(np.datetime64("2023-11-01"), "A", "C")]
        expected = df.loc[(dt.date(2023, 11, 1), "A", "C")]
        tm.assert_series_equal(result, expected)

    def test_loc_np_datetime64_full_key_nonexistent_raises(self, mi_df_datetime_date):
        # GH#55969
        df = mi_df_datetime_date
        with pytest.raises(KeyError):
            df.loc[(np.datetime64("2023-11-01"), "A", "Z")]

    def test_loc_np_datetime64_stability_across_random_values(self):
        # GH#55969 - the original bug produced different results depending
        # on random values in the DataFrame due to hash/comparison issues
        for seed in range(20):
            rng = np.random.default_rng(seed)
            dates = [
                dt.date(2023, 11, 1),
                dt.date(2023, 11, 1),
                dt.date(2023, 11, 2),
            ]
            vals = rng.uniform(size=3)
            df = DataFrame(
                data=np.array([dates, ["A", "B", "C"], ["C", "D", "E"], vals]).T,
                columns=["dates", "t1", "t2", "vals"],
            )
            df.set_index(["dates", "t1", "t2"], inplace=True)

            date = np.datetime64("2023-11-01")

            result_a = df.loc[(date, "A")]
            assert len(result_a) == 1, f"seed={seed}: (date, 'A') got {len(result_a)}"

            result_b = df.loc[(date, "B")]
            assert len(result_b) == 1, f"seed={seed}: (date, 'B') got {len(result_b)}"

            with pytest.raises(KeyError):
                df.loc[(date, "C")]
