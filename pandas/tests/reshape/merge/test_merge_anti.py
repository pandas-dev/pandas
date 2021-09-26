import numpy as np
import pytest

import pandas as pd
from pandas import (
    Categorical,
    DataFrame,
)
import pandas._testing as tm
from pandas.core.reshape.merge import merge


class Test_AntiJoin:
    @pytest.mark.parametrize(
        "how, exp_index, exp_values",
        [
            ("anti_left", ["c"], [[3, 30, np.nan, np.nan]]),
            ("anti_right", ["d"], [[np.nan, np.nan, 4, 40]]),
            (
                "anti_full",
                ["c", "d"],
                [[3, 30, np.nan, np.nan], [np.nan, np.nan, 4, 40]],
            ),
        ],
    )
    def test_basic_anti_index(self, how, exp_index, exp_values):
        # basic test containing NaNs w/o on param
        left = DataFrame({"A": [1, 2, 3], "C": [10, 20, 30]}, index=["a", "b", "c"])
        right = DataFrame({"B": [1, 2, 4], "C": [10, 20, 40]}, index=["a", "b", "d"])
        expected = DataFrame(
            exp_values, index=exp_index, columns=["A", "C_x", "B", "C_y"]
        )
        result = merge(left, right, how=how, left_index=True, right_index=True)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "on, how, data",
        [
            (
                ["C"],
                "anti_left",
                [[1, 5, np.nan], [2, 6, np.nan]],
            ),
            (
                ["C"],
                "anti_right",
                [[np.nan, 8, 2], [np.nan, 9, 4]],
            ),
            (
                ["C"],
                "anti_full",
                [[1, 5, np.nan], [2, 6, np.nan], [np.nan, 8, 2], [np.nan, 9, 4]],
            ),
            (
                None,
                "anti_left",
                [[1, 5, np.nan], [2, 6, np.nan]],
            ),
            (
                None,
                "anti_right",
                [[np.nan, 8, 2], [np.nan, 9, 4]],
            ),
            (
                None,
                "anti_full",
                [[1, 5, np.nan], [2, 6, np.nan], [np.nan, 8, 2], [np.nan, 9, 4]],
            ),
        ],
    )
    def test_basic_anti_on(self, on, how, data):
        # basic test containing NaNs with on param
        left = DataFrame({"A": [1, 2, 3], "C": [5, 6, 7]}, index=["a", "b", "c"])
        right = DataFrame({"B": [1, 2, 4], "C": [7, 8, 9]}, index=["a", "b", "d"])
        expected = DataFrame(data, columns=["A", "C", "B"])
        result = merge(left, right, how=how, on=on)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "expected, how, left_on, right_on",
        [
            (
                DataFrame({"A": 3, "C_x": 7, "B": np.nan, "C_y": np.nan}, index=[0]),
                "anti_left",
                ["A"],
                ["B"],
            ),
            (
                DataFrame({"A": np.nan, "C_x": np.nan, "B": 4, "C_y": 9}, index=[0]),
                "anti_right",
                ["A"],
                ["B"],
            ),
            (
                DataFrame(
                    {
                        "A": [3, np.nan],
                        "C_x": [7, np.nan],
                        "B": [np.nan, 4],
                        "C_y": [np.nan, 9],
                    },
                ),
                "anti_full",
                ["A"],
                ["B"],
            ),
        ],
    )
    def test_basic_anti_lefton_righton(self, expected, how, left_on, right_on):
        # basic test containing NaNs with left_on / right_on params
        left = DataFrame({"A": [1, 2, 3], "C": [5, 6, 7]}, index=["a", "b", "c"])
        right = DataFrame({"B": [1, 2, 4], "C": [7, 8, 9]}, index=["a", "b", "d"])
        result = merge(
            left,
            right,
            how=how,
            left_on=left_on,
            right_on=right_on,
        )
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "expected, how",
        [
            (
                DataFrame(
                    {
                        "A": [3],
                        "B_x": [6],
                        "C_x": [7],
                        "B_y": [np.nan],
                        "C_y": [np.nan],
                        "D": ["c"],
                    },
                    index=[np.nan],
                ),
                "anti_left",
            ),
            (
                DataFrame(
                    {
                        "A": [np.nan],
                        "B_x": [np.nan],
                        "C_x": [np.nan],
                        "B_y": [9],
                        "C_y": [7],
                        "D": ["d"],
                    },
                    index=[2],
                ),
                "anti_right",
            ),
            (
                DataFrame(
                    {
                        "A": [3, np.nan],
                        "B_x": [6, np.nan],
                        "C_x": [7, np.nan],
                        "B_y": [np.nan, 9],
                        "C_y": [np.nan, 7],
                        "D": ["c", "d"],
                    },
                    index=[np.nan, 2],
                ),
                "anti_full",
            ),
        ],
    )
    def test_anti_index_with_col(self, expected, how):
        # basic test containing NaNs with left_index and right_on params
        left = DataFrame(
            {"A": [1, 2, 3], "B": [4, 5, 6], "C": [5, 6, 7]}, index=["a", "b", "c"]
        )
        right = DataFrame({"B": [5, 5, 9], "C": [4, 6, 7], "D": ["a", "b", "d"]})
        result = merge(left, right, how=how, left_index=True, right_on=["D"])
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "expected, how, on",
        [
            (
                DataFrame(
                    {"A": [1, 3], "B": [4, 6], "C": [5, 7], "D": [np.nan, np.nan]}
                ).astype({"D": object}),
                "anti_left",
                None,
            ),
            (
                DataFrame(
                    {"A": [np.nan, np.nan], "B": [5, 9], "C": [4, 7], "D": ["a", "d"]}
                ),
                "anti_right",
                None,
            ),
            (
                DataFrame(
                    {
                        "A": [1, 3, np.nan, np.nan],
                        "B": [4, 6, 5, 9],
                        "C": [5, 7, 4, 7],
                        "D": [np.nan, np.nan, "a", "d"],
                    }
                ),
                "anti_full",
                None,
            ),
            (
                DataFrame(
                    {"A": [1, 3], "B": [4, 6], "C": [5, 7], "D": [np.nan, np.nan]}
                ).astype({"D": object}),
                "anti_left",
                ["B", "C"],
            ),
            (
                DataFrame(
                    {"A": [np.nan, np.nan], "B": [5, 9], "C": [4, 7], "D": ["a", "d"]}
                ),
                "anti_right",
                ["B", "C"],
            ),
            (
                DataFrame(
                    {
                        "A": [1, 3, np.nan, np.nan],
                        "B": [4, 6, 5, 9],
                        "C": [5, 7, 4, 7],
                        "D": [np.nan, np.nan, "a", "d"],
                    }
                ),
                "anti_full",
                ["B", "C"],
            ),
        ],
    )
    def test_anti_multicol(self, expected, how, on):
        # test with multicol with and w/o on param
        right = DataFrame({"B": [5, 5, 9], "C": [4, 6, 7], "D": ["a", "b", "d"]})
        left = DataFrame(
            {"A": [1, 2, 3], "B": [4, 5, 6], "C": [5, 6, 7]}, index=["a", "b", "c"]
        )
        result = merge(left, right, how=how, on=on)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "expected, how, on, left_on, right_on",
        [
            (
                DataFrame({"A": np.nan, "B": "c", "C": np.nan}, index=[0]),
                "anti_right",
                None,
                None,
                None,
            ),
            (
                DataFrame({"A": np.nan, "B": 3, "C": np.nan}, index=[0]).astype(
                    {"B": object}
                ),
                "anti_left",
                None,
                None,
                None,
            ),
            (
                DataFrame({"A": np.nan, "B": "c", "C": np.nan}, index=[0]),
                "anti_right",
                ["B"],
                None,
                None,
            ),
            (
                DataFrame({"A": np.nan, "B": 3, "C": np.nan}, index=[0]).astype(
                    {"B": object}
                ),
                "anti_left",
                ["B"],
                None,
                None,
            ),
            (
                DataFrame(
                    {"A": [2.0], "B_x": [2], "C": [np.nan], "B_y": [np.nan]}
                ).astype({"B_x": object, "B_y": object}),
                "anti_left",
                None,
                ["A"],
                ["C"],
            ),
            (
                DataFrame(
                    {
                        "A": [np.nan, np.nan],
                        "B_x": [np.nan, np.nan],
                        "C": [1.0, 3.0],
                        "B_y": ["a", 2],
                    }
                ).astype({"B_x": object}),
                "anti_right",
                None,
                ["A"],
                ["C"],
            ),
        ],
    )
    def test_anti_with_nan(self, expected, how, on, left_on, right_on):
        # basic anti_joins with mixed dtypes
        left = DataFrame({"A": [np.nan, 2, np.nan], "B": ["a", 2, 3]})
        right = DataFrame({"C": [1, 3, np.nan], "B": ["a", 2, "c"]})
        result = merge(left, right, on=on, how=how, left_on=left_on, right_on=right_on)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "expected, how, left_on, right_on",
        [
            (
                DataFrame({"A": [np.nan, pd.NA], "B": ["a", 3], "C": [np.nan, np.nan]}),
                "anti_left",
                "B",
                "B",
            ),
            (
                DataFrame(
                    {"A": [np.nan, np.nan], "B": [pd.NA, "c"], "C": [1, np.nan]}
                ).astype({"A": object}),
                "anti_right",
                "B",
                "B",
            ),
            (
                DataFrame(
                    {
                        "A": [np.nan, pd.NA, np.nan, np.nan],
                        "B": ["a", 3, pd.NA, "c"],
                        "C": [np.nan, np.nan, 1, np.nan],
                    }
                ),
                "anti_full",
                "B",
                "B",
            ),
            (
                DataFrame(
                    {
                        "A": [2, pd.NA],
                        "B_x": [2, 3],
                        "C": [np.nan, np.nan],
                        "B_y": [np.nan, np.nan],
                    }
                ).astype({"B_x": object, "B_y": object}),
                "anti_left",
                "A",
                "C",
            ),
            (
                DataFrame(
                    {
                        "A": [np.nan, np.nan],
                        "B_x": [np.nan, np.nan],
                        "C": [1.0, 3],
                        "B_y": [pd.NA, 2],
                    }
                ).astype(
                    {
                        "A": object,
                        "B_x": object,
                    }
                ),
                "anti_right",
                "A",
                "C",
            ),
            (
                DataFrame(
                    {
                        "A": [2, pd.NA, np.nan, np.nan],
                        "B_x": [2, 3, np.nan, np.nan],
                        "C": [np.nan, np.nan, 1, 3],
                        "B_y": [np.nan, np.nan, pd.NA, 2],
                    }
                ).astype({"B_x": object}),
                "anti_full",
                "A",
                "C",
            ),
        ],
    )
    def test_anti_with_nan_and_NA(self, expected, how, left_on, right_on):
        # test to check np.nan isn't matched with pd.NA
        left = DataFrame({"A": [np.nan, 2, pd.NA], "B": ["a", 2, 3]})
        right = DataFrame({"C": [1, 3, np.nan], "B": [pd.NA, 2, "c"]})
        result = merge(left, right, how=how, left_on=left_on, right_on=right_on)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "how, expected",
        [
            (
                "anti_left",
                DataFrame(
                    {"vals_x": [20, 17], "vals_y": [np.nan] * 2},
                    index=pd.date_range("1/2/2010", periods=2, freq="2d"),
                ),
            ),
            (
                "anti_right",
                DataFrame(
                    {"vals_x": [np.nan] * 2, "vals_y": [17, 21]},
                    index=pd.date_range("1/7/2010", periods=2, freq="2d"),
                ),
            ),
            (
                "anti_full",
                DataFrame(
                    {
                        "vals_x": [20, 17, np.nan, np.nan],
                        "vals_y": [np.nan, np.nan, 17, 21],
                    },
                    index=pd.date_range("1/7/2010", periods=2, freq="2d").union(
                        pd.date_range("1/2/2010", periods=2, freq="2d")
                    ),
                ),
            ),
        ],
    )
    def test_anti_datetime(self, how, expected):
        left = DataFrame(
            {"vals": [10, 20, 15, 17, 21]},
            index=pd.date_range("1/1/2010", periods=5, freq="D"),
        )
        right = DataFrame(
            {"vals": [10, 20, 15, 17, 21]},
            index=pd.date_range("1/1/2010", periods=5, freq="2D"),
        )
        result = merge(left, right, left_index=True, right_index=True, how=how)
        tm.assert_frame_equal(result, expected)

    def test_anti_datetime_tz(self):
        expected = DataFrame(
            {
                "Date": pd.date_range(
                    "10-20-2021", periods=2, freq="6D", tz="Asia/Kolkata"
                ),
                "a_x": [3, 4],
                "a_y": [np.nan, np.nan],
            }
        )
        left = DataFrame(
            {
                "Date": pd.date_range(
                    "10-02-2021", periods=5, freq="6D", tz="Asia/Kolkata"
                ),
                "a": range(5),
            }
        )
        right = DataFrame(
            {
                "Date": pd.date_range(
                    "10-02-2021", periods=5, freq="3D", tz="Asia/Kolkata"
                ),
                "a": range(5),
            }
        )
        result = merge(left, right, how="anti_left", on="Date")
        tm.assert_frame_equal(result, expected)

    def test_anti_categorical(self):
        left = DataFrame({"A": list("abca"), "B": list("bccd")}, dtype="category")
        right = DataFrame({"A": list("dad"), "C": list("gap")}, dtype="category")
        expected = DataFrame(
            {
                "A": ["b", "c"],
                "B": Categorical(["c", "c"], categories=list("bcd")),
                "C": Categorical([np.nan, np.nan], categories=list("agp")),
            }
        )
        result = merge(left, right, how="anti_left")
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "dtype", ["Int64", "Int32", "UInt32", "UInt64", "Float32", "Float64"]
    )
    def test_anti_EA_dtypes(self, dtype):
        left = DataFrame({"A": [1, 2, 3], "B": [4, 5, 6]}, dtype=dtype)
        right = DataFrame({"A": [1, 4, 5], "C": [7, 6, 8]}, dtype=dtype)
        result = merge(left, right, how="anti_right")
        expected = DataFrame(
            {"A": [4, 5], "B": [pd.NA, pd.NA], "C": [6, 8]}, dtype=dtype
        )
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "dtype", ["Int64", "Int32", "UInt32", "UInt64", "Float32", "Float64"]
    )
    def test_anti_EA_dtypes_with_multicol(self, dtype):
        left = DataFrame(
            {"A": [1, 2, 3], "B": [4, 5, 6], "C": [5, 6, 7]},
            index=["a", "b", "c"],
            dtype=dtype,
        )
        right = DataFrame({"B": [5, 5, 9], "C": [4, 6, 7], "D": [1, 0, 0]}, dtype=dtype)
        expected = DataFrame(
            columns=list("ABCD"), data=[[1, 4, 5, pd.NA], [3, 6, 7, pd.NA]], dtype=dtype
        )
        result = merge(left, right, how="anti_left")
        tm.assert_frame_equal(result, expected)
