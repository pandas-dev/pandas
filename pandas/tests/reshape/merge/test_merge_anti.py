import numpy as np
import pytest

from pandas import DataFrame
import pandas._testing as tm
from pandas.core.reshape.merge import merge


class Test_AntiJoin:
    @pytest.mark.parametrize(
        "how, expected",
        [
            (
                "anti_left",
                DataFrame({"A": 3, "C_x": 7, "B": np.nan, "C_y": np.nan}, index=["c"]),
            ),
            (
                "anti_right",
                DataFrame({"A": np.nan, "C_x": np.nan, "B": 4, "C_y": 9}, index=["d"]),
            ),
            (
                "anti_full",
                DataFrame(
                    {
                        "A": [3, np.nan],
                        "C_x": [7, np.nan],
                        "B": [np.nan, 4],
                        "C_y": [np.nan, 9],
                    },
                    index=["c", "d"],
                ),
            ),
        ],
    )
    def test_basic_anti_index(self, how, expected):
        left = DataFrame({"A": [1, 2, 3], "C": [5, 6, 7]}, index=["a", "b", "c"])
        right = DataFrame({"B": [1, 2, 4], "C": [7, 8, 9]}, index=["a", "b", "d"])
        result = merge(left, right, how=how, left_index=True, right_index=True)
        tm.assert_frame_equal(result, expected)

    @pytest.mark.parametrize(
        "on, how, expected",
        [
            (
                ["C"],
                "anti_left",
                DataFrame(
                    {"A": [1, 2], "C": [5, 6], "B": [np.nan, np.nan]}, index=[0, 1]
                ),
            ),
            (
                ["C"],
                "anti_right",
                DataFrame(
                    {"A": [np.nan, np.nan], "C": [8, 9], "B": [2, 4]}, index=[0, 1]
                ),
            ),
            (
                ["C"],
                "anti_full",
                DataFrame(
                    {
                        "A": [1, 2, np.nan, np.nan],
                        "C": [5, 6, 8, 9],
                        "B": [np.nan, np.nan, 2, 4],
                    },
                    index=[0, 1, 2, 3],
                ),
            ),
            (
                None,
                "anti_left",
                DataFrame({"A": [1, 2], "C": [5, 6], "B": [np.nan, np.nan]}),
            ),
            (
                None,
                "anti_right",
                DataFrame({"A": [np.nan, np.nan], "C": [8, 9], "B": [2, 4]}),
            ),
            (
                None,
                "anti_full",
                DataFrame(
                    {
                        "A": [1, 2, np.nan, np.nan],
                        "C": [5, 6, 8, 9],
                        "B": [np.nan, np.nan, 2, 4],
                    },
                ),
            ),
        ],
    )
    def test_basic_anti_on(self, on, how, expected):
        left = DataFrame({"A": [1, 2, 3], "C": [5, 6, 7]}, index=["a", "b", "c"])
        right = DataFrame({"B": [1, 2, 4], "C": [7, 8, 9]}, index=["a", "b", "d"])
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
        df_2 = DataFrame({"B": [5, 5, 9], "C": [4, 6, 7], "D": ["a", "b", "d"]})
        df_1 = DataFrame(
            {"A": [1, 2, 3], "B": [4, 5, 6], "C": [5, 6, 7]}, index=["a", "b", "c"]
        )
        result = merge(df_1, df_2, how=how, on=on)
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
        df_1 = DataFrame({"A": [np.nan, 2, np.nan], "B": ["a", 2, 3]})
        df_2 = DataFrame({"C": [1, 3, np.nan], "B": ["a", 2, "c"]})
        result = merge(df_1, df_2, on=on, how=how, left_on=left_on, right_on=right_on)
        tm.assert_frame_equal(result, expected)
