"""
Tests for core/methods/corr.py
"""

import numpy as np
import pytest

from pandas import (
    Categorical,
    DataFrame,
    Series,
)
import pandas._testing as tm
from pandas.core.methods.corr import transform_ord_cat_cols_to_coded_cols


@pytest.mark.parametrize(
    ("input_df_dict", "expected_df_dict"),
    [
        pytest.param(
            # 1) Simple: two ordered categorical columns (with and without None)
            {
                "ord_cat": Categorical(
                    ["low", "m", "h", "vh"],
                    categories=["low", "m", "h", "vh"],
                    ordered=True,
                ),
                "ord_cat_none": Categorical(
                    ["low", "m", "h", None],
                    categories=["low", "m", "h"],
                    ordered=True,
                ),
            },
            {
                # codes: low=0, m=1, h=2, vh=3
                "ord_cat": Series([0, 1, 2, 3], dtype="int8"),
                # codes: low=0, m=1, h=2, None -> NaN
                "ord_cat_none": [0, 1.0, 2.0, np.nan],
            },
            id="ordered-categoricals-basic",
        ),
        pytest.param(
            # 2) Mixed dtypes: only the ordered categorical should change
            {
                "ordered": Categorical(
                    ["a", "c", "b"],
                    categories=["a", "b", "c"],
                    ordered=True,
                ),
                "unordered": Categorical(["x", "y", "x"], ordered=False),
                "num": [10, 20, 30],
                "text": ["u", "v", "w"],
            },
            {
                # codes: a=0, c=2, b=1
                "ordered": Series([0, 2, 1], dtype="int8"),
                # unordered categorical should be untouched (still categorical)
                "unordered": Categorical(["x", "y", "x"], ordered=False),
                "num": [10, 20, 30],
                "text": ["u", "v", "w"],
            },
            id="mixed-types-only-ordered-changes",
        ),
    ],
)
def test_transform_ord_cat_cols_to_coded_cols(
    input_df_dict: dict, expected_df_dict: dict
) -> None:
    # GH #60306
    input_df = DataFrame(input_df_dict)
    expected_df = DataFrame(expected_df_dict)
    out_df = transform_ord_cat_cols_to_coded_cols(input_df)
    assert list(out_df.columns) == list(expected_df.columns)
    tm.assert_frame_equal(out_df, expected_df)


@pytest.mark.parametrize(
    ("input_df_dict", "expected_df_dict"),
    [
        pytest.param(
            {
                "dup_1": Categorical(
                    ["low", "m", "h"],
                    categories=["low", "m", "h"],
                    ordered=True,
                ),
                "dup_2": [5, 6, 7],
            },
            {
                # After transform: position 0 (ordered cat) becomes codes [0,1,2],
                # position 1 remains untouched numbers [5,6,7].
                "dup_1": Series([0, 1, 2], dtype="int8"),
                "dup_2": [5, 6, 7],
            },
            id="duplicate-names-ordered-first",
        ),
        pytest.param(
            {
                "dup_1": ["a", "b", "c"],  # non-categorical
                "dup_2": Categorical(
                    ["p", "q", None],
                    categories=["p", "q"],
                    ordered=True,
                ),
                "dup_3": Categorical(
                    ["low", "m", "h"],
                    categories=["low", "m", "h"],
                    ordered=True,
                ),
            },
            {
                # First stays object; second turns into codes [0, 1, NaN]
                # and third changes into codes [0, 1, 2]
                "dup_1": ["a", "b", "c"],
                "dup_2": [0.0, 1.0, np.nan],
                "dup_3": Series([0, 1, 2], dtype="int8"),
            },
            id="duplicate-names-ordered-and-non-categorical-and-none",
        ),
    ],
)
def test_transform_ord_cat_cols_to_coded_cols_duplicated_col(
    input_df_dict: dict, expected_df_dict: dict
) -> None:
    # GH #60306
    input_df = DataFrame(input_df_dict)
    expected_df = DataFrame(expected_df_dict)
    input_df.columns = ["dup" for _ in input_df.columns]
    expected_df.columns = ["dup" for _ in expected_df.columns]

    out_df = transform_ord_cat_cols_to_coded_cols(input_df)
    tm.assert_frame_equal(out_df, expected_df)
