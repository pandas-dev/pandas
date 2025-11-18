"""
Tests for core/methods/corr.py
"""

import pytest
import numpy as np
from pandas import DataFrame, Series, Categorical
import pandas._testing as tm
from pandas.core.methods.corr import transform_ord_cat_cols_to_coded_cols


@pytest.mark.parametrize(
    ("input_df", "expected_df"),
    [
        pytest.param(
            # 1) Simple: two ordered categorical columns (with and without None)
            DataFrame(
                {
                    "ord_cat": Series(
                        Categorical(
                            ["low", "m", "h", "vh"],
                            categories=["low", "m", "h", "vh"],
                            ordered=True,
                        )
                    ),
                    "ord_cat_none": Series(
                        Categorical(
                            ["low", "m", "h", None],
                            categories=["low", "m", "h"],
                            ordered=True,
                        )
                    ),
                }
            ),
            DataFrame(
                {
                    # codes: low=0, m=1, h=2, vh=3
                    "ord_cat": Series([0, 1, 2, 3], dtype="int8"),
                    # codes: low=0, m=1, h=2, None -> NaN
                    "ord_cat_none": Series([0, 1.0, 2.0, np.nan]),
                }
            ),
            id="ordered-categoricals-basic",
        ),
        pytest.param(
            # 2) Mixed dtypes: only the ordered categorical should change
            DataFrame(
                {
                    "ordered": Series(
                        Categorical(
                            ["a", "c", "b"],
                            categories=["a", "b", "c"],
                            ordered=True,
                        )
                    ),
                    "unordered": Series(Categorical(["x", "y", "x"], ordered=False)),
                    "num": Series([10, 20, 30]),
                    "text": Series(["u", "v", "w"]),
                }
            ),
            DataFrame(
                {
                    # codes: a=0, c=2, b=1
                    "ordered": Series([0, 2, 1], dtype="int8"),
                    # unordered categorical should be untouched (still categorical)
                    "unordered": Series(Categorical(["x", "y", "x"], ordered=False)),
                    "num": Series([10, 20, 30]),
                    "text": Series(["u", "v", "w"]),
                }
            ),
            id="mixed-types-only-ordered-changes",
        ),
        pytest.param(
            # 3 Duplicate column names: first 'dup' is ordered categorical,
            # second 'dup' is non-categorical
            DataFrame(
                {
                    "dup": Series(
                        Categorical(
                            ["low", "m", "h"],
                            categories=["low", "m", "h"],
                            ordered=True,
                        )
                    ),
                    "dup": Series([5, 6, 7]),  # duplicate name, later column
                }
            ),
            DataFrame(
                {
                    # After transform: position 0 (ordered cat) becomes codes [0,1,2],
                    # position 1 remains untouched numbers [5,6,7].
                    "dup": Series([0, 1, 2], dtype="int8"),
                    "dup": Series([5, 6, 7]),
                }
            ),
            id="duplicate-names-ordered-first",
        ),
        pytest.param(
            # 4 Duplicate column names: first 'dup' is non-categorical,
            # second 'dup' is ordered categorical, third 'dup' is ordered categorical
            DataFrame(
                {
                    "dup": Series(["a", "b", "c"]),  # non-categorical (object)
                    "dup": Series(
                        Categorical(
                            ["p", "q", None],
                            categories=["p", "q"],
                            ordered=True,
                        )
                    ),
                    "dup": Series(
                        Categorical(
                            ["low", "m", "h"],
                            categories=["low", "m", "h"],
                            ordered=True,
                        )
                    ),
                }
            ),
            DataFrame(
                {
                    # First stays object; second turns into codes [0, 1, NaN]
                    # and third changes into codes [0, 1, 2]
                    "dup": Series(["a", "b", "c"]),
                    "dup": Series([0.0, 1.0, np.nan]),
                    "dup": Series([0, 1, 2], dtype="int8"),
                }
            ),
            id="duplicate-names-ordered-and-non-categorical-and-none",
        ),
    ],
)
def test_transform_ord_cat_cols_to_coded_cols(input_df, expected_df):
    out_df = transform_ord_cat_cols_to_coded_cols(input_df)
    assert list(out_df.columns) == list(expected_df.columns)
    for i, col in enumerate(out_df.columns):
        tm.assert_series_equal(out_df.iloc[:, i], expected_df.iloc[:, i])
