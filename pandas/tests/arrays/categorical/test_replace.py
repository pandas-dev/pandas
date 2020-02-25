import pytest

import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize(
    "to_replace,value,expected,check_types,check_categorical",
    [
        # one-to-one
        (1, 2, [2, 2, 3], True, True),
        (1, 4, [4, 2, 3], True, True),
        (4, 1, [1, 2, 3], True, True),
        (5, 6, [1, 2, 3], True, True),
        # many-to-one
        ([1], 2, [2, 2, 3], True, True),
        ([1, 2], 3, [3, 3, 3], True, True),
        ([1, 2], 4, [4, 4, 3], True, True),
        ((1, 2, 4), 5, [5, 5, 3], True, True),
        ((5, 6), 2, [1, 2, 3], True, True),
        # many-to-many, handled outside of Categorical and results in separate dtype
        ([1], [2], [2, 2, 3], False, False),
        ([1, 4], [5, 2], [5, 2, 3], False, False),
        # check_categorical sorts categories, which crashes on mixed dtypes
        (3, "4", [1, 2, "4"], True, False),
        ([1, 2, "3"], "5", ["5", "5", 3], True, False),
    ],
)
def test_replace(to_replace, value, expected, check_types, check_categorical):
    # GH 31720
    s = pd.Series([1, 2, 3], dtype="category")
    result = s.replace(to_replace, value)
    expected = pd.Series(expected, dtype="category")
    s.replace(to_replace, value, inplace=True)
    tm.assert_series_equal(
        expected,
        result,
        check_dtype=check_types,
        check_categorical=check_categorical,
        check_category_order=False,
    )
    tm.assert_series_equal(
        expected,
        s,
        check_dtype=check_types,
        check_categorical=check_categorical,
        check_category_order=False,
    )
