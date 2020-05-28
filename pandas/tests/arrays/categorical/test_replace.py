import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize(
    "to_replace,value,expected,flip_categories",
    [
        # one-to-one
        (1, 2, [2, 2, 3], False),
        (1, 4, [4, 2, 3], False),
        (4, 1, [1, 2, 3], False),
        (5, 6, [1, 2, 3], False),
        # many-to-one
        ([1], 2, [2, 2, 3], False),
        ([1, 2], 3, [3, 3, 3], False),
        ([1, 2], 4, [4, 4, 3], False),
        ((1, 2, 4), 5, [5, 5, 3], False),
        ((5, 6), 2, [1, 2, 3], False),
        # many-to-many, handled outside of Categorical and results in separate dtype
        ([1], [2], [2, 2, 3], True),
        ([1, 4], [5, 2], [5, 2, 3], True),
        # check_categorical sorts categories, which crashes on mixed dtypes
        (3, "4", [1, 2, "4"], False),
        ([1, 2, "3"], "5", ["5", "5", 3], True),
    ],
)
def test_replace(to_replace, value, expected, flip_categories):
    # GH 31720
    stays_categorical = not isinstance(value, list)

    s = pd.Series([1, 2, 3], dtype="category")
    result = s.replace(to_replace, value)
    expected = pd.Series(expected, dtype="category")
    s.replace(to_replace, value, inplace=True)

    if flip_categories:
        expected = expected.cat.set_categories(expected.cat.categories[::-1])

    if not stays_categorical:
        # the replace call loses categorical dtype
        expected = pd.Series(np.asarray(expected))

    tm.assert_series_equal(
        expected, result, check_category_order=False,
    )
    tm.assert_series_equal(
        expected, s, check_category_order=False,
    )
