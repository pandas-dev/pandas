import pytest

import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize(
    "index, to_replace, value, expected",
    [
        ([1, 2, 3], 3, "a", [1, 2, "a"]),
        (
            [1, None, 2],
            [1, 2],
            "a",
            ["a", None, "a"],
        ),
    ],
)
def test_categorical_index_replace(index, to_replace, value, expected):
    index = pd.CategoricalIndex(index)
    expected = pd.CategoricalIndex(expected)

    result = index.replace(to_replace=to_replace, value=value)

    tm.assert_equal(result, expected)


def test_categorical_index_replace_dict_and_value():
    index = pd.CategoricalIndex([1, 2, 3])

    msg = "Series.replace cannot use dict-like to_replace and non-None value"
    with pytest.raises(ValueError, match=msg):
        index.replace({1: "a", 3: "c"}, "x")


@pytest.mark.parametrize(
    "index, to_replace, value, expected",
    [
        ([1, 2, 3], [2, 3], ["b", "c"], [1, "b", "c"]),
        ([1, 2, 3], 3, "c", [1, 2, "c"]),
        (
            [1, None, 2],
            [1, 2],
            "a",
            ["a", None, "a"],
        ),
    ],
)
def test_index_replace(index, to_replace, value, expected):
    index = pd.CategoricalIndex(index)
    expected = pd.CategoricalIndex(expected)

    result = index.replace(to_replace=to_replace, value=value)

    tm.assert_equal(result, expected)
