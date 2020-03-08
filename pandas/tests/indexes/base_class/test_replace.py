import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize(
    "index, to_replace, value, expected",
    [
        ([1, 2, 3], [1, 3], ["a", "c"], ["a", 2, "c"]),
        ([1, 2, 3], 1, "a", ["a", 2, 3]),
        ([1, None, 2], [1, 2], "a", ["a", None, "a"],),
    ],
)
def test_index_replace(index, to_replace, value, expected):
    index = pd.Index(index)
    expected = pd.Index(expected)

    result = index.replace(to_replace=to_replace, value=value)

    tm.assert_equal(result, expected)


@pytest.mark.parametrize(
    "index, to_replace, value, regex, expected",
    [
        (
            ["bat", "foo", "baait", "bar"],
            r"^ba.$",
            "new",
            True,
            ["new", "foo", "baait", "new"],
        ),
        (
            ["bat", "foo", "baait", "bar"],
            None,
            None,
            {r"^ba.$": "new", "foo": "xyz"},
            ["new", "xyz", "baait", "new"],
        ),
    ],
)
def test_index_replace_regex(index, to_replace, value, regex, expected):
    index = pd.Index(index)
    expected = pd.Index(expected)

    result = index.replace(to_replace=to_replace, value=value, regex=regex)
    tm.assert_equal(expected, result)


def test_index_replace_dict_and_value():
    index = pd.Index([1, 2, 3])

    msg = "If 'to_replace' is a dict, 'value' should be None."
    with pytest.raises(TypeError, match=msg):
        index.replace({1: "a", 3: "c"}, "x")


def test_index_replace_scalar_only():
    index = pd.Index([1, None, 2])
    with pytest.raises(NotImplementedError):
        index.replace(np.nan)


def test_index_replace_bfill():
    index = pd.Index([0, 1, 2, 3, 4])
    expected = pd.Index([0, 3, 3, 3, 4])

    result = index.replace([1, 2], method="bfill")
    tm.assert_equal(expected, result)
