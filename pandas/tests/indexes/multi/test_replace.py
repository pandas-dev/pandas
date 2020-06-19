import pytest

import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize(
    "names, arrays, to_replace, value, expected_arrays",
    [
        (
            [None, None],
            [[1, 1, 2, 2], ["red", "blue", "red", "blue"]],
            [1, "red"],
            [0, "black"],
            [[0, 0, 2, 2], ["black", "blue", "black", "blue"]],
        ),
        # names should be preserved
        (
            ["digits", "colors"],
            [[1, 1, 2, 2], ["red", "blue", "red", "blue"]],
            1,
            0,
            [[0, 0, 2, 2], ["red", "blue", "red", "blue"]],
        ),
        (
            [None, None],
            [[1, 1, 2, 2], ["red", "blue", "red", "blue"]],
            1,
            0,
            [[0, 0, 2, 2], ["red", "blue", "red", "blue"]],
        ),
        (
            [None, None],
            [[1, 1, 2, 2], ["red", "blue", "red", "blue"]],
            [1, 2],
            0,
            [[0, 0, 0, 0], ["red", "blue", "red", "blue"]],
        ),
        (
            [None, None],
            [[1, 1, 2, 2], ["red", "blue", "red", "blue"]],
            [1, 2],
            0,
            [[0, 0, 0, 0], ["red", "blue", "red", "blue"]],
        ),
        # nested dicts
        (
            ["digits", "colors"],
            [[1, 1, 2, 2], ["red", "blue", "red", "blue"]],
            {"digits": {1: 0}, "colors": {"red": "black"}},
            None,
            [[0, 0, 2, 2], ["black", "blue", "black", "blue"]],
        ),
        # dicts and value
        (
            ["digits", "colors"],
            [[1, 1, 2, 2], ["red", "blue", "red", "blue"]],
            {"digits": [1], "colors": ["red", "blue"]},
            "x",
            [["x", "x", 2, 2], ["x", "x", "x", "x"]],
        ),
    ],
)
def test_multi_index_replace(names, arrays, to_replace, value, expected_arrays):
    multi_index = pd.MultiIndex.from_arrays(arrays, names=names)
    expected = pd.MultiIndex.from_arrays(expected_arrays, names=names)

    result = multi_index.replace(to_replace=to_replace, value=value)

    tm.assert_equal(result, expected)
