import pytest

import pandas.util._test_decorators as td

import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize(
    "vals_left, vals_right, dtype",
    [
        ([1, 2, 3], [1, 2], "int64"),
        ([1.5, 2.5, 3.5], [1.5, 2.5], "float64"),
        ([True, True, False], [True, True], "bool"),
        (["a", "b", "c"], ["a", "b"], "object"),
        pytest.param(
            ["a", "b", "c"],
            ["a", "b"],
            "string[pyarrow]",
            marks=td.skip_if_no("pyarrow"),
        ),
        pytest.param(
            ["a", "b", "c"],
            ["a", "b"],
            "str",
            marks=td.skip_if_no("pyarrow"),
        ),
    ],
)
def test_left_semi(vals_left, vals_right, dtype):
    vals_left = pd.Series(vals_left, dtype=dtype)
    vals_right = pd.Series(vals_right, dtype=dtype)
    left = pd.DataFrame({"a": vals_left, "b": [1, 2, 3]})
    right = pd.DataFrame({"a": vals_right, "c": 1})
    expected = pd.DataFrame({"a": vals_right, "b": [1, 2]})
    result = left.merge(right, how="left_semi")
    tm.assert_frame_equal(result, expected)

    result = left.set_index("a").merge(
        right.set_index("a"), how="left_semi", left_index=True, right_index=True
    )
    tm.assert_frame_equal(result, expected.set_index("a"))

    result = left.set_index("a").merge(
        right, how="left_semi", left_index=True, right_on="a"
    )
    tm.assert_frame_equal(result, expected.set_index("a"))

    result = left.merge(
        right.set_index("a"), how="left_semi", right_index=True, left_on="a"
    )
    tm.assert_frame_equal(result, expected)

    right = pd.DataFrame({"d": vals_right, "c": 1})
    result = left.merge(right, how="left_semi", left_on="a", right_on="d")
    tm.assert_frame_equal(result, expected)

    right = pd.DataFrame({"d": vals_right, "c": 1})
    result = left.merge(right, how="left_semi", left_on=["a", "b"], right_on=["d", "c"])
    tm.assert_frame_equal(result, expected.head(1))


def test_left_semi_invalid():
    left = pd.DataFrame({"a": [1, 2, 3], "b": [1, 2, 3]})
    right = pd.DataFrame({"a": [1, 2], "c": 1})
    msg = "indicator is not supported for semi-join."
    with pytest.raises(NotImplementedError, match=msg):
        left.merge(right, how="left_semi", indicator=True)

    msg = "sort is not supported for semi-join. Sort your DataFrame afterwards."
    with pytest.raises(NotImplementedError, match=msg):
        left.merge(right, how="left_semi", sort=True)
