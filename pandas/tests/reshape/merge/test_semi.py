import pytest

import pandas.util._test_decorators as td

import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize(
    "vals_left, vals_right, dtype",
    [
        ([1, 2, 3], [1, 2], "int64"),
        (["a", "b", "c"], ["a", "b"], "object"),
        pytest.param(
            ["a", "b", "c"],
            ["a", "b"],
            "string[pyarrow]",
            marks=td.skip_if_no("pyarrow"),
        ),
    ],
)
def test_leftsemi(vals_left, vals_right, dtype):
    vals_left = pd.Series(vals_left, dtype=dtype)
    vals_right = pd.Series(vals_right, dtype=dtype)
    left = pd.DataFrame({"a": vals_left, "b": [1, 2, 3]})
    right = pd.DataFrame({"a": vals_right, "c": 1})
    expected = pd.DataFrame({"a": vals_right, "b": [1, 2]})
    result = left.merge(right, how="leftsemi")
    tm.assert_frame_equal(result, expected)

    right = pd.DataFrame({"d": vals_right, "c": 1})
    result = left.merge(right, how="leftsemi", left_on="a", right_on="d")
    tm.assert_frame_equal(result, expected)

    right = pd.DataFrame({"d": vals_right, "c": 1})
    result = left.merge(right, how="leftsemi", left_on=["a", "b"], right_on=["d", "c"])
    tm.assert_frame_equal(result, expected.head(1))


def test_leftsemi_invalid():
    left = pd.DataFrame({"a": [1, 2, 3], "b": [1, 2, 3]})
    right = pd.DataFrame({"a": [1, 2], "c": 1})

    msg = "left_index or right_index are not supported for semi-join."
    with pytest.raises(NotImplementedError, match=msg):
        left.merge(right, how="leftsemi", left_index=True, right_on="a")
    with pytest.raises(NotImplementedError, match=msg):
        left.merge(right, how="leftsemi", right_index=True, left_on="a")

    msg = "validate is not supported for semi-join."
    with pytest.raises(NotImplementedError, match=msg):
        left.merge(right, how="leftsemi", validate="one_to_one")

    msg = "indicator is not supported for semi-join."
    with pytest.raises(NotImplementedError, match=msg):
        left.merge(right, how="leftsemi", indicator=True)

    msg = "sort is not supported for semi-join. Sort your DataFrame afterwards."
    with pytest.raises(NotImplementedError, match=msg):
        left.merge(right, how="leftsemi", sort=True)
