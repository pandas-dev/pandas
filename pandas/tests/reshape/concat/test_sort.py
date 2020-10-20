import pytest

import pandas as pd
import pandas._testing as tm


@pytest.fixture(params=[True, False])
def sort(request):
    """Boolean sort keyword for concat and DataFrame.append."""
    return request.param


def test_concat_sorts_columns(sort):
    # GH-4588
    df1 = pd.DataFrame({"a": [1, 2], "b": [1, 2]}, columns=["b", "a"])
    df2 = pd.DataFrame({"a": [3, 4], "c": [5, 6]})

    # for sort=True/None
    expected = pd.DataFrame(
        {"a": [1, 2, 3, 4], "b": [1, 2, None, None], "c": [None, None, 5, 6]},
        columns=["a", "b", "c"],
    )

    if sort is False:
        expected = expected[["b", "a", "c"]]

    # default
    with tm.assert_produces_warning(None):
        result = pd.concat([df1, df2], ignore_index=True, sort=sort)
    tm.assert_frame_equal(result, expected)


def test_concat_sorts_index(sort):
    df1 = pd.DataFrame({"a": [1, 2, 3]}, index=["c", "a", "b"])
    df2 = pd.DataFrame({"b": [1, 2]}, index=["a", "b"])

    # For True/None
    expected = pd.DataFrame(
        {"a": [2, 3, 1], "b": [1, 2, None]}, index=["a", "b", "c"], columns=["a", "b"]
    )
    if sort is False:
        expected = expected.loc[["c", "a", "b"]]

    # Warn and sort by default
    with tm.assert_produces_warning(None):
        result = pd.concat([df1, df2], axis=1, sort=sort)
    tm.assert_frame_equal(result, expected)


def test_concat_inner_sort(sort):
    # https://github.com/pandas-dev/pandas/pull/20613
    df1 = pd.DataFrame({"a": [1, 2], "b": [1, 2], "c": [1, 2]}, columns=["b", "a", "c"])
    df2 = pd.DataFrame({"a": [1, 2], "b": [3, 4]}, index=[3, 4])

    with tm.assert_produces_warning(None):
        # unset sort should *not* warn for inner join
        # since that never sorted
        result = pd.concat([df1, df2], sort=sort, join="inner", ignore_index=True)

    expected = pd.DataFrame({"b": [1, 2, 3, 4], "a": [1, 2, 1, 2]}, columns=["b", "a"])
    if sort is True:
        expected = expected[["a", "b"]]
    tm.assert_frame_equal(result, expected)


def test_concat_aligned_sort():
    # GH-4588
    df = pd.DataFrame({"c": [1, 2], "b": [3, 4], "a": [5, 6]}, columns=["c", "b", "a"])
    result = pd.concat([df, df], sort=True, ignore_index=True)
    expected = pd.DataFrame(
        {"a": [5, 6, 5, 6], "b": [3, 4, 3, 4], "c": [1, 2, 1, 2]},
        columns=["a", "b", "c"],
    )
    tm.assert_frame_equal(result, expected)

    result = pd.concat([df, df[["c", "b"]]], join="inner", sort=True, ignore_index=True)
    expected = expected[["b", "c"]]
    tm.assert_frame_equal(result, expected)


def test_concat_aligned_sort_does_not_raise():
    # GH-4588
    # We catch TypeErrors from sorting internally and do not re-raise.
    df = pd.DataFrame({1: [1, 2], "a": [3, 4]}, columns=[1, "a"])
    expected = pd.DataFrame({1: [1, 2, 1, 2], "a": [3, 4, 3, 4]}, columns=[1, "a"])
    result = pd.concat([df, df], ignore_index=True, sort=True)
    tm.assert_frame_equal(result, expected)
