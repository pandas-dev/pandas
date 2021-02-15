import pytest

import pandas as pd
import pandas._testing as tm


def test_only_unique():
    df = pd.DataFrame(
        {"A": list("abbacc"), "B": list("122133"), "C": [1, 2, 2, 1, 3, 3]}
    )

    expected = pd.DataFrame({"A": list("abc"), "B": list("123"), "C": [1, 2, 3]})
    result = df.groupby("A", as_index=False).only_unique()
    tm.assert_frame_equal(result, expected)


def test_only_unique_as_index():
    df = pd.DataFrame({"A": list("abbacc"), "B": list("122133")})
    expected = pd.DataFrame({"A": list("abc"), "B": list("123")})
    result = df.groupby("A").only_unique()
    tm.assert_frame_equal(result, expected.set_index("A"))


def test_only_unique_aggregate():
    df = pd.DataFrame({"A": list("abbacc"), "B": list("122133")})
    expected = pd.DataFrame({"A": list("abc"), "B": list("123")})
    result = df.groupby("A").aggregate("only_unique")
    tm.assert_frame_equal(result, expected.set_index("A"))


def test_not_only_unique():
    df = pd.DataFrame({"A": list("abbacc"), "B": list("122333")})

    with pytest.raises(ValueError):
        df.groupby("A").only_unique()


def test_with_missing():
    df = pd.DataFrame({"A": list("abbcc"), "B": [None, 1, 1, 3, 3]})

    expected = pd.DataFrame({"A": list("abc"), "B": [None, 1, 3]}).set_index("A")
    result = df.groupby("A").only_unique()
    tm.assert_frame_equal(result, expected)


@pytest.mark.xfail(reason="missing values not yet handled properly")
def test_with_multiple_missing():
    df = pd.DataFrame({"A": list("abbcc"), "B": [None, 1, 1, None, 3, 3]})

    expected = pd.DataFrame({"A": list("abc"), "B": [None, 1, 3]}).set_index("A")
    result = df.groupby("A").only_unique()
    tm.assert_frame_equal(result, expected)
