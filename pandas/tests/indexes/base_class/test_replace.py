import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


def test_index_replace():
    index = pd.Index([1, 2, 3])
    expected = pd.Index(["a", 2, "c"])

    result = index.replace([1, 3], ["a", "c"])

    tm.assert_equal(result, expected)


def test_index_replace_1():
    index = pd.Index([1, 2, 3])
    expected = pd.Index(["a", 2, 3])

    result = index.replace(1, "a")

    tm.assert_equal(result, expected)


def test_index_replace_2():
    index = pd.Index([1, 2, 3])
    expected = pd.Index(["a", 2, "c"])

    result = index.replace({1: "a", 3: "c"})
    tm.assert_equal(result, expected)


def test_index_replace_3():
    index = pd.Index([1, None, 2])
    with pytest.raises(NotImplementedError):
        index.replace(np.nan)


def test_index_replace_4():
    index = pd.Index([1, None, 2])
    expected = pd.Index(["a", None, "a"])

    result = index.replace([1, 2], "a")
    tm.assert_equal(expected, result)


def test_index_replace_5():
    index = pd.Index([0, 1, 2, 3, 4])
    expected = pd.Index([0, 3, 3, 3, 4])

    result = index.replace([1, 2], method="bfill")
    tm.assert_equal(expected, result)


def test_index_replace_6():
    index = pd.Index(["bat", "foo", "baait", "bar"])
    expected = pd.Index(["new", "foo", "baait", "new"])

    result = index.replace(to_replace=r"^ba.$", value="new", regex=True)
    tm.assert_equal(expected, result)


def test_index_replace_7():
    index = pd.Index(["bat", "foo", "baait", "bar"])
    expected = pd.Index(["new", "foo", "baait", "new"])

    result = index.replace(regex=r"^ba.$", value="new")
    tm.assert_equal(expected, result)


def test_index_replace_8():
    index = pd.Index(["bat", "foo", "baait", "bar"])
    expected = pd.Index(["new", "xyz", "baait", "new"])

    result = index.replace(regex={r"^ba.$": "new", "foo": "xyz"})
    tm.assert_equal(expected, result)


if __name__ == "__main__":
    # %load_ext autoreload
    # %autoreload 2
    index = pd.Index([1, 2, 3])
    index.replace([1, 2], ["a", "b"])
