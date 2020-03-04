import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


def test_index_replace():
    index = pd.Index([1, 2, 3])
    expected = pd.Index(["a", 2, "c"])

    result = index.replace([1, 3], ["a", "c"])

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


if __name__ == "__main__":
    # %load_ext autoreload
    # %autoreload 2
    index = pd.Index([1, 2, 3])
    index.replace([1, 2], ["a", "b"])
