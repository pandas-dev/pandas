import pytest

import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize("keep", ["first", "last", "all", "invalid_parameter"])
def test_argmax(keep):
    # GH#35257
    s = pd.Series(
        [10, 11, 22, 13, 14, 22, 10, 17, 10, 20, 21, 22],
        index=["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"],
    )
    if keep == "first":
        expected = 2
    elif keep == "last":
        expected = 11
    elif keep == "all":
        expected = [2, 5, 11]
        result = s.argmax(keep=keep)
        return tm.equalContents(result, expected)
    else:
        try:
            result = s.argmax(keep=keep)
        except ValueError:
            return True
    result = s.argmax(keep=keep)
    assert result == expected


@pytest.mark.parametrize("keep", ["first", "last", "all", "invalid_parameter"])
def test_argmin(keep):
    # GH#35257
    s = pd.Series(
        [10, 11, 22, 13, 14, 22, 10, 17, 10, 20, 21, 22],
        index=["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"],
    )
    if keep == "first":
        expected = 0
    elif keep == "last":
        expected = 8
    elif keep == "all":
        expected = [0, 6, 8]
        result = s.argmin(keep=keep)
        return tm.equalContents(result, expected)
    else:
        try:
            result = s.argmin(keep=keep)
        except ValueError:
            return True
    result = s.argmin(keep=keep)
    assert result == expected


@pytest.mark.parametrize("keep", ["first", "last", "all", "invalid_parameter"])
def test_idxmax(keep):
    # GH#35257
    s = pd.Series(
        [10, 11, 22, 13, 14, 22, 10, 17, 10, 20, 21, 22],
        index=["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"],
    )
    if keep == "first":
        expected = "c"
    elif keep == "last":
        expected = "l"
    elif keep == "all":
        expected = ["c", "f", "l"]
        result = s.idxmax(keep=keep)
        return tm.equalContents(result, expected)
    else:
        try:
            result = s.idxmax(keep=keep)
        except ValueError:
            return True
    result = s.idxmax(keep=keep)
    assert result == expected


@pytest.mark.parametrize("keep", ["first", "last", "all", "invalid_parameter"])
def test_idxmin(keep):
    # GH#35257
    s = pd.Series(
        [10, 11, 22, 13, 14, 22, 10, 17, 10, 20, 21, 22],
        index=["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"],
    )
    if keep == "first":
        expected = "a"
    elif keep == "last":
        expected = "i"
    elif keep == "all":
        expected = ["a", "g", "i"]
        result = s.idxmin(keep=keep)
        return tm.equalContents(result, expected)
    else:
        try:
            result = s.idxmin(keep=keep)
        except ValueError:
            return True
    result = s.idxmin(keep=keep)
    assert result == expected
