import pandas as pd


def test_argmax():
    # GH#35257
    s = pd.Series(
        [10, 11, 22, 13, 14, 22, 10, 17, 10, 20, 21, 22],
        index=["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"],
    )
    assert s.argmax(keep="last") == 11
    assert all(s.argmax(keep="all") == [2, 5, 11])
    assert s.argmax() == 2


def test_argmin():
    # GH#35257
    s = pd.Series(
        [10, 11, 22, 13, 14, 22, 10, 17, 10, 20, 21, 22],
        index=["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"],
    )
    assert s.argmin(keep="last") == 8
    assert all(s.argmin(keep="all") == [0, 6, 8])
    assert s.argmin() == 0


def test_idxmax():
    # GH#35257
    s = pd.Series(
        [10, 11, 22, 13, 14, 22, 10, 17, 10, 20, 21, 22],
        index=["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"],
    )
    assert s.idxmax(keep="last") == "l"
    assert all(s.idxmax(keep="all") == ["c", "f", "l"])
    assert s.idxmax() == "c"


def test_idxmin():
    # GH#35257
    s = pd.Series(
        [10, 11, 22, 13, 14, 22, 10, 17, 10, 20, 21, 22],
        index=["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"],
    )
    assert s.idxmin(keep="last") == "i"
    assert all(s.idxmin(keep="all") == ["a", "g", "i"])
    assert s.idxmin() == "a"
