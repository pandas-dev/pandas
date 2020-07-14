import pandas as pd
import pandas._testing as tm


def test_idxmax():
    # GH#35257
    df = pd.DataFrame(
        {
            "a": [1, 5, 2, 4, 3, 5, 4, 2],
            "b": [9, 4, 3, 4, 2, 2, 4, 9],
            "c": [1, 8, 2, 4, 1, 5, 8, 3],
        },
        index=["A", "B", "C", "D", "E", "F", "G", "H"]
    )
    assert all(
        df.idxmax(keep="last") == pd.Series(["F", "H", "G"], index=["a", "b", "c"])
    )
    tm.assert_series_equal(
        df.idxmax(keep="all"),
        pd.Series([["B", "F"], ["A", "H"], ["B", "G"]], index=["a", "b", "c"]),
    )
    assert all(
        df.idxmax(keep="first") == pd.Series(["B", "A", "B"], index=["a", "b", "c"])
    )

def test_idxmin():
    # GH#35257
    df = pd.DataFrame(
        {
            "a": [1, 5, 2, 4, 3, 5, 4, 2],
            "b": [9, 4, 3, 4, 2, 2, 4, 9],
            "c": [1, 8, 2, 4, 1, 5, 8, 3],
        },
        index=["A", "B", "C", "D", "E", "F", "G", "H"]
    )
    assert all(
        df.idxmin(keep="last") == pd.Series(["A", "F", "E"], index=["a", "b", "c"])
    )
    tm.assert_series_equal(
        df.idxmin(keep="all"),
        pd.Series([["A"], ["E", "F"], ["A", "E"]], index=["a", "b", "c"]),
    )
    assert all(
        df.idxmin(keep="first") == pd.Series(["A", "E", "A"], index=["a", "b", "c"])
    )