import pytest

import pandas as pd


def test_idxmax():
    # GH#35257
    df = pd.DataFrame(
        {
            "a": [1, 5, 2, 4, 3, 5, 4, 2],
            "b": [9, 4, 3, 4, 2, 2, 4, 9],
            "c": [1, 8, 2, 4, 1, 5, 8, 3],
        },
        index=["A", "B", "C", "D", "E", "F", "G", "H"],
    )
    assert all(
        df.idxmax(keep="last") == pd.Series(["F", "H", "G"], index=["a", "b", "c"])
    )
    all_solution = pd.Series(
        [["B", "F"], ["A", "H"], ["B", "G"]], index=["a", "b", "c"]
    )
    for i in range(0, df.idxmax(keep="all").shape[0]):
        assert all(df.idxmax(keep="all")[i] == all_solution[i])
    assert all(
        df.idxmax(keep="first") == pd.Series(["B", "A", "B"], index=["a", "b", "c"])
    )
