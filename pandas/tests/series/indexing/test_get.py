import numpy as np

import pandas as pd
from pandas import Series


def test_get():
    # GH 6383
    s = Series(
        np.array(
            [
                43,
                48,
                60,
                48,
                50,
                51,
                50,
                45,
                57,
                48,
                56,
                45,
                51,
                39,
                55,
                43,
                54,
                52,
                51,
                54,
            ]
        )
    )

    result = s.get(25, 0)
    expected = 0
    assert result == expected

    s = Series(
        np.array(
            [
                43,
                48,
                60,
                48,
                50,
                51,
                50,
                45,
                57,
                48,
                56,
                45,
                51,
                39,
                55,
                43,
                54,
                52,
                51,
                54,
            ]
        ),
        index=pd.Float64Index(
            [
                25.0,
                36.0,
                49.0,
                64.0,
                81.0,
                100.0,
                121.0,
                144.0,
                169.0,
                196.0,
                1225.0,
                1296.0,
                1369.0,
                1444.0,
                1521.0,
                1600.0,
                1681.0,
                1764.0,
                1849.0,
                1936.0,
            ]
        ),
    )

    result = s.get(25, 0)
    expected = 43
    assert result == expected

    # GH 7407
    # with a boolean accessor
    df = pd.DataFrame({"i": [0] * 3, "b": [False] * 3})
    vc = df.i.value_counts()
    result = vc.get(99, default="Missing")
    assert result == "Missing"

    vc = df.b.value_counts()
    result = vc.get(False, default="Missing")
    assert result == 3

    result = vc.get(True, default="Missing")
    assert result == "Missing"


def test_get_nan():
    # GH 8569
    s = pd.Float64Index(range(10)).to_series()
    assert s.get(np.nan) is None
    assert s.get(np.nan, default="Missing") == "Missing"


def test_get_nan_multiple():
    # GH 8569
    # ensure that fixing "test_get_nan" above hasn't broken get
    # with multiple elements
    s = pd.Float64Index(range(10)).to_series()

    idx = [2, 30]
    assert s.get(idx) is None

    idx = [2, np.nan]
    assert s.get(idx) is None

    # GH 17295 - all missing keys
    idx = [20, 30]
    assert s.get(idx) is None

    idx = [np.nan, np.nan]
    assert s.get(idx) is None


def test_get_with_default():
    # GH#7725
    d0 = ["a", "b", "c", "d"]
    d1 = np.arange(4, dtype="int64")
    others = ["e", 10]

    for data, index in ((d0, d1), (d1, d0)):
        s = Series(data, index=index)
        for i, d in zip(index, data):
            assert s.get(i) == d
            assert s.get(i, d) == d
            assert s.get(i, "z") == d
            for other in others:
                assert s.get(other, "z") == "z"
                assert s.get(other, other) == other
