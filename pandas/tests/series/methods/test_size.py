import pytest

from pandas import Series


@pytest.mark.parametrize(
    "data, index, expected",
    [
        ([1, 2, 3], None, 3),
        ({"a": 1, "b": 2, "c": 3}, None, 3),
        ([1, 2, 3], ["x", "y", "z"], 3),
        ([1, 2, 3, 4, 5], ["x", "y", "z", "w", "n"], 5),
        ([1, 2, 3], None, 3),
        ([1, 2, 3], ["x", "y", "z"], 3),
        ([1, 2, 3, 4], ["x", "y", "z", "w"], 4),
    ],
)
def test_series(data, index, expected):
    s = Series(data, index=index)
    assert s.size == expected
    assert isinstance(s.size, int)
