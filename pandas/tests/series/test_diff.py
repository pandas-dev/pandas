from pandas import (
    Series
)

from numpy import nan


def test_diff():
    data = Series(
        [0, -1, -2, -3, -4, -3, -2, -1, 0, -1, -1, 0, -1, -2, -3, -2, 0]
    )
    filtered = data.between(-2, 0, inclusive=True)
    diff_boolean = filtered.diff()
    expected_boolean = Series(
        [nan, False, False, True, False, False, True, False, False,
         False, False, False, False, False, True, True, False]
    )
    assert diff_boolean.equals(expected_boolean)
    diff_data = data.diff()
    expected_data = Series(
        [nan, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0,
         -1.0, 0.0, 1.0, -1.0, -1.0, -1.0, 1.0, 2.0]
    )
    assert diff_data.equals(expected_data)
