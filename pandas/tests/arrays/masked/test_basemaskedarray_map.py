from pandas import (
    NA,
    Series,
    isna,
)
from pandas.testing import assert_series_equal


def test_basemaskedarray_map():
    s = Series([1, 2, None, 4], dtype="Int32")

    def transform(x):
        if isna(x):
            return x
        return x + 1

    result = s.map(transform)
    expected = Series([2, 3, NA, 5], dtype=result.dtype)

    assert_series_equal(result, expected)
