import pytest

import pandas as pd
from pandas.tests.arrays.masked_shared import (
    ComparisonOps,
    NumericOps,
)


class TestComparisonOps(NumericOps, ComparisonOps):
    @pytest.mark.parametrize("other", [True, False, pd.NA, -1, 0, 1])
    def test_scalar(self, other, all_compare_operators, dtype):
        ComparisonOps.test_scalar(self, other, all_compare_operators, dtype)

    def test_compare_to_int(self, dtype, all_compare_operators):
        # GH 28930
        s1 = pd.Series([1, None, 3], dtype=dtype)
        s2 = pd.Series([1, None, 3], dtype="float")

        method = getattr(s1, all_compare_operators)
        result = method(2)

        method = getattr(s2, all_compare_operators)
        expected = method(2).astype("boolean")
        expected[s2.isna()] = pd.NA

        self.assert_series_equal(result, expected)


def test_equals():
    # GH-30652
    # equals is generally tested in /tests/extension/base/methods, but this
    # specifically tests that two arrays of the same class but different dtype
    # do not evaluate equal
    a1 = pd.array([1, 2, None], dtype="Int64")
    a2 = pd.array([1, 2, None], dtype="Int32")
    assert a1.equals(a2) is False
