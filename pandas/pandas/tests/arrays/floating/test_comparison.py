import pytest

import pandas as pd
import pandas._testing as tm
from pandas.tests.arrays.masked_shared import (
    ComparisonOps,
    NumericOps,
)


class TestComparisonOps(NumericOps, ComparisonOps):
    @pytest.mark.parametrize("other", [True, False, pd.NA, -1.0, 0.0, 1])
    def test_scalar(self, other, all_compare_operators, dtype):
        ComparisonOps.test_scalar(self, other, all_compare_operators, dtype)

    def test_compare_with_integerarray(self, all_compare_operators):
        op = self.get_op_from_name(all_compare_operators)
        a = pd.array([0, 1, None] * 3, dtype="Int64")
        b = pd.array([0] * 3 + [1] * 3 + [None] * 3, dtype="Float64")
        other = b.astype("Int64")
        expected = op(a, other)
        result = op(a, b)
        tm.assert_extension_array_equal(result, expected)
        expected = op(other, a)
        result = op(b, a)
        tm.assert_extension_array_equal(result, expected)


def test_equals():
    # GH-30652
    # equals is generally tested in /tests/extension/base/methods, but this
    # specifically tests that two arrays of the same class but different dtype
    # do not evaluate equal
    a1 = pd.array([1, 2, None], dtype="Float64")
    a2 = pd.array([1, 2, None], dtype="Float32")
    assert a1.equals(a2) is False
