"""
Shareable tests for &, |, ^
"""
import operator

import numpy as np
import pytest

from pandas import DataFrame, Series
import pandas._testing as tm


class TestLogicalOps:
    @pytest.mark.parametrize(
        "left, right, op, expected",
        [
            (
                [True, False, np.nan],
                [True, False, True],
                operator.and_,
                [True, False, False],
            ),
            (
                [True, False, True],
                [True, False, np.nan],
                operator.and_,
                [True, False, False],
            ),
            (
                [True, False, np.nan],
                [True, False, True],
                operator.or_,
                [True, False, False],
            ),
            (
                [True, False, True],
                [True, False, np.nan],
                operator.or_,
                [True, False, True],
            ),
        ],
    )
    @pytest.mark.parametrize("klass", [Series, DataFrame])
    def test_logical_operators_nans(self, left, right, op, expected, klass):
        # GH#13896
        result = op(klass(left), klass(right))
        expected = klass(expected)

        tm.assert_equal(result, expected)
