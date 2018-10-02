import operator
from decimal import Decimal

import pandas as pd
import pandas.util.testing as tm
from pandas.tests.extension.decimal_array import DecimalArray


def test_combine_from_sequence_raises():
    # https://github.com/pandas-dev/pandas/issues/22850
    class BadDecimalArray(DecimalArray):
        def _from_sequence(cls, scalars, dtype=None, copy=False):
            raise KeyError("For the test")

    ser = pd.Series(BadDecimalArray([Decimal("1.0"), Decimal("2.0")]))
    result = ser.combine(ser, operator.add)

    # note: object dtype
    expected = pd.Series([Decimal("2.0"), Decimal("4.0")], dtype="object")
    tm.assert_series_equal(result, expected)
