import pytest
import numpy as np

import pandas as pd
import pandas.util.testing as tm
from pandas._libs.tslibs import iNaT
from pandas._libs.tslibs.period import IncompatibleFrequency
from pandas.core.arrays import PeriodArray


@pytest.mark.parametrize('key, value, expected', [
    ([0], pd.Period("2000", "D"), [10957, 1, 2]),
    ([0], None, [iNaT, 1, 2]),
    ([0, 1, 2], pd.Period("2000", "D"), [10957] * 3),
    ([0, 1, 2], [pd.Period("2000", "D"),
                 pd.Period("2001", "D"),
                 pd.Period("2002", "D")],
     [10957, 11323, 11688]),
])
def test_setitem(key, value, expected):
    arr = PeriodArray(np.arange(3), freq="D")
    expected = PeriodArray(expected, freq="D")
    arr[key] = value
    tm.assert_period_array_equal(arr, expected)


def test_setitem_raises():
    arr = PeriodArray(np.arange(3), freq="D")
    with tm.assert_raises_regex(IncompatibleFrequency, "freq"):
        arr[0] = pd.Period("2000", freq="A")

    with tm.assert_raises_regex(ValueError, "length"):
        arr[[0, 1]] = [pd.Period("2000", freq="D")]

    with tm.assert_raises_regex(TypeError, "int"):
        arr[0] = 1

