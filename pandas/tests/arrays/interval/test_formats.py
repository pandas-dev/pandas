import numpy as np

import pandas as pd
from pandas.core.arrays import IntervalArray


def test_repr():
    # GH#25022
    arr = IntervalArray.from_tuples([(0, 1), (1, 2)])
    result = repr(arr)
    expected = (
        "<IntervalArray>\n[(0, 1], (1, 2]]\nLength: 2, dtype: interval[int64, right]"
    )
    assert result == expected


def test_series_repr_float_precision():
    # GH#25920
    arr = IntervalArray.from_arrays([0.55555555, np.nan], [1.55555555, np.nan])
    ser = pd.Series(arr)
    with pd.option_context("display.precision", 3):
        result = repr(ser)
    expected = "0    (0.556, 1.556]\n1               NaN\ndtype: interval"
    assert result == expected
