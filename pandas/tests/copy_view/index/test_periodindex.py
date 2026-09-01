import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm
from pandas.tests.copy_view.util import get_array


@pytest.mark.parametrize("box", [lambda x: x, pd.PeriodIndex])
def test_periodindex(box):
    dt = pd.period_range("2019-12-31", periods=3, freq="D")
    ser = pd.Series(dt)
    idx = box(pd.PeriodIndex(ser))
    expected = idx.copy(deep=True)
    ser.iloc[0] = pd.Period("2020-12-31")
    tm.assert_index_equal(idx, expected)


def test_constructor_copy_input_period_ea_default():
    # GH 63388
    arr = pd.array(["2020-01-01", "2020-01-02"], dtype="period[D]")
    idx = pd.PeriodIndex(arr)
    assert not tm.shares_memory(arr, idx.array)


def test_series_from_temporary_periodindex_readonly_data():
    # GH 63388
    arr = pd.array(["2020-01-01", "2020-01-02"], dtype="period[D]")
    arr._ndarray.flags.writeable = False
    ser = pd.Series(pd.PeriodIndex(arr))
    assert not np.shares_memory(arr._ndarray, get_array(ser))
    ser.iloc[0] = pd.Period("2022-01-01", freq="D")
    expected = pd.Series(
        [pd.Period("2022-01-01", freq="D"), pd.Period("2020-01-02", freq="D")],
        dtype="period[D]",
    )
    tm.assert_series_equal(ser, expected)
