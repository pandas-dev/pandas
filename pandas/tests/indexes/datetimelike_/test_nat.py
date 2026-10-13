import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm


@pytest.mark.parametrize(
    "index_without_na",
    [
        pd.TimedeltaIndex(["1 days", "2 days"]),
        pd.PeriodIndex(["2011-01-01", "2011-01-02"], freq="D"),
        pd.DatetimeIndex(["2011-01-01", "2011-01-02"]),
        pd.DatetimeIndex(["2011-01-01", "2011-01-02"], tz="UTC"),
    ],
)
def test_nat(index_without_na):
    empty_index = index_without_na[:0]

    index_with_na = index_without_na.copy(deep=True)
    index_with_na._data[1] = pd.NaT

    assert empty_index._na_value is pd.NaT
    assert index_with_na._na_value is pd.NaT
    assert index_without_na._na_value is pd.NaT

    idx = index_without_na
    assert idx._can_hold_na

    tm.assert_numpy_array_equal(idx._isnan, np.array([False, False]))
    assert idx.hasnans is False

    idx = index_with_na
    assert idx._can_hold_na

    tm.assert_numpy_array_equal(idx._isnan, np.array([False, True]))
    assert idx.hasnans is True
