import numpy as np
import pytest

import pandas as pd


def test_contains_nan(using_nan_is_na):
    # GH#52840
    arr = pd.array(range(5)) / 0

    assert np.isnan(arr._data[0])
    if using_nan_is_na:
        assert arr.isna()[0]
    else:
        assert not arr.isna()[0]
    assert np.nan in arr


def test_hasnans_from_arrow(using_nan_is_na):
    # GH#49818
    pyarrow = pytest.importorskip("pyarrow")
    arr = pyarrow.array([np.nan, 1, 2])
    ser = pd.Series(pd.Float64Dtype().__from_arrow__(arr))
    if using_nan_is_na:
        # NaN folded into mask -> hasnans is True
        assert ser.hasnans is True
    else:
        # NaN is a value, not in the mask -> hasnans is False
        assert ser.hasnans is False
        # But NaN still exists in the underlying data
        assert np.isnan(ser._values._data[0])
