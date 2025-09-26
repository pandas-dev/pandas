import numpy as np

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
