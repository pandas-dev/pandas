import numpy as np

from pandas.core.dtypes.cast import dict_compat

import pandas as pd


def test_dict_compat():
    data_datetime64 = {np.datetime64("1990-03-15"): 1, np.datetime64("2015-03-15"): 2}
    data_unchanged = {1: 2, 3: 4, 5: 6}
    expected = {pd.Timestamp("1990-3-15"): 1, pd.Timestamp("2015-03-15"): 2}
    assert dict_compat(data_datetime64) == expected
    assert dict_compat(expected) == expected
    assert dict_compat(data_unchanged) == data_unchanged
