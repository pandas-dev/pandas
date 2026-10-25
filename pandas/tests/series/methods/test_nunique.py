import numpy as np

import pandas as pd


def test_nunique():
    # basics.rst doc example
    series = pd.Series(np.random.default_rng(2).standard_normal(500))
    series[20:500] = np.nan
    series[10:20] = 5000
    result = series.nunique()
    assert result == 11


def test_nunique_categorical():
    # GH#18051
    ser = pd.Series(pd.Categorical([]))
    assert ser.nunique() == 0

    ser = pd.Series(pd.Categorical([np.nan]))
    assert ser.nunique() == 0
