# coding=utf-8

import pytest
import numpy as np
import pandas as pd
import pandas.util.testing as tm

def test_histogram():
    np.random.seed(3)
    s = pd.Series(np.random.normal(0, 1, 100))
    h, b = s.histogram(20)
    _h, _b = np.histogram(s, 20)
    assert np.all(h == _h)
    assert np.all(b == _b)
