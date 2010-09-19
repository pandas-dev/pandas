from pandas.core.common import notnull, isnull
import pandas.core.common as common

import numpy as np

def test_notnull():
    assert notnull(1.)
    assert not notnull(None)
    assert not notnull(np.NaN)
    assert not notnull(np.inf)
    assert not notnull(-np.inf)

def test_isnull():
    assert not isnull(1.)
    assert isnull(None)
    assert isnull(np.NaN)
    assert isnull(np.inf)
    assert isnull(-np.inf)

