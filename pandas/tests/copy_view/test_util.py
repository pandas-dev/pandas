import numpy as np

import pandas as pd
from pandas.tests.copy_view.util import get_array


def test_get_array_numpy():
    df = pd.DataFrame({"a": [1, 2, 3]})
    assert np.shares_memory(get_array(df, "a"), get_array(df, "a"))


def test_get_array_masked():
    df = pd.DataFrame({"a": [1, 2, 3]}, dtype="Int64")
    assert np.shares_memory(get_array(df, "a"), get_array(df, "a"))
