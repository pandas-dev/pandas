import numpy as np

from pandas import DataFrame, Series
import pandas._testing as tm
from pandas.arrays import SparseArray


def test_is_sparse():
    # GH26706 Add is_sparse method to check for sparse columns in a DataFrame
    df = DataFrame({"A": SparseArray([1, np.nan, 1]), "B": [1, 2, 3]})
    result = df.is_sparse

    expected = Series([True, False], index=["A", "B"])
    tm.assert_series_equal(result, expected)
