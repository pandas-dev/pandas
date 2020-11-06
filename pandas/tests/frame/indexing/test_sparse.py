import pandas as pd
import pandas._testing as tm
from pandas.arrays import SparseArray


class TestSparseDataFrameIndexing:
    def test_getitem_sparse_column(self):
        # https://github.com/pandas-dev/pandas/issues/23559
        data = SparseArray([0, 1])
        df = pd.DataFrame({"A": data})
        expected = pd.Series(data, name="A")
        result = df["A"]
        tm.assert_series_equal(result, expected)

        result = df.iloc[:, 0]
        tm.assert_series_equal(result, expected)

        result = df.loc[:, "A"]
        tm.assert_series_equal(result, expected)
