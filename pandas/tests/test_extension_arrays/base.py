import pandas as pd
import pandas.util.testing as tm
from pandas.core.internals import ExtensionBlock


class BaseArrayTests:

    def test_series_constructor(self, test_data):
        result = pd.Series(test_data)
        assert result.dtype == test_data.dtype
        assert len(result) == len(test_data)
        assert isinstance(result._data.blocks[0], ExtensionBlock)

    def test_dataframe_constructor(self, test_data):
        result = pd.DataFrame({"A": test_data})
        assert result.dtypes['A'] == test_data.dtype
        assert result.shape == (len(test_data), 1)
        assert isinstance(result._data.blocks[0], ExtensionBlock)

    def test_concat(self, test_data):
        result = pd.concat([
            pd.Series(test_data),
            pd.Series(test_data),
        ], ignore_index=True)
        assert len(result) == len(test_data) * 2

    def test_iloc(self, test_data):
        ser = pd.Series(test_data)
        result = ser.iloc[:4]
        expected = pd.Series(test_data[:4])
        tm.assert_series_equal(result, expected)

    def test_loc(self, test_data):
        ser = pd.Series(test_data)
        result = ser.loc[[0, 1, 2, 3]]
        expected = pd.Series(test_data[:4])
        tm.assert_series_equal(result, expected)
