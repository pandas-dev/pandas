import pytest
import pandas as pd
import numpy as np
import pandas.testing as tm
from pandas.core.outliers import (
    is_outlier, 
    has_outliers, 
    drop_outliers, 
    fill_outliers
    )


# Data to test with
list_a = np.linspace(1, 5, num = 1000).tolist() + [100]
list_b = np.linspace(2, 4, num = 1000).tolist() + [75]

class TestIsOutlier:
    @pytest.mark.parametrize("method", ["iqr", "zscore"])
    def test_basic_detection(self, method):
        s = pd.Series(list_a)
        mask = is_outlier(s, method=method)
        assert mask.iloc[-1]

    def test_all_missing(self):
        s = pd.Series([np.nan, np.nan])
        mask = is_outlier(s)
        tm.assert_series_equal(mask, pd.Series([False, False]))
        
class TestHasOutliers:
    def test_series_has_outliers_true(self):
        s = pd.Series(list_a)
        result = has_outliers(s)
        assert result is True

    def test_series_has_outliers_false(self):
        s = pd.Series(list_a[:-1])
        result = has_outliers(s)
        assert result is False

    def test_dataframe_has_outliers_true(self):
        df = pd.DataFrame({"A": list_a, "B": list_b})
        result = has_outliers(df)
        assert result is True

    def test_dataframe_has_outliers_false(self):
        df = pd.DataFrame({"A": [1, 2, 3], "B": [10, 20, 30]})
        result = has_outliers(df)
        assert result is False

    def test_dataframe_with_non_numeric_column(self):
        df = pd.DataFrame({"A": list_a, "B": ["b"] * len(list_a)})
        result = has_outliers(df)
        assert result is True

    def test_empty_series(self):
        s = pd.Series([], dtype=float)
        result = has_outliers(s)
        assert result is False

    def test_empty_dataframe(self):
        df = pd.DataFrame({"A": []})
        result = has_outliers(df)
        assert result is False


class TestDropOutliers:
    def test_series_drop(self):
        s = pd.Series(list_a)
        result = drop_outliers(s)
        expected = pd.Series(list_a[:-1])
        tm.assert_series_equal(result, expected)

    def test_dataframe_drop_rows(self):
        df = pd.DataFrame({"A": list_a, "B": ["b"] * len(list_a)})
        result = drop_outliers(df, axis=0)
        expected = pd.DataFrame({"A": list_a[:-1], "B": ["b"] * len(list_a[:-1])})
        tm.assert_frame_equal(result, expected)

class TestFillOutliers:
    def test_fill_constant(self):
        s = pd.Series([1, 2, 3, 4, 100])
        result = fill_outliers(s, value=5)
        expected = pd.Series([1, 2, 3, 4, 5])
        tm.assert_series_equal(result, expected)

    def test_fill_median(self):
        s = pd.Series([1, 2, 3, 4, 100])
        result = fill_outliers(s, method={"filling": "median", "detection": "iqr"})
        expected = pd.Series([1.0, 2.0, 3.0, 4.0, 2.5])
        tm.assert_series_equal(result, expected)