from numpy import nan
from pandas import Series
import pandas.util.testing as tm


def test_diff():
    '''
    Tests the pd.Series diff function on boolean series.
    '''
    data = Series(
        [0, -1, -2, -3, -4, -3, -2, -1, 0]
    )
    filtered = data.between(-2, 0, inclusive=True)
    diff_boolean = filtered.diff()
    expected_boolean = Series(
        [nan, False, False, True, False, False, True, False, False]
    )
    tm.assert_series_equal(diff_boolean, expected_boolean)
    diff_data = data.diff()
    expected_data = Series(
        [nan, -1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0]
    )
    tm.assert_series_equal(diff_data, expected_data)
