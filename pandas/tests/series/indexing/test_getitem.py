"""
Series.__getitem__ test classes are organized by the type of key passed.
"""

import numpy as np

import pandas as pd
from pandas import Series, period_range
import pandas._testing as tm


class TestSeriesGetitemScalars:
    pass


class TestSeriesGetitemSlices:
    def test_getitem_slice_2d(self, datetime_series):

        # This is currently failing because the test was relying on
        # the DeprecationWarning coming through Index.__getitem__.
        # We want to implement a warning specifically for Series.__getitem__
        # at which point this will become a Deprecation/FutureWarning
        with tm.assert_produces_warning(None):
            # GH#30588 multi-dimensional indexing deprecated
            result = datetime_series[:, np.newaxis]
        expected = datetime_series.values[:, np.newaxis]
        tm.assert_almost_equal(result, expected)


class TestGetitemListLike:
    def test_getitem_intlist_intindex_periodvalues(self):
        ser = Series(period_range("2000-01-01", periods=10, freq="D"))

        result = ser[[2, 4]]
        exp = pd.Series(
            [pd.Period("2000-01-03", freq="D"), pd.Period("2000-01-05", freq="D")],
            index=[2, 4],
            dtype="Period[D]",
        )
        tm.assert_series_equal(result, exp)
        assert result.dtype == "Period[D]"
