from __future__ import division
# pylint: disable-msg=W0612,E1101
import numpy as np
import pandas as pan
import pandas.util.testing as tm
from pandas.core.api import (DataFrame, Index, Series, notnull, isnull,
                             MultiIndex, DatetimeIndex, Timestamp, Period)
from pandas.util.testing import (assert_almost_equal,
                                 assert_series_equal,
                                 assert_frame_equal,
                                 makeCustomDataframe as mkdf,
                                 ensure_clean)

class TestDivUnderTruediv(object):
    def test_frame_div_dtype(self):
        p = DataFrame({"A": np.arange(10)})
        result = p.div(5)
        assert result.A.dtype.kind == "f", "Expected float dtype, instead saw %r" % result.A.dtype

    def test_series_div_dtype(self):
        p = Series(np.arange(10))
        result = p.div(4)
        assert result.dtype.kind == "f", "Expected float dtype, instead saw %r" % result.dtype

    def test_frame_div(self):
        p = DataFrame(tm.getIntegerSeriesData())
        result = p.div(3)
        expected = p.truediv(3)
        assert_frame_equal(result, expected)

        result = p.div(p.irow(0), axis=1)
        expected = p.truediv(p.irow(0), axis=1)
        assert_frame_equal(result, expected)

    def test_series_div(self):
        p = DataFrame(tm.getIntegerSeriesData())
        series = p.icol(0)
        result = series.div(5)
        expected = series.truediv(5)
        assert_series_equal(result, expected)

