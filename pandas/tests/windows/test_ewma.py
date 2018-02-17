from datetime import datetime

import numpy as np
import pytest
from numpy.random import randn

import pandas.core.window as rwindow
import pandas.util.testing as tm
from pandas import DataFrame, Series, bdate_range
from pandas.errors import UnsupportedFunctionCall

N, K = 100, 10


def assert_equal(left, right):
    if isinstance(left, Series):
        tm.assert_series_equal(left, right)
    else:
        tm.assert_frame_equal(left, right)


class Base(object):

    _nan_locs = np.arange(20, 40)
    _inf_locs = np.array([])

    def _create_data(self):
        arr = randn(N)
        arr[self._nan_locs] = np.NaN

        self.arr = arr
        self.rng = bdate_range(datetime(2009, 1, 1), periods=N)
        self.series = Series(arr.copy(), index=self.rng)
        self.frame = DataFrame(randn(N, K), index=self.rng,
                               columns=np.arange(K))


class TestEWM(Base):

    def setup_method(self, method):
        self._create_data()

    def test_doc_string(self):

        df = DataFrame({'B': [0, 1, 2, np.nan, 4]})
        df
        df.ewm(com=0.5).mean()

    def test_constructor(self):
        for o in [self.series, self.frame]:
            c = o.ewm

            # valid
            c(com=0.5)
            c(span=1.5)
            c(alpha=0.5)
            c(halflife=0.75)
            c(com=0.5, span=None)
            c(alpha=0.5, com=None)
            c(halflife=0.75, alpha=None)

            # not valid: mutually exclusive
            with pytest.raises(ValueError):
                c(com=0.5, alpha=0.5)
            with pytest.raises(ValueError):
                c(span=1.5, halflife=0.75)
            with pytest.raises(ValueError):
                c(alpha=0.5, span=1.5)

            # not valid: com < 0
            with pytest.raises(ValueError):
                c(com=-0.5)

            # not valid: span < 1
            with pytest.raises(ValueError):
                c(span=0.5)

            # not valid: halflife <= 0
            with pytest.raises(ValueError):
                c(halflife=0)

            # not valid: alpha <= 0 or alpha > 1
            for alpha in (-0.5, 1.5):
                with pytest.raises(ValueError):
                    c(alpha=alpha)

    def test_numpy_compat(self):
        # see gh-12811
        e = rwindow.EWM(Series([2, 4, 6]), alpha=0.5)

        msg = "numpy operations are not valid with window objects"

        for func in ('std', 'mean', 'var'):
            tm.assert_raises_regex(UnsupportedFunctionCall, msg,
                                   getattr(e, func), 1, 2, 3)
            tm.assert_raises_regex(UnsupportedFunctionCall, msg,
                                   getattr(e, func), dtype=np.float64)
