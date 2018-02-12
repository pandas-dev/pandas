from datetime import datetime

import numpy as np
import pytest
from numpy.random import randn

import pandas as pd
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


class TestExpanding(Base):

    def setup_method(self, method):
        self._create_data()

    def test_doc_string(self):

        df = DataFrame({'B': [0, 1, 2, np.nan, 4]})
        df
        df.expanding(2).sum()

    def test_constructor(self):
        # GH 12669

        for o in [self.series, self.frame]:
            c = o.expanding

            # valid
            c(min_periods=1)
            c(min_periods=1, center=True)
            c(min_periods=1, center=False)

            # not valid
            for w in [2., 'foo', np.array([2])]:
                with pytest.raises(ValueError):
                    c(min_periods=w)
                with pytest.raises(ValueError):
                    c(min_periods=1, center=w)

    def test_numpy_compat(self):
        # see gh-12811
        e = rwindow.Expanding(Series([2, 4, 6]), window=2)

        msg = "numpy operations are not valid with window objects"

        for func in ('std', 'mean', 'sum', 'max', 'min', 'var'):
            tm.assert_raises_regex(UnsupportedFunctionCall, msg,
                                   getattr(e, func), 1, 2, 3)
            tm.assert_raises_regex(UnsupportedFunctionCall, msg,
                                   getattr(e, func), dtype=np.float64)

    @pytest.mark.parametrize(
        'expander',
        [1, pytest.param('ls', marks=pytest.mark.xfail(
                         reason='GH 16425 expanding with '
                                'offset not supported'))])
    def test_empty_df_expanding(self, expander):
        # GH 15819 Verifies that datetime and integer expanding windows can be
        # applied to empty DataFrames

        expected = DataFrame()
        result = DataFrame().expanding(expander).sum()
        tm.assert_frame_equal(result, expected)

        # Verifies that datetime and integer expanding windows can be applied
        # to empty DataFrames with datetime index
        expected = DataFrame(index=pd.DatetimeIndex([]))
        result = DataFrame(
            index=pd.DatetimeIndex([])).expanding(expander).sum()
        tm.assert_frame_equal(result, expected)

    def test_missing_minp_zero(self):
        # https://github.com/pandas-dev/pandas/pull/18921
        # minp=0
        x = pd.Series([np.nan])
        result = x.expanding(min_periods=0).sum()
        expected = pd.Series([0.0])
        tm.assert_series_equal(result, expected)

        # minp=1
        result = x.expanding(min_periods=1).sum()
        expected = pd.Series([np.nan])
        tm.assert_series_equal(result, expected)
