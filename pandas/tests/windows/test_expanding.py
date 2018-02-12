import pytest

from numpy.random import randn
import numpy as np
from datetime import datetime
import pandas as pd
from pandas import (Series, DataFrame, bdate_range,
                    notna)
import pandas.core.window as rwindow
from pandas.errors import UnsupportedFunctionCall
import pandas.util.testing as tm

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


# create the data only once as we are not setting it
def _create_consistency_data():
    def create_series():
        return [Series(),
                Series([np.nan]),
                Series([np.nan, np.nan]),
                Series([3.]),
                Series([np.nan, 3.]),
                Series([3., np.nan]),
                Series([1., 3.]),
                Series([2., 2.]),
                Series([3., 1.]),
                Series([5., 5., 5., 5., np.nan, np.nan, np.nan, 5., 5., np.nan,
                        np.nan]),
                Series([np.nan, 5., 5., 5., np.nan, np.nan, np.nan, 5., 5.,
                        np.nan, np.nan]),
                Series([np.nan, np.nan, 5., 5., np.nan, np.nan, np.nan, 5., 5.,
                        np.nan, np.nan]),
                Series([np.nan, 3., np.nan, 3., 4., 5., 6., np.nan, np.nan, 7.,
                        12., 13., 14., 15.]),
                Series([np.nan, 5., np.nan, 2., 4., 0., 9., np.nan, np.nan, 3.,
                        12., 13., 14., 15.]),
                Series([2., 3., np.nan, 3., 4., 5., 6., np.nan, np.nan, 7.,
                        12., 13., 14., 15.]),
                Series([2., 5., np.nan, 2., 4., 0., 9., np.nan, np.nan, 3.,
                        12., 13., 14., 15.]),
                Series(range(10)),
                Series(range(20, 0, -2)), ]

    def create_dataframes():
        return ([DataFrame(),
                 DataFrame(columns=['a']),
                 DataFrame(columns=['a', 'a']),
                 DataFrame(columns=['a', 'b']),
                 DataFrame(np.arange(10).reshape((5, 2))),
                 DataFrame(np.arange(25).reshape((5, 5))),
                 DataFrame(np.arange(25).reshape((5, 5)),
                           columns=['a', 'b', 99, 'd', 'd'])] +
                [DataFrame(s) for s in create_series()])

    def is_constant(x):
        values = x.values.ravel()
        return len(set(values[notna(values)])) == 1

    def no_nans(x):
        return x.notna().all().all()

    # data is a tuple(object, is_contant, no_nans)
    data = create_series() + create_dataframes()

    return [(x, is_constant(x), no_nans(x)) for x in data]


_consistency_data = _create_consistency_data()


class TestGrouperGrouping(object):

    def setup_method(self, method):
        self.series = Series(np.arange(10))
        self.frame = DataFrame({'A': [1] * 20 + [2] * 12 + [3] * 8,
                                'B': np.arange(40)})

    def test_mutated(self):

        def f():
            self.frame.groupby('A', foo=1)
        pytest.raises(TypeError, f)

        g = self.frame.groupby('A')
        assert not g.mutated
        g = self.frame.groupby('A', mutated=True)
        assert g.mutated

    def test_getitem(self):
        g = self.frame.groupby('A')
        g_mutated = self.frame.groupby('A', mutated=True)

        expected = g_mutated.B.apply(lambda x: x.rolling(2).mean())

        result = g.rolling(2).mean().B
        tm.assert_series_equal(result, expected)

        result = g.rolling(2).B.mean()
        tm.assert_series_equal(result, expected)

        result = g.B.rolling(2).mean()
        tm.assert_series_equal(result, expected)

        result = self.frame.B.groupby(self.frame.A).rolling(2).mean()
        tm.assert_series_equal(result, expected)

    def test_getitem_multiple(self):

        # GH 13174
        g = self.frame.groupby('A')
        r = g.rolling(2)
        g_mutated = self.frame.groupby('A', mutated=True)
        expected = g_mutated.B.apply(lambda x: x.rolling(2).count())

        result = r.B.count()
        tm.assert_series_equal(result, expected)

        result = r.B.count()
        tm.assert_series_equal(result, expected)

    def test_expanding(self):
        g = self.frame.groupby('A')
        r = g.expanding()

        for f in ['sum', 'mean', 'min', 'max', 'count', 'kurt', 'skew']:

            result = getattr(r, f)()
            expected = g.apply(lambda x: getattr(x.expanding(), f)())
            tm.assert_frame_equal(result, expected)

        for f in ['std', 'var']:
            result = getattr(r, f)(ddof=0)
            expected = g.apply(lambda x: getattr(x.expanding(), f)(ddof=0))
            tm.assert_frame_equal(result, expected)

        result = r.quantile(0.5)
        expected = g.apply(lambda x: x.expanding().quantile(0.5))
        tm.assert_frame_equal(result, expected)

    def test_expanding_corr_cov(self):
        g = self.frame.groupby('A')
        r = g.expanding()

        for f in ['corr', 'cov']:
            result = getattr(r, f)(self.frame)

            def func(x):
                return getattr(x.expanding(), f)(self.frame)
            expected = g.apply(func)
            tm.assert_frame_equal(result, expected)

            result = getattr(r.B, f)(pairwise=True)

            def func(x):
                return getattr(x.B.expanding(), f)(pairwise=True)
            expected = g.apply(func)
            tm.assert_series_equal(result, expected)

    def test_expanding_apply(self):
        g = self.frame.groupby('A')
        r = g.expanding()

        # reduction
        result = r.apply(lambda x: x.sum())
        expected = g.apply(lambda x: x.expanding().apply(lambda y: y.sum()))
        tm.assert_frame_equal(result, expected)
