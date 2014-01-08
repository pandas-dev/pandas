from pandas import compat
import nose

from numpy import nan
import numpy as np

from pandas import Series, DataFrame

from pandas.compat import product
from pandas.util.testing import (assert_frame_equal,
                                 assert_series_equal,
                                 assert_almost_equal)
import pandas.util.testing as tm


class TestRank(tm.TestCase):
    _multiprocess_can_split_ = True
    s = Series([1, 3, 4, 2, nan, 2, 1, 5, nan, 3])
    df = DataFrame({'A': s, 'B': s})

    results = {
        'average': np.array([1.5, 5.5, 7.0, 3.5, nan,
                             3.5, 1.5, 8.0, nan, 5.5]),
        'min': np.array([1, 5, 7, 3, nan, 3, 1, 8, nan, 5]),
        'max': np.array([2, 6, 7, 4, nan, 4, 2, 8, nan, 6]),
        'first': np.array([1, 5, 7, 3, nan, 4, 2, 8, nan, 6])
    }

    def test_rank_tie_methods(self):
        s = self.s

        def _check(s, expected, method='average'):
            result = s.rank(method=method)
            assert_almost_equal(result, expected)

        dtypes = [None, object]
        disabled = set([(object, 'first')])
        results = self.results

        for method, dtype in product(results, dtypes):
            if (dtype, method) in disabled:
                continue
            series = s if dtype is None else s.astype(dtype)
            _check(series, results[method], method=method)

    def test_rank_descending(self):
        dtypes = ['O', 'f8', 'i8']

        for dtype, method in product(dtypes, self.results):
            if 'i' in dtype:
                s = self.s.dropna()
                df = self.df.dropna()
            else:
                s = self.s.astype(dtype)
                df = self.df.astype(dtype)

            res = s.rank(ascending=False)
            expected = (s.max() - s).rank()
            assert_series_equal(res, expected)

            res = df.rank(ascending=False)
            expected = (df.max() - df).rank()
            assert_frame_equal(res, expected)

            if method == 'first' and dtype == 'O':
                continue

            expected = (s.max() - s).rank(method=method)
            res2 = s.rank(method=method, ascending=False)
            assert_series_equal(res2, expected)

            expected = (df.max() - df).rank(method=method)

            if dtype != 'O':
                res2 = df.rank(method=method, ascending=False,
                               numeric_only=True)
                assert_frame_equal(res2, expected)

            res3 = df.rank(method=method, ascending=False,
                           numeric_only=False)
            assert_frame_equal(res3, expected)

    def test_rank_2d_tie_methods(self):
        s = self.s
        df = self.df

        def _check2d(df, expected, method='average', axis=0):
            exp_df = DataFrame({'A': expected, 'B': expected})

            if axis == 1:
                df = df.T
                exp_df = exp_df.T

            result = df.rank(method=method, axis=axis)
            assert_frame_equal(result, exp_df)

        dtypes = [None, object]
        disabled = set([(object, 'first')])
        results = self.results

        for method, axis, dtype in product(results, [0, 1], dtypes):
            if (dtype, method) in disabled:
                continue
            frame = df if dtype is None else df.astype(dtype)
            _check2d(frame, results[method], method=method, axis=axis)

    def test_rank_int(self):
        s = self.s.dropna().astype('i8')

        for method, res in compat.iteritems(self.results):
            result = s.rank(method=method)
            expected = Series(res).dropna()
            expected.index = result.index
            assert_series_equal(result, expected)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
