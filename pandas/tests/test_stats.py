# -*- coding: utf-8 -*-
from pandas import compat
import nose

from distutils.version import LooseVersion
from numpy import nan
import numpy as np

from pandas import Series, DataFrame

from pandas.compat import product
from pandas.util.testing import (assert_frame_equal, assert_series_equal)
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
        'first': np.array([1, 5, 7, 3, nan, 4, 2, 8, nan, 6]),
        'dense': np.array([1, 3, 4, 2, nan, 2, 1, 5, nan, 3]),
    }

    def test_rank_tie_methods(self):
        s = self.s

        def _check(s, expected, method='average'):
            result = s.rank(method=method)
            tm.assert_series_equal(result, Series(expected))

        dtypes = [None, object]
        disabled = set([(object, 'first')])
        results = self.results

        for method, dtype in product(results, dtypes):
            if (dtype, method) in disabled:
                continue
            series = s if dtype is None else s.astype(dtype)
            _check(series, results[method], method=method)

    def test_rank_methods_series(self):
        tm.skip_if_no_package('scipy', '0.13', 'scipy.stats.rankdata')
        import scipy
        from scipy.stats import rankdata

        xs = np.random.randn(9)
        xs = np.concatenate([xs[i:] for i in range(0, 9, 2)])  # add duplicates
        np.random.shuffle(xs)

        index = [chr(ord('a') + i) for i in range(len(xs))]

        for vals in [xs, xs + 1e6, xs * 1e-6]:
            ts = Series(vals, index=index)

            for m in ['average', 'min', 'max', 'first', 'dense']:
                result = ts.rank(method=m)
                sprank = rankdata(vals, m if m != 'first' else 'ordinal')
                expected = Series(sprank, index=index)

                if LooseVersion(scipy.__version__) >= '0.17.0':
                    expected = expected.astype('float64')
                tm.assert_series_equal(result, expected)

    def test_rank_methods_frame(self):
        tm.skip_if_no_package('scipy', '0.13', 'scipy.stats.rankdata')
        import scipy
        from scipy.stats import rankdata

        xs = np.random.randint(0, 21, (100, 26))
        xs = (xs - 10.0) / 10.0
        cols = [chr(ord('z') - i) for i in range(xs.shape[1])]

        for vals in [xs, xs + 1e6, xs * 1e-6]:
            df = DataFrame(vals, columns=cols)

            for ax in [0, 1]:
                for m in ['average', 'min', 'max', 'first', 'dense']:
                    result = df.rank(axis=ax, method=m)
                    sprank = np.apply_along_axis(
                        rankdata, ax, vals,
                        m if m != 'first' else 'ordinal')
                    sprank = sprank.astype(np.float64)
                    expected = DataFrame(sprank, columns=cols)

                    if LooseVersion(scipy.__version__) >= '0.17.0':
                        expected = expected.astype('float64')
                    tm.assert_frame_equal(result, expected)

    def test_rank_dense_method(self):
        dtypes = ['O', 'f8', 'i8']
        in_out = [([1], [1]),
                  ([2], [1]),
                  ([0], [1]),
                  ([2, 2], [1, 1]),
                  ([1, 2, 3], [1, 2, 3]),
                  ([4, 2, 1], [3, 2, 1],),
                  ([1, 1, 5, 5, 3], [1, 1, 3, 3, 2]),
                  ([-5, -4, -3, -2, -1], [1, 2, 3, 4, 5])]

        for ser, exp in in_out:
            for dtype in dtypes:
                s = Series(ser).astype(dtype)
                result = s.rank(method='dense')
                expected = Series(exp).astype(result.dtype)
                assert_series_equal(result, expected)

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

    def test_rank_object_bug(self):
        # GH 13445

        # smoke tests
        Series([np.nan] * 32).astype(object).rank(ascending=True)
        Series([np.nan] * 32).astype(object).rank(ascending=False)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
