import nose
import unittest

from numpy import nan

from pandas import *

from pandas.util.compat import product
from pandas.util.testing import (assert_frame_equal,
                                 assert_almost_equal)

class TestStats(unittest.TestCase):

    def test_rank_tie_methods(self):
        s = Series([1, 3, 4, 2, nan, 2, 1, 5, nan, 3])

        # c(1, 3, 4, 2, NA, 2, 1, 5, NA, 3)

        ranks = s.rank()
        expected = np.array([1.5, 5.5, 7.0, 3.5, nan,
                             3.5, 1.5, 8.0, nan, 5.5])
        assert_almost_equal(ranks, expected)

        ranks = s.rank(method='min')
        expected = np.array([1, 5, 7, 3, nan, 3, 1, 8, nan, 5])
        assert_almost_equal(ranks, expected)

        ranks = s.rank(method='max')
        expected = np.array([2, 6, 7, 4, nan, 4, 2, 8, nan, 6])
        assert_almost_equal(ranks, expected)

        ranks = s.rank(method='first')
        expected = np.array([1, 5, 7, 3, nan, 4, 2, 8, nan, 6])
        assert_almost_equal(ranks, expected)

    def test_rank_2d_tie_methods(self):
        s = Series([1, 3, 4, 2, nan, 2, 1, 5, nan, 3])
        df = DataFrame({'A': s, 'B': s})

        def _check2d(df, expected, method='average', axis=0):
            exp_df = DataFrame({'A': expected, 'B': expected})

            if axis == 1:
                df = df.T
                exp_df = exp_df.T

            result = df.rank(method=method, axis=axis)
            assert_frame_equal(result, exp_df)

        results = {
            'average': np.array([1.5, 5.5, 7.0, 3.5, nan,
                                 3.5, 1.5, 8.0, nan, 5.5]),
            'min': np.array([1, 5, 7, 3, nan, 3, 1, 8, nan, 5]),
            'max': np.array([2, 6, 7, 4, nan, 4, 2, 8, nan, 6]),
            'first': np.array([1, 5, 7, 3, nan, 4, 2, 8, nan, 6])
        }

        dtypes = [None]

        for method, axis, dtype in product(results, [0, 1], dtypes):
            frame = df if dtype is None else df.astype(dtype)
            _check2d(frame, results[method], method=method, axis=axis)

if __name__ == '__main__':
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

