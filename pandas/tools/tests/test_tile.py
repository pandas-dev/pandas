import nose
import unittest

import numpy as np

from pandas import DataFrame, Series, unique
import pandas.util.testing as tm
import pandas.core.common as com

from pandas.core.algorithms import quantile
from pandas.tools.tile import cut, qcut

from numpy.testing import assert_equal, assert_almost_equal

class TestCut(unittest.TestCase):

    def test_simple(self):
        data = np.ones(5)
        result = cut(data, 4, labels=False)
        desired = [1, 1, 1, 1, 1]
        assert_equal(result, desired)

    def test_bins(self):
        data = np.array([.2, 1.4, 2.5, 6.2, 9.7, 2.1])
        result, bins = cut(data, 3, retbins=True)
        assert_equal(result.labels, [0, 0, 0, 1, 2, 0])
        assert_almost_equal(bins, [ 0.1905, 3.36666667, 6.53333333, 9.7])

    def test_right(self):
        data = np.array([.2, 1.4, 2.5, 6.2, 9.7, 2.1, 2.575])
        result, bins = cut(data, 4, right=True, retbins=True)
        assert_equal(result.labels, [0, 0, 0, 2, 3, 0, 0])
        assert_almost_equal(bins, [0.1905, 2.575, 4.95, 7.325, 9.7])

    def test_noright(self):
        data = np.array([.2, 1.4, 2.5, 6.2, 9.7, 2.1, 2.575])
        result, bins = cut(data, 4, right=False, retbins=True)
        assert_equal(result.labels, [0, 0, 0, 2, 3, 0, 1])
        assert_almost_equal(bins, [ 0.2, 2.575, 4.95, 7.325, 9.7095])

    def test_arraylike(self):
        data = [.2, 1.4, 2.5, 6.2, 9.7, 2.1]
        result, bins = cut(data, 3, retbins=True)
        assert_equal(result.labels, [0, 0, 0, 1, 2, 0])
        assert_almost_equal(bins, [ 0.1905, 3.36666667, 6.53333333, 9.7])

    def test_bins_not_monotonic(self):
        data = [.2, 1.4, 2.5, 6.2, 9.7, 2.1]
        self.assertRaises(ValueError, cut, data, [0.1, 1.5, 1, 10])

    def test_labels(self):
        arr = np.tile(np.arange(0, 1.01, 0.1), 4)

        result, bins = cut(arr, 4, retbins=True)
        ex_levels = ['(-0.001, 0.25]', '(0.25, 0.5]', '(0.5, 0.75]',
                     '(0.75, 1]']
        self.assert_(np.array_equal(result.levels, ex_levels))

        result, bins = cut(arr, 4, retbins=True, right=False)
        ex_levels = ['[0, 0.25)', '[0.25, 0.5)', '[0.5, 0.75)',
                     '[0.75, 1.001)']
        self.assert_(np.array_equal(result.levels, ex_levels))

    def test_cut_pass_series_name_to_factor(self):
        s = Series(np.random.randn(100), name='foo')

        factor = cut(s, 4)
        self.assertEquals(factor.name, 'foo')

    def test_label_precision(self):
        arr = np.arange(0, 0.73, 0.01)

        result = cut(arr, 4, precision=2)
        ex_levels = ['(-0.00072, 0.18]', '(0.18, 0.36]', '(0.36, 0.54]',
                     '(0.54, 0.72]']
        self.assert_(np.array_equal(result.levels, ex_levels))

    def test_na_handling(self):
        arr = np.arange(0, 0.75, 0.01)
        arr[::3] = np.nan

        result = cut(arr, 4)

        result_arr = np.asarray(result)

        ex_arr = np.where(com.isnull(arr), np.nan, result_arr)

        tm.assert_almost_equal(result_arr, ex_arr)

        result = cut(arr, 4, labels=False)
        ex_result = np.where(com.isnull(arr), np.nan, result)
        tm.assert_almost_equal(result, ex_result)

    def test_qcut(self):
        arr = np.random.randn(1000)

        labels, bins = qcut(arr, 4, retbins=True)
        ex_bins = quantile(arr, [0, .25, .5, .75, 1.])
        ex_bins[0] -= (arr.max() - arr.min()) * 0.001
        assert_almost_equal(bins, ex_bins)

        ex_levels = cut(arr, ex_bins)
        self.assert_(np.array_equal(labels, ex_levels))

    def test_qcut_bounds(self):
        np.random.seed(12345)
        arr = np.random.randn(1000)

        factor = qcut(arr, 10, labels=False)
        self.assert_(len(np.unique(factor)) == 10)

    def test_qcut_specify_quantiles(self):
        arr = np.random.randn(100)

        factor = qcut(arr, [0, .25, .5, .75, 1.])
        expected = qcut(arr, 4)
        self.assert_(factor.equals(expected))

    def test_cut_out_of_bounds(self):
        np.random.seed(12345)

        arr = np.random.randn(100)
        self.assertRaises(ValueError, cut, arr, [-1, 0, 1])

        arr = np.where(arr < -1, 0, arr)
        self.assertRaises(ValueError, cut, arr, [-1, 0, 1])

    def test_cut_pass_labels(self):
        arr = [50, 5, 10, 15, 20, 30, 70]
        bins = [0, 25, 50, 100]
        labels = ['Small', 'Medium', 'Large']

        result = cut(arr, bins, labels=labels)

        exp = cut(arr, bins)
        exp.levels = labels

        self.assert_(result.equals(exp))

if __name__ == '__main__':
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)


