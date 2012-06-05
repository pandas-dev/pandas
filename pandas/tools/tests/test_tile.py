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
        desired = [2, 2, 2, 2, 2]
        assert_equal(result, desired)

    def test_bins(self):
        data = np.array([.2, 1.4, 2.5, 6.2, 9.7, 2.1])
        result, bins = cut(data, 3, labels=False, retbins=True)
        assert_equal(result, [1, 1, 1, 2, 3, 1])
        assert_almost_equal(bins, [ 0.1905, 3.36666667, 6.53333333, 9.7])

    def test_right(self):
        data = np.array([.2, 1.4, 2.5, 6.2, 9.7, 2.1, 2.575])
        result, bins = cut(data, 4, right=True, labels=False, retbins=True)
        assert_equal(result, [1, 1, 1, 3, 4, 1, 1])
        assert_almost_equal(bins, [0.1905, 2.575, 4.95, 7.325, 9.7])

    def test_noright(self):
        data = np.array([.2, 1.4, 2.5, 6.2, 9.7, 2.1, 2.575])
        result, bins = cut(data, 4, right=False, labels=False, retbins=True)
        assert_equal(result, [1, 1, 1, 3, 4, 1, 2])
        assert_almost_equal(bins, [ 0.2, 2.575, 4.95, 7.325, 9.7095])

    def test_arraylike(self):
        data = [.2, 1.4, 2.5, 6.2, 9.7, 2.1]
        result, bins = cut(data, 3, labels=False, retbins=True)
        assert_equal(result, [1, 1, 1, 2, 3, 1])
        assert_almost_equal(bins, [ 0.1905, 3.36666667, 6.53333333, 9.7])

    def test_bins_not_monotonic(self):
        data = [.2, 1.4, 2.5, 6.2, 9.7, 2.1]
        self.assertRaises(ValueError, cut, data, [0.1, 1.5, 1, 10])

    def test_labels(self):
        arr = np.tile(np.arange(0, 1.01, 0.1), 4)

        labels, bins = cut(arr, 4, retbins=True)
        distinct_labels = sorted(unique(labels))
        ex_labels = ['(-0.001, 0.25]', '(0.25, 0.5]', '(0.5, 0.75]',
                     '(0.75, 1]']
        self.assertEqual(distinct_labels, ex_labels)

        labels, bins = cut(arr, 4, retbins=True, right=False)
        distinct_labels = sorted(unique(labels))
        ex_labels = ['[0, 0.25)', '[0.25, 0.5)', '[0.5, 0.75)',
                     '[0.75, 1.001)']
        self.assertEqual(distinct_labels, ex_labels)

    def test_label_precision(self):
        arr = np.arange(0, 0.73, 0.01)

        labels = cut(arr, 4, precision=2)
        distinct_labels = sorted(unique(labels))
        ex_labels = ['(-0.00072, 0.18]', '(0.18, 0.36]', '(0.36, 0.54]',
                     '(0.54, 0.72]']
        self.assertEqual(distinct_labels, ex_labels)

    def test_na_handling(self):
        arr = np.arange(0, 0.75, 0.01)
        arr[::3] = np.nan

        labels = cut(arr, 4)
        ex_labels = np.where(com.isnull(arr), np.nan, labels)

        tm.assert_almost_equal(labels, ex_labels)

        labels = cut(arr, 4, labels=False)
        ex_labels = np.where(com.isnull(arr), np.nan, labels)
        tm.assert_almost_equal(labels, ex_labels)

    def test_qcut(self):
        arr = np.random.randn(1000)

        labels, bins = qcut(arr, 4, retbins=True)

        ex_bins = quantile(arr, [0, .25, .5, .75, 1.])

        assert_almost_equal(bins, ex_bins)

        ex_labels = cut(arr, ex_bins)

        self.assert_(np.array_equal(labels, ex_labels))


if __name__ == '__main__':
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)


