import os
import nose
import unittest

import numpy as np

from pandas import DataFrame, Series, unique
import pandas.util.testing as tm
from pandas.util.testing import assertRaisesRegexp
import pandas.core.common as com

from pandas.core.algorithms import quantile
from pandas.tools.tile import cut, qcut
import pandas.tools.tile as tmod

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
        assert_almost_equal(bins, [0.1905, 3.36666667, 6.53333333, 9.7])

    def test_right(self):
        data = np.array([.2, 1.4, 2.5, 6.2, 9.7, 2.1, 2.575])
        result, bins = cut(data, 4, right=True, retbins=True)
        assert_equal(result.labels, [0, 0, 0, 2, 3, 0, 0])
        assert_almost_equal(bins, [0.1905, 2.575, 4.95, 7.325, 9.7])

    def test_noright(self):
        data = np.array([.2, 1.4, 2.5, 6.2, 9.7, 2.1, 2.575])
        result, bins = cut(data, 4, right=False, retbins=True)
        assert_equal(result.labels, [0, 0, 0, 2, 3, 0, 1])
        assert_almost_equal(bins, [0.2, 2.575, 4.95, 7.325, 9.7095])

    def test_arraylike(self):
        data = [.2, 1.4, 2.5, 6.2, 9.7, 2.1]
        result, bins = cut(data, 3, retbins=True)
        assert_equal(result.labels, [0, 0, 0, 1, 2, 0])
        assert_almost_equal(bins, [0.1905, 3.36666667, 6.53333333, 9.7])

    def test_bins_not_monotonic(self):
        data = [.2, 1.4, 2.5, 6.2, 9.7, 2.1]
        self.assertRaises(ValueError, cut, data, [0.1, 1.5, 1, 10])

    def test_wrong_num_labels(self):
        data = [.2, 1.4, 2.5, 6.2, 9.7, 2.1]
        self.assertRaises(ValueError, cut, data, [0, 1, 10],
                          labels=['foo', 'bar', 'baz'])

    def test_cut_corner(self):
        # h3h
        self.assertRaises(ValueError, cut, [], 2)

        self.assertRaises(ValueError, cut, [1, 2, 3], 0.5)

    def test_cut_out_of_range_more(self):
        # #1511
        s = Series([0, -1, 0, 1, -3])
        ind = cut(s, [0, 1], labels=False)
        exp = [np.nan, np.nan, np.nan, 0, np.nan]
        assert_almost_equal(ind, exp)

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
        assert_almost_equal(bins, ex_bins)

        ex_levels = cut(arr, ex_bins, include_lowest=True)
        self.assert_(np.array_equal(labels, ex_levels))

    def test_qcut_bounds(self):
        arr = np.random.randn(1000)

        factor = qcut(arr, 10, labels=False)
        self.assert_(len(np.unique(factor)) == 10)

    def test_qcut_specify_quantiles(self):
        arr = np.random.randn(100)

        factor = qcut(arr, [0, .25, .5, .75, 1.])
        expected = qcut(arr, 4)
        self.assert_(factor.equals(expected))

    def test_qcut_all_bins_same(self):
        assertRaisesRegexp(ValueError, "edges.*unique", qcut, [0,0,0,0,0,0,0,0,0,0], 3)

    def test_cut_out_of_bounds(self):
        arr = np.random.randn(100)

        result = cut(arr, [-1, 0, 1])

        mask = result.labels == -1
        ex_mask = (arr < -1) | (arr > 1)
        self.assert_(np.array_equal(mask, ex_mask))

    def test_cut_pass_labels(self):
        arr = [50, 5, 10, 15, 20, 30, 70]
        bins = [0, 25, 50, 100]
        labels = ['Small', 'Medium', 'Large']

        result = cut(arr, bins, labels=labels)

        exp = cut(arr, bins)
        exp.levels = labels

        self.assert_(result.equals(exp))

    def test_qcut_include_lowest(self):
        values = np.arange(10)

        cats = qcut(values, 4)

        ex_levels = ['[0, 2.25]', '(2.25, 4.5]', '(4.5, 6.75]', '(6.75, 9]']
        self.assert_((cats.levels == ex_levels).all())

    def test_qcut_nas(self):
        arr = np.random.randn(100)
        arr[:20] = np.nan

        result = qcut(arr, 4)
        self.assert_(com.isnull(result[:20]).all())

    def test_label_formatting(self):
        self.assertEquals(tmod._trim_zeros('1.000'), '1')

        # it works
        result = cut(np.arange(11.), 2)

        result = cut(np.arange(11.) / 1e10, 2)

        # #1979, negative numbers

        result = tmod._format_label(-117.9998, precision=3)
        self.assertEquals(result, '-118')
        result = tmod._format_label(117.9998, precision=3)
        self.assertEquals(result, '118')

    def test_qcut_binning_issues(self):
        # #1978, 1979
        path = os.path.join(curpath(), 'cut_data.csv')

        arr = np.loadtxt(path)

        result = qcut(arr, 20)

        starts = []
        ends = []
        for lev in result.levels:
            s, e = lev[1:-1].split(',')

            self.assertTrue(s != e)

            starts.append(float(s))
            ends.append(float(e))

        for (sp, sn), (ep, en) in zip(zip(starts[:-1], starts[1:]),
                                      zip(ends[:-1], ends[1:])):
            self.assertTrue(sp < sn)
            self.assertTrue(ep < en)
            self.assertTrue(ep <= sn)


def curpath():
    pth, _ = os.path.split(os.path.abspath(__file__))
    return pth


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
