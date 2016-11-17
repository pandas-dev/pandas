import os
import nose

import numpy as np
from pandas.compat import zip

from pandas import Series, Index
import pandas.util.testing as tm
from pandas.util.testing import assertRaisesRegexp
import pandas.core.common as com

from pandas.core.algorithms import quantile
from pandas.tools.tile import cut, qcut
import pandas.tools.tile as tmod


class TestCut(tm.TestCase):

    def test_simple(self):
        data = np.ones(5)
        result = cut(data, 4, labels=False)
        desired = np.array([1, 1, 1, 1, 1])
        tm.assert_numpy_array_equal(result, desired,
                                    check_dtype=False)

    def test_bins(self):
        data = np.array([.2, 1.4, 2.5, 6.2, 9.7, 2.1])
        result, bins = cut(data, 3, retbins=True)

        exp_codes = np.array([0, 0, 0, 1, 2, 0], dtype=np.int8)
        tm.assert_numpy_array_equal(result.codes, exp_codes)
        exp = np.array([0.1905, 3.36666667, 6.53333333, 9.7])
        tm.assert_almost_equal(bins, exp)

    def test_right(self):
        data = np.array([.2, 1.4, 2.5, 6.2, 9.7, 2.1, 2.575])
        result, bins = cut(data, 4, right=True, retbins=True)
        exp_codes = np.array([0, 0, 0, 2, 3, 0, 0], dtype=np.int8)
        tm.assert_numpy_array_equal(result.codes, exp_codes)
        exp = np.array([0.1905, 2.575, 4.95, 7.325, 9.7])
        tm.assert_numpy_array_equal(bins, exp)

    def test_noright(self):
        data = np.array([.2, 1.4, 2.5, 6.2, 9.7, 2.1, 2.575])
        result, bins = cut(data, 4, right=False, retbins=True)
        exp_codes = np.array([0, 0, 0, 2, 3, 0, 1], dtype=np.int8)
        tm.assert_numpy_array_equal(result.codes, exp_codes)
        exp = np.array([0.2, 2.575, 4.95, 7.325, 9.7095])
        tm.assert_almost_equal(bins, exp)

    def test_arraylike(self):
        data = [.2, 1.4, 2.5, 6.2, 9.7, 2.1]
        result, bins = cut(data, 3, retbins=True)
        exp_codes = np.array([0, 0, 0, 1, 2, 0], dtype=np.int8)
        tm.assert_numpy_array_equal(result.codes, exp_codes)
        exp = np.array([0.1905, 3.36666667, 6.53333333, 9.7])
        tm.assert_almost_equal(bins, exp)

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
        s = Series([0, -1, 0, 1, -3], name='x')
        ind = cut(s, [0, 1], labels=False)
        exp = Series([np.nan, np.nan, np.nan, 0, np.nan], name='x')
        tm.assert_series_equal(ind, exp)

    def test_labels(self):
        arr = np.tile(np.arange(0, 1.01, 0.1), 4)

        result, bins = cut(arr, 4, retbins=True)
        ex_levels = Index(['(-0.001, 0.25]', '(0.25, 0.5]', '(0.5, 0.75]',
                           '(0.75, 1]'])
        self.assert_index_equal(result.categories, ex_levels)

        result, bins = cut(arr, 4, retbins=True, right=False)
        ex_levels = Index(['[0, 0.25)', '[0.25, 0.5)', '[0.5, 0.75)',
                           '[0.75, 1.001)'])
        self.assert_index_equal(result.categories, ex_levels)

    def test_cut_pass_series_name_to_factor(self):
        s = Series(np.random.randn(100), name='foo')

        factor = cut(s, 4)
        self.assertEqual(factor.name, 'foo')

    def test_label_precision(self):
        arr = np.arange(0, 0.73, 0.01)

        result = cut(arr, 4, precision=2)
        ex_levels = Index(['(-0.00072, 0.18]', '(0.18, 0.36]',
                           '(0.36, 0.54]', '(0.54, 0.72]'])
        self.assert_index_equal(result.categories, ex_levels)

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

    def test_inf_handling(self):
        data = np.arange(6)
        data_ser = Series(data, dtype='int64')

        result = cut(data, [-np.inf, 2, 4, np.inf])
        result_ser = cut(data_ser, [-np.inf, 2, 4, np.inf])

        ex_categories = Index(['(-inf, 2]', '(2, 4]', '(4, inf]'])

        tm.assert_index_equal(result.categories, ex_categories)
        tm.assert_index_equal(result_ser.cat.categories, ex_categories)
        self.assertEqual(result[5], '(4, inf]')
        self.assertEqual(result[0], '(-inf, 2]')
        self.assertEqual(result_ser[5], '(4, inf]')
        self.assertEqual(result_ser[0], '(-inf, 2]')

    def test_qcut(self):
        arr = np.random.randn(1000)

        labels, bins = qcut(arr, 4, retbins=True)
        ex_bins = quantile(arr, [0, .25, .5, .75, 1.])
        tm.assert_almost_equal(bins, ex_bins)

        ex_levels = cut(arr, ex_bins, include_lowest=True)
        self.assert_categorical_equal(labels, ex_levels)

    def test_qcut_bounds(self):
        arr = np.random.randn(1000)

        factor = qcut(arr, 10, labels=False)
        self.assertEqual(len(np.unique(factor)), 10)

    def test_qcut_specify_quantiles(self):
        arr = np.random.randn(100)

        factor = qcut(arr, [0, .25, .5, .75, 1.])
        expected = qcut(arr, 4)
        tm.assert_categorical_equal(factor, expected)

    def test_qcut_all_bins_same(self):
        assertRaisesRegexp(ValueError, "edges.*unique", qcut,
                           [0, 0, 0, 0, 0, 0, 0, 0, 0, 0], 3)

    def test_cut_out_of_bounds(self):
        arr = np.random.randn(100)

        result = cut(arr, [-1, 0, 1])

        mask = result.codes == -1
        ex_mask = (arr < -1) | (arr > 1)
        self.assert_numpy_array_equal(mask, ex_mask)

    def test_cut_pass_labels(self):
        arr = [50, 5, 10, 15, 20, 30, 70]
        bins = [0, 25, 50, 100]
        labels = ['Small', 'Medium', 'Large']

        result = cut(arr, bins, labels=labels)

        exp = cut(arr, bins)
        exp.categories = labels

        tm.assert_categorical_equal(result, exp)

    def test_qcut_include_lowest(self):
        values = np.arange(10)

        cats = qcut(values, 4)

        ex_levels = ['[0, 2.25]', '(2.25, 4.5]', '(4.5, 6.75]', '(6.75, 9]']
        self.assertTrue((cats.categories == ex_levels).all())

    def test_qcut_nas(self):
        arr = np.random.randn(100)
        arr[:20] = np.nan

        result = qcut(arr, 4)
        self.assertTrue(com.isnull(result[:20]).all())

    def test_label_formatting(self):
        self.assertEqual(tmod._trim_zeros('1.000'), '1')

        # it works
        result = cut(np.arange(11.), 2)

        result = cut(np.arange(11.) / 1e10, 2)

        # #1979, negative numbers

        result = tmod._format_label(-117.9998, precision=3)
        self.assertEqual(result, '-118')
        result = tmod._format_label(117.9998, precision=3)
        self.assertEqual(result, '118')

    def test_qcut_binning_issues(self):
        # #1978, 1979
        path = os.path.join(tm.get_data_path(), 'cut_data.csv')
        arr = np.loadtxt(path)

        result = qcut(arr, 20)

        starts = []
        ends = []
        for lev in result.categories:
            s, e = lev[1:-1].split(',')

            self.assertTrue(s != e)

            starts.append(float(s))
            ends.append(float(e))

        for (sp, sn), (ep, en) in zip(zip(starts[:-1], starts[1:]),
                                      zip(ends[:-1], ends[1:])):
            self.assertTrue(sp < sn)
            self.assertTrue(ep < en)
            self.assertTrue(ep <= sn)

    def test_cut_return_categorical(self):
        from pandas import Categorical
        s = Series([0, 1, 2, 3, 4, 5, 6, 7, 8])
        res = cut(s, 3)
        exp = Series(Categorical.from_codes([0, 0, 0, 1, 1, 1, 2, 2, 2],
                                            ["(-0.008, 2.667]",
                                             "(2.667, 5.333]", "(5.333, 8]"],
                                            ordered=True))
        tm.assert_series_equal(res, exp)

    def test_qcut_return_categorical(self):
        from pandas import Categorical
        s = Series([0, 1, 2, 3, 4, 5, 6, 7, 8])
        res = qcut(s, [0, 0.333, 0.666, 1])
        exp = Series(Categorical.from_codes([0, 0, 0, 1, 1, 1, 2, 2, 2],
                                            ["[0, 2.664]",
                                             "(2.664, 5.328]", "(5.328, 8]"],
                                            ordered=True))
        tm.assert_series_equal(res, exp)

    def test_series_retbins(self):
        # GH 8589
        s = Series(np.arange(4))
        result, bins = cut(s, 2, retbins=True)
        tm.assert_numpy_array_equal(result.cat.codes.values,
                                    np.array([0, 0, 1, 1], dtype=np.int8))
        tm.assert_numpy_array_equal(bins, np.array([-0.003, 1.5, 3]))

        result, bins = qcut(s, 2, retbins=True)
        tm.assert_numpy_array_equal(result.cat.codes.values,
                                    np.array([0, 0, 1, 1], dtype=np.int8))
        tm.assert_numpy_array_equal(bins, np.array([0, 1.5, 3]))

    def test_single_bin(self):
        # issue 14652
        expected = Series([0, 0])

        s = Series([9., 9.])
        result = cut(s, 1, labels=False)
        tm.assert_series_equal(result, expected)

        s = Series([-9., -9.])
        result = cut(s, 1, labels=False)
        tm.assert_series_equal(result, expected)


def curpath():
    pth, _ = os.path.split(os.path.abspath(__file__))
    return pth


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
