# pylint: disable-msg=E1101,W0612

import nose  # noqa
import numpy as np
import pandas as pd
import pandas.util.testing as tm


class TestSparseSeriesIndexing(tm.TestCase):

    _multiprocess_can_split_ = True

    def test_loc(self):
        orig = pd.Series([1, np.nan, np.nan, 3, np.nan])
        sparse = orig.to_sparse()

        self.assertEqual(sparse.loc[0], 1)
        self.assertTrue(np.isnan(sparse.loc[1]))

        result = sparse.loc[[1, 3, 4]]
        exp = orig.loc[[1, 3, 4]].to_sparse()
        tm.assert_sp_series_equal(result, exp)

        # exceeds the bounds
        result = sparse.loc[[1, 3, 4, 5]]
        exp = orig.loc[[1, 3, 4, 5]].to_sparse()
        tm.assert_sp_series_equal(result, exp)
        # padded with NaN
        self.assertTrue(np.isnan(result[-1]))

        # dense array
        result = sparse.loc[orig % 2 == 1]
        exp = orig.loc[orig % 2 == 1].to_sparse()
        tm.assert_sp_series_equal(result, exp)

        # sparse array (actuary it coerces to normal Series)
        result = sparse.loc[sparse % 2 == 1]
        exp = orig.loc[orig % 2 == 1].to_sparse()
        tm.assert_sp_series_equal(result, exp)

    def test_loc_index(self):
        orig = pd.Series([1, np.nan, np.nan, 3, np.nan], index=list('ABCDE'))
        sparse = orig.to_sparse()

        self.assertEqual(sparse.loc['A'], 1)
        self.assertTrue(np.isnan(sparse.loc['B']))

        result = sparse.loc[['A', 'C', 'D']]
        exp = orig.loc[['A', 'C', 'D']].to_sparse()
        tm.assert_sp_series_equal(result, exp)

        # dense array
        result = sparse.loc[orig % 2 == 1]
        exp = orig.loc[orig % 2 == 1].to_sparse()
        tm.assert_sp_series_equal(result, exp)

        # sparse array (actuary it coerces to normal Series)
        result = sparse.loc[sparse % 2 == 1]
        exp = orig.loc[orig % 2 == 1].to_sparse()
        tm.assert_sp_series_equal(result, exp)

    def test_loc_slice(self):
        orig = pd.Series([1, np.nan, np.nan, 3, np.nan])
        sparse = orig.to_sparse()
        tm.assert_sp_series_equal(sparse.loc[2:], orig.loc[2:].to_sparse())

    def test_iloc(self):
        orig = pd.Series([1, np.nan, np.nan, 3, np.nan])
        sparse = orig.to_sparse()

        self.assertEqual(sparse.iloc[3], 3)
        self.assertTrue(np.isnan(sparse.iloc[2]))

        result = sparse.iloc[[1, 3, 4]]
        exp = orig.iloc[[1, 3, 4]].to_sparse()
        tm.assert_sp_series_equal(result, exp)

        with tm.assertRaises(IndexError):
            sparse.iloc[[1, 3, 5]]

    def test_iloc_slice(self):
        orig = pd.Series([1, np.nan, np.nan, 3, np.nan])
        sparse = orig.to_sparse()
        tm.assert_sp_series_equal(sparse.iloc[2:], orig.iloc[2:].to_sparse())
