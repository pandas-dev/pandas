# pylint: disable=E1101,E1103,W0232

from datetime import datetime
import unittest
import nose

import numpy as np

from pandas.core.factor import Factor
from pandas.core.index import Index, Int64Index, MultiIndex
from pandas.util.testing import assert_almost_equal

import pandas.util.testing as tm


class TestFactor(unittest.TestCase):

    def setUp(self):
        self.factor = Factor.from_array(['a', 'b', 'b', 'a',
                                         'a', 'c', 'c', 'c'])

    def test_getitem(self):
        self.assertEqual(self.factor[0], 'a')
        self.assertEqual(self.factor[-1], 'c')

        subf = self.factor[[0, 1, 2]]
        tm.assert_almost_equal(subf.labels, [0, 1, 1])

        subf = self.factor[np.asarray(self.factor) == 'c']
        tm.assert_almost_equal(subf.labels, [2, 2, 2])

    def test_constructor_unsortable(self):
        arr = np.array([1, 2, 3, datetime.now()], dtype='O')

        # it works!
        factor = Factor.from_array(arr)

    def test_factor_agg(self):
        import pandas.core.frame as frame

        arr = np.arange(len(self.factor))

        f = np.sum
        agged = frame.factor_agg(self.factor, arr, f)
        labels = self.factor.labels
        for i, idx in enumerate(self.factor.levels):
            self.assertEqual(f(arr[labels == i]), agged[i])

    def test_comparisons(self):
        result = self.factor[self.factor == 'a']
        expected = self.factor[np.asarray(self.factor) == 'a']
        self.assert_(result.equals(expected))

        result = self.factor[self.factor != 'a']
        expected = self.factor[np.asarray(self.factor) != 'a']
        self.assert_(result.equals(expected))

        result = self.factor[self.factor < 'c']
        expected = self.factor[np.asarray(self.factor) < 'c']
        self.assert_(result.equals(expected))

        result = self.factor[self.factor > 'a']
        expected = self.factor[np.asarray(self.factor) > 'a']
        self.assert_(result.equals(expected))

        result = self.factor[self.factor >= 'b']
        expected = self.factor[np.asarray(self.factor) >= 'b']
        self.assert_(result.equals(expected))

        result = self.factor[self.factor <= 'b']
        expected = self.factor[np.asarray(self.factor) <= 'b']
        self.assert_(result.equals(expected))

        n = len(self.factor)

        other = self.factor[np.random.permutation(n)]
        result = self.factor == other
        expected = np.asarray(self.factor) == np.asarray(other)
        self.assert_(np.array_equal(result, expected))

        result = self.factor == 'd'
        expected = np.repeat(False, len(self.factor))
        self.assert_(np.array_equal(result, expected))

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                         # '--with-coverage', '--cover-package=pandas.core'],
                   exit=False)


