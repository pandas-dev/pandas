import unittest

import numpy as np

from pandas.core.api import Series
import pandas as pd

import pandas.core.algorithms as algos
import pandas.util.testing as tm


class TestMatch(unittest.TestCase):
    _multiprocess_can_split_ = True

    def test_ints(self):
        values = np.array([0, 2, 1])
        to_match = np.array([0, 1, 2, 2, 0, 1, 3, 0])

        result = algos.match(to_match, values)
        expected = np.array([0, 2, 1, 1, 0, 2, -1, 0])
        self.assert_(np.array_equal(result, expected))

    def test_strings(self):
        values = ['foo', 'bar', 'baz']
        to_match = ['bar', 'foo', 'qux', 'foo', 'bar', 'baz', 'qux']

        result = algos.match(to_match, values)
        expected = np.array([1, 0, -1, 0, 1, 2, -1])
        self.assert_(np.array_equal(result, expected))


class TestUnique(unittest.TestCase):
    _multiprocess_can_split_ = True

    def test_ints(self):
        arr = np.random.randint(0, 100, size=50)

        result = algos.unique(arr)
        self.assert_(isinstance(result, np.ndarray))

    def test_objects(self):
        arr = np.random.randint(0, 100, size=50).astype('O')

        result = algos.unique(arr)
        self.assert_(isinstance(result, np.ndarray))

    def test_object_refcount_bug(self):
        lst = ['A', 'B', 'C', 'D', 'E']
        for i in xrange(1000):
            len(algos.unique(lst))

    def test_on_index_object(self):
        mindex = pd.MultiIndex.from_arrays([np.arange(5).repeat(5),
                                            np.tile(np.arange(5), 5)])
        mindex = mindex.repeat(2)

        result = pd.unique(mindex)
        result.sort()

        expected = mindex.values
        expected.sort()

        tm.assert_almost_equal(result, expected)


def test_quantile():
    s = Series(np.random.randn(100))

    result = algos.quantile(s, [0, .25, .5, .75, 1.])
    expected = algos.quantile(s.values, [0, .25, .5, .75, 1.])
    tm.assert_almost_equal(result, expected)

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
