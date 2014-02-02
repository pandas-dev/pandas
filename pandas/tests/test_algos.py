from pandas.compat import range

import numpy as np

from pandas.core.api import Series, Categorical
import pandas as pd

import pandas.core.algorithms as algos
import pandas.util.testing as tm

class TestMatch(tm.TestCase):
    _multiprocess_can_split_ = True

    def test_ints(self):
        values = np.array([0, 2, 1])
        to_match = np.array([0, 1, 2, 2, 0, 1, 3, 0])

        result = algos.match(to_match, values)
        expected = np.array([0, 2, 1, 1, 0, 2, -1, 0])
        self.assert_(np.array_equal(result, expected))

        result = Series(algos.match(to_match, values, np.nan))
        expected = Series(np.array([0, 2, 1, 1, 0, 2, np.nan, 0]))
        tm.assert_series_equal(result,expected)

        s = pd.Series(np.arange(5),dtype=np.float32)
        result = algos.match(s, [2,4])
        expected = np.array([-1, -1, 0, -1, 1])
        self.assert_(np.array_equal(result, expected))

        result = Series(algos.match(s, [2,4], np.nan))
        expected = Series(np.array([np.nan, np.nan, 0, np.nan, 1]))
        tm.assert_series_equal(result,expected)

    def test_strings(self):
        values = ['foo', 'bar', 'baz']
        to_match = ['bar', 'foo', 'qux', 'foo', 'bar', 'baz', 'qux']

        result = algos.match(to_match, values)
        expected = np.array([1, 0, -1, 0, 1, 2, -1])
        self.assert_(np.array_equal(result, expected))

        result = Series(algos.match(to_match, values, np.nan))
        expected = Series(np.array([1, 0, np.nan, 0, 1, 2, np.nan]))
        tm.assert_series_equal(result,expected)

class TestFactorize(tm.TestCase):
    _multiprocess_can_split_ = True

    def test_basic(self):

        labels, uniques = algos.factorize(['a', 'b', 'b', 'a',
                                           'a', 'c', 'c', 'c'])
        self.assert_(np.array_equal(labels, np.array([ 0, 1, 1, 0, 0, 2, 2, 2],dtype=np.int64)))
        self.assert_(np.array_equal(uniques, np.array(['a','b','c'], dtype=object)))

        labels, uniques = algos.factorize(['a', 'b', 'b', 'a',
                                           'a', 'c', 'c', 'c'], sort=True)
        self.assert_(np.array_equal(labels, np.array([ 0, 1, 1, 0, 0, 2, 2, 2],dtype=np.int64)))
        self.assert_(np.array_equal(uniques, np.array(['a','b','c'], dtype=object)))

        labels, uniques = algos.factorize(list(reversed(range(5))))
        self.assert_(np.array_equal(labels, np.array([0, 1, 2, 3, 4], dtype=np.int64)))
        self.assert_(np.array_equal(uniques, np.array([ 4, 3, 2, 1, 0],dtype=np.int64)))

        labels, uniques = algos.factorize(list(reversed(range(5))), sort=True)
        self.assert_(np.array_equal(labels, np.array([ 4, 3, 2, 1, 0],dtype=np.int64)))
        self.assert_(np.array_equal(uniques, np.array([0, 1, 2, 3, 4], dtype=np.int64)))

        labels, uniques = algos.factorize(list(reversed(np.arange(5.))))
        self.assert_(np.array_equal(labels, np.array([0., 1., 2., 3., 4.], dtype=np.float64)))
        self.assert_(np.array_equal(uniques, np.array([ 4, 3, 2, 1, 0],dtype=np.int64)))

        labels, uniques = algos.factorize(list(reversed(np.arange(5.))), sort=True)
        self.assert_(np.array_equal(labels, np.array([ 4, 3, 2, 1, 0],dtype=np.int64)))
        self.assert_(np.array_equal(uniques, np.array([0., 1., 2., 3., 4.], dtype=np.float64)))

    def test_mixed(self):

        # doc example reshaping.rst
        x = Series(['A', 'A', np.nan, 'B', 3.14, np.inf])
        labels, uniques = algos.factorize(x)

        self.assert_(np.array_equal(labels, np.array([ 0,  0, -1,  1,  2,  3],dtype=np.int64)))
        self.assert_(np.array_equal(uniques, np.array(['A', 'B', 3.14, np.inf], dtype=object)))

        labels, uniques = algos.factorize(x, sort=True)
        self.assert_(np.array_equal(labels, np.array([ 2,  2, -1,  3,  0,  1],dtype=np.int64)))
        self.assert_(np.array_equal(uniques, np.array([3.14, np.inf, 'A', 'B'], dtype=object)))

    def test_datelike(self):

        # M8
        v1 = pd.Timestamp('20130101 09:00:00.00004')
        v2 = pd.Timestamp('20130101')
        x = Series([v1,v1,v1,v2,v2,v1])
        labels, uniques = algos.factorize(x)
        self.assert_(np.array_equal(labels, np.array([ 0,0,0,1,1,0],dtype=np.int64)))
        self.assert_(np.array_equal(uniques, np.array([v1.value,v2.value],dtype='M8[ns]')))

        labels, uniques = algos.factorize(x, sort=True)
        self.assert_(np.array_equal(labels, np.array([ 1,1,1,0,0,1],dtype=np.int64)))
        self.assert_(np.array_equal(uniques, np.array([v2.value,v1.value],dtype='M8[ns]')))

        # period
        v1 = pd.Period('201302',freq='M')
        v2 = pd.Period('201303',freq='M')
        x = Series([v1,v1,v1,v2,v2,v1])

        # periods are not 'sorted' as they are converted back into an index
        labels, uniques = algos.factorize(x)
        self.assert_(np.array_equal(labels, np.array([ 0,0,0,1,1,0],dtype=np.int64)))
        self.assert_(np.array_equal(uniques, np.array([v1, v2],dtype=object)))

        labels, uniques = algos.factorize(x,sort=True)
        self.assert_(np.array_equal(labels, np.array([ 0,0,0,1,1,0],dtype=np.int64)))
        self.assert_(np.array_equal(uniques, np.array([v1, v2],dtype=object)))

class TestUnique(tm.TestCase):
    _multiprocess_can_split_ = True

    def test_ints(self):
        arr = np.random.randint(0, 100, size=50)

        result = algos.unique(arr)
        tm.assert_isinstance(result, np.ndarray)

    def test_objects(self):
        arr = np.random.randint(0, 100, size=50).astype('O')

        result = algos.unique(arr)
        tm.assert_isinstance(result, np.ndarray)

    def test_object_refcount_bug(self):
        lst = ['A', 'B', 'C', 'D', 'E']
        for i in range(1000):
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

class TestValueCounts(tm.TestCase):
    _multiprocess_can_split_ = True

    def test_value_counts(self):
        from pandas.tools.tile import cut

        arr = np.random.randn(4)
        factor = cut(arr, 4)

        tm.assert_isinstance(factor, Categorical)

        result = algos.value_counts(factor)
        expected = algos.value_counts(np.asarray(factor))
        tm.assert_series_equal(result, expected)

    def test_value_counts_bins(self):
        s = [1, 2, 3, 4]
        result = algos.value_counts(s, bins=1)
        self.assertEqual(result.tolist(), [4])
        self.assertEqual(result.index[0], 0.997)

        result = algos.value_counts(s, bins=2, sort=False)
        self.assertEqual(result.tolist(), [2, 2])
        self.assertEqual(result.index[0], 0.997)
        self.assertEqual(result.index[1], 2.5)

    def test_value_counts_dtypes(self):
        result = algos.value_counts([1, 1.])
        self.assertEqual(len(result), 1)

        result = algos.value_counts([1, 1.], bins=1)
        self.assertEqual(len(result), 1)

        result = algos.value_counts(Series([1, 1., '1']))  # object
        self.assertEqual(len(result), 2)

        self.assertRaises(TypeError, lambda s: algos.value_counts(s, bins=1), ['1', 1])


def test_quantile():
    s = Series(np.random.randn(100))

    result = algos.quantile(s, [0, .25, .5, .75, 1.])
    expected = algos.quantile(s.values, [0, .25, .5, .75, 1.])
    tm.assert_almost_equal(result, expected)

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
