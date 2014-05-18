import numpy as np
import unittest
import nose
import pandas.hashtable as _hash
import pandas as pd

class TestFactorizer(unittest.TestCase):

    def test_factorize_nan(self):
        # nan should map to na_sentinel, not reverse_indexer[na_sentinel]
        # rizer.factorize should not raise an exception if na_sentinel indexes
        # outside of reverse_indexer
        key = np.array([1, 2, 1, np.nan], dtype='O')
        rizer = _hash.Factorizer(len(key))
        for na_sentinel in (-1, 20):
            ids = rizer.factorize(key, sort=True, na_sentinel=na_sentinel)
            expected = np.array([0, 1, 0, na_sentinel], dtype='int32')
            self.assertEqual(len(set(key)), len(set(expected)))
            self.assert_(np.array_equal(pd.isnull(key), expected == na_sentinel))

        # nan still maps to na_sentinel when sort=False
        key = np.array([0, np.nan, 1], dtype='O')
        na_sentinel = -1        
        ids = rizer.factorize(key, sort=False, na_sentinel=na_sentinel)
        expected = np.array([ 2, -1,  0], dtype='int32')
        self.assertEqual(len(set(key)), len(set(expected)))
        self.assert_(np.array_equal(pd.isnull(key), expected == na_sentinel))        

    def test_vector_resize(self):
        # Test for memory errors after internal vector
        # reallocations (pull request #7157)

        def _test_vector_resize(htable, uniques, dtype, nvals):
            vals = np.array(np.random.randn(1000), dtype=dtype)
            # get_labels appends to the vector
            htable.get_labels(vals[:nvals], uniques, 0, -1)
            # to_array resizes the vector
            uniques.to_array()
            htable.get_labels(vals, uniques, 0, -1)

        test_cases = [
            (_hash.PyObjectHashTable, _hash.ObjectVector, 'object'),
            (_hash.Float64HashTable, _hash.Float64Vector, 'float64'),
            (_hash.Int64HashTable, _hash.Int64Vector, 'int64')]

        for (tbl, vect, dtype) in test_cases:
            # resizing to empty is a special case
            _test_vector_resize(tbl(), vect(), dtype, 0)
            _test_vector_resize(tbl(), vect(), dtype, 10)


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
