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
            self.assertTrue(np.array_equal(pd.isnull(key), expected == na_sentinel))

        # nan still maps to na_sentinel when sort=False
        key = np.array([0, np.nan, 1], dtype='O')
        na_sentinel = -1
        ids = rizer.factorize(key, sort=False, na_sentinel=na_sentinel)
        expected = np.array([ 2, -1,  0], dtype='int32')
        self.assertEqual(len(set(key)), len(set(expected)))
        self.assertTrue(np.array_equal(pd.isnull(key), expected == na_sentinel))

if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
