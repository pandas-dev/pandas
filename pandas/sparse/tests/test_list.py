from pandas.compat import range
import unittest

from numpy import nan
import numpy as np

from pandas.sparse.api import SparseList, SparseArray
from pandas.util.testing import assert_almost_equal

from .test_sparse import assert_sp_array_equal


def assert_sp_list_equal(left, right):
    assert_sp_array_equal(left.to_array(), right.to_array())


class TestSparseList(unittest.TestCase):

    _multiprocess_can_split_ = True

    def setUp(self):
        self.na_data = np.array([nan, nan, 1, 2, 3, nan, 4, 5, nan, 6])
        self.zero_data = np.array([0, 0, 1, 2, 3, 0, 4, 5, 0, 6])

    def test_constructor(self):
        lst1 = SparseList(self.na_data[:5])
        exp = SparseList()
        exp.append(self.na_data[:5])
        assert_sp_list_equal(lst1, exp)

    def test_len(self):
        arr = self.na_data
        splist = SparseList()
        splist.append(arr[:5])
        self.assertEquals(len(splist), 5)
        splist.append(arr[5])
        self.assertEquals(len(splist), 6)
        splist.append(arr[6:])
        self.assertEquals(len(splist), 10)

    def test_append_na(self):
        arr = self.na_data
        splist = SparseList()
        splist.append(arr[:5])
        splist.append(arr[5])
        splist.append(arr[6:])

        sparr = splist.to_array()
        assert_sp_array_equal(sparr, SparseArray(arr))

    def test_append_zero(self):
        arr = self.zero_data
        splist = SparseList(fill_value=0)
        splist.append(arr[:5])
        splist.append(arr[5])
        splist.append(arr[6:])

        sparr = splist.to_array()
        assert_sp_array_equal(sparr, SparseArray(arr, fill_value=0))

    def test_consolidate(self):
        arr = self.na_data
        exp_sparr = SparseArray(arr)

        splist = SparseList()
        splist.append(arr[:5])
        splist.append(arr[5])
        splist.append(arr[6:])

        consol = splist.consolidate(inplace=False)
        self.assert_(consol.nchunks == 1)
        self.assert_(splist.nchunks == 3)
        assert_sp_array_equal(consol.to_array(), exp_sparr)

        splist.consolidate()
        self.assert_(splist.nchunks == 1)
        assert_sp_array_equal(splist.to_array(), exp_sparr)

    def test_copy(self):
        arr = self.na_data
        exp_sparr = SparseArray(arr)

        splist = SparseList()
        splist.append(arr[:5])
        splist.append(arr[5])

        cp = splist.copy()
        cp.append(arr[6:])
        self.assertEquals(splist.nchunks, 2)
        assert_sp_array_equal(cp.to_array(), exp_sparr)

    def test_getitem(self):
        arr = self.na_data
        splist = SparseList()
        splist.append(arr[:5])
        splist.append(arr[5])
        splist.append(arr[6:])

        for i in range(len(arr)):
            assert_almost_equal(splist[i], arr[i])
            assert_almost_equal(splist[-i], arr[-i])


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
