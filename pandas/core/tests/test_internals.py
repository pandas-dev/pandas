import unittest

import numpy as np

from pandas import Index
from pandas.core.internals import *

from pandas.util.testing import (assert_almost_equal, randn)

class TestBlock(unittest.TestCase):

    def setUp(self):
        n = 10
        floats = np.repeat(np.atleast_2d(np.arange(3.)), n, axis=0)
        self.fblock = make_block(floats, [0, 2, 4],
                                 ['a', 'b', 'c', 'd', 'e'])

    def test_merge(self):
        pass

    def test_copy(self):
        pass

    def test_columns(self):
        cols = self.fblock.columns
        self.assert_(np.array_equal(cols, ['a', 'c', 'e']))

        cols2 = self.fblock.columns
        self.assert_(cols is cols2)

    def test_assign_ref_columns(self):
        self.fblock.ref_columns = ['foo', 'bar', 'baz', 'quux', 'hi']
        self.assert_(np.array_equal(self.fblock.columns,
                                    ['foo', 'baz', 'hi']))

    def test_reindex_index(self):
        pass

    def test_reindex_columns_from(self):
        new_cols = Index(['e', 'b', 'c', 'f'])
        reindexed = self.fblock.reindex_columns_from(new_cols)
        assert_almost_equal(reindexed.ref_locs, [0, 2])
        self.assertEquals(reindexed.values.shape[1], 2)
        self.assert_((reindexed.values[:, 0] == 2).all())
        self.assert_((reindexed.values[:, 1] == 1).all())

    def test_insert(self):
        pass

    def test_delete(self):
        newb = self.fblock.delete('a')
        assert_almost_equal(newb.ref_locs, [2, 4])
        self.assert_((newb.values[:, 0] == 1).all())

        newb = self.fblock.delete('c')
        assert_almost_equal(newb.ref_locs, [0, 4])
        self.assert_((newb.values[:, 1] == 2).all())

        newb = self.fblock.delete('e')
        assert_almost_equal(newb.ref_locs, [0, 2])
        self.assert_((newb.values[:, 1] == 1).all())

        self.assertRaises(Exception, self.fblock.delete, 'b')

    def test_get(self):
        pass

    def test_set(self):
        pass

    def test_fillna(self):
        pass

    def test_repr(self):
        pass


class TestBlockManager(unittest.TestCase):

    def test_interleave(self):
        pass

    def test_consolidate(self):
        pass

    def test_xs(self):
        pass


if __name__ == '__main__':
    # unittest.main()
    import nose
    # nose.runmodule(argv=[__file__,'-vvs','-x', '--pdb-failure'],
    #                exit=False)
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

