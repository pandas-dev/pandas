import unittest

import numpy as np

from pandas import Index
from pandas.core.internals import *

from pandas.util.testing import (assert_almost_equal, randn)

def assert_block_equal(left, right):
    assert_almost_equal(left.values, right.values)
    assert(left.columns.equals(right.columns))
    assert(left.ref_columns.equals(right.ref_columns))

def get_float_mat(n, k):
    return np.repeat(np.atleast_2d(np.arange(k, dtype=float)), n, axis=0)

TEST_COLS = ['a', 'b', 'c', 'd', 'e', 'f']
N = 10

def get_float_ex():
    floats = get_float_mat(N, 3)
    return make_block(floats, ['a', 'c', 'e'], TEST_COLS)

def get_obj_ex():
    mat = np.empty((N, 2), dtype=object)
    mat[:, 0] = 'foo'
    mat[:, 1] = 'bar'
    return make_block(mat, ['b', 'd'], TEST_COLS)

def get_bool_ex():
    mat = np.ones((N, 1), dtype=bool)
    return make_block(mat, ['f'], TEST_COLS)

class TestBlock(unittest.TestCase):

    def setUp(self):
        self.fblock = get_float_ex()
        self.oblock = get_obj_ex()
        self.bool_block = get_bool_ex()

    def test_constructor(self):
        pass

    def test_ref_locs(self):
        assert_almost_equal(self.fblock.ref_locs, [0, 2, 4])

    def test_attrs(self):
        self.assert_(self.fblock.shape == self.fblock.values.shape)
        self.assert_(self.fblock.dtype == self.fblock.values.dtype)

    def test_merge(self):
        avals = randn(10, 2)
        bvals = randn(10, 2)

        ref_cols = ['e', 'a', 'b', 'd', 'f']

        ablock = make_block(avals, ['e', 'b'], ref_cols)
        bblock = make_block(bvals, ['a', 'd'], ref_cols)
        merged = ablock.merge(bblock)
        exvals = np.hstack((avals, bvals))
        excols = ['e', 'b', 'a', 'd']
        eblock = make_block(exvals, excols, ref_cols)
        eblock = eblock.reindex_columns_from(ref_cols)
        assert_block_equal(merged, eblock)

        # TODO: merge with mixed type?

    def test_copy(self):
        cop = self.fblock.copy()
        self.assert_(cop is not self.fblock)
        assert_block_equal(self.fblock, cop)

    def test_columns(self):
        cols = self.fblock.columns
        self.assert_(np.array_equal(cols, ['a', 'c', 'e']))

        cols2 = self.fblock.columns
        self.assert_(cols is cols2)

    def test_assign_ref_columns(self):
        new_cols = Index(['foo', 'bar', 'baz', 'quux', 'hi'])
        self.fblock.set_ref_columns(new_cols)
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

    def test_reindex_cast(self):
        pass

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

    def setUp(self):
        self.blocks = [get_float_ex(),
                       get_obj_ex(),
                       get_bool_ex()]
        self.mgr = BlockManager.from_blocks(self.blocks, np.arange(N))

    def test_attrs(self):
        self.assertEquals(self.mgr.nblocks, len(self.mgr.blocks))
        self.assertEquals(len(self.mgr), len(self.mgr.index))

    def test_contains(self):
        self.assert_('a' in self.mgr)
        self.assert_('g' not in self.mgr)

    def test_get(self):
        pass

    def test_set(self):
        pass

    def test_set_change_dtype(self):
        self.mgr.set('g', np.zeros(N, dtype=bool))

        self.mgr.set('g', np.repeat('foo', N))
        self.assert_(self.mgr.get('g').dtype == np.object_)

        mgr2 = self.mgr.consolidate()
        mgr2.set('g', np.repeat('foo', N))
        self.assert_(mgr2.get('g').dtype == np.object_)

    def test_copy(self):
        pass

    def test_as_matrix(self):
        pass

    def test_xs(self):
        pass

    def test_from_blocks(self):
        self.assert_(np.array_equal(self.mgr.columns, TEST_COLS))

    def test_interleave(self):
        pass

    def test_consolidate(self):
        pass

    def test_consolidate_ordering_issues(self):
        self.mgr.set('f', randn(N))
        self.mgr.set('d', randn(N))
        self.mgr.set('b', randn(N))

        cons = self.mgr.consolidate()
        self.assertEquals(cons.nblocks, 1)
        self.assert_(cons.blocks[0].columns.equals(cons.columns))

    def test_reindex_index(self):
        pass

    def test_reindex_columns(self):
        def _check_cols(before, after, cols):
            for col in cols:
                assert_almost_equal(after.get(col), before.get(col))

        # not consolidated
        vals = randn(N)
        self.mgr.set('g', vals)
        reindexed = self.mgr.reindex_columns(['g', 'c', 'a', 'd'])
        self.assertEquals(reindexed.nblocks, 2)
        assert_almost_equal(reindexed.get('g'), vals)
        _check_cols(self.mgr, reindexed, ['c', 'a', 'd'])

    def test_xs(self):
        pass


if __name__ == '__main__':
    # unittest.main()
    import nose
    # nose.runmodule(argv=[__file__,'-vvs','-x', '--pdb-failure'],
    #                exit=False)
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

