# pylint: disable=W0102

import unittest

import numpy as np

from pandas import Index, MultiIndex, DataFrame
from pandas.core.internals import *
import pandas.core.internals as internals

from pandas.util.testing import (assert_almost_equal, assert_frame_equal, randn)

def assert_block_equal(left, right):
    assert_almost_equal(left.values, right.values)
    assert(left.dtype == right.dtype)
    assert(left.items.equals(right.items))
    assert(left.ref_items.equals(right.ref_items))

def get_float_mat(n, k):
    return np.repeat(np.atleast_2d(np.arange(k, dtype=float)), n, axis=0)

TEST_COLS = ['a', 'b', 'c', 'd', 'e', 'f', 'g']
N = 10

def get_float_ex(cols=['a', 'c', 'e']):
    floats = get_float_mat(N, 3).T
    return make_block(floats, cols, TEST_COLS)

def get_obj_ex(cols=['b', 'd']):
    mat = np.empty((N, 2), dtype=object)
    mat[:, 0] = 'foo'
    mat[:, 1] = 'bar'
    return make_block(mat.T, cols, TEST_COLS)

def get_bool_ex(cols=['f']):
    mat = np.ones((N, 1), dtype=bool)
    return make_block(mat.T, cols, TEST_COLS)

def get_int_ex(cols=['g']):
    mat = randn(N, 1).astype(int)
    return make_block(mat.T, cols, TEST_COLS)

class TestBlock(unittest.TestCase):

    def setUp(self):
        self.fblock = get_float_ex()
        self.oblock = get_obj_ex()
        self.bool_block = get_bool_ex()
        self.int_block = get_int_ex()

    def test_constructor(self):
        pass

    def test_pickle(self):
        import pickle

        def _check(blk):
            pickled = pickle.dumps(blk)
            unpickled = pickle.loads(pickled)
            assert_block_equal(blk, unpickled)

        _check(self.fblock)
        _check(self.oblock)
        _check(self.bool_block)

    def test_ref_locs(self):
        assert_almost_equal(self.fblock.ref_locs, [0, 2, 4])

    def test_attrs(self):
        self.assert_(self.fblock.shape == self.fblock.values.shape)
        self.assert_(self.fblock.dtype == self.fblock.values.dtype)
        self.assert_(len(self.fblock) == len(self.fblock.values))

    def test_merge(self):
        avals = randn(2, 10)
        bvals = randn(2, 10)

        ref_cols = ['e', 'a', 'b', 'd', 'f']

        ablock = make_block(avals, ['e', 'b'], ref_cols)
        bblock = make_block(bvals, ['a', 'd'], ref_cols)
        merged = ablock.merge(bblock)
        exvals = np.vstack((avals, bvals))
        excols = ['e', 'b', 'a', 'd']
        eblock = make_block(exvals, excols, ref_cols)
        eblock = eblock.reindex_items_from(ref_cols)
        assert_block_equal(merged, eblock)

        # TODO: merge with mixed type?

    def test_copy(self):
        cop = self.fblock.copy()
        self.assert_(cop is not self.fblock)
        assert_block_equal(self.fblock, cop)

    def test_items(self):
        cols = self.fblock.items
        self.assert_(np.array_equal(cols, ['a', 'c', 'e']))

        cols2 = self.fblock.items
        self.assert_(cols is cols2)

    def test_assign_ref_items(self):
        new_cols = Index(['foo', 'bar', 'baz', 'quux', 'hi'])
        self.fblock.set_ref_items(new_cols)
        self.assert_(np.array_equal(self.fblock.items,
                                    ['foo', 'baz', 'hi']))

    def test_reindex_index(self):
        pass

    def test_reindex_items_from(self):
        new_cols = Index(['e', 'b', 'c', 'f'])
        reindexed = self.fblock.reindex_items_from(new_cols)
        assert_almost_equal(reindexed.ref_locs, [0, 2])
        self.assertEquals(reindexed.values.shape[0], 2)
        self.assert_((reindexed.values[0] == 2).all())
        self.assert_((reindexed.values[1] == 1).all())

    def test_reindex_cast(self):
        pass

    def test_insert(self):
        pass

    def test_delete(self):
        newb = self.fblock.delete('a')
        assert_almost_equal(newb.ref_locs, [2, 4])
        self.assert_((newb.values[0] == 1).all())

        newb = self.fblock.delete('c')
        assert_almost_equal(newb.ref_locs, [0, 4])
        self.assert_((newb.values[1] == 2).all())

        newb = self.fblock.delete('e')
        assert_almost_equal(newb.ref_locs, [0, 2])
        self.assert_((newb.values[1] == 1).all())

        self.assertRaises(Exception, self.fblock.delete, 'b')

    def test_split_block_at(self):
        left, right = self.fblock.split_block_at('a')
        self.assert_(left is None)
        self.assert_(np.array_equal(right.items, ['c', 'e']))

        left, right = self.fblock.split_block_at('c')
        self.assert_(np.array_equal(left.items, ['a']))
        self.assert_(np.array_equal(right.items, ['e']))

        left, right = self.fblock.split_block_at('e')
        self.assert_(np.array_equal(left.items, ['a', 'c']))
        self.assert_(right is None)

        bblock = get_bool_ex(['f'])
        left, right = bblock.split_block_at('f')
        self.assert_(left is None)
        self.assert_(right is None)

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
                       get_bool_ex(),
                       get_int_ex()]
        self.mgr = BlockManager.from_blocks(self.blocks, np.arange(N))

    def test_constructor_corner(self):
        pass

    def test_attrs(self):
        self.assertEquals(self.mgr.nblocks, len(self.mgr.blocks))
        self.assertEquals(len(self.mgr), len(self.mgr.items))

    def test_is_mixed_dtype(self):
        self.assert_(self.mgr.is_mixed_dtype())

        blocks = [get_bool_ex(['a']), get_bool_ex(['b'])]
        mgr = BlockManager.from_blocks(blocks, np.arange(N))
        self.assert_(not mgr.is_mixed_dtype())

    def test_is_indexed_like(self):
        self.assert_(self.mgr._is_indexed_like(self.mgr))
        mgr2 = self.mgr.reindex_axis(np.arange(N - 1), axis=1)
        self.assert_(not self.mgr._is_indexed_like(mgr2))

    def test_block_id_vector_item_dtypes(self):
        expected = [0, 1, 0, 1, 0, 2, 3]
        result = self.mgr.block_id_vector
        assert_almost_equal(expected, result)

        result = self.mgr.item_dtypes
        expected = ['float64', 'object', 'float64', 'object', 'float64',
                    'bool', 'int64']
        self.assert_(np.array_equal(result, expected))

    def test_union_block_items(self):
        blocks = [get_float_ex(['a', 'b', 'c']),
                  get_float_ex(['c', 'd', 'e'])]
        self.assertRaises(Exception, internals._union_block_items, blocks)

        blocks = [get_float_ex(['a', 'b', 'c']),
                  get_float_ex(['f', 'e', 'd'])]
        self.assert_(np.array_equal(internals._union_block_items(blocks),
                                    ['a', 'b', 'c', 'd', 'e', 'f']))

    def test_contains(self):
        self.assert_('a' in self.mgr)
        self.assert_('baz' not in self.mgr)

    def test_pickle(self):
        import pickle

        pickled = pickle.dumps(self.mgr)
        mgr2 = pickle.loads(pickled)

        # same result
        assert_frame_equal(DataFrame(self.mgr), DataFrame(mgr2))

        # share ref_items
        self.assert_(mgr2.blocks[0].ref_items is mgr2.blocks[1].ref_items)

    def test_get(self):
        pass

    def test_get_scalar(self):
        for item in self.mgr.items:
            for i, index in enumerate(self.mgr.axes[1]):
                res = self.mgr.get_scalar((item, index))
                exp = self.mgr.get(item)[i]
                assert_almost_equal(res, exp)

    def test_set(self):
        pass

    def test_set_change_dtype(self):
        self.mgr.set('baz', np.zeros(N, dtype=bool))

        self.mgr.set('baz', np.repeat('foo', N))
        self.assert_(self.mgr.get('baz').dtype == np.object_)

        mgr2 = self.mgr.consolidate()
        mgr2.set('baz', np.repeat('foo', N))
        self.assert_(mgr2.get('baz').dtype == np.object_)

        mgr2.set('quux', randn(N).astype(int))
        self.assert_(mgr2.get('quux').dtype == np.int64)

        mgr2.set('quux', randn(N))
        self.assert_(mgr2.get('quux').dtype == np.float_)

    def test_copy(self):
        shallow = self.mgr.copy(deep=False)

        for cp_blk, blk in zip(shallow.blocks, self.mgr.blocks):
            self.assert_(cp_blk.values is blk.values)

    def test_as_matrix(self):
        pass

    def test_as_matrix_int_bool(self):
        blocks = [get_bool_ex(['a']), get_bool_ex(['b'])]
        index_sz = blocks[0].values.shape[1]
        mgr = BlockManager.from_blocks(blocks, np.arange(index_sz))
        self.assert_(mgr.as_matrix().dtype == np.bool_)

        blocks = [get_int_ex(['a']), get_int_ex(['b'])]
        mgr = BlockManager.from_blocks(blocks, np.arange(index_sz))
        self.assert_(mgr.as_matrix().dtype == np.int64)

    def test_xs(self):
        pass

    def test_from_blocks(self):
        self.assert_(np.array_equal(self.mgr.items, TEST_COLS))

    def test_interleave(self):
        pass

    def test_consolidate(self):
        pass

    def test_consolidate_ordering_issues(self):
        self.mgr.set('f', randn(N))
        self.mgr.set('d', randn(N))
        self.mgr.set('b', randn(N))
        self.mgr.set('g', randn(N))

        cons = self.mgr.consolidate()
        self.assertEquals(cons.nblocks, 1)
        self.assert_(cons.blocks[0].items.equals(cons.items))

    def test_reindex_index(self):
        pass

    def test_reindex_items(self):
        def _check_cols(before, after, cols):
            for col in cols:
                assert_almost_equal(after.get(col), before.get(col))

        # not consolidated
        vals = randn(N)
        self.mgr.set('g', vals)
        reindexed = self.mgr.reindex_items(['g', 'c', 'a', 'd'])
        self.assertEquals(reindexed.nblocks, 2)
        assert_almost_equal(reindexed.get('g'), vals.squeeze())
        _check_cols(self.mgr, reindexed, ['c', 'a', 'd'])

    def test_xs(self):
        index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'],
                                   ['one', 'two', 'three']],
                           labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                                   [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                           names=['first', 'second'])

        self.mgr.set_axis(1, index)

        result = self.mgr.xs('bar', axis=1)
        expected = self.mgr.get_slice(slice(3, 5), axis=1)

        assert_frame_equal(DataFrame(result), DataFrame(expected))

if __name__ == '__main__':
    # unittest.main()
    import nose
    # nose.runmodule(argv=[__file__,'-vvs','-x', '--pdb-failure'],
    #                exit=False)
    nose.runmodule(argv=[__file__,'-vvs','-x','--pdb', '--pdb-failure'],
                   exit=False)

