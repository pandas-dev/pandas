# pylint: disable=W0102

import nose
import numpy as np

from pandas import Index, MultiIndex, DataFrame, Series
from pandas.sparse.array import SparseArray
from pandas.core.internals import *
import pandas.core.internals as internals
import pandas.util.testing as tm

from pandas.util.testing import (
    assert_almost_equal, assert_frame_equal, randn)
from pandas.compat import zip, u


def assert_block_equal(left, right):
    assert_almost_equal(left.values, right.values)
    assert(left.dtype == right.dtype)
    assert(left.items.equals(right.items))
    assert(left.ref_items.equals(right.ref_items))


def get_float_mat(n, k, dtype):
    return np.repeat(np.atleast_2d(np.arange(k, dtype=dtype)), n, axis=0)

TEST_COLS = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 's1', 's2']
N = 10


def get_float_ex(cols=['a', 'c', 'e'], dtype = np.float_):
    floats = get_float_mat(N, len(cols), dtype = dtype).T
    return make_block(floats, cols, TEST_COLS)


def get_complex_ex(cols=['h']):
    complexes = (get_float_mat(N, 1, dtype = np.float_).T * 1j).astype(np.complex128)
    return make_block(complexes, cols, TEST_COLS)


def get_obj_ex(cols=['b', 'd']):
    mat = np.empty((N, 2), dtype=object)
    mat[:, 0] = 'foo'
    mat[:, 1] = 'bar'
    return make_block(mat.T, cols, TEST_COLS)

def get_bool_ex(cols=['f']):
    mat = np.ones((N, 1), dtype=bool)
    return make_block(mat.T, cols, TEST_COLS)


def get_int_ex(cols=['g'], dtype = np.int_):
    mat = randn(N, 1).astype(dtype)
    return make_block(mat.T, cols, TEST_COLS)


def get_dt_ex(cols=['h']):
    mat = randn(N, 1).astype(int).astype('M8[ns]')
    return make_block(mat.T, cols, TEST_COLS)

def get_sparse_ex1():
    sa1 = SparseArray([0, 0, 1, 2, 3, 0, 4, 5, 0, 6], fill_value=0)
    return make_block(sa1, ['s1'], TEST_COLS)

def get_sparse_ex2():
    sa2 = SparseArray([0, 0, 2, 3, 4, 0, 6, 7, 0, 8], fill_value=0)
    return make_block(sa2, ['s2'], TEST_COLS)

def create_blockmanager(blocks):
    l = []
    for b in blocks:
        l.extend(b.items)
    items = Index(l)
    for b in blocks:
        b.ref_items = items

    index_sz = blocks[0].shape[1]
    return BlockManager(blocks, [items, np.arange(index_sz)])

def create_singleblockmanager(blocks):
    l = []
    for b in blocks:
        l.extend(b.items)
    items = Index(l)
    for b in blocks:
        b.ref_items = items

    return SingleBlockManager(blocks, [items])

class TestBlock(tm.TestCase):

    _multiprocess_can_split_ = True

    def setUp(self):
        self.fblock = get_float_ex()
        self.cblock = get_complex_ex()
        self.oblock = get_obj_ex()
        self.bool_block = get_bool_ex()
        self.int_block = get_int_ex()

    def test_constructor(self):
        int32block = get_int_ex(['a'],dtype = np.int32)
        self.assert_(int32block.dtype == np.int32)

    def test_pickle(self):
        import pickle

        def _check(blk):
            pickled = pickle.dumps(blk)
            unpickled = pickle.loads(pickled)
            assert_block_equal(blk, unpickled)

        _check(self.fblock)
        _check(self.cblock)
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

        # with dup column support this method was taken out
        # GH3679
        raise nose.SkipTest("skipping for now")

        bs = list(self.fblock.split_block_at('a'))
        self.assertEqual(len(bs), 1)
        self.assertTrue(np.array_equal(bs[0].items, ['c', 'e']))

        bs = list(self.fblock.split_block_at('c'))
        self.assertEqual(len(bs), 2)
        self.assertTrue(np.array_equal(bs[0].items, ['a']))
        self.assertTrue(np.array_equal(bs[1].items, ['e']))

        bs = list(self.fblock.split_block_at('e'))
        self.assertEqual(len(bs), 1)
        self.assertTrue(np.array_equal(bs[0].items, ['a', 'c']))

        bblock = get_bool_ex(['f'])
        bs = list(bblock.split_block_at('f'))
        self.assertEqual(len(bs), 0)

    def test_unicode_repr(self):
        mat = np.empty((N, 2), dtype=object)
        mat[:, 0] = 'foo'
        mat[:, 1] = 'bar'
        cols = ['b', u("\u05d0")]
        str_repr = repr(make_block(mat.T, cols, TEST_COLS))

    def test_get(self):
        pass

    def test_set(self):
        pass

    def test_fillna(self):
        pass

    def test_repr(self):
        pass


class TestBlockManager(tm.TestCase):

    _multiprocess_can_split_ = True

    def setUp(self):
        self.blocks = [get_float_ex(),
                       get_obj_ex(),
                       get_bool_ex(),
                       get_int_ex(),
                       get_complex_ex()]

        all_items = [b.items for b in self.blocks]

        items = sorted(all_items[0].append(all_items[1:]))
        items = Index(items)
        for b in self.blocks:
            b.ref_items = items

        self.mgr = BlockManager(self.blocks, [items, np.arange(N)])

    def test_constructor_corner(self):
        pass

    def test_attrs(self):
        self.assertEquals(self.mgr.nblocks, len(self.mgr.blocks))
        self.assertEquals(len(self.mgr), len(self.mgr.items))

    def test_is_mixed_dtype(self):
        self.assert_(self.mgr.is_mixed_type)

        mgr = create_blockmanager([get_bool_ex(['a']), get_bool_ex(['b'])])
        self.assert_(not mgr.is_mixed_type)

    def test_is_indexed_like(self):
        self.assert_(self.mgr._is_indexed_like(self.mgr))
        mgr2 = self.mgr.reindex_axis(np.arange(N - 1), axis=1)
        self.assert_(not self.mgr._is_indexed_like(mgr2))

    def test_block_id_vector_item_dtypes(self):
        expected = [0, 1, 0, 1, 0, 2, 3, 4]
        result = self.mgr.block_id_vector
        assert_almost_equal(expected, result)

        result = self.mgr.item_dtypes

        # as the platform may not exactly match this, pseudo match
        expected = ['float64', 'object', 'float64', 'object', 'float64',
                    'bool', 'int64', 'complex128']
        for e, r in zip(expected, result):
            np.dtype(e).kind == np.dtype(r).kind

    def test_duplicate_item_failure(self):
        items = Index(['a', 'a'])
        blocks = [get_bool_ex(['a']), get_float_ex(['a'])]
        for b in blocks:
            b.ref_items = items

        # test trying to create _ref_locs with/o ref_locs set on the blocks
        self.assertRaises(AssertionError, BlockManager, blocks, [items, np.arange(N)])

        blocks[0].set_ref_locs([0])
        blocks[1].set_ref_locs([1])
        mgr = BlockManager(blocks, [items, np.arange(N)])
        mgr.iget(1)

        # invalidate the _ref_locs
        for b in blocks:
            b._ref_locs = None
        mgr._ref_locs = None
        mgr._items_map = None
        self.assertRaises(AssertionError, mgr._set_ref_locs, do_refs=True)

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

        # GH2431
        self.assertTrue(hasattr(mgr2, "_is_consolidated"))
        self.assertTrue(hasattr(mgr2, "_known_consolidated"))

        # reset to False on load
        self.assertFalse(mgr2._is_consolidated)
        self.assertFalse(mgr2._known_consolidated)

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
        self.assert_(mgr2.get('quux').dtype == np.int_)

        mgr2.set('quux', randn(N))
        self.assert_(mgr2.get('quux').dtype == np.float_)

    def test_copy(self):
        shallow = self.mgr.copy(deep=False)

        # we don't guaranteee block ordering
        for blk in self.mgr.blocks:
            found = False
            for cp_blk in shallow.blocks:
                if cp_blk.values is blk.values:
                    found = True
                    break
            self.assert_(found == True)

    def test_sparse(self):
        mgr = create_blockmanager([get_sparse_ex1(),get_sparse_ex2()])

        # what to test here?
        self.assert_(mgr.as_matrix().dtype == np.float64)

    def test_sparse_mixed(self):
        mgr = create_blockmanager([get_sparse_ex1(),get_sparse_ex2(),get_float_ex()])
        self.assert_(len(mgr.blocks) == 3)
        self.assert_(isinstance(mgr,BlockManager))

        # what to test here?

    def test_as_matrix_float(self):

        mgr = create_blockmanager([get_float_ex(['c'],np.float32), get_float_ex(['d'],np.float16), get_float_ex(['e'],np.float64)])
        self.assert_(mgr.as_matrix().dtype == np.float64)

        mgr = create_blockmanager([get_float_ex(['c'],np.float32), get_float_ex(['d'],np.float16)])
        self.assert_(mgr.as_matrix().dtype == np.float32)

    def test_as_matrix_int_bool(self):

        mgr = create_blockmanager([get_bool_ex(['a']), get_bool_ex(['b'])])
        self.assert_(mgr.as_matrix().dtype == np.bool_)

        mgr = create_blockmanager([get_int_ex(['a'],np.int64), get_int_ex(['b'],np.int64), get_int_ex(['c'],np.int32), get_int_ex(['d'],np.int16), get_int_ex(['e'],np.uint8) ])
        self.assert_(mgr.as_matrix().dtype == np.int64)

        mgr = create_blockmanager([get_int_ex(['c'],np.int32), get_int_ex(['d'],np.int16), get_int_ex(['e'],np.uint8) ])
        self.assert_(mgr.as_matrix().dtype == np.int32)

    def test_as_matrix_datetime(self):
        mgr = create_blockmanager([get_dt_ex(['h']), get_dt_ex(['g'])])
        self.assert_(mgr.as_matrix().dtype == 'M8[ns]')

    def test_astype(self):

        # coerce all
        mgr = create_blockmanager([get_float_ex(['c'],np.float32), get_float_ex(['d'],np.float16), get_float_ex(['e'],np.float64)])

        for t in ['float16','float32','float64','int32','int64']:
            tmgr = mgr.astype(t)
            self.assert_(tmgr.as_matrix().dtype == np.dtype(t))

        # mixed
        mgr = create_blockmanager([get_obj_ex(['a','b']),get_bool_ex(['c']),get_dt_ex(['d']),get_float_ex(['e'],np.float32), get_float_ex(['f'],np.float16), get_float_ex(['g'],np.float64)])
        for t in ['float16','float32','float64','int32','int64']:
            tmgr = mgr.astype(t, raise_on_error = False).get_numeric_data()
            self.assert_(tmgr.as_matrix().dtype == np.dtype(t))

    def test_convert(self):

        def _compare(old_mgr, new_mgr):
            """ compare the blocks, numeric compare ==, object don't """
            old_blocks = set(old_mgr.blocks)
            new_blocks = set(new_mgr.blocks)
            self.assert_(len(old_blocks) == len(new_blocks))

            # compare non-numeric
            for b in old_blocks:
                found = False
                for nb in new_blocks:
                    if (b.values == nb.values).all():
                        found = True
                        break
                self.assert_(found == True)

            for b in new_blocks:
                found = False
                for ob in old_blocks:
                    if (b.values == ob.values).all():
                        found = True
                        break
                self.assert_(found == True)

        # noops
        mgr = create_blockmanager([get_int_ex(['f']), get_float_ex(['g'])])
        new_mgr = mgr.convert()
        _compare(mgr,new_mgr)

        mgr = create_blockmanager([get_obj_ex(['a','b']), get_int_ex(['f']), get_float_ex(['g'])])
        new_mgr = mgr.convert()
        _compare(mgr,new_mgr)

        # there could atcually be multiple dtypes resulting
        def _check(new_mgr,block_type, citems):
            items = set()
            for b in new_mgr.blocks:
                if isinstance(b,block_type):
                    for i in list(b.items):
                        items.add(i)
            self.assert_(items == set(citems))

        # convert
        mat = np.empty((N, 3), dtype=object)
        mat[:, 0] = '1'
        mat[:, 1] = '2.'
        mat[:, 2] = 'foo'
        b = make_block(mat.T, ['a','b','foo'], TEST_COLS)

        mgr = create_blockmanager([b, get_int_ex(['f']), get_float_ex(['g'])])
        new_mgr = mgr.convert(convert_numeric = True)

        _check(new_mgr,FloatBlock,['b','g'])
        _check(new_mgr,IntBlock,['a','f'])

        mgr = create_blockmanager([b, get_int_ex(['f'],np.int32), get_bool_ex(['bool']), get_dt_ex(['dt']),
                                   get_int_ex(['i'],np.int64), get_float_ex(['g'],np.float64), get_float_ex(['h'],np.float16)])
        new_mgr = mgr.convert(convert_numeric = True)

        _check(new_mgr,FloatBlock,['b','g','h'])
        _check(new_mgr,IntBlock,['a','f','i'])
        _check(new_mgr,ObjectBlock,['foo'])
        _check(new_mgr,BoolBlock,['bool'])
        _check(new_mgr,DatetimeBlock,['dt'])

    def test_interleave(self):
        pass

    def test_interleave_non_unique_cols(self):
        df = DataFrame([
            [Timestamp('20130101'), 3.5],
            [Timestamp('20130102'), 4.5]],
            columns=['x', 'x'],
            index=[1, 2])

        df_unique = df.copy()
        df_unique.columns = ['x', 'y']
        np.testing.assert_array_equal(df_unique.values, df.values)

    def test_consolidate(self):
        pass

    def test_consolidate_ordering_issues(self):
        self.mgr.set('f', randn(N))
        self.mgr.set('d', randn(N))
        self.mgr.set('b', randn(N))
        self.mgr.set('g', randn(N))
        self.mgr.set('h', randn(N))

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

    def test_get_numeric_data(self):
        int_ser = Series(np.array([0, 1, 2]))
        float_ser = Series(np.array([0., 1., 2.]))
        complex_ser = Series(np.array([0j, 1j, 2j]))
        str_ser = Series(np.array(['a', 'b', 'c']))
        bool_ser = Series(np.array([True, False, True]))
        obj_ser = Series(np.array([1, 'a', 5]))
        dt_ser = Series(tm.makeDateIndex(3))
        # check types
        df = DataFrame({'int': int_ser, 'float': float_ser,
                        'complex': complex_ser, 'str': str_ser,
                        'bool': bool_ser, 'obj': obj_ser,
                        'dt': dt_ser})
        xp = DataFrame({'int': int_ser, 'float': float_ser,
                        'complex': complex_ser, 'bool': bool_ser})
        rs = DataFrame(df._data.get_numeric_data())
        assert_frame_equal(xp, rs)

        xp = DataFrame({'bool': bool_ser})
        rs = DataFrame(df._data.get_bool_data())
        assert_frame_equal(xp, rs)

        rs = DataFrame(df._data.get_bool_data())
        df.ix[0, 'bool'] = not df.ix[0, 'bool']

        self.assertEqual(rs.ix[0, 'bool'], df.ix[0, 'bool'])

        rs = DataFrame(df._data.get_bool_data(copy=True))
        df.ix[0, 'bool'] = not df.ix[0, 'bool']

        self.assertEqual(rs.ix[0, 'bool'], not df.ix[0, 'bool'])

    def test_missing_unicode_key(self):
        df = DataFrame({"a": [1]})
        try:
            df.ix[:, u("\u05d0")]  # should not raise UnicodeEncodeError
        except KeyError:
            pass  # this is the expected exception

    def test_equals(self):
        # unique items
        index = Index(list('abcdef'))
        block1 = make_block(np.arange(12).reshape(3,4), list('abc'), index)
        block2 = make_block(np.arange(12).reshape(3,4)*10, list('def'), index)
        block1.ref_items = block2.ref_items = index
        bm1 = BlockManager([block1, block2], [index, np.arange(block1.shape[1])])
        bm2 = BlockManager([block2, block1], [index, np.arange(block1.shape[1])])
        self.assert_(bm1.equals(bm2))

        # non-unique items
        index = Index(list('aaabbb'))
        block1 = make_block(np.arange(12).reshape(3,4), list('aaa'), index,
                            placement=[0,1,2])
        block2 = make_block(np.arange(12).reshape(3,4)*10, list('bbb'), index,
                            placement=[3,4,5])
        block1.ref_items = block2.ref_items = index
        bm1 = BlockManager([block1, block2], [index, np.arange(block1.shape[1])])
        bm2 = BlockManager([block2, block1], [index, np.arange(block1.shape[1])])
        self.assert_(bm1.equals(bm2))
        
if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
