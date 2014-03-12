# pylint: disable=W0102

import nose
import numpy as np

from pandas import Index, MultiIndex, DataFrame, Series
from pandas.compat import OrderedDict
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
    assert_almost_equal(left.ref_locs, right.ref_locs)


def get_numeric_mat(n, k, dtype):
    return np.repeat(np.atleast_2d(np.arange(k, dtype=dtype)), n, axis=0)


N = 10


def create_block(typestr, placement, num_rows=None, num_offset=None):
    placement = np.asanyarray(placement)

    if num_offset is None:
        num_offset = 0

    if num_rows is None:
        num_rows = N

    if typestr in ('float', 'f8', 'f4', 'f2',
                   'int', 'i8', 'i4', 'i2', 'i1',
                   'uint', 'u8', 'u4', 'u2', 'u1'):
        values = get_numeric_mat(num_rows, len(placement),
                                 dtype=np.dtype(typestr)).T + num_offset
    elif typestr in ('complex', 'c16', 'c8'):
        values = get_numeric_mat(num_rows, len(placement),
                                 dtype=np.dtype(typestr)).T + num_offset
        values *= 1.j
    elif typestr in ('object', 'string', 'O'):
        values = np.repeat(
            np.array([['A%s' % i
                       for i in np.arange(len(placement)) + num_offset]]),
            num_rows, axis=0).T
    elif typestr in ('bool'):
        values = np.ones((num_rows, len(placement)), dtype=np.bool_).T
    elif typestr in ('datetime', 'dt'):
        values = (randn(num_rows, len(placement)).astype(int)
                  .astype('M8[ns]')).T
    elif typestr in ('sparse',):
        # FIXME: doesn't support num_rows != 10
        assert len(placement) == 1
        assert num_rows == 10
        values = SparseArray([0, 0, 1, 2, 3, 0, 4, 5, 0, 6], fill_value=0)
        arr = values.sp_values.view()
        arr += (num_offset - 1)
    else:
        raise ValueError('Unsupported typestr: "%s"' % typestr)

    return make_block(values, placement=placement)


def create_mgr(descr, num_rows=None):
    """
    Construct BlockManager from string description.

    String description syntax looks similar to np.matrix initializer.  It looks
    like this::

        a,b,c: f8; d,e,f: i8

    Rules are rather simple:

    * supported datatypes:

      * float, f8, f4, f2
      * int, i8, i4, i2, i1
      * uint, u8, u4, u2, u1
      * complex, c16, c8
      * bool
      * object, string, O
      * datetime, dt
      * sparse

    * components are semicolon-separated
    * each component is `NAME,NAME,NAME: DTYPE_ID`
    * whitespace around colons & semicolons are removed
    * components with same DTYPE_ID are combined into single block
    * to force multiple blocks with same dtype, use '-SUFFIX'::

        'a:f8-1; b:f8-2; c:f8-foobar'

    """
    if num_rows is None:
        num_rows = N

    offset = 0
    mgr_items = []
    block_placements = OrderedDict()
    for d in descr.split(';'):
        d = d.strip()
        names, blockstr = d.partition(':')[::2]
        blockstr = blockstr.strip()
        names = names.strip().split(',')

        mgr_items.extend(names)
        placement = list(np.arange(len(names)) + offset)
        try:
            block_placements[blockstr].extend(placement)
        except KeyError:
            block_placements[blockstr] = placement
        offset += len(names)

    mgr_items = Index(mgr_items)

    blocks = []
    num_offset = 0
    for blockstr, placement in block_placements.items():
        typestr = blockstr.split('-')[0]
        blocks.append(create_block(typestr, placement, num_rows=num_rows,
                                   num_offset=num_offset,))
        num_offset += len(placement)

    return BlockManager(sorted(blocks, key=lambda b: b.ref_locs[0]),
                        [mgr_items, np.arange(num_rows)])



class TestBlock(tm.TestCase):

    _multiprocess_can_split_ = True

    def setUp(self):
        # self.fblock = get_float_ex()  # a,c,e
        # self.cblock = get_complex_ex() #
        # self.oblock = get_obj_ex()
        # self.bool_block = get_bool_ex()
        # self.int_block = get_int_ex()

        self.fblock = create_block('float', [0, 2, 4])
        self.cblock = create_block('complex', [7])
        self.oblock = create_block('object', [1, 3])
        self.bool_block = create_block('bool', [5])
        self.int_block = create_block('int', [6])

    def test_constructor(self):
        int32block = create_block('i4', [0])
        self.assertEqual(int32block.dtype, np.int32)

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
        self.assertEqual(self.fblock.shape, self.fblock.values.shape)
        self.assertEqual(self.fblock.dtype, self.fblock.values.dtype)
        self.assertEqual(len(self.fblock), len(self.fblock.values))

    def test_merge(self):
        avals = randn(2, 10)
        bvals = randn(2, 10)

        ref_cols = Index(['e', 'a', 'b', 'd', 'f'])

        ablock = make_block(avals,
                            ref_cols.get_indexer(['e', 'b']))
        bblock = make_block(bvals,
                            ref_cols.get_indexer(['a', 'd']))
        merged = ablock.merge(bblock)
        assert_almost_equal(merged.ref_locs, [0, 1, 2, 3])
        assert_almost_equal(merged.values[[0, 2]], avals)
        assert_almost_equal(merged.values[[1, 3]], bvals)

        # TODO: merge with mixed type?

    def test_copy(self):
        cop = self.fblock.copy()
        self.assertIsNot(cop, self.fblock)
        assert_block_equal(self.fblock, cop)

    def test_items(self):
        raise nose.SkipTest('items are removed from Block')
        cols = self.fblock.items
        self.assert_numpy_array_equal(cols, ['a', 'c', 'e'])

        cols2 = self.fblock.items
        # disabled: items are generated
        # self.assertIs(cols, cols2)

    def test_assign_ref_items(self):
        raise nose.SkipTest('ref_items are removed from Block')
        new_cols = Index(['foo', 'bar', 'baz', 'quux', 'hi'])
        self.fblock.set_ref_items(new_cols)
        self.assert_numpy_array_equal(self.fblock.items, ['foo', 'baz', 'hi'])

    def test_reindex_index(self):
        pass

    def test_reindex_items_from(self):
        raise nose.SkipTest('reindex_items_from is removed from Block')
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
        newb = self.fblock.delete(0)
        assert_almost_equal(newb.ref_locs, [2, 4])
        self.assert_((newb.values[0] == 1).all())

        newb = self.fblock.delete(1)
        assert_almost_equal(newb.ref_locs, [0, 4])
        self.assert_((newb.values[1] == 2).all())

        newb = self.fblock.delete(2)
        assert_almost_equal(newb.ref_locs, [0, 2])
        self.assert_((newb.values[1] == 1).all())

        self.assertRaises(Exception, self.fblock.delete, 3)

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
        raise nose.SkipTest('No items to test unicode on...')
        str_repr = repr(create_block('object', [0, 1],
                                     ref_items=['b', u("\u05d0")]))

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
        self.mgr = create_mgr('a: f8; b: object; c: f8; d: object; e: f8;'
                              'f: bool; g: i8; h: complex')

    def test_constructor_corner(self):
        pass

    def test_attrs(self):
        mgr = create_mgr('a,b,c: f8-1; d,e,f: f8-2')
        self.assertEquals(mgr.nblocks, 2)
        self.assertEquals(len(mgr), 6)

    def test_is_mixed_dtype(self):
        self.assertFalse(create_mgr('a,b:f8').is_mixed_type)
        self.assertFalse(create_mgr('a:f8-1; b:f8-2').is_mixed_type)

        self.assertTrue(create_mgr('a,b:f8; c,d: f4').is_mixed_type)
        self.assertTrue(create_mgr('a,b:f8; c,d: object').is_mixed_type)

    def test_is_indexed_like(self):
        mgr1 = create_mgr('a,b: f8')
        mgr2 = create_mgr('a:i8; b:bool')
        mgr3 = create_mgr('a,b,c: f8')
        self.assertTrue(mgr1._is_indexed_like(mgr1))
        self.assertTrue(mgr1._is_indexed_like(mgr2))
        self.assertTrue(mgr1._is_indexed_like(mgr3))

        self.assertFalse(mgr1._is_indexed_like(
            mgr1.get_slice(slice(-1), axis=1)))

    def test_duplicate_ref_loc_failure(self):
        tmp_mgr = create_mgr('a:bool; a: f8')

        axes, blocks = tmp_mgr.axes, tmp_mgr.blocks

        blocks[0]._ref_locs = np.array([0])
        blocks[1]._ref_locs = np.array([0])
        # test trying to create block manager with overlapping ref locs
        self.assertRaises(AssertionError, BlockManager, blocks, axes)

        blocks[0]._ref_locs = np.array([0])
        blocks[1]._ref_locs = np.array([1])
        mgr = BlockManager(blocks, axes)
        mgr.iget(1)

        # invalidate the _ref_locs
        for b in blocks:
            b._ref_locs = None
        mgr._ref_locs = None
        mgr._items_map = None
        self.assertRaises(Exception, mgr._rebuild_ref_locs)

    def test_contains(self):
        self.assertIn('a', self.mgr)
        self.assertNotIn('baz', self.mgr)

    def test_pickle(self):
        import pickle

        pickled = pickle.dumps(self.mgr)
        mgr2 = pickle.loads(pickled)

        # same result
        assert_frame_equal(DataFrame(self.mgr), DataFrame(mgr2))

        # share ref_items
        # self.assertIs(mgr2.blocks[0].ref_items, mgr2.blocks[1].ref_items)

        # GH2431
        self.assertTrue(hasattr(mgr2, "_is_consolidated"))
        self.assertTrue(hasattr(mgr2, "_known_consolidated"))

        # reset to False on load
        self.assertFalse(mgr2._is_consolidated)
        self.assertFalse(mgr2._known_consolidated)

    def test_get_scalar(self):
        for item in self.mgr.items:
            for i, index in enumerate(self.mgr.axes[1]):
                res = self.mgr.get_scalar((item, index))
                exp = self.mgr.get(item)[i]
                assert_almost_equal(res, exp)

    def test_get(self):
        cols = Index(list('abc'))
        values = np.random.rand(3, 3)
        block = make_block(values=values.copy(),
                           placement=np.arange(3))
        mgr = BlockManager(blocks=[block], axes=[cols, np.arange(3)])

        assert_almost_equal(mgr.get('a'), values[0])
        assert_almost_equal(mgr.get('b'), values[1])
        assert_almost_equal(mgr.get('c'), values[2])

    def test_set(self):
        mgr = create_mgr('a,b,c: int', num_rows=3)

        mgr.set('d', np.array(['foo'] * 3))
        mgr.set('b', np.array(['bar'] * 3))
        assert_almost_equal(mgr.get('a'), [0] * 3)
        assert_almost_equal(mgr.get('b'), ['bar'] * 3)
        assert_almost_equal(mgr.get('c'), [2] * 3)
        assert_almost_equal(mgr.get('d'), ['foo'] * 3)

    def test_insert(self):
        self.mgr.insert(0, 'inserted', np.arange(N))

        self.assertEqual(self.mgr.items[0], 'inserted')
        assert_almost_equal(self.mgr.get('inserted'), np.arange(N))

        for blk in self.mgr.blocks:
            yield self.assertIs, self.mgr.items, blk.ref_items

    def test_set_change_dtype(self):
        self.mgr.set('baz', np.zeros(N, dtype=bool))

        self.mgr.set('baz', np.repeat('foo', N))
        self.assertEqual(self.mgr.get('baz').dtype, np.object_)

        mgr2 = self.mgr.consolidate()
        mgr2.set('baz', np.repeat('foo', N))
        self.assertEqual(mgr2.get('baz').dtype, np.object_)

        mgr2.set('quux', randn(N).astype(int))
        self.assertEqual(mgr2.get('quux').dtype, np.int_)

        mgr2.set('quux', randn(N))
        self.assertEqual(mgr2.get('quux').dtype, np.float_)

    def test_copy(self):
        shallow = self.mgr.copy(deep=False)

        # we don't guaranteee block ordering
        for blk in self.mgr.blocks:
            found = False
            for cp_blk in shallow.blocks:
                if cp_blk.values is blk.values:
                    found = True
                    break
            self.assertTrue(found)

    def test_sparse(self):
        mgr = create_mgr('a: sparse-1; b: sparse-2')

        # what to test here?
        self.assertEqual(mgr.as_matrix().dtype, np.float64)

    def test_sparse_mixed(self):
        mgr = create_mgr('a: sparse-1; b: sparse-2; c: f8')
        self.assertEqual(len(mgr.blocks), 3)
        self.assertIsInstance(mgr, BlockManager)

        # what to test here?

    def test_as_matrix_float(self):
        mgr = create_mgr('c: f4; d: f2; e: f8')
        self.assertEqual(mgr.as_matrix().dtype, np.float64)

        mgr = create_mgr('c: f4; d: f2')
        self.assertEqual(mgr.as_matrix().dtype, np.float32)

    def test_as_matrix_int_bool(self):
        mgr = create_mgr('a: bool-1; b: bool-2')
        self.assertEqual(mgr.as_matrix().dtype, np.bool_)

        mgr = create_mgr('a: i8-1; b: i8-2; c: i4; d: i2; e: u1')
        self.assertEqual(mgr.as_matrix().dtype, np.int64)

        mgr = create_mgr('c: i4; d: i2; e: u1')
        self.assertEqual(mgr.as_matrix().dtype, np.int32)

    def test_as_matrix_datetime(self):
        mgr = create_mgr('h: datetime-1; g: datetime-2')
        self.assertEqual(mgr.as_matrix().dtype, 'M8[ns]')

    def test_astype(self):
        # coerce all
        mgr = create_mgr('c: f4; d: f2; e: f8')
        for t in ['float16', 'float32', 'float64', 'int32', 'int64']:
            t = np.dtype(t)
            tmgr = mgr.astype(t)
            self.assertEqual(tmgr.get('c').dtype.type, t)
            self.assertEqual(tmgr.get('d').dtype.type, t)
            self.assertEqual(tmgr.get('e').dtype.type, t)

        # mixed
        mgr = create_mgr('a,b: object; c: bool; d: datetime;'
                         'e: f4; f: f2; g: f8')
        for t in ['float16', 'float32', 'float64', 'int32', 'int64']:
            t = np.dtype(t)
            tmgr = mgr.astype(t, raise_on_error=False)
            self.assertEqual(tmgr.get('c').dtype.type, t)
            self.assertEqual(tmgr.get('e').dtype.type, t)
            self.assertEqual(tmgr.get('f').dtype.type, t)
            self.assertEqual(tmgr.get('g').dtype.type, t)

            self.assertEqual(tmgr.get('a').dtype.type, np.object_)
            self.assertEqual(tmgr.get('b').dtype.type, np.object_)
            if t != np.int64:
                self.assertEqual(tmgr.get('d').dtype.type, np.datetime64)
            else:
                self.assertEqual(tmgr.get('d').dtype.type, t)

    def test_convert(self):
        def _compare(old_mgr, new_mgr):
            """ compare the blocks, numeric compare ==, object don't """
            old_blocks = set(old_mgr.blocks)
            new_blocks = set(new_mgr.blocks)
            self.assertEqual(len(old_blocks), len(new_blocks))

            # compare non-numeric
            for b in old_blocks:
                found = False
                for nb in new_blocks:
                    if (b.values == nb.values).all():
                        found = True
                        break
                self.assertTrue(found)

            for b in new_blocks:
                found = False
                for ob in old_blocks:
                    if (b.values == ob.values).all():
                        found = True
                        break
                self.assertTrue(found)

        # noops
        mgr = create_mgr('f: i8; g: f8')
        new_mgr = mgr.convert()
        _compare(mgr,new_mgr)

        mgr = create_mgr('a, b: object; f: i8; g: f8')
        new_mgr = mgr.convert()
        _compare(mgr,new_mgr)

        # convert
        mgr = create_mgr('a,b,foo: object; f: i8; g: f8')
        mgr.set('a', np.array(['1'] * N, dtype=np.object_))
        mgr.set('b', np.array(['2.'] * N, dtype=np.object_))
        mgr.set('foo', np.array(['foo.'] * N, dtype=np.object_))
        new_mgr = mgr.convert(convert_numeric=True)
        self.assertEquals(new_mgr.get('a').dtype.type, np.int64)
        self.assertEquals(new_mgr.get('b').dtype.type, np.float64)
        self.assertEquals(new_mgr.get('foo').dtype.type, np.object_)
        self.assertEquals(new_mgr.get('f').dtype.type, np.int64)
        self.assertEquals(new_mgr.get('g').dtype.type, np.float64)

        mgr = create_mgr('a,b,foo: object; f: i4; bool: bool; dt: datetime;'
                         'i: i8; g: f8; h: f2')
        mgr.set('a', np.array(['1'] * N, dtype=np.object_))
        mgr.set('b', np.array(['2.'] * N, dtype=np.object_))
        mgr.set('foo', np.array(['foo.'] * N, dtype=np.object_))
        new_mgr = mgr.convert(convert_numeric=True)
        self.assertEquals(new_mgr.get('a').dtype.type, np.int64)
        self.assertEquals(new_mgr.get('b').dtype.type, np.float64)
        self.assertEquals(new_mgr.get('foo').dtype.type, np.object_)
        self.assertEquals(new_mgr.get('f').dtype.type, np.int32)
        self.assertEquals(new_mgr.get('bool').dtype.type, np.bool_)
        self.assertEquals(new_mgr.get('dt').dtype.type, np.datetime64)
        self.assertEquals(new_mgr.get('i').dtype.type, np.int64)
        self.assertEquals(new_mgr.get('g').dtype.type, np.float64)
        self.assertEquals(new_mgr.get('h').dtype.type, np.float16)

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
        assert_almost_equal(cons.blocks[0].ref_locs,
                            np.arange(len(cons.items)))

    def test_reindex_index(self):
        pass

    def test_reindex_items(self):
        # mgr is not consolidated, f8 & f8-2 blocks
        mgr = create_mgr('a: f8; b: i8; c: f8; d: i8; e: f8;'
                         'f: bool; g: f8-2')

        reindexed = mgr.reindex_axis(['g', 'c', 'a', 'd'], axis=0)
        self.assertEquals(reindexed.nblocks, 2)
        assert_almost_equal(reindexed.items, ['g', 'c', 'a', 'd'])
        assert_almost_equal(mgr.get('g'), reindexed.get('g'))
        assert_almost_equal(mgr.get('c'), reindexed.get('c'))
        assert_almost_equal(mgr.get('a'), reindexed.get('a'))
        assert_almost_equal(mgr.get('d'), reindexed.get('d'))

    def test_multiindex_xs(self):
        mgr = create_mgr('a,b,c: f8; d,e,f: i8')

        index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'],
                                   ['one', 'two', 'three']],
                           labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                                   [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                           names=['first', 'second'])

        mgr.set_axis(1, index)
        result = mgr.xs('bar', axis=1)
        self.assertEqual(result.shape, (6, 2))
        self.assertEqual(result.axes[1][0], ('bar', 'one'))
        self.assertEqual(result.axes[1][1], ('bar', 'two'))

    def test_get_numeric_data(self):
        mgr = create_mgr('int: int; float: float; complex: complex;'
                         'str: object; bool: bool; obj: object; dt: datetime',
                         num_rows=3)
        mgr.set('obj', np.array([1, 2, 3], dtype=np.object_))

        numeric = mgr.get_numeric_data()
        assert_almost_equal(numeric.items, ['int', 'float', 'complex', 'bool'])
        assert_almost_equal(mgr.get('float'), numeric.get('float'))

        # Check sharing
        numeric.set('float', np.array([100., 200., 300.]))
        assert_almost_equal(mgr.get('float'), np.array([100., 200., 300.]))

        numeric2 = mgr.get_numeric_data(copy=True)
        assert_almost_equal(numeric.items, ['int', 'float', 'complex', 'bool'])
        numeric2.set('float', np.array([1000., 2000., 3000.]))
        assert_almost_equal(mgr.get('float'), np.array([100., 200., 300.]))

    def test_get_bool_data(self):
        mgr = create_mgr('int: int; float: float; complex: complex;'
                         'str: object; bool: bool; obj: object; dt: datetime',
                         num_rows=3)
        mgr.set('obj', np.array([True, False, True], dtype=np.object_))

        bools = mgr.get_bool_data()
        assert_almost_equal(bools.items, ['bool'])
        assert_almost_equal(mgr.get('bool'), bools.get('bool'))

        bools.set('bool', np.array([True, False, True]))
        assert_almost_equal(mgr.get('bool'), [True, False, True])

        # Check sharing
        bools2 = mgr.get_bool_data(copy=True)
        bools2.set('bool', np.array([False, True, False]))
        assert_almost_equal(mgr.get('bool'), [True, False, True])

    def test_unicode_repr_doesnt_raise(self):
        str_repr = repr(create_mgr(u('b,\u05d0: object')))

    def test_missing_unicode_key(self):
        df = DataFrame({"a": [1]})
        try:
            df.ix[:, u("\u05d0")]  # should not raise UnicodeEncodeError
        except KeyError:
            pass  # this is the expected exception

    def test_equals(self):
        # unique items
        bm1 = create_mgr('a,b,c: i8-1; d,e,f: i8-2')
        bm2 = BlockManager(bm1.blocks[::-1], bm1.axes)
        self.assertTrue(bm1.equals(bm2))

        bm1 = create_mgr('a,a,a: i8-1; b,b,b: i8-2')
        bm2 = BlockManager(bm1.blocks[::-1], bm1.axes)
        self.assertTrue(bm1.equals(bm2))

if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
