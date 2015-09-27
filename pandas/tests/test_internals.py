# -*- coding: utf-8 -*-
# pylint: disable=W0102

from datetime import datetime, date

import nose
import numpy as np

import re
import itertools
from pandas import Index, MultiIndex, DataFrame, DatetimeIndex, Series, Categorical
from pandas.compat import OrderedDict, lrange
from pandas.sparse.array import SparseArray
from pandas.core.internals import (BlockPlacement, SingleBlockManager, make_block,
                                   BlockManager)
import pandas.core.common as com
import pandas.core.internals as internals
import pandas.util.testing as tm
import pandas as pd
from pandas.util.testing import (
    assert_almost_equal, assert_frame_equal, randn, assert_series_equal)
from pandas.compat import zip, u


def assert_block_equal(left, right):
    assert_almost_equal(left.values, right.values)
    assert(left.dtype == right.dtype)
    assert_almost_equal(left.mgr_locs, right.mgr_locs)


def get_numeric_mat(shape):
    arr = np.arange(shape[0])
    return np.lib.stride_tricks.as_strided(
        x=arr, shape=shape,
        strides=(arr.itemsize,) + (0,) * (len(shape) - 1)).copy()


N = 10


def create_block(typestr, placement, item_shape=None, num_offset=0):
    """
    Supported typestr:

        * float, f8, f4, f2
        * int, i8, i4, i2, i1
        * uint, u8, u4, u2, u1
        * complex, c16, c8
        * bool
        * object, string, O
        * datetime, dt, M8[ns], M8[ns, tz]
        * timedelta, td, m8[ns]
        * sparse (SparseArray with fill_value=0.0)
        * sparse_na (SparseArray with fill_value=np.nan)
        * category, category2

    """
    placement = BlockPlacement(placement)
    num_items = len(placement)

    if item_shape is None:
        item_shape = (N,)

    shape = (num_items,) + item_shape

    mat = get_numeric_mat(shape)

    if typestr in ('float', 'f8', 'f4', 'f2',
                   'int', 'i8', 'i4', 'i2', 'i1',
                   'uint', 'u8', 'u4', 'u2', 'u1'):
        values = mat.astype(typestr) + num_offset
    elif typestr in ('complex', 'c16', 'c8'):
        values = 1.j * (mat.astype(typestr) + num_offset)
    elif typestr in ('object', 'string', 'O'):
        values = np.reshape(['A%d' % i for i in mat.ravel() + num_offset],
                            shape)
    elif typestr in ('b','bool',):
        values = np.ones(shape, dtype=np.bool_)
    elif typestr in ('datetime', 'dt', 'M8[ns]'):
        values = (mat * 1e9).astype('M8[ns]')
    elif typestr.startswith('M8[ns'):
        # datetime with tz
        m = re.search('M8\[ns,\s*(\w+\/?\w*)\]', typestr)
        assert m is not None, "incompatible typestr -> {0}".format(typestr)
        tz = m.groups()[0]
        assert num_items == 1, "must have only 1 num items for a tz-aware"
        values = DatetimeIndex(np.arange(N) * 1e9, tz=tz)
    elif typestr in ('timedelta', 'td', 'm8[ns]'):
        values = (mat * 1).astype('m8[ns]')
    elif typestr in ('category',):
        values = Categorical([1,1,2,2,3,3,3,3,4,4])
    elif typestr in ('category2',):
        values = Categorical(['a','a','a','a','b','b','c','c','c','d'])
    elif typestr in ('sparse', 'sparse_na'):
        # FIXME: doesn't support num_rows != 10
        assert shape[-1] == 10
        assert all(s == 1 for s in shape[:-1])
        if typestr.endswith('_na'):
            fill_value = np.nan
        else:
            fill_value = 0.0
        values = SparseArray([fill_value, fill_value, 1, 2, 3, fill_value,
                              4, 5, fill_value, 6], fill_value=fill_value)
        arr = values.sp_values.view()
        arr += (num_offset - 1)
    else:
        raise ValueError('Unsupported typestr: "%s"' % typestr)

    return make_block(values, placement=placement, ndim=len(shape))


def create_single_mgr(typestr, num_rows=None):
    if num_rows is None:
        num_rows = N

    return SingleBlockManager(
        create_block(typestr, placement=slice(0, num_rows), item_shape=()),
        np.arange(num_rows))


def create_mgr(descr, item_shape=None):
    """
    Construct BlockManager from string description.

    String description syntax looks similar to np.matrix initializer.  It looks
    like this::

        a,b,c: f8; d,e,f: i8

    Rules are rather simple:

    * see list of supported datatypes in `create_block` method
    * components are semicolon-separated
    * each component is `NAME,NAME,NAME: DTYPE_ID`
    * whitespace around colons & semicolons are removed
    * components with same DTYPE_ID are combined into single block
    * to force multiple blocks with same dtype, use '-SUFFIX'::

        'a:f8-1; b:f8-2; c:f8-foobar'

    """
    if item_shape is None:
        item_shape = (N,)

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
        blocks.append(create_block(typestr, placement, item_shape=item_shape,
                                   num_offset=num_offset,))
        num_offset += len(placement)

    return BlockManager(sorted(blocks, key=lambda b: b.mgr_locs[0]),
                        [mgr_items] + [np.arange(n) for n in item_shape])



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

        def _check(blk):
            assert_block_equal(self.round_trip_pickle(blk), blk)

        _check(self.fblock)
        _check(self.cblock)
        _check(self.oblock)
        _check(self.bool_block)

    def test_mgr_locs(self):
        assert_almost_equal(self.fblock.mgr_locs, [0, 2, 4])

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
        assert_almost_equal(merged.mgr_locs, [0, 1, 2, 3])
        assert_almost_equal(merged.values[[0, 2]], avals)
        assert_almost_equal(merged.values[[1, 3]], bvals)

        # TODO: merge with mixed type?

    def test_copy(self):
        cop = self.fblock.copy()
        self.assertIsNot(cop, self.fblock)
        assert_block_equal(self.fblock, cop)

    def test_reindex_index(self):
        pass

    def test_reindex_cast(self):
        pass

    def test_insert(self):
        pass

    def test_delete(self):
        newb = self.fblock.copy()
        newb.delete(0)
        assert_almost_equal(newb.mgr_locs, [2, 4])
        self.assertTrue((newb.values[0] == 1).all())

        newb = self.fblock.copy()
        newb.delete(1)
        assert_almost_equal(newb.mgr_locs, [0, 4])
        self.assertTrue((newb.values[1] == 2).all())

        newb = self.fblock.copy()
        newb.delete(2)
        assert_almost_equal(newb.mgr_locs, [0, 2])
        self.assertTrue((newb.values[1] == 1).all())

        newb = self.fblock.copy()
        self.assertRaises(Exception, newb.delete, 3)

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

    def test_get(self):
        pass

    def test_set(self):
        pass

    def test_fillna(self):
        pass

    def test_repr(self):
        pass


class TestDatetimeBlock(tm.TestCase):
    _multiprocess_can_split_ = True

    def test_try_coerce_arg(self):
        block = create_block('datetime', [0])

        # coerce None
        none_coerced = block._try_coerce_args(block.values, None)[1]
        self.assertTrue(pd.Timestamp(none_coerced) is pd.NaT)

        # coerce different types of date bojects
        vals = (np.datetime64('2010-10-10'),
                datetime(2010, 10, 10),
                date(2010, 10, 10))
        for val in vals:
            coerced = block._try_coerce_args(block.values, val)[1]
            self.assertEqual(np.int64, type(coerced))
            self.assertEqual(pd.Timestamp('2010-10-10'), pd.Timestamp(coerced))


class TestBlockManager(tm.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        self.mgr = create_mgr('a: f8; b: object; c: f8; d: object; e: f8;'
                              'f: bool; g: i8; h: complex')

    def test_constructor_corner(self):
        pass

    def test_attrs(self):
        mgr = create_mgr('a,b,c: f8-1; d,e,f: f8-2')
        self.assertEqual(mgr.nblocks, 2)
        self.assertEqual(len(mgr), 6)

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

        blocks[0].mgr_locs = np.array([0])
        blocks[1].mgr_locs = np.array([0])
        # test trying to create block manager with overlapping ref locs
        self.assertRaises(AssertionError, BlockManager, blocks, axes)

        blocks[0].mgr_locs = np.array([0])
        blocks[1].mgr_locs = np.array([1])
        mgr = BlockManager(blocks, axes)
        mgr.iget(1)

    def test_contains(self):
        self.assertIn('a', self.mgr)
        self.assertNotIn('baz', self.mgr)

    def test_pickle(self):

        mgr2 = self.round_trip_pickle(self.mgr)
        assert_frame_equal(DataFrame(self.mgr), DataFrame(mgr2))

        # share ref_items
        # self.assertIs(mgr2.blocks[0].ref_items, mgr2.blocks[1].ref_items)

        # GH2431
        self.assertTrue(hasattr(mgr2, "_is_consolidated"))
        self.assertTrue(hasattr(mgr2, "_known_consolidated"))

        # reset to False on load
        self.assertFalse(mgr2._is_consolidated)
        self.assertFalse(mgr2._known_consolidated)

    def test_non_unique_pickle(self):

        mgr = create_mgr('a,a,a:f8')
        mgr2 = self.round_trip_pickle(mgr)
        assert_frame_equal(DataFrame(mgr), DataFrame(mgr2))

        mgr = create_mgr('a: f8; a: i8')
        mgr2 = self.round_trip_pickle(mgr)
        assert_frame_equal(DataFrame(mgr), DataFrame(mgr2))

    def test_categorical_block_pickle(self):
        mgr = create_mgr('a: category')
        mgr2 = self.round_trip_pickle(mgr)
        assert_frame_equal(DataFrame(mgr), DataFrame(mgr2))

        smgr = create_single_mgr('category')
        smgr2 = self.round_trip_pickle(smgr)
        assert_series_equal(Series(smgr), Series(smgr2))

    def test_get_scalar(self):
        for item in self.mgr.items:
            for i, index in enumerate(self.mgr.axes[1]):
                res = self.mgr.get_scalar((item, index))
                exp = self.mgr.get(item, fastpath=False)[i]
                assert_almost_equal(res, exp)
                exp = self.mgr.get(item).internal_values()[i]
                assert_almost_equal(res, exp)

    def test_get(self):
        cols = Index(list('abc'))
        values = np.random.rand(3, 3)
        block = make_block(values=values.copy(),
                           placement=np.arange(3))
        mgr = BlockManager(blocks=[block], axes=[cols, np.arange(3)])

        assert_almost_equal(mgr.get('a', fastpath=False), values[0])
        assert_almost_equal(mgr.get('b', fastpath=False), values[1])
        assert_almost_equal(mgr.get('c', fastpath=False), values[2])
        assert_almost_equal(mgr.get('a').internal_values(), values[0])
        assert_almost_equal(mgr.get('b').internal_values(), values[1])
        assert_almost_equal(mgr.get('c').internal_values(), values[2])

    def test_set(self):
        mgr = create_mgr('a,b,c: int', item_shape=(3,))

        mgr.set('d', np.array(['foo'] * 3))
        mgr.set('b', np.array(['bar'] * 3))
        assert_almost_equal(mgr.get('a').internal_values(), [0] * 3)
        assert_almost_equal(mgr.get('b').internal_values(), ['bar'] * 3)
        assert_almost_equal(mgr.get('c').internal_values(), [2] * 3)
        assert_almost_equal(mgr.get('d').internal_values(), ['foo'] * 3)

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

    def test_set_change_dtype_slice(self): # GH8850
        cols = MultiIndex.from_tuples([('1st','a'), ('2nd','b'), ('3rd','c')])
        df = DataFrame([[1.0, 2, 3], [4.0, 5, 6]], columns=cols)
        df['2nd'] = df['2nd'] * 2.0

        self.assertEqual(sorted(df.blocks.keys()), ['float64', 'int64'])
        assert_frame_equal(df.blocks['float64'],
                DataFrame([[1.0, 4.0], [4.0, 10.0]], columns=cols[:2]))
        assert_frame_equal(df.blocks['int64'],
                DataFrame([[3], [6]], columns=cols[2:]))

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

    def test_as_matrix_datetime_tz(self):
        mgr = create_mgr('h: M8[ns, US/Eastern]; g: M8[ns, CET]')
        self.assertEqual(mgr.get('h').dtype, 'datetime64[ns, US/Eastern]')
        self.assertEqual(mgr.get('g').dtype, 'datetime64[ns, CET]')
        self.assertEqual(mgr.as_matrix().dtype, 'object')

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
        new_mgr = mgr.convert(numeric=True)
        self.assertEqual(new_mgr.get('a').dtype, np.int64)
        self.assertEqual(new_mgr.get('b').dtype, np.float64)
        self.assertEqual(new_mgr.get('foo').dtype, np.object_)
        self.assertEqual(new_mgr.get('f').dtype, np.int64)
        self.assertEqual(new_mgr.get('g').dtype, np.float64)

        mgr = create_mgr('a,b,foo: object; f: i4; bool: bool; dt: datetime;'
                         'i: i8; g: f8; h: f2')
        mgr.set('a', np.array(['1'] * N, dtype=np.object_))
        mgr.set('b', np.array(['2.'] * N, dtype=np.object_))
        mgr.set('foo', np.array(['foo.'] * N, dtype=np.object_))
        new_mgr = mgr.convert(numeric=True)
        self.assertEqual(new_mgr.get('a').dtype, np.int64)
        self.assertEqual(new_mgr.get('b').dtype, np.float64)
        self.assertEqual(new_mgr.get('foo').dtype, np.object_)
        self.assertEqual(new_mgr.get('f').dtype, np.int32)
        self.assertEqual(new_mgr.get('bool').dtype, np.bool_)
        self.assertEqual(new_mgr.get('dt').dtype.type, np.datetime64)
        self.assertEqual(new_mgr.get('i').dtype, np.int64)
        self.assertEqual(new_mgr.get('g').dtype, np.float64)
        self.assertEqual(new_mgr.get('h').dtype, np.float16)

    def test_interleave(self):


        # self
        for dtype in ['f8','i8','object','bool','complex','M8[ns]','m8[ns]']:
            mgr = create_mgr('a: {0}'.format(dtype))
            self.assertEqual(mgr.as_matrix().dtype,dtype)
            mgr = create_mgr('a: {0}; b: {0}'.format(dtype))
            self.assertEqual(mgr.as_matrix().dtype,dtype)

        # will be converted according the actual dtype of the underlying
        mgr = create_mgr('a: category')
        self.assertEqual(mgr.as_matrix().dtype,'i8')
        mgr = create_mgr('a: category; b: category')
        self.assertEqual(mgr.as_matrix().dtype,'i8'),
        mgr = create_mgr('a: category; b: category2')
        self.assertEqual(mgr.as_matrix().dtype,'object')
        mgr = create_mgr('a: category2')
        self.assertEqual(mgr.as_matrix().dtype,'object')
        mgr = create_mgr('a: category2; b: category2')
        self.assertEqual(mgr.as_matrix().dtype,'object')

        # combinations
        mgr = create_mgr('a: f8')
        self.assertEqual(mgr.as_matrix().dtype,'f8')
        mgr = create_mgr('a: f8; b: i8')
        self.assertEqual(mgr.as_matrix().dtype,'f8')
        mgr = create_mgr('a: f4; b: i8')
        self.assertEqual(mgr.as_matrix().dtype,'f4')
        mgr = create_mgr('a: f4; b: i8; d: object')
        self.assertEqual(mgr.as_matrix().dtype,'object')
        mgr = create_mgr('a: bool; b: i8')
        self.assertEqual(mgr.as_matrix().dtype,'object')
        mgr = create_mgr('a: complex')
        self.assertEqual(mgr.as_matrix().dtype,'complex')
        mgr = create_mgr('a: f8; b: category')
        self.assertEqual(mgr.as_matrix().dtype,'object')
        mgr = create_mgr('a: M8[ns]; b: category')
        self.assertEqual(mgr.as_matrix().dtype,'object')
        mgr = create_mgr('a: M8[ns]; b: bool')
        self.assertEqual(mgr.as_matrix().dtype,'object')
        mgr = create_mgr('a: M8[ns]; b: i8')
        self.assertEqual(mgr.as_matrix().dtype,'object')
        mgr = create_mgr('a: m8[ns]; b: bool')
        self.assertEqual(mgr.as_matrix().dtype,'object')
        mgr = create_mgr('a: m8[ns]; b: i8')
        self.assertEqual(mgr.as_matrix().dtype,'object')
        mgr = create_mgr('a: M8[ns]; b: m8[ns]')
        self.assertEqual(mgr.as_matrix().dtype,'object')

    def test_interleave_non_unique_cols(self):
        df = DataFrame([
            [pd.Timestamp('20130101'), 3.5],
            [pd.Timestamp('20130102'), 4.5]],
            columns=['x', 'x'],
            index=[1, 2])

        df_unique = df.copy()
        df_unique.columns = ['x', 'y']
        self.assertEqual(df_unique.values.shape, df.values.shape)
        tm.assert_numpy_array_equal(df_unique.values[0], df.values[0])
        tm.assert_numpy_array_equal(df_unique.values[1], df.values[1])

    def test_consolidate(self):
        pass

    def test_consolidate_ordering_issues(self):
        self.mgr.set('f', randn(N))
        self.mgr.set('d', randn(N))
        self.mgr.set('b', randn(N))
        self.mgr.set('g', randn(N))
        self.mgr.set('h', randn(N))

        cons = self.mgr.consolidate()
        self.assertEqual(cons.nblocks, 1)
        assert_almost_equal(cons.blocks[0].mgr_locs,
                            np.arange(len(cons.items)))

    def test_reindex_index(self):
        pass

    def test_reindex_items(self):
        # mgr is not consolidated, f8 & f8-2 blocks
        mgr = create_mgr('a: f8; b: i8; c: f8; d: i8; e: f8;'
                         'f: bool; g: f8-2')

        reindexed = mgr.reindex_axis(['g', 'c', 'a', 'd'], axis=0)
        self.assertEqual(reindexed.nblocks, 2)
        assert_almost_equal(reindexed.items, ['g', 'c', 'a', 'd'])
        assert_almost_equal(mgr.get('g',fastpath=False), reindexed.get('g',fastpath=False))
        assert_almost_equal(mgr.get('c',fastpath=False), reindexed.get('c',fastpath=False))
        assert_almost_equal(mgr.get('a',fastpath=False), reindexed.get('a',fastpath=False))
        assert_almost_equal(mgr.get('d',fastpath=False), reindexed.get('d',fastpath=False))
        assert_almost_equal(mgr.get('g').internal_values(), reindexed.get('g').internal_values())
        assert_almost_equal(mgr.get('c').internal_values(), reindexed.get('c').internal_values())
        assert_almost_equal(mgr.get('a').internal_values(), reindexed.get('a').internal_values())
        assert_almost_equal(mgr.get('d').internal_values(), reindexed.get('d').internal_values())

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
                         item_shape=(3,))
        mgr.set('obj', np.array([1, 2, 3], dtype=np.object_))

        numeric = mgr.get_numeric_data()
        assert_almost_equal(numeric.items, ['int', 'float', 'complex', 'bool'])
        assert_almost_equal(mgr.get('float',fastpath=False), numeric.get('float',fastpath=False))
        assert_almost_equal(mgr.get('float').internal_values(), numeric.get('float').internal_values())

        # Check sharing
        numeric.set('float', np.array([100., 200., 300.]))
        assert_almost_equal(mgr.get('float',fastpath=False), np.array([100., 200., 300.]))
        assert_almost_equal(mgr.get('float').internal_values(), np.array([100., 200., 300.]))

        numeric2 = mgr.get_numeric_data(copy=True)
        assert_almost_equal(numeric.items, ['int', 'float', 'complex', 'bool'])
        numeric2.set('float', np.array([1000., 2000., 3000.]))
        assert_almost_equal(mgr.get('float',fastpath=False), np.array([100., 200., 300.]))
        assert_almost_equal(mgr.get('float').internal_values(), np.array([100., 200., 300.]))

    def test_get_bool_data(self):
        mgr = create_mgr('int: int; float: float; complex: complex;'
                         'str: object; bool: bool; obj: object; dt: datetime',
                         item_shape=(3,))
        mgr.set('obj', np.array([True, False, True], dtype=np.object_))

        bools = mgr.get_bool_data()
        assert_almost_equal(bools.items, ['bool'])
        assert_almost_equal(mgr.get('bool',fastpath=False), bools.get('bool',fastpath=False))
        assert_almost_equal(mgr.get('bool').internal_values(), bools.get('bool').internal_values())

        bools.set('bool', np.array([True, False, True]))
        assert_almost_equal(mgr.get('bool',fastpath=False), [True, False, True])
        assert_almost_equal(mgr.get('bool').internal_values(), [True, False, True])

        # Check sharing
        bools2 = mgr.get_bool_data(copy=True)
        bools2.set('bool', np.array([False, True, False]))
        assert_almost_equal(mgr.get('bool',fastpath=False), [True, False, True])
        assert_almost_equal(mgr.get('bool').internal_values(), [True, False, True])

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

    def test_equals_block_order_different_dtypes(self):
        # GH 9330

        mgr_strings = [
            "a:i8;b:f8", # basic case
            "a:i8;b:f8;c:c8;d:b", # many types
            "a:i8;e:dt;f:td;g:string", # more types
            "a:i8;b:category;c:category2;d:category2", # categories
            "c:sparse;d:sparse_na;b:f8", # sparse
            ]

        for mgr_string in mgr_strings:
            bm = create_mgr(mgr_string)
            block_perms = itertools.permutations(bm.blocks)
            for bm_perm in block_perms:
                bm_this = BlockManager(bm_perm, bm.axes)
                self.assertTrue(bm.equals(bm_this))
                self.assertTrue(bm_this.equals(bm))

    def test_single_mgr_ctor(self):
        mgr = create_single_mgr('f8', num_rows=5)
        self.assertEqual(mgr.as_matrix().tolist(), [0., 1., 2., 3., 4.])


class TestIndexing(object):
    # Nosetests-style data-driven tests.
    #
    # This test applies different indexing routines to block managers and
    # compares the outcome to the result of same operations on np.ndarray.
    #
    # NOTE: sparse (SparseBlock with fill_value != np.nan) fail a lot of tests
    #       and are disabled.

    MANAGERS = [
        create_single_mgr('f8', N),
        create_single_mgr('i8', N),
        #create_single_mgr('sparse', N),
        create_single_mgr('sparse_na', N),

        # 2-dim
        create_mgr('a,b,c,d,e,f: f8', item_shape=(N,)),
        create_mgr('a,b,c,d,e,f: i8', item_shape=(N,)),
        create_mgr('a,b: f8; c,d: i8; e,f: string', item_shape=(N,)),
        create_mgr('a,b: f8; c,d: i8; e,f: f8', item_shape=(N,)),
        #create_mgr('a: sparse', item_shape=(N,)),
        create_mgr('a: sparse_na', item_shape=(N,)),

        # 3-dim
        create_mgr('a,b,c,d,e,f: f8', item_shape=(N, N)),
        create_mgr('a,b,c,d,e,f: i8', item_shape=(N, N)),
        create_mgr('a,b: f8; c,d: i8; e,f: string', item_shape=(N, N)),
        create_mgr('a,b: f8; c,d: i8; e,f: f8', item_shape=(N, N)),
        # create_mgr('a: sparse', item_shape=(1, N)),
    ]

    # MANAGERS = [MANAGERS[6]]

    def test_get_slice(self):
        def assert_slice_ok(mgr, axis, slobj):
            # import pudb; pudb.set_trace()
            mat = mgr.as_matrix()

            # we maybe using an ndarray to test slicing and
            # might not be the full length of the axis
            if isinstance(slobj, np.ndarray):
                ax = mgr.axes[axis]
                if len(ax) and len(slobj) and len(slobj) != len(ax):
                    slobj = np.concatenate([slobj, np.zeros(len(ax)-len(slobj),dtype=bool)])
            sliced = mgr.get_slice(slobj, axis=axis)
            mat_slobj = (slice(None),) * axis + (slobj,)
            assert_almost_equal(mat[mat_slobj], sliced.as_matrix())
            assert_almost_equal(mgr.axes[axis][slobj], sliced.axes[axis])

        for mgr in self.MANAGERS:
            for ax in range(mgr.ndim):
                # slice
                yield assert_slice_ok, mgr, ax, slice(None)
                yield assert_slice_ok, mgr, ax, slice(3)
                yield assert_slice_ok, mgr, ax, slice(100)
                yield assert_slice_ok, mgr, ax, slice(1, 4)
                yield assert_slice_ok, mgr, ax, slice(3, 0, -2)

                # boolean mask
                yield assert_slice_ok, mgr, ax, np.array([], dtype=np.bool_)
                yield (assert_slice_ok, mgr, ax,
                       np.ones(mgr.shape[ax], dtype=np.bool_))
                yield (assert_slice_ok, mgr, ax,
                       np.zeros(mgr.shape[ax], dtype=np.bool_))

                if mgr.shape[ax] >= 3:
                    yield (assert_slice_ok, mgr, ax,
                           np.arange(mgr.shape[ax]) % 3 == 0)
                    yield (assert_slice_ok, mgr, ax,
                           np.array([True, True, False], dtype=np.bool_))

                # fancy indexer
                yield assert_slice_ok, mgr, ax, []
                yield assert_slice_ok, mgr, ax, lrange(mgr.shape[ax])

                if mgr.shape[ax] >= 3:
                    yield assert_slice_ok, mgr, ax, [0, 1, 2]
                    yield assert_slice_ok, mgr, ax, [-1, -2, -3]

    def test_take(self):
        def assert_take_ok(mgr, axis, indexer):
            mat = mgr.as_matrix()
            taken = mgr.take(indexer, axis)
            assert_almost_equal(np.take(mat, indexer, axis),
                                taken.as_matrix())
            assert_almost_equal(mgr.axes[axis].take(indexer),
                                taken.axes[axis])

        for mgr in self.MANAGERS:
            for ax in range(mgr.ndim):
                # take/fancy indexer
                yield assert_take_ok, mgr, ax, []
                yield assert_take_ok, mgr, ax, [0, 0, 0]
                yield assert_take_ok, mgr, ax, lrange(mgr.shape[ax])

                if mgr.shape[ax] >= 3:
                    yield assert_take_ok, mgr, ax, [0, 1, 2]
                    yield assert_take_ok, mgr, ax, [-1, -2, -3]

    def test_reindex_axis(self):
        def assert_reindex_axis_is_ok(mgr, axis, new_labels,
                                      fill_value):
            mat = mgr.as_matrix()
            indexer = mgr.axes[axis].get_indexer_for(new_labels)

            reindexed = mgr.reindex_axis(new_labels, axis,
                                         fill_value=fill_value)
            assert_almost_equal(com.take_nd(mat, indexer, axis,
                                            fill_value=fill_value),
                                reindexed.as_matrix())
            assert_almost_equal(reindexed.axes[axis], new_labels)

        for mgr in self.MANAGERS:
            for ax in range(mgr.ndim):
                for fill_value in (None, np.nan, 100.):
                    yield assert_reindex_axis_is_ok, mgr, ax, [], fill_value
                    yield (assert_reindex_axis_is_ok, mgr, ax,
                           mgr.axes[ax], fill_value)
                    yield (assert_reindex_axis_is_ok, mgr, ax,
                           mgr.axes[ax][[0, 0, 0]], fill_value)
                    yield (assert_reindex_axis_is_ok, mgr, ax,
                           ['foo', 'bar', 'baz'], fill_value)
                    yield (assert_reindex_axis_is_ok, mgr, ax,
                           ['foo', mgr.axes[ax][0], 'baz'], fill_value)

                    if mgr.shape[ax] >= 3:
                        yield (assert_reindex_axis_is_ok, mgr, ax,
                               mgr.axes[ax][:-3], fill_value)
                        yield (assert_reindex_axis_is_ok, mgr, ax,
                               mgr.axes[ax][-3::-1], fill_value)
                        yield (assert_reindex_axis_is_ok, mgr, ax,
                               mgr.axes[ax][[0, 1, 2, 0, 1, 2]], fill_value)

    def test_reindex_indexer(self):
        def assert_reindex_indexer_is_ok(mgr, axis, new_labels, indexer,
                                         fill_value):
            mat = mgr.as_matrix()
            reindexed_mat = com.take_nd(mat, indexer, axis,
                                        fill_value=fill_value)
            reindexed = mgr.reindex_indexer(new_labels, indexer, axis,
                                            fill_value=fill_value)
            assert_almost_equal(reindexed_mat, reindexed.as_matrix())
            assert_almost_equal(reindexed.axes[axis], new_labels)

        for mgr in self.MANAGERS:
            for ax in range(mgr.ndim):
                for fill_value in (None, np.nan, 100.):
                    yield (assert_reindex_indexer_is_ok, mgr, ax,
                           [], [], fill_value)
                    yield (assert_reindex_indexer_is_ok, mgr, ax,
                           mgr.axes[ax], np.arange(mgr.shape[ax]), fill_value)
                    yield (assert_reindex_indexer_is_ok, mgr, ax,
                           ['foo'] * mgr.shape[ax], np.arange(mgr.shape[ax]),
                           fill_value)

                    yield (assert_reindex_indexer_is_ok, mgr, ax,
                           mgr.axes[ax][::-1], np.arange(mgr.shape[ax]),
                           fill_value)
                    yield (assert_reindex_indexer_is_ok, mgr, ax,
                           mgr.axes[ax], np.arange(mgr.shape[ax])[::-1],
                           fill_value)
                    yield (assert_reindex_indexer_is_ok, mgr, ax,
                           ['foo', 'bar', 'baz'], [0, 0, 0], fill_value)
                    yield (assert_reindex_indexer_is_ok, mgr, ax,
                           ['foo', 'bar', 'baz'], [-1, 0, -1], fill_value)
                    yield (assert_reindex_indexer_is_ok, mgr, ax,
                           ['foo', mgr.axes[ax][0], 'baz'], [-1, -1, -1],
                           fill_value)

                    if mgr.shape[ax] >= 3:
                        yield (assert_reindex_indexer_is_ok, mgr, ax,
                               ['foo', 'bar', 'baz'], [0, 1, 2], fill_value)


    # test_get_slice(slice_like, axis)
    # take(indexer, axis)
    # reindex_axis(new_labels, axis)
    # reindex_indexer(new_labels, indexer, axis)


class TestBlockPlacement(tm.TestCase):
    _multiprocess_can_split_ = True

    def test_slice_len(self):
        self.assertEqual(len(BlockPlacement(slice(0, 4))), 4)
        self.assertEqual(len(BlockPlacement(slice(0, 4, 2))), 2)
        self.assertEqual(len(BlockPlacement(slice(0, 3, 2))), 2)

        self.assertEqual(len(BlockPlacement(slice(0, 1, 2))), 1)
        self.assertEqual(len(BlockPlacement(slice(1, 0, -1))), 1)

    def test_zero_step_raises(self):
        self.assertRaises(ValueError, BlockPlacement, slice(1, 1, 0))
        self.assertRaises(ValueError, BlockPlacement, slice(1, 2, 0))

    def test_unbounded_slice_raises(self):
        def assert_unbounded_slice_error(slc):
            # assertRaisesRegexp is not available in py2.6
            # self.assertRaisesRegexp(ValueError, "unbounded slice",
            #                         lambda: BlockPlacement(slc))
            self.assertRaises(ValueError, BlockPlacement, slc)

        assert_unbounded_slice_error(slice(None, None))
        assert_unbounded_slice_error(slice(10, None))
        assert_unbounded_slice_error(slice(None, None, -1))
        assert_unbounded_slice_error(slice(None, 10, -1))

        # These are "unbounded" because negative index will change depending on
        # container shape.
        assert_unbounded_slice_error(slice(-1, None))
        assert_unbounded_slice_error(slice(None, -1))
        assert_unbounded_slice_error(slice(-1, -1))
        assert_unbounded_slice_error(slice(-1, None, -1))
        assert_unbounded_slice_error(slice(None, -1, -1))
        assert_unbounded_slice_error(slice(-1, -1, -1))

    def test_not_slice_like_slices(self):
        def assert_not_slice_like(slc):
            self.assertTrue(not BlockPlacement(slc).is_slice_like)

        assert_not_slice_like(slice(0, 0))
        assert_not_slice_like(slice(100, 0))

        assert_not_slice_like(slice(100, 100, -1))
        assert_not_slice_like(slice(0, 100, -1))

        self.assertTrue(not BlockPlacement(slice(0, 0)).is_slice_like)
        self.assertTrue(not BlockPlacement(slice(100, 100)).is_slice_like)

    def test_array_to_slice_conversion(self):
        def assert_as_slice_equals(arr, slc):
            self.assertEqual(BlockPlacement(arr).as_slice, slc)

        assert_as_slice_equals([0], slice(0, 1, 1))
        assert_as_slice_equals([100], slice(100, 101, 1))

        assert_as_slice_equals([0, 1, 2], slice(0, 3, 1))
        assert_as_slice_equals([0, 5, 10], slice(0, 15, 5))
        assert_as_slice_equals([0, 100], slice(0, 200, 100))

        assert_as_slice_equals([2, 1], slice(2, 0, -1))
        assert_as_slice_equals([2, 1, 0], slice(2, None, -1))
        assert_as_slice_equals([100, 0], slice(100, None, -100))

    def test_not_slice_like_arrays(self):
        def assert_not_slice_like(arr):
            self.assertTrue(not BlockPlacement(arr).is_slice_like)

        assert_not_slice_like([])
        assert_not_slice_like([-1])
        assert_not_slice_like([-1, -2, -3])
        assert_not_slice_like([-10])
        assert_not_slice_like([-1])
        assert_not_slice_like([-1, 0, 1, 2])
        assert_not_slice_like([-2, 0, 2, 4])
        assert_not_slice_like([1, 0, -1])
        assert_not_slice_like([1, 1, 1])

    def test_slice_iter(self):
        self.assertEqual(list(BlockPlacement(slice(0, 3))), [0, 1, 2])
        self.assertEqual(list(BlockPlacement(slice(0, 0))), [])
        self.assertEqual(list(BlockPlacement(slice(3, 0))), [])

        self.assertEqual(list(BlockPlacement(slice(3, 0, -1))), [3, 2, 1])
        self.assertEqual(list(BlockPlacement(slice(3, None, -1))),
                          [3, 2, 1, 0])

    def test_slice_to_array_conversion(self):
        def assert_as_array_equals(slc, asarray):
            tm.assert_numpy_array_equal(
                BlockPlacement(slc).as_array,
                np.asarray(asarray))

        assert_as_array_equals(slice(0, 3), [0, 1, 2])
        assert_as_array_equals(slice(0, 0), [])
        assert_as_array_equals(slice(3, 0), [])

        assert_as_array_equals(slice(3, 0, -1), [3, 2, 1])
        assert_as_array_equals(slice(3, None, -1), [3, 2, 1, 0])
        assert_as_array_equals(slice(31, None, -10), [31, 21, 11, 1])

    def test_blockplacement_add(self):
        bpl = BlockPlacement(slice(0, 5))
        self.assertEqual(bpl.add(1).as_slice, slice(1, 6, 1))
        self.assertEqual(bpl.add(np.arange(5)).as_slice,
                          slice(0, 10, 2))
        self.assertEqual(list(bpl.add(np.arange(5, 0, -1))),
                          [5, 5, 5, 5, 5])

    def test_blockplacement_add_int(self):
        def assert_add_equals(val, inc, result):
            self.assertEqual(list(BlockPlacement(val).add(inc)),
                              result)

        assert_add_equals(slice(0, 0), 0, [])
        assert_add_equals(slice(1, 4), 0, [1, 2, 3])
        assert_add_equals(slice(3, 0, -1), 0, [3, 2, 1])
        assert_add_equals(slice(2, None, -1), 0, [2, 1, 0])
        assert_add_equals([1, 2, 4], 0, [1, 2, 4])

        assert_add_equals(slice(0, 0), 10, [])
        assert_add_equals(slice(1, 4), 10, [11, 12, 13])
        assert_add_equals(slice(3, 0, -1), 10, [13, 12, 11])
        assert_add_equals(slice(2, None, -1), 10, [12, 11, 10])
        assert_add_equals([1, 2, 4], 10, [11, 12, 14])

        assert_add_equals(slice(0, 0), -1, [])
        assert_add_equals(slice(1, 4), -1, [0, 1, 2])
        assert_add_equals(slice(3, 0, -1), -1, [2, 1, 0])
        assert_add_equals([1, 2, 4], -1, [0, 1, 3])

        self.assertRaises(ValueError,
                          lambda: BlockPlacement(slice(1, 4)).add(-10))
        self.assertRaises(ValueError,
                          lambda: BlockPlacement([1, 2, 4]).add(-10))
        self.assertRaises(ValueError,
                          lambda: BlockPlacement(slice(2, None, -1)).add(-1))

    # def test_blockplacement_array_add(self):

    #     assert_add_equals(slice(0, 2), [0, 1, 1], [0, 2, 3])
    #     assert_add_equals(slice(2, None, -1), [1, 1, 0], [3, 2, 0])


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
