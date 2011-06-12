import itertools

from numpy import nan
import numpy as np

from pandas.core.index import Index, NULL_INDEX
from pandas.core.common import _ensure_index
from pandas.core.series import Series
import pandas.core.common as common

class Block(object):
    """
    Canonical unit of homogeneous dtype contained in DataMatrix

    Index-ignorant; let the container take care of that
    """
    def __init__(self, values, columns):
        values = _convert_if_1d(values)
        if issubclass(values.dtype.type, basestring):
            values = np.array(values, dtype=object)

        self.values = values
        self.columns = _ensure_index(columns)
        assert(len(self.columns) == values.shape[1])

    def __repr__(self):
        x, y = self.shape
        name = type(self).__name__
        return '%s: %s, %d x %d, dtype %s' % (name, self.columns,
                                              x, y, self.dtype)

    def __contains__(self, col):
        return col in self.columns

    def __len__(self):
        return len(self.values)

    def __getstate__(self):
        return (np.asarray(self.columns),
                self.values)

    def __setstate__(self, state):
        columns, values = state
        self.columns = _ensure_index(columns)
        self.values = values

    @property
    def shape(self):
        return self.values.shape

    @property
    def dtype(self):
        return self.values.dtype

    def copy(self):
        return make_block(self.values.copy(), self.columns)

    def merge(self, other):
        return _merge_blocks([self, other])

    def reindex_index(self, indexer, notmask, needs_masking):
        """
        Reindex using pre-computed indexer information
        """
        new_values = self.values.take(indexer, axis=0)
        if needs_masking:
            new_values = _cast_if_bool_int(new_values)
            common.null_out_axis(new_values, notmask, 0)
        return make_block(new_values, self.columns)

    def reindex_columns(self, new_columns):
        """

        """
        indexer, mask = self.columns.get_indexer(new_columns)
        new_values = self.values.take(indexer, axis=1)

        notmask = -mask
        if len(mask) > 0 and notmask.any():
            new_values = _cast_if_bool_int(new_values)
            common.null_out_axis(new_values, notmask, 1)

        return make_block(new_values, new_columns)

    def reindex_columns_from(self, columns):
        """
        Reindex to only those columns contained in the input set of columns

        E.g. if you have ['a', 'b'], and the input columns is ['b', 'c', 'd'],
        then the resulting columns will be ['b']

        Returns
        -------
        reindexed : Block
        """
        indexer, mask = self.columns.get_indexer(columns)
        masked_idx = indexer[mask]
        new_values = self.values.take(masked_idx, axis=1)
        new_columns = self.columns.take(masked_idx)
        return make_block(new_values, new_columns)

    def insert(self, col, value, loc=None):
        """
        Insert new column into Block, return new Block

        Returns
        -------
        y : Block (new object)
        """
        assert(col not in self.columns)
        if loc is None:
            loc = len(self.columns)

        new_columns = _insert_into_columns(self.columns, col, loc)
        new_values = _insert_into_values(self.values, value, loc)
        return make_block(new_values, new_columns)

    def get(self, col):
        loc = self.columns.get_loc(col)
        return self.values[:, loc]

    def set(self, col, value):
        """
        Modify Block in-place with new column value

        Returns
        -------
        None
        """
        loc = self.columns.get_loc(col)
        self.values[:, loc] = value

    def delete(self, col):
        """
        Returns
        -------
        y : Block (new object)
        """
        loc = self.columns.get_loc(col)
        new_columns = _delete_from_columns(self.columns, loc)
        new_values = _delete_from_values(self.values, loc)
        return make_block(new_values, new_columns)

    def fillna(self, value):
        new_values = self.values.copy()
        mask = common.isnull(new_values.ravel())
        new_values.flat[mask] = value
        return make_block(new_values, self.columns)

def _cast_if_bool_int(values):
    if issubclass(values.dtype.type, np.int_):
        values = values.astype(float)
    elif issubclass(values.dtype.type, np.bool_):
        values = values.astype(object)
    return values

def _insert_into_columns(columns, col, loc):
    columns = np.asarray(columns)
    new_columns = np.insert(columns, loc, col)
    return Index(new_columns)

def _insert_into_values(values, new_vec, loc):
    return np.insert(values, loc, new_vec, 1)

def _delete_from_columns(columns, loc):
    columns = np.asarray(columns)
    new_columns = np.delete(columns, loc)
    return Index(new_columns)

def _delete_from_values(values, loc):
    return np.delete(values, loc, 1)

def _convert_if_1d(values):
    if values.ndim == 1:
        values = np.atleast_2d(values).T

    return values

#-------------------------------------------------------------------------------
# Is this even possible?

class FloatBlock(Block):
    pass

class IntBlock(Block):
    pass

class BoolBlock(Block):
    pass

class ObjectBlock(Block):
    pass

def make_block(values, columns):
    dtype = values.dtype
    vtype = dtype.type

    if issubclass(vtype, np.floating):
        klass = FloatBlock
    elif issubclass(vtype, np.integer):
        klass = IntBlock
    elif dtype == np.bool_:
        klass = BoolBlock
    else:
        klass = ObjectBlock

    return klass(values, columns)

# TODO: flexible with index=None and/or columns=None

class BlockManager(object):
    """
    Manage a bunch of labeled 2D mixed-type ndarrays. Essentially it's a
    lightweight blocked set of labeled data to be manipulated by the DataFrame
    public API class

    This is *not* a public API class
    """
    def __init__(self, blocks, index=None, columns=None,
                 skip_integrity_check=False):
        self.index = _ensure_index(index)
        self.columns = _ensure_index(columns)
        self.blocks = blocks

        if not skip_integrity_check:
            self._verify_integrity()

    def __getstate__(self):
        return (np.asarray(self.index),
                np.asarray(self.columns),
                self.blocks)

    def __setstate__(self, state):
        index, columns, blocks = state
        self.index = _ensure_index(index)
        self.columns = _ensure_index(columns)
        self.blocks = blocks

    def __repr__(self):
        output = 'BlockManager'
        for block in self.blocks:
            output += '\n%s' % repr(block)
        return output

    def _verify_integrity(self):
        _union_block_columns(self.blocks)
        length = len(self)
        for block in self.blocks:
            assert(len(block) == length)

        tot_cols = sum(len(x.columns) for x in self.blocks)
        assert(len(self.columns) == tot_cols)

    def cast(self, dtype):
        new_blocks = []
        for block in self.blocks:
            newb = make_block(block.values.astype(dtype), block.columns)
            new_blocks.append(newb)

        new_mgr = BlockManager(new_blocks, self.index, self.columns)
        return new_mgr.consolidate()

    def is_consolidated(self):
        """
        Return True if more than one block with the same dtype
        """
        dtypes = [blk.dtype for blk in self.blocks]
        return len(dtypes) == len(set(dtypes))

    def get_slice(self, slice_obj):
        new_blocks = _slice_blocks(self.blocks, slice_obj)
        new_index = self.index[slice_obj]
        return BlockManager(new_blocks, index=new_index, columns=self.columns)

    def get_series_dict(self, index):
        return _blocks_to_series_dict(self.blocks, index)

    def __len__(self):
        # number of blocks
        return len(self.index)

    @classmethod
    def from_blocks(cls, blocks, index):
        # also checks for overlap
        columns = _union_block_columns(blocks)
        return BlockManager(blocks, index, columns)

    def __contains__(self, column):
        return column in self.columns

    @property
    def nblocks(self):
        return len(self.blocks)

    def copy(self):
        copy_blocks = [block.copy() for block in self.blocks]
        return BlockManager(copy_blocks, self.index, self.columns)

    def as_matrix(self, columns=None):
        if columns is None:
            columns = self.columns

        if len(self.blocks) == 0:
            return np.empty((len(self.index), 0), dtype=float)
        elif len(self.blocks) == 1:
            blk = self.blocks[0]
            if blk.columns.equals(columns):
                # if not, then just call interleave per below
                return blk.values
        return _interleave(self.blocks, columns)

    def xs(self, i, copy=True):
        # TODO: fix this awful mess

        if len(self.blocks) > 1:
            if not copy:
                raise Exception('cannot get view of mixed-type DataFrame')
            vals = np.concatenate([b.values[i] for b in self.blocks])
            cols = np.concatenate([b.columns for b in self.blocks])
            xs = Series(vals, index=cols).reindex(self.columns)
        else:
            vals = self.blocks[0].values[i]
            cols = self.blocks[0].columns
            xs = Series(vals, cols)
            if copy:
                xs = xs.copy()
        return xs

    def consolidate(self):
        """
        Join together blocks having same dtype

        Returns
        -------
        y : BlockManager
        """
        if self.is_consolidated():
            return self

        new_blocks = _consolidate(self.blocks)
        return BlockManager(new_blocks, self.index, self.columns)

    def get(self, col):
        _, block = self._find_block(col)
        return block.get(col)

    def delete(self, col):
        i, block = self._find_block(col)
        loc = self.columns.get_loc(col)
        self.columns = Index(np.delete(np.asarray(self.columns), loc))
        self._delete_from_block(i, col)

    def set(self, col, value):
        """
        Set new column in-place. Does not consolidate. Adds new Block if not
        contained in the current set of columns
        """
        assert(len(value) == len(self))
        if col in self.columns:
            i, block = self._find_block(col)
            if _needs_other_dtype(block, value):
                # delete from block, create and append new block
                self._delete_from_block(i, col)
                self._add_new_block(col, value)
            else:
                block.set(col, value)
        else:
            # new block
            self._add_new_block(col, value)

            # TODO: where to insert?
            self.columns = _insert_into_columns(self.columns, col,
                                                len(self.columns))

    def _delete_from_block(self, i, col):
        """
        Delete and maybe remove the whole block
        """
        block = self.blocks[i]
        assert(col in block.columns)
        if len(block.columns) == 1:
            self.blocks.pop(i)
        else:
            new_block = block.delete(col)
            self.blocks[i] = new_block

    def _add_new_block(self, col, value):
        # Do we care about dtype at the moment?
        new_block = make_block(value, [col])
        self._push_new_block(new_block)

    def _find_block(self, col):
        self._check_have(col)
        for i, block in enumerate(self.blocks):
            if col in block:
                return i, block

        raise Exception('technically unreachable code')

    def _push_new_block(self, block):
        self.blocks.append(block)

    def _check_have(self, col):
        if col not in self.columns:
            raise KeyError('no column named %s' % col)

    def _chunk_index(self, col):
        pass

    def reindex_index(self, new_index, method=None):
        new_index = _ensure_index(new_index)
        indexer, mask = self.index.get_indexer(new_index, method)

        # TODO: deal with length-0 case? or does it fall out?
        notmask = -mask
        needs_masking = len(new_index) > 0 and notmask.any()

        new_blocks = []
        for block in self.blocks:

            newb = block.reindex_index(indexer, notmask, needs_masking)
            new_blocks.append(newb)

        return BlockManager(new_blocks, new_index, self.columns)

    def merge(self, other):
        # TODO
        assert(self.index.equals(other.index))
        consolidated = _consolidate(self.blocks + other.blocks)
        cons_columns = _union_block_columns(consolidated)
        return BlockManager(consolidated, self.index, cons_columns)

    def join_on(self, other, on):
        reindexed = other.reindex_index(on)
        reindexed.index = self.index
        return self.merge(reindexed)

    def reindex_columns(self, new_columns):
        """

        """
        new_columns = _ensure_index(new_columns)
        data = self
        if not data.is_consolidated():
            data = data.consolidate()
            return data.reindex_columns(new_columns)

        new_blocks = []
        for block in self.blocks:
            newb = block.reindex_columns_from(new_columns)
            if len(newb.columns) > 0:
                new_blocks.append(newb)

        # will put these in the float bucket
        extra_columns = new_columns - self.columns
        if len(extra_columns):
            new_blocks = add_na_columns(new_blocks, extra_columns,
                                        self.index, new_columns)

        return BlockManager(new_blocks, self.index, new_columns)

    def rename_index(self, mapper):
        new_index = [mapper(x) for x in self.index]
        return BlockManager(self.blocks, new_index, self.columns)

    def rename_columns(self, mapper):
        def _rename_block(blk):
            new_cols = [mapper(x) for x in blk.columns]
            return make_block(blk.values, new_cols)

        new_blocks = [_rename_block(x) for x in self.blocks]
        new_columns = [mapper(x) for x in self.columns]
        return BlockManager(new_blocks, self.index, new_columns)

    def fillna(self, value):
        new_blocks = [b.fillna(value) for b in self.blocks]
        return BlockManager(new_blocks, self.index, self.columns)


def add_na_columns(blocks, new_columns, index, columns):
    # create new block, then consolidate
    indexer, mask = new_columns.get_indexer(columns)

    # reorder to match relative order of new columns
    new_columns = new_columns.take(indexer[mask])

    values = _nan_array(index, new_columns)
    newb = make_block(values, new_columns)
    blocks.append(newb)
    return _consolidate(blocks)

def _slice_blocks(blocks, slice_obj):
    new_blocks = []
    for block in blocks:
        newb = make_block(block.values[slice_obj], block.columns)
        new_blocks.append(newb)
    return new_blocks

# TODO!
def _needs_other_dtype(block, to_insert):
    if block.dtype == np.float64:
        return not issubclass(to_insert.dtype.type, (np.integer, np.floating))
    elif block.dtype == np.object_:
        return issubclass(to_insert.dtype.type, (np.integer, np.floating))
    else:
        raise Exception('have not handled this case yet')

def _blocks_to_series_dict(blocks, index=None):
    series_dict = {}

    if index is None:
        index = Index(np.arange(len(blocks[0])))

    for block in blocks:
        for col, vec in zip(block.columns, block.values.T):
            series_dict[col] = Series(vec, index=index)
    return series_dict

def _insert_column(blocks, column, value):
    """
    Default: new block
    """
    pass

def _set_column(blocks, column, value):
    pass

def _delete_column(blocks, columns):
    pass

def _reindex_blocks(blocks, old_columns, new_columns):
    pass

def _interleave(blocks, columns):
    """
    Return ndarray from blocks with specified column order
    """
    dtype = _interleaved_dtype(blocks)
    columns = _ensure_index(columns)

    result = np.empty((len(blocks[0]), len(columns)), dtype=dtype)
    result.fill(nan)

    for block in blocks:
        indexer, mask = columns.get_indexer(block.columns)

        if mask.all():
            result[:, indexer] = block.values
        else:
            result[:, indexer[mask]] = block.values[:, mask]

    return result

def _interleaved_dtype(blocks):
    for block in blocks:
        if not issubclass(block.dtype.type, np.floating):
            return object
    return np.float64

def _consolidate(blocks):
    """
    Merge blocks having same dtype
    """
    get_dtype = lambda x: x.dtype

    # sort by dtype
    grouper = itertools.groupby(sorted(blocks, key=get_dtype),
                                lambda x: x.dtype)

    new_blocks = []
    for dtype, group_blocks in grouper:
        new_block = _merge_blocks(list(group_blocks))
        new_blocks.append(new_block)

    return new_blocks

def _merge_blocks(blocks):
    new_values = np.hstack([b.values for b in blocks])
    new_columns = np.concatenate([b.columns for b in blocks])
    new_block = make_block(new_values, new_columns)
    return new_block

def _xs(blocks, i, copy=True):
    if copy:
        return np.concatenate([b[i] for b in blocks])
    else:
        if len(blocks) == 1:
            return blocks[0].values[i]
        else:
            raise Exception('cannot get view with mixed-type data')

def _union_block_columns(blocks):
    seen = Index([])

    for block in blocks:
        len_before = len(seen)
        seen = seen.union(block.columns)
        if len(seen) != len_before + len(block.columns):
            raise Exception('column names overlap')

    return seen

def _nan_manager_matching(index, columns):
    # what if one of these is empty?
    values = _nan_array(index, columns)
    block = make_block(values, columns)
    return BlockManager([block], columns)

def _nan_array(index, columns, dtype=np.float64):
    if index is None:
        index = NULL_INDEX
    if columns is None:
        columns = NULL_INDEX

    values = np.empty((len(index), len(columns)), dtype=dtype)
    values.fill(nan)
    return values

import unittest
class TestBlockOperations(unittest.TestCase):

    def test_interleave(self):
        pass

    def test_consolidate(self):
        pass

    def test_xs(self):
        pass

if __name__ == '__main__':
    n = 10
    floats = np.repeat(np.atleast_2d(np.arange(3.)), n, axis=0)
    objects = np.empty((n, 2), dtype=object)
    objects[:, 0] = 'foo'
    objects[:, 1] = 'bar'

    float_cols = Index(['a', 'c', 'e'])
    object_cols = Index(['b', 'd'])
    columns = Index(sorted(float_cols + object_cols))
    index = np.arange(n)
    new_columns = Index(['a', 'c', 'e', 'b', 'd'])

    fblock = make_block(floats, float_cols)
    oblock = make_block(objects, object_cols)

    blocks = [fblock, oblock]

    interleaved = _interleave(blocks, columns)

    mgr = BlockManager(blocks, index, columns)
