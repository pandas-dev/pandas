from operator import attrgetter
import itertools

from numpy import nan
import numpy as np

from pandas.core.index import Index, NULL_INDEX
from pandas.core.common import _ensure_index, _try_sort
from pandas.core.series import Series
import pandas.core.common as common

class Block(object):
    """
    Canonical unit of homogeneous dtype contained in DataFrame

    Index-ignorant; let the container take care of that
    """
    def __init__(self, values, columns, ref_columns):
        values = _convert_if_1d(values)
        if issubclass(values.dtype.type, basestring):
            values = np.array(values, dtype=object)

        self.values = values
        assert(len(columns) == values.shape[1])
        self.columns = _ensure_index(columns)
        self.ref_columns = _ensure_index(ref_columns)

    _ref_locs = None
    @property
    def ref_locs(self):
        if self._ref_locs is None:
            indexer, mask = self.ref_columns.get_indexer(self.columns)
            assert(mask.all())
            self._ref_locs = indexer
        return self._ref_locs

    def set_ref_columns(self, ref_columns, maybe_rename=True):
        """
        If maybe_rename=True, need to set the columns for this guy
        """
        assert(isinstance(ref_columns, Index))
        if maybe_rename:
            self.columns = ref_columns.take(self.ref_locs)
        self.ref_columns = ref_columns

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
        # should not pickle generally (want to share ref_columns), but here for
        # completeness
        return (np.asarray(self.columns),
                np.asarray(self.ref_columns),
                self.values)

    def __setstate__(self, state):
        columns, ref_columns, values = state
        self.columns = Index(columns)
        self.ref_columns = Index(ref_columns)
        self.values = values

    @property
    def shape(self):
        return self.values.shape

    @property
    def dtype(self):
        return self.values.dtype

    def copy(self):
        return make_block(self.values.copy(), self.columns,
                          self.ref_columns)

    def merge(self, other):
        assert(self.ref_columns.equals(other.ref_columns))

        # Not sure whether to allow this or not
        # if not union_ref.equals(other.ref_columns):
        #     union_ref = self.ref_columns + other.ref_columns
        return _merge_blocks([self, other], self.ref_columns)

    def reindex_index(self, indexer, notmask, needs_masking):
        """
        Reindex using pre-computed indexer information
        """
        new_values = self.values.take(indexer, axis=0)
        if needs_masking:
            new_values = _cast_if_bool_int(new_values)
            common.null_out_axis(new_values, notmask, 0)
        return make_block(new_values, self.columns, self.ref_columns)

    def reindex_columns_from(self, new_columns):
        """
        Reindex to only those columns contained in the input set of columns

        E.g. if you have ['a', 'b'], and the input columns is ['b', 'c', 'd'],
        then the resulting columns will be ['b']

        Returns
        -------
        reindexed : Block
        """
        indexer, mask = self.columns.get_indexer(new_columns)
        masked_idx = indexer[mask]
        new_values = self.values.take(masked_idx, axis=1)
        new_cols = self.columns.take(masked_idx)
        return make_block(new_values, new_cols, new_columns)

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
        new_cols = np.delete(np.asarray(self.columns), loc)
        new_values = np.delete(self.values, loc, 1)
        return make_block(new_values, new_cols, self.ref_columns)

    def fillna(self, value):
        new_values = self.values.copy()
        mask = common.isnull(new_values.ravel())
        new_values.flat[mask] = value
        return make_block(new_values, self.columns, self.ref_columns)

def _insert_into_columns(columns, col, loc):
    columns = np.asarray(columns)
    new_columns = np.insert(columns, loc, col)
    return Index(new_columns)

def _cast_if_bool_int(values):
    if issubclass(values.dtype.type, np.int_):
        values = values.astype(float)
    elif issubclass(values.dtype.type, np.bool_):
        values = values.astype(object)
    return values

def _convert_if_1d(values):
    if values.ndim == 1:
        values = np.atleast_2d(values).T

    return values

#-------------------------------------------------------------------------------
# Is this even possible?

class FloatBlock(Block):

    def can_store(self, value):
        return issubclass(value.dtype.type, (np.integer, np.floating))

class IntBlock(Block):

    def can_store(self, value):
        return issubclass(value.dtype.type, np.integer)

class BoolBlock(Block):

    def can_store(self, value):
        return issubclass(value.dtype.type, np.bool_)

class ObjectBlock(Block):

    def can_store(self, value):
        return not issubclass(value.dtype.type,
                              (np.integer, np.floating, np.bool_))

def make_block(values, columns, ref_columns):
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

    return klass(values, columns, ref_columns)

# TODO: flexible with index=None and/or columns=None

class BlockManager(object):
    """
    Core internal data structure to implement DataFrame

    Manage a bunch of labeled 2D mixed-type ndarrays. Essentially it's a
    lightweight blocked set of labeled data to be manipulated by the DataFrame
    public API class

    Parameters
    ----------


    Notes
    -----
    This is *not* a public API class
    """
    def __init__(self, blocks, index=None, columns=None,
                 skip_integrity_check=False):
        self.index = _ensure_index(index)
        self._columns = _ensure_index(columns)
        self.blocks = blocks

        if not skip_integrity_check:
            self._verify_integrity()

    _columns = None
    def _set_columns(self, value):
        self._columns = _ensure_index(value)

        for block in self.blocks:
            block.set_ref_columns(self._columns, maybe_rename=True)

    columns = property(fget=attrgetter('_columns'),
                       fset=_set_columns)

    def set_columns_norename(self, value):
        self._columns = _ensure_index(value)

        for block in self.blocks:
            block.set_ref_columns(self._columns, maybe_rename=False)

    def __getstate__(self):
        block_values = [b.values for b in self.blocks]
        block_columns = [np.asarray(b.columns) for b in self.blocks]
        return (np.asarray(self.index),
                np.asarray(self.columns),
                block_values, block_columns)

    def __setstate__(self, state):
        index, ref_columns, bvalues, bcolumns = state
        self.index = _ensure_index(index)
        self._columns = _ensure_index(ref_columns)

        blocks = []
        for values, columns in zip(bvalues, bcolumns):
            blk = make_block(values, columns, self.columns)
            blocks.append(blk)
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
            newb = make_block(block.values.astype(dtype), block.columns,
                              block.ref_columns)
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
        # TODO: fix this mess

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

        new_blocks = _consolidate(self.blocks, self.columns)
        return BlockManager(new_blocks, self.index, self.columns)

    def get(self, col):
        _, block = self._find_block(col)
        return block.get(col)

    def delete(self, col):
        i, _ = self._find_block(col)
        loc = self.columns.get_loc(col)
        new_cols = Index(np.delete(np.asarray(self.columns), loc))

        self._delete_from_block(i, col)
        self.set_columns_norename(new_cols)

    def set(self, col, value):
        """
        Set new column in-place. Does not consolidate. Adds new Block if not
        contained in the current set of columns
        """
        assert(len(value) == len(self))
        if col in self.columns:
            i, block = self._find_block(col)
            if not block.can_store(value):
                # delete from block, create and append new block
                self._delete_from_block(i, col)
                self._add_new_block(col, value)
            else:
                block.set(col, value)
        else:
            # TODO: where to insert?
            new_cols = _insert_into_columns(self.columns, col,
                                            len(self.columns))
            self.set_columns_norename(new_cols)
            # new block
            self._add_new_block(col, value)

    def _delete_from_block(self, i, col):
        """
        Delete and maybe remove the whole block
        """
        block = self.blocks[i]
        newb = block.delete(col)

        if len(newb.ref_locs) == 0:
            self.blocks.pop(i)
        else:
            self.blocks[i] = newb

    def _add_new_block(self, col, value):
        # Do we care about dtype at the moment?
        new_block = make_block(value, [col], self.columns)
        self.blocks.append(new_block)

    def _find_block(self, col):
        self._check_have(col)
        for i, block in enumerate(self.blocks):
            if col in block:
                return i, block

    def _check_have(self, col):
        if col not in self.columns:
            raise KeyError('no column named %s' % col)

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

        intersection = self.columns.intersection(other.columns)
        try:
            assert(len(intersection) == 0)
        except AssertionError:
            raise Exception('columns overlap: %s' % intersection)

        cons_cols = self.columns + other.columns
        consolidated = _consolidate(self.blocks + other.blocks, cons_cols)
        return BlockManager(consolidated, self.index, cons_cols)

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

        # TODO: this part could be faster (!)
        _, mask = self.columns.get_indexer(new_columns)
        notmask = -mask

        if notmask.any():
            extra_columns = new_columns[notmask]
            na_block = add_na_columns(extra_columns, self.index, new_columns)
            new_blocks.append(na_block)
            new_blocks = _consolidate(new_blocks, new_columns)

        return BlockManager(new_blocks, self.index, new_columns)

    def rename_index(self, mapper):
        new_index = [mapper(x) for x in self.index]
        return BlockManager(self.blocks, new_index, self.columns)

    def rename_columns(self, mapper):
        new_columns = Index([mapper(x) for x in self.columns])
        new_blocks = []
        for block in self.blocks:
            newb = block.copy()
            newb.set_ref_columns(new_columns, maybe_rename=True)
            new_blocks.append(newb)
        return BlockManager(new_blocks, self.index, new_columns)

    def fillna(self, value):
        """

        """
        new_blocks = [b.fillna(value) for b in self.blocks]
        return BlockManager(new_blocks, self.index, self.columns)

    @property
    def block_id_vector(self):
        # TODO
        result = np.empty(len(self.columns), dtype=int)
        result.fill(-1)

        for i, blk in enumerate(self.blocks):
            indexer, mask = self.columns.get_indexer(blk.columns)
            assert(mask.all())
            result.put(indexer, i)

        assert((result >= 0).all())
        return result

_data_types = [np.float_, np.int_]
def form_blocks(data, index, columns):
    from pandas.core.internals import add_na_columns

    # pre-filter out columns if we passed it
    if columns is None:
        columns = Index(_try_sort(data.keys()))
        extra_columns = NULL_INDEX
    else:
        columns = _ensure_index(columns)
        extra_columns = columns - Index(data.keys())

    # put "leftover" columns in float bucket, where else?
    # generalize?
    num_dict = {}
    object_dict = {}
    for k, v in data.iteritems():
        if issubclass(v.dtype.type, (np.floating, np.integer)):
            num_dict[k] = v
        else:
            object_dict[k] = v

    blocks = []

    if len(num_dict) > 0:
        num_dtypes = set(v.dtype for v in num_dict.values())
        if len(num_dtypes) > 1:
            num_dtype = np.float_
        else:
            num_dtype = list(num_dtypes)[0]

        # TODO: find corner cases
        # TODO: check type inference
        num_block = _simple_blockify(num_dict, columns, num_dtype)
        blocks.append(num_block)

    if len(object_dict) > 0:
        object_block = _simple_blockify(object_dict, columns, np.object_)
        blocks.append(object_block)

    if len(extra_columns):
        na_block = add_na_columns(extra_columns, index, columns)
        blocks.append(na_block)
        blocks = _consolidate(blocks, columns)

    return blocks, columns

def _simple_blockify(dct, ref_columns, dtype):
    block_columns, values = _stack_dict(dct)
    # CHECK DTYPE?
    if values.dtype != dtype:
        values = values.astype(dtype)

    return make_block(values, block_columns, ref_columns)

def _stack_dict(dct):
    columns = Index(_try_sort(dct))
    stacked = np.vstack([dct[k].values for k in columns]).T
    return columns, stacked

# def _float_blockify(dct, index, columns):
#     n = len(index)
#     k = len(columns)
#     values = np.empty((n, k), dtype=np.float64)
#     values.fill(nan)

#     if len(dct) > 0:
#         dict_columns, stacked = _stack_dict(dct)
#         indexer, mask = columns.get_indexer(dict_columns)
#         assert(mask.all())
#         values[:, indexer] = stacked

#     # do something with dtype?
#     return make_block(values, columns)

def add_na_columns(new_columns, index, ref_columns):
    # create new block, then consolidate
    values = _nan_array(index, new_columns)
    return make_block(values, new_columns, ref_columns)

def _slice_blocks(blocks, slice_obj):
    new_blocks = []
    for block in blocks:
        newb = make_block(block.values[slice_obj], block.columns,
                          block.ref_columns)
        new_blocks.append(newb)
    return new_blocks

def _blocks_to_series_dict(blocks, index=None):
    series_dict = {}

    # if index is None:
    #     index = Index(np.arange(len(blocks[0])))

    for block in blocks:
        for col, vec in zip(block.columns, block.values.T):
            series_dict[col] = Series(vec, index=index)
    return series_dict

def _interleave(blocks, columns):
    """
    Return ndarray from blocks with specified column order
    """
    dtype = _interleaved_dtype(blocks)
    columns = _ensure_index(columns)

    result = np.empty((len(blocks[0]), len(columns)), dtype=dtype)
    colmask = np.zeros(len(columns), dtype=bool)

    for block in blocks:
        indexer, mask = columns.get_indexer(block.columns)

        if mask.all():
            result[:, indexer] = block.values
        else:
            indexer = indexer[mask]
            result[:, indexer] = block.values[:, mask]

        colmask[indexer] = 1

    # By construction, all of the column should be covered by one of the blocks
    assert(colmask.all())
    return result

def _interleaved_dtype(blocks):
    have_int = False
    have_bool = False
    have_object = False
    have_float = False

    for block in blocks:
        if isinstance(block, FloatBlock):
            have_float = True
        elif isinstance(block, IntBlock):
            have_int = True
        elif isinstance(block, BoolBlock):
            have_bool = True
        elif isinstance(block, ObjectBlock):
            have_object = True
        else: # pragma: no cover
            raise Exception('Unrecognized block type')

    have_numeric = have_float or have_int

    if have_object:
        return np.object_
    elif have_bool and have_numeric:
        return np.object_
    elif have_bool:
        return np.bool_
    elif have_int and not have_float:
        return np.int_
    else:
        return np.float64

def _consolidate(blocks, columns):
    """
    Merge blocks having same dtype
    """
    get_dtype = lambda x: x.dtype

    # sort by dtype
    grouper = itertools.groupby(sorted(blocks, key=get_dtype),
                                lambda x: x.dtype)

    new_blocks = []
    for dtype, group_blocks in grouper:
        new_block = _merge_blocks(list(group_blocks), columns)
        new_blocks.append(new_block)

    return new_blocks

def _merge_blocks(blocks, columns):
    new_values = np.hstack([b.values for b in blocks])
    new_cols = np.concatenate([b.columns for b in blocks])
    new_block = make_block(new_values, new_cols, columns)
    return new_block.reindex_columns_from(columns)

def _union_block_columns(blocks):
    seen = None
    for block in blocks:
        len_before = 0 if seen is None else len(seen)
        if seen is None:
            seen = block.columns
        else:
            seen = seen.union(block.columns)
        if len(seen) != len_before + len(block.columns):
            raise Exception('column names overlap')

    return seen

def _nan_array(index, columns, dtype=np.float64):
    if index is None:
        index = NULL_INDEX
    if columns is None:
        columns = NULL_INDEX

    values = np.empty((len(index), len(columns)), dtype=dtype)
    values.fill(nan)
    return values

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

    float_locs = new_columns.get_indexer(float_cols)[0]
    obj_locs = new_columns.get_indexer(object_cols)[0]

    fblock = make_block(floats, float_cols, float_cols)
    oblock = make_block(objects, object_cols, object_cols)

    # blocks = [fblock, oblock]

    # interleaved = _interleave(blocks, columns)

    # mgr = BlockManager(blocks, index, columns)
