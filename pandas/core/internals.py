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
    Canonical n-dimensional unit of homogeneous dtype contained in a pandas data
    structure

    Index-ignorant; let the container take care of that
    """
    def __init__(self, values, items, ref_items, ndim=2):
        # values = _convert_if_1d(values)

        if issubclass(values.dtype.type, basestring):
            values = np.array(values, dtype=object)

        assert(values.ndim == ndim)
        assert(len(items) == len(values))

        self.values = values
        self.ndim = 2
        self.items = _ensure_index(items)
        self.ref_items = _ensure_index(ref_items)

    _ref_locs = None
    @property
    def ref_locs(self):
        if self._ref_locs is None:
            indexer, mask = self.ref_items.get_indexer(self.items)
            assert(mask.all())
            self._ref_locs = indexer
        return self._ref_locs

    def set_ref_items(self, ref_items, maybe_rename=True):
        """
        If maybe_rename=True, need to set the items for this guy
        """
        assert(isinstance(ref_items, Index))
        if maybe_rename:
            self.items = ref_items.take(self.ref_locs)
        self.ref_items = ref_items

    def __repr__(self):
        shape = ' x '.join([str(s) for s in self.shape])
        name = type(self).__name__
        return '%s: %s, %s, dtype %s' % (name, self.items, shape, self.dtype)

    def __contains__(self, item):
        return item in self.items

    def __len__(self):
        return len(self.values)

    def __getstate__(self):
        # should not pickle generally (want to share ref_items), but here for
        # completeness
        return (np.asarray(self.items), np.asarray(self.ref_items),
                self.values)

    def __setstate__(self, state):
        items, ref_items, values = state
        self.items = Index(items)
        self.ref_items = Index(ref_items)
        self.values = values

    @property
    def shape(self):
        return self.values.shape

    @property
    def dtype(self):
        return self.values.dtype

    def copy(self):
        return make_block(self.values.copy(), self.items,
                          self.ref_items, ndim=self.ndim)

    def merge(self, other):
        assert(self.ref_items.equals(other.ref_items))

        # Not sure whether to allow this or not
        # if not union_ref.equals(other.ref_items):
        #     union_ref = self.ref_items + other.ref_items
        return _merge_blocks([self, other], self.ref_items)

    def reindex_axis(self, indexer, notmask, needs_masking, axis=0):
        """
        Reindex using pre-computed indexer information
        """
        new_values = self.values.take(indexer, axis=axis)
        if needs_masking:
            new_values = _cast_if_bool_int(new_values)
            common.null_out_axis(new_values, notmask, axis)
        return make_block(new_values, self.items, self.ref_items,
                          ndim=self.ndim)

    def reindex_items_from(self, new_ref_items):
        """
        Reindex to only those items contained in the input set of items

        E.g. if you have ['a', 'b'], and the input items is ['b', 'c', 'd'],
        then the resulting items will be ['b']

        Returns
        -------
        reindexed : Block
        """
        indexer, mask = self.items.get_indexer(new_ref_items)
        masked_idx = indexer[mask]
        new_values = self.values.take(masked_idx, axis=0)
        new_items = self.items.take(masked_idx)
        return make_block(new_values, new_items, new_ref_items, ndim=self.ndim)

    def get(self, item):
        loc = self.items.get_loc(item)
        return self.values[loc]

    def set(self, item, value):
        """
        Modify Block in-place with new item value

        Returns
        -------
        None
        """
        loc = self.items.get_loc(item)
        self.values[loc] = value

    def delete(self, item):
        """
        Returns
        -------
        y : Block (new object)
        """
        loc = self.items.get_loc(item)
        new_items = np.delete(np.asarray(self.items), loc)
        new_values = np.delete(self.values, loc, 0)
        return make_block(new_values, new_items, self.ref_items, ndim=self.ndim)

    def fillna(self, value):
        new_values = self.values.copy()
        mask = common.isnull(new_values.ravel())
        new_values.flat[mask] = value
        return make_block(new_values, self.items, self.ref_items,
                          ndim=self.ndim)

def _insert_into_items(items, item, loc):
    items = np.asarray(items)
    new_items = np.insert(items, loc, item)
    return Index(new_items)

def _cast_if_bool_int(values):
    if issubclass(values.dtype.type, np.int_):
        values = values.astype(float)
    elif issubclass(values.dtype.type, np.bool_):
        values = values.astype(object)
    return values

def _convert_if_1d(values):
    if values.ndim == 1:
        values = np.atleast_2d(values)

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

def make_block(values, items, ref_items, ndim=2):
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

    return klass(values, items, ref_items, ndim=ndim)

# TODO: flexible with index=None and/or items=None

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
    def __init__(self, blocks, axes, skip_integrity_check=False, ndim=2):
        self.axes = [_ensure_index(ax) for ax in axes]
        self.blocks = blocks

        if not skip_integrity_check:
            self._verify_integrity()

    def is_mixed_dtype(self):
        counts = set()
        for block in self.blocks:
            counts.add(block.dtype)
            if len(counts) > 1:
                return True
        return False

    def set_axis(self, axis, value):
        cur_axis = self.axes[axis]
        if len(value) != len(cur_axis):
            raise Exception('Length mismatch (%d vs %d)'
                            % (len(index), len(cur_axis)))
        self.axes[axis] = _ensure_index(value)

    def _set_items(self, value):
        self.set_axis(0, value)
        for block in self.blocks:
            block.set_ref_items(self.items, maybe_rename=True)

    def _get_items(self):
        return self.axes[0]

    items = property(fget=_get_items, fset=_set_items)

    def set_items_norename(self, value):
        value = _ensure_index(value)
        self.axes[0] = value

        for block in self.blocks:
            block.set_ref_items(value, maybe_rename=False)

    def __getstate__(self):
        block_values = [b.values for b in self.blocks]
        block_items = [np.asarray(b.items) for b in self.blocks]
        axes_array = [np.asarray(ax) for ax in self.axes]
        return axes_array, block_values, block_items, self.ndim

    def __setstate__(self, state):
        ax_arrays, bvalues, bitems, ndim = state

        self.axes = [_ensure_index(ax) for ax in ax_arrays]

        blocks = []
        for values, items in zip(bvalues, bitems):
            blk = make_block(values, items, self.axes[0])
            blocks.append(blk)
        self.blocks = blocks

    def __repr__(self):
        output = 'BlockManager'
        for i, ax in enumerate(self.axes):
            if i == 0:
                output += 'Items: \n%s' % ax
            else:
                output += 'Axis %d: \n%s' % (i, ax)

        for block in self.blocks:
            output += '\n%s' % repr(block)
        return output

    @property
    def shape(self):
        return tuple(len(ax) for ax in self.axes)

    def _verify_integrity(self):
        _union_block_items(self.blocks)

        mgr_shape = self.shape
        for block in self.blocks:
            assert(block.values.shape[1:] == mgr_shape[1:])
        tot_items = sum(len(x.items) for x in self.blocks)
        assert(len(self.items) == tot_items)

    def cast(self, dtype):
        new_blocks = []
        for block in self.blocks:
            newb = make_block(block.values.astype(dtype), block.items,
                              block.ref_items, ndim=self.ndim)
            new_blocks.append(newb)

        new_mgr = BlockManager(new_blocks, [self.index, self.items],
                               ndim=self.ndim)
        return new_mgr.consolidate()

    def is_consolidated(self):
        """
        Return True if more than one block with the same dtype
        """
        dtypes = [blk.dtype for blk in self.blocks]
        return len(dtypes) == len(set(dtypes))

    def get_slice(self, slice_obj, axis=0):
        new_blocks = _slice_blocks(self.blocks, slice_obj, axis)

        new_axes = list(self.axes)
        new_axes[axis] = new_axes[axis][slice_obj]
        return BlockManager(new_blocks, new_axes, ndim=self.ndim)

    def get_series_dict(self, index):
        return _blocks_to_series_dict(self.blocks, index)

    def __len__(self):
        # number of blocks
        return len(self.index)

    @classmethod
    def from_blocks(cls, blocks, index):
        # also checks for overlap
        items = _union_block_items(blocks)
        return BlockManager(blocks, index, items)

    def __contains__(self, item):
        return item in self.items

    @property
    def nblocks(self):
        return len(self.blocks)

    def copy(self):
        copy_blocks = [block.copy() for block in self.blocks]
        return BlockManager(copy_blocks, self.index, self.items)

    def as_matrix(self, items=None):
        if len(self.blocks) == 0:
            mat = np.empty((len(self.index), 0), dtype=float)
        elif len(self.blocks) == 1:
            blk = self.blocks[0]
            if items is None or blk.items.equals(items):
                # if not, then just call interleave per below
                mat = blk.values
        else:
            if items is None:
                mat = _interleave(self.blocks, self.items)
            else:
                mat = self.reindex_items(items).as_matrix()

        return mat

    def xs(self, i, copy=True):
        # TODO: fix this mess

        if len(self.blocks) > 1:
            if not copy:
                raise Exception('cannot get view of mixed-type or '
                                'non-consolidated DataFrame')
            vals = np.concatenate([b.values[i] for b in self.blocks])
            items = np.concatenate([b.items for b in self.blocks])
            xs = Series(vals, index=items).reindex(self.items)
        else:
            vals = self.blocks[0].values[i]
            items = self.blocks[0].items
            xs = Series(vals, items)
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

        new_blocks = _consolidate(self.blocks, self.items)
        return BlockManager(new_blocks, self.index, self.items)

    def get(self, item):
        _, block = self._find_block(item)
        return block.get(item)

    def delete(self, item):
        i, _ = self._find_block(item)
        loc = self.items.get_loc(item)
        new_items = Index(np.delete(np.asarray(self.items), loc))

        self._delete_from_block(i, item)
        self.set_items_norename(new_items)

    def set(self, item, value):
        """
        Set new item in-place. Does not consolidate. Adds new Block if not
        contained in the current set of items
        """
        assert(len(value) == len(self))
        if item in self.items:
            i, block = self._find_block(item)
            if not block.can_store(value):
                # delete from block, create and append new block
                self._delete_from_block(i, item)
                self._add_new_block(item, value)
            else:
                block.set(item, value)
        else:
            # TODO: where to insert?
            new_items = _insert_into_items(self.items, item,
                                            len(self.items))
            self.set_items_norename(new_items)
            # new block
            self._add_new_block(item, value)

    def _delete_from_block(self, i, item):
        """
        Delete and maybe remove the whole block
        """
        block = self.blocks[i]
        newb = block.delete(item)

        if len(newb.ref_locs) == 0:
            self.blocks.pop(i)
        else:
            self.blocks[i] = newb

    def _add_new_block(self, item, value):
        # Do we care about dtype at the moment?
        new_block = make_block(value, [item], self.items)
        self.blocks.append(new_block)

    def _find_block(self, item):
        self._check_have(item)
        for i, block in enumerate(self.blocks):
            if item in block:
                return i, block

    def _check_have(self, item):
        if item not in self.items:
            raise KeyError('no item named %s' % item)

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

        return BlockManager(new_blocks, new_index, self.items)

    def merge(self, other):
        # TODO
        assert(self.index.equals(other.index))

        intersection = self.items.intersection(other.items)
        try:
            assert(len(intersection) == 0)
        except AssertionError:
            raise Exception('items overlap: %s' % intersection)

        cons_items = self.items + other.items
        consolidated = _consolidate(self.blocks + other.blocks, cons_items)
        return BlockManager(consolidated, self.index, cons_items)

    def join_on(self, other, on):
        reindexed = other.reindex_index(on)
        reindexed.index = self.index
        return self.merge(reindexed)

    def reindex_items(self, new_items):
        """

        """
        new_items = _ensure_index(new_items)
        data = self
        if not data.is_consolidated():
            data = data.consolidate()
            return data.reindex_items(new_items)

        new_blocks = []
        for block in self.blocks:
            newb = block.reindex_items_from(new_items)
            if len(newb.items) > 0:
                new_blocks.append(newb)

        # TODO: this part could be faster (!)
        _, mask = self.items.get_indexer(new_items)
        notmask = -mask

        if notmask.any():
            extra_items = new_items[notmask]
            na_block = add_na_items(extra_items, self.index, new_items)
            new_blocks.append(na_block)
            new_blocks = _consolidate(new_blocks, new_items)

        return BlockManager(new_blocks, self.index, new_items)

    def rename_index(self, mapper):
        new_index = [mapper(x) for x in self.index]
        return BlockManager(self.blocks, new_index, self.items)

    def rename_items(self, mapper):
        new_items = Index([mapper(x) for x in self.items])
        new_blocks = []
        for block in self.blocks:
            newb = block.copy()
            newb.set_ref_items(new_items, maybe_rename=True)
            new_blocks.append(newb)
        return BlockManager(new_blocks, self.index, new_items)

    def fillna(self, value):
        """

        """
        new_blocks = [b.fillna(value) for b in self.blocks]
        return BlockManager(new_blocks, self.index, self.items)

    @property
    def block_id_vector(self):
        # TODO
        result = np.empty(len(self.items), dtype=int)
        result.fill(-1)

        for i, blk in enumerate(self.blocks):
            indexer, mask = self.items.get_indexer(blk.items)
            assert(mask.all())
            result.put(indexer, i)

        assert((result >= 0).all())
        return result

_data_types = [np.float_, np.int_]
def form_blocks(data, index, items):
    from pandas.core.internals import add_na_items

    # pre-filter out items if we passed it
    if items is None:
        items = Index(_try_sort(data.keys()))
        extra_items = NULL_INDEX
    else:
        items = _ensure_index(items)
        extra_items = items - Index(data.keys())

    # put "leftover" items in float bucket, where else?
    # generalize?
    num_dict = {}
    bool_dict = {}
    object_dict = {}
    for k, v in data.iteritems():
        if issubclass(v.dtype.type, (np.floating, np.integer)):
            num_dict[k] = v
        elif v.dtype == np.bool_:
            bool_dict[k] = v
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
        num_block = _simple_blockify(num_dict, items, num_dtype)
        blocks.append(num_block)

    if len(bool_dict):
        bool_block = _simple_blockify(bool_dict, items, np.bool_)
        blocks.append(bool_block)

    if len(object_dict) > 0:
        object_block = _simple_blockify(object_dict, items, np.object_)
        blocks.append(object_block)

    if len(extra_items):
        na_block = add_na_items(extra_items, index, items)
        blocks.append(na_block)
        blocks = _consolidate(blocks, items)

    return blocks, items

def _simple_blockify(dct, ref_items, dtype):
    block_items, values = _stack_dict(dct)
    # CHECK DTYPE?
    if values.dtype != dtype:
        values = values.astype(dtype)

    return make_block(values, block_items, ref_items)

def _stack_dict(dct):
    items = Index(_try_sort(dct))
    stacked = np.vstack([dct[k].values for k in items]).T
    return items, stacked

def add_na_items(new_items, index, ref_items):
    # create new block, then consolidate
    values = _nan_array(index, new_items)
    return make_block(values, new_items, ref_items)

def _slice_blocks(blocks, slice_obj, axis):
    new_blocks = []
    for block in blocks:
        newb = make_block(block.values[slice_obj], block.items,
                          block.ref_items)
        new_blocks.append(newb)
    return new_blocks

def _blocks_to_series_dict(blocks, index=None):
    series_dict = {}

    for block in blocks:
        for item, vec in zip(block.items, block.values.T):
            series_dict[item] = Series(vec, index=index)
    return series_dict

def _interleave(blocks, items):
    """
    Return ndarray from blocks with specified item order
    Items must be contained in the blocks
    """
    dtype = _interleaved_dtype(blocks)
    items = _ensure_index(items)

    result = np.empty((len(blocks[0]), len(items)), dtype=dtype)
    itemmask = np.zeros(len(items), dtype=bool)

    # By construction, all of the item should be covered by one of the blocks
    for block in blocks:
        indexer, mask = items.get_indexer(block.items)
        assert(mask.all())
        result[:, indexer] = block.values
        itemmask[indexer] = 1
    assert(itemmask.all())
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

def _consolidate(blocks, items):
    """
    Merge blocks having same dtype
    """
    get_dtype = lambda x: x.dtype

    # sort by dtype
    grouper = itertools.groupby(sorted(blocks, key=get_dtype),
                                lambda x: x.dtype)

    new_blocks = []
    for dtype, group_blocks in grouper:
        new_block = _merge_blocks(list(group_blocks), items)
        new_blocks.append(new_block)

    return new_blocks

def _merge_blocks(blocks, items):
    new_values = np.hstack([b.values for b in blocks])
    new_items = np.concatenate([b.items for b in blocks])
    new_block = make_block(new_values, new_items, items)
    return new_block.reindex_items_from(items)

def _merge_blocks2(blocks, items):
    new_values = np.vstack([b.values for b in blocks])
    new_items = np.concatenate([b.items for b in blocks])
    new_block = make_block(new_values, new_items, items)
    return new_block.reindex_items_from(items)

def _union_block_items(blocks):
    seen = None
    for block in blocks:
        len_before = 0 if seen is None else len(seen)
        if seen is None:
            seen = block.items
        else:
            seen = seen.union(block.items)
        if len(seen) != len_before + len(block.items):
            raise Exception('item names overlap')

    return seen

def _nan_array(index, items, dtype=np.float64):
    if index is None:
        index = NULL_INDEX
    if items is None:
        items = NULL_INDEX

    values = np.empty((len(index), len(items)), dtype=dtype)
    values.fill(nan)
    return values

if __name__ == '__main__':
    n = 10
    floats = np.repeat(np.atleast_2d(np.arange(3.)), n, axis=0)
    objects = np.empty((n, 2), dtype=object)
    objects[:, 0] = 'foo'
    objects[:, 1] = 'bar'

    float_items = Index(['a', 'c', 'e'])
    object_items = Index(['b', 'd'])
    items = Index(sorted(float_items + object_items))
    index = np.arange(n)
    new_items = Index(['a', 'c', 'e', 'b', 'd'])

    float_locs = new_items.get_indexer(float_items)[0]
    obj_locs = new_items.get_indexer(object_items)[0]

    fblock = make_block(floats, float_items, float_items)
    oblock = make_block(objects, object_items, object_items)

    # blocks = [fblock, oblock]

    # interleaved = _interleave(blocks, items)

    # mgr = BlockManager(blocks, index, items)
