import itertools

from numpy import nan
import numpy as np

from pandas.core.index import Index, _ensure_index
import pandas.core.common as common
import pandas._tseries as lib

class Block(object):
    """
    Canonical n-dimensional unit of homogeneous dtype contained in a pandas data
    structure

    Index-ignorant; let the container take care of that
    """
    __slots__ = ['items', 'ref_items', '_ref_locs', 'values', 'ndim']

    def __init__(self, values, items, ref_items, ndim=2,
                 do_integrity_check=False):
        if issubclass(values.dtype.type, basestring):
            values = np.array(values, dtype=object)

        assert(values.ndim == ndim)
        assert(len(items) == len(values))

        self.values = values
        self.ndim = ndim
        self.items = _ensure_index(items)
        self.ref_items = _ensure_index(ref_items)

        if do_integrity_check:
            self._check_integrity()

    def _check_integrity(self):
        if len(self.items) < 2:
            return
        # monotonicity
        return (self.ref_locs[1:] > self.ref_locs[:-1]).all()

    _ref_locs = None
    @property
    def ref_locs(self):
        if self._ref_locs is None:
            indexer = self.ref_items.get_indexer(self.items)
            assert((indexer != -1).all())
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
        return (self.items, self.ref_items, self.values)

    def __setstate__(self, state):
        items, ref_items, values = state
        self.items = _ensure_index(items)
        self.ref_items = _ensure_index(ref_items)
        self.values = values
        self.ndim = values.ndim

    @property
    def shape(self):
        return self.values.shape

    @property
    def dtype(self):
        return self.values.dtype

    def copy(self, deep=True):
        values = self.values
        if deep:
            values = values.copy()
        return make_block(values, self.items, self.ref_items)

    def merge(self, other):
        assert(self.ref_items.equals(other.ref_items))

        # Not sure whether to allow this or not
        # if not union_ref.equals(other.ref_items):
        #     union_ref = self.ref_items + other.ref_items
        return _merge_blocks([self, other], self.ref_items)

    def reindex_axis(self, indexer, mask, needs_masking, axis=0):
        """
        Reindex using pre-computed indexer information
        """
        if self.values.size > 0:
            new_values = common.take_fast(self.values, indexer, mask,
                                          needs_masking, axis=axis)
        else:
            shape = list(self.shape)
            shape[axis] = len(indexer)
            new_values = np.empty(shape)
            new_values.fill(np.nan)
        return make_block(new_values, self.items, self.ref_items)

    def reindex_items_from(self, new_ref_items):
        """
        Reindex to only those items contained in the input set of items

        E.g. if you have ['a', 'b'], and the input items is ['b', 'c', 'd'],
        then the resulting items will be ['b']

        Returns
        -------
        reindexed : Block
        """
        new_ref_items, indexer = self.items.reindex(new_ref_items)
        mask = indexer != -1
        masked_idx = indexer[mask]
        new_values = self.values.take(masked_idx, axis=0)
        new_items = self.items.take(masked_idx)
        return make_block(new_values, new_items, new_ref_items)

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
        new_items = self.items.delete(loc)
        new_values = np.delete(self.values, loc, 0)
        return make_block(new_values, new_items, self.ref_items)

    def fillna(self, value):
        new_values = self.values.copy()
        mask = common.isnull(new_values.ravel())
        new_values.flat[mask] = value
        return make_block(new_values, self.items, self.ref_items)

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

def make_block(values, items, ref_items, do_integrity_check=False):
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

    return klass(values, items, ref_items, ndim=values.ndim,
                 do_integrity_check=do_integrity_check)

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
    __slots__ = ['axes', 'blocks', 'ndim']

    def __init__(self, blocks, axes, do_integrity_check=True):
        self.axes = [_ensure_index(ax) for ax in axes]
        self.blocks = blocks

        ndim = len(axes)
        for block in blocks:
            assert(ndim == block.values.ndim)

        if do_integrity_check:
            self._verify_integrity()

    @property
    def ndim(self):
        return len(self.axes)

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
                            % (len(value), len(cur_axis)))
        self.axes[axis] = _ensure_index(value)

        if axis == 0:
            for block in self.blocks:
                block.set_ref_items(self.items, maybe_rename=True)

    # make items read only for now
    def _get_items(self):
        return self.axes[0]
    items = property(fget=_get_items)

    def set_items_norename(self, value):
        value = _ensure_index(value)
        self.axes[0] = value

        for block in self.blocks:
            block.set_ref_items(value, maybe_rename=False)

    def __getstate__(self):
        block_values = [b.values for b in self.blocks]
        block_items = [b.items for b in self.blocks]
        axes_array = [ax for ax in self.axes]
        return axes_array, block_values, block_items

    def __setstate__(self, state):
        # discard anything after 3rd, support beta pickling format for a little
        # while longer
        ax_arrays, bvalues, bitems = state[:3]

        self.axes = [_ensure_index(ax) for ax in ax_arrays]
        blocks = []
        for values, items in zip(bvalues, bitems):
            blk = make_block(values, items, self.axes[0],
                             do_integrity_check=True)
            blocks.append(blk)
        self.blocks = blocks

    def __len__(self):
        return len(self.items)

    def __repr__(self):
        output = 'BlockManager'
        for i, ax in enumerate(self.axes):
            if i == 0:
                output += '\nItems: %s' % ax
            else:
                output += '\nAxis %d: %s' % (i, ax)

        for block in self.blocks:
            output += '\n%s' % repr(block)
        return output

    @property
    def shape(self):
        return tuple(len(ax) for ax in self.axes)

    def _verify_integrity(self):
        union_items = _union_block_items(self.blocks)

        mgr_shape = self.shape
        for block in self.blocks:
            assert(block.values.shape[1:] == mgr_shape[1:])
        tot_items = sum(len(x.items) for x in self.blocks)
        assert(len(self.items) == tot_items)

    def astype(self, dtype):
        new_blocks = []
        for block in self.blocks:
            newb = make_block(block.values.astype(dtype), block.items,
                              block.ref_items)
            new_blocks.append(newb)

        new_mgr = BlockManager(new_blocks, self.axes)
        return new_mgr.consolidate()

    def is_consolidated(self):
        """
        Return True if more than one block with the same dtype
        """
        dtypes = [blk.dtype for blk in self.blocks]
        return len(dtypes) == len(set(dtypes))

    def get_slice(self, slobj, axis=0):
        new_axes = list(self.axes)
        new_axes[axis] = new_axes[axis][slobj]

        if axis == 0:
            new_items = new_axes[0]
            if len(self.blocks) == 1:
                blk = self.blocks[0]
                newb = make_block(blk.values[slobj], new_items,
                                  new_items)
                new_blocks = [newb]
            else:
                return self.reindex_items(new_items)
        else:
            new_blocks = self._slice_blocks(slobj, axis)

        return BlockManager(new_blocks, new_axes, do_integrity_check=False)

    def _slice_blocks(self, slobj, axis):
        new_blocks = []

        slicer = [slice(None, None) for _ in range(self.ndim)]
        slicer[axis] = slobj
        slicer = tuple(slicer)

        for block in self.blocks:
            newb = make_block(block.values[slicer], block.items,
                              block.ref_items)
            new_blocks.append(newb)
        return new_blocks

    def get_series_dict(self):
        # For DataFrame
        return _blocks_to_series_dict(self.blocks, self.axes[1])

    @classmethod
    def from_blocks(cls, blocks, index):
        # also checks for overlap
        items = _union_block_items(blocks)
        return BlockManager(blocks, [items, index])

    def __contains__(self, item):
        return item in self.items

    @property
    def nblocks(self):
        return len(self.blocks)

    def copy(self, deep=True):
        """
        Make deep or shallow copy of BlockManager

        Parameters
        ----------
        deep : boolean, default True
            If False, return shallow copy (do not copy data)

        Returns
        -------
        copy : BlockManager
        """
        copy_blocks = [block.copy(deep=deep) for block in self.blocks]
        return BlockManager(copy_blocks, self.axes)

    def as_matrix(self, items=None):
        if len(self.blocks) == 0:
            mat = np.empty(self.shape, dtype=float)
        elif len(self.blocks) == 1:
            blk = self.blocks[0]
            if items is None or blk.items.equals(items):
                # if not, then just call interleave per below
                mat = blk.values
            else:
                mat = self.reindex_items(items).as_matrix()
        else:
            if items is None:
                mat = self._interleave(self.items)
            else:
                mat = self.reindex_items(items).as_matrix()

        return mat

    def _interleave(self, items):
        """
        Return ndarray from blocks with specified item order
        Items must be contained in the blocks
        """
        dtype = _interleaved_dtype(self.blocks)
        items = _ensure_index(items)

        result = np.empty(self.shape, dtype=dtype)
        itemmask = np.zeros(len(items), dtype=bool)

        # By construction, all of the item should be covered by one of the
        # blocks
        for block in self.blocks:
            indexer = items.get_indexer(block.items)
            assert((indexer != -1).all())
            result[indexer] = block.values
            itemmask[indexer] = 1
        assert(itemmask.all())
        return result

    def xs(self, key, axis=1, copy=True):
        assert(axis >= 1)

        loc = self.axes[axis].get_loc(key)
        slicer = [slice(None, None) for _ in range(self.ndim)]
        slicer[axis] = loc
        slicer = tuple(slicer)

        new_axes = list(self.axes)

        # could be an array indexer!
        if isinstance(loc, (slice, np.ndarray)):
            new_axes[axis] = new_axes[axis][loc]
        else:
            new_axes.pop(axis)

        new_blocks = []
        if len(self.blocks) > 1:
            if not copy:
                raise Exception('cannot get view of mixed-type or '
                                'non-consolidated DataFrame')
            for blk in self.blocks:
                newb = make_block(blk.values[slicer], blk.items, blk.ref_items)
                new_blocks.append(newb)
        elif len(self.blocks) == 1:
            vals = self.blocks[0].values[slicer]
            if copy:
                vals = vals.copy()
            new_blocks = [make_block(vals, self.items, self.items)]

        return BlockManager(new_blocks, new_axes)

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
        return BlockManager(new_blocks, self.axes)

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
        if value.ndim == self.ndim - 1:
            value = value.reshape((1,) + value.shape)
        assert(value.shape[1:] == self.shape[1:])
        if item in self.items:
            i, block = self._find_block(item)
            if not block.can_store(value):
                # delete from block, create and append new block
                self._delete_from_block(i, item)
                self._add_new_block(item, value)
            else:
                block.set(item, value)
        else:
            # insert at end
            self.insert(len(self.items), item, value)

    def insert(self, loc, item, value):
        if item in self.items:
            raise Exception('cannot insert %s, already exists' % item)

        new_items = self.items.insert(loc, item)
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

        # hm, elaborate hack?
        loc = self.items.get_loc(item)
        new_block = make_block(value, self.items[loc:loc+1], self.items)
        self.blocks.append(new_block)

    def _find_block(self, item):
        self._check_have(item)
        for i, block in enumerate(self.blocks):
            if item in block:
                return i, block

    def _check_have(self, item):
        if item not in self.items:
            raise KeyError('no item named %s' % str(item))

    def reindex_axis(self, new_axis, method=None, axis=0):
        if axis == 0:
            assert(method is None)
            return self.reindex_items(new_axis)

        new_axis = _ensure_index(new_axis)
        cur_axis = self.axes[axis]

        new_axis, indexer = cur_axis.reindex(new_axis, method)
        mask = indexer == -1

        # TODO: deal with length-0 case? or does it fall out?
        needs_masking = len(new_axis) > 0 and mask.any()

        new_blocks = []
        for block in self.blocks:
            newb = block.reindex_axis(indexer, mask, needs_masking,
                                      axis=axis)
            new_blocks.append(newb)

        new_axes = list(self.axes)
        new_axes[axis] = new_axis
        return BlockManager(new_blocks, new_axes)

    def reindex_items(self, new_items):
        """

        """
        new_items = _ensure_index(new_items)
        data = self
        if not data.is_consolidated():
            data = data.consolidate()
            return data.reindex_items(new_items)

        # TODO: this part could be faster (!)
        new_items, indexer = self.items.reindex(new_items)
        mask = indexer == -1

        new_blocks = []
        for block in self.blocks:
            newb = block.reindex_items_from(new_items)
            if len(newb.items) > 0:
                new_blocks.append(newb)

        if mask.any():
            extra_items = new_items[mask]

            block_shape = list(self.shape)
            block_shape[0] = len(extra_items)
            block_values = np.empty(block_shape, dtype=np.float64)
            block_values.fill(nan)
            na_block = make_block(block_values, extra_items, new_items,
                                  do_integrity_check=True)
            new_blocks.append(na_block)
            new_blocks = _consolidate(new_blocks, new_items)

        new_axes = list(self.axes)
        new_axes[0] = new_items

        return BlockManager(new_blocks, new_axes)

    def take(self, indices, axis=1):
        if axis == 0:
            raise NotImplementedError

        new_axes = list(self.axes)
        new_axes[axis] = self.axes[axis].take(indices)
        new_blocks = []
        for blk in self.blocks:
            newb = make_block(blk.values.take(indices, axis=axis), blk.items,
                              self.items)
            new_blocks.append(newb)

        return BlockManager(new_blocks, new_axes)

    def merge(self, other, lsuffix=None, rsuffix=None):
        assert(self._is_indexed_like(other))

        this, other = self._maybe_rename_join(other, lsuffix, rsuffix)

        cons_items = this.items + other.items
        consolidated = _consolidate(this.blocks + other.blocks, cons_items)

        new_axes = list(this.axes)
        new_axes[0] = cons_items

        return BlockManager(consolidated, new_axes)

    def _maybe_rename_join(self, other, lsuffix, rsuffix, copydata=True):
        intersection = self.items.intersection(other.items)

        if len(intersection) > 0:
            if not lsuffix and not rsuffix:
                raise Exception('columns overlap: %s' % intersection)

            def lrenamer(x):
                if x in intersection:
                    return '%s%s' % (x, lsuffix)
                return x

            def rrenamer(x):
                if x in intersection:
                    return '%s%s' % (x, rsuffix)
                return x

            # XXX: COPIES DATA!
            this = self.rename_items(lrenamer, copydata=copydata)
            other = other.rename_items(rrenamer, copydata=copydata)
        else:
            this = self

        return this, other

    def _is_indexed_like(self, other):
        """
        Check all axes except items
        """
        assert(self.ndim == other.ndim)
        for ax, oax in zip(self.axes[1:], other.axes[1:]):
            if not ax.equals(oax):
                return False
        return True

    def join_on(self, other, on, axis=1, lsuffix=None, rsuffix=None):
        this, other = self._maybe_rename_join(other, lsuffix, rsuffix)

        other_axis = other.axes[axis]
        indexer = lib.merge_indexer_object(on.astype(object),
                                           other_axis.indexMap)

        # TODO: deal with length-0 case? or does it fall out?
        mask = indexer == -1
        needs_masking = len(on) > 0 and mask.any()
        other_blocks = []
        for block in other.blocks:
            newb = block.reindex_axis(indexer, mask, needs_masking, axis=axis)
            other_blocks.append(newb)

        cons_items = this.items + other.items
        consolidated = _consolidate(this.blocks + other_blocks, cons_items)

        new_axes = list(this.axes)
        new_axes[0] = cons_items
        return BlockManager(consolidated, new_axes)

    def rename_axis(self, mapper, axis=1):
        new_axis = Index([mapper(x) for x in self.axes[axis]])
        new_axis._verify_integrity()

        new_axes = list(self.axes)
        new_axes[axis] = new_axis
        return BlockManager(self.blocks, new_axes)

    def rename_items(self, mapper, copydata=True):
        new_items = Index([mapper(x) for x in self.items])
        new_items._verify_integrity()

        new_blocks = []
        for block in self.blocks:
            newb = block.copy(deep=copydata)
            newb.set_ref_items(new_items, maybe_rename=True)
            new_blocks.append(newb)
        new_axes = list(self.axes)
        new_axes[0] = new_items
        return BlockManager(new_blocks, new_axes)

    def add_prefix(self, prefix):
        f = (('%s' % prefix) + '%s').__mod__
        return self.rename_items(f)

    def add_suffix(self, suffix):
        f = ('%s' + ('%s' % suffix)).__mod__
        return self.rename_items(f)

    def fillna(self, value):
        """

        """
        new_blocks = [b.fillna(value) for b in self.blocks]
        return BlockManager(new_blocks, self.axes)

    @property
    def block_id_vector(self):
        # TODO
        result = np.empty(len(self.items), dtype=int)
        result.fill(-1)

        for i, blk in enumerate(self.blocks):
            indexer = self.items.get_indexer(blk.items)
            assert((indexer != -1).all())
            result.put(indexer, i)

        assert((result >= 0).all())
        return result

_data_types = [np.float_, np.int_]
def form_blocks(data, axes):
    # pre-filter out items if we passed it
    items = axes[0]

    if len(data) < len(items):
        extra_items = items - Index(data.keys())
    else:
        extra_items = []

    # put "leftover" items in float bucket, where else?
    # generalize?
    float_dict = {}
    int_dict = {}
    bool_dict = {}
    object_dict = {}
    for k, v in data.iteritems():
        if issubclass(v.dtype.type, np.floating):
            float_dict[k] = v
        elif issubclass(v.dtype.type, np.integer):
            int_dict[k] = v
        elif v.dtype == np.bool_:
            bool_dict[k] = v
        else:
            object_dict[k] = v

    blocks = []
    if len(float_dict):
        float_block = _simple_blockify(float_dict, items, np.float64)
        blocks.append(float_block)

    if len(int_dict):
        int_block = _simple_blockify(int_dict, items, np.int_)
        blocks.append(int_block)

    if len(bool_dict):
        bool_block = _simple_blockify(bool_dict, items, np.bool_)
        blocks.append(bool_block)

    if len(object_dict) > 0:
        object_block = _simple_blockify(object_dict, items, np.object_)
        blocks.append(object_block)

    if len(extra_items):
        shape = (len(extra_items),) + tuple(len(x) for x in axes[1:])
        block_values = np.empty(shape, dtype=float)
        block_values.fill(nan)

        na_block = make_block(block_values, extra_items, items,
                              do_integrity_check=True)
        blocks.append(na_block)
        blocks = _consolidate(blocks, items)

    return blocks

def _simple_blockify(dct, ref_items, dtype):
    block_items, values = _stack_dict(dct, ref_items)
    # CHECK DTYPE?
    if values.dtype != dtype: # pragma: no cover
        values = values.astype(dtype)

    return make_block(values, block_items, ref_items, do_integrity_check=True)

def _stack_dict(dct, ref_items):
    items = [x for x in ref_items if x in dct]
    stacked = np.vstack([np.asarray(dct[k]) for k in items])
    return items, stacked

def _blocks_to_series_dict(blocks, index=None):
    from pandas.core.series import Series

    series_dict = {}

    for block in blocks:
        for item, vec in zip(block.items, block.values):
            series_dict[item] = Series(vec, index=index)
    return series_dict

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

# TODO: this could be much optimized

def _merge_blocks(blocks, items):
    if len(blocks) == 1:
        return blocks[0]
    new_values = np.vstack([b.values for b in blocks])
    new_items = np.concatenate([b.items for b in blocks])
    new_block = make_block(new_values, new_items, items,
                           do_integrity_check=True)
    return new_block.reindex_items_from(items)

def _union_block_items(blocks):
    tot_len = 0
    all_items = []
    slow = False
    for b in blocks:
        tot_len += len(b.items)
        if type(b.items) != Index:
            slow = True
        all_items.append(b.items)

    if slow:
        the_union = _union_items_slow(all_items)
    else:
        the_union = Index(lib.fast_unique_multiple(all_items))

    if tot_len > len(the_union):
        raise Exception('item names overlap')
    return the_union

def _union_items_slow(all_items):
    seen = None
    for items in all_items:
        if seen is None:
            seen = items
        else:
            seen = seen.union(items)
    return seen


def join_managers(left, right, axis=1, how='left'):
    """
    Parameters
    ----------
    other
    lindexer
    lmask
    rindexer
    rmask

    Returns
    -------
    merged : BlockManager
    """
    assert(left.is_consolidated())
    assert(right.is_consolidated())

    laxis = left.axes[axis]
    raxis = right.axes[axis]

    join_index, lindexer, rindexer = laxis.join(raxis, how=how,
                                                return_indexers=True)

    N = len(join_index)

    if lindexer is None:
        lmask = None
        lneed_masking = None
    else:
        lmask = lindexer == -1
        lneed_masking = lmask.any()

    if rindexer is None:
        rmask = None
        rneed_masking = None
    else:
        rmask = rindexer == -1
        rneed_masking = rmask.any()

    lblocks = _maybe_upcast_blocks(left.blocks, lneed_masking)
    rblocks = _maybe_upcast_blocks(right.blocks, rneed_masking)

    left_blockmap = dict((type(blk), blk) for blk in lblocks)
    right_blockmap = dict((type(blk), blk) for blk in rblocks)

    # do NOT sort
    result_items = left.items.append(right.items)

    result_axes = list(left.axes)
    result_axes[0] = result_items
    result_axes[axis] = join_index

    result_blocks = []

    # copies all data by definition

    kinds = set(left_blockmap) | set(right_blockmap)
    for klass in kinds:
        if klass in left_blockmap and klass in right_blockmap:
            # true merge, do not produce intermediate copy
            lblk = left_blockmap[klass]
            rblk = right_blockmap[klass]
            new_values = _merge_blocks_fast(lblk, rblk,
                                            lindexer, lmask, lneed_masking,
                                            rindexer, rmask, rneed_masking,
                                            axis=axis)
            new_items = lblk.items.append(rblk.items)
            res_blk = make_block(new_values, new_items, result_items)
        elif klass in left_blockmap:
            # only take necessary
            blk = left_blockmap[klass]
            if lindexer is None:
                res_blk = blk.copy()
            else:
                res_blk = blk.reindex_axis(lindexer, lmask, lneed_masking,
                                           axis=axis)
            res_blk.ref_items = result_items
        elif klass in right_blockmap:
            # only take necessary
            blk = right_blockmap[klass]
            if rindexer is None:
                res_blk = blk.copy()
            else:
                res_blk = blk.reindex_axis(rindexer, rmask, rneed_masking,
                                           axis=axis)
            res_blk.ref_items = result_items

        result_blocks.append(res_blk)

    return BlockManager(result_blocks, result_axes)

def _maybe_upcast_blocks(blocks, needs_masking):
    """
    Upcast and consolidate if necessary
    """
    if not needs_masking:
        return blocks
    new_blocks = []
    for block in blocks:
        if isinstance(block, IntBlock):
            newb = make_block(block.values.astype(float), block.items,
                              block.ref_items)
        elif isinstance(block, BoolBlock):
            newb = make_block(block.values.astype(object), block.items,
                              block.ref_items)
        else:
            newb = block
        new_blocks.append(newb)

    # use any ref_items
    return _consolidate(new_blocks, newb.ref_items)

def _merge_blocks_fast(left, right, lindexer, lmask, lneed_masking,
                       rindexer, rmask, rneed_masking, axis=1):

    n = left.values.shape[axis] if lindexer is None else len(lindexer)
    lk = len(left.items)
    rk = len(right.items)

    out_shape = list(left.shape)
    out_shape[0] = lk + rk
    out_shape[axis] = n

    out = np.empty(out_shape, dtype=left.values.dtype)

    if lindexer is None:
        common.take_fast(left.values, np.arange(n, dtype=np.int32),
                         None, False, axis=axis, out=out[:lk])
    else:
        common.take_fast(left.values, lindexer, lmask, lneed_masking,
                         axis=axis, out=out[:lk])

    if rindexer is None:
        common.take_fast(right.values, np.arange(n, dtype=np.int32),
                         None, False, axis=axis, out=out[lk:])
    else:
        common.take_fast(right.values, rindexer, rmask, rneed_masking,
                         axis=axis, out=out[lk:])
    return out
