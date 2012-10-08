import itertools
from datetime import datetime

from numpy import nan
import numpy as np

from pandas.core.index import Index, _ensure_index, _handle_legacy_indexes
import pandas.core.common as com
import pandas.lib as lib

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
            indexer = com._ensure_platform_int(indexer)
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

    def reindex_axis(self, indexer, mask, needs_masking, axis=0,
                     fill_value=np.nan):
        """
        Reindex using pre-computed indexer information
        """
        if self.values.size > 0:
            new_values = com.take_fast(self.values, indexer, mask,
                                       needs_masking, axis=axis,
                                       fill_value=fill_value)
        else:
            shape = list(self.shape)
            shape[axis] = len(indexer)
            new_values = np.empty(shape)
            new_values.fill(fill_value)
        return make_block(new_values, self.items, self.ref_items)

    def reindex_items_from(self, new_ref_items, copy=True):
        """
        Reindex to only those items contained in the input set of items

        E.g. if you have ['a', 'b'], and the input items is ['b', 'c', 'd'],
        then the resulting items will be ['b']

        Returns
        -------
        reindexed : Block
        """
        new_ref_items, indexer = self.items.reindex(new_ref_items)
        if indexer is None:
            new_items = new_ref_items
            new_values = self.values.copy() if copy else self.values
        else:
            mask = indexer != -1
            masked_idx = indexer[mask]

            if self.values.ndim == 2:
                new_values = com.take_2d(self.values, masked_idx, axis=0,
                                         needs_masking=False)
            else:
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

    def split_block_at(self, item):
        """
        Split block around given column, for "deleting" a column without
        having to copy data by returning views on the original array

        Returns
        -------
        leftb, rightb : (Block or None, Block or None)
        """
        loc = self.items.get_loc(item)

        if len(self.items) == 1:
            # no blocks left
            return None, None

        if loc == 0:
            # at front
            left_block = None
            right_block = make_block(self.values[1:], self.items[1:].copy(),
                                      self.ref_items)
        elif loc == len(self.values) - 1:
            # at back
            left_block = make_block(self.values[:-1], self.items[:-1].copy(),
                                    self.ref_items)
            right_block = None
        else:
            # in the middle
            left_block = make_block(self.values[:loc],
                                    self.items[:loc].copy(), self.ref_items)
            right_block = make_block(self.values[loc + 1:],
                                     self.items[loc + 1:].copy(),
                                     self.ref_items)

        return left_block, right_block

    def fillna(self, value, inplace=False):
        new_values = self.values if inplace else self.values.copy()

        mask = com.isnull(new_values)
        np.putmask(new_values, mask, value)

        if inplace:
            return self
        else:
            return make_block(new_values, self.items, self.ref_items)

    def _can_hold_element(self, value):
        raise NotImplementedError()

    def _try_cast(self, value):
        raise NotImplementedError()

    def replace(self, to_replace, value, inplace=False):
        new_values = self.values if inplace else self.values.copy()
        if self._can_hold_element(value):
            value = self._try_cast(value)

        if not isinstance(to_replace, (list, np.ndarray)):
            if self._can_hold_element(to_replace):
                to_replace = self._try_cast(to_replace)
                np.putmask(new_values, com.mask_missing(new_values, to_replace),
                           value)
        else:
            try:
                to_replace = np.array(to_replace, dtype=self.dtype)
                np.putmask(new_values, com.mask_missing(new_values, to_replace),
                           value)
            except:
                to_replace = np.array(to_replace, dtype=object)
                for r in to_replace:
                    if self._can_hold_element(r):
                        r = self._try_cast(r)
                np.putmask(new_values, com.mask_missing(new_values, to_replace),
                           value)

        if inplace:
            return self
        else:
            return make_block(new_values, self.items, self.ref_items)

    def putmask(self, mask, new, inplace=False):
        new_values = self.values if inplace else self.values.copy()
        if self._can_hold_element(new):
            new = self._try_cast(new)
            np.putmask(new_values, mask, new)
        if inplace:
            return self
        else:
            return make_block(new_values, self.items, self.ref_items)

    def interpolate(self, method='pad', axis=0, inplace=False,
                    limit=None, missing=None):
        values = self.values if inplace else self.values.copy()

        if values.ndim != 2:
            raise NotImplementedError

        transf = (lambda x: x) if axis == 0 else (lambda x: x.T)

        if missing is None:
            mask = None
        else: # todo create faster fill func without masking
            mask = _mask_missing(transf(values), missing)

        if method == 'pad':
            com.pad_2d(transf(values), limit=limit, mask=mask)
        else:
            com.backfill_2d(transf(values), limit=limit, mask=mask)

        return make_block(values, self.items, self.ref_items)

    def take(self, indexer, axis=1, fill_value=np.nan):
        assert(axis >= 1)
        new_values = com.take_fast(self.values, indexer, None,
                                   None, axis=axis,
                                   fill_value=fill_value)
        return make_block(new_values, self.items, self.ref_items)

    def get_values(self, dtype):
        return self.values

def _mask_missing(array, missing_values):
    if not isinstance(missing_values, (list, np.ndarray)):
        missing_values = [missing_values]

    mask = None
    missing_values = np.array(missing_values, dtype=object)
    if com.isnull(missing_values).any():
        mask = com.isnull(array)
        missing_values = missing_values[com.notnull(missing_values)]

    for v in missing_values:
        if mask is None:
            mask = array == missing_values
        else:
            mask |= array == missing_values
    return mask

#-------------------------------------------------------------------------------
# Is this even possible?

class FloatBlock(Block):
    _can_hold_na = True

    def _can_hold_element(self, element):
        return isinstance(element, (float, int))

    def _try_cast(self, element):
        try:
            return float(element)
        except: # pragma: no cover
            return element

    def should_store(self, value):
        # when inserting a column should not coerce integers to floats
        # unnecessarily
        return issubclass(value.dtype.type, np.floating)

class ComplexBlock(Block):
    _can_hold_na = True

    def _can_hold_element(self, element):
        return isinstance(element, complex)

    def _try_cast(self, element):
        try:
            return complex(element)
        except: # pragma: no cover
            return element

    def should_store(self, value):
        return issubclass(value.dtype.type, np.complexfloating)

class IntBlock(Block):
    _can_hold_na = False

    def _can_hold_element(self, element):
        return com.is_integer(element)

    def _try_cast(self, element):
        try:
            return int(element)
        except: # pragma: no cover
            return element

    def should_store(self, value):
        return issubclass(value.dtype.type, np.integer)

class BoolBlock(Block):
    _can_hold_na = False

    def _can_hold_element(self, element):
        return isinstance(element, (int, bool))

    def _try_cast(self, element):
        try:
            return bool(element)
        except: # pragma: no cover
            return element

    def should_store(self, value):
        return issubclass(value.dtype.type, np.bool_)

class ObjectBlock(Block):
    _can_hold_na = True

    def _can_hold_element(self, element):
        return True

    def _try_cast(self, element):
        return element

    def should_store(self, value):
        return not issubclass(value.dtype.type,
                              (np.integer, np.floating, np.complexfloating,
                               np.datetime64, np.bool_))

_NS_DTYPE = np.dtype('M8[ns]')

class DatetimeBlock(Block):
    _can_hold_na = True

    def __init__(self, values, items, ref_items, ndim=2,
                 do_integrity_check=False):
        if values.dtype != _NS_DTYPE:
            values = lib.cast_to_nanoseconds(values)

        Block.__init__(self, values, items, ref_items, ndim=ndim,
                       do_integrity_check=do_integrity_check)

    def _can_hold_element(self, element):
        return com.is_integer(element) or isinstance(element, datetime)

    def _try_cast(self, element):
        try:
            return int(element)
        except:
            return element

    def should_store(self, value):
        return issubclass(value.dtype.type, np.datetime64)

    def set(self, item, value):
        """
        Modify Block in-place with new item value

        Returns
        -------
        None
        """
        loc = self.items.get_loc(item)

        if value.dtype != _NS_DTYPE:
            value = lib.cast_to_nanoseconds(value)

        self.values[loc] = value

    def get_values(self, dtype):
        if dtype == object:
            flat_i8 = self.values.ravel().view(np.int64)
            res = lib.ints_to_pydatetime(flat_i8)
            return res.reshape(self.values.shape)
        return self.values


def make_block(values, items, ref_items, do_integrity_check=False):
    dtype = values.dtype
    vtype = dtype.type

    if issubclass(vtype, np.floating):
        klass = FloatBlock
    elif issubclass(vtype, np.complexfloating):
        klass = ComplexBlock
    elif issubclass(vtype, np.datetime64):
        klass = DatetimeBlock
    elif issubclass(vtype, np.integer):
        if vtype != np.int64:
            values = values.astype('i8')
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

    @classmethod
    def make_empty(self):
        return BlockManager([], [[], []])

    def __nonzero__(self):
        return True

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
        self.axes = _handle_legacy_indexes(self.axes)

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
        # _union_block_items(self.blocks)
        mgr_shape = self.shape
        for block in self.blocks:
            assert(block.ref_items is self.items)
            assert(block.values.shape[1:] == mgr_shape[1:])
        tot_items = sum(len(x.items) for x in self.blocks)
        assert(len(self.items) == tot_items)

    def astype(self, dtype):
        new_blocks = []
        for block in self.blocks:
            newb = make_block(com._astype_nansafe(block.values, dtype),
                              block.items, block.ref_items)
            new_blocks.append(newb)

        new_mgr = BlockManager(new_blocks, self.axes)
        return new_mgr.consolidate()

    def is_consolidated(self):
        """
        Return True if more than one block with the same dtype
        """
        dtypes = [blk.dtype.type for blk in self.blocks]
        return len(dtypes) == len(set(dtypes))

    def get_numeric_data(self, copy=False, type_list=None):
        """
        Parameters
        ----------
        copy : boolean, default False
            Whether to copy the blocks
        type_list : tuple of type, default None
            Numeric types by default (Float/Complex/Int but not Datetime)
        """
        if type_list is None:
            def filter_blocks(block):
                return (isinstance(block, (IntBlock, FloatBlock, ComplexBlock))
                        and not isinstance(block, DatetimeBlock))
        else:
            type_list = self._get_clean_block_types(type_list)
            filter_blocks = lambda block: isinstance(block, type_list)

        maybe_copy = lambda b: b.copy() if copy else b
        num_blocks = [maybe_copy(b) for b in self.blocks if filter_blocks(b)]

        if len(num_blocks) == 0:
            return BlockManager.make_empty()

        indexer = np.sort(np.concatenate([b.ref_locs for b in num_blocks]))
        new_items = self.items.take(indexer)

        new_blocks = []
        for b in num_blocks:
            b = b.copy(deep=False)
            b.ref_items = new_items
            new_blocks.append(b)
        new_axes = list(self.axes)
        new_axes[0] = new_items
        return BlockManager(new_blocks, new_axes, do_integrity_check=False)

    def _get_clean_block_types(self, type_list):
        if not isinstance(type_list, tuple):
            try:
                type_list = tuple(type_list)
            except TypeError:
                type_list = (type_list,)

        type_map = {int : IntBlock, float : FloatBlock,
                    complex : ComplexBlock,
                    np.datetime64 : DatetimeBlock,
                    datetime : DatetimeBlock,
                    bool : BoolBlock,
                    object : ObjectBlock}

        type_list = tuple([type_map.get(t, t) for t in type_list])
        return type_list

    def get_bool_data(self, copy=False):
        return self.get_numeric_data(copy=copy, type_list=(BoolBlock,))

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
        # copy_axes = [ax.copy() for ax in self.axes]
        copy_axes = list(self.axes)
        return BlockManager(copy_blocks, copy_axes, do_integrity_check=False)

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
            result[indexer] = block.get_values(dtype)
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

    def fast_2d_xs(self, loc, copy=False):
        """

        """
        if len(self.blocks) == 1:
            result = self.blocks[0].values[:, loc]
            if copy:
                result = result.copy()
            return result

        if not copy:
            raise Exception('cannot get view of mixed-type or '
                            'non-consolidated DataFrame')

        dtype = _interleaved_dtype(self.blocks)

        items = self.items
        n = len(items)
        result = np.empty(n, dtype=dtype)
        for blk in self.blocks:
            values = blk.values
            for j, item in enumerate(blk.items):
                i = items.get_loc(item)
                result[i] = values[j, loc]

        return result

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

    def _consolidate_inplace(self):
        self.blocks = _consolidate(self.blocks, self.items)

    def get(self, item):
        _, block = self._find_block(item)
        return block.get(item)

    def iget(self, i):
        item = self.items[i]
        if self.items.is_unique:
            return self.get(item)
        else:
            # ugh
            try:
                inds, = (self.items == item).nonzero()
            except AttributeError: #MultiIndex
                inds, = self.items.map(lambda x: x == item).nonzero()

            _, block = self._find_block(item)

            try:
                binds, = (block.items == item).nonzero()
            except AttributeError: #MultiIndex
                binds, = block.items.map(lambda x: x == item).nonzero()

            for j, (k, b) in enumerate(zip(inds, binds)):
                if i == k:
                    return block.values[b]

            raise Exception('Cannot have duplicate column names '
                            'split across dtypes')

    def get_scalar(self, tup):
        """
        Retrieve single item
        """
        item = tup[0]
        _, blk = self._find_block(item)

        # this could obviously be seriously sped up in cython
        item_loc = blk.items.get_loc(item),
        full_loc = item_loc + tuple(ax.get_loc(x)
                                    for ax, x in zip(self.axes[1:], tup[1:]))
        return blk.values[full_loc]

    def delete(self, item):
        i, _ = self._find_block(item)
        loc = self.items.get_loc(item)

        new_items = self.items.delete(loc)

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
            if not block.should_store(value):
                # delete from block, create and append new block
                self._delete_from_block(i, item)
                self._add_new_block(item, value, loc=None)
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
        self._add_new_block(item, value, loc=loc)

        if len(self.blocks) > 100:
            self._consolidate_inplace()

    def set_items_norename(self, value):
        value = _ensure_index(value)
        self.axes[0] = value

        for block in self.blocks:
            block.set_ref_items(value, maybe_rename=False)

    def _delete_from_block(self, i, item):
        """
        Delete and maybe remove the whole block
        """
        block = self.blocks.pop(i)
        new_left, new_right = block.split_block_at(item)

        if new_left is not None:
            self.blocks.append(new_left)

        if new_right is not None:
            self.blocks.append(new_right)

    def _add_new_block(self, item, value, loc=None):
        # Do we care about dtype at the moment?

        # hm, elaborate hack?
        if loc is None:
            loc = self.items.get_loc(item)
        new_block = make_block(value, self.items[loc:loc+1].copy(),
                               self.items)
        self.blocks.append(new_block)

    def _find_block(self, item):
        self._check_have(item)
        for i, block in enumerate(self.blocks):
            if item in block:
                return i, block

    def _check_have(self, item):
        if item not in self.items:
            raise KeyError('no item named %s' % str(item))

    def reindex_axis(self, new_axis, method=None, axis=0, copy=True):
        new_axis = _ensure_index(new_axis)
        cur_axis = self.axes[axis]

        if new_axis.equals(cur_axis):
            if copy:
                result = self.copy(deep=True)
                result.axes[axis] = new_axis

                if axis == 0:
                    # patch ref_items, #1823
                    for blk in result.blocks:
                        blk.ref_items = new_axis

                return result
            else:
                return self

        if axis == 0:
            assert(method is None)
            return self.reindex_items(new_axis)

        new_axis, indexer = cur_axis.reindex(new_axis, method)
        return self.reindex_indexer(new_axis, indexer, axis=axis)

    def reindex_indexer(self, new_axis, indexer, axis=1, fill_value=np.nan):
        """
        pandas-indexer with -1's only.
        """
        if axis == 0:
            return self._reindex_indexer_items(new_axis, indexer, fill_value)

        mask = indexer == -1

        # TODO: deal with length-0 case? or does it fall out?
        needs_masking = len(new_axis) > 0 and mask.any()

        new_blocks = []
        for block in self.blocks:
            newb = block.reindex_axis(indexer, mask, needs_masking,
                                      axis=axis, fill_value=fill_value)
            new_blocks.append(newb)

        new_axes = list(self.axes)
        new_axes[axis] = new_axis
        return BlockManager(new_blocks, new_axes)

    def _reindex_indexer_items(self, new_items, indexer, fill_value):
        # TODO: less efficient than I'd like

        item_order = com.take_1d(self.items.values, indexer)

        # keep track of what items aren't found anywhere
        mask = np.zeros(len(item_order), dtype=bool)

        new_blocks = []
        for blk in self.blocks:
            blk_indexer = blk.items.get_indexer(item_order)
            selector = blk_indexer != -1
            # update with observed items
            mask |= selector

            if not selector.any():
                continue

            new_block_items = new_items.take(selector.nonzero()[0])
            new_values = com.take_fast(blk.values, blk_indexer[selector],
                                       None, False, axis=0)
            new_blocks.append(make_block(new_values, new_block_items,
                                         new_items))

        if not mask.all():
            na_items = new_items[-mask]
            na_block = self._make_na_block(na_items, new_items,
                                           fill_value=fill_value)
            new_blocks.append(na_block)
            new_blocks = _consolidate(new_blocks, new_items)

        return BlockManager(new_blocks, [new_items] + self.axes[1:])

    def reindex_items(self, new_items, copy=True, fill_value=np.nan):
        """

        """
        new_items = _ensure_index(new_items)
        data = self
        if not data.is_consolidated():
            data = data.consolidate()
            return data.reindex_items(new_items)

        # TODO: this part could be faster (!)
        new_items, indexer = self.items.reindex(new_items)

        # could have some pathological (MultiIndex) issues here
        new_blocks = []
        if indexer is None:
            for blk in self.blocks:
                if copy:
                    new_blocks.append(blk.reindex_items_from(new_items))
                else:
                    blk.ref_items = new_items
                    new_blocks.append(blk)
        else:
            for block in self.blocks:
                newb = block.reindex_items_from(new_items, copy=copy)
                if len(newb.items) > 0:
                    new_blocks.append(newb)

            mask = indexer == -1
            if mask.any():
                extra_items = new_items[mask]
                na_block = self._make_na_block(extra_items, new_items,
                                               fill_value=fill_value)
                new_blocks.append(na_block)
                new_blocks = _consolidate(new_blocks, new_items)

        return BlockManager(new_blocks, [new_items] + self.axes[1:])

    def _make_na_block(self, items, ref_items, fill_value=np.nan):
        # TODO: infer dtypes other than float64 from fill_value

        block_shape = list(self.shape)
        block_shape[0] = len(items)

        dtype = com._infer_dtype(fill_value)
        block_values = np.empty(block_shape, dtype=dtype)
        block_values.fill(fill_value)
        na_block = make_block(block_values, items, ref_items,
                              do_integrity_check=True)
        return na_block

    def take(self, indexer, axis=1):
        if axis == 0:
            raise NotImplementedError

        indexer = np.asarray(indexer, dtype='i4')

        n = len(self.axes[axis])
        if ((indexer == -1) | (indexer >= n)).any():
            raise Exception('Indices must be nonzero and less than '
                            'the axis length')

        new_axes = list(self.axes)
        new_axes[axis] = self.axes[axis].take(indexer)
        new_blocks = []
        for blk in self.blocks:
            new_values = com.take_fast(blk.values, indexer,
                                       None, False, axis=axis)
            newb = make_block(new_values, blk.items, self.items)
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
        to_rename = self.items.intersection(other.items)
        if len(to_rename) > 0:
            if not lsuffix and not rsuffix:
                raise Exception('columns overlap: %s' % to_rename)

            def lrenamer(x):
                if x in to_rename:
                    return '%s%s' % (x, lsuffix)
                return x

            def rrenamer(x):
                if x in to_rename:
                    return '%s%s' % (x, rsuffix)
                return x

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

    def rename_axis(self, mapper, axis=1):
        new_axis = Index([mapper(x) for x in self.axes[axis]])
        assert(new_axis.is_unique)

        new_axes = list(self.axes)
        new_axes[axis] = new_axis
        return BlockManager(self.blocks, new_axes)

    def rename_items(self, mapper, copydata=True):
        new_items = Index([mapper(x) for x in self.items])
        new_items.is_unique

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

    def fillna(self, value, inplace=False):
        new_blocks = [b.fillna(value, inplace=inplace)
                      if b._can_hold_na else b
                      for b in self.blocks]
        if inplace:
            return self
        return BlockManager(new_blocks, self.axes)

    def replace(self, to_replace, value, inplace=False):
        new_blocks = [b.replace(to_replace, value, inplace=inplace)
                      for b in self.blocks]
        if inplace:
            return self
        return BlockManager(new_blocks, self.axes)

    def _replace_list(self, src_lst, dest_lst):
        sset = set(src_lst)
        if any([k in sset for k in dest_lst]):
            masks = {}
            for s in src_lst:
                masks[s] = [b.values == s for b in self.blocks]

            for s, d in zip(src_lst, dest_lst):
                [b.putmask(masks[s][i], d, inplace=True) for i, b in
                 enumerate(self.blocks)]
        else:
            for s, d in zip(src_lst, dest_lst):
                self.replace(s, d, inplace=True)

        return self

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

    @property
    def item_dtypes(self):
        result = np.empty(len(self.items), dtype='O')
        mask = np.zeros(len(self.items), dtype=bool)
        for i, blk in enumerate(self.blocks):
            indexer = self.items.get_indexer(blk.items)
            result.put(indexer, blk.values.dtype.name)
            mask.put(indexer, 1)
        assert(mask.all())
        return result

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
    complex_dict = {}
    int_dict = {}
    bool_dict = {}
    object_dict = {}
    datetime_dict = {}
    for k, v in data.iteritems():
        if issubclass(v.dtype.type, np.floating):
            float_dict[k] = v
        elif issubclass(v.dtype.type, np.complexfloating):
            complex_dict[k] = v
        elif issubclass(v.dtype.type, np.datetime64):
            datetime_dict[k] = v
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

    if len(complex_dict):
        complex_block = _simple_blockify(complex_dict, items, np.complex128)
        blocks.append(complex_block)

    if len(int_dict):
        int_block = _simple_blockify(int_dict, items, np.int64)
        blocks.append(int_block)

    for k, v in list(datetime_dict.items()):
        # hackeroo
        if hasattr(v, 'tz') and v.tz is not None:
            del datetime_dict[k]
            object_dict[k] = v.asobject

    if len(datetime_dict):
        datetime_block = _simple_blockify(datetime_dict, items,
                                          np.dtype('M8[ns]'))
        blocks.append(datetime_block)

    if len(bool_dict):
        bool_block = _simple_blockify(bool_dict, items, np.bool_)
        blocks.append(bool_block)

    if len(object_dict) > 0:
        object_block = _simple_blockify(object_dict, items, np.object_)
        blocks.append(object_block)

    if len(extra_items):
        shape = (len(extra_items),) + tuple(len(x) for x in axes[1:])

        # empty items -> dtype object
        block_values = np.empty(shape, dtype=object)

        block_values.fill(nan)

        na_block = make_block(block_values, extra_items, items,
                              do_integrity_check=True)
        blocks.append(na_block)
        blocks = _consolidate(blocks, items)

    return blocks

def _simple_blockify(dct, ref_items, dtype):
    block_items, values = _stack_dict(dct, ref_items, dtype)
    # CHECK DTYPE?
    if values.dtype != dtype: # pragma: no cover
        values = values.astype(dtype)

    return make_block(values, block_items, ref_items, do_integrity_check=True)

def _stack_dict(dct, ref_items, dtype):
    from pandas.core.series import Series

    # fml
    def _asarray_compat(x):
        # asarray shouldn't be called on SparseSeries
        if isinstance(x, Series):
            return x.values
        else:
            return np.asarray(x)

    def _shape_compat(x):
        # sparseseries
        if isinstance(x, Series):
            return len(x),
        else:
            return x.shape

    # index may box values
    items = ref_items[[x in dct for x in ref_items]]

    first = dct[items[0]]
    shape = (len(dct),) + _shape_compat(first)

    stacked = np.empty(shape, dtype=dtype)
    for i, item in enumerate(items):
        stacked[i] = _asarray_compat(dct[item])

    # stacked = np.vstack([_asarray_compat(dct[k]) for k in items])
    return items, stacked

def _blocks_to_series_dict(blocks, index=None):
    from pandas.core.series import Series

    series_dict = {}

    for block in blocks:
        for item, vec in zip(block.items, block.values):
            series_dict[item] = Series(vec, index=index, name=item)
    return series_dict

def _interleaved_dtype(blocks):
    from collections import defaultdict
    counts = defaultdict(lambda: 0)
    for x in blocks:
        counts[type(x)] += 1

    have_int = counts[IntBlock] > 0
    have_bool = counts[BoolBlock] > 0
    have_object = counts[ObjectBlock] > 0
    have_float = counts[FloatBlock] > 0
    have_complex = counts[ComplexBlock] > 0
    have_dt64 = counts[DatetimeBlock] > 0
    have_numeric = have_float or have_complex or have_int

    if (have_object or
        (have_bool and have_numeric) or
        (have_numeric and have_dt64)):
        return np.dtype(object)
    elif have_bool:
        return np.dtype(bool)
    elif have_int and not have_float and not have_complex:
        return np.dtype('i8')
    elif have_dt64 and not have_float and not have_complex:
        return np.dtype('M8[ns]')
    elif have_complex:
        return np.dtype('c16')
    else:
        return np.dtype('f8')

def _consolidate(blocks, items):
    """
    Merge blocks having same dtype
    """
    get_dtype = lambda x: x.dtype.name

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
    new_values = _vstack([b.values for b in blocks])
    new_items = blocks[0].items.append([b.items for b in blocks[1:]])
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

def _vstack(to_stack):
    if all(x.dtype == _NS_DTYPE for x in to_stack):
        # work around NumPy 1.6 bug
        new_values = np.vstack([x.view('i8') for x in to_stack])
        return new_values.view(_NS_DTYPE)
    else:
        return np.vstack(to_stack)
