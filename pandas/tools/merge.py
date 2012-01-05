"""
SQL-style merge routines
"""

import numpy as np

from pandas.core.frame import DataFrame, _merge_doc
from pandas.core.groupby import get_group_index
from pandas.core.index import Index, MultiIndex, _get_combined_index
from pandas.core.internals import (IntBlock, BoolBlock, BlockManager,
                                   make_block, _consolidate)
from pandas.sparse.frame import SparseDataFrame
from pandas.util.decorators import cache_readonly
import pandas.core.common as com

import pandas._tseries as lib

def merge(left, right, how='inner', on=None, left_on=None, right_on=None,
          left_index=False, right_index=False, sort=True,
          suffixes=('.x', '.y'), copy=True):
    op = _MergeOperation(left, right, how=how, on=on, left_on=left_on,
                         right_on=right_on, left_index=left_index,
                         right_index=right_index, sort=sort, suffixes=suffixes,
                         copy=copy)
    return op.get_result()
if __debug__: merge.__doc__ = _merge_doc % '\nleft : DataFrame'


def concat(objs, axis=0, join='outer', join_axes=None,
           ignore_index=False, verify_integrity=False):
    """
    Concatenate DataFrame objects row or column wise

    Parameters
    ----------
    objs : list of DataFrame objects
    axis : {0, 1}, default 0
        The axis to concatenate along
    join : {'inner', 'outer'}, default 'outer'
        How to handle indexes on other axis(es)
    join_index : index-like
    verify_integrity : boolean, default False
        Check whether the new concatenated axis contains duplicates. This can
        be very expensive relative to the actual data concatenation

    Returns
    -------
    concatenated : DataFrame
    """
    op = Concatenator(objs, axis=axis, join_axes=join_axes,
                      ignore_index=ignore_index, join=join,
                      verify_integrity=verify_integrity)
    return op.get_result()



# TODO: NA group handling
# TODO: transformations??
# TODO: only copy DataFrames when modification necessary

class _MergeOperation(object):
    """
    Perform a database (SQL) merge operation between two DataFrame objects
    using either columns as keys or their row indexes
    """

    def __init__(self, left, right, how='inner', on=None,
                 left_on=None, right_on=None, axis=1,
                 left_index=False, right_index=False, sort=True,
                 suffixes=('.x', '.y'), copy=True):
        self.left = self.orig_left = left
        self.right = self.orig_right = right
        self.how = how
        self.axis = axis

        self.on = _maybe_make_list(on)
        self.left_on = _maybe_make_list(left_on)
        self.right_on = _maybe_make_list(right_on)

        self.drop_keys = False # set this later...kludge

        self.copy = copy

        self.suffixes = suffixes

        self.sort = sort

        self.left_index = left_index
        self.right_index = right_index

        # note this function has side effects
        (self.left_join_keys,
         self.right_join_keys,
         self.join_names) = self._get_merge_keys()

    def get_result(self):
        join_index, left_indexer, right_indexer = self._get_join_info()

        # this is a bit kludgy
        ldata, rdata = self._get_merge_data(self.join_names)

        # TODO: more efficiently handle group keys to avoid extra consolidation!
        join_op = _BlockJoinOperation([ldata, rdata], join_index,
                                      [left_indexer, right_indexer], axis=1,
                                      copy=self.copy)

        result_data = join_op.get_result()
        result = DataFrame(result_data)

        self._maybe_add_join_keys(result, left_indexer, right_indexer)

        return result

    def _maybe_add_join_keys(self, result, left_indexer, right_indexer):
        if not self.drop_keys:
            # do nothing, already found in one of the DataFrames
            return

        # insert group keys
        for i, name in enumerate(self.join_names):
            # a faster way?
            key_col = com.take_1d(self.left_join_keys[i], left_indexer)
            na_indexer = (left_indexer == -1).nonzero()[0]
            right_na_indexer = right_indexer.take(na_indexer)
            key_col.put(na_indexer, com.take_1d(self.right_join_keys[i],
                                                right_na_indexer))
            result.insert(i, name, key_col)

    def _get_join_info(self):
        left_ax = self.left._data.axes[self.axis]
        right_ax = self.right._data.axes[self.axis]
        if self.left_index and self.right_index:
            join_index, left_indexer, right_indexer = \
                left_ax.join(right_ax, how=self.how, return_indexers=True)
        elif self.right_index and self.how == 'left':
            join_index = left_ax
            left_indexer = None

            if len(self.left_join_keys) > 1:
                assert(isinstance(right_ax, MultiIndex) and
                       len(self.left_join_keys) == right_ax.nlevels)

                right_indexer = _get_multiindex_indexer(self.left_join_keys,
                                                        right_ax, sort=False)
            else:
                right_indexer = right_ax.get_indexer(self.left_join_keys[0])

        elif self.left_index and self.how == 'right':
            join_index = right_ax
            right_indexer = None

            if len(self.right_join_keys) > 1:
                assert(isinstance(left_ax, MultiIndex) and
                       len(self.right_join_keys) == left_ax.nlevels)
                left_indexer = _get_multiindex_indexer(self.right_join_keys,
                                                       left_ax, sort=False)
            else:
                left_indexer = left_ax.get_indexer(self.right_join_keys[0])
        else:
            # max groups = largest possible number of distinct groups
            left_key, right_key, max_groups = self._get_group_keys()

            join_func = _join_functions[self.how]
            left_indexer, right_indexer = join_func(left_key.astype('i4'),
                                                    right_key.astype('i4'),
                                                    max_groups)

            if self.right_index:
                join_index = self.left.index.take(left_indexer)
            elif self.left_index:
                join_index = self.right.index.take(right_indexer)
            else:
                join_index = Index(np.arange(len(left_indexer)))

        return join_index, left_indexer, right_indexer

    def _get_merge_data(self, join_names):
        """
        Handles overlapping column names etc.
        """
        ldata, rdata = self.left._data, self.right._data
        lsuf, rsuf = self.suffixes
        exclude_names = [x for x in join_names if x is not None]
        ldata, rdata = ldata._maybe_rename_join(rdata, lsuf, rsuf,
                                                exclude=exclude_names,
                                                copydata=False)

        return ldata, rdata

    def _get_merge_keys(self):
        """
        Note: has side effects (copy/delete key columns)

        Parameters
        ----------
        left
        right
        on

        Returns
        -------
        left_keys, right_keys
        """
        # Hm, any way to make this logic less complicated??
        join_names = []

        if (self.on is None and self.left_on is None
            and self.right_on is None):

            if self.left_index and self.right_index:
                pass
            elif self.left_index:
                if self.right_on is None:
                    raise Exception('Must pass right_on or right_index=True')
            elif self.right_index:
                if self.left_on is None:
                    raise Exception('Must pass left_on or left_index=True')
            else:
                # use the common columns
                common_cols = self.left.columns.intersection(self.right.columns)
                self.left_on = self.right_on = common_cols
                self.drop_keys = True

        elif self.on is not None:
            if self.left_on is not None or self.right_on is not None:
                raise Exception('Can only pass on OR left_on and '
                                'right_on')
            self.left_on = self.right_on = self.on
            self.drop_keys = True

        # this is a touch kludgy, but accomplishes the goal
        left_keys = None
        if self.left_on is not None:
            self.left, left_keys, left_names = \
                _get_keys(self.left, self.left_on, drop=self.drop_keys)
            join_names = left_names

        right_keys = None
        if self.right_on is not None:
            self.right, right_keys, right_names = \
                _get_keys(self.right, self.right_on, drop=self.drop_keys)
            join_names = right_names

        return left_keys, right_keys, join_names

    def _get_group_keys(self):
        """

        Parameters
        ----------

        Returns
        -------

        """
        if self.left_index:
            if isinstance(self.left.index, MultiIndex):
                left_keys = [lev.values.take(lab)
                             for lev, lab in zip(self.left.index.levels,
                                                 self.left.index.labels)]
            else:
                left_keys = [self.left.index.values]
        else:
            left_keys = self.left_join_keys

        if self.right_index:
            if isinstance(self.right.index, MultiIndex):
                right_keys = [lev.values.take(lab)
                              for lev, lab in zip(self.right.index.levels,
                                                  self.right.index.labels)]
            else:
                right_keys = [self.right.index.values]
        else:
            right_keys = self.right_join_keys

        assert(len(left_keys) == len(right_keys))

        left_labels = []
        right_labels = []
        group_sizes = []

        for lk, rk in zip(left_keys, right_keys):
            llab, rlab, count = _factorize_objects(lk, rk, sort=self.sort)

            left_labels.append(llab)
            right_labels.append(rlab)
            group_sizes.append(count)

        left_group_key = get_group_index(left_labels, group_sizes)
        right_group_key = get_group_index(right_labels, group_sizes)

        max_groups = 1L
        for x in group_sizes:
            max_groups *= long(x)

        if max_groups > 2**63:  # pragma: no cover
            raise Exception('Combinatorial explosion! (boom)')

        left_group_key, right_group_key, max_groups = \
            _factorize_int64(left_group_key, right_group_key,
                             sort=self.sort)
        return left_group_key, right_group_key, max_groups

def _get_keys(frame, on, drop=False):
    to_drop = []
    keys = []
    names = []
    for k in on:
        if isinstance(k, np.ndarray) and len(k) == len(frame):
            keys.append(k)
            names.append(None) # super kludge-tastic
        else:
            to_drop.append(k)
            keys.append(frame[k].values)
            names.append(k)


    if drop:
        frame = frame.copy()
        for k in to_drop:
            del frame[k]

        # this is a bit too expensive...
        # frame = frame.drop(to_drop, axis=1)

    return frame, keys, names


def _get_multiindex_indexer(join_keys, index, sort=True):
    shape = []
    labels = []
    for level, key in zip(index.levels, join_keys):
        llab, rlab, count = _factorize_objects(level, key, sort=False)
        labels.append(rlab)
        shape.append(count)

    left_group_key = get_group_index(labels, shape) #.astype('i4')
    right_group_key = get_group_index(index.labels, shape) #.astype('i4')

    left_group_key, right_group_key, max_groups = \
        _factorize_int64(left_group_key, right_group_key,
                         sort=sort)

    left_indexer, right_indexer = \
        lib.left_outer_join(left_group_key.astype('i4'),
                            right_group_key.astype('i4'),
                            max_groups)

    return right_indexer

    # after refactorizing, I don't think reordering is necessary

    # NOW! reorder
    #right_indexer.take(left_indexer.argsort())

def _maybe_make_list(obj):
    if obj is not None and not isinstance(obj, (tuple, list)):
        return [obj]
    return obj

def _right_outer_join(x, y, max_groups):
    right_indexer, left_indexer = lib.left_outer_join(y, x, max_groups)
    return left_indexer, right_indexer

_join_functions = {
    'inner' : lib.inner_join,
    'left' : lib.left_outer_join,
    'right' : _right_outer_join,
    'outer' : lib.full_outer_join,
}

def _factorize_int64(left_index, right_index, sort=True):
    rizer = lib.Int64Factorizer(max(len(left_index), len(right_index)))

    # 32-bit compatibility
    if left_index.dtype != np.int64:  # pragma: no cover
        left_index = left_index.astype('i8')

    if right_index.dtype != np.int64:  # pragma: no cover
        right_index = right_index.astype('i8')

    llab, _ = rizer.factorize(left_index)
    rlab, _ = rizer.factorize(right_index)

    if sort:
        llab, rlab = _sort_labels(np.array(rizer.uniques), llab, rlab)

    return llab, rlab, rizer.get_count()

def _factorize_objects(left_index, right_index, sort=True):
    rizer = lib.Factorizer(max(len(left_index), len(right_index)))

    llab, _ = rizer.factorize(left_index.astype('O'))
    rlab, _ = rizer.factorize(right_index.astype('O'))

    count = rizer.get_count()

    if sort:
        llab, rlab = _sort_labels(rizer.uniques, llab, rlab)

        # TODO: na handling

    return llab, rlab, count

def _sort_labels(uniques, left, right):
    if not isinstance(uniques, np.ndarray):
        # tuplesafe
        uniques = Index(uniques).values

    sorter = uniques.argsort()

    reverse_indexer = np.empty(len(sorter), dtype=np.int32)
    reverse_indexer.put(sorter, np.arange(len(sorter)))
    return reverse_indexer.take(left), reverse_indexer.take(right)


class _BlockJoinOperation(object):
    """
    BlockJoinOperation made generic for N DataFrames

    Object responsible for orchestrating efficient join operation between two
    BlockManager data structures
    """
    def __init__(self, data_list, join_index, indexers, axis=1, copy=True):
        if axis <= 0:  # pragma: no cover
            raise Exception('Only axis >= 1 supported for this operation')

        assert(len(data_list) == len(indexers))

        self.units = []
        for data, indexer in zip(data_list, indexers):
            if not data.is_consolidated():
                data = data.consolidate()
            self.units.append(_JoinUnit(data.blocks, indexer))

        self.join_index = join_index
        self.axis = axis
        self.copy = copy

        # do NOT sort
        self.result_items = _concat_indexes([d.items for d in data_list])
        self.result_axes = list(data_list[0].axes)
        self.result_axes[0] = self.result_items
        self.result_axes[axis] = self.join_index

    def _prepare_blocks(self):
        blockmaps = []

        for unit in self.units:
            join_blocks = unit.get_upcasted_blocks()
            type_map = dict((type(blk), blk) for blk in join_blocks)
            blockmaps.append(type_map)

        return blockmaps

    def get_result(self):
        """
        Returns
        -------
        merged : BlockManager
        """
        blockmaps = self._prepare_blocks()
        kinds = _get_all_block_kinds(blockmaps)

        result_blocks = []

        # maybe want to enable flexible copying <-- what did I mean?
        for klass in kinds:
            klass_blocks = [mapping.get(klass) for mapping in blockmaps]
            res_blk = self._get_merged_block(klass_blocks)
            result_blocks.append(res_blk)

        return BlockManager(result_blocks, self.result_axes)

    def _get_merged_block(self, blocks):

        to_merge = []

        for unit, block in zip(self.units, blocks):
            if block is not None:
                to_merge.append((unit, block))

        if len(to_merge) > 1:
            return self._merge_blocks(to_merge)
        else:
            unit, block = to_merge[0]
            return unit.reindex_block(block, self.axis,
                                      self.result_items, copy=self.copy)

    def _merge_blocks(self, merge_chunks):
        """
        merge_chunks -> [(_JoinUnit, Block)]
        """
        funit, fblock = merge_chunks[0]
        fidx = funit.indexer

        out_shape = list(fblock.values.shape)

        n = len(fidx) if fidx is not None else out_shape[self.axis]

        out_shape[0] = sum(len(blk) for unit, blk in merge_chunks)
        out_shape[self.axis] = n

        # Should use Fortran order??
        out = np.empty(out_shape, dtype=fblock.values.dtype)

        sofar = 0
        for unit, blk in merge_chunks:
            out_chunk = out[sofar : sofar + len(blk)]

            if unit.indexer is None:
            # is this really faster than assigning to arr.flat?
                com.take_fast(blk.values, np.arange(n, dtype='i4'),
                              None, False,
                              axis=self.axis, out=out_chunk)
            else:
                # write out the values to the result array
                com.take_fast(blk.values, unit.indexer,
                              None, False,
                              axis=self.axis, out=out_chunk)

            sofar += len(blk)

        # does not sort
        new_block_items = _concat_indexes([b.items for _, b in merge_chunks])
        return make_block(out, new_block_items, self.result_items)



class _JoinUnit(object):
    """
    Blocks plus indexer
    """

    def __init__(self, blocks, indexer):
        self.blocks = blocks
        self.indexer = indexer

    @cache_readonly
    def mask_info(self):
        if self.indexer is None or not _may_need_upcasting(self.blocks):
            mask = None
            need_masking = False
        else:
            mask = self.indexer == -1
            need_masking = mask.any()

        return mask, need_masking

    @property
    def need_masking(self):
        return self.mask_info[1]

    def get_upcasted_blocks(self):
        # will short-circuit and not compute lneed_masking if indexer is None
        if self.need_masking:
            return _upcast_blocks(self.blocks)
        return self.blocks

    def reindex_block(self, block, axis, ref_items, copy=True):
        # still some inefficiency here for bool/int64 because in the case where
        # no masking is needed, take_fast will recompute the mask

        mask, need_masking = self.mask_info

        if self.indexer is None:
            if copy:
                result = block.copy()
            else:
                result = block
        else:
            result = block.reindex_axis(self.indexer, mask, need_masking,
                                        axis=axis)

        result.ref_items = ref_items
        return result

def _may_need_upcasting(blocks):
    for block in blocks:
        if isinstance(block, (IntBlock, BoolBlock)):
            return True
    return False


def _upcast_blocks(blocks):
    """
    Upcast and consolidate if necessary
    """
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

def _get_all_block_kinds(blockmaps):
    kinds = set()
    for mapping in blockmaps:
        kinds |= set(mapping)
    return kinds

#----------------------------------------------------------------------
# Concatenate DataFrame objects


class Concatenator(object):
    """
    Orchestrates a concatenation operation for BlockManagers, with little hacks
    to support sparse data structures, etc.
    """

    def __init__(self, objs, axis=0, join='outer', join_axes=None,
                 ignore_index=False, verify_integrity=False):
        if join == 'outer':
            self.intersect = False
        elif join == 'inner':
            self.intersect = True
        else:  # pragma: no cover
            raise ValueError('Only can inner (intersect) or outer (union) join '
                             'the other axis')

        # consolidate data
        for obj in objs:
            obj.consolidate(inplace=True)
        self.objs = objs

        # Need to flip BlockManager axis in the DataFrame special case
        if isinstance(objs[0], DataFrame):
            axis = 1 if axis == 0 else 0

        # note: this is the BlockManager axis (since DataFrame is transposed)
        self.axis = axis

        self.join_axes = join_axes

        self.ignore_index = ignore_index
        self.verify_integrity = verify_integrity

        self.new_axes = self._get_new_axes()

    def get_result(self):
        first = self.objs[0]

        if len(self.objs) == 1:
            return first

        new_data = self._get_concatenated_data()
        return first._from_axes(new_data, self.new_axes)

    def _get_concatenated_data(self):
        try:
            # need to conform to same other (joined) axes for block join
            reindexed_data = self._get_reindexed_data()

            blockmaps = []
            for data in reindexed_data:
                type_map = dict((type(blk), blk) for blk in data.blocks)
                blockmaps.append(type_map)
            kinds = _get_all_block_kinds(blockmaps)

            new_blocks = []
            for kind in kinds:
                klass_blocks = [mapping.get(kind) for mapping in blockmaps]
                stacked_block = self._concat_blocks(klass_blocks)
                new_blocks.append(stacked_block)
            new_data = BlockManager(new_blocks, self.new_axes)
        except Exception:  # EAFP
            # should not be possible to fail here for the expected reason with
            # axis = 0
            if self.axis == 0:  # pragma: no cover
                raise

            new_data = {}
            for item in self.new_axes[0]:
                new_data[item] = self._concat_single_item(item)

        return new_data

    def _get_reindexed_data(self):
        # HACK: ugh

        reindexed_data = []
        if isinstance(self.objs[0], SparseDataFrame):
            pass
        else:
            axes_to_reindex = list(enumerate(self.new_axes))
            axes_to_reindex.pop(self.axis)

            for obj in self.objs:
                data = obj._data
                for i, ax in axes_to_reindex:
                    data = data.reindex_axis(ax, axis=i, copy=False)
                reindexed_data.append(data)

        return reindexed_data

    def _concat_blocks(self, blocks):
        concat_values = np.concatenate([b.values for b in blocks],
                                       axis=self.axis)

        if self.axis > 0:
            # Not safe to remove this check, need to profile
            if not _all_indexes_same([b.items for b in blocks]):
                raise Exception('dtypes are not consistent throughout '
                                'DataFrames')
            return make_block(concat_values, blocks[0].items, self.new_axes[0])
        else:
            concat_items = _concat_indexes([b.items for b in blocks])
            # TODO: maybe want to "take" from the new columns?
            return make_block(concat_values, concat_items, self.new_axes[0])

    def _concat_single_item(self, item):
        all_values = []
        dtypes = set()
        for obj in self.objs:
            try:
                values = obj._data.get(item)
                dtypes.add(values.dtype)
                all_values.append(values)
            except KeyError:
                all_values.append(None)

        # this stinks
        have_object = False
        for dtype in dtypes:
            if issubclass(dtype.type, (np.object_, np.bool_)):
                have_object = True
        if have_object:
            empty_dtype = np.object_
        else:
            empty_dtype = np.float64

        to_concat = []
        for df, col_values in zip(self.objs, all_values):
            if col_values is None:
                missing_arr = np.empty(len(df), dtype=empty_dtype)
                missing_arr.fill(np.nan)
                to_concat.append(missing_arr)
            else:
                to_concat.append(col_values)

        return np.concatenate(to_concat)

    def _get_new_axes(self):
        ndim = self.objs[0].ndim
        new_axes = [None] * ndim

        if self.ignore_index:
            concat_axis = None
        else:
            concat_axis = _concat_indexes([x._data.axes[self.axis]
                                           for x in self.objs])
            self._maybe_check_integrity(concat_axis)

        new_axes[self.axis] = concat_axis

        if self.join_axes is None:
            for i in range(ndim):
                if i == self.axis:
                    continue
                all_indexes = [x._data.axes[i] for x in self.objs]
                comb_axis = _get_combined_index(all_indexes,
                                                intersect=self.intersect)
                new_axes[i] = comb_axis

        else:
            assert(len(self.join_axes) == ndim - 1)

            # ufff...
            indices = range(ndim)
            indices.remove(self.axis)

            for i, ax in zip(indices, self.join_axes):
                new_axes[i] = ax

        return new_axes

    def _maybe_check_integrity(self, concat_index):
        if self.verify_integrity:
            if not concat_index._verify_integrity():
                overlap = concat_index.get_duplicates()
                raise Exception('Indexes have overlapping values: %s'
                                % str(overlap))


def _concat_frames_hierarchical(frames, keys, groupings, axis=0):
    names = [ping.name for ping in groupings]
    levels = [ping.group_index for ping in groupings]

    if axis == 0:
        indexes = [x.index for x in frames]
        new_index = _make_concat_multiindex(indexes, keys, levels, names)
        new_columns = frames[0].columns
    else:
        all_columns = [x.columns for x in frames]
        new_columns = _make_concat_multiindex(all_columns, keys,
                                              levels, names)
        new_index = frames[0].index

    if frames[0]._is_mixed_type:
        new_data = {}
        for col in new_columns:
            new_data[col] = np.concatenate([x[col].values for x in frames])
        return DataFrame(new_data, index=new_index, columns=new_columns)
    else:
        new_values = np.concatenate([x.values for x in frames], axis=axis)
        return DataFrame(new_values, index=new_index, columns=new_columns)

def _concat_indexes(indexes):
    return indexes[0].append(indexes[1:])

def _make_concat_multiindex(indexes, keys, levels, names):
    single_level = len(levels) == 1

    if not _all_indexes_same(indexes):
        label_list = []

        # things are potentially different sizes, so compute the exact labels
        # for each level and pass those to MultiIndex.from_arrays
        if single_level:
            zipped = [keys]
        else:
            zipped = zip(*keys)

        for hlevel in zipped:
            to_concat = []
            for k, index in zip(hlevel, indexes):
                to_concat.append(np.repeat(k, len(index)))
            label_list.append(np.concatenate(to_concat))

        concat_index = _concat_indexes(indexes)

        # these go at the end
        if isinstance(concat_index, MultiIndex):
            for level in range(concat_index.nlevels):
                label_list.append(concat_index.get_level_values(level))
        else:
            label_list.append(concat_index.values)

        consensus_name = indexes[0].names
        for index in indexes[1:]:
            if index.names != consensus_name:
                consensus_name = [None] * index.nlevels
                break
        names.extend(consensus_name)

        return MultiIndex.from_arrays(label_list, names=names)

    new_index = indexes[0]
    n = len(new_index)

    names.append(indexes[0].name)

    new_levels = list(levels)

    # do something a bit more speedy
    new_levels.append(new_index)

    # construct labels
    labels = []

    if single_level:
        zipped = [keys]
    else:
        zipped = zip(*keys)

    for hlevel, level in zip(zipped, levels):
        mapped = level.get_indexer(hlevel)
        labels.append(np.repeat(mapped, n))

    # last labels for the new level
    labels.append(np.tile(np.arange(n), len(indexes)))
    return MultiIndex(levels=new_levels, labels=labels, names=names)

def _all_indexes_same(indexes):
    first = indexes[0]
    for index in indexes[1:]:
        if not first.equals(index):
            return False
    return True

