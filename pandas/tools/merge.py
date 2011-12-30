"""
SQL-style merge routines
"""

import numpy as np

from pandas.core.frame import DataFrame
from pandas.core.index import Index
from pandas.core.internals import (IntBlock, BoolBlock, BlockManager,
                                   make_block, _consolidate)
from pandas.util.decorators import cache_readonly
import pandas.core.common as com

import pandas._tseries as lib

def merge(left, right, how='left', on=None, left_on=None, right_on=None,
          left_index=False, right_index=False, sort=True,
          suffixes=('.x', '.y'), copy=True):
    """
    Merge DataFrame objects by performing a database-style join operation by
    columns or indexes

    Parameters
    ----------
    left : DataFrame
    right : DataFrame
    how : {'left', 'right', 'outer', 'inner'}, default 'left'
        * left: use only keys from left frame (SQL: left outer join)
        * right: use only keys from right frame (SQL: right outer join)
        * outer: use union of keys from both frames (SQL: full outer join)
        * inner: use intersection of keys from both frames (SQL: inner join)
    on : label or list
        Field names to join on. Must be found in both DataFrames
    left_on : label or list
        Field names to join on in left DataFrame
    right_on : label or list
        Field names to join on in right DataFrame
    left_index : boolean, default True
        Use the index from the left DataFrame as the join key
    right_index : boolean, default True
        Use the index from the right DataFrame as the join key
    sort : boolean, default True
        Sort the join keys lexicographically in the result DataFrame
    suffixes : 2-length sequence (tuple, list, ...)
        Suffix to apply to overlapping column names in the left and right
        side, respectively
    copy : boolean, default True
        If False, do not copy data unnecessarily

    Returns
    -------
    merged : DataFrame
    """
    op = _MergeOperation(left, right, how=how, on=on, left_on=left_on,
                         right_on=right_on, left_index=left_index,
                         right_index=right_index, sort=sort, suffixes=suffixes,
                         copy=copy)
    return op.get_result()


# TODO: shortcuts with MultiIndex labels already computed
# TODO: NA group handling
# TODO: DONE group column names in result
# TODO: transformations??
# TODO: only copy DataFrames when modification necessary

def join_managers(left, right, axis=1, how='left', copy=True):
    join_index, left_indexer, right_indexer = \
        left.axes[axis].join(right.axes[axis], how=how, return_indexers=True)
    op = _JoinOperation(left, right, join_index, left_indexer,
                        right_indexer, axis=axis)
    return op.get_result(copy=copy)

class _MergeOperation(object):

    def __init__(self, left, right, how='inner', on=None,
                 left_on=None, right_on=None,
                 left_index=False, right_index=False, sort=True,
                 suffixes=('.x', '.y'), copy=True):
        self.left = self.orig_left = left
        self.right = self.orig_right = right
        self.how = how

        self.on = _maybe_make_list(on)
        self.left_on = _maybe_make_list(left_on)
        self.right_on = _maybe_make_list(right_on)

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
        left_indexer, right_indexer = self._get_join_indexers()
        new_axis = self._get_new_axis(left_indexer)

        # this is a bit kludgy
        ldata, rdata = self._get_merge_data(self.join_names)

        # TODO: more efficiently handle group keys to avoid extra consolidation!
        join_op = _JoinOperation(ldata, rdata, new_axis,
                                 left_indexer, right_indexer, axis=1)

        result_data = join_op.get_result(copy=self.copy)
        result = DataFrame(result_data)

        self._maybe_add_join_keys(result, left_indexer, right_indexer)

        return result

    def _maybe_add_join_keys(self, result, left_indexer, right_indexer):
        # insert group keys
        for i, name in enumerate(self.join_names):
            # a faster way?
            key_col = com.take_1d(self.left_join_keys[i], left_indexer)
            na_indexer = (left_indexer == -1).nonzero()[0]
            right_na_indexer = right_indexer.take(na_indexer)
            key_col.put(na_indexer, com.take_1d(self.right_join_keys[i],
                                                right_na_indexer))
            result.insert(i, name, key_col)

    def _get_join_indexers(self):
        # max groups = largest possible number of distinct groups
        left_key, right_key, max_groups = \
            _get_group_keys(self.left_join_keys, self.right_join_keys,
                            sort=self.sort)

        join_func = _join_functions[self.how]
        left_indexer, right_indexer = join_func(left_key.astype('i4'),
                                                right_key.astype('i4'),
                                                max_groups)

        return left_indexer, right_indexer

    def _get_new_axis(self, left_indexer):
        if left_indexer is None:
            new_axis = self.left.index
        else:
            new_axis = Index(np.arange(len(left_indexer)))
        return new_axis

    def _get_merge_data(self, join_names):
        """
        Handles overlapping column names etc.
        """
        ldata, rdata = self.left._data, self.right._data
        lsuf, rsuf = self.suffixes

        # basically by construction the column names are stored in
        # left_on...for now
        ldata, rdata = ldata._maybe_rename_join(rdata, lsuf, rsuf,
                                                exclude=join_names,
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
        left_keys = []
        right_keys = []
        join_names = []

        # need_set_names = False
        # pop_right = False

        if (self.on is None and self.left_on is None
            and self.right_on is None):

            if self.left_index and self.right_index:
                left_keys.append(self.left.index.values)
                right_keys.append(self.right.index.values)

                # need_set_names = True

                # XXX something better than this
                join_names.append('join_key')
            elif self.left_index:
                left_keys.append(self.left.index.values)
                if self.right_on is None:
                    raise Exception('Must pass right_on or right_index=True')
            elif self.right_index:
                right_keys.append(self.right.index.values)
                if self.left_on is None:
                    raise Exception('Must pass left_on or left_index=True')
            else:
                # use the common columns
                common_cols = self.left.columns.intersection(self.right.columns)
                self.left_on = self.right_on = common_cols

                # pop_right = True

        elif self.on is not None:
            if self.left_on is not None or self.right_on is not None:
                raise Exception('Can only pass on OR left_on and '
                                'right_on')
            self.left_on = self.right_on = self.on

            # pop_right = True

        # this is a touch kludgy, but accomplishes the goal
        if self.right_on is not None:
            right = self.right.copy()
            right_keys.extend([right.pop(k) for k in self.right_on])
            self.right = right

        if self.left_on is not None:
            left = self.left.copy()
            left_keys.extend([left.pop(k) for k in self.left_on])
            self.left = left

            # TODO: something else?
            join_names = self.left_on

        return left_keys, right_keys, join_names

def _get_group_keys(left_keys, right_keys, sort=True):
    """

    Parameters
    ----------

    Returns
    -------

    """
    from pandas.core.groupby import get_group_index

    assert(len(left_keys) == len(right_keys))

    left_labels = []
    right_labels = []
    group_sizes = []

    for lk, rk in zip(left_keys, right_keys):
        llab, rlab, count = _factorize_objects(lk, rk, sort=sort)

        left_labels.append(llab)
        right_labels.append(rlab)
        group_sizes.append(count)

    left_group_key = get_group_index(left_labels, group_sizes)
    right_group_key = get_group_index(right_labels, group_sizes)
    max_groups = np.prod(group_sizes)

    if max_groups > 1000000:
        # compress
        left_group_key, right_group_key, max_groups = \
            _factorize_int64(left_group_key, right_group_key,
                             sort=sort)

    return left_group_key, right_group_key, max_groups

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


class _JoinOperation(object):
    """
    Object responsible for orchestrating efficient join operation between two
    BlockManager data structures
    """
    def __init__(self, left, right, join_index, left_indexer, right_indexer,
                 axis=1):
        assert(axis > 0)

        if not left.is_consolidated():
            left = left.consolidate()
        if not right.is_consolidated():
            right = right.consolidate()

        self.left = left
        self.right = right
        self.axis = axis

        self.join_index = join_index
        self.lindexer = left_indexer
        self.rindexer = right_indexer

        # do NOT sort
        self.result_items = left.items.append(right.items)
        self.result_axes = list(left.axes)
        self.result_axes[0] = self.result_items
        self.result_axes[axis] = self.join_index

    def get_result(self, copy=False):
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
        left_blockmap, right_blockmap = self._prepare_blocks()

        result_blocks = []

        # maybe want to enable flexible copying

        kinds = set(left_blockmap) | set(right_blockmap)
        for klass in kinds:
            lblk = left_blockmap.get(klass)
            rblk = right_blockmap.get(klass)

            if lblk and rblk:
                # true merge, do not produce intermediate copy
                res_blk = self._merge_blocks(lblk, rblk)
            elif lblk:
                res_blk = self._reindex_block(lblk, side='left')
            else:
                res_blk = self._reindex_block(rblk, side='right')

            result_blocks.append(res_blk)

        return BlockManager(result_blocks, self.result_axes)

    def _prepare_blocks(self):
        lblocks = self.left.blocks
        rblocks = self.right.blocks

        # will short-circuit and not compute lneed_masking
        if self.lneed_masking:
            lblocks = self._upcast_blocks(lblocks)

        if self.rneed_masking:
            rblocks = self._upcast_blocks(rblocks)

        left_blockmap = dict((type(blk), blk) for blk in lblocks)
        right_blockmap = dict((type(blk), blk) for blk in rblocks)

        return left_blockmap, right_blockmap

    def _reindex_block(self, block, side='left', copy=True):
        if side == 'left':
            indexer = self.lindexer
            mask, need_masking = self.lmask_info
        else:
            indexer = self.rindexer
            mask, need_masking = self.rmask_info

        # still some inefficiency here for bool/int64 because in the case where
        # no masking is needed, take_fast will recompute the mask

        if indexer is None and copy:
            result = block.copy()
        else:
            result = block.reindex_axis(indexer, mask, need_masking,
                                        axis=self.axis)

        result.ref_items = self.result_items
        return result

    @cache_readonly
    def lmask_info(self):
        if (self.lindexer is None or
            not self._may_need_upcasting(self.left.blocks)):
            lmask = None
            lneed_masking = False
        else:
            lmask = self.lindexer == -1
            lneed_masking = lmask.any()

        return lmask, lneed_masking

    @cache_readonly
    def rmask_info(self):
        if (self.rindexer is None or
            not self._may_need_upcasting(self.right.blocks)):
            rmask = None
            rneed_masking = False
        else:
            rmask = self.rindexer == -1
            rneed_masking = rmask.any()

        return rmask, rneed_masking

    @property
    def lneed_masking(self):
        return self.lmask_info[1]

    @property
    def rneed_masking(self):
        return self.rmask_info[1]

    @staticmethod
    def _may_need_upcasting(blocks):
        for block in blocks:
            if isinstance(block, (IntBlock, BoolBlock)):
                return True
        return False

    def _merge_blocks(self, lblk, rblk):
        lidx = self.lindexer
        ridx = self.rindexer

        n = lblk.values.shape[self.axis] if lidx is None else len(lidx)
        lk = len(lblk.items)
        rk = len(rblk.items)

        out_shape = list(lblk.shape)
        out_shape[0] = lk + rk
        out_shape[self.axis] = n

        out = np.empty(out_shape, dtype=lblk.values.dtype)

        # is this really faster than assigning to arr.flat?
        if lidx is None:
            # out[:lk] = lblk.values
            com.take_fast(lblk.values, np.arange(n, dtype='i4'),
                          None, False,
                          axis=self.axis, out=out[:lk])
        else:
            # write out the values to the result array
            com.take_fast(lblk.values, lidx, None, False,
                             axis=self.axis, out=out[:lk])
        if ridx is None:
            # out[lk:] = lblk.values
            com.take_fast(rblk.values, np.arange(n, dtype='i4'),
                          None, False,
                          axis=self.axis, out=out[lk:])
        else:
            com.take_fast(rblk.values, ridx, None, False,
                          axis=self.axis, out=out[lk:])

        # does not sort
        new_items = lblk.items.append(rblk.items)
        return make_block(out, new_items, self.result_items)

    @staticmethod
    def _upcast_blocks(blocks):
        """
        Upcast and consolidate if necessary
        """
        # if not need_masking:
        #     return blocks

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
