"""
SQL-style merge routines
"""
import types

import numpy as np
from pandas.compat import range, long, lrange, lzip, zip
import pandas.compat as compat
from pandas.core.categorical import Categorical
from pandas.core.frame import DataFrame, _merge_doc
from pandas.core.generic import NDFrame
from pandas.core.groupby import get_group_index
from pandas.core.series import Series
from pandas.core.index import (Index, MultiIndex, _get_combined_index,
                               _ensure_index, _get_consensus_names,
                               _all_indexes_same)
from pandas.core.internals import (TimeDeltaBlock, IntBlock, BoolBlock, BlockManager,
                                   make_block, _consolidate)
from pandas.util.decorators import cache_readonly, Appender, Substitution
from pandas.core.common import (PandasError, ABCSeries,
                                is_timedelta64_dtype, is_datetime64_dtype,
                                is_integer_dtype, isnull)

import pandas.core.common as com

import pandas.lib as lib
import pandas.algos as algos
import pandas.hashtable as _hash
import pandas.tslib as tslib

@Substitution('\nleft : DataFrame')
@Appender(_merge_doc, indents=0)
def merge(left, right, how='inner', on=None, left_on=None, right_on=None,
          left_index=False, right_index=False, sort=False,
          suffixes=('_x', '_y'), copy=True):
    op = _MergeOperation(left, right, how=how, on=on, left_on=left_on,
                         right_on=right_on, left_index=left_index,
                         right_index=right_index, sort=sort, suffixes=suffixes,
                         copy=copy)
    return op.get_result()
if __debug__:
    merge.__doc__ = _merge_doc % '\nleft : DataFrame'


class MergeError(Exception):
    pass


def ordered_merge(left, right, on=None, left_by=None, right_by=None,
                  left_on=None, right_on=None,
                  fill_method=None, suffixes=('_x', '_y')):
    """Perform merge with optional filling/interpolation designed for ordered
    data like time series data. Optionally perform group-wise merge (see
    examples)

    Parameters
    ----------
    left : DataFrame
    right : DataFrame
    fill_method : {'ffill', None}, default None
        Interpolation method for data
    on : label or list
        Field names to join on. Must be found in both DataFrames.
    left_on : label or list, or array-like
        Field names to join on in left DataFrame. Can be a vector or list of
        vectors of the length of the DataFrame to use a particular vector as
        the join key instead of columns
    right_on : label or list, or array-like
        Field names to join on in right DataFrame or vector/list of vectors per
        left_on docs
    left_by : column name or list of column names
        Group left DataFrame by group columns and merge piece by piece with
        right DataFrame
    right_by : column name or list of column names
        Group right DataFrame by group columns and merge piece by piece with
        left DataFrame
    suffixes : 2-length sequence (tuple, list, ...)
        Suffix to apply to overlapping column names in the left and right
        side, respectively

    Examples
    --------
    >>> A                      >>> B
          key  lvalue group        key  rvalue
    0   a       1     a        0     b       1
    1   c       2     a        1     c       2
    2   e       3     a        2     d       3
    3   a       1     b
    4   c       2     b
    5   e       3     b

    >>> ordered_merge(A, B, fill_method='ffill', left_by='group')
       key  lvalue group  rvalue
    0    a       1     a     NaN
    1    b       1     a       1
    2    c       2     a       2
    3    d       2     a       3
    4    e       3     a       3
    5    f       3     a       4
    6    a       1     b     NaN
    7    b       1     b       1
    8    c       2     b       2
    9    d       2     b       3
    10   e       3     b       3
    11   f       3     b       4

    Returns
    -------
    merged : DataFrame
    """
    def _merger(x, y):
        op = _OrderedMerge(x, y, on=on, left_on=left_on, right_on=right_on,
                           # left_index=left_index, right_index=right_index,
                           suffixes=suffixes, fill_method=fill_method)
        return op.get_result()

    if left_by is not None and right_by is not None:
        raise ValueError('Can only group either left or right frames')
    elif left_by is not None:
        if not isinstance(left_by, (list, tuple)):
            left_by = [left_by]
        pieces = []
        for key, xpiece in left.groupby(left_by):
            merged = _merger(xpiece, right)
            for k in left_by:
                # May have passed ndarray
                try:
                    if k in merged:
                        merged[k] = key
                except:
                    pass
            pieces.append(merged)
        return concat(pieces, ignore_index=True)
    elif right_by is not None:
        if not isinstance(right_by, (list, tuple)):
            right_by = [right_by]
        pieces = []
        for key, ypiece in right.groupby(right_by):
            merged = _merger(left, ypiece)
            for k in right_by:
                try:
                    if k in merged:
                        merged[k] = key
                except:
                    pass
            pieces.append(merged)
        return concat(pieces, ignore_index=True)
    else:
        return _merger(left, right)


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
                 suffixes=('_x', '_y'), copy=True):
        self.left = self.orig_left = left
        self.right = self.orig_right = right
        self.how = how
        self.axis = axis

        self.on = com._maybe_make_list(on)
        self.left_on = com._maybe_make_list(left_on)
        self.right_on = com._maybe_make_list(right_on)

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
        ldata, rdata = self._get_merge_data()

        # TODO: more efficiently handle group keys to avoid extra
        #       consolidation!
        join_op = _BlockJoinOperation([ldata, rdata], join_index,
                                      [left_indexer, right_indexer], axis=1,
                                      copy=self.copy)

        result_data = join_op.get_result()
        result = DataFrame(result_data)

        self._maybe_add_join_keys(result, left_indexer, right_indexer)

        return result

    def _maybe_add_join_keys(self, result, left_indexer, right_indexer):
        # insert group keys

        keys = zip(self.join_names, self.left_on, self.right_on)
        for i, (name, lname, rname) in enumerate(keys):
            if not _should_fill(lname, rname):
                continue

            if name in result:
                key_col = result[name]

                if left_indexer is not None and right_indexer is not None:

                    if name in self.left:
                        na_indexer = (left_indexer == -1).nonzero()[0]
                        if len(na_indexer) == 0:
                            continue

                        right_na_indexer = right_indexer.take(na_indexer)
                        key_col.put(
                            na_indexer, com.take_1d(self.right_join_keys[i],
                                                    right_na_indexer))
                    elif name in self.right:
                        na_indexer = (right_indexer == -1).nonzero()[0]
                        if len(na_indexer) == 0:
                            continue

                        left_na_indexer = left_indexer.take(na_indexer)
                        key_col.put(na_indexer, com.take_1d(self.left_join_keys[i],
                                                            left_na_indexer))

            elif left_indexer is not None:
                if name is None:
                    name = 'key_%d' % i

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
            join_index, left_indexer, right_indexer = \
                _left_join_on_index(left_ax, right_ax, self.left_join_keys,
                                    sort=self.sort)

        elif self.left_index and self.how == 'right':
            join_index, right_indexer, left_indexer = \
                _left_join_on_index(right_ax, left_ax, self.right_join_keys,
                                    sort=self.sort)
        else:
            (left_indexer,
             right_indexer) = _get_join_indexers(self.left_join_keys,
                                                 self.right_join_keys,
                                                 sort=self.sort, how=self.how)

            if self.right_index:
                join_index = self.left.index.take(left_indexer)
            elif self.left_index:
                join_index = self.right.index.take(right_indexer)
            else:
                join_index = Index(np.arange(len(left_indexer)))

        return join_index, left_indexer, right_indexer

    def _get_merge_data(self):
        """
        Handles overlapping column names etc.
        """
        ldata, rdata = self.left._data, self.right._data
        lsuf, rsuf = self.suffixes
        ldata, rdata = ldata._maybe_rename_join(rdata, lsuf, rsuf,
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
        self._validate_specification()

        left_keys = []
        right_keys = []
        join_names = []
        right_drop = []
        left_drop = []
        left, right = self.left, self.right

        is_lkey = lambda x: isinstance(x, (np.ndarray, ABCSeries)) and len(x) == len(left)
        is_rkey = lambda x: isinstance(x, (np.ndarray, ABCSeries)) and len(x) == len(right)

        # ugh, spaghetti re #733
        if _any(self.left_on) and _any(self.right_on):
            for lk, rk in zip(self.left_on, self.right_on):
                if is_lkey(lk):
                    left_keys.append(lk)
                    if is_rkey(rk):
                        right_keys.append(rk)
                        join_names.append(None)  # what to do?
                    else:
                        right_keys.append(right[rk].values)
                        join_names.append(rk)
                else:
                    if not is_rkey(rk):
                        right_keys.append(right[rk].values)
                        if lk == rk:
                            # avoid key upcast in corner case (length-0)
                            if len(left) > 0:
                                right_drop.append(rk)
                            else:
                                left_drop.append(lk)
                    else:
                        right_keys.append(rk)
                    left_keys.append(left[lk].values)
                    join_names.append(lk)
        elif _any(self.left_on):
            for k in self.left_on:
                if is_lkey(k):
                    left_keys.append(k)
                    join_names.append(None)
                else:
                    left_keys.append(left[k].values)
                    join_names.append(k)
            if isinstance(self.right.index, MultiIndex):
                right_keys = [lev.values.take(lab)
                              for lev, lab in zip(self.right.index.levels,
                                                  self.right.index.labels)]
            else:
                right_keys = [self.right.index.values]
        elif _any(self.right_on):
            for k in self.right_on:
                if is_rkey(k):
                    right_keys.append(k)
                    join_names.append(None)
                else:
                    right_keys.append(right[k].values)
                    join_names.append(k)
            if isinstance(self.left.index, MultiIndex):
                left_keys = [lev.values.take(lab)
                             for lev, lab in zip(self.left.index.levels,
                                                 self.left.index.labels)]
            else:
                left_keys = [self.left.index.values]

        if left_drop:
            self.left = self.left.drop(left_drop, axis=1)

        if right_drop:
            self.right = self.right.drop(right_drop, axis=1)

        return left_keys, right_keys, join_names

    def _validate_specification(self):
        # Hm, any way to make this logic less complicated??
        if (self.on is None and self.left_on is None
                and self.right_on is None):

            if self.left_index and self.right_index:
                self.left_on, self.right_on = (), ()
            elif self.left_index:
                if self.right_on is None:
                    raise MergeError('Must pass right_on or right_index=True')
            elif self.right_index:
                if self.left_on is None:
                    raise MergeError('Must pass left_on or left_index=True')
            else:
                if not self.left.columns.is_unique:
                    raise MergeError("Left data columns not unique: %s"
                                     % repr(self.left.columns))

                if not self.right.columns.is_unique:
                    raise MergeError("Right data columns not unique: %s"
                                     % repr(self.right.columns))

                # use the common columns
                common_cols = self.left.columns.intersection(
                    self.right.columns)
                if len(common_cols) == 0:
                    raise MergeError('No common columns to perform merge on')
                self.left_on = self.right_on = common_cols
        elif self.on is not None:
            if self.left_on is not None or self.right_on is not None:
                raise MergeError('Can only pass on OR left_on and '
                                 'right_on')
            self.left_on = self.right_on = self.on
        elif self.left_on is not None:
            n = len(self.left_on)
            if self.right_index:
                if len(self.left_on) != self.right.index.nlevels:
                    raise ValueError('len(left_on) must equal the number '
                                         'of levels in the index of "right"')
                self.right_on = [None] * n
        elif self.right_on is not None:
            n = len(self.right_on)
            if self.left_index:
                if len(self.right_on) != self.left.index.nlevels:
                    raise ValueError('len(right_on) must equal the number '
                                         'of levels in the index of "left"')
                self.left_on = [None] * n
        if len(self.right_on) != len(self.left_on):
            raise ValueError("len(right_on) must equal len(left_on)")


def _get_join_indexers(left_keys, right_keys, sort=False, how='inner'):
    """

    Parameters
    ----------

    Returns
    -------

    """
    if len(left_keys) != len(right_keys):
        raise AssertionError('left_key and right_keys must be the same length')

    left_labels = []
    right_labels = []
    group_sizes = []

    for lk, rk in zip(left_keys, right_keys):
        llab, rlab, count = _factorize_keys(lk, rk, sort=sort)

        left_labels.append(llab)
        right_labels.append(rlab)
        group_sizes.append(count)

    max_groups = long(1)
    for x in group_sizes:
        max_groups *= long(x)

    if max_groups > 2 ** 63:  # pragma: no cover
        left_group_key, right_group_key, max_groups = \
            _factorize_keys(lib.fast_zip(left_labels),
                            lib.fast_zip(right_labels))
    else:
        left_group_key = get_group_index(left_labels, group_sizes)
        right_group_key = get_group_index(right_labels, group_sizes)

        left_group_key, right_group_key, max_groups = \
            _factorize_keys(left_group_key, right_group_key, sort=sort)

    join_func = _join_functions[how]
    return join_func(left_group_key, right_group_key, max_groups)


class _OrderedMerge(_MergeOperation):

    def __init__(self, left, right, on=None, by=None, left_on=None,
                 right_on=None, axis=1, left_index=False, right_index=False,
                 suffixes=('_x', '_y'), copy=True,
                 fill_method=None):

        self.fill_method = fill_method

        _MergeOperation.__init__(self, left, right, on=on, left_on=left_on,
                                 right_on=right_on, axis=axis,
                                 left_index=left_index,
                                 right_index=right_index,
                                 how='outer', suffixes=suffixes,
                                 sort=True  # sorts when factorizing
                                 )

    def get_result(self):
        join_index, left_indexer, right_indexer = self._get_join_info()

        # this is a bit kludgy
        ldata, rdata = self._get_merge_data()

        if self.fill_method == 'ffill':
            left_join_indexer = algos.ffill_indexer(left_indexer)
            right_join_indexer = algos.ffill_indexer(right_indexer)
        else:
            left_join_indexer = left_indexer
            right_join_indexer = right_indexer

        join_op = _BlockJoinOperation([ldata, rdata], join_index,
                                      [left_join_indexer, right_join_indexer],
                                      axis=1, copy=self.copy)

        result_data = join_op.get_result()
        result = DataFrame(result_data)

        self._maybe_add_join_keys(result, left_indexer, right_indexer)

        return result


def _get_multiindex_indexer(join_keys, index, sort=False):
    shape = []
    labels = []
    for level, key in zip(index.levels, join_keys):
        llab, rlab, count = _factorize_keys(level, key, sort=False)
        labels.append(rlab)
        shape.append(count)

    left_group_key = get_group_index(labels, shape)
    right_group_key = get_group_index(index.labels, shape)

    left_group_key, right_group_key, max_groups = \
        _factorize_keys(left_group_key, right_group_key,
                        sort=False)

    left_indexer, right_indexer = \
        algos.left_outer_join(com._ensure_int64(left_group_key),
                              com._ensure_int64(right_group_key),
                              max_groups, sort=False)

    return left_indexer, right_indexer


def _get_single_indexer(join_key, index, sort=False):
    left_key, right_key, count = _factorize_keys(join_key, index, sort=sort)

    left_indexer, right_indexer = \
        algos.left_outer_join(com._ensure_int64(left_key),
                              com._ensure_int64(right_key),
                              count, sort=sort)

    return left_indexer, right_indexer


def _left_join_on_index(left_ax, right_ax, join_keys, sort=False):
    join_index = left_ax
    left_indexer = None

    if len(join_keys) > 1:
        if not ((isinstance(right_ax, MultiIndex) and
                 len(join_keys) == right_ax.nlevels)):
            raise AssertionError("If more than one join key is given then "
                                 "'right_ax' must be a MultiIndex and the "
                                 "number of join keys must be the number of "
                                 "levels in right_ax")

        left_tmp, right_indexer = \
            _get_multiindex_indexer(join_keys, right_ax,
                                    sort=sort)
        if sort:
            left_indexer = left_tmp
            join_index = left_ax.take(left_indexer)
    else:
        jkey = join_keys[0]
        if sort:
            left_indexer, right_indexer = \
                _get_single_indexer(jkey, right_ax, sort=sort)
            join_index = left_ax.take(left_indexer)
        else:
            right_indexer = right_ax.get_indexer(jkey)

    return join_index, left_indexer, right_indexer


def _right_outer_join(x, y, max_groups):
    right_indexer, left_indexer = algos.left_outer_join(y, x, max_groups)
    return left_indexer, right_indexer

_join_functions = {
    'inner': algos.inner_join,
    'left': algos.left_outer_join,
    'right': _right_outer_join,
    'outer': algos.full_outer_join,
}


def _factorize_keys(lk, rk, sort=True):
    if com._is_int_or_datetime_dtype(lk) and com._is_int_or_datetime_dtype(rk):
        klass = _hash.Int64Factorizer
        lk = com._ensure_int64(lk)
        rk = com._ensure_int64(rk)
    else:
        klass = _hash.Factorizer
        lk = com._ensure_object(lk)
        rk = com._ensure_object(rk)

    rizer = klass(max(len(lk), len(rk)))

    llab = rizer.factorize(lk)
    rlab = rizer.factorize(rk)

    count = rizer.get_count()

    if sort:
        uniques = rizer.uniques.to_array()
        llab, rlab = _sort_labels(uniques, llab, rlab)

    # NA group
    lmask = llab == -1
    lany = lmask.any()
    rmask = rlab == -1
    rany = rmask.any()

    if lany or rany:
        if lany:
            np.putmask(llab, lmask, count)
        if rany:
            np.putmask(rlab, rmask, count)
        count += 1

    return llab, rlab, count


def _sort_labels(uniques, left, right):
    if not isinstance(uniques, np.ndarray):
        # tuplesafe
        uniques = Index(uniques).values

    sorter = uniques.argsort()

    reverse_indexer = np.empty(len(sorter), dtype=np.int64)
    reverse_indexer.put(sorter, np.arange(len(sorter)))

    new_left = reverse_indexer.take(com._ensure_platform_int(left))
    np.putmask(new_left, left == -1, -1)

    new_right = reverse_indexer.take(com._ensure_platform_int(right))
    np.putmask(new_right, right == -1, -1)

    return new_left, new_right


class _BlockJoinOperation(object):
    """
    BlockJoinOperation made generic for N DataFrames

    Object responsible for orchestrating efficient join operation between two
    BlockManager data structures
    """
    def __init__(self, data_list, join_index, indexers, axis=1, copy=True):
        if axis <= 0:  # pragma: no cover
            raise MergeError('Only axis >= 1 supported for this operation')

        if len(data_list) != len(indexers):
            raise AssertionError("data_list and indexers must have the same "
                                 "length")

        self.units = []
        for data, indexer in zip(data_list, indexers):
            if not data.is_consolidated():
                data = data.consolidate()
            data._set_ref_locs()
            self.units.append(_JoinUnit(data.blocks, indexer))

        self.join_index = join_index
        self.axis = axis
        self.copy = copy
        self.offsets = None

        # do NOT sort
        self.result_items = _concat_indexes([d.items for d in data_list])
        self.result_axes = list(data_list[0].axes)
        self.result_axes[0] = self.result_items
        self.result_axes[axis] = self.join_index

    def _prepare_blocks(self):
        blockmaps = []

        for unit in self.units:
            join_blocks = unit.get_upcasted_blocks()
            type_map = {}
            for blk in join_blocks:
                type_map.setdefault(blk.ftype, []).append(blk)
            blockmaps.append((unit, type_map))

        return blockmaps

    def get_result(self):
        """
        Returns
        -------
        merged : BlockManager
        """
        blockmaps = self._prepare_blocks()
        kinds = _get_merge_block_kinds(blockmaps)

        # maybe want to enable flexible copying <-- what did I mean?
        kind_blocks = []
        for klass in kinds:
            klass_blocks = []
            for unit, mapping in blockmaps:
                if klass in mapping:
                    klass_blocks.extend((unit, b) for b in mapping[klass])

            # blocks that we are going to merge
            kind_blocks.append(klass_blocks)

        # create the merge offsets, essentially where the resultant blocks go in the result
        if not self.result_items.is_unique:

            # length of the merges for each of the klass blocks
            self.offsets = np.zeros(len(blockmaps))
            for kb in kind_blocks:
                kl = list(b.get_merge_length() for unit, b in kb)
                self.offsets += np.array(kl)

        # merge the blocks to create the result blocks
        result_blocks = []
        for klass_blocks in kind_blocks:
            res_blk = self._get_merged_block(klass_blocks)
            result_blocks.append(res_blk)

        return BlockManager(result_blocks, self.result_axes)

    def _get_merged_block(self, to_merge):
        if len(to_merge) > 1:

            # placement set here
            return self._merge_blocks(to_merge)
        else:
            unit, block = to_merge[0]
            blk = unit.reindex_block(block, self.axis,
                                     self.result_items, copy=self.copy)

            # set placement / invalidate on a unique result
            if self.result_items.is_unique and blk._ref_locs is not None:
                if not self.copy:
                    blk = blk.copy()
                blk.set_ref_locs(None)

            return blk


    def _merge_blocks(self, merge_chunks):
        """
        merge_chunks -> [(_JoinUnit, Block)]
        """
        funit, fblock = merge_chunks[0]
        fidx = funit.indexer

        out_shape = list(fblock.get_values().shape)

        n = len(fidx) if fidx is not None else out_shape[self.axis]

        merge_lengths = list(blk.get_merge_length() for unit, blk in merge_chunks)
        out_shape[0] = sum(merge_lengths)
        out_shape[self.axis] = n

        # Should use Fortran order??
        block_dtype = _get_block_dtype([x[1] for x in merge_chunks])
        out = np.empty(out_shape, dtype=block_dtype)

        sofar = 0
        for unit, blk in merge_chunks:
            out_chunk = out[sofar: sofar + len(blk)]
            com.take_nd(blk.get_values(), unit.indexer, self.axis, out=out_chunk)
            sofar += len(blk)

        # does not sort
        new_block_items = _concat_indexes([b.items for _, b in merge_chunks])

        # need to set placement if we have a non-unique result
        # calculate by the existing placement plus the offset in the result set
        placement = None
        if not self.result_items.is_unique:
            placement = []
            offsets = np.append(np.array([0]),self.offsets.cumsum()[:-1])
            for (unit, blk), offset in zip(merge_chunks,offsets):
                placement.extend(blk.ref_locs+offset)

        return make_block(out, new_block_items, self.result_items, placement=placement)


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
            return None
        else:
            mask = self.indexer == -1
            needs_masking = mask.any()
            return (mask, needs_masking)

    def get_upcasted_blocks(self):
        # will short-circuit and not compute needs_masking if indexer is None
        if self.mask_info is not None and self.mask_info[1]:
            return _upcast_blocks(self.blocks)
        return self.blocks

    def reindex_block(self, block, axis, ref_items, copy=True):
        if self.indexer is None:
            result = block.copy() if copy else block
        else:
            result = block.reindex_axis(self.indexer, axis=axis,
                                        mask_info=self.mask_info)
        result.ref_items = ref_items
        return result


def _may_need_upcasting(blocks):
    for block in blocks:
        if isinstance(block, (IntBlock, BoolBlock)) and not isinstance(block, TimeDeltaBlock):
            return True
    return False


def _upcast_blocks(blocks):
    """
    Upcast and consolidate if necessary
    """
    new_blocks = []
    for block in blocks:
        if isinstance(block, TimeDeltaBlock):
            # these are int blocks underlying, but are ok
            newb = block
        elif isinstance(block, IntBlock):
            newb = make_block(block.values.astype(float), block.items,
                              block.ref_items, placement=block._ref_locs)
        elif isinstance(block, BoolBlock):
            newb = make_block(block.values.astype(object), block.items,
                              block.ref_items, placement=block._ref_locs)
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


def _get_merge_block_kinds(blockmaps):
    kinds = set()
    for _, mapping in blockmaps:
        kinds |= set(mapping)
    return kinds


def _get_block_dtype(blocks):
    if len(blocks) == 0:
        return object
    blk1 = blocks[0]
    dtype = blk1.dtype

    if issubclass(dtype.type, np.floating):
        for blk in blocks:
            if blk.dtype.type == np.float64:
                return blk.dtype

    return dtype

#----------------------------------------------------------------------
# Concatenate DataFrame objects


def concat(objs, axis=0, join='outer', join_axes=None, ignore_index=False,
           keys=None, levels=None, names=None, verify_integrity=False):
    """
    Concatenate pandas objects along a particular axis with optional set logic
    along the other axes. Can also add a layer of hierarchical indexing on the
    concatenation axis, which may be useful if the labels are the same (or
    overlapping) on the passed axis number

    Parameters
    ----------
    objs : list or dict of Series, DataFrame, or Panel objects
        If a dict is passed, the sorted keys will be used as the `keys`
        argument, unless it is passed, in which case the values will be
        selected (see below). Any None objects will be dropped silently unless
        they are all None in which case an Exception will be raised
    axis : {0, 1, ...}, default 0
        The axis to concatenate along
    join : {'inner', 'outer'}, default 'outer'
        How to handle indexes on other axis(es)
    join_axes : list of Index objects
        Specific indexes to use for the other n - 1 axes instead of performing
        inner/outer set logic
    verify_integrity : boolean, default False
        Check whether the new concatenated axis contains duplicates. This can
        be very expensive relative to the actual data concatenation
    keys : sequence, default None
        If multiple levels passed, should contain tuples. Construct
        hierarchical index using the passed keys as the outermost level
    levels : list of sequences, default None
        Specific levels (unique values) to use for constructing a
        MultiIndex. Otherwise they will be inferred from the keys
    names : list, default None
        Names for the levels in the resulting hierarchical index
    ignore_index : boolean, default False
        If True, do not use the index values along the concatenation axis. The
        resulting axis will be labeled 0, ..., n - 1. This is useful if you are
        concatenating objects where the concatenation axis does not have
        meaningful indexing information. Note the the index values on the other
        axes are still respected in the join.

    Notes
    -----
    The keys, levels, and names arguments are all optional

    Returns
    -------
    concatenated : type of objects
    """
    op = _Concatenator(objs, axis=axis, join_axes=join_axes,
                       ignore_index=ignore_index, join=join,
                       keys=keys, levels=levels, names=names,
                       verify_integrity=verify_integrity)
    return op.get_result()


class _Concatenator(object):
    """
    Orchestrates a concatenation operation for BlockManagers
    """

    def __init__(self, objs, axis=0, join='outer', join_axes=None,
                 keys=None, levels=None, names=None,
                 ignore_index=False, verify_integrity=False):
        if not isinstance(objs, (list,tuple,types.GeneratorType,dict)):
            raise AssertionError('first argument must be a list-like of pandas '
                                 'objects, you passed an object of type '
                                 '"{0}"'.format(type(objs).__name__))

        if join == 'outer':
            self.intersect = False
        elif join == 'inner':
            self.intersect = True
        else:  # pragma: no cover
            raise ValueError('Only can inner (intersect) or outer (union) '
                             'join the other axis')

        if isinstance(objs, dict):
            if keys is None:
                keys = sorted(objs)
            objs = [objs[k] for k in keys]

        if keys is None:
            objs = [obj for obj in objs if obj is not None]
        else:
            # #1649
            clean_keys = []
            clean_objs = []
            for k, v in zip(keys, objs):
                if v is None:
                    continue
                clean_keys.append(k)
                clean_objs.append(v)
            objs = clean_objs
            keys = clean_keys

        if len(objs) == 0:
            raise Exception('All objects passed were None')

        # consolidate data
        for obj in objs:
            if isinstance(obj, NDFrame):
                obj.consolidate(inplace=True)
        self.objs = objs

        sample = objs[0]

        # Need to flip BlockManager axis in the DataFrame special case
        if isinstance(sample, DataFrame):
            axis = 1 if axis == 0 else 0

        self._is_series = isinstance(sample, ABCSeries)
        if not 0 <= axis <= sample.ndim:
            raise AssertionError("axis must be between 0 and {0}, "
                                 "input was {1}".format(sample.ndim, axis))

        # note: this is the BlockManager axis (since DataFrame is transposed)
        self.axis = axis

        self.join_axes = join_axes

        self.keys = keys
        self.names = names
        self.levels = levels

        self.ignore_index = ignore_index
        self.verify_integrity = verify_integrity

        self.new_axes = self._get_new_axes()

    def get_result(self):
        if self._is_series and self.axis == 0:
            new_data = com._concat_compat([x.get_values() for x in self.objs])
            name = com._consensus_name_attr(self.objs)
            new_data = self._post_merge(new_data)
            return Series(new_data, index=self.new_axes[0], name=name)
        elif self._is_series:
            data = dict(zip(range(len(self.objs)), self.objs))
            index, columns = self.new_axes
            tmpdf = DataFrame(data, index=index)
            if columns is not None:
                tmpdf.columns = columns
            return tmpdf
        else:
            new_data = self._get_concatenated_data()
            new_data = self._post_merge(new_data)
            return self.objs[0]._from_axes(new_data, self.new_axes)

    def _post_merge(self, data):
        if isinstance(data, BlockManager):
            data = data.post_merge(self.objs)
        return data

    def _get_fresh_axis(self):
        return Index(np.arange(len(self._get_concat_axis())))

    def _prepare_blocks(self):
        reindexed_data = self._get_reindexed_data()

        # we are consolidating as we go, so just add the blocks, no-need for dtype mapping
        blockmaps = []
        for data in reindexed_data:
            data = data.consolidate()
            data._set_ref_locs()
            blockmaps.append(data.get_block_map(typ='dict'))
        return blockmaps, reindexed_data

    def _get_concatenated_data(self):
        # need to conform to same other (joined) axes for block join
        blockmaps, rdata = self._prepare_blocks()
        kinds = _get_all_block_kinds(blockmaps)

        try:
            # need to conform to same other (joined) axes for block join
            new_blocks = []
            for kind in kinds:
                klass_blocks = []
                for mapping in blockmaps:
                    l = mapping.get(kind)
                    if l is None:
                        l = [ None ]
                    klass_blocks.extend(l)
                stacked_block = self._concat_blocks(klass_blocks)
                new_blocks.append(stacked_block)

            if self.axis == 0 and self.ignore_index:
                self.new_axes[0] = self._get_fresh_axis()

            for blk in new_blocks:
                blk.ref_items = self.new_axes[0]

            new_data = BlockManager(new_blocks, self.new_axes)

        # Eventual goal would be to move everything to PandasError or other explicit error
        except (Exception, PandasError):  # EAFP

            # should not be possible to fail here for the expected reason with
            # axis = 0
            if self.axis == 0:  # pragma: no cover
                raise

            new_data = {}
            for item in self.new_axes[0]:
                new_data[item] = self._concat_single_item(rdata, item)

        return new_data

    def _get_reindexed_data(self):
        # HACK: ugh

        reindexed_data = []
        axes_to_reindex = list(enumerate(self.new_axes))
        axes_to_reindex.pop(self.axis)

        for obj in self.objs:
            data = obj._data.prepare_for_merge()
            for i, ax in axes_to_reindex:
                data = data.reindex_axis(ax, axis=i, copy=False)
            reindexed_data.append(data)

        return reindexed_data

    def _concat_blocks(self, blocks):

        values_list = [b.get_values() for b in blocks if b is not None]
        concat_values = com._concat_compat(values_list, axis=self.axis)

        if self.axis > 0:
            # Not safe to remove this check, need to profile
            if not _all_indexes_same([b.items for b in blocks]):
                # TODO: Either profile this piece or remove.
                # FIXME: Need to figure out how to test whether this line exists or does not...(unclear if even possible
                #        or maybe would require performance test)
                raise PandasError('dtypes are not consistent throughout '
                                  'DataFrames')
            return make_block(concat_values,
                              blocks[0].items,
                              self.new_axes[0],
                              placement=blocks[0]._ref_locs)
        else:

            offsets = np.r_[0, np.cumsum([len(x._data.axes[0]) for
                                          x in self.objs])]
            indexer = np.concatenate([offsets[i] + b.ref_locs
                                      for i, b in enumerate(blocks)
                                      if b is not None])
            if self.ignore_index:
                concat_items = indexer
            else:
                concat_items = self.new_axes[0].take(indexer)

            if self.ignore_index:
                ref_items = self._get_fresh_axis()
                return make_block(concat_values, concat_items, ref_items)

            block = make_block(concat_values, concat_items, self.new_axes[0])

            # we need to set the ref_locs in this block so we have the mapping
            # as we now have a non-unique index across dtypes, and we need to
            # map the column location to the block location
            # GH3602
            if not self.new_axes[0].is_unique:
                block.set_ref_locs(indexer)

            return block

    def _concat_single_item(self, objs, item):
        # this is called if we don't have consistent dtypes in a row-wise append
        all_values = []
        dtypes = []
        alls = set()

        # figure out the resulting dtype of the combination
        for data, orig in zip(objs, self.objs):
            d = dict([ (t,False) for t in ['object','datetime','timedelta','other'] ])
            if item in orig:
                values = data.get(item)
                if hasattr(values,'to_dense'):
                    values = values.to_dense()
                all_values.append(values)

                dtype = values.dtype

                if issubclass(dtype.type, (np.object_, np.bool_)):
                    d['object'] = True
                    alls.add('object')
                elif is_datetime64_dtype(dtype):
                    d['datetime'] = True
                    alls.add('datetime')
                elif is_timedelta64_dtype(dtype):
                    d['timedelta'] = True
                    alls.add('timedelta')
                else:
                    d['other'] = True
                    alls.add('other')

            else:
                all_values.append(None)
                d['other'] = True
                alls.add('other')

            dtypes.append(d)

        if 'datetime' in alls or 'timedelta' in alls:

            if 'object' in alls or 'other' in alls:

                for v, d in zip(all_values,dtypes):
                    if d.get('datetime') or d.get('timedelta'):
                        pass

                    # if we have all null, then leave a date/time like type
                    # if we have only that type left
                    elif v is None or isnull(v).all():

                        alls.discard('other')
                        alls.discard('object')

        # create the result
        if 'object' in alls:
            empty_dtype, fill_value = np.object_, np.nan
        elif 'other' in alls:
            empty_dtype, fill_value = np.float64, np.nan
        elif 'datetime' in alls:
            empty_dtype, fill_value = 'M8[ns]', tslib.iNaT
        elif 'timedelta' in alls:
            empty_dtype, fill_value = 'm8[ns]', tslib.iNaT
        else: # pragma
            raise AssertionError("invalid dtype determination in concat_single_item")

        to_concat = []
        for obj, item_values in zip(objs, all_values):
            if item_values is None or isnull(item_values).all():
                shape = obj.shape[1:]
                missing_arr = np.empty(shape, dtype=empty_dtype)
                missing_arr.fill(fill_value)
                to_concat.append(missing_arr)
            else:
                to_concat.append(item_values)

        # this method only gets called with axis >= 1
        if self.axis < 1:
            raise AssertionError("axis must be >= 1, input was"
                                 " {0}".format(self.axis))
        return com._concat_compat(to_concat, axis=self.axis - 1)

    def _get_result_dim(self):
        if self._is_series and self.axis == 1:
            return 2
        else:
            return self.objs[0].ndim

    def _get_new_axes(self):
        ndim = self._get_result_dim()
        new_axes = [None] * ndim

        if self.join_axes is None:
            for i in range(ndim):
                if i == self.axis:
                    continue
                new_axes[i] = self._get_comb_axis(i)
        else:
            if len(self.join_axes) != ndim - 1:
                raise AssertionError("length of join_axes must not be "
                                     "equal to {0}".format(ndim - 1))

            # ufff...
            indices = lrange(ndim)
            indices.remove(self.axis)

            for i, ax in zip(indices, self.join_axes):
                new_axes[i] = ax

        if self.ignore_index:
            concat_axis = None
        else:
            concat_axis = self._get_concat_axis()

        new_axes[self.axis] = concat_axis

        return new_axes

    def _get_comb_axis(self, i):
        if self._is_series:
            all_indexes = [x.index for x in self.objs]
        else:
            try:
                all_indexes = [x._data.axes[i] for x in self.objs]
            except IndexError:
                types = [type(x).__name__ for x in self.objs]
                raise TypeError("Cannot concatenate list of %s" % types)

        return _get_combined_index(all_indexes, intersect=self.intersect)

    def _get_concat_axis(self):
        if self._is_series:
            if self.axis == 0:
                indexes = [x.index for x in self.objs]
            elif self.keys is None:
                names = []
                for x in self.objs:
                    if not isinstance(x, Series):
                        raise TypeError("Cannot concatenate type 'Series' "
                                        "with object of type "
                                        "%r" % type(x).__name__)
                    if x.name is not None:
                        names.append(x.name)
                    else:
                        return Index(np.arange(len(self.objs)))
                return Index(names)
            else:
                return _ensure_index(self.keys)
        else:
            indexes = [x._data.axes[self.axis] for x in self.objs]

        if self.keys is None:
            concat_axis = _concat_indexes(indexes)
        else:
            concat_axis = _make_concat_multiindex(indexes, self.keys,
                                                  self.levels, self.names)

        self._maybe_check_integrity(concat_axis)

        return concat_axis

    def _maybe_check_integrity(self, concat_index):
        if self.verify_integrity:
            if not concat_index.is_unique:
                overlap = concat_index.get_duplicates()
                raise ValueError('Indexes have overlapping values: %s'
                                % str(overlap))


def _concat_indexes(indexes):
    return indexes[0].append(indexes[1:])


def _make_concat_multiindex(indexes, keys, levels=None, names=None):
    if ((levels is None and isinstance(keys[0], tuple)) or
            (levels is not None and len(levels) > 1)):
        zipped = lzip(*keys)
        if names is None:
            names = [None] * len(zipped)

        if levels is None:
            levels = [Categorical.from_array(zp).levels for zp in zipped]
        else:
            levels = [_ensure_index(x) for x in levels]
    else:
        zipped = [keys]
        if names is None:
            names = [None]

        if levels is None:
            levels = [_ensure_index(keys)]
        else:
            levels = [_ensure_index(x) for x in levels]

    if not _all_indexes_same(indexes):
        label_list = []

        # things are potentially different sizes, so compute the exact labels
        # for each level and pass those to MultiIndex.from_arrays

        for hlevel, level in zip(zipped, levels):
            to_concat = []
            for key, index in zip(hlevel, indexes):
                try:
                    i = level.get_loc(key)
                except KeyError:
                    raise ValueError('Key %s not in level %s'
                                     % (str(key), str(level)))

                to_concat.append(np.repeat(i, len(index)))
            label_list.append(np.concatenate(to_concat))

        concat_index = _concat_indexes(indexes)

        # these go at the end
        if isinstance(concat_index, MultiIndex):
            levels.extend(concat_index.levels)
            label_list.extend(concat_index.labels)
        else:
            factor = Categorical.from_array(concat_index)
            levels.append(factor.levels)
            label_list.append(factor.labels)

        if len(names) == len(levels):
            names = list(names)
        else:
            # make sure that all of the passed indices have the same nlevels
            if not len(set([ i.nlevels for i in indexes ])) == 1:
                raise AssertionError("Cannot concat indices that do"
                                     " not have the same number of levels")

            # also copies
            names = names + _get_consensus_names(indexes)

        return MultiIndex(levels=levels, labels=label_list, names=names,
                          verify_integrity=False)

    new_index = indexes[0]
    n = len(new_index)
    kpieces = len(indexes)

    # also copies
    new_names = list(names)
    new_levels = list(levels)

    # construct labels
    new_labels = []

    # do something a bit more speedy

    for hlevel, level in zip(zipped, levels):
        hlevel = _ensure_index(hlevel)
        mapped = level.get_indexer(hlevel)

        mask = mapped == -1
        if mask.any():
            raise ValueError('Values not found in passed level: %s'
                             % str(hlevel[mask]))

        new_labels.append(np.repeat(mapped, n))

    if isinstance(new_index, MultiIndex):
        new_levels.extend(new_index.levels)
        new_labels.extend([np.tile(lab, kpieces) for lab in new_index.labels])
    else:
        new_levels.append(new_index)
        new_labels.append(np.tile(np.arange(n), kpieces))

    if len(new_names) < len(new_levels):
        new_names.extend(new_index.names)

    return MultiIndex(levels=new_levels, labels=new_labels, names=new_names,
                      verify_integrity=False)


def _should_fill(lname, rname):
    if not isinstance(lname, compat.string_types) or not isinstance(rname, compat.string_types):
        return True
    return lname == rname


def _any(x):
    return x is not None and len(x) > 0 and any([y is not None for y in x])
