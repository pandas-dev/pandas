"""
SQL-style merge routines
"""

import numpy as np
from pandas.compat import range, lrange, lzip, zip, map, filter
import pandas.compat as compat
from pandas.core.categorical import Categorical
from pandas.core.frame import DataFrame, _merge_doc
from pandas.core.generic import NDFrame
from pandas.core.series import Series
from pandas.core.index import (Index, MultiIndex, _get_combined_index,
                               _ensure_index, _get_consensus_names,
                               _all_indexes_same)
from pandas.core.internals import (items_overlap_with_suffix,
                                   concatenate_block_managers)
from pandas.util.decorators import Appender, Substitution
from pandas.core.common import ABCSeries, isnull

import pandas.core.common as com

import pandas.algos as algos
import pandas.hashtable as _hash


@Substitution('\nleft : DataFrame')
@Appender(_merge_doc, indents=0)
def merge(left, right, how='inner', on=None, left_on=None, right_on=None,
          left_index=False, right_index=False, sort=False,
          suffixes=('_x', '_y'), copy=True, indicator=False):
    op = _MergeOperation(left, right, how=how, on=on, left_on=left_on,
                         right_on=right_on, left_index=left_index,
                         right_index=right_index, sort=sort, suffixes=suffixes,
                         copy=copy, indicator=indicator)
    return op.get_result()
if __debug__:
    merge.__doc__ = _merge_doc % '\nleft : DataFrame'


class MergeError(ValueError):
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
        The output type will the be same as 'left', if it is a subclass
        of DataFrame.
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
                 suffixes=('_x', '_y'), copy=True, indicator=False):
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

        self.indicator = indicator

        if isinstance(self.indicator, compat.string_types):
            self.indicator_name = self.indicator
        elif isinstance(self.indicator, bool):
            self.indicator_name = '_merge' if self.indicator else None
        else:
            raise ValueError('indicator option can only accept boolean or string arguments')


        # note this function has side effects
        (self.left_join_keys,
         self.right_join_keys,
         self.join_names) = self._get_merge_keys()

    def get_result(self):
        if self.indicator:
            self.left, self.right = self._indicator_pre_merge(self.left, self.right)

        join_index, left_indexer, right_indexer = self._get_join_info()

        ldata, rdata = self.left._data, self.right._data
        lsuf, rsuf = self.suffixes

        llabels, rlabels = items_overlap_with_suffix(ldata.items, lsuf,
                                                     rdata.items, rsuf)

        lindexers = {1: left_indexer} if left_indexer is not None else {}
        rindexers = {1: right_indexer} if right_indexer is not None else {}

        result_data = concatenate_block_managers(
            [(ldata, lindexers), (rdata, rindexers)],
            axes=[llabels.append(rlabels), join_index],
            concat_axis=0, copy=self.copy)

        typ = self.left._constructor
        result = typ(result_data).__finalize__(self, method='merge')

        if self.indicator:
            result = self._indicator_post_merge(result)

        self._maybe_add_join_keys(result, left_indexer, right_indexer)

        return result

    def _indicator_pre_merge(self, left, right):

        columns = left.columns.union(right.columns)

        for i in ['_left_indicator', '_right_indicator']:
            if i in columns:
                raise ValueError("Cannot use `indicator=True` option when data contains a column named {}".format(i))
        if self.indicator_name in columns:
            raise ValueError("Cannot use name of an existing column for indicator column")

        left = left.copy()
        right = right.copy()

        left['_left_indicator'] = 1
        left['_left_indicator'] = left['_left_indicator'].astype('int8')

        right['_right_indicator'] = 2
        right['_right_indicator'] = right['_right_indicator'].astype('int8')

        return left, right

    def _indicator_post_merge(self, result):

        result['_left_indicator'] = result['_left_indicator'].fillna(0)
        result['_right_indicator'] = result['_right_indicator'].fillna(0)

        result[self.indicator_name] = Categorical((result['_left_indicator'] + result['_right_indicator']), categories=[1,2,3])
        result[self.indicator_name] = result[self.indicator_name].cat.rename_categories(['left_only', 'right_only', 'both'])

        result = result.drop(labels=['_left_indicator', '_right_indicator'], axis=1)

        return result

    def _maybe_add_join_keys(self, result, left_indexer, right_indexer):
        # insert group keys

        keys = zip(self.join_names, self.left_on, self.right_on)
        for i, (name, lname, rname) in enumerate(keys):
            if not _should_fill(lname, rname):
                continue

            if name in result:
                key_indexer = result.columns.get_loc(name)

                if left_indexer is not None and right_indexer is not None:

                    if name in self.left:
                        if len(self.left) == 0:
                            continue

                        na_indexer = (left_indexer == -1).nonzero()[0]
                        if len(na_indexer) == 0:
                            continue

                        right_na_indexer = right_indexer.take(na_indexer)
                        result.iloc[na_indexer,key_indexer] = com.take_1d(self.right_join_keys[i],
                                                                          right_na_indexer)
                    elif name in self.right:
                        if len(self.right) == 0:
                            continue

                        na_indexer = (right_indexer == -1).nonzero()[0]
                        if len(na_indexer) == 0:
                            continue

                        left_na_indexer = left_indexer.take(na_indexer)
                        result.iloc[na_indexer,key_indexer] = com.take_1d(self.left_join_keys[i],
                                                                          left_na_indexer)
            elif left_indexer is not None \
                    and isinstance(self.left_join_keys[i], np.ndarray):

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
                if len(self.left) > 0:
                    join_index = self.left.index.take(left_indexer)
                else:
                    join_index = self.right.index.take(right_indexer)
                    left_indexer = np.array([-1] * len(join_index))
            elif self.left_index:
                if len(self.right) > 0:
                    join_index = self.right.index.take(right_indexer)
                else:
                    join_index = self.left.index.take(left_indexer)
                    right_indexer = np.array([-1] * len(join_index))
            else:
                join_index = Index(np.arange(len(left_indexer)))

        if len(join_index) == 0:
            join_index = join_index.astype(object)
        return join_index, left_indexer, right_indexer

    def _get_merge_data(self):
        """
        Handles overlapping column names etc.
        """
        ldata, rdata = self.left._data, self.right._data
        lsuf, rsuf = self.suffixes

        llabels, rlabels = items_overlap_with_suffix(
            ldata.items, lsuf, rdata.items, rsuf)

        if not llabels.equals(ldata.items):
            ldata = ldata.copy(deep=False)
            ldata.set_axis(0, llabels)

        if not rlabels.equals(rdata.items):
            rdata = rdata.copy(deep=False)
            rdata.set_axis(0, rlabels)

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
                        right_keys.append(right[rk]._values)
                        join_names.append(rk)
                else:
                    if not is_rkey(rk):
                        right_keys.append(right[rk]._values)
                        if lk == rk:
                            # avoid key upcast in corner case (length-0)
                            if len(left) > 0:
                                right_drop.append(rk)
                            else:
                                left_drop.append(lk)
                    else:
                        right_keys.append(rk)
                    left_keys.append(left[lk]._values)
                    join_names.append(lk)
        elif _any(self.left_on):
            for k in self.left_on:
                if is_lkey(k):
                    left_keys.append(k)
                    join_names.append(None)
                else:
                    left_keys.append(left[k]._values)
                    join_names.append(k)
            if isinstance(self.right.index, MultiIndex):
                right_keys = [lev._values.take(lab)
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
                    right_keys.append(right[k]._values)
                    join_names.append(k)
            if isinstance(self.left.index, MultiIndex):
                left_keys = [lev._values.take(lab)
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
                # use the common columns
                common_cols = self.left.columns.intersection(
                    self.right.columns)
                if len(common_cols) == 0:
                    raise MergeError('No common columns to perform merge on')
                if not common_cols.is_unique:
                    raise MergeError("Data columns not unique: %s"
                                     % repr(common_cols))
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
    from functools import partial

    assert len(left_keys) == len(right_keys), \
            'left_key and right_keys must be the same length'

    # bind `sort` arg. of _factorize_keys
    fkeys = partial(_factorize_keys, sort=sort)

    # get left & right join labels and num. of levels at each location
    llab, rlab, shape = map(list, zip( * map(fkeys, left_keys, right_keys)))

    # get flat i8 keys from label lists
    lkey, rkey = _get_join_keys(llab, rlab, shape, sort)

    # factorize keys to a dense i8 space
    # `count` is the num. of unique keys
    # set(lkey) | set(rkey) == range(count)
    lkey, rkey, count = fkeys(lkey, rkey)

    # preserve left frame order if how == 'left' and sort == False
    kwargs = {'sort':sort} if how == 'left' else {}
    join_func = _join_functions[how]
    return join_func(lkey, rkey, count, **kwargs)


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
        ldata, rdata = self.left._data, self.right._data
        lsuf, rsuf = self.suffixes

        llabels, rlabels = items_overlap_with_suffix(ldata.items, lsuf,
                                                     rdata.items, rsuf)

        if self.fill_method == 'ffill':
            left_join_indexer = algos.ffill_indexer(left_indexer)
            right_join_indexer = algos.ffill_indexer(right_indexer)
        else:
            left_join_indexer = left_indexer
            right_join_indexer = right_indexer

        lindexers = {1: left_join_indexer} if left_join_indexer is not None else {}
        rindexers = {1: right_join_indexer} if right_join_indexer is not None else {}

        result_data = concatenate_block_managers(
            [(ldata, lindexers), (rdata, rindexers)],
            axes=[llabels.append(rlabels), join_index],
            concat_axis=0, copy=self.copy)

        typ = self.left._constructor
        result = typ(result_data).__finalize__(self, method='ordered_merge')

        self._maybe_add_join_keys(result, left_indexer, right_indexer)

        return result


def _get_multiindex_indexer(join_keys, index, sort):
    from functools import partial

    # bind `sort` argument
    fkeys = partial(_factorize_keys, sort=sort)

    # left & right join labels and num. of levels at each location
    rlab, llab, shape = map(list, zip( * map(fkeys, index.levels, join_keys)))
    if sort:
        rlab = list(map(np.take, rlab, index.labels))
    else:
        i8copy = lambda a: a.astype('i8', subok=False, copy=True)
        rlab = list(map(i8copy, index.labels))

    # fix right labels if there were any nulls
    for i in range(len(join_keys)):
        mask = index.labels[i] == -1
        if mask.any():
            # check if there already was any nulls at this location
            # if there was, it is factorized to `shape[i] - 1`
            a = join_keys[i][llab[i] == shape[i] - 1]
            if a.size == 0 or not a[0] != a[0]:
                shape[i] += 1

            rlab[i][mask] = shape[i] - 1

    # get flat i8 join keys
    lkey, rkey = _get_join_keys(llab, rlab, shape, sort)

    # factorize keys to a dense i8 space
    lkey, rkey, count = fkeys(lkey, rkey)

    return algos.left_outer_join(lkey, rkey, count, sort=sort)


def _get_single_indexer(join_key, index, sort=False):
    left_key, right_key, count = _factorize_keys(join_key, index, sort=sort)

    left_indexer, right_indexer = \
        algos.left_outer_join(com._ensure_int64(left_key),
                              com._ensure_int64(right_key),
                              count, sort=sort)

    return left_indexer, right_indexer


def _left_join_on_index(left_ax, right_ax, join_keys, sort=False):
    if len(join_keys) > 1:
        if not ((isinstance(right_ax, MultiIndex) and
                 len(join_keys) == right_ax.nlevels)):
            raise AssertionError("If more than one join key is given then "
                                 "'right_ax' must be a MultiIndex and the "
                                 "number of join keys must be the number of "
                                 "levels in right_ax")

        left_indexer, right_indexer = \
            _get_multiindex_indexer(join_keys, right_ax, sort=sort)
    else:
        jkey = join_keys[0]

        left_indexer, right_indexer = \
            _get_single_indexer(jkey, right_ax, sort=sort)

    if sort or len(left_ax) != len(left_indexer):
        # if asked to sort or there are 1-to-many matches
        join_index = left_ax.take(left_indexer)
        return join_index, left_indexer, right_indexer

    # left frame preserves order & length of its index
    return left_ax, None, right_indexer


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
    if com.is_datetime64tz_dtype(lk) and com.is_datetime64tz_dtype(rk):
        lk = lk.values
        rk = rk.values
    if com.is_int_or_datetime_dtype(lk) and com.is_int_or_datetime_dtype(rk):
        klass = _hash.Int64Factorizer
        lk = com._ensure_int64(com._values_from_object(lk))
        rk = com._ensure_int64(com._values_from_object(rk))
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


def _get_join_keys(llab, rlab, shape, sort):
    from pandas.core.groupby import _int64_overflow_possible

    # how many levels can be done without overflow
    pred = lambda i: not _int64_overflow_possible(shape[:i])
    nlev = next(filter(pred, range(len(shape), 0, -1)))

    # get keys for the first `nlev` levels
    stride = np.prod(shape[1:nlev], dtype='i8')
    lkey = stride * llab[0].astype('i8', subok=False, copy=False)
    rkey = stride * rlab[0].astype('i8', subok=False, copy=False)

    for i in range(1, nlev):
        stride //= shape[i]
        lkey += llab[i] * stride
        rkey += rlab[i] * stride

    if nlev == len(shape):  # all done!
        return lkey, rkey

    # densify current keys to avoid overflow
    lkey, rkey, count = _factorize_keys(lkey, rkey, sort=sort)

    llab = [lkey] + llab[nlev:]
    rlab = [rkey] + rlab[nlev:]
    shape = [count] + shape[nlev:]

    return _get_join_keys(llab, rlab, shape, sort)

#----------------------------------------------------------------------
# Concatenate DataFrame objects


def concat(objs, axis=0, join='outer', join_axes=None, ignore_index=False,
           keys=None, levels=None, names=None, verify_integrity=False, copy=True):
    """
    Concatenate pandas objects along a particular axis with optional set logic
    along the other axes. Can also add a layer of hierarchical indexing on the
    concatenation axis, which may be useful if the labels are the same (or
    overlapping) on the passed axis number

    Parameters
    ----------
    objs : a sequence or mapping of Series, DataFrame, or Panel objects
        If a dict is passed, the sorted keys will be used as the `keys`
        argument, unless it is passed, in which case the values will be
        selected (see below). Any None objects will be dropped silently unless
        they are all None in which case a ValueError will be raised
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
    copy : boolean, default True
        If False, do not copy data unnecessarily

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
                       verify_integrity=verify_integrity,
                       copy=copy)
    return op.get_result()


class _Concatenator(object):
    """
    Orchestrates a concatenation operation for BlockManagers
    """

    def __init__(self, objs, axis=0, join='outer', join_axes=None,
                 keys=None, levels=None, names=None,
                 ignore_index=False, verify_integrity=False, copy=True):
        if isinstance(objs, (NDFrame, compat.string_types)):
            raise TypeError('first argument must be an iterable of pandas '
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
        else:
            objs = list(objs)

        if len(objs) == 0:
            raise ValueError('No objects to concatenate')

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
            raise ValueError('All objects passed were None')

        # consolidate data & figure out what our result ndim is going to be
        ndims = set()
        for obj in objs:
            if not isinstance(obj, NDFrame):
                raise TypeError("cannot concatenate a non-NDFrame object")

            # consolidate
            obj.consolidate(inplace=True)
            ndims.add(obj.ndim)

        # get the sample
        # want the higest ndim that we have, and must be non-empty
        # unless all objs are empty
        sample = None
        if len(ndims) > 1:
            max_ndim = max(ndims)
            for obj in objs:
                if obj.ndim == max_ndim and np.sum(obj.shape):
                    sample = obj
                    break

        else:
            # filter out the empties
            # if we have not multi-index possibiltes
            df = DataFrame([ obj.shape for obj in objs ]).sum(1)
            non_empties = df[df!=0]
            if len(non_empties) and (keys is None and names is None and levels is None and join_axes is None):
                objs = [ objs[i] for i in non_empties.index ]
                sample = objs[0]

        if sample is None:
            sample = objs[0]
        self.objs = objs

        # Need to flip BlockManager axis in the DataFrame special case
        self._is_frame = isinstance(sample, DataFrame)
        if self._is_frame:
            axis = 1 if axis == 0 else 0

        self._is_series = isinstance(sample, ABCSeries)
        if not 0 <= axis <= sample.ndim:
            raise AssertionError("axis must be between 0 and {0}, "
                                 "input was {1}".format(sample.ndim, axis))

        # if we have mixed ndims, then convert to highest ndim
        # creating column numbers as needed
        if len(ndims) > 1:
            current_column = 0
            max_ndim = sample.ndim
            self.objs, objs = [], self.objs
            for obj in objs:

                ndim = obj.ndim
                if ndim == max_ndim:
                    pass

                elif ndim != max_ndim-1:
                    raise ValueError("cannot concatenate unaligned mixed "
                                     "dimensional NDFrame objects")

                else:
                    name = getattr(obj,'name',None)
                    if ignore_index or name is None:
                        name = current_column
                        current_column += 1

                    # doing a row-wise concatenation so need everything
                    # to line up
                    if self._is_frame and axis == 1:
                        name = 0
                    obj = sample._constructor({ name : obj })

                self.objs.append(obj)

        # note: this is the BlockManager axis (since DataFrame is transposed)
        self.axis = axis
        self.join_axes = join_axes
        self.keys = keys
        self.names = names
        self.levels = levels

        self.ignore_index = ignore_index
        self.verify_integrity = verify_integrity
        self.copy = copy

        self.new_axes = self._get_new_axes()

    def get_result(self):

        # series only
        if self._is_series:

            # stack blocks
            if self.axis == 0:
                new_data = com._concat_compat([x._values for x in self.objs])
                name = com._consensus_name_attr(self.objs)
                return Series(new_data, index=self.new_axes[0], name=name).__finalize__(self, method='concat')

            # combine as columns in a frame
            else:
                data = dict(zip(range(len(self.objs)), self.objs))
                index, columns = self.new_axes
                tmpdf = DataFrame(data, index=index)
                # checks if the column variable already stores valid column names (because set via the 'key' argument
                # in the 'concat' function call. If that's not the case, use the series names as column names
                if columns.equals(Index(np.arange(len(self.objs)))) and not self.ignore_index:
                    columns = np.array([ data[i].name for i in range(len(data)) ], dtype='object')
                    indexer = isnull(columns)
                    if indexer.any():
                        columns[indexer] = np.arange(len(indexer[indexer]))
                tmpdf.columns = columns
                return tmpdf.__finalize__(self, method='concat')

        # combine block managers
        else:
            mgrs_indexers = []
            for obj in self.objs:
                mgr = obj._data
                indexers = {}
                for ax, new_labels in enumerate(self.new_axes):
                    if ax == self.axis:
                        # Suppress reindexing on concat axis
                        continue

                    obj_labels = mgr.axes[ax]
                    if not new_labels.equals(obj_labels):
                        indexers[ax] = obj_labels.reindex(new_labels)[1]

                mgrs_indexers.append((obj._data, indexers))

            new_data = concatenate_block_managers(
                mgrs_indexers, self.new_axes, concat_axis=self.axis, copy=self.copy)
            if not self.copy:
                new_data._consolidate_inplace()

            return self.objs[0]._from_axes(new_data, self.new_axes).__finalize__(self, method='concat')

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

        new_axes[self.axis] = self._get_concat_axis()
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
        """
        Return index to be used along concatenation axis.
        """
        if self._is_series:
            if self.axis == 0:
                indexes = [x.index for x in self.objs]
            elif self.ignore_index:
                idx = Index(np.arange(len(self.objs)))
                idx.is_unique = True  # arange is always unique
                return idx
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
                        idx = Index(np.arange(len(self.objs)))
                        idx.is_unique = True
                        return idx

                return Index(names)
            else:
                return _ensure_index(self.keys)
        else:
            indexes = [x._data.axes[self.axis] for x in self.objs]

        if self.ignore_index:
            idx = Index(np.arange(sum(len(i) for i in indexes)))
            idx.is_unique = True
            return idx

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
            levels = [Categorical.from_array(zp, ordered=True).categories for zp in zipped]
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
            factor = Categorical.from_array(concat_index, ordered=True)
            levels.append(factor.categories)
            label_list.append(factor.codes)

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
