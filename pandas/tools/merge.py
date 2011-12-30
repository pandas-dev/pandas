"""
SQL-style merge routines
"""

import numpy as np

from pandas.core.frame import DataFrame
from pandas.core.index import Index
from pandas.core.internals import _JoinOperation
import pandas.core.common as com

import pandas._tseries as lib
from pandas._sandbox import Factorizer

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
    how : {'left', 'right', 'outer', 'inner'}
        How to handle indexes of the two objects. Default: 'left'
        for joining on index, None otherwise
        * left: use only keys from left frame (SQL: left outer join)
        * right: use only keys from right frame (SQL: right outer join)
        * outer: use union of keys from both frames (SQL: full outer join)
        * inner: use intersection of keys from both frames (SQL: inner join)
    on : label or list

    left_on : label or list

    right_on : label or list

    left_index : boolean, default True

    right_index : boolean, default True

    sort : boolean, default True

    suffixes : 2-length sequence (tuple, list, ...)
        Suffix to apply to overlapping column names in the left and right
        side, respectively
    copy : boolean, default True
        If False, do not copy data unnecessarily

    Examples
    --------

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

class _MergeOperation(object):

    def __init__(self, left, right, how='inner', on=None,
                 left_on=None, right_on=None,
                 left_index=False, right_index=False, sort=True,
                 suffixes=('.x', '.y'), copy=True):
        self.left = left
        self.right = right
        self.how = how

        self.on = _maybe_make_list(on)
        self.left_on = _maybe_make_list(left_on)
        self.right_on = _maybe_make_list(right_on)

        self.copy = copy

        self.suffixes = suffixes

        self.sort = sort

        self.left_index = left_index
        self.right_index = right_index

    def get_result(self):
        # note this function has side effects
        left_join_keys, right_join_keys, join_names = self._get_merge_keys()

        # this is a bit kludgy
        ldata, rdata = self._get_merge_data(join_names)

        # max groups = largest possible number of distinct groups
        left_key, right_key, max_groups = \
            _get_group_keys(left_join_keys, right_join_keys, sort=self.sort)

        join_func = _join_functions[self.how]
        left_indexer, right_indexer = join_func(left_key.astype('i4'),
                                                right_key.astype('i4'),
                                                max_groups)

        new_axis = Index(np.arange(len(left_indexer)))

        # TODO: more efficiently handle group keys to avoid extra consolidation!

        join_op = _JoinOperation(ldata, rdata, new_axis,
                                 left_indexer, right_indexer, axis=1)

        result_data = join_op.get_result(copy=self.copy)
        result = DataFrame(result_data)

        # insert group keys
        for i, name in enumerate(join_names):
            # a faster way?
            key_col = com.take_1d(left_join_keys[i], left_indexer)
            na_indexer = (left_indexer == -1).nonzero()[0]
            right_na_indexer = right_indexer.take(na_indexer)
            key_col.put(na_indexer, com.take_1d(right_join_keys[i],
                                                right_na_indexer))
            result.insert(i, name, key_col)

        return result

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
        rizer = Factorizer(max(len(lk), len(rk)))

        llab, _ = rizer.factorize(lk.astype('O'))
        rlab, _ = rizer.factorize(rk.astype('O'))

        count = rizer.get_count()

        if sort:
            sorter = Index(rizer.uniques).argsort()
            reverse_indexer = np.empty(len(sorter), dtype=np.int32)
            reverse_indexer.put(sorter, np.arange(len(sorter)))

            llab = reverse_indexer.take(llab)
            rlab = reverse_indexer.take(rlab)

            # TODO: na handling

        left_labels.append(llab)
        right_labels.append(rlab)
        group_sizes.append(count)

    left_group_key = get_group_index(left_labels, group_sizes)
    right_group_key = get_group_index(right_labels, group_sizes)
    max_groups = np.prod(group_sizes)

    return left_group_key, right_group_key, max_groups

import pandas._sandbox as sbx

def _maybe_make_list(obj):
    if obj is not None and not isinstance(obj, (tuple, list)):
        return [obj]
    return obj

def _right_outer_join(x, y, max_groups):
    right_indexer, left_indexer = sbx.left_outer_join(y, x, max_groups)
    return left_indexer, right_indexer

_join_functions = {
    'inner' : sbx.inner_join,
    'left' : sbx.left_outer_join,
    'right' : _right_outer_join,
    'outer' : sbx.full_outer_join,
}
