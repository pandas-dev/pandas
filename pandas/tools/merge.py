"""
SQL-style merge routines
"""

import numpy as np

from pandas.core.frame import DataFrame
from pandas.core.index import Index
from pandas.core.internals import _JoinOperation

import pandas._tseries as lib
from pandas._sandbox import Factorizer

def merge(left, right, how='inner', cols=None, left_cols=None, right_cols=None,
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
        * left: use only keys from left frame
        * right: use only keys from right frame
        * outer: use union of keys from both frames
        * inner: use intersection of keys from both frames
    cols
    left_cols
    right_cols
    left_index
    right_index
    sort
    suffixes
    copy : boolean, default True
        If False, do not copy data unnecessarily

    Examples
    --------

    Returns
    -------
    merged : DataFrame
    """
    left_join_keys, right_join_keys = _get_merge_keys(left, right, cols,
                                                      left_cols, right_cols,
                                                      left_index, right_index)

    # max groups = largest possible number of distinct groups
    left_key, right_key, max_groups = _get_group_keys(left_join_keys,
                                                      right_join_keys)

    join_func = _join_functions[how]
    left_indexer, right_indexer = join_func(left_key, right_key, max_groups)
    new_axis = Index(np.arange(len(left_indexer)))

    join_op = _JoinOperation(left, right, new_axis, left_indexer,
                             right_indexer, axis=1)
    result_data = join_op.get_result(copy=copy)
    return DataFrame(result_data)

class _MergeOperation(object):

    def __init__(self, left, right, how='inner', cols=None,
                 left_cols=None, right_cols=None,
                 left_index=False, right_index=False, sort=True,
                 suffixes=('.x', '.y'), copy=True):
        pass

def _get_merge_keys(left, right, cols, left_cols, right_cols,
                    left_index=False, right_index=False):
    """

    Parameters
    ----------

    Returns
    -------

    """
    if on is None:
        pass
    else:
        pass

def _get_group_keys(left_keys, right_keys):
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

        left_labels.append(llab)
        right_labels.append(rlab)
        group_sizes.append(rizer.get_count())

    left_group_key = get_group_index(left_labels, group_sizes)
    right_group_key = get_group_index(right_labels, group_sizes)
    max_groups = np.prod(group_sizes)

    return left_group_key, right_group_key, max_groups

import pandas._sandbox as sbx

def _right_outer_join(x, y):
    right_indexer, left_indexer = sbx.left_outer_join(y, x)
    return left_indexer, right_indexer

_join_functions = {
    'inner' : sbx.inner_join,
    'left' : sbx.left_outer_join,
    'right' : _right_outer_join,
    'outer' : sbx.full_outer_join,
}
