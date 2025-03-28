"""
Algorithms that Involve Multiple DataFrames
===========================================

The pandas operations ``concat``, ``join``, and ``merge`` combine multiple
DataFrames.  This module contains analogous algorithms in the parallel case.

There are two important cases:

1.  We combine along a partitioned index
2.  We combine along an unpartitioned index or other column

In the first case we know which partitions of each dataframe interact with
which others.  This lets us be significantly more clever and efficient.

In the second case each partition from one dataset interacts with all
partitions from the other.  We handle this through a shuffle operation.

Partitioned Joins
-----------------

In the first case where we join along a partitioned index we proceed in the
following stages.

1.  Align the partitions of all inputs to be the same.  This involves a call
    to ``dd.repartition`` which will split up and concat existing partitions as
    necessary.  After this step all inputs have partitions that align with
    each other.  This step is relatively cheap.
    See the function ``align_partitions``.
2.  Remove unnecessary partitions based on the type of join we perform (left,
    right, inner, outer).  We can do this at the partition level before any
    computation happens.  We'll do it again on each partition when we call the
    in-memory function.  See the function ``require``.
3.  Embarrassingly parallel calls to ``pd.concat``, ``pd.join``, or
    ``pd.merge``.  Now that the data is aligned and unnecessary blocks have
    been removed we can rely on the fast in-memory Pandas join machinery to
    execute joins per-partition.  We know that all intersecting records exist
    within the same partition


Hash Joins via Shuffle
----------------------

When we join along an unpartitioned index or along an arbitrary column any
partition from one input might interact with any partition in another.  In
this case we perform a hash-join by shuffling data in each input by that
column.  This results in new inputs with the same partition structure cleanly
separated along that column.

We proceed with hash joins in the following stages:

1.  Shuffle each input on the specified column.  See the function
    ``dd.shuffle``.
2.  Perform embarrassingly parallel join across shuffled inputs.
"""

from __future__ import annotations

import pickle
import warnings

import numpy as np
import pandas as pd
from pandas.api.types import is_dtype_equal

from dask.dataframe import methods
from dask.dataframe.core import _concat
from dask.dataframe.dispatch import group_split_dispatch, hash_object_dispatch
from dask.dataframe.shuffle import partitioning_index, shuffle_group
from dask.dataframe.utils import asciitable

###############################################################
# Join / Merge
###############################################################


def merge_chunk(
    lhs,
    *args,
    result_meta,
    **kwargs,
):
    rhs, *args = args
    left_index = kwargs.get("left_index", False)
    right_index = kwargs.get("right_index", False)
    empty_index_dtype = result_meta.index.dtype
    categorical_columns = result_meta.select_dtypes(include="category").columns

    if categorical_columns is not None:
        for col in categorical_columns:
            left = None
            right = None

            if col in lhs:
                left = lhs[col]
            elif col == kwargs.get("right_on", None) and left_index:
                if isinstance(lhs.index.dtype, pd.CategoricalDtype):
                    left = lhs.index

            if col in rhs:
                right = rhs[col]
            elif col == kwargs.get("left_on", None) and right_index:
                if isinstance(rhs.index.dtype, pd.CategoricalDtype):
                    right = rhs.index

            dtype = "category"
            if left is not None and right is not None:
                dtype = methods.union_categoricals(
                    [left.astype("category"), right.astype("category")]
                ).dtype

            if left is not None:
                if isinstance(left, pd.Index):
                    lhs.index = left.astype(dtype)
                else:
                    lhs = lhs.assign(**{col: left.astype(dtype)})
            if right is not None:
                if isinstance(right, pd.Index):
                    rhs.index = right.astype(dtype)
                else:
                    rhs = rhs.assign(**{col: right.astype(dtype)})

    if len(args) and args[0] == "leftsemi" or kwargs.get("how", None) == "leftsemi":
        if isinstance(rhs, (pd.DataFrame, pd.Series)):
            # otherwise it's cudf
            rhs = rhs.drop_duplicates()
            if len(args):
                args[0] = "inner"
            else:
                kwargs["how"] = "inner"
    out = lhs.merge(rhs, *args, **kwargs)

    # Workaround for pandas bug where if the left frame of a merge operation is
    # empty, the resulting dataframe can have columns in the wrong order.
    # https://github.com/pandas-dev/pandas/issues/9937
    if len(lhs) == 0:
        out = out[result_meta.columns]

    # Workaround pandas bug where if the output result of a merge operation is
    # an empty dataframe, the output index is `int64` in all cases, regardless
    # of input dtypes.
    if len(out) == 0 and empty_index_dtype is not None:
        out.index = out.index.astype(empty_index_dtype)
    return out


def warn_dtype_mismatch(left, right, left_on, right_on):
    """Checks for merge column dtype mismatches and throws a warning (#4574)"""

    if not isinstance(left_on, list):
        left_on = [left_on]
    if not isinstance(right_on, list):
        right_on = [right_on]

    if all(col in left.columns for col in left_on) and all(
        col in right.columns for col in right_on
    ):
        dtype_mism = [
            ((lo, ro), left.dtypes[lo], right.dtypes[ro])
            for lo, ro in zip(left_on, right_on)
            if not is_dtype_equal(left.dtypes[lo], right.dtypes[ro])
        ]

        if dtype_mism:
            col_tb = asciitable(
                ("Merge columns", "left dtype", "right dtype"), dtype_mism
            )

            warnings.warn(
                (
                    "Merging dataframes with merge column data "
                    "type mismatches: \n{}\nCast dtypes explicitly to "
                    "avoid unexpected results."
                ).format(col_tb)
            )


###############################################################
# ASOF Join
###############################################################


def pair_partitions(L, R):
    """Returns which partitions to pair for the merge_asof algorithm and the
    bounds on which to split them up
    """
    result = []

    n, m = len(L) - 1, len(R) - 1
    i, j = 0, -1
    while j + 1 < m and R[j + 1] <= L[i]:
        j += 1
    J = []
    while i < n:
        partition = max(0, min(m - 1, j))
        lower = R[j] if j >= 0 and R[j] > L[i] else None
        upper = (
            R[j + 1]
            if j + 1 < m
            and (R[j + 1] < L[i + 1] or R[j + 1] == L[i + 1] and i == n - 1)
            else None
        )

        J.append((partition, lower, upper))

        i1 = i + 1 if j + 1 == m or (i + 1 < n and R[j + 1] >= L[i + 1]) else i
        j1 = j + 1 if i + 1 == n or (j + 1 < m and L[i + 1] >= R[j + 1]) else j
        if i1 > i:
            result.append(J)
            J = []
        elif i == n - 1 and R[j1] > L[n]:
            result.append(J)
            break
        i, j = i1, j1

    return result


def merge_asof_padded(left, right, prev=None, next=None, **kwargs):
    """merge_asof but potentially adding rows to the beginning/end of right"""
    frames = []
    if prev is not None:
        frames.append(prev)
    frames.append(right)
    if next is not None:
        frames.append(next)

    frame = pd.concat(frames)
    result = pd.merge_asof(left, frame, **kwargs)
    # pd.merge_asof() resets index name (and dtype) if left is empty df
    if result.index.name != left.index.name:
        result.index.name = left.index.name
    return result


###############################################################
# Concat
###############################################################


def concat_and_check(dfs, ignore_order=False):
    if len(set(map(len, dfs))) != 1:
        raise ValueError("Concatenated DataFrames of different lengths")
    return methods.concat(dfs, axis=1, ignore_order=ignore_order)


def _contains_index_name(df, columns_or_index):
    """
    Test whether ``columns_or_index`` contains a reference
    to the index of ``df

    This is the local (non-collection) version of
    ``dask.core.DataFrame._contains_index_name``.
    """

    def _is_index_level_reference(x, key):
        return (
            x.index.name is not None
            and (np.isscalar(key) or isinstance(key, tuple))
            and key == x.index.name
            and key not in getattr(x, "columns", ())
        )

    if isinstance(columns_or_index, list):
        return any(_is_index_level_reference(df, n) for n in columns_or_index)
    else:
        return _is_index_level_reference(df, columns_or_index)


def _select_columns_or_index(df, columns_or_index):
    """
    Returns a DataFrame with columns corresponding to each
    column or index level in columns_or_index.  If included,
    the column corresponding to the index level is named _index.

    This is the local (non-collection) version of
    ``dask.core.DataFrame._select_columns_or_index``.
    """

    def _is_column_label_reference(df, key):
        return (np.isscalar(key) or isinstance(key, tuple)) and key in df.columns

    # Ensure columns_or_index is a list
    columns_or_index = (
        columns_or_index if isinstance(columns_or_index, list) else [columns_or_index]
    )

    column_names = [n for n in columns_or_index if _is_column_label_reference(df, n)]

    selected_df = df[column_names]
    if _contains_index_name(df, columns_or_index):
        # Index name was included
        selected_df = selected_df.assign(_index=df.index)

    return selected_df


def _split_partition(df, on, nsplits):
    """
    Split-by-hash a DataFrame into `nsplits` groups.

    Hashing will be performed on the columns or index specified by `on`.
    """

    if isinstance(on, bytes):
        on = pickle.loads(on)

    if isinstance(on, str) or pd.api.types.is_list_like(on):
        # If `on` is a column name or list of column names, we
        # can hash/split by those columns.
        on = [on] if isinstance(on, str) else list(on)
        nset = set(on)
        if nset.intersection(set(df.columns)) == nset:
            o = df[on]
            dtypes = {}
            for col, dtype in o.dtypes.items():
                if pd.api.types.is_numeric_dtype(dtype):
                    dtypes[col] = np.float64
            if not dtypes:
                ind = hash_object_dispatch(df[on], index=False)
            else:
                ind = hash_object_dispatch(df[on].astype(dtypes), index=False)

            ind = ind % nsplits
            return group_split_dispatch(df, ind, nsplits, ignore_index=False)

    # We are not joining (purely) on columns.  Need to
    # add a "_partitions" column to perform the split.
    from dask.dataframe.dask_expr._collection import FrameBase

    if not isinstance(on, FrameBase):
        on = _select_columns_or_index(df, on)

    dtypes = {}
    for col, dtype in on.dtypes.items():
        if pd.api.types.is_numeric_dtype(dtype):
            dtypes[col] = np.float64
    if not dtypes:
        dtypes = None

    partitions = partitioning_index(on, nsplits, cast_dtype=dtypes)
    df2 = df.assign(_partitions=partitions)
    return shuffle_group(
        df2,
        ["_partitions"],
        0,
        nsplits,
        nsplits,
        False,
        nsplits,
    )


def _concat_wrapper(dfs):
    """Concat and remove temporary "_partitions" column"""
    df = _concat(dfs, False)
    if "_partitions" in df.columns:
        del df["_partitions"]
    return df


def _merge_chunk_wrapper(*args, **kwargs):
    return merge_chunk(
        *args,
        **{
            k: pickle.loads(v) if isinstance(v, bytes) else v for k, v in kwargs.items()
        },
    )
