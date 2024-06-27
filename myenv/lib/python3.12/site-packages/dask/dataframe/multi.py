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
    ``dask.dataframe.shuffle.shuffle``.
2.  Perform embarrassingly parallel join across shuffled inputs.
"""
from __future__ import annotations

import math
import pickle
import warnings
from functools import partial, wraps

import numpy as np
import pandas as pd
from pandas.api.types import is_dtype_equal
from tlz import merge_sorted, unique

from dask.base import is_dask_collection, tokenize
from dask.dataframe import methods
from dask.dataframe.core import (
    DataFrame,
    Index,
    Series,
    _concat,
    _deprecated_kwarg,
    _Frame,
    _maybe_from_pandas,
    is_broadcastable,
    map_partitions,
    new_dd_object,
    prefix_reduction,
    suffix_reduction,
)
from dask.dataframe.dispatch import group_split_dispatch, hash_object_dispatch
from dask.dataframe.io import from_pandas
from dask.dataframe.shuffle import (
    partitioning_index,
    rearrange_by_divisions,
    shuffle,
    shuffle_group,
)
from dask.dataframe.utils import (
    asciitable,
    check_meta,
    is_dataframe_like,
    is_series_like,
    make_meta,
    strip_unknown_categories,
)
from dask.highlevelgraph import HighLevelGraph
from dask.layers import BroadcastJoinLayer
from dask.utils import M, apply, get_default_shuffle_method


def align_partitions(*dfs):
    """Mutually partition and align DataFrame blocks

    This serves as precursor to multi-dataframe operations like join, concat,
    or merge.

    Parameters
    ----------
    dfs: sequence of dd.DataFrame, dd.Series and dd.base.Scalar
        Sequence of dataframes to be aligned on their index

    Returns
    -------
    dfs: sequence of dd.DataFrame, dd.Series and dd.base.Scalar
        These must have consistent divisions with each other
    divisions: tuple
        Full divisions sequence of the entire result
    result: list
        A list of lists of keys that show which data exist on which
        divisions
    """
    _is_broadcastable = partial(is_broadcastable, dfs)
    dfs1 = [df for df in dfs if isinstance(df, _Frame) and not _is_broadcastable(df)]
    if len(dfs) == 0:
        raise ValueError("dfs contains no DataFrame and Series")
    if not all(df.known_divisions for df in dfs1):
        raise ValueError(
            "Not all divisions are known, can't align "
            "partitions. Please use `set_index` "
            "to set the index."
        )

    divisions = list(unique(merge_sorted(*[df.divisions for df in dfs1])))
    if len(divisions) == 1:  # single value for index
        divisions = (divisions[0], divisions[0])
    dfs2 = [
        df.repartition(divisions, force=True) if isinstance(df, _Frame) else df
        for df in dfs
    ]

    result = list()
    inds = [0 for df in dfs]
    for d in divisions[:-1]:
        L = list()
        for i, df in enumerate(dfs2):
            if isinstance(df, _Frame):
                j = inds[i]
                divs = df.divisions
                if j < len(divs) - 1 and divs[j] == d:
                    L.append((df._name, inds[i]))
                    inds[i] += 1
                else:
                    L.append(None)
            else:  # Scalar has no divisions
                L.append(None)
        result.append(L)
    return dfs2, tuple(divisions), result


def _maybe_align_partitions(args):
    """Align DataFrame blocks if divisions are different.

    Note that if all divisions are unknown, but have equal npartitions, then
    they will be passed through unchanged. This is different than
    `align_partitions`, which will fail if divisions aren't all known"""
    _is_broadcastable = partial(is_broadcastable, args)
    dfs = [df for df in args if isinstance(df, _Frame) and not _is_broadcastable(df)]
    if not dfs:
        return args

    divisions = dfs[0].divisions
    if not all(df.divisions == divisions for df in dfs):
        dfs2 = iter(align_partitions(*dfs)[0])
        return [a if not isinstance(a, _Frame) else next(dfs2) for a in args]
    return args


def require(divisions, parts, required=None):
    """Clear out divisions where required components are not present

    In left, right, or inner joins we exclude portions of the dataset if one
    side or the other is not present.  We can achieve this at the partition
    level as well

    >>> divisions = [1, 3, 5, 7, 9]
    >>> parts = [(('a', 0), None),
    ...          (('a', 1), ('b', 0)),
    ...          (('a', 2), ('b', 1)),
    ...          (None, ('b', 2))]

    >>> divisions2, parts2 = require(divisions, parts, required=[0])
    >>> divisions2
    (1, 3, 5, 7)
    >>> parts2  # doctest: +NORMALIZE_WHITESPACE
    ((('a', 0), None),
     (('a', 1), ('b', 0)),
     (('a', 2), ('b', 1)))

    >>> divisions2, parts2 = require(divisions, parts, required=[1])
    >>> divisions2
    (3, 5, 7, 9)
    >>> parts2  # doctest: +NORMALIZE_WHITESPACE
    ((('a', 1), ('b', 0)),
     (('a', 2), ('b', 1)),
     (None, ('b', 2)))

    >>> divisions2, parts2 = require(divisions, parts, required=[0, 1])
    >>> divisions2
    (3, 5, 7)
    >>> parts2  # doctest: +NORMALIZE_WHITESPACE
    ((('a', 1), ('b', 0)),
     (('a', 2), ('b', 1)))
    """
    if not required:
        return divisions, parts
    for i in required:
        present = [j for j, p in enumerate(parts) if p[i] is not None]
        divisions = tuple(divisions[min(present) : max(present) + 2])
        parts = tuple(parts[min(present) : max(present) + 1])
    return divisions, parts


###############################################################
# Join / Merge
###############################################################


required = {
    "left": [0],
    "leftsemi": [0],
    "leftanti": [0],
    "right": [1],
    "inner": [0, 1],
    "outer": [],
}
allowed_left = ("inner", "left", "leftsemi", "leftanti")
allowed_right = ("inner", "right")


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


def merge_indexed_dataframes(lhs, rhs, left_index=True, right_index=True, **kwargs):
    """Join two partitioned dataframes along their index"""
    how = kwargs.get("how", "left")
    kwargs["left_index"] = left_index
    kwargs["right_index"] = right_index

    (lhs, rhs), divisions, parts = align_partitions(lhs, rhs)
    divisions, parts = require(divisions, parts, required[how])

    name = "join-indexed-" + tokenize(lhs, rhs, **kwargs)

    meta = lhs._meta_nonempty.merge(rhs._meta_nonempty, **kwargs)
    kwargs["result_meta"] = meta

    dsk = dict()
    for i, (a, b) in enumerate(parts):
        dsk[(name, i)] = (apply, merge_chunk, [a, b], kwargs)

    graph = HighLevelGraph.from_collections(name, dsk, dependencies=[lhs, rhs])
    return new_dd_object(graph, name, meta, divisions)


shuffle_func = shuffle  # name sometimes conflicts with keyword argument


@_deprecated_kwarg("shuffle", "shuffle_method")
def hash_join(
    lhs,
    left_on,
    rhs,
    right_on,
    how="inner",
    npartitions=None,
    suffixes=("_x", "_y"),
    shuffle_method=None,
    indicator=False,
    max_branch=None,
):
    """Join two DataFrames on particular columns with hash join

    This shuffles both datasets on the joined column and then performs an
    embarrassingly parallel join partition-by-partition

    >>> hash_join(lhs, 'id', rhs, 'id', how='left', npartitions=10)  # doctest: +SKIP
    """
    if shuffle_method is None:
        shuffle_method = get_default_shuffle_method()
    if shuffle_method == "p2p":
        from distributed.shuffle import hash_join_p2p

        return hash_join_p2p(
            lhs=lhs,
            left_on=left_on,
            rhs=rhs,
            right_on=right_on,
            how=how,
            npartitions=npartitions,
            suffixes=suffixes,
            indicator=indicator,
        )
    if npartitions is None:
        npartitions = max(lhs.npartitions, rhs.npartitions)

    lhs2 = shuffle_func(
        lhs,
        left_on,
        npartitions=npartitions,
        shuffle_method=shuffle_method,
        max_branch=max_branch,
    )
    rhs2 = shuffle_func(
        rhs,
        right_on,
        npartitions=npartitions,
        shuffle_method=shuffle_method,
        max_branch=max_branch,
    )

    if isinstance(left_on, Index):
        left_on = None
        left_index = True
    else:
        left_index = False

    if isinstance(right_on, Index):
        right_on = None
        right_index = True
    else:
        right_index = False

    kwargs = dict(
        how=how,
        left_on=left_on,
        right_on=right_on,
        left_index=left_index,
        right_index=right_index,
        suffixes=suffixes,
        indicator=indicator,
    )

    # dummy result
    # Avoid using dummy data for a collection it is empty
    _lhs_meta = lhs._meta_nonempty if len(lhs.columns) else lhs._meta
    _rhs_meta = rhs._meta_nonempty if len(rhs.columns) else rhs._meta
    meta = _lhs_meta.merge(_rhs_meta, **kwargs)

    if isinstance(left_on, list):
        left_on = (list, tuple(left_on))
    if isinstance(right_on, list):
        right_on = (list, tuple(right_on))

    kwargs["result_meta"] = meta

    joined = map_partitions(
        merge_chunk,
        lhs2,
        rhs2,
        meta=meta,
        enforce_metadata=False,
        transform_divisions=False,
        align_dataframes=False,
        **kwargs,
    )

    return joined


def single_partition_join(left, right, **kwargs):
    # if the merge is performed on_index, divisions can be kept, otherwise the
    # new index will not necessarily correspond with the current divisions

    meta = left._meta_nonempty.merge(right._meta_nonempty, **kwargs)

    use_left = kwargs.get("right_index") or right._contains_index_name(
        kwargs.get("right_on")
    )
    use_right = kwargs.get("left_index") or left._contains_index_name(
        kwargs.get("left_on")
    )

    if len(meta) == 0:
        if use_left:
            meta.index = meta.index.astype(left.index.dtype)
        elif use_right:
            meta.index = meta.index.astype(right.index.dtype)
        else:
            meta.index = meta.index.astype("int64")

    kwargs["result_meta"] = meta

    if right.npartitions == 1 and kwargs["how"] in allowed_left:
        if use_left:
            divisions = left.divisions
        elif use_right and len(right.divisions) == len(left.divisions):
            divisions = right.divisions
        else:
            divisions = [None for _ in left.divisions]

    elif left.npartitions == 1 and kwargs["how"] in allowed_right:
        if use_right:
            divisions = right.divisions
        elif use_left and len(left.divisions) == len(right.divisions):
            divisions = left.divisions
        else:
            divisions = [None for _ in right.divisions]
    else:
        raise NotImplementedError(
            "single_partition_join has no fallback for invalid calls"
        )

    joined = map_partitions(
        merge_chunk,
        left,
        right,
        meta=meta,
        enforce_metadata=False,
        transform_divisions=False,
        align_dataframes=False,
        **kwargs,
    )
    joined.divisions = tuple(divisions)
    return joined


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


@_deprecated_kwarg("shuffle", "shuffle_method")
@wraps(pd.merge)
def merge(
    left,
    right,
    how="inner",
    on=None,
    left_on=None,
    right_on=None,
    left_index=False,
    right_index=False,
    suffixes=("_x", "_y"),
    indicator=False,
    npartitions=None,
    shuffle_method=None,
    max_branch=None,
    broadcast=None,
):
    for o in [on, left_on, right_on]:
        if isinstance(o, _Frame):
            raise NotImplementedError(
                "Dask collections not currently allowed in merge columns"
            )
    if not on and not left_on and not right_on and not left_index and not right_index:
        on = [c for c in left.columns if c in right.columns]
        if not on:
            left_index = right_index = True

    if on and not left_on and not right_on:
        left_on = right_on = on
        on = None

    supported_how = ("left", "right", "outer", "inner", "leftanti", "leftsemi")
    if how not in supported_how:
        raise ValueError(
            f"dask.dataframe.merge does not support how='{how}'. Options are: {supported_how}."
            f" Note that 'leftanti' and 'leftsemi' are only dask_cudf options."
        )

    if isinstance(left, (pd.Series, pd.DataFrame)) and isinstance(
        right, (pd.Series, pd.DataFrame)
    ):
        return pd.merge(
            left,
            right,
            how=how,
            on=on,
            left_on=left_on,
            right_on=right_on,
            left_index=left_index,
            right_index=right_index,
            suffixes=suffixes,
            indicator=indicator,
        )

    # Transform pandas objects into dask.dataframe objects
    if not is_dask_collection(left):
        if right_index and left_on:  # change to join on index
            left = left.set_index(left[left_on])
            left_on = None
            left_index = True
        left = from_pandas(left, npartitions=1)  # turn into DataFrame

    if not is_dask_collection(right):
        if left_index and right_on:  # change to join on index
            right = right.set_index(right[right_on])
            right_on = None
            right_index = True
        right = from_pandas(right, npartitions=1)  # turn into DataFrame

    # Both sides are now dd.DataFrame or dd.Series objects
    merge_indexed_left = (
        left_index or left._contains_index_name(left_on)
    ) and left.known_divisions

    merge_indexed_right = (
        right_index or right._contains_index_name(right_on)
    ) and right.known_divisions

    # Both sides indexed
    if merge_indexed_left and merge_indexed_right:  # Do indexed join
        return merge_indexed_dataframes(
            left,
            right,
            how=how,
            suffixes=suffixes,
            indicator=indicator,
            left_on=left_on,
            right_on=right_on,
            left_index=left_index,
            right_index=right_index,
        )

    # Single partition on one side
    # Note that cudf supports "leftsemi" and "leftanti" joins
    elif (
        left.npartitions == 1
        and how in allowed_right
        or right.npartitions == 1
        and how in allowed_left
    ):
        return single_partition_join(
            left,
            right,
            how=how,
            right_on=right_on,
            left_on=left_on,
            left_index=left_index,
            right_index=right_index,
            suffixes=suffixes,
            indicator=indicator,
        )

    # One side is indexed, the other not
    elif (
        left_index
        and left.known_divisions
        and not right_index
        or right_index
        and right.known_divisions
        and not left_index
    ):
        left_empty = left._meta_nonempty
        right_empty = right._meta_nonempty
        meta = left_empty.merge(
            right_empty,
            how=how,
            on=on,
            left_on=left_on,
            right_on=right_on,
            left_index=left_index,
            right_index=right_index,
            suffixes=suffixes,
            indicator=indicator,
        )

        if merge_indexed_left and left.known_divisions:
            right = rearrange_by_divisions(
                right,
                right_on,
                left.divisions,
                max_branch,
                shuffle_method=shuffle_method,
            )
            left = left.clear_divisions()
        elif merge_indexed_right and right.known_divisions:
            left = rearrange_by_divisions(
                left,
                left_on,
                right.divisions,
                max_branch,
                shuffle_method=shuffle_method,
            )
            right = right.clear_divisions()
        return map_partitions(
            merge_chunk,
            left,
            right,
            meta=meta,
            how=how,
            on=on,
            left_on=left_on,
            right_on=right_on,
            left_index=left_index,
            right_index=right_index,
            suffixes=suffixes,
            indicator=indicator,
            result_meta=meta,
        )
    # Catch all hash join
    else:
        if left_on and right_on:
            warn_dtype_mismatch(left, right, left_on, right_on)

        # Check if we should use a broadcast_join
        # See note on `broadcast_bias` below.
        broadcast_bias = 0.5
        if isinstance(broadcast, float):
            broadcast_bias = float(broadcast)
            broadcast = None
        elif not isinstance(broadcast, bool) and broadcast is not None:
            # Let's be strict about the `broadcast` type to
            # avoid arbitrarily casting int to float or bool.
            raise ValueError(
                f"Optional `broadcast` argument must be float or bool."
                f"Type={type(broadcast)} is not supported."
            )
        bcast_side = "left" if left.npartitions < right.npartitions else "right"
        n_small = min(left.npartitions, right.npartitions)
        n_big = max(left.npartitions, right.npartitions)
        if (
            shuffle_method in ("tasks", "p2p", None)
            and how in ("inner", "left", "right")
            and how != bcast_side
            and broadcast is not False
        ):
            # Note on `broadcast_bias`:
            # We can expect the broadcast merge to be competitive with
            # the shuffle merge when the number of partitions in the
            # smaller collection is less than the logarithm of the number
            # of partitions in the larger collection.  By default, we add
            # a small preference for the shuffle-based merge by multiplying
            # the log result by a 0.5 scaling factor.  We call this factor
            # the `broadcast_bias`, because a larger number will make Dask
            # more likely to select the `broadcast_join` code path.  If
            # the user specifies a floating-point value for the `broadcast`
            # kwarg, that value will be used as the `broadcast_bias`.

            # FIXME: We never evaluated how P2P compares against broadcast and
            # where a suitable cutoff point is. While scaling likely still
            # depends on number of partitions, the broadcast bias should likely
            # be different
            if broadcast or (n_small < math.log2(n_big) * broadcast_bias):
                return broadcast_join(
                    left,
                    left.index if left_index else left_on,
                    right,
                    right.index if right_index else right_on,
                    how,
                    npartitions,
                    suffixes,
                    indicator=indicator,
                )

        return hash_join(
            left,
            left.index if left_index else left_on,
            right,
            right.index if right_index else right_on,
            how,
            npartitions,
            suffixes,
            shuffle_method=shuffle_method,
            indicator=indicator,
            max_branch=max_branch,
        )


###############################################################
# ASOF Join
###############################################################


def most_recent_tail(left, right):
    if len(right.index) == 0:
        return left
    return right.tail(1)


def most_recent_tail_summary(left, right, by=None):
    return pd.concat([left, right]).drop_duplicates(subset=by, keep="last")


def compute_tails(ddf, by=None):
    """For each partition, returns the last row of the most recent nonempty
    partition.
    """
    empty = ddf._meta.iloc[0:0]

    if by is None:
        return prefix_reduction(most_recent_tail, ddf, empty)
    else:
        kwargs = {"by": by}
        return prefix_reduction(most_recent_tail_summary, ddf, empty, **kwargs)


def most_recent_head(left, right):
    if len(left.index) == 0:
        return right
    return left.head(1)


def most_recent_head_summary(left, right, by=None):
    return pd.concat([left, right]).drop_duplicates(subset=by, keep="first")


def compute_heads(ddf, by=None):
    """For each partition, returns the first row of the next nonempty
    partition.
    """
    empty = ddf._meta.iloc[0:0]

    if by is None:
        return suffix_reduction(most_recent_head, ddf, empty)
    else:
        kwargs = {"by": by}
        return suffix_reduction(most_recent_head_summary, ddf, empty, **kwargs)


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


def get_unsorted_columns(frames):
    """
    Determine the unsorted column order.

    This should match the output of concat([frames], sort=False)
    """
    new_columns = pd.concat([frame._meta for frame in frames]).columns
    order = []
    for frame in frames:
        order.append(new_columns.get_indexer_for(frame.columns))

    order = np.concatenate(order)
    order = pd.unique(order)
    order = new_columns.take(order)
    return order


def merge_asof_indexed(left, right, **kwargs):
    dsk = dict()
    name = "asof-join-indexed-" + tokenize(left, right, **kwargs)
    meta = pd.merge_asof(left._meta_nonempty, right._meta_nonempty, **kwargs)

    if all(map(pd.isnull, left.divisions)):
        # results in an empty df that looks like ``meta``
        return from_pandas(meta.iloc[len(meta) :], npartitions=left.npartitions)

    if all(map(pd.isnull, right.divisions)):
        # results in an df that looks like ``left`` with nulls for
        # all ``right.columns``
        return map_partitions(
            pd.merge_asof,
            left,
            right=right,
            left_index=True,
            right_index=True,
            meta=meta,
        )

    dependencies = [left, right]
    tails = heads = None
    if kwargs["direction"] in ["backward", "nearest"]:
        tails = compute_tails(right, by=kwargs["right_by"])
        dependencies.append(tails)
    if kwargs["direction"] in ["forward", "nearest"]:
        heads = compute_heads(right, by=kwargs["right_by"])
        dependencies.append(heads)

    for i, J in enumerate(pair_partitions(left.divisions, right.divisions)):
        frames = []
        for j, lower, upper in J:
            slice = (methods.boundary_slice, (left._name, i), lower, upper, False)
            tail = (tails._name, j) if tails is not None else None
            head = (heads._name, j) if heads is not None else None
            frames.append(
                (
                    apply,
                    merge_asof_padded,
                    [slice, (right._name, j), tail, head],
                    kwargs,
                )
            )
        dsk[(name, i)] = (methods.concat, frames)

    graph = HighLevelGraph.from_collections(name, dsk, dependencies=dependencies)
    result = new_dd_object(graph, name, meta, left.divisions)
    return result


@wraps(pd.merge_asof)
def merge_asof(
    left,
    right,
    on=None,
    left_on=None,
    right_on=None,
    left_index=False,
    right_index=False,
    by=None,
    left_by=None,
    right_by=None,
    suffixes=("_x", "_y"),
    tolerance=None,
    allow_exact_matches=True,
    direction="backward",
):
    if direction not in ["backward", "forward", "nearest"]:
        raise ValueError(
            "Invalid merge_asof direction. Choose from 'backward'"
            " 'forward', or 'nearest'"
        )

    kwargs = {
        "on": on,
        "left_on": left_on,
        "right_on": right_on,
        "left_index": left_index,
        "right_index": right_index,
        "by": by,
        "left_by": left_by,
        "right_by": right_by,
        "suffixes": suffixes,
        "tolerance": tolerance,
        "allow_exact_matches": allow_exact_matches,
        "direction": direction,
    }

    if left is None or right is None:
        raise ValueError("Cannot merge_asof on None")

    # if is_dataframe_like(left) and is_dataframe_like(right):
    if isinstance(left, pd.DataFrame) and isinstance(right, pd.DataFrame):
        return pd.merge_asof(left, right, **kwargs)

    if on is not None:
        if left_on is not None or right_on is not None:
            raise ValueError(
                "Can only pass argument 'on' OR 'left_on' and 'right_on', not a "
                "combination of both."
            )
        left_on = right_on = on

    for o in [left_on, right_on]:
        if isinstance(o, _Frame):
            raise NotImplementedError(
                "Dask collections not currently allowed in merge columns"
            )

    if not is_dask_collection(left):
        left = from_pandas(left, npartitions=1)
    ixname = ixcol = divs = None
    if left_on is not None:
        if right_index:
            divs = left.divisions if left.known_divisions else None
            ixname = left.index.name
            left = left.reset_index()
            ixcol = left.columns[0]
        left = left.set_index(left_on, sorted=True)

    if not is_dask_collection(right):
        right = from_pandas(right, npartitions=1)
    if right_on is not None:
        right = right.set_index(right_on, drop=(left_on == right_on), sorted=True)

    if by is not None:
        if left_by is not None or right_by is not None:
            raise ValueError(
                "Can only pass argument 'by' OR 'left_by' and 'right_by', not a combination of both."
            )
        kwargs["left_by"] = kwargs["right_by"] = by
    if left_by is None and right_by is not None:
        raise ValueError("Must specify both left_on and right_on if one is specified.")
    if left_by is not None and right_by is None:
        raise ValueError("Must specify both left_on and right_on if one is specified.")

    del kwargs["on"], kwargs["left_on"], kwargs["right_on"], kwargs["by"]
    kwargs["left_index"] = kwargs["right_index"] = True

    if not left.known_divisions or not right.known_divisions:
        raise ValueError("merge_asof input must be sorted!")

    result = merge_asof_indexed(left, right, **kwargs)
    if left_on or right_on:
        result = result.reset_index()
        if ixcol is not None:
            if divs is not None:
                result = result.set_index(ixcol, sorted=True, divisions=divs)
            else:
                result = result.map_partitions(M.set_index, ixcol)
            result = result.map_partitions(M.rename_axis, ixname)

    return result


###############################################################
# Concat
###############################################################


def concat_and_check(dfs, ignore_order=False):
    if len(set(map(len, dfs))) != 1:
        raise ValueError("Concatenated DataFrames of different lengths")
    return methods.concat(dfs, axis=1, ignore_order=ignore_order)


def concat_unindexed_dataframes(dfs, ignore_order=False, **kwargs):
    name = "concat-" + tokenize(*dfs)

    dsk = {
        (name, i): (concat_and_check, [(df._name, i) for df in dfs], ignore_order)
        for i in range(dfs[0].npartitions)
    }
    kwargs.update({"ignore_order": ignore_order})
    meta = methods.concat([df._meta for df in dfs], axis=1, **kwargs)

    graph = HighLevelGraph.from_collections(name, dsk, dependencies=dfs)
    return new_dd_object(graph, name, meta, dfs[0].divisions)


def concat_indexed_dataframes(dfs, axis=0, join="outer", ignore_order=False, **kwargs):
    """Concatenate indexed dataframes together along the index"""
    warn = axis != 0
    kwargs.update({"ignore_order": ignore_order})
    meta = methods.concat(
        [df._meta for df in dfs],
        axis=axis,
        join=join,
        filter_warning=warn,
        **kwargs,
    )
    empties = [strip_unknown_categories(df._meta) for df in dfs]

    dfs2, divisions, parts = align_partitions(*dfs)

    name = "concat-indexed-" + tokenize(join, *dfs)

    parts2 = [
        [df if df is not None else empty for df, empty in zip(part, empties)]
        for part in parts
    ]

    filter_warning = True
    uniform = False

    dsk = {
        (name, i): (methods.concat, part, axis, join, uniform, filter_warning, kwargs)
        for i, part in enumerate(parts2)
    }
    for df in dfs2:
        dsk.update(df.dask)

    return new_dd_object(dsk, name, meta, divisions)


def stack_partitions(dfs, divisions, join="outer", ignore_order=False, **kwargs):
    """Concatenate partitions on axis=0 by doing a simple stack"""
    # Use _meta_nonempty as pandas.concat will incorrectly cast float to datetime
    # for empty data frames. See https://github.com/pandas-dev/pandas/issues/32934.

    kwargs.update({"ignore_order": ignore_order})

    meta = make_meta(
        methods.concat(
            [
                df._meta_nonempty
                for df in dfs
                if not is_dataframe_like(df) or len(df._meta_nonempty.columns) > 0
            ],
            join=join,
            filter_warning=False,
            **kwargs,
        )
    )
    empty = strip_unknown_categories(meta)

    name = f"concat-{tokenize(*dfs)}"
    dsk = {}
    i = 0
    astyped_dfs = []
    for df in dfs:
        # dtypes of all dfs need to be coherent
        # refer to https://github.com/dask/dask/issues/4685
        # and https://github.com/dask/dask/issues/5968.
        if is_dataframe_like(df):
            shared_columns = df.columns.intersection(meta.columns)
            needs_astype = [
                col
                for col in shared_columns
                if df[col].dtype != meta[col].dtype
                and not isinstance(df[col].dtype, pd.CategoricalDtype)
            ]

            if needs_astype:
                # Copy to avoid mutating the caller inplace
                df = df.copy()
                df[needs_astype] = df[needs_astype].astype(meta[needs_astype].dtypes)

        if is_series_like(df) and is_series_like(meta):
            if not df.dtype == meta.dtype and not isinstance(
                df.dtype, pd.CategoricalDtype
            ):
                df = df.astype(meta.dtype)
        else:
            pass  # TODO: there are other non-covered cases here

        astyped_dfs.append(df)

        # An error will be raised if the schemas or categories don't match. In
        # this case we need to pass along the meta object to transform each
        # partition, so they're all equivalent.
        try:
            check_meta(df._meta, meta)
            match = True
        except (ValueError, TypeError):
            match = False

        filter_warning = True
        uniform = False

        for key in df.__dask_keys__():
            if match:
                dsk[(name, i)] = key
            else:
                dsk[(name, i)] = (
                    apply,
                    methods.concat,
                    [[empty, key], 0, join, uniform, filter_warning],
                    kwargs,
                )
            i += 1

    graph = HighLevelGraph.from_collections(name, dsk, dependencies=astyped_dfs)

    return new_dd_object(graph, name, meta, divisions)


def concat(
    dfs,
    axis=0,
    join="outer",
    interleave_partitions=False,
    ignore_unknown_divisions=False,
    ignore_order=False,
    **kwargs,
):
    """Concatenate DataFrames along rows.

    - When axis=0 (default), concatenate DataFrames row-wise:

      - If all divisions are known and ordered, concatenate DataFrames keeping
        divisions. When divisions are not ordered, specifying
        interleave_partition=True allows concatenate divisions each by each.

      - If any of division is unknown, concatenate DataFrames resetting its
        division to unknown (None)

    - When axis=1, concatenate DataFrames column-wise:

      - Allowed if all divisions are known.

      - If any of division is unknown, it raises ValueError.

    Parameters
    ----------
    dfs : list
        List of dask.DataFrames to be concatenated
    axis : {0, 1, 'index', 'columns'}, default 0
        The axis to concatenate along
    join : {'inner', 'outer'}, default 'outer'
        How to handle indexes on other axis
    interleave_partitions : bool, default False
        Whether to concatenate DataFrames ignoring its order. If True, every
        divisions are concatenated each by each.
    ignore_unknown_divisions : bool, default False
        By default a warning is raised if any input has unknown divisions.
        Set to True to disable this warning.
    ignore_order : bool, default False
        Whether to ignore order when doing the union of categoricals.

    Notes
    -----
    This differs in from ``pd.concat`` in the when concatenating Categoricals
    with different categories. Pandas currently coerces those to objects
    before concatenating. Coercing to objects is very expensive for large
    arrays, so dask preserves the Categoricals by taking the union of
    the categories.

    Examples
    --------
    If all divisions are known and ordered, divisions are kept.

    >>> import dask.dataframe as dd
    >>> a                                               # doctest: +SKIP
    dd.DataFrame<x, divisions=(1, 3, 5)>
    >>> b                                               # doctest: +SKIP
    dd.DataFrame<y, divisions=(6, 8, 10)>
    >>> dd.concat([a, b])                               # doctest: +SKIP
    dd.DataFrame<concat-..., divisions=(1, 3, 6, 8, 10)>

    Unable to concatenate if divisions are not ordered.

    >>> a                                               # doctest: +SKIP
    dd.DataFrame<x, divisions=(1, 3, 5)>
    >>> b                                               # doctest: +SKIP
    dd.DataFrame<y, divisions=(2, 3, 6)>
    >>> dd.concat([a, b])                               # doctest: +SKIP
    ValueError: All inputs have known divisions which cannot be concatenated
    in order. Specify interleave_partitions=True to ignore order

    Specify interleave_partitions=True to ignore the division order.

    >>> dd.concat([a, b], interleave_partitions=True)   # doctest: +SKIP
    dd.DataFrame<concat-..., divisions=(1, 2, 3, 5, 6)>

    If any of division is unknown, the result division will be unknown

    >>> a                                               # doctest: +SKIP
    dd.DataFrame<x, divisions=(None, None)>
    >>> b                                               # doctest: +SKIP
    dd.DataFrame<y, divisions=(1, 4, 10)>
    >>> dd.concat([a, b])                               # doctest: +SKIP
    dd.DataFrame<concat-..., divisions=(None, None, None, None)>

    By default concatenating with unknown divisions will raise a warning.
    Set ``ignore_unknown_divisions=True`` to disable this:

    >>> dd.concat([a, b], ignore_unknown_divisions=True)# doctest: +SKIP
    dd.DataFrame<concat-..., divisions=(None, None, None, None)>

    Different categoricals are unioned

    >>> dd.concat([
    ...     dd.from_pandas(pd.Series(['a', 'b'], dtype='category'), 1),
    ...     dd.from_pandas(pd.Series(['a', 'c'], dtype='category'), 1),
    ... ], interleave_partitions=True).dtype
    CategoricalDtype(categories=['a', 'b', 'c'], ordered=False, categories_dtype=object)
    """

    if not isinstance(dfs, list):
        raise TypeError("dfs must be a list of DataFrames/Series objects")
    if len(dfs) == 0:
        raise ValueError("No objects to concatenate")
    if len(dfs) == 1:
        if axis == 1 and isinstance(dfs[0], Series):
            return dfs[0].to_frame()
        else:
            return dfs[0]

    if join not in ("inner", "outer"):
        raise ValueError("'join' must be 'inner' or 'outer'")

    axis = DataFrame._validate_axis(axis)

    if axis == 1:
        try:
            # remove any empty DataFrames
            dfs = [df for df in dfs if bool(len(df.columns))]
        except AttributeError:
            # 'Series' object has no attribute 'columns'
            pass
    dasks = [df for df in dfs if isinstance(df, _Frame)]
    dfs = _maybe_from_pandas(dfs)

    if axis == 1:
        if all(df.known_divisions for df in dasks):
            return concat_indexed_dataframes(
                dfs, axis=axis, join=join, ignore_order=ignore_order, **kwargs
            )
        elif (
            len(dasks) == len(dfs)
            and all(not df.known_divisions for df in dfs)
            and len({df.npartitions for df in dasks}) == 1
        ):
            if not ignore_unknown_divisions:
                warnings.warn(
                    "Concatenating dataframes with unknown divisions.\n"
                    "We're assuming that the indices of each dataframes"
                    " are \n aligned. This assumption is not generally "
                    "safe."
                )
            return concat_unindexed_dataframes(dfs, ignore_order=ignore_order, **kwargs)
        else:
            raise ValueError(
                "Unable to concatenate DataFrame with unknown "
                "division specifying axis=1"
            )
    else:
        if all(df.known_divisions for df in dasks):
            # each DataFrame's division must be greater than previous one
            if all(
                dfs[i].divisions[-1] < dfs[i + 1].divisions[0]
                for i in range(len(dfs) - 1)
            ):
                divisions = []
                for df in dfs[:-1]:
                    # remove last to concatenate with next
                    divisions += df.divisions[:-1]
                divisions += dfs[-1].divisions
                return stack_partitions(
                    dfs, divisions, join=join, ignore_order=ignore_order, **kwargs
                )
            elif interleave_partitions:
                return concat_indexed_dataframes(
                    dfs, join=join, ignore_order=ignore_order, **kwargs
                )
            else:
                divisions = [None] * (sum(df.npartitions for df in dfs) + 1)
                return stack_partitions(
                    dfs, divisions, join=join, ignore_order=ignore_order, **kwargs
                )
        else:
            divisions = [None] * (sum(df.npartitions for df in dfs) + 1)
            return stack_partitions(
                dfs, divisions, join=join, ignore_order=ignore_order, **kwargs
            )


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
    if not isinstance(on, _Frame):
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


@_deprecated_kwarg("shuffle", "shuffle_method")
def broadcast_join(
    lhs,
    left_on,
    rhs,
    right_on,
    how="inner",
    npartitions=None,
    suffixes=("_x", "_y"),
    shuffle_method=None,
    indicator=False,
    parts_out=None,
):
    """Join two DataFrames on particular columns by broadcasting

    This broadcasts the partitions of the smaller DataFrame to each
    partition of the larger DataFrame, joins each partition pair,
    and then concatenates the new data for each output partition.
    """

    if npartitions:
        # Repartition the larger collection before the merge
        if lhs.npartitions < rhs.npartitions:
            rhs = rhs.repartition(npartitions=npartitions)
        else:
            lhs = lhs.repartition(npartitions=npartitions)

    if how not in ("inner", "left", "right"):
        # Broadcast algorithm cannot handle an "outer" join
        raise ValueError(
            "Only 'inner', 'left' and 'right' broadcast joins are supported."
        )

    if how == "left" and lhs.npartitions < rhs.npartitions:
        # Must broadcast rhs for a "left" broadcast join
        raise ValueError("'left' broadcast join requires rhs broadcast.")

    if how == "right" and rhs.npartitions <= lhs.npartitions:
        # Must broadcast lhs for a "right" broadcast join
        raise ValueError("'right' broadcast join requires lhs broadcast.")

    # TODO: It *may* be beneficial to perform the hash
    # split for "inner" join as well (even if it is not
    # technically needed for correctness).  More testing
    # is needed here.
    if how != "inner":
        # Shuffle to-be-broadcasted side by hash. This
        # means that we will need to perform a local
        # shuffle and split on each partition of the
        # "other" collection (with the same hashing
        # approach) to ensure the correct rows are
        # joined by `merge_chunk`.  The local hash and
        # split of lhs is in `_split_partition`.
        if lhs.npartitions < rhs.npartitions:
            lhs2 = shuffle_func(
                lhs,
                left_on,
                shuffle_method="tasks",
            )
            lhs_name = lhs2._name
            lhs_dep = lhs2
            rhs_name = rhs._name
            rhs_dep = rhs
        else:
            rhs2 = shuffle_func(
                rhs,
                right_on,
                shuffle_method="tasks",
            )
            lhs_name = lhs._name
            lhs_dep = lhs
            rhs_name = rhs2._name
            rhs_dep = rhs2
    else:
        lhs_name = lhs._name
        lhs_dep = lhs
        rhs_name = rhs._name
        rhs_dep = rhs

    if isinstance(left_on, Index):
        left_on = None
        left_index = True
    else:
        left_index = False

    if isinstance(right_on, Index):
        right_on = None
        right_index = True
    else:
        right_index = False

    merge_kwargs = dict(
        how=how,
        left_on=left_on,
        right_on=right_on,
        left_index=left_index,
        right_index=right_index,
        suffixes=suffixes,
        indicator=indicator,
    )

    # dummy result
    meta = lhs._meta_nonempty.merge(rhs._meta_nonempty, **merge_kwargs)
    merge_kwargs["result_meta"] = meta

    # Assume the output partitions/divisions
    # should correspond to the collection that
    # is NOT broadcasted.
    if lhs.npartitions < rhs.npartitions:
        npartitions = rhs.npartitions
        divisions = rhs.divisions
        _index_names = set(rhs._meta_nonempty.index.names)
    else:
        npartitions = lhs.npartitions
        divisions = lhs.divisions
        _index_names = set(lhs._meta_nonempty.index.names)

    # Cannot preserve divisions if the index is lost
    if _index_names != set(meta.index.names):
        divisions = [None] * (npartitions + 1)

    token = tokenize(lhs, rhs, npartitions, **merge_kwargs)
    name = "bcast-join-" + token
    broadcast_join_layer = BroadcastJoinLayer(
        name,
        npartitions,
        lhs_name,
        lhs.npartitions,
        rhs_name,
        rhs.npartitions,
        parts_out=parts_out,
        **merge_kwargs,
    )

    graph = HighLevelGraph.from_collections(
        name,
        broadcast_join_layer,
        dependencies=[lhs_dep, rhs_dep],
    )

    return new_dd_object(graph, name, meta, divisions)


def _recursive_pairwise_outer_join(
    dataframes_to_merge, on, lsuffix, rsuffix, npartitions, shuffle_method
):
    """
    Schedule the merging of a list of dataframes in a pairwise method. This is a recursive function that results
    in a much more efficient scheduling of merges than a simple loop
    from:
    [A] [B] [C] [D] -> [AB] [C] [D] -> [ABC] [D] -> [ABCD]
    to:
    [A] [B] [C] [D] -> [AB] [CD] -> [ABCD]
    Note that either way, n-1 merges are still required, but using a pairwise reduction it can be completed in parallel.
    :param dataframes_to_merge: A list of Dask dataframes to be merged together on their index
    :return: A single Dask Dataframe, comprised of the pairwise-merges of all provided dataframes
    """
    number_of_dataframes_to_merge = len(dataframes_to_merge)

    merge_options = {
        "on": on,
        "lsuffix": lsuffix,
        "rsuffix": rsuffix,
        "npartitions": npartitions,
        "shuffle_method": shuffle_method,
    }

    # Base case 1: just return the provided dataframe and merge with `left`
    if number_of_dataframes_to_merge == 1:
        return dataframes_to_merge[0]

    # Base case 2: merge the two provided dataframe to be merged with `left`
    if number_of_dataframes_to_merge == 2:
        merged_ddf = dataframes_to_merge[0].join(
            dataframes_to_merge[1], how="outer", **merge_options
        )
        return merged_ddf

    # Recursive case: split the list of dfs into two ~even sizes and continue down
    else:
        middle_index = number_of_dataframes_to_merge // 2
        merged_ddf = _recursive_pairwise_outer_join(
            [
                _recursive_pairwise_outer_join(
                    dataframes_to_merge[:middle_index], **merge_options
                ),
                _recursive_pairwise_outer_join(
                    dataframes_to_merge[middle_index:], **merge_options
                ),
            ],
            **merge_options,
        )
        return merged_ddf
