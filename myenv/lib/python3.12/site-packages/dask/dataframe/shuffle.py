from __future__ import annotations

import contextlib
import logging
import math
import shutil
import tempfile
import uuid
import warnings
from collections.abc import Callable, Mapping, Sequence
from typing import Any, Literal

import numpy as np
import pandas as pd
import tlz as toolz
from pandas.api.types import is_numeric_dtype

from dask import config
from dask.base import compute, compute_as_if_collection, is_dask_collection, tokenize
from dask.dataframe import methods
from dask.dataframe._compat import PANDAS_GE_300
from dask.dataframe.core import (
    DataFrame,
    Series,
    _deprecated_kwarg,
    _Frame,
    map_partitions,
    new_dd_object,
)
from dask.dataframe.dispatch import (
    group_split_dispatch,
    hash_object_dispatch,
    partd_encode_dispatch,
)
from dask.dataframe.utils import UNKNOWN_CATEGORIES
from dask.highlevelgraph import HighLevelGraph
from dask.layers import ShuffleLayer, SimpleShuffleLayer
from dask.sizeof import sizeof
from dask.utils import M, digit, get_default_shuffle_method

logger = logging.getLogger(__name__)


def _calculate_divisions(
    df: DataFrame,
    partition_col: Series,
    repartition: bool,
    npartitions: int,
    upsample: float = 1.0,
    partition_size: float = 128e6,
    ascending: bool = True,
) -> tuple[list, list, list, bool]:
    """
    Utility function to calculate divisions for calls to `map_partitions`
    """
    sizes = df.map_partitions(sizeof) if repartition else []
    divisions = partition_col._repartition_quantiles(npartitions, upsample=upsample)
    mins = partition_col.map_partitions(M.min)
    maxes = partition_col.map_partitions(M.max)

    try:
        divisions, sizes, mins, maxes = compute(divisions, sizes, mins, maxes)
    except TypeError as e:
        # When there are nulls and a column is non-numeric, a TypeError is sometimes raised as a result of
        # 1) computing mins/maxes above, 2) every null being switched to NaN, and 3) NaN being a float.
        # Also, Pandas ExtensionDtypes may cause TypeErrors when dealing with special nulls such as pd.NaT or pd.NA.
        # If this happens, we hint the user about eliminating nulls beforehand.
        if not is_numeric_dtype(partition_col.dtype):
            obj, suggested_method = (
                ("column", f"`.dropna(subset=['{partition_col.name}'])`")
                if any(partition_col._name == df[c]._name for c in df)
                else ("series", "`.loc[series[~series.isna()]]`")
            )
            raise NotImplementedError(
                f"Divisions calculation failed for non-numeric {obj} '{partition_col.name}'.\n"
                f"This is probably due to the presence of nulls, which Dask does not entirely support in the index.\n"
                f"We suggest you try with {suggested_method}."
            ) from e
        # For numeric types there shouldn't be problems with nulls, so we raise as-it-is this particular TypeError
        else:
            raise e

    empty_dataframe_detected = pd.isna(divisions).all()
    if repartition or empty_dataframe_detected:
        total = sum(sizes)
        npartitions = max(math.ceil(total / partition_size), 1)
        npartitions = min(npartitions, df.npartitions)
        n = divisions.size
        try:
            divisions = np.interp(
                x=np.linspace(0, n - 1, npartitions + 1),
                xp=np.linspace(0, n - 1, n),
                fp=divisions.tolist(),
            ).tolist()
        except (TypeError, ValueError):  # str type
            indexes = np.linspace(0, n - 1, npartitions + 1).astype(int)
            divisions = divisions.iloc[indexes].tolist()
    else:
        # Drop duplicate divisions returned by partition quantiles
        n = divisions.size
        divisions = (
            list(divisions.iloc[: n - 1].unique()) + divisions.iloc[n - 1 :].tolist()
        )

    mins = mins.bfill()
    maxes = maxes.bfill()
    if isinstance(partition_col.dtype, pd.CategoricalDtype):
        dtype = partition_col.dtype
        mins = mins.astype(dtype)
        maxes = maxes.astype(dtype)

    if mins.isna().any() or maxes.isna().any():
        presorted = False
    else:
        n = mins.size
        maxes2 = (maxes.iloc[: n - 1] if ascending else maxes.iloc[1:]).reset_index(
            drop=True
        )
        mins2 = (mins.iloc[1:] if ascending else mins.iloc[: n - 1]).reset_index(
            drop=True
        )
        presorted = (
            mins.tolist() == mins.sort_values(ascending=ascending).tolist()
            and maxes.tolist() == maxes.sort_values(ascending=ascending).tolist()
            and (maxes2 < mins2).all()
        )

    return divisions, mins.tolist(), maxes.tolist(), presorted


@_deprecated_kwarg("shuffle", "shuffle_method")
def sort_values(
    df: DataFrame,
    by: str | list[str],
    npartitions: int | Literal["auto"] | None = None,
    shuffle_method: str | None = None,
    ascending: bool | list[bool] = True,
    na_position: Literal["first"] | Literal["last"] = "last",
    upsample: float = 1.0,
    partition_size: float = 128e6,
    sort_function: Callable[[pd.DataFrame], pd.DataFrame] | None = None,
    sort_function_kwargs: Mapping[str, Any] | None = None,
) -> DataFrame:
    """See DataFrame.sort_values for docstring"""
    if na_position not in ("first", "last"):
        raise ValueError("na_position must be either 'first' or 'last'")
    if not isinstance(by, list):
        by = [by]
    if any(not isinstance(b, str) for b in by):
        raise NotImplementedError(
            "Dataframes only support sorting by named columns which must be passed as a "
            "string or a list of strings.\n"
            "You passed %s" % str(by)
        )

    if (
        ascending is not None
        and not isinstance(ascending, bool)
        and not len(ascending) == len(by)
    ):
        raise ValueError(f"Length of {ascending=} != length of {by=}")

    sort_kwargs = {
        "by": by,
        "ascending": ascending,
        "na_position": na_position,
    }
    if sort_function is None:
        sort_function = M.sort_values
    if sort_function_kwargs is not None:
        sort_kwargs.update(sort_function_kwargs)

    if df.npartitions == 1:
        return df.map_partitions(sort_function, **sort_kwargs)

    if npartitions == "auto":
        warnings.warn(
            "`npartitions='auto'` is deprecated, either set it as an integer or leave as `None`.",
            FutureWarning,
            2,
        )
        repartition = True
        npartitions = max(100, df.npartitions)
    else:
        if npartitions is None:
            npartitions = df.npartitions
        repartition = False

    sort_by_col = df[by[0]]
    divisions_ascending = ascending
    if divisions_ascending and not isinstance(divisions_ascending, bool):
        divisions_ascending = divisions_ascending[0]
    assert divisions_ascending is None or isinstance(divisions_ascending, bool)
    divisions, _, _, presorted = _calculate_divisions(
        df,
        sort_by_col,
        repartition,
        npartitions,
        upsample,
        partition_size,
        divisions_ascending,
    )

    if len(divisions) == 2:
        return df.repartition(npartitions=1).map_partitions(
            sort_function, **sort_kwargs
        )

    if presorted and npartitions == df.npartitions:
        # divisions are in the right place
        return df.map_partitions(sort_function, **sort_kwargs)

    df = rearrange_by_divisions(
        df,
        by[0],
        divisions,
        shuffle_method=shuffle_method,
        ascending=divisions_ascending,
        na_position=na_position,
        duplicates=False,
    )
    df = df.map_partitions(sort_function, **sort_kwargs)
    return df


@_deprecated_kwarg("compute")
@_deprecated_kwarg("shuffle", "shuffle_method")
def set_index(
    df: DataFrame,
    index: str | Series,
    npartitions: int | Literal["auto"] | None = None,
    shuffle_method: str | None = None,
    compute: bool = False,
    drop: bool = True,
    upsample: float = 1.0,
    divisions: Sequence | None = None,
    partition_size: float = 128e6,
    sort: bool = True,
    **kwargs,
) -> DataFrame:
    """See _Frame.set_index for docstring"""
    if not sort:
        return df.map_partitions(
            M.set_index, index, align_dataframes=False, drop=drop, **kwargs
        ).clear_divisions()

    if npartitions == "auto":
        warnings.warn(
            "`npartitions='auto'` is deprecated, either set it as an integer or leave as `None`.",
            FutureWarning,
            2,
        )
        repartition = True
        npartitions = max(100, df.npartitions)
    else:
        if npartitions is None:
            npartitions = df.npartitions
        repartition = False

    if not isinstance(index, Series):
        index2 = df[index]
    else:
        index2 = index

    if divisions is None:
        divisions, mins, maxes, presorted = _calculate_divisions(
            df, index2, repartition, npartitions, upsample, partition_size
        )

        if presorted and npartitions == df.npartitions:
            divisions = mins + [maxes[-1]]
            result = set_sorted_index(df, index, drop=drop, divisions=divisions)
            return result.map_partitions(M.sort_index)

    return set_partition(
        df,
        index,
        divisions,
        shuffle_method=shuffle_method,
        drop=drop,
        compute=compute,
        **kwargs,
    )


@_deprecated_kwarg("shuffle", "shuffle_method")
def set_partition(
    df: DataFrame,
    index: str | Series,
    divisions: Sequence,
    max_branch: int = 32,
    drop: bool = True,
    shuffle_method: str | None = None,
    compute: bool | None = None,
) -> DataFrame:
    """Group DataFrame by index

    Sets a new index and partitions data along that index according to
    divisions.  Divisions are often found by computing approximate quantiles.
    The function ``set_index`` will do both of these steps.

    Parameters
    ----------
    df: DataFrame/Series
        Data that we want to re-partition
    index: string or Series
        Column to become the new index
    divisions: list
        Values to form new divisions between partitions
    drop: bool, default True
        Whether to delete columns to be used as the new index
    shuffle_method: str (optional)
        Either 'disk' for an on-disk shuffle or 'tasks' to use the task
        scheduling framework.  Use 'disk' if you are on a single machine
        and 'tasks' if you are on a distributed cluster.
    max_branch: int (optional)
        If using the task-based shuffle, the amount of splitting each
        partition undergoes.  Increase this for fewer copies but more
        scheduler overhead.

    See Also
    --------
    set_index
    shuffle
    partd
    """
    if isinstance(divisions, tuple):
        # pd.isna considers tuples to be scalars. Convert to a list.
        divisions = list(divisions)

    if not isinstance(index, Series):
        dtype = df[index].dtype
    else:
        dtype = index.dtype

    if pd.isna(divisions).any() and pd.api.types.is_integer_dtype(dtype):
        # Can't construct a Series[int64] when any / all of the divisions are NaN.
        divisions = df._meta._constructor_sliced(divisions)
    elif (
        isinstance(dtype, pd.CategoricalDtype)
        and UNKNOWN_CATEGORIES in dtype.categories
    ):
        # If categories are unknown, leave as a string dtype instead.
        divisions = df._meta._constructor_sliced(divisions)
    else:
        divisions = df._meta._constructor_sliced(divisions, dtype=dtype)

    meta = df._meta._constructor_sliced([0])
    # Ensure that we have the same index as before to avoid alignment
    # when calculating meta dtypes later on
    meta.index = df._meta_nonempty.index[:1]

    if not isinstance(index, Series):
        partitions = df[index].map_partitions(
            set_partitions_pre, divisions=divisions, meta=meta
        )
        df2 = df.assign(_partitions=partitions)
    else:
        partitions = index.map_partitions(
            set_partitions_pre, divisions=divisions, meta=meta
        )
        df2 = df.assign(_partitions=partitions, _index=index)

    df3 = rearrange_by_column(
        df2,
        "_partitions",
        max_branch=max_branch,
        npartitions=len(divisions) - 1,
        shuffle_method=shuffle_method,
        compute=compute,
        ignore_index=True,
    )

    if not isinstance(index, Series):
        df4 = df3.map_partitions(
            set_index_post_scalar,
            index_name=index,
            drop=drop,
            column_dtype=df.columns.dtype,
        )
    else:
        df4 = df3.map_partitions(
            set_index_post_series,
            index_name=index.name,
            drop=drop,
            column_dtype=df.columns.dtype,
        )

    divisions = methods.tolist(divisions)
    # None and pd.NA values are not sortable
    df4.divisions = tuple(i if not pd.isna(i) else np.nan for i in divisions)

    return df4.map_partitions(M.sort_index)


@_deprecated_kwarg("shuffle", "shuffle_method")
def shuffle(
    df,
    index,
    shuffle_method=None,
    npartitions=None,
    max_branch=32,
    ignore_index=False,
    compute=None,
):
    """Group DataFrame by index

    Hash grouping of elements. After this operation all elements that have
    the same index will be in the same partition. Note that this requires
    full dataset read, serialization and shuffle. This is expensive. If
    possible you should avoid shuffles.

    This does not preserve a meaningful index/partitioning scheme. This is not
    deterministic if done in parallel.

    See Also
    --------
    set_index
    set_partition
    shuffle_disk
    """
    list_like = pd.api.types.is_list_like(index) and not is_dask_collection(index)
    shuffle_method = shuffle_method or get_default_shuffle_method()

    if not isinstance(index, _Frame):
        if list_like:
            # Make sure we don't try to select with pd.Series/pd.Index
            index = list(index)
        index = df._select_columns_or_index(index)
    elif hasattr(index, "to_frame"):
        # If this is an index, we should still convert to a
        # DataFrame. Otherwise, the hashed values of a column
        # selection will not match (important when merging).
        index = index.to_frame()

    dtypes = {}
    for col, dtype in index.dtypes.items():
        if pd.api.types.is_numeric_dtype(dtype):
            dtypes[col] = np.float64
    if not dtypes:
        dtypes = None

    meta = df._meta._constructor_sliced([0])
    # Ensure that we have the same index as before to avoid alignment
    # when calculating meta dtypes later on
    meta.index = df._meta_nonempty.index[:1]
    partitions = index.map_partitions(
        partitioning_index,
        npartitions=npartitions or df.npartitions,
        meta=meta,
        transform_divisions=False,
        cast_dtype=dtypes,
    )
    df2 = df.assign(_partitions=partitions)
    df2._meta.index.name = df._meta.index.name
    df3 = rearrange_by_column(
        df2,
        "_partitions",
        npartitions=npartitions,
        max_branch=max_branch,
        shuffle_method=shuffle_method,
        compute=compute,
        ignore_index=ignore_index,
    )
    del df3["_partitions"]
    return df3


@_deprecated_kwarg("shuffle", "shuffle_method")
def rearrange_by_divisions(
    df,
    column,
    divisions,
    max_branch=None,
    shuffle_method=None,
    ascending=True,
    na_position="last",
    duplicates=True,
):
    """Shuffle dataframe so that column separates along divisions"""
    divisions = df._meta._constructor_sliced(divisions)
    # duplicates need to be removed sometimes to properly sort null dataframes
    if not duplicates:
        divisions = divisions.drop_duplicates()
    meta = df._meta._constructor_sliced([0])
    # Ensure that we have the same index as before to avoid alignment
    # when calculating meta dtypes later on
    meta.index = df._meta_nonempty.index[:1]
    # Assign target output partitions to every row
    partitions = df[column].map_partitions(
        set_partitions_pre,
        divisions=divisions,
        ascending=ascending,
        na_position=na_position,
        meta=meta,
    )
    df2 = df.assign(_partitions=partitions)

    # Perform shuffle
    df3 = rearrange_by_column(
        df2,
        "_partitions",
        max_branch=max_branch,
        npartitions=len(divisions) - 1,
        shuffle_method=shuffle_method,
    )
    del df3["_partitions"]
    return df3


@_deprecated_kwarg("shuffle", "shuffle_method")
def rearrange_by_column(
    df,
    col,
    npartitions=None,
    max_branch=None,
    shuffle_method=None,
    compute=None,
    ignore_index=False,
):
    shuffle_method = shuffle_method or get_default_shuffle_method()

    # if the requested output partitions < input partitions
    # we repartition first as shuffling overhead is
    # proportionate to the number of input partitions

    if (
        shuffle_method != "p2p"
        and npartitions is not None
        and npartitions < df.npartitions
    ):
        df = df.repartition(npartitions=npartitions)

    if shuffle_method == "disk":
        return rearrange_by_column_disk(df, col, npartitions, compute=compute)
    elif shuffle_method == "tasks":
        df2 = rearrange_by_column_tasks(
            df, col, max_branch, npartitions, ignore_index=ignore_index
        )
        if ignore_index:
            df2._meta = df2._meta.reset_index(drop=True)
        return df2
    elif shuffle_method == "p2p":
        from distributed.shuffle import rearrange_by_column_p2p

        return rearrange_by_column_p2p(df, col, npartitions)
    else:
        raise NotImplementedError("Unknown shuffle method %s" % shuffle_method)


class maybe_buffered_partd:
    """
    If serialized, will return non-buffered partd. Otherwise returns a buffered partd
    """

    def __init__(self, encode_cls=None, buffer=True, tempdir=None):
        self.tempdir = tempdir or config.get("temporary_directory", None)
        self.buffer = buffer
        self.compression = config.get("dataframe.shuffle.compression", None)
        self.encode_cls = encode_cls
        if encode_cls is None:
            import partd

            self.encode_cls = partd.PandasBlocks

    def __reduce__(self):
        if self.tempdir:
            return (maybe_buffered_partd, (self.encode_cls, False, self.tempdir))
        else:
            return (maybe_buffered_partd, (self.encode_cls, False))

    def __call__(self, *args, **kwargs):
        import partd

        path = tempfile.mkdtemp(suffix=".partd", dir=self.tempdir)

        try:
            partd_compression = (
                getattr(partd.compressed, self.compression)
                if self.compression
                else None
            )
        except AttributeError as e:
            raise ImportError(
                "Not able to import and load {} as compression algorithm."
                "Please check if the library is installed and supported by Partd.".format(
                    self.compression
                )
            ) from e
        file = partd.File(path)
        partd.file.cleanup_files.append(path)
        # Envelope partd file with compression, if set and available
        if partd_compression:
            file = partd_compression(file)
        if self.buffer:
            return self.encode_cls(partd.Buffer(partd.Dict(), file))
        else:
            return self.encode_cls(file)


def rearrange_by_column_disk(df, column, npartitions=None, compute=False):
    """Shuffle using local disk

    See Also
    --------
    rearrange_by_column_tasks:
        Same function, but using tasks rather than partd
        Has a more informative docstring
    """
    if npartitions is None:
        npartitions = df.npartitions

    token = tokenize(df, column, npartitions)
    always_new_token = uuid.uuid1().hex

    p = ("zpartd-" + always_new_token,)
    encode_cls = partd_encode_dispatch(df._meta)
    dsk1 = {p: (maybe_buffered_partd(encode_cls=encode_cls),)}

    # Partition data on disk
    name = "shuffle-partition-" + always_new_token
    dsk2 = {
        (name, i): (shuffle_group_3, key, column, npartitions, p)
        for i, key in enumerate(df.__dask_keys__())
    }

    dependencies = []
    if compute:
        graph = HighLevelGraph.merge(df.dask, dsk1, dsk2)
        graph = HighLevelGraph.from_collections(name, graph, dependencies=[df])
        keys = [p, sorted(dsk2)]
        pp, values = compute_as_if_collection(DataFrame, graph, keys)
        dsk1 = {p: pp}
        dsk2 = dict(zip(sorted(dsk2), values))
    else:
        dependencies.append(df)

    # Barrier
    barrier_token = "barrier-" + always_new_token
    dsk3 = {barrier_token: (barrier, list(dsk2))}

    # Collect groups
    name = "shuffle-collect-" + token
    dsk4 = {
        (name, i): (collect, p, i, df._meta, barrier_token) for i in range(npartitions)
    }

    divisions = (None,) * (npartitions + 1)

    layer = toolz.merge(dsk1, dsk2, dsk3, dsk4)
    graph = HighLevelGraph.from_collections(name, layer, dependencies=dependencies)
    return new_dd_object(graph, name, df._meta, divisions)


def _noop(x, cleanup_token):
    """
    A task that does nothing.
    """
    return x


def rearrange_by_column_tasks(
    df, column, max_branch=32, npartitions=None, ignore_index=False
):
    """Order divisions of DataFrame so that all values within column(s) align

    This enacts a task-based shuffle.  It contains most of the tricky logic
    around the complex network of tasks.  Typically before this function is
    called a new column, ``"_partitions"`` has been added to the dataframe,
    containing the output partition number of every row.  This function
    produces a new dataframe where every row is in the proper partition.  It
    accomplishes this by splitting each input partition into several pieces,
    and then concatenating pieces from different input partitions into output
    partitions.  If there are enough partitions then it does this work in
    stages to avoid scheduling overhead.

    Lets explain the motivation for this further.  Imagine that we have 1000
    input partitions and 1000 output partitions. In theory we could split each
    input into 1000 pieces, and then move the 1 000 000 resulting pieces
    around, and then concatenate them all into 1000 output groups.  This would
    be fine, but the central scheduling overhead of 1 000 000 tasks would
    become a bottleneck.  Instead we do this in stages so that we split each of
    the 1000 inputs into 30 pieces (we now have 30 000 pieces) move those
    around, concatenate back down to 1000, and then do the same process again.
    This has the same result as the full transfer, but now we've moved data
    twice (expensive) but done so with only 60 000 tasks (cheap).

    Note that the `column` input may correspond to a list of columns (rather
    than just a single column name).  In this case, the `shuffle_group` and
    `shuffle_group_2` functions will use hashing to map each row to an output
    partition. This approach may require the same rows to be hased multiple
    times, but avoids the need to assign a new "_partitions" column.

    Parameters
    ----------
    df: dask.dataframe.DataFrame
    column: str or list
        A column name on which we want to split, commonly ``"_partitions"``
        which is assigned by functions upstream.  This could also be a list of
        columns (in which case shuffle_group will create a hash array/column).
    max_branch: int
        The maximum number of splits per input partition.  Defaults to 32.
        If there are more partitions than this then the shuffling will occur in
        stages in order to avoid creating npartitions**2 tasks
        Increasing this number increases scheduling overhead but decreases the
        number of full-dataset transfers that we have to make.
    npartitions: Optional[int]
        The desired number of output partitions

    Returns
    -------
    df3: dask.dataframe.DataFrame

    See also
    --------
    rearrange_by_column_disk: same operation, but uses partd
    rearrange_by_column: parent function that calls this or rearrange_by_column_disk
    shuffle_group: does the actual splitting per-partition
    """

    max_branch = max_branch or 32

    if (npartitions or df.npartitions) <= max_branch:
        # We are creating a small number of output partitions.
        # No need for staged shuffling. Staged shuffling will
        # sometimes require extra work/communication in this case.
        token = tokenize(df, column, npartitions)
        shuffle_name = f"simple-shuffle-{token}"
        npartitions = npartitions or df.npartitions
        shuffle_layer = SimpleShuffleLayer(
            shuffle_name,
            column,
            npartitions,
            df.npartitions,
            ignore_index,
            df._name,
            df._meta,
        )
        graph = HighLevelGraph.from_collections(
            shuffle_name, shuffle_layer, dependencies=[df]
        )
        return new_dd_object(graph, shuffle_name, df._meta, [None] * (npartitions + 1))

    n = df.npartitions
    stages = int(math.ceil(math.log(n) / math.log(max_branch)))
    if stages > 1:
        k = int(math.ceil(n ** (1 / stages)))
    else:
        k = n

    inputs = [tuple(digit(i, j, k) for j in range(stages)) for i in range(k**stages)]

    npartitions_orig = df.npartitions
    token = tokenize(df, stages, column, n, k)
    for stage in range(stages):
        stage_name = f"shuffle-{stage}-{token}"
        stage_layer = ShuffleLayer(
            stage_name,
            column,
            inputs,
            stage,
            npartitions,
            n,
            k,
            ignore_index,
            df._name,
            df._meta,
        )
        graph = HighLevelGraph.from_collections(
            stage_name, stage_layer, dependencies=[df]
        )
        df = new_dd_object(graph, stage_name, df._meta, df.divisions)

    if npartitions is not None and npartitions != npartitions_orig:
        token = tokenize(df, npartitions)
        repartition_group_token = "repartition-group-" + token

        dsk = {
            (repartition_group_token, i): (
                shuffle_group_2,
                k,
                column,
                ignore_index,
                npartitions,
            )
            for i, k in enumerate(df.__dask_keys__())
        }

        repartition_get_name = "repartition-get-" + token

        for p in range(npartitions):
            dsk[(repartition_get_name, p)] = (
                shuffle_group_get,
                (repartition_group_token, p % npartitions_orig),
                p,
            )

        graph2 = HighLevelGraph.from_collections(
            repartition_get_name, dsk, dependencies=[df]
        )
        df2 = new_dd_object(
            graph2, repartition_get_name, df._meta, [None] * (npartitions + 1)
        )
    else:
        df2 = df
        df2.divisions = (None,) * (npartitions_orig + 1)

    return df2


########################################################
# Various convenience functions to be run by the above #
########################################################


def partitioning_index(df, npartitions, cast_dtype=None):
    """
    Computes a deterministic index mapping each record to a partition.

    Identical rows are mapped to the same partition.

    Parameters
    ----------
    df : DataFrame/Series/Index
    npartitions : int
        The number of partitions to group into.
    cast_dtype : dtype, optional
        The dtype to cast to to avoid nullability issues

    Returns
    -------
    partitions : ndarray
        An array of int64 values mapping each record to a partition.
    """
    if cast_dtype is not None:
        # Fixme: astype raises with strings in numeric columns, but raising
        # here might be very noisy
        df = df.astype(cast_dtype, errors="ignore")
    res = hash_object_dispatch(df, index=False) % int(npartitions)
    # Note: Use a signed integer since pandas is more efficient at handling
    # this since there is not always a fastpath for uints
    return res.astype(np.min_scalar_type(-(npartitions - 1)))


def barrier(args):
    list(args)
    return 0


def cleanup_partd_files(p, keys):
    """
    Cleanup the files in a partd.File dataset.

    Parameters
    ----------
    p : partd.Interface
        File or Encode wrapping a file should be OK.
    keys: List
        Just for scheduling purposes, not actually used.
    """
    import partd

    if isinstance(p, partd.Encode):
        maybe_file = p.partd
    else:
        maybe_file = None

    if isinstance(maybe_file, partd.File):
        path = maybe_file.path
    else:
        path = None

    if path:
        shutil.rmtree(path, ignore_errors=True)


def collect(p, part, meta, barrier_token):
    """Collect partitions from partd, yield dataframes"""
    with ensure_cleanup_on_exception(p):
        res = p.get(part)
        return res if len(res) > 0 else meta


def set_partitions_pre(s, divisions, ascending=True, na_position="last"):
    try:
        if ascending:
            partitions = divisions.searchsorted(s, side="right") - 1
        else:
            partitions = len(divisions) - divisions.searchsorted(s, side="right") - 1
    except (TypeError, ValueError):
        # `searchsorted` fails if either `divisions` or `s` contains nulls and strings
        partitions = np.empty(len(s), dtype="int32")
        not_null = s.notna()
        divisions_notna = divisions[divisions.notna()]
        if ascending:
            partitions[not_null] = (
                divisions_notna.searchsorted(s[not_null], side="right") - 1
            )
        else:
            partitions[not_null] = (
                len(divisions)
                - divisions_notna.searchsorted(s[not_null], side="right")
                - 1
            )
    partitions[(partitions < 0) | (partitions >= len(divisions) - 1)] = (
        len(divisions) - 2 if ascending else 0
    )
    nas = s.isna()
    # We could be a ndarray already (datetime dtype)
    nas = getattr(nas, "values", nas)
    partitions[nas] = len(divisions) - 2 if na_position == "last" else 0
    return partitions


def shuffle_group_2(df, cols, ignore_index, nparts):
    if not len(df):
        return {}, df

    if isinstance(cols, str):
        cols = [cols]

    if cols and cols[0] == "_partitions":
        ind = df[cols[0]].astype(np.int32)
    else:
        ind = (
            hash_object_dispatch(df[cols] if cols else df, index=False) % int(nparts)
        ).astype(np.int32)

    n = ind.max() + 1
    result2 = group_split_dispatch(df, ind, n, ignore_index=ignore_index)
    return result2, df.iloc[:0]


def shuffle_group_get(g_head, i):
    g, head = g_head
    if i in g:
        return g[i]
    else:
        return head


def shuffle_group(df, cols, stage, k, npartitions, ignore_index, nfinal):
    """Splits dataframe into groups

    The group is determined by their final partition, and which stage we are in
    in the shuffle

    Parameters
    ----------
    df: DataFrame
    cols: str or list
        Column name(s) on which to split the dataframe. If ``cols`` is not
        "_partitions", hashing will be used to determine target partition
    stage: int
        We shuffle dataframes with many partitions we in a few stages to avoid
        a quadratic number of tasks.  This number corresponds to which stage
        we're in, starting from zero up to some small integer
    k: int
        Desired number of splits from this dataframe
    npartition: int
        Total number of output partitions for the full dataframe
    nfinal: int
        Total number of output partitions after repartitioning

    Returns
    -------
    out: Dict[int, DataFrame]
        A dictionary mapping integers in {0..k} to dataframes such that the
        hash values of ``df[col]`` are well partitioned.
    """
    if isinstance(cols, str):
        cols = [cols]

    if cols and cols[0] == "_partitions":
        ind = df[cols[0]]
    else:
        ind = hash_object_dispatch(df[cols] if cols else df, index=False)
        if nfinal and nfinal != npartitions:
            ind = ind % int(nfinal)

    typ = np.min_scalar_type(npartitions * 2)
    # Here we convert the final output index `ind` into the output index
    # for the current stage.
    kwargs = {} if PANDAS_GE_300 else {"copy": False}
    ind = (ind % npartitions).astype(typ, **kwargs) // k**stage % k
    return group_split_dispatch(df, ind, k, ignore_index=ignore_index)


@contextlib.contextmanager
def ensure_cleanup_on_exception(p):
    """Ensure a partd.File is cleaned up.

    We have several tasks referring to a `partd.File` instance. We want to
    ensure that the file is cleaned up if and only if there's an exception
    in the tasks using the `partd.File`.
    """
    try:
        yield
    except Exception:
        # the function (e.g. shuffle_group_3) had an internal exception.
        # We'll cleanup our temporary files and re-raise.
        try:
            p.drop()
        except Exception:
            logger.exception("ignoring exception in ensure_cleanup_on_exception")
        raise


def shuffle_group_3(df, col, npartitions, p):
    with ensure_cleanup_on_exception(p):
        g = df.groupby(col)
        d = {i: g.get_group(i) for i in g.groups}
        p.append(d, fsync=True)


def set_index_post_scalar(df, index_name, drop, column_dtype):
    df2 = df.drop("_partitions", axis=1).set_index(index_name, drop=drop)
    df2.columns = df2.columns.astype(column_dtype)
    return df2


def set_index_post_series(df, index_name, drop, column_dtype):
    df2 = df.drop("_partitions", axis=1).set_index("_index", drop=True)
    df2.index.name = index_name
    df2.columns = df2.columns.astype(column_dtype)
    return df2


def drop_overlap(df, index):
    return df.drop(index) if index in df.index else df


def get_overlap(df, index):
    return df.loc[[index]] if index in df.index else df._constructor()


def fix_overlap(ddf, mins, maxes, lens):
    """Ensures that the upper bound on each partition of ddf (except the last) is exclusive

    This is accomplished by first removing empty partitions, then altering existing
    partitions as needed to include all the values for a particular index value in
    one partition.
    """
    name = "fix-overlap-" + tokenize(ddf, mins, maxes, lens)

    non_empties = [i for i, length in enumerate(lens) if length != 0]
    # If all empty, collapse into one partition
    if len(non_empties) == 0:
        divisions = (None, None)
        dsk = {(name, 0): (ddf._name, 0)}
        graph = HighLevelGraph.from_collections(name, dsk, dependencies=[ddf])
        return new_dd_object(graph, name, ddf._meta, divisions)

    # drop empty partitions by mapping each partition in a new graph to a particular
    # partition on the old graph.
    dsk = {(name, i): (ddf._name, div) for i, div in enumerate(non_empties)}
    ddf_keys = list(dsk.values())
    divisions = tuple(mins) + (maxes[-1],)

    overlap = [i for i in range(1, len(mins)) if mins[i] >= maxes[i - 1]]

    frames = []
    for i in overlap:
        # `frames` is a list of data from previous partitions that we may want to
        # move to partition i.  Here, we add "overlap" from the previous partition
        # (i-1) to this list.
        frames.append((get_overlap, ddf_keys[i - 1], divisions[i]))

        # Make sure that any data added from partition i-1 to `frames` is removed
        # from partition i-1.
        dsk[(name, i - 1)] = (drop_overlap, dsk[(name, i - 1)], divisions[i])

        # We do not want to move "overlap" from the previous partition (i-1) into
        # this partition (i) if the data from this partition will need to be moved
        # to the next partition (i+1) anyway.  If we concatenate data too early,
        # we may lose rows (https://github.com/dask/dask/issues/6972).
        if divisions[i] == divisions[i + 1] and i + 1 in overlap:
            continue

        frames.append(ddf_keys[i])
        dsk[(name, i)] = (methods.concat, frames)
        frames = []

    graph = HighLevelGraph.from_collections(name, dsk, dependencies=[ddf])
    return new_dd_object(graph, name, ddf._meta, divisions)


def _compute_partition_stats(
    column: Series, allow_overlap: bool = False, **kwargs
) -> tuple[list, list, list[int]]:
    """For a given column, compute the min, max, and len of each partition.

    And make sure that the partitions are sorted relative to each other.
    NOTE: this does not guarantee that every partition is internally sorted.
    """
    mins = column.map_partitions(M.min, meta=column)
    maxes = column.map_partitions(M.max, meta=column)
    lens = column.map_partitions(len, meta=column)
    mins, maxes, lens = compute(mins, maxes, lens, **kwargs)
    mins = mins.bfill().tolist()
    maxes = maxes.bfill().tolist()
    non_empty_mins = [m for m, length in zip(mins, lens) if length != 0]
    non_empty_maxes = [m for m, length in zip(maxes, lens) if length != 0]
    if (
        sorted(non_empty_mins) != non_empty_mins
        or sorted(non_empty_maxes) != non_empty_maxes
    ):
        raise ValueError(
            f"Partitions are not sorted ascending by {column.name or 'the index'}",
            f"In your dataset the (min, max, len) values of {column.name or 'the index'} "
            f"for each partition are : {list(zip(mins, maxes, lens))}",
        )
    if not allow_overlap and any(
        a <= b for a, b in zip(non_empty_mins[1:], non_empty_maxes[:-1])
    ):
        warnings.warn(
            "Partitions have overlapping values, so divisions are non-unique."
            "Use `set_index(sorted=True)` with no `divisions` to allow dask to fix the overlap. "
            f"In your dataset the (min, max, len) values of {column.name or 'the index'} "
            f"for each partition are : {list(zip(mins, maxes, lens))}",
            UserWarning,
        )
    lens = methods.tolist(lens)
    if not allow_overlap:
        return (mins, maxes, lens)
    else:
        return (non_empty_mins, non_empty_maxes, lens)


def compute_divisions(df: DataFrame, col: Any | None = None, **kwargs) -> tuple:
    column = df.index if col is None else df[col]
    mins, maxes, _ = _compute_partition_stats(column, allow_overlap=False, **kwargs)

    return tuple(mins) + (maxes[-1],)


def compute_and_set_divisions(df: DataFrame, **kwargs) -> DataFrame:
    mins, maxes, lens = _compute_partition_stats(df.index, allow_overlap=True, **kwargs)
    if len(mins) == len(df.divisions) - 1:
        df._divisions = tuple(mins) + (maxes[-1],)
        if not any(mins[i] >= maxes[i - 1] for i in range(1, len(mins))):
            return df

    return fix_overlap(df, mins, maxes, lens)


def set_sorted_index(
    df: DataFrame,
    index: str | Series,
    drop: bool = True,
    divisions: Sequence | None = None,
    **kwargs,
) -> DataFrame:
    if isinstance(index, Series):
        meta = df._meta.set_index(index._meta, drop=drop)
    else:
        meta = df._meta.set_index(index, drop=drop)

    result = map_partitions(
        M.set_index,
        df,
        index,
        drop=drop,
        meta=meta,
        align_dataframes=False,
        transform_divisions=False,
    )

    if not divisions:
        return compute_and_set_divisions(result, **kwargs)
    elif len(divisions) != len(df.divisions):
        msg = (
            "When doing `df.set_index(col, sorted=True, divisions=...)`, "
            "divisions indicates known splits in the index column. In this "
            "case divisions must be the same length as the existing "
            "divisions in `df`\n\n"
            "If the intent is to repartition into new divisions after "
            "setting the index, you probably want:\n\n"
            "`df.set_index(col, sorted=True).repartition(divisions=divisions)`"
        )
        raise ValueError(msg)

    result.divisions = tuple(divisions)
    return result
