from __future__ import annotations

import datetime
import warnings
from numbers import Integral

import pandas as pd
from pandas.api.types import is_datetime64_any_dtype
from pandas.core.window import Rolling as pd_Rolling

from dask.array.core import normalize_arg
from dask.base import tokenize
from dask.blockwise import BlockwiseDepDict
from dask.dataframe import methods
from dask.dataframe._compat import check_axis_keyword_deprecation
from dask.dataframe.core import (
    Scalar,
    _Frame,
    _get_divisions_map_partitions,
    _get_meta_map_partitions,
    _maybe_from_pandas,
    apply_and_enforce,
    new_dd_object,
    partitionwise_graph,
)
from dask.dataframe.io import from_pandas
from dask.dataframe.multi import _maybe_align_partitions
from dask.dataframe.utils import (
    insert_meta_param_description,
    is_dask_collection,
    is_dataframe_like,
    is_series_like,
)
from dask.delayed import unpack_collections
from dask.highlevelgraph import HighLevelGraph
from dask.typing import no_default
from dask.utils import M, apply, derived_from, funcname, has_keyword

CombinedOutput = type("CombinedOutput", (tuple,), {})


def _combined_parts(prev_part, current_part, next_part, before, after):
    msg = (
        "Partition size is less than overlapping "
        "window size. Try using ``df.repartition`` "
        "to increase the partition size."
    )

    if prev_part is not None and isinstance(before, Integral):
        if prev_part.shape[0] != before:
            raise NotImplementedError(msg)

    if next_part is not None and isinstance(after, Integral):
        if next_part.shape[0] != after:
            raise NotImplementedError(msg)

    parts = [p for p in (prev_part, current_part, next_part) if p is not None]
    combined = methods.concat(parts)

    return CombinedOutput(
        (
            combined,
            len(prev_part) if prev_part is not None else None,
            len(next_part) if next_part is not None else None,
        )
    )


def overlap_chunk(func, before, after, *args, **kwargs):
    dfs = [df for df in args if isinstance(df, CombinedOutput)]
    combined, prev_part_length, next_part_length = dfs[0]

    args = [arg[0] if isinstance(arg, CombinedOutput) else arg for arg in args]

    out = func(*args, **kwargs)

    if prev_part_length is None:
        before = None
    if isinstance(before, datetime.timedelta):
        before = prev_part_length

    expansion = None
    if combined.shape[0] != 0:
        expansion = out.shape[0] // combined.shape[0]
    if before and expansion:
        before *= expansion
    if next_part_length is None:
        return out.iloc[before:]
    if isinstance(after, datetime.timedelta):
        after = next_part_length
    if after and expansion:
        after *= expansion
    return out.iloc[before:-after]


@insert_meta_param_description
def map_overlap(
    func,
    df,
    before,
    after,
    *args,
    meta=no_default,
    enforce_metadata=True,
    transform_divisions=True,
    align_dataframes=True,
    **kwargs,
):
    """Apply a function to each partition, sharing rows with adjacent partitions.

    Parameters
    ----------
    func : function
        The function applied to each partition. If this function accepts
        the special ``partition_info`` keyword argument, it will receive
        information on the partition's relative location within the
        dataframe.
    df: dd.DataFrame, dd.Series
    args, kwargs :
        Positional and keyword arguments to pass to the function.
        Positional arguments are computed on a per-partition basis, while
        keyword arguments are shared across all partitions. The partition
        itself will be the first positional argument, with all other
        arguments passed *after*. Arguments can be ``Scalar``, ``Delayed``,
        or regular Python objects. DataFrame-like args (both dask and
        pandas) will be repartitioned to align (if necessary) before
        applying the function; see ``align_dataframes`` to control this
        behavior.
    enforce_metadata : bool, default True
        Whether to enforce at runtime that the structure of the DataFrame
        produced by ``func`` actually matches the structure of ``meta``.
        This will rename and reorder columns for each partition,
        and will raise an error if this doesn't work,
        but it won't raise if dtypes don't match.
    before : int, timedelta or string timedelta
        The rows to prepend to partition ``i`` from the end of
        partition ``i - 1``.
    after : int, timedelta or string timedelta
        The rows to append to partition ``i`` from the beginning
        of partition ``i + 1``.
    transform_divisions : bool, default True
        Whether to apply the function onto the divisions and apply those
        transformed divisions to the output.
    align_dataframes : bool, default True
        Whether to repartition DataFrame- or Series-like args
        (both dask and pandas) so their divisions align before applying
        the function. This requires all inputs to have known divisions.
        Single-partition inputs will be split into multiple partitions.

        If False, all inputs must have either the same number of partitions
        or a single partition. Single-partition inputs will be broadcast to
        every partition of multi-partition inputs.
    $META

    See Also
    --------
    dd.DataFrame.map_overlap
    """

    df = (
        from_pandas(df, 1)
        if (is_series_like(df) or is_dataframe_like(df)) and not is_dask_collection(df)
        else df
    )

    args = (df,) + args

    if isinstance(before, str):
        before = pd.to_timedelta(before)
    if isinstance(after, str):
        after = pd.to_timedelta(after)

    if isinstance(before, datetime.timedelta) or isinstance(after, datetime.timedelta):
        if not is_datetime64_any_dtype(df.index._meta_nonempty.inferred_type):
            raise TypeError(
                "Must have a `DatetimeIndex` when using string offset "
                "for `before` and `after`"
            )
    else:
        if not (
            isinstance(before, Integral)
            and before >= 0
            and isinstance(after, Integral)
            and after >= 0
        ):
            raise ValueError("before and after must be positive integers")

    name = kwargs.pop("token", None)
    parent_meta = kwargs.pop("parent_meta", None)

    assert callable(func)
    if name is not None:
        token = tokenize(meta, before, after, *args, **kwargs)
    else:
        name = "overlap-" + funcname(func)
        token = tokenize(func, meta, before, after, *args, **kwargs)
    name = f"{name}-{token}"

    if align_dataframes:
        args = _maybe_from_pandas(args)
        try:
            args = _maybe_align_partitions(args)
        except ValueError as e:
            raise ValueError(
                f"{e}. If you don't want the partitions to be aligned, and are "
                "calling `map_overlap` directly, pass `align_dataframes=False`."
            ) from e

    dfs = [df for df in args if isinstance(df, _Frame)]

    meta = _get_meta_map_partitions(args, dfs, func, kwargs, meta, parent_meta)

    if all(isinstance(arg, Scalar) for arg in args):
        layer = {
            (name, 0): (
                apply,
                func,
                (tuple, [(arg._name, 0) for arg in args]),
                kwargs,
            )
        }
        graph = HighLevelGraph.from_collections(name, layer, dependencies=args)
        return Scalar(graph, name, meta)

    args2 = []
    dependencies = []

    divisions = _get_divisions_map_partitions(
        align_dataframes, transform_divisions, dfs, func, args, kwargs
    )

    def _handle_frame_argument(arg):
        dsk = {}

        prevs_parts_dsk, prevs = _get_previous_partitions(arg, before)
        dsk.update(prevs_parts_dsk)

        nexts_parts_dsk, nexts = _get_nexts_partitions(arg, after)
        dsk.update(nexts_parts_dsk)

        name_a = "overlap-concat-" + tokenize(arg)
        for i, (prev, current, next) in enumerate(
            zip(prevs, arg.__dask_keys__(), nexts)
        ):
            key = (name_a, i)
            dsk[key] = (_combined_parts, prev, current, next, before, after)

        graph = HighLevelGraph.from_collections(name_a, dsk, dependencies=[arg])
        return new_dd_object(graph, name_a, meta, divisions)

    for arg in args:
        if isinstance(arg, _Frame):
            arg = _handle_frame_argument(arg)
            args2.append(arg)
            dependencies.append(arg)
            continue
        arg = normalize_arg(arg)
        arg2, collections = unpack_collections(arg)
        if collections:
            args2.append(arg2)
            dependencies.extend(collections)
        else:
            args2.append(arg)

    kwargs3 = {}
    simple = True
    for k, v in kwargs.items():
        v = normalize_arg(v)
        v, collections = unpack_collections(v)
        dependencies.extend(collections)
        kwargs3[k] = v
        if collections:
            simple = False

    if has_keyword(func, "partition_info"):
        partition_info = {
            (i,): {"number": i, "division": division}
            for i, division in enumerate(divisions[:-1])
        }

        args2.insert(0, BlockwiseDepDict(partition_info))
        orig_func = func

        def func(partition_info, *args, **kwargs):
            return orig_func(*args, **kwargs, partition_info=partition_info)

    if enforce_metadata:
        dsk = partitionwise_graph(
            apply_and_enforce,
            name,
            func,
            before,
            after,
            *args2,
            dependencies=dependencies,
            _func=overlap_chunk,
            _meta=meta,
            **kwargs3,
        )
    else:
        kwargs4 = kwargs if simple else kwargs3
        dsk = partitionwise_graph(
            overlap_chunk,
            name,
            func,
            before,
            after,
            *args2,
            **kwargs4,
            dependencies=dependencies,
        )

    graph = HighLevelGraph.from_collections(name, dsk, dependencies=dependencies)
    return new_dd_object(graph, name, meta, divisions)


def _get_nexts_partitions(df, after):
    """
    Helper to get the nexts partitions required for the overlap
    """
    dsk = {}
    df_name = df._name

    timedelta_partition_message = (
        "Partition size is less than specified window. "
        "Try using ``df.repartition`` to increase the partition size"
    )

    name_b = "overlap-append-" + tokenize(df, after)
    if after and isinstance(after, Integral):
        nexts = []
        for i in range(1, df.npartitions):
            key = (name_b, i)
            dsk[key] = (M.head, (df_name, i), after)
            nexts.append(key)
        nexts.append(None)
    elif isinstance(after, datetime.timedelta):
        # TODO: Do we have a use-case for this? Pandas doesn't allow negative rolling windows
        deltas = pd.Series(df.divisions).diff().iloc[1:-1]
        if (after > deltas).any():
            raise ValueError(timedelta_partition_message)

        nexts = []
        for i in range(1, df.npartitions):
            key = (name_b, i)
            dsk[key] = (_head_timedelta, (df_name, i - 0), (df_name, i), after)
            nexts.append(key)
        nexts.append(None)
    else:
        nexts = [None] * df.npartitions
    return dsk, nexts


def _get_previous_partitions(df, before):
    """
    Helper to get the previous partitions required for the overlap
    """
    dsk = {}
    df_name = df._name

    name_a = "overlap-prepend-" + tokenize(df, before)
    if before and isinstance(before, Integral):
        prevs = [None]
        for i in range(df.npartitions - 1):
            key = (name_a, i)
            dsk[key] = (M.tail, (df_name, i), before)
            prevs.append(key)

    elif isinstance(before, datetime.timedelta):
        # Assumes monotonic (increasing?) index
        divs = pd.Series(df.divisions)
        deltas = divs.diff().iloc[1:-1]

        # In the first case window-size is larger than at least one partition, thus it is
        # necessary to calculate how many partitions must be used for each rolling task.
        # Otherwise, these calculations can be skipped (faster)

        if (before > deltas).any():
            pt_z = divs[0]
            prevs = [None]
            for i in range(df.npartitions - 1):
                # Select all indexes of relevant partitions between the current partition and
                # the partition with the highest division outside the rolling window (before)
                pt_i = divs[i + 1]

                # lower-bound the search to the first division
                lb = max(pt_i - before, pt_z)

                first, j = divs[i], i
                while first > lb and j > 0:
                    first = first - deltas[j]
                    j = j - 1

                key = (name_a, i)
                dsk[key] = (
                    _tail_timedelta,
                    [(df_name, k) for k in range(j, i + 1)],
                    (df_name, i + 1),
                    before,
                )
                prevs.append(key)

        else:
            prevs = [None]
            for i in range(df.npartitions - 1):
                key = (name_a, i)
                dsk[key] = (
                    _tail_timedelta,
                    [(df_name, i)],
                    (df_name, i + 1),
                    before,
                )
                prevs.append(key)
    else:
        prevs = [None] * df.npartitions
    return dsk, prevs


def _head_timedelta(current, next_, after):
    """Return rows of ``next_`` whose index is before the last
    observation in ``current`` + ``after``.

    Parameters
    ----------
    current : DataFrame
    next_ : DataFrame
    after : timedelta

    Returns
    -------
    overlapped : DataFrame
    """
    return next_[next_.index < (current.index.max() + after)]


def _tail_timedelta(prevs, current, before):
    """Return the concatenated rows of each dataframe in ``prevs`` whose
    index is after the first observation in ``current`` - ``before``.

    Parameters
    ----------
    current : DataFrame
    prevs : list of DataFrame objects
    before : timedelta

    Returns
    -------
    overlapped : DataFrame
    """
    selected = methods.concat(
        [prev[prev.index > (current.index.min() - before)] for prev in prevs]
    )
    return selected


class Rolling:
    """Provides rolling window calculations."""

    def __init__(
        self,
        obj,
        window=None,
        min_periods=None,
        center=False,
        win_type=None,
        axis=no_default,
    ):
        self.obj = obj  # dataframe or series
        self.window = window
        self.min_periods = min_periods
        self.center = center
        self.axis = axis
        self.win_type = win_type
        # Allow pandas to raise if appropriate
        obj._meta.rolling(**self._rolling_kwargs())
        # Using .rolling(window='2s'), pandas will convert the
        # offset str to a window in nanoseconds. But pandas doesn't
        # accept the integer window with win_type='freq', so we store
        # that information here.
        # See https://github.com/pandas-dev/pandas/issues/15969
        self._win_type = None if isinstance(self.window, int) else "freq"

        if self.axis in ("index", 1, "rows"):
            warnings.warn(
                "Using axis=1 in Rolling has been deprecated and will be removed "
                "in a future version.",
                FutureWarning,
            )

    def _rolling_kwargs(self):
        kwargs = {
            "window": self.window,
            "min_periods": self.min_periods,
            "center": self.center,
            "win_type": self.win_type,
        }
        if self.axis is not no_default:
            kwargs["axis"] = self.axis
        return kwargs

    @property
    def _has_single_partition(self):
        """
        Indicator for whether the object has a single partition (True)
        or multiple (False).
        """
        return (
            self.axis in (1, "columns")
            or (isinstance(self.window, Integral) and self.window <= 1)
            or self.obj.npartitions == 1
        )

    @staticmethod
    def pandas_rolling_method(df, rolling_kwargs, name, *args, **kwargs):
        with check_axis_keyword_deprecation():
            rolling = df.rolling(**rolling_kwargs)
        return getattr(rolling, name)(*args, **kwargs)

    def _call_method(self, method_name, *args, **kwargs):
        rolling_kwargs = self._rolling_kwargs()
        meta = self.pandas_rolling_method(
            self.obj._meta_nonempty, rolling_kwargs, method_name, *args, **kwargs
        )

        if self._has_single_partition:
            # There's no overlap just use map_partitions
            return self.obj.map_partitions(
                self.pandas_rolling_method,
                rolling_kwargs,
                method_name,
                *args,
                token=method_name,
                meta=meta,
                **kwargs,
            )
        # Convert window to overlap
        if self.center:
            before = self.window // 2
            after = self.window - before - 1
        elif self._win_type == "freq":
            before = pd.Timedelta(self.window)
            after = 0
        else:
            before = self.window - 1
            after = 0
        return map_overlap(
            self.pandas_rolling_method,
            self.obj,
            before,
            after,
            rolling_kwargs,
            method_name,
            *args,
            token=method_name,
            meta=meta,
            **kwargs,
        )

    @derived_from(pd_Rolling)
    def count(self):
        return self._call_method("count")

    @derived_from(pd_Rolling)
    def cov(self):
        return self._call_method("cov")

    @derived_from(pd_Rolling)
    def sum(self):
        return self._call_method("sum")

    @derived_from(pd_Rolling)
    def mean(self):
        return self._call_method("mean")

    @derived_from(pd_Rolling)
    def median(self):
        return self._call_method("median")

    @derived_from(pd_Rolling)
    def min(self):
        return self._call_method("min")

    @derived_from(pd_Rolling)
    def max(self):
        return self._call_method("max")

    @derived_from(pd_Rolling)
    def std(self, ddof=1):
        return self._call_method("std", ddof=1)

    @derived_from(pd_Rolling)
    def var(self, ddof=1):
        return self._call_method("var", ddof=1)

    @derived_from(pd_Rolling)
    def skew(self):
        return self._call_method("skew")

    @derived_from(pd_Rolling)
    def kurt(self):
        return self._call_method("kurt")

    @derived_from(pd_Rolling)
    def quantile(self, quantile):
        return self._call_method("quantile", quantile)

    @derived_from(pd_Rolling)
    def apply(
        self,
        func,
        raw=False,
        engine="cython",
        engine_kwargs=None,
        args=None,
        kwargs=None,
    ):
        kwargs = kwargs or {}
        args = args or ()
        return self._call_method(
            "apply",
            func,
            raw=raw,
            engine=engine,
            engine_kwargs=engine_kwargs,
            args=args,
            kwargs=kwargs,
        )

    @derived_from(pd_Rolling)
    def aggregate(self, func, *args, **kwargs):
        return self._call_method("agg", func, *args, **kwargs)

    agg = aggregate

    def __repr__(self):
        def order(item):
            k, v = item
            _order = {
                "window": 0,
                "min_periods": 1,
                "center": 2,
                "win_type": 3,
                "axis": 4,
            }
            return _order[k]

        rolling_kwargs = self._rolling_kwargs()
        rolling_kwargs["window"] = self.window
        rolling_kwargs["win_type"] = self._win_type
        return "Rolling [{}]".format(
            ",".join(
                f"{k}={v}"
                for k, v in sorted(rolling_kwargs.items(), key=order)
                if v is not None
            )
        )


class RollingGroupby(Rolling):
    def __init__(
        self,
        groupby,
        window=None,
        min_periods=None,
        center=False,
        win_type=None,
        axis=0,
    ):
        self._groupby_kwargs = groupby._groupby_kwargs
        self._groupby_slice = groupby._slice

        obj = groupby.obj
        if self._groupby_slice is not None:
            if isinstance(self._groupby_slice, str):
                sliced_plus = [self._groupby_slice]
            else:
                sliced_plus = list(self._groupby_slice)
            if isinstance(groupby.by, str):
                sliced_plus.append(groupby.by)
            else:
                sliced_plus.extend(groupby.by)
            obj = obj[sliced_plus]

        super().__init__(
            obj,
            window=window,
            min_periods=min_periods,
            center=center,
            win_type=win_type,
            axis=axis,
        )

    def _rolling_kwargs(self):
        kwargs = super()._rolling_kwargs()
        if kwargs.get("axis", None) in (0, "index"):
            # it's a default, no need to pass
            kwargs.pop("axis")
        return kwargs

    @staticmethod
    def pandas_rolling_method(
        df,
        rolling_kwargs,
        name,
        *args,
        groupby_kwargs=None,
        groupby_slice=None,
        **kwargs,
    ):
        groupby = df.groupby(**groupby_kwargs)
        if groupby_slice:
            groupby = groupby[groupby_slice]
        rolling = groupby.rolling(**rolling_kwargs)
        return getattr(rolling, name)(*args, **kwargs).sort_index(level=-1)

    def _call_method(self, method_name, *args, **kwargs):
        return super()._call_method(
            method_name,
            *args,
            groupby_kwargs=self._groupby_kwargs,
            groupby_slice=self._groupby_slice,
            **kwargs,
        )
