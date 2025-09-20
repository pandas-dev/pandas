from __future__ import annotations

import functools
import math
import warnings
from collections.abc import Callable

import numpy as np
import pandas as pd
from pandas.core.apply import reconstruct_func, validate_func_kwargs

from dask import is_dask_collection
from dask._task_spec import Task
from dask.core import flatten
from dask.dataframe.core import (
    _concat,
    apply_and_enforce,
    is_dataframe_like,
    is_series_like,
)
from dask.dataframe.dask_expr._collection import (
    FrameBase,
    Index,
    Series,
    new_collection,
)
from dask.dataframe.dask_expr._expr import (
    Assign,
    Blockwise,
    Expr,
    MapPartitions,
    Projection,
    RenameFrame,
    RenameSeries,
    ToFrame,
    _DeepCopy,
    _extract_meta,
    are_co_aligned,
    determine_column_projection,
    no_default,
)
from dask.dataframe.dask_expr._reductions import ApplyConcatApply, Chunk, Reduction
from dask.dataframe.dask_expr._shuffle import RearrangeByColumn
from dask.dataframe.dask_expr._util import (
    PANDAS_GE_300,
    _convert_to_list,
    get_specified_shuffle,
)
from dask.dataframe.dispatch import concat, make_meta, meta_nonempty
from dask.dataframe.groupby import (
    GROUP_KEYS_DEFAULT,
    _agg_finalize,
    _aggregate_docstring,
    _apply_chunk,
    _build_agg_args,
    _cov_agg,
    _cov_chunk,
    _cum_agg_aligned,
    _cum_agg_filled,
    _cumcount_aggregate,
    _determine_levels,
    _groupby_aggregate,
    _groupby_aggregate_spec,
    _groupby_apply_funcs,
    _groupby_get_group,
    _groupby_raise_unaligned,
    _groupby_slice_apply,
    _groupby_slice_shift,
    _groupby_slice_transform,
    _head_aggregate,
    _head_chunk,
    _non_agg_chunk,
    _normalize_spec,
    _nunique_df_chunk,
    _nunique_df_combine,
    _tail_aggregate,
    _tail_chunk,
    _unique_aggregate,
    _value_counts,
    _value_counts_aggregate,
    _var_agg,
    _var_chunk,
)
from dask.dataframe.utils import insert_meta_param_description, is_scalar
from dask.typing import Key
from dask.utils import M, derived_from, is_index_like


def _as_dict(key, value):
    # Utility to convert a single kwarg to a dict.
    # The dict will be empty if the value is None
    return {} if value is None else {key: value}


def _adjust_split_out_for_group_keys(npartitions, by):
    if len(by) == 1:
        return math.ceil(npartitions / 15)
    return math.ceil(npartitions / (10 / (len(by) - 1)))


class Aggregation:
    """User defined groupby-aggregation.

    This class allows users to define their own custom aggregation in terms of
    operations on Pandas dataframes in a map-reduce style. You need to specify
    what operation to do on each chunk of data, how to combine those chunks of
    data together, and then how to finalize the result.

    See :ref:`dataframe.groupby.aggregate` for more.

    Parameters
    ----------
    name : str
        the name of the aggregation. It should be unique, since intermediate
        result will be identified by this name.
    chunk : callable
        a function that will be called with the grouped column of each
        partition, takes a Pandas SeriesGroupBy in input.
        It can either return a single series or a tuple of series.
        The index has to be equal to the groups.
    agg : callable
        a function that will be called to aggregate the results of each chunk.
        Again the argument(s) will be a Pandas SeriesGroupBy. If ``chunk``
        returned a tuple, ``agg`` will be called with all of them as
        individual positional arguments.
    finalize : callable
        an optional finalizer that will be called with the results from the
        aggregation.

    Examples
    --------
    We could implement ``sum`` as follows:

    >>> custom_sum = dd.Aggregation(
    ...     name='custom_sum',
    ...     chunk=lambda s: s.sum(),
    ...     agg=lambda s0: s0.sum()
    ... )  # doctest: +SKIP
    >>> df.groupby('g').agg(custom_sum)  # doctest: +SKIP

    We can implement ``mean`` as follows:

    >>> custom_mean = dd.Aggregation(
    ...     name='custom_mean',
    ...     chunk=lambda s: (s.count(), s.sum()),
    ...     agg=lambda count, sum: (count.sum(), sum.sum()),
    ...     finalize=lambda count, sum: sum / count,
    ... )  # doctest: +SKIP
    >>> df.groupby('g').agg(custom_mean)  # doctest: +SKIP

    Though of course, both of these are built-in and so you don't need to
    implement them yourself.
    """

    def __init__(self, name, chunk, agg, finalize=None):
        self.chunk = chunk
        self.agg = agg
        self.finalize = finalize
        self.__name__ = name


###
### Groupby-aggregation expressions
###


class GroupByBase:
    @functools.cached_property
    def _by_meta(self):
        return [meta_nonempty(x._meta) if isinstance(x, Expr) else x for x in self.by]

    @functools.cached_property
    def _by_columns(self):
        return [x for x in self.by if not isinstance(x, Expr)]

    @property
    def split_by(self):
        return list(
            flatten(
                [[x] if not isinstance(x, Expr) else x.columns for x in self.by],
                container=list,
            )
        )

    @functools.cached_property
    def by(self):
        return self.operands[len(self._parameters) :]

    @functools.cached_property
    def levels(self):
        return _determine_levels(self.by)

    @property
    def shuffle_by_index(self):
        return True


class GroupByChunk(Chunk, GroupByBase):
    _preserves_partitioning_information = True

    @functools.cached_property
    def _args(self) -> list:
        return [self.frame] + self.by

    @functools.cached_property
    def _meta(self):
        args = [
            meta_nonempty(op._meta) if isinstance(op, Expr) else op for op in self._args
        ]
        return make_meta(self.operation(*args, **self._kwargs))


class GroupByApplyConcatApply(ApplyConcatApply, GroupByBase):
    _chunk_cls = GroupByChunk

    @functools.cached_property
    def _meta_chunk(self):
        meta = meta_nonempty(self.frame._meta)
        return self.chunk(meta, *self._by_meta, **self.chunk_kwargs)

    @property
    def _chunk_cls_args(self):
        return self.by

    @property
    def split_out(self):
        if self.operand("split_out") is None:
            return 1
        return super().split_out

    @property
    def _projection_columns(self):
        return self.frame.columns

    def _tune_down(self):
        if self.operand("split_out") is None:
            return self.substitute_parameters(
                {
                    "split_out": functools.partial(
                        _adjust_split_out_for_group_keys, by=self.by
                    )
                }
            )


class SingleAggregation(GroupByApplyConcatApply, GroupByBase):
    """Single groupby aggregation

    This is an abstract class. Sub-classes must implement
    the following methods:

    -   `groupby_chunk`: Applied to each group within
        the `chunk` method of `GroupByApplyConcatApply`
    -   `groupby_aggregate`: Applied to each group within
        the `aggregate` method of `GroupByApplyConcatApply`

    Parameters
    ----------
    frame: Expr
        Dataframe- or series-like expression to group.
    by: str, list or Series
        The key for grouping
    observed:
        Passed through to dataframe backend.
    dropna:
        Whether rows with NA values should be dropped.
    chunk_kwargs:
        Key-word arguments to pass to `groupby_chunk`.
    aggregate_kwargs:
        Key-word arguments to pass to `aggregate_chunk`.
    """

    _parameters = [
        "frame",
        "observed",
        "dropna",
        "chunk_kwargs",
        "aggregate_kwargs",
        "_slice",
        "split_every",
        "split_out",
        "sort",
        "shuffle_method",
    ]
    _defaults = {
        "observed": None,
        "dropna": None,
        "chunk_kwargs": None,
        "aggregate_kwargs": None,
        "_slice": None,
        "split_every": 8,
        "split_out": None,
        "sort": None,
        "shuffle_method": None,
    }

    groupby_chunk: Callable | None = None
    groupby_aggregate: Callable | None = None

    @classmethod
    def chunk(cls, df, *by, **kwargs):
        return _apply_chunk(df, *by, **kwargs)

    @classmethod
    def aggregate(cls, inputs, **kwargs):
        return _groupby_aggregate(_concat(inputs), **kwargs)

    @property
    def chunk_kwargs(self) -> dict:  # type: ignore
        chunk_kwargs = self.operand("chunk_kwargs") or {}
        columns = self._slice
        return {
            "chunk": self.groupby_chunk,
            "columns": columns,
            **_as_dict("observed", self.observed),
            **_as_dict("dropna", self.dropna),
            **chunk_kwargs,
        }

    @property
    def aggregate_kwargs(self) -> dict:  # type: ignore
        aggregate_kwargs = self.operand("aggregate_kwargs") or {}
        groupby_aggregate = self.groupby_aggregate or self.groupby_chunk
        return {
            "aggfunc": groupby_aggregate,
            "levels": self.levels,
            "sort": self.sort,
            **_as_dict("observed", self.observed),
            **_as_dict("dropna", self.dropna),
            **aggregate_kwargs,
        }

    def _simplify_up(self, parent, dependents):
        return groupby_projection(self, parent, dependents)


class GroupbyAggregationBase(GroupByApplyConcatApply, GroupByBase):
    """Base class for groupby aggregation

    This class can be subclassed to perform a general
    groupby aggregation by passing in a `str`, `list` or
    `dict`-based specification using the `arg` operand.

    Parameters
    ----------
    frame: Expr
        Dataframe- or series-like expression to group.
    by: str, list or Series
        The key for grouping
    arg: str, list or dict
        Aggregation spec defining the specific aggregations
        to perform.
    observed:
        Passed through to dataframe backend.
    dropna:
        Whether rows with NA values should be dropped.
    """

    _parameters = [
        "frame",
        "arg",
        "observed",
        "dropna",
        "split_every",
        "split_out",
        "sort",
        "shuffle_method",
        "_slice",
    ]
    _defaults = {
        "observed": None,
        "dropna": None,
        "split_every": 8,
        "split_out": None,
        "sort": None,
        "shuffle_method": None,
        "_slice": None,
    }

    @functools.cached_property
    def spec(self):
        # Converts the `arg` operand into specific
        # chunk, aggregate, and finalizer functions
        if is_dataframe_like(self.frame._meta):
            group_columns = self._by_columns
            if self._slice:
                non_group_columns = self._slice
                if is_scalar(non_group_columns):
                    non_group_columns = [non_group_columns]
            else:
                non_group_columns = [
                    col for col in self.frame.columns if col not in group_columns
                ]
            spec = _normalize_spec(self.arg, non_group_columns)
        elif is_series_like(self.frame._meta):
            if isinstance(self.arg, (list, tuple, dict)):
                spec = _normalize_spec({None: self.arg}, [])
                spec = [
                    (result_column, func, input_column)
                    for ((_, result_column), func, input_column) in spec
                ]

            else:
                spec = _normalize_spec({None: self.arg}, [])
                spec = [
                    (self.frame.columns[0], func, input_column)
                    for (_, func, input_column) in spec
                ]
        else:
            raise ValueError(f"aggregate on unknown object {self.frame._meta}")

        return spec

    @functools.cached_property
    def agg_args(self):
        keys = ["chunk_funcs", "aggregate_funcs", "finalizers"]
        return dict(zip(keys, _build_agg_args(self.spec)))

    def _simplify_down(self):
        if not isinstance(self.arg, dict):
            return

        # Use agg-spec information to add column projection
        required_columns = (
            set(self._by_columns)
            .union(self.arg.keys())
            .intersection(self.frame.columns)
        )
        column_projection = [
            column for column in self.frame.columns if column in required_columns
        ]
        if column_projection != self.frame.columns:
            return type(self)(self.frame[column_projection], *self.operands[1:])


class GroupbyAggregation(GroupbyAggregationBase):
    """Logical groupby aggregation class

    This class lowers itself to concrete implementations for decomposable
    or holistic aggregations.
    """

    @functools.cached_property
    def _meta(self):
        return self._lower()._meta

    @functools.cached_property
    def _is_decomposable(self):
        return not any(s[1] in ("median", np.median) for s in self.spec)

    def _lower(self):
        cls = (
            DecomposableGroupbyAggregation
            if self._is_decomposable
            else HolisticGroupbyAggregation
        )
        return cls(
            self.frame,
            self.arg,
            self.observed,
            self.dropna,
            self.split_every,
            self.split_out,
            self.sort,
            self.shuffle_method,
            self._slice,
            *self.by,
        )


class HolisticGroupbyAggregation(GroupbyAggregationBase):
    """Groupby aggregation for both decomposable and non-decomposable aggregates

    This class always calculates the aggregates by first collecting all the data for
    the groups and then aggregating at once.

    We are always shuffling, so we will never call combine
    """

    @functools.cached_property
    def _meta(self):
        meta = self._meta_chunk
        aggregate = self.aggregate or (lambda x: x)
        aggregate_kwargs = self.aggregate_kwargs
        meta = aggregate([meta], **aggregate_kwargs)
        return make_meta(meta)

    chunk = staticmethod(_non_agg_chunk)

    @property
    def should_shuffle(self):
        return True

    @classmethod
    def aggregate(cls, inputs, **kwargs):
        return _groupby_aggregate_spec(_concat(inputs), **kwargs)

    @property
    def chunk_kwargs(self) -> dict:  # type: ignore
        return {
            "by": self._by_columns,
            "key": [col for col in self.frame.columns if col not in self._by_columns],
            **_as_dict("observed", self.observed),
            **_as_dict("dropna", self.dropna),
        }

    @property
    def aggregate_kwargs(self) -> dict:  # type: ignore
        return {
            "spec": self.arg,
            "levels": _determine_levels(self.by),
            **_as_dict("observed", self.observed),
            **_as_dict("dropna", self.dropna),
        }


class DecomposableGroupbyAggregation(GroupbyAggregationBase):
    """Groupby aggregation for decomposable aggregates

    The results may be calculated via tree or shuffle reduction.
    """

    chunk = staticmethod(_groupby_apply_funcs)

    @classmethod
    def combine(cls, inputs, **kwargs):
        return _groupby_apply_funcs(_concat(inputs), **kwargs)

    @classmethod
    def aggregate(cls, inputs, **kwargs):
        return _agg_finalize(_concat(inputs), **kwargs)

    @property
    def chunk_kwargs(self) -> dict:  # type: ignore
        return {
            "funcs": self.agg_args["chunk_funcs"],
            "sort": self.sort,
            **_as_dict("observed", self.observed),
            **_as_dict("dropna", self.dropna),
        }

    @property
    def combine_kwargs(self) -> dict:  # type: ignore
        return {
            "funcs": self.agg_args["aggregate_funcs"],
            "level": self.levels,
            "sort": self.sort,
            **_as_dict("observed", self.observed),
            **_as_dict("dropna", self.dropna),
        }

    @property
    def aggregate_kwargs(self) -> dict:  # type: ignore
        return {
            "aggregate_funcs": self.agg_args["aggregate_funcs"],
            "arg": self.arg,
            "columns": self._slice,
            "finalize_funcs": self.agg_args["finalizers"],
            "is_series": self.frame._meta.ndim == 1,
            "level": self.levels,
            "sort": self.sort,
            **_as_dict("observed", self.observed),
            **_as_dict("dropna", self.dropna),
        }


class Sum(SingleAggregation):
    groupby_chunk = M.sum


class Prod(SingleAggregation):
    groupby_chunk = M.prod


class Min(SingleAggregation):
    groupby_chunk = M.min


class Max(SingleAggregation):
    groupby_chunk = M.max


class First(SingleAggregation):
    groupby_chunk = M.first


class Last(SingleAggregation):
    groupby_chunk = M.last


class Count(SingleAggregation):
    groupby_chunk = M.count
    groupby_aggregate = M.sum


class Size(SingleAggregation):
    groupby_chunk = M.size
    groupby_aggregate = M.sum

    def _simplify_down(self):
        if (
            self._slice is not None
            and not isinstance(self._slice, list)
            or self.frame.ndim == 1
        ):
            # Scalar slices influence the result and are allowed, i.e., the name of
            # the series is different
            return

        # We can remove every column since pandas reduces to a Series anyway
        by_columns = self._by_columns
        by_columns = [c for c in by_columns if c in self.frame.columns]
        if set(by_columns) == set(self.frame.columns):
            return

        slice_idx = self._parameters.index("_slice")
        ops = [op if i != slice_idx else None for i, op in enumerate(self.operands)]
        return type(self)(self.frame[by_columns], *ops[1:])


class IdxMin(SingleAggregation):
    groupby_chunk = M.idxmin
    groupby_aggregate = M.first


class IdxMax(IdxMin):
    groupby_chunk = M.idxmax
    groupby_aggregate = M.first


class ValueCounts(SingleAggregation):
    groupby_chunk = staticmethod(_value_counts)
    groupby_aggregate = staticmethod(_value_counts_aggregate)


class Unique(SingleAggregation):
    groupby_chunk = M.unique
    groupby_aggregate = staticmethod(_unique_aggregate)

    @functools.cached_property
    def aggregate_kwargs(self) -> dict:  # type: ignore
        kwargs = super().aggregate_kwargs
        meta = self.frame._meta
        if meta.ndim == 1:
            name = meta.name
        else:
            name = meta[self._slice].name
        return {**kwargs, "name": name}


class Cov(SingleAggregation):
    chunk = staticmethod(_cov_chunk)
    std = False

    @classmethod
    def combine(cls, g, levels):
        return _concat(g)

    @classmethod
    def aggregate(cls, inputs, **kwargs):
        return _cov_agg(_concat(inputs), **kwargs)

    @property
    def chunk_kwargs(self) -> dict:  # type: ignore
        return self.operand("chunk_kwargs")

    @property
    def aggregate_kwargs(self) -> dict:  # type: ignore
        kwargs = self.operand("aggregate_kwargs").copy()
        kwargs["sort"] = self.sort
        kwargs["std"] = self.std
        kwargs["levels"] = self.levels
        return kwargs

    @property
    def combine_kwargs(self) -> dict:  # type: ignore
        return {"levels": self.levels}


class Corr(Cov):
    std = True


class GroupByReduction(Reduction, GroupByBase):
    _chunk_cls = GroupByChunk

    def _tune_down(self):
        if self.operand("split_out") is None:
            return self.substitute_parameters(
                {
                    "split_out": functools.partial(
                        _adjust_split_out_for_group_keys, by=self.by
                    )
                }
            )

    @property
    def split_out(self):
        if self.operand("split_out") is None:
            return 1
        return super().split_out

    @property
    def _chunk_cls_args(self):
        return self.by

    @functools.cached_property
    def _meta_chunk(self):
        meta = meta_nonempty(self.frame._meta)
        return self.chunk(meta, *self._by_meta, **self.chunk_kwargs)

    def _divisions(self):
        if self.sort:
            return (None, None)
        split_out = self.split_out
        if split_out is True:
            split_out = self.frame.npartitions
        return (None,) * (split_out + 1)

    def _simplify_up(self, parent, dependents):
        return groupby_projection(self, parent, dependents)

    @functools.cached_property
    def combine_kwargs(self):
        return {"levels": self.levels, "observed": self.observed, "dropna": self.dropna}

    @functools.cached_property
    def chunk_kwargs(self):
        return {"observed": self.observed, "dropna": self.dropna}

    @functools.cached_property
    def aggregate_kwargs(self):
        return {
            "levels": self.levels,
            "sort": self.sort,
            "observed": self.observed,
            "dropna": self.dropna,
        }


def _var_combine(g, levels, sort=False, observed=False, dropna=True):
    return g.groupby(level=levels, sort=sort, observed=observed, dropna=dropna).sum()


class Var(GroupByReduction):
    _parameters = [
        "frame",
        "ddof",
        "numeric_only",
        "split_out",
        "split_every",
        "sort",
        "dropna",
        "observed",
        "shuffle_method",
    ]
    _defaults = {
        "split_out": 1,
        "sort": None,
        "observed": None,
        "dropna": None,
        "split_every": None,
        "shuffle_method": None,
    }
    reduction_aggregate = staticmethod(_var_agg)
    reduction_combine = staticmethod(_var_combine)
    chunk = staticmethod(_var_chunk)

    @functools.cached_property
    def aggregate_kwargs(self):
        return {
            "ddof": self.ddof,
            "numeric_only": self.numeric_only,
            **super().aggregate_kwargs,
        }

    @functools.cached_property
    def chunk_kwargs(self):
        return {"numeric_only": self.numeric_only, **super().chunk_kwargs}


class Std(Var):
    def _lower(self):
        v = Var(*self.operands)
        return MapPartitions(
            v,
            func=np.sqrt,
            meta=v._meta,
            enforce_metadata=True,
            transform_divisions=True,
            clear_divisions=True,
        )


def _mean_chunk(df, *by, observed=None, dropna=None):
    if is_series_like(df):
        df = df.to_frame()

    g = _groupby_raise_unaligned(df, by=by, observed=observed, dropna=dropna)
    x = g.sum(numeric_only=True)
    n = g[x.columns].count().rename(columns=lambda c: c + "-count")
    return concat([x, n], axis=1)


def _mean_combine(g, levels, sort=False, observed=None, dropna=None):
    return g.groupby(level=levels, sort=sort, observed=observed, dropna=dropna).sum()


def _mean_agg(g, levels, sort=False, observed=False, dropna=True):
    result = g.groupby(level=levels, sort=sort, observed=observed, dropna=dropna).sum()
    s = result[result.columns[: len(result.columns) // 2]]
    c = result[result.columns[len(result.columns) // 2 :]]
    c.columns = s.columns
    return s / c


class Mean(GroupByReduction):
    _parameters = SingleAggregation._parameters
    _defaults = SingleAggregation._defaults
    reduction_aggregate = staticmethod(_mean_agg)
    reduction_combine = staticmethod(_mean_combine)
    chunk = staticmethod(_mean_chunk)


def nunique_df_combine(dfs, *args, **kwargs):
    return _nunique_df_combine(concat(dfs), *args, **kwargs)


def nunique_df_aggregate(dfs, levels, name, sort=False):
    df = concat(dfs)
    if df.ndim == 1:
        # split out reduces to a Series
        return df.groupby(level=levels, sort=sort, observed=True).nunique()
    else:
        return df.groupby(level=levels, sort=sort, observed=True)[name].nunique()


class NUnique(SingleAggregation):
    aggregate = staticmethod(nunique_df_aggregate)
    combine = staticmethod(nunique_df_combine)

    @staticmethod
    def chunk(df, *by, **kwargs):
        if df.ndim == 1:
            df = df.to_frame()
            kwargs = dict(name=df.columns[0], levels=_determine_levels(by))
        return _nunique_df_chunk(df, *by, **kwargs)

    @functools.cached_property
    def chunk_kwargs(self) -> dict:  # type: ignore
        kwargs = super().chunk_kwargs
        kwargs["name"] = self._slice
        return kwargs

    @functools.cached_property
    def aggregate_kwargs(self) -> dict:  # type: ignore
        return {"levels": self.levels, "name": self._slice}

    @functools.cached_property
    def combine_kwargs(self):
        return {"levels": self.levels}


class Head(SingleAggregation):
    groupby_chunk = staticmethod(_head_chunk)
    groupby_aggregate = staticmethod(_head_aggregate)

    @classmethod
    def combine(cls, inputs, **kwargs):
        return _concat(inputs)


class Tail(Head):
    groupby_chunk = staticmethod(_tail_chunk)
    groupby_aggregate = staticmethod(_tail_aggregate)


class GroupByApply(Expr, GroupByBase):
    _parameters = [
        "frame",
        "observed",
        "dropna",
        "_slice",
        "group_keys",
        "func",
        "meta",
        "args",
        "kwargs",
        "shuffle_method",
    ]
    _defaults = {
        "observed": None,
        "dropna": None,
        "_slice": None,
        "group_keys": True,
        "shuffle_method": None,
    }

    @functools.cached_property
    def grp_func(self):
        return functools.partial(groupby_slice_apply, func=self.func)

    @functools.cached_property
    def _meta(self):
        if self.operand("meta") is not no_default:
            return make_meta(self.operand("meta"), parent_meta=self.frame._meta)
        return _meta_apply_transform(self, self.grp_func)

    def _divisions(self):
        if self.need_to_shuffle:
            return (None,) * (self.frame.npartitions + 1)
        return self.frame.divisions

    def _shuffle_grp_func(self, shuffled=False):
        return self.grp_func

    @functools.cached_property
    def unique_partition_mapping_columns_from_shuffle(self):
        if not self.need_to_shuffle:
            return self.frame.unique_partition_mapping_columns_from_shuffle
        elif not any(isinstance(b, Expr) for b in self.by):
            return {tuple(self._by_columns)}
        else:
            return set()

    @functools.cached_property
    def need_to_shuffle(self):
        if not any(isinstance(b, Expr) for b in self.by):
            if any(
                set(self._by_columns) >= set(cols)
                for cols in self.frame.unique_partition_mapping_columns_from_shuffle
            ):
                return False

        return any(div is None for div in self.frame.divisions) or not any(
            _contains_index_name(self.frame._meta.index.name, b) for b in self.by
        )

    def _lower(self):
        df = self.frame
        by = self.by

        if self.need_to_shuffle:

            def get_map_columns(df):
                map_columns = {col: str(col) for col in df.columns if col != str(col)}
                unmap_columns = {v: k for k, v in map_columns.items()}
                return map_columns, unmap_columns

            # Map Tuple[str] column names to str before the shuffle

            if any(isinstance(b, Expr) for b in self.by):
                is_series = df.ndim == 1
                if is_series:
                    df = ToFrame(df)
                cols, assign_exprs = [], []
                for i, b in enumerate(self.by):
                    if isinstance(b, Expr):
                        assign_exprs.extend([f"_by_{i}", b])
                        cols.append(f"_by_{i}")
                if len(assign_exprs):
                    df = Assign(df, *assign_exprs)

                map_columns, unmap_columns = get_map_columns(df)
                if map_columns:
                    df = RenameFrame(df, map_columns)

                df = RearrangeByColumn(
                    df,
                    [map_columns.get(c, c) for c in cols],
                    df.npartitions,
                    method=self.shuffle_method,
                )

                if unmap_columns:
                    df = RenameFrame(df, unmap_columns)

                # pandas checks if the group keys are part of the initial
                # DataFrame through reference tracking. This is fine as long
                # as we don't trigger a copy after the Assign above, since the
                # blocks stay separate normally, disk shuffle triggers a copy
                # though, and we shouldn't rely on those internals in pandas,
                # so we can trigger a deep copy here to clear the references
                # since we know more about the query than pandas does.
                by = [
                    (
                        b
                        if not isinstance(b, Expr)
                        else _DeepCopy(
                            RenameSeries(
                                Projection(df, f"_by_{i}"), index=self.by[i].columns[0]
                            )
                        )
                    )
                    for i, b in enumerate(self.by)
                ]
                cols = [col for col in df.columns if col not in cols]
                if is_series:
                    cols = cols[0]
                df = Projection(df, cols)
            else:
                map_columns, unmap_columns = get_map_columns(df)
                if map_columns:
                    df = RenameFrame(df, map_columns)
                df = RearrangeByColumn(
                    df,
                    map_columns.get(self.by[0], self.by[0]),
                    self.npartitions,
                    method=self.shuffle_method,
                )

                if unmap_columns:
                    df = RenameFrame(df, unmap_columns)

            grp_func = self._shuffle_grp_func(True)
        else:
            grp_func = self._shuffle_grp_func(False)

        return GroupByUDFBlockwise(
            df,
            self._slice,
            self.group_keys,
            self.observed,
            self.dropna,
            self.operand("args"),
            self.operand("kwargs"),
            grp_func,
            self.operand("meta"),
            *by,
        )


class GroupByTransform(GroupByApply):
    @functools.cached_property
    def grp_func(self):
        return functools.partial(groupby_slice_transform, func=self.func)


def _fillna(group, *, what, **kwargs):
    return getattr(group, what)(**kwargs)


class GroupByBFill(GroupByTransform):
    func = staticmethod(functools.partial(_fillna, what="bfill"))

    def _simplify_up(self, parent, dependents):
        if isinstance(parent, Projection):
            return groupby_projection(self, parent, dependents)


class GroupByFFill(GroupByBFill):
    func = staticmethod(functools.partial(_fillna, what="ffill"))


class GroupByShift(GroupByApply):
    _defaults = {
        "observed": None,
        "dropna": None,
        "_slice": None,
        "func": None,
        "group_keys": True,
    }

    @functools.cached_property
    def grp_func(self):
        return functools.partial(groupby_slice_shift, shuffled=False)

    def _shuffle_grp_func(self, shuffled=False):
        return functools.partial(groupby_slice_shift, shuffled=shuffled)


class Median(GroupByShift):
    _parameters = GroupByApply._parameters + ["split_every"]
    default = {**GroupByShift._defaults, "split_every": None}

    @functools.cached_property
    def grp_func(self):
        return functools.partial(_median_groupby_aggregate)

    def _shuffle_grp_func(self, shuffled=False):
        return self.grp_func

    def _simplify_up(self, parent, dependents):
        if isinstance(parent, Projection):
            return groupby_projection(self, parent, dependents)

    @functools.cached_property
    def npartitions(self):
        npartitions = self.frame.npartitions
        if self.split_every is not None:
            npartitions = npartitions // self.split_every
        return npartitions


def groupby_get_group(df, *by_key, get_key=None, columns=None):
    if PANDAS_GE_300 and is_scalar(get_key):
        get_key = (get_key,)
    return _groupby_get_group(df, list(by_key), get_key, columns)


class GetGroup(Blockwise, GroupByBase):
    operation = staticmethod(groupby_get_group)
    _parameters = ["frame", "get_key", "columns"]
    _keyword_only = ["get_key", "columns"]

    @property
    def _args(self) -> list:
        return [self.frame] + self.by

    @property
    def _kwargs(self) -> dict:
        cols = self.operand("columns")
        return {
            "get_key": self.get_key,
            "columns": cols if cols is not None else self.frame.columns,
        }


def _median_groupby_aggregate(
    df,
    by=None,
    key=None,
    group_keys=True,  # not used
    dropna=None,
    observed=None,
    numeric_only=False,
    args=None,
    **kwargs,
):
    dropna = {"dropna": dropna} if dropna is not None else {}
    observed = {"observed": observed} if observed is not None else {}

    g = df.groupby(by=by, **observed, **dropna)
    if key is not None:
        g = g[key]
    return g.median(numeric_only=numeric_only)


class GroupByUDFBlockwise(Blockwise, GroupByBase):
    _parameters = [
        "frame",
        "_slice",
        "group_keys",
        "observed",
        "dropna",
        "args",
        "kwargs",
        "dask_func",
        "meta",
    ]
    _defaults = {"observed": None, "dropna": None}
    _keyword_only = [
        "_slice",
        "group_keys",
        "observed",
        "dropna",
        "args",
        "kwargs",
        "dask_func",
    ]

    @property
    def _args(self) -> list:
        return [self.frame] + self.by

    @functools.cached_property
    def _meta(self):
        if self.operand("meta") is not no_default:
            return make_meta(self.operand("meta"), parent_meta=self.frame._meta)
        return _meta_apply_transform(self, self.dask_func)

    def _task(self, name: Key, index: int) -> Task:
        args = [self._blockwise_arg(op, index) for op in self._args]
        kwargs = self._kwargs.copy()
        kwargs.update(
            {
                "_func": self.operation,
                "_meta": self._meta,
            }
        )
        return Task(name, apply_and_enforce, *args, **kwargs)

    @staticmethod
    def operation(
        frame,
        *by,
        _slice,
        group_keys=None,
        observed=None,
        dropna=None,
        args=None,
        kwargs=None,
        dask_func=None,
    ):
        if args is None:
            args = ()
        if kwargs is None:
            kwargs = {}
        return dask_func(
            frame,
            list(by),
            key=_slice,
            group_keys=group_keys,
            args=args,
            **_as_dict("observed", observed),
            **_as_dict("dropna", dropna),
            **kwargs,
        )


def _contains_index_name(index_name, by):
    if index_name is None:
        return False

    if isinstance(by, Expr):
        return False

    if not is_scalar(by):
        return False

    return index_name == by


def groupby_slice_apply(
    df,
    grouper,
    key,
    func,
    args,
    group_keys=GROUP_KEYS_DEFAULT,
    dropna=None,
    observed=None,
    **kwargs,
):
    return _groupby_slice_apply(
        df,
        grouper,
        key,
        func,
        *args,
        group_keys=group_keys,
        dropna=dropna,
        observed=observed,
        **kwargs,
    )


def groupby_slice_shift(
    df,
    grouper,
    key,
    args,
    shuffled,
    group_keys=GROUP_KEYS_DEFAULT,
    dropna=None,
    observed=None,
    **kwargs,
):
    return _groupby_slice_shift(
        df, grouper, key, shuffled, group_keys, dropna, observed, **kwargs
    )


def groupby_slice_transform(
    df,
    grouper,
    key,
    func,
    args,
    group_keys=GROUP_KEYS_DEFAULT,
    dropna=None,
    observed=None,
    **kwargs,
):
    return _groupby_slice_transform(
        df,
        grouper,
        key,
        func,
        *args,
        group_keys=group_keys,
        dropna=dropna,
        observed=observed,
        **kwargs,
    )


def _meta_apply_transform(obj, grp_func):
    kwargs = obj.operand("kwargs")
    by_meta = obj._by_meta
    by_meta = [x if is_scalar(x) else meta_nonempty(x) for x in by_meta]
    meta_args, meta_kwargs = _extract_meta((obj.operand("args"), kwargs), nonempty=True)
    return make_meta(
        grp_func(
            meta_nonempty(obj.frame._meta),
            by_meta,
            key=obj._slice,
            args=meta_args,
            **_as_dict("observed", obj.observed),
            **_as_dict("dropna", obj.dropna),
            **_as_dict("group_keys", obj.group_keys),
            **meta_kwargs,
        )
    )


def groupby_projection(expr, parent, dependents):
    if isinstance(parent, Projection):
        columns = determine_column_projection(
            expr, parent, dependents, additional_columns=expr._by_columns
        )
        columns = _convert_to_list(columns)
        columns = [col for col in expr.frame.columns if col in columns]
        if columns == expr.frame.columns:
            return
        return type(parent)(
            type(expr)(expr.frame[columns], *expr.operands[1:]),
            *parent.operands[1:],
        )
    return


###
### Groupby Collection API
###


def _clean_by_expr(obj, by):
    if (
        isinstance(by, Series)
        and by.name in obj.columns
        and obj.ndim == 2
        and by._name == obj[by.name]._name
    ):
        return by.name
    elif isinstance(by, Index) and by._name == obj.index._name:
        return by.expr
    elif isinstance(by, (Series, Index)):
        if not are_co_aligned(obj.expr, by.expr):
            raise NotImplementedError(
                "by must be in the DataFrames columns or aligned with the DataFrame."
            )
        if isinstance(by, Index):
            by = by.to_series()
            by.index = obj.index
        return by.expr

    # By is a column name, e.g. str or int
    return by


class GroupByCumulative(Expr, GroupByBase):
    _parameters = ["frame", "dropna", "_slice", "numeric_only"]
    _defaults = {"numeric_only": None, "dropna": None, "_slice": None}
    chunk = None
    aggregate: Callable | None = None
    initial = 0

    @functools.cached_property
    def _meta(self):
        cols = None if self._slice is None else self._slice
        return _apply_chunk(
            self.frame._meta,
            *self._by_meta,
            chunk=self.chunk,
            columns=cols,
            **self.numeric_only,
        )

    def _divisions(self):
        return self.frame.divisions

    @property
    def numeric_only(self):
        no = self.operand("numeric_only")
        return {} if no is None else {"numeric_only": no}

    def _lower(self):
        meta = self._meta
        dropna = {} if self.dropna is None else {"dropna": self.dropna}
        columns = meta.name if is_series_like(meta) else meta.columns

        frame = MapPartitions(
            self.frame,
            _apply_chunk,
            meta,
            True,
            True,
            False,
            True,
            None,
            None,
            None,
            {"chunk": self.chunk, "columns": columns, **dropna, **self.numeric_only},
            len(self.by),
            *self.by,
        )
        cum_raw = frame

        if frame.ndim == 1:
            frame = ToFrame(frame)

        by = self.by.copy()
        for i, b in enumerate(by):
            if not isinstance(b, Expr):
                if b in self.frame.columns:
                    frame = Assign(frame, f"_by_{b}", self.frame[b])
                else:
                    frame = Assign(frame, f"_by_{b}", self.frame.index)

                by[i] = f"_by_{b}"

        columns = 0 if columns is None else columns
        cum_last = MapPartitions(
            frame,
            _apply_chunk,
            no_default,
            True,
            True,
            False,
            True,
            None,
            None,
            None,
            {"chunk": M.last, "columns": columns, **dropna},
            len(by),
            *by,
        )
        return GroupByCumulativeFinalizer(
            frame,
            cum_raw,
            cum_last,
            meta,
            self.aggregate,
            self.initial,
            columns,
            *by,
        )


class GroupByCumulativeFinalizer(Expr, GroupByBase):
    _parameters = [
        "frame",
        "cum_raw",
        "cum_last",
        "meta",
        "aggregate",
        "initial",
        "columns",
    ]

    @functools.cached_property
    def _meta(self):
        return self.meta

    def _divisions(self):
        return self.frame.divisions

    def _layer(self) -> dict:
        dsk = {(self._name, 0): (self.cum_raw._name, 0)}
        name_cum = "cum-last" + self._name

        for i in range(1, self.frame.npartitions):
            # store each cumulative step to graph to reduce computation
            if i == 1:
                dsk[(name_cum, i)] = (self.cum_last._name, i - 1)
            else:
                # aggregate with previous cumulation results
                dsk[(name_cum, i)] = (  # type: ignore
                    _cum_agg_filled,
                    (name_cum, i - 1),
                    (self.cum_last._name, i - 1),
                    self.aggregate,
                    self.initial,
                )
            dsk[(self._name, i)] = (  # type: ignore
                _cum_agg_aligned,
                (self.frame._name, i),
                (name_cum, i),
                self.by,
                self.operand("columns"),
                self.aggregate,
                self.initial,
            )
        return dsk


class GroupByCumsum(GroupByCumulative):
    chunk = M.cumsum
    aggregate = M.add
    initial = 0


class GroupByCumprod(GroupByCumulative):
    chunk = M.cumprod
    aggregate = M.mul
    initial = 1


class GroupByCumcount(GroupByCumulative):
    chunk = M.cumcount
    aggregate = staticmethod(_cumcount_aggregate)
    initial = -1


class GroupBy:
    """Collection container for groupby aggregations

    The purpose of this class is to expose an API similar
    to Pandas' `Groupby` for dask-expr collections.

    See Also
    --------
    SingleAggregation
    """

    def __init__(
        self,
        obj,
        by,
        group_keys=True,
        sort=None,
        observed=None,
        dropna=None,
        slice=None,
    ):
        if isinstance(by, (tuple, list)):
            by = [_clean_by_expr(obj, x) for x in by]
        else:
            by = _clean_by_expr(obj, by)

        by_ = by if isinstance(by, (tuple, list)) else [by]
        if any(isinstance(key, pd.Grouper) for key in by_):
            raise NotImplementedError("pd.Grouper is currently not supported by Dask.")
        self._slice = slice
        # Check if we can project columns
        projection = None
        if (
            np.isscalar(slice)
            or isinstance(slice, (str, list, tuple))
            or (
                (is_index_like(slice) or is_series_like(slice))
                and not is_dask_collection(slice)
            )
        ):
            projection = set(by_).union(
                {slice} if (np.isscalar(slice) or isinstance(slice, str)) else slice
            )
            projection = [c for c in obj.columns if c in projection]

        self.obj = obj[projection] if projection is not None else obj
        self.sort = sort
        self.observed = (
            observed if observed is not None else False if not PANDAS_GE_300 else True
        )
        self.dropna = dropna
        self.group_keys = group_keys
        self.by = (
            [by] if np.isscalar(by) or isinstance(by, (Expr, Callable)) else list(by)
        )
        # surface pandas errors
        self._meta = self.obj._meta.groupby(
            by,
            group_keys=group_keys,
            sort=sort,
            **_as_dict("observed", observed),
            **_as_dict("dropna", dropna),
        )
        if slice is not None:
            if isinstance(slice, tuple):
                slice = list(slice)
            self._meta = self._meta[slice]

    def _numeric_only_kwargs(self, numeric_only):
        kwargs = {"numeric_only": numeric_only}
        return {"chunk_kwargs": kwargs.copy(), "aggregate_kwargs": kwargs.copy()}

    def _single_agg(
        self,
        expr_cls,
        split_every=None,
        split_out=None,
        chunk_kwargs=None,
        aggregate_kwargs=None,
        shuffle_method=None,
    ):
        if split_every is None:
            split_every = 8
        return new_collection(
            expr_cls(
                self.obj.expr,
                self.observed,
                self.dropna,
                chunk_kwargs,
                aggregate_kwargs,
                self._slice,
                split_every,
                split_out,
                self.sort,
                get_specified_shuffle(shuffle_method),
                *self.by,
            )
        )

    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError as e:
            raise AttributeError(e) from e

    def __dir__(self):
        return sorted(
            set(
                dir(type(self))
                + list(self.__dict__)
                + list(filter(M.isidentifier, self.obj.columns))
            )
        )

    def compute(self, **kwargs):
        raise NotImplementedError(
            "DataFrameGroupBy does not allow compute method."
            "Please chain it with an aggregation method (like ``.mean()``) or get a "
            "specific group using ``.get_group()`` before calling ``compute()``"
        )

    def __getitem__(self, key):
        if is_scalar(key):
            return SeriesGroupBy(
                self.obj,
                by=self.by,
                group_keys=self.group_keys,
                slice=key,
                sort=self.sort,
                dropna=self.dropna,
                observed=self.observed,
            )
        g = GroupBy(
            self.obj,
            by=self.by,
            slice=key,
            sort=self.sort,
            dropna=self.dropna,
            observed=self.observed,
            group_keys=self.group_keys,
        )
        return g

    @derived_from(
        pd.core.groupby.GroupBy,
        inconsistencies="If the group is not present, Dask will return an empty Series/DataFrame.",
    )
    def get_group(self, key):
        return new_collection(GetGroup(self.obj.expr, key, self._slice, *self.by))

    @derived_from(pd.core.groupby.GroupBy)
    def count(self, **kwargs):
        return self._single_agg(Count, **kwargs)

    @derived_from(pd.core.groupby.GroupBy)
    def sum(self, numeric_only=False, min_count=None, **kwargs):
        numeric_kwargs = self._numeric_only_kwargs(numeric_only)
        result = self._single_agg(Sum, **kwargs, **numeric_kwargs)
        if min_count:
            return result.where(self.count() >= min_count, other=np.nan)
        return result

    @derived_from(pd.core.groupby.GroupBy)
    def prod(self, numeric_only=False, min_count=None, **kwargs):
        numeric_kwargs = self._numeric_only_kwargs(numeric_only)
        result = self._single_agg(Prod, **kwargs, **numeric_kwargs)
        if min_count:
            return result.where(self.count() >= min_count, other=np.nan)
        return result

    def _cum_agg(self, cls, numeric_only=None):
        return new_collection(
            cls(
                self.obj.expr,
                self.dropna,
                self._slice,
                numeric_only,
                *self.by,
            )
        )

    @derived_from(pd.core.groupby.GroupBy)
    def cumsum(self, numeric_only=False):
        return self._cum_agg(GroupByCumsum, numeric_only)

    @derived_from(pd.core.groupby.GroupBy)
    def cumprod(self, numeric_only=False):
        return self._cum_agg(GroupByCumprod, numeric_only)

    @derived_from(pd.core.groupby.GroupBy)
    def cumcount(self):
        return self._cum_agg(GroupByCumcount)

    def _all_numeric(self):
        """Are all columns that we're not grouping on numeric?"""
        numerics = self.obj._meta._get_numeric_data()
        # This computes a groupby but only on the empty meta
        post_group_columns = self._meta.count().columns
        return len(set(post_group_columns) - set(numerics.columns)) == 0

    @derived_from(pd.core.groupby.GroupBy)
    def mean(self, numeric_only=False, split_out=None, **kwargs):
        if not numeric_only and not self._all_numeric():
            raise NotImplementedError(
                "'numeric_only=False' is not implemented in Dask."
            )
        numeric_kwargs = self._numeric_only_kwargs(numeric_only)
        result = self._single_agg(Mean, split_out=split_out, **kwargs, **numeric_kwargs)
        return self._postprocess_series_squeeze(result)

    def _postprocess_series_squeeze(self, result):
        if (
            isinstance(self.obj, Series)
            or is_scalar(self._slice)
            and self._slice is not None
        ):
            if len(result.columns) < 1:
                raise NotImplementedError(
                    "Cannot call `SeriesGroupBy.var` or `SeriesGroupBy.mean` on the key "
                    "column. Please use `aggregate` if you really need to do this."
                )
            result = result[result.columns[0]]
        return result

    @derived_from(pd.core.groupby.GroupBy)
    def min(self, numeric_only=False, **kwargs):
        numeric_kwargs = self._numeric_only_kwargs(numeric_only)
        return self._single_agg(Min, **kwargs, **numeric_kwargs)

    @derived_from(pd.core.groupby.GroupBy)
    def max(self, numeric_only=False, **kwargs):
        numeric_kwargs = self._numeric_only_kwargs(numeric_only)
        return self._single_agg(Max, **kwargs, **numeric_kwargs)

    @derived_from(pd.core.groupby.GroupBy)
    def first(self, numeric_only=False, sort=None, **kwargs):
        if sort:
            raise NotImplementedError()
        numeric_kwargs = self._numeric_only_kwargs(numeric_only)
        return self._single_agg(First, **kwargs, **numeric_kwargs)

    @derived_from(pd.DataFrame)
    def cov(
        self,
        ddof=1,
        split_every=None,
        split_out=None,
        numeric_only=False,
        shuffle_method=None,
    ):
        numeric_kwargs = self._numeric_only_kwargs(numeric_only)
        return self._single_agg(
            Cov,
            split_every,
            split_out,
            chunk_kwargs=numeric_kwargs["chunk_kwargs"],
            aggregate_kwargs={"ddof": ddof},
        )

    @derived_from(pd.DataFrame)
    def corr(
        self, split_every=None, split_out=None, numeric_only=False, shuffle_method=None
    ):
        numeric_kwargs = self._numeric_only_kwargs(numeric_only)
        return self._single_agg(
            Corr,
            split_every,
            split_out,
            chunk_kwargs=numeric_kwargs["chunk_kwargs"],
            aggregate_kwargs={"ddof": 1},
        )

    @derived_from(pd.core.groupby.GroupBy)
    def last(self, numeric_only=False, sort=None, **kwargs):
        if sort:
            raise NotImplementedError()
        numeric_kwargs = self._numeric_only_kwargs(numeric_only)
        return self._single_agg(Last, **kwargs, **numeric_kwargs)

    @derived_from(pd.core.groupby.GroupBy)
    def ffill(self, limit=None, shuffle_method=None):
        return self._transform_like_op(
            GroupByFFill, None, limit=limit, shuffle_method=shuffle_method
        )

    @derived_from(pd.core.groupby.GroupBy)
    def bfill(self, limit=None, shuffle_method=None):
        return self._transform_like_op(
            GroupByBFill, None, limit=limit, shuffle_method=shuffle_method
        )

    @derived_from(pd.core.groupby.GroupBy)
    def size(self, **kwargs):
        return self._single_agg(Size, **kwargs)

    @derived_from(pd.DataFrame)
    def idxmin(
        self,
        split_every=None,
        split_out=None,
        skipna=True,
        numeric_only=False,
        shuffle_method=None,
    ):
        numeric_kwargs = self._numeric_only_kwargs(numeric_only)
        numeric_kwargs["chunk_kwargs"]["skipna"] = skipna
        return self._single_agg(
            IdxMin,
            split_every=split_every,
            split_out=split_out,
            shuffle_method=shuffle_method,
            **numeric_kwargs,
        )

    @derived_from(pd.DataFrame)
    def idxmax(
        self,
        split_every=None,
        split_out=None,
        skipna=True,
        numeric_only=False,
        shuffle_method=None,
    ):
        numeric_kwargs = self._numeric_only_kwargs(numeric_only)
        numeric_kwargs["chunk_kwargs"]["skipna"] = skipna
        return self._single_agg(
            IdxMax,
            split_every=split_every,
            split_out=split_out,
            shuffle_method=shuffle_method,
            **numeric_kwargs,
        )

    @derived_from(pd.core.groupby.SeriesGroupBy)
    def head(self, n=5, split_every=None, split_out=None):
        chunk_kwargs = {"n": n}
        aggregate_kwargs = {
            "n": n,
            "index_levels": len(self.by) if not isinstance(self.by, Expr) else 1,
        }
        return self._single_agg(
            Head,
            split_every=split_every,
            split_out=split_out,
            chunk_kwargs=chunk_kwargs,
            aggregate_kwargs=aggregate_kwargs,
        )

    @derived_from(pd.core.groupby.SeriesGroupBy)
    def tail(self, n=5, split_every=None, split_out=None):
        chunk_kwargs = {"n": n}
        aggregate_kwargs = {
            "n": n,
            "index_levels": len(self.by) if not isinstance(self.by, Expr) else 1,
        }
        return self._single_agg(
            Tail,
            split_every=split_every,
            split_out=split_out,
            chunk_kwargs=chunk_kwargs,
            aggregate_kwargs=aggregate_kwargs,
        )

    @derived_from(pd.core.groupby.GroupBy)
    def var(
        self,
        ddof=1,
        split_every=None,
        split_out=None,
        numeric_only=False,
        shuffle_method=None,
    ):
        if not numeric_only and not self._all_numeric():
            raise NotImplementedError(
                "'numeric_only=False' is not implemented in Dask."
            )
        result = new_collection(
            Var(
                self.obj.expr,
                ddof,
                numeric_only,
                split_out,
                split_every,
                self.sort,
                self.dropna,
                self.observed,
                shuffle_method,
                *self.by,
            )
        )
        return self._postprocess_series_squeeze(result)

    @derived_from(pd.core.groupby.GroupBy)
    def std(
        self,
        ddof=1,
        split_every=None,
        split_out=None,
        numeric_only=False,
        shuffle_method=None,
    ):
        if not numeric_only and not self._all_numeric():
            raise NotImplementedError(
                "'numeric_only=False' is not implemented in Dask."
            )
        result = new_collection(
            Std(
                self.obj.expr,
                ddof,
                numeric_only,
                split_out,
                split_every,
                self.sort,
                self.dropna,
                self.observed,
                shuffle_method,
                *self.by,
            )
        )
        return self._postprocess_series_squeeze(result)

    @_aggregate_docstring(based_on="pd.core.groupby.DataFrameGroupBy.agg")
    def aggregate(
        self, arg=None, split_every=8, split_out=None, shuffle_method=None, **kwargs
    ):
        relabeling, order, columns = None, None, None
        if arg is None:
            if not isinstance(self, SeriesGroupBy):
                relabeling, arg, columns, order = reconstruct_func(arg, **kwargs)
            elif isinstance(self, SeriesGroupBy):
                columns, arg = validate_func_kwargs(kwargs)
                relabeling = True

        if arg == "size":
            return self.size()

        result = new_collection(
            GroupbyAggregation(
                self.obj.expr,
                arg,
                self.observed,
                self.dropna,
                split_every,
                split_out,
                self.sort,
                shuffle_method,
                self._slice,
                *self.by,
            )
        )
        if relabeling and result is not None:
            if order is not None:
                result = result.iloc[:, order]
            result.columns = columns
        return result

    def agg(self, *args, **kwargs):
        return self.aggregate(*args, **kwargs)

    def _warn_if_no_meta(self, meta, method="apply"):
        if meta is no_default:
            msg = f"""`meta` is not specified, inferred from partial data.
Please provide `meta` if the result is unexpected.
  Before: .{method}(func)
  After:  .{method}(func, meta={{'x': 'f8', 'y': 'f8'}}) for dataframe result
  or:     .{method}(func, meta=('x', 'f8'))            for series result
"""
            warnings.warn(msg, stacklevel=3)

    @insert_meta_param_description(pad=12)
    def apply(self, func, *args, meta=no_default, shuffle_method=None, **kwargs):
        """Parallel version of pandas GroupBy.apply

        This mimics the pandas version except for the following:

        1.  If the grouper does not align with the index then this causes a full
            shuffle.  The order of rows within each group may not be preserved.
        2.  Dask's GroupBy.apply is not appropriate for aggregations. For custom
            aggregations, use :class:`dask.dataframe.groupby.Aggregation`.

        .. warning::

           Pandas' groupby-apply can be used to to apply arbitrary functions,
           including aggregations that result in one row per group. Dask's
           groupby-apply will apply ``func`` once on each group, doing a shuffle
           if needed, such that each group is contained in one partition.
           When ``func`` is a reduction, e.g., you'll end up with one row
           per group. To apply a custom aggregation with Dask,
           use :class:`dask.dataframe.groupby.Aggregation`.

        Parameters
        ----------
        func: function
            Function to apply
        args, kwargs : Scalar, Delayed or object
            Arguments and keywords to pass to the function.
        $META

        Returns
        -------
        applied : Series or DataFrame depending on columns keyword
        """
        self._warn_if_no_meta(meta)
        return new_collection(
            GroupByApply(
                self.obj.expr,
                self.observed,
                self.dropna,
                self._slice,
                self.group_keys,
                func,
                meta,
                args,
                kwargs,
                get_specified_shuffle(shuffle_method),
                *self.by,
            )
        )

    def _transform_like_op(
        self, expr_cls, func, meta=no_default, shuffle_method=None, *args, **kwargs
    ):
        return new_collection(
            expr_cls(
                self.obj.expr,
                self.observed,
                self.dropna,
                self._slice,
                self.group_keys,
                func,
                meta,
                args,
                kwargs,
                get_specified_shuffle(shuffle_method),
                *self.by,
            )
        )

    @insert_meta_param_description(pad=12)
    def transform(self, func, meta=no_default, shuffle_method=None, *args, **kwargs):
        """Parallel version of pandas GroupBy.transform

        This mimics the pandas version except for the following:

        1.  If the grouper does not align with the index then this causes a full
            shuffle.  The order of rows within each group may not be preserved.
        2.  Dask's GroupBy.transform is not appropriate for aggregations. For custom
            aggregations, use :class:`dask.dataframe.groupby.Aggregation`.

        .. warning::

           Pandas' groupby-transform can be used to apply arbitrary functions,
           including aggregations that result in one row per group. Dask's
           groupby-transform will apply ``func`` once on each group, doing a shuffle
           if needed, such that each group is contained in one partition.
           When ``func`` is a reduction, e.g., you'll end up with one row
           per group. To apply a custom aggregation with Dask,
           use :class:`dask.dataframe.groupby.Aggregation`.

        Parameters
        ----------
        func: function
            Function to apply
        args, kwargs : Scalar, Delayed or object
            Arguments and keywords to pass to the function.
        $META

        Returns
        -------
        applied : Series or DataFrame depending on columns keyword
        """
        self._warn_if_no_meta(meta, method="transform")
        return self._transform_like_op(
            GroupByTransform, func, meta, shuffle_method, *args, **kwargs
        )

    @insert_meta_param_description(pad=12)
    def shift(self, periods=1, meta=no_default, shuffle_method=None, *args, **kwargs):
        """Parallel version of pandas GroupBy.shift

        This mimics the pandas version except for the following:

        If the grouper does not align with the index then this causes a full
        shuffle.  The order of rows within each group may not be preserved.

        Parameters
        ----------
        periods : Delayed, Scalar or int, default 1
            Number of periods to shift.
        freq : Delayed, Scalar or str, optional
            Frequency string.
        fill_value : Scalar, Delayed or object, optional
            The scalar value to use for newly introduced missing values.
        $META

        Returns
        -------
        shifted : Series or DataFrame shifted within each group.

        Examples
        --------
        >>> import dask
        >>> ddf = dask.datasets.timeseries(freq="1h")
        >>> result = ddf.groupby("name").shift(1, meta={"id": int, "x": float, "y": float})
        """
        if "axis" in kwargs:
            raise TypeError("axis is not supported in shift.")
        self._warn_if_no_meta(meta, method="shift")
        kwargs = {"periods": periods, **kwargs}
        return self._transform_like_op(
            GroupByShift, None, meta, shuffle_method, *args, **kwargs
        )

    @derived_from(pd.core.groupby.GroupBy)
    def median(
        self, split_every=None, split_out=True, shuffle_method=None, numeric_only=False
    ):
        result = new_collection(
            Median(
                self.obj.expr,
                self.observed,
                self.dropna,
                self._slice,
                self.group_keys,
                None,
                no_default,
                (),
                {"numeric_only": numeric_only},
                get_specified_shuffle(shuffle_method),
                split_every,
                *self.by,
            )
        )
        if split_out is not True:
            result = result.repartition(npartitions=split_out)
        return result

    def rolling(self, window, min_periods=None, center=False, win_type=None, axis=0):
        """Provides rolling transformations.

        .. note::

            Since MultiIndexes are not well supported in Dask, this method returns a
            dataframe with the same index as the original data. The groupby column is
            not added as the first level of the index like pandas does.

            This method works differently from other groupby methods. It does a groupby
            on each partition (plus some overlap). This means that the output has the
            same shape and number of partitions as the original.

        Parameters
        ----------
        window : str, offset
           Size of the moving window. This is the number of observations used
           for calculating the statistic. Data must have a ``DatetimeIndex``
        min_periods : int, default None
            Minimum number of observations in window required to have a value
            (otherwise result is NA).
        center : boolean, default False
            Set the labels at the center of the window.
        win_type : string, default None
            Provide a window type. The recognized window types are identical
            to pandas.
        axis : int, default 0

        Returns
        -------
        a Rolling object on which to call a method to compute a statistic

        Examples
        --------
        >>> import dask
        >>> ddf = dask.datasets.timeseries(freq="1h")
        >>> result = ddf.groupby("name").x.rolling('1D').max()
        """
        from dask.dataframe.dask_expr._rolling import Rolling

        return Rolling(
            self.obj,
            window,
            min_periods=min_periods,
            center=center,
            win_type=win_type,
            groupby_kwargs={
                "by": self.by,
                "sort": self.sort,
                "observed": self.observed,
                "dropna": self.dropna,
                "group_keys": self.group_keys,
            },
            groupby_slice=self._slice,
        )


class SeriesGroupBy(GroupBy):
    def __init__(
        self,
        obj,
        by,
        group_keys=True,
        sort=None,
        observed=None,
        dropna=None,
        slice=None,
    ):
        # Raise pandas errors if applicable
        if isinstance(obj, Series):
            if isinstance(by, FrameBase):
                obj._meta.groupby(by._meta, **_as_dict("observed", observed))
            elif isinstance(by, (list, tuple)) and any(
                isinstance(x, FrameBase) for x in by
            ):
                metas = [x._meta if isinstance(x, FrameBase) else x for x in by]
                obj._meta.groupby(metas, **_as_dict("observed", observed))
            elif isinstance(by, list):
                if len(by) == 0:
                    raise ValueError("No group keys passed!")

                non_series_items = [item for item in by if not isinstance(item, Series)]
                obj._meta.groupby(non_series_items, **_as_dict("observed", observed))
            else:
                obj._meta.groupby(by, **_as_dict("observed", observed))

        super().__init__(
            obj,
            by=by,
            group_keys=group_keys,
            slice=slice,
            observed=observed,
            dropna=dropna,
            sort=sort,
        )

    @derived_from(pd.core.groupby.SeriesGroupBy)
    def value_counts(self, **kwargs):
        return self._single_agg(ValueCounts, **kwargs)

    @derived_from(pd.core.groupby.SeriesGroupBy)
    def unique(self, **kwargs):
        return self._single_agg(Unique, **kwargs)

    def idxmin(
        self,
        split_every=None,
        split_out=None,
        skipna=True,
        numeric_only=False,
        **kwargs,
    ):
        # pandas doesn't support numeric_only here, which is odd
        return self._single_agg(
            IdxMin,
            split_every=None,
            split_out=split_out,
            chunk_kwargs=dict(skipna=skipna),
        )

    def idxmax(
        self,
        split_every=None,
        split_out=None,
        skipna=True,
        numeric_only=False,
        **kwargs,
    ):
        # pandas doesn't support numeric_only here, which is odd
        return self._single_agg(
            IdxMax,
            split_every=split_every,
            split_out=split_out,
            chunk_kwargs=dict(skipna=skipna),
        )

    @derived_from(pd.core.groupby.SeriesGroupBy)
    def nunique(self, split_every=None, split_out=True, shuffle_method=None):
        """
        Examples
        --------
        >>> import pandas as pd
        >>> import dask.dataframe as dd
        >>> d = {'col1': [1, 2, 3, 4], 'col2': [5, 6, 7, 8]}
        >>> df = pd.DataFrame(data=d)
        >>> ddf = dd.from_pandas(df, 2)
        >>> ddf.groupby(['col1']).col2.nunique().compute()
        """
        slice = self._slice or self.obj.name
        return new_collection(
            NUnique(
                self.obj.expr,
                self.observed,
                self.dropna,
                None,
                None,
                slice,
                split_every,
                split_out,
                self.sort,
                get_specified_shuffle(shuffle_method),
                *self.by,
            )
        )

    def cov(self, *args, **kwargs):
        raise NotImplementedError("cov is not implemented for SeriesGroupBy objects.")

    def corr(self, *args, **kwargs):
        raise NotImplementedError("cov is not implemented for SeriesGroupBy objects.")

    def _all_numeric(self):
        return True
