from __future__ import annotations

import datetime
import functools
import numbers
import operator
import warnings
import weakref
from collections import defaultdict
from collections.abc import Callable, Collection, Mapping
from functools import partial
from typing import Any as AnyType

import numpy as np
import pandas as pd
from pandas.errors import PerformanceWarning
from tlz import merge_sorted, partition, unique

from dask import _expr as core
from dask._expr import Expr as BaseExpr
from dask._expr import FinalizeCompute
from dask._task_spec import Alias, DataNode, Task, TaskRef, execute_graph
from dask.array import Array
from dask.base import collections_to_expr
from dask.core import flatten
from dask.dataframe import methods
from dask.dataframe._pyarrow import to_pyarrow_string
from dask.dataframe.core import (
    _concat,
    _get_divisions_map_partitions,
    _rename,
    apply_and_enforce,
    has_parallel_type,
    is_dataframe_like,
    is_index_like,
    is_series_like,
    safe_head,
    total_mem_usage,
)
from dask.dataframe.dask_expr._util import (
    _BackendData,
    _calc_maybe_new_divisions,
    _convert_to_list,
    _tokenize_partial,
)
from dask.dataframe.dispatch import make_meta, meta_nonempty
from dask.dataframe.rolling import CombinedOutput, _head_timedelta, overlap_chunk
from dask.dataframe.shuffle import drop_overlap, get_overlap
from dask.dataframe.utils import (
    clear_known_categories,
    drop_by_shallow_copy,
    is_scalar,
    raise_on_meta_error,
    valid_divisions,
)
from dask.typing import Key, no_default
from dask.utils import (
    M,
    funcname,
    get_meta_library,
    has_keyword,
    is_arraylike,
    partial_by_order,
    pseudorandom,
    random_state_data,
)

optimize = core.optimize


class Expr(core.SingletonExpr):
    """Primary class for all Expressions

    This mostly includes Dask protocols and various Pandas-like method
    definitions to make us look more like a DataFrame.
    """

    _is_length_preserving = False
    _filter_passthrough = False

    def _filter_passthrough_available(self, parent, dependents):
        return self._filter_passthrough and is_filter_pushdown_available(
            self, parent, dependents
        )

    @functools.cached_property
    def ndim(self):
        meta = self._meta
        try:
            return meta.ndim
        except AttributeError:
            return 0

    def __dask_keys__(self):
        return [(self._name, i) for i in range(self.npartitions)]

    def optimize(self, **kwargs):
        return optimize(self, **kwargs)

    def __hash__(self):
        return hash(self._name)

    @property
    def index(self):
        return Index(self)

    @property
    def size(self):
        return Size(self)

    @property
    def nbytes(self):
        return NBytes(self)

    def _tree_repr_lines(self, indent=0, recursive=True):
        header = funcname(type(self)) + ":"
        lines = []
        for i, op in enumerate(self.operands):
            if isinstance(op, Expr):
                if recursive:
                    lines.extend(op._tree_repr_lines(2))
            else:
                if isinstance(op, _BackendData):
                    op = op._data

                # TODO: this stuff is pandas-specific
                if isinstance(op, pd.core.base.PandasObject):
                    op = "<pandas>"
                elif is_dataframe_like(op):
                    op = "<dataframe>"
                elif is_index_like(op):
                    op = "<index>"
                elif is_series_like(op):
                    op = "<series>"
                elif is_arraylike(op):
                    op = "<array>"
                header = self._tree_repr_argument_construction(i, op, header)

        lines = [header] + lines
        lines = [" " * indent + line for line in lines]

        return lines

    def _operands_for_repr(self):
        to_include = []
        for param, operand in zip(self._parameters, self.operands):
            if isinstance(operand, Expr) or (
                not isinstance(operand, (pd.Series, pd.DataFrame))
                and operand != self._defaults.get(param)
            ):
                to_include.append(f"{param}={operand!r}")
        return to_include

    def __getattr__(self, key):
        try:
            return super().__getattr__(key)
        except AttributeError:

            if is_dataframe_like(self._meta) and key in self._meta.columns:
                return self[key]
            raise

    def __getitem__(self, other):
        if isinstance(other, Expr):
            return Filter(self, other)
        else:
            return Projection(self, other)  # df[["a", "b", "c"]]

    def __bool__(self):
        raise ValueError(
            f"The truth value of a {self.__class__.__name__} is ambiguous. "
            "Use a.any() or a.all()."
        )

    def __add__(self, other):
        return Add(self, other)

    def __radd__(self, other):
        return Add(other, self)

    def __sub__(self, other):
        return Sub(self, other)

    def __rsub__(self, other):
        return Sub(other, self)

    def __mul__(self, other):
        return Mul(self, other)

    def __rmul__(self, other):
        return Mul(other, self)

    def __pow__(self, power):
        return Pow(self, power)

    def __rpow__(self, power):
        return Pow(power, self)

    def __truediv__(self, other):
        return Div(self, other)

    def __rtruediv__(self, other):
        return Div(other, self)

    def __lt__(self, other):
        return LT(self, other)

    def __rlt__(self, other):
        return LT(other, self)

    def __gt__(self, other):
        return GT(self, other)

    def __rgt__(self, other):
        return GT(other, self)

    def __le__(self, other):
        return LE(self, other)

    def __rle__(self, other):
        return LE(other, self)

    def __ge__(self, other):
        return GE(self, other)

    def __rge__(self, other):
        return GE(other, self)

    def __eq__(self, other):
        return EQ(self, other)

    def __ne__(self, other):
        return NE(self, other)

    def __and__(self, other):
        return And(self, other)

    def __rand__(self, other):
        return And(other, self)

    def __or__(self, other):
        return Or(self, other)

    def __ror__(self, other):
        return Or(other, self)

    def __xor__(self, other):
        return XOr(self, other)

    def __rxor__(self, other):
        return XOr(other, self)

    def __invert__(self):
        return Invert(self)

    def __neg__(self):
        return Neg(self)

    def __pos__(self):
        return Pos(self)

    def __mod__(self, other):
        return Mod(self, other)

    def __rmod__(self, other):
        return Mod(other, self)

    def __floordiv__(self, other):
        return FloorDiv(self, other)

    def __rfloordiv__(self, other):
        return FloorDiv(other, self)

    def __divmod__(self, other):
        res1 = self // other
        res2 = self % other
        return res1, res2

    def __rdivmod__(self, other):
        res1 = other // self
        res2 = other % self
        return res1, res2

    def sum(self, skipna=True, numeric_only=False, split_every=False, axis=0):
        return Sum(self, skipna, numeric_only, split_every, axis)

    def prod(self, skipna=True, numeric_only=False, split_every=False, axis=0):
        return Prod(self, skipna, numeric_only, split_every, axis)

    def var(self, axis=0, skipna=True, ddof=1, numeric_only=False, split_every=False):
        if axis == 0:
            return Var(self, skipna, ddof, numeric_only, split_every)
        elif axis == 1:
            return VarColumns(self, skipna, ddof, numeric_only)
        else:
            raise ValueError(f"axis={axis} not supported. Please specify 0 or 1")

    def std(self, axis=0, skipna=True, ddof=1, numeric_only=False, split_every=False):
        return Sqrt(self.var(axis, skipna, ddof, numeric_only, split_every=split_every))

    def mean(self, skipna=True, numeric_only=False, split_every=False, axis=0):
        return Mean(
            self,
            skipna=skipna,
            numeric_only=numeric_only,
            split_every=split_every,
            axis=axis,
        )

    def max(self, skipna=True, numeric_only=False, split_every=False, axis=0):
        return Max(self, skipna, numeric_only, split_every, axis)

    def any(self, skipna=True, split_every=False):
        return Any(self, skipna=skipna, split_every=split_every)

    def all(self, skipna=True, split_every=False):
        return All(self, skipna=skipna, split_every=split_every)

    def idxmin(self, skipna=True, numeric_only=False, split_every=False):
        return IdxMin(
            self, skipna=skipna, numeric_only=numeric_only, split_every=split_every
        )

    def idxmax(self, skipna=True, numeric_only=False, split_every=False):
        return IdxMax(
            self, skipna=skipna, numeric_only=numeric_only, split_every=split_every
        )

    def mode(self, dropna=True, split_every=False):
        return Mode(self, dropna=dropna, split_every=split_every)

    def min(self, skipna=True, numeric_only=False, split_every=False, axis=0):
        return Min(self, skipna, numeric_only, split_every=split_every, axis=axis)

    def count(self, numeric_only=False, split_every=False):
        return Count(self, numeric_only, split_every)

    def cumsum(self, skipna=True):
        from dask.dataframe.dask_expr._cumulative import CumSum

        return CumSum(self, skipna=skipna)

    def cumprod(self, skipna=True):
        from dask.dataframe.dask_expr._cumulative import CumProd

        return CumProd(self, skipna=skipna)

    def cummax(self, skipna=True):
        from dask.dataframe.dask_expr._cumulative import CumMax

        return CumMax(self, skipna=skipna)

    def cummin(self, skipna=True):
        from dask.dataframe.dask_expr._cumulative import CumMin

        return CumMin(self, skipna=skipna)

    def abs(self):
        return Abs(self)

    def astype(self, dtypes):
        return AsType(self, dtypes)

    def clip(self, lower=None, upper=None, axis=None):
        return Clip(self, lower=lower, upper=upper, axis=axis)

    def combine_first(self, other):
        if are_co_aligned(self, other):
            return CombineFirst(self, other=other)
        else:
            return CombineFirstAlign(self, other)

    def to_timestamp(self, freq=None, how="start"):
        return ToTimestamp(self, freq=freq, how=how)

    def isna(self):
        return IsNa(self)

    def isnull(self):
        # These are the same anyway
        return IsNa(self)

    def round(self, decimals=0):
        return Round(self, decimals=decimals)

    def where(self, cond, other=np.nan):
        if not are_co_aligned(self, *[c for c in [cond, other] if isinstance(c, Expr)]):
            return WhereAlign(self, cond=cond, other=other)
        return Where(self, cond=cond, other=other)

    def mask(self, cond, other=np.nan):
        if not are_co_aligned(self, *[c for c in [cond, other] if isinstance(c, Expr)]):
            return MaskAlign(self, cond=cond, other=other)
        return Mask(self, cond=cond, other=other)

    def apply(self, function, *args, meta=None, **kwargs):
        return Apply(self, function, args, meta, kwargs)

    def replace(self, to_replace=None, value=no_default, regex=False):
        return Replace(self, to_replace=to_replace, value=value, regex=regex)

    def fillna(self, value=None):
        if isinstance(value, Expr) and not are_co_aligned(self, value):
            return FillnaAlign(self, value=value)
        return Fillna(self, value=value)

    def rename_axis(
        self, mapper=no_default, index=no_default, columns=no_default, axis=0
    ):
        return RenameAxis(self, mapper=mapper, index=index, columns=columns, axis=axis)

    def align(self, other, join="outer", axis=None, fill_value=None):
        from dask.dataframe.dask_expr._collection import new_collection

        if not are_co_aligned(self, other):
            aligned = AlignAlignPartitions(self, other, join, axis, fill_value)
        else:
            aligned = _Align(self, other, join, axis=axis, fill_value=fill_value)

        return new_collection(AlignGetitem(aligned, position=0)), new_collection(
            AlignGetitem(aligned, position=1)
        )

    def nunique_approx(self, split_every=None):
        return NuniqueApprox(self, b=16, split_every=split_every)

    def memory_usage_per_partition(self, index=True, deep=False):
        return MemoryUsagePerPartition(self, index, deep)

    @functools.cached_property
    def divisions(self):
        return tuple(self._divisions())

    def _divisions(self):
        raise NotImplementedError()

    @property
    def known_divisions(self):
        """Whether divisions are already known"""
        return len(self.divisions) > 0 and self.divisions[0] is not None

    @property
    def npartitions(self):
        if "npartitions" in self._parameters:
            return self.operand("npartitions")
        else:
            return len(self.divisions) - 1

    @property
    def columns(self) -> list:
        try:
            return list(self._meta.columns)
        except AttributeError:
            if self.ndim == 1:
                return [self.name]
            return []
        except Exception:
            raise

    @functools.cached_property
    def unique_partition_mapping_columns_from_shuffle(self) -> set:
        """Preserves the columns defining the partition mapping from shuffles.

        This property specifies if a column or a set of columns have a unique
        partition mapping that was defined by a shuffle operation. The mapping
        is created by hashing the values and the separating them onto partitions.
        It is important that this property is only propagated if the values
        in those columns did not change in this expression. The property is
        only populated if the mapping was created by the ``partitioning_index``
        function.

        Simply knowing that every value is in only one partition is not a
        satisfying condition, because we also use this property on merge
        operations, where we need these values to be in matching partitions.

        This is also the reason why set_index or sort_values can't set the
        property, they fulfill a weaker condition than what this property enforces.

        Normally, this set contains one tuple of either one or multiple columns.
        It can contain 2, when the operation shuffles multiple columns of the
        result, i.e. a merge operation and the left and right join columns.


        Returns
        -------
            A set of column groups that have a unique partition mapping as
            defined by a shuffle.
        """
        return set()

    @property
    def _projection_columns(self):
        return self.columns

    @property
    def name(self):
        return self._meta.name

    @property
    def dtypes(self):
        return self._meta.dtypes

    def _filter_simplification(self, parent, predicate=None):
        if predicate is None:
            predicate = parent.predicate.substitute(self, self.frame)
        if are_co_aligned(self.frame, predicate):
            # Only do this if we are aligned
            return type(self)(self.frame[predicate], *self.operands[1:])

    def fuse(self):
        return optimize_blockwise_fusion(self)

    def finalize_compute(self):
        return FinalizeComputeDF(self)


class FinalizeComputeDF(FinalizeCompute, Expr):
    _parameters = ["frame"]

    def _simplify_down(self):
        from dask.dataframe.dask_expr._repartition import Repartition

        return Repartition(self.frame, 1)


class Literal(Expr):
    """Represent a literal (known) value as an `Expr`"""

    _parameters = ["value"]

    def _divisions(self):
        return (None, None)

    @functools.cached_property
    def _meta(self):
        return make_meta(self.value)

    def _task(self, name: Key, index: int) -> Task:
        assert index == 0
        return DataNode(name, self.value)  # type: ignore


class Blockwise(Expr):
    """Super-class for block-wise operations

    This is fairly generic, and includes definitions for `_meta`, `divisions`,
    `_layer` that are often (but not always) correct.  Mostly this helps us
    avoid duplication in the future.

    Note that `Fused` expressions rely on every `Blockwise`
    expression defining a proper `_task` method.
    """

    operation: Callable | None = None
    _keyword_only: list[str] = []
    _projection_passthrough = False
    _preserves_partitioning_information = False

    @functools.cached_property
    def _meta(self):
        args = [op._meta if isinstance(op, Expr) else op for op in self._args]
        return self.operation(*args, **self._kwargs)

    @functools.cached_property
    def _kwargs(self) -> dict:
        if self._keyword_only:
            return {
                p: self.operand(p)
                for p in self._parameters
                if p in self._keyword_only and self.operand(p) is not no_default
            }
        return {}

    @functools.cached_property
    def _args(self) -> list:
        if self._keyword_only:
            args = [
                self.operand(p) for p in self._parameters if p not in self._keyword_only
            ] + self.operands[len(self._parameters) :]
            return args
        return self.operands

    def _broadcast_dep(self, dep: Expr):
        # Checks if a dependency should be broadcasted to
        # all partitions of this `Blockwise` operation
        return dep.npartitions == 1 and dep.ndim < self.ndim

    def _divisions(self):
        # This is an issue.  In normal Dask we re-divide everything in a step
        # which combines divisions and graph.
        # We either have to create a new Align layer (ok) or combine divisions
        # and graph into a single operation.
        dependencies = self.dependencies()
        for arg in dependencies:
            if not self._broadcast_dep(arg):
                assert arg.divisions == dependencies[0].divisions
        return dependencies[0].divisions

    @functools.cached_property
    def _name(self):
        if self.operation:
            head = funcname(self.operation)
        else:
            head = funcname(type(self)).lower()
        return head + "-" + self.deterministic_token

    def _blockwise_arg(self, arg, i):
        """Return a Blockwise-task argument"""
        if isinstance(arg, Expr):
            # Make key for Expr-based argument
            if self._broadcast_dep(arg):
                return TaskRef((arg._name, 0))
            else:
                return TaskRef((arg._name, i))

        else:
            return arg

    def _task(self, name: Key, index: int) -> Task:
        """Produce the task for a specific partition

        Parameters
        ----------
        index:
            Partition index for this task.

        Returns
        -------
        task: tuple
        """
        args = [self._blockwise_arg(op, index) for op in self._args]
        if self._kwargs:
            return Task(name, self.operation, *args, **self._kwargs)  # type: ignore
        else:
            return Task(name, self.operation, *args)  # type: ignore

    def _simplify_up(self, parent, dependents):
        if self._projection_passthrough and isinstance(parent, Projection):
            return plain_column_projection(self, parent, dependents)

    @functools.cached_property
    def unique_partition_mapping_columns_from_shuffle(self):
        if self._preserves_partitioning_information:
            return self.frame.unique_partition_mapping_columns_from_shuffle
        return set()


class MapPartitions(Blockwise):
    _parameters = [
        "frame",
        "func",
        "meta",
        "enforce_metadata",
        "transform_divisions",
        "clear_divisions",
        "align_dataframes",
        "parent_meta",
        "required_columns",
        "token",
        "kwargs",
        "nargs",
    ]
    _defaults: dict = {
        "kwargs": None,
        "align_dataframes": True,
        "parent_meta": None,
        "required_columns": None,
        "token": None,
        "nargs": 0,
    }

    @functools.cached_property
    def token(self):
        if "token" in self._parameters:
            return self.operand("token")
        return None

    def __str__(self):
        return f"MapPartitions({funcname(self.func)})"

    @functools.cached_property
    def _name(self):
        if self.token is not None:
            head = self.token
        else:
            head = funcname(self.func).lower()
        return head + "-" + self.deterministic_token

    def _broadcast_dep(self, dep: Expr):
        # Always broadcast single-partition dependencies in MapPartitions
        return dep.npartitions == 1

    @functools.cached_property
    def args(self):
        return [self.frame] + self.operands[
            len(self._parameters) : len(self._parameters) + self.nargs
        ]

    @functools.cached_property
    def _meta(self):
        meta = self.operand("meta")
        return _get_meta_map_partitions(
            self.args,
            [
                e
                for e in self.args
                if isinstance(e, Expr) and not isinstance(e, _DelayedExpr)
            ],
            self.func,
            self.kwargs,
            meta,
            self.parent_meta,
        )

    def _divisions(self):
        # Unknown divisions
        dfs = [arg for arg in self.args if isinstance(arg, Expr)]

        if self.clear_divisions:
            max_partitions = max(df.npartitions for df in dfs)
            return (None,) * (max_partitions + 1)

        # (Possibly) known divisions
        return _get_divisions_map_partitions(
            self.align_dataframes,
            self.transform_divisions,
            dfs,
            self.func,
            self.args,
            self.kwargs,
        )

    @functools.cached_property
    def _has_partition_info(self):
        return has_keyword(self.func, "partition_info")

    def _task(self, name: Key, index: int) -> Task:
        args = [self._blockwise_arg(op, index) for op in self.args]
        kwargs = dict(self.kwargs if self.kwargs is not None else {})
        if self._has_partition_info:
            kwargs["partition_info"] = {
                "number": index,
                "division": self.divisions[index],
            }

        if self.enforce_metadata:
            kwargs.update(
                {
                    "_func": self.func,
                    "_meta": self._meta,
                }
            )
            return Task(name, apply_and_enforce, *args, **kwargs)
        else:
            return Task(
                name,
                self.func,
                *args,
                **kwargs,
            )

    @staticmethod
    def projected_operation(mapped_func, post_projection, *args, **kwargs):
        # Apply a mapped function and then project columns.
        # Used by `_simplify_up` to apply column projection.
        return mapped_func(*args, **kwargs)[post_projection]

    def _simplify_up(self, parent, dependents):
        if isinstance(parent, Projection) and self.required_columns is not None:
            if missing := set(self.required_columns) - set(self.frame.columns):
                raise KeyError(
                    f"Some elements of `required_columns` are missing: {missing}"
                )

            columns = determine_column_projection(
                self, parent, dependents, additional_columns=self.required_columns
            )
            columns = [col for col in self.frame.columns if col in columns]

            if columns == self.frame.columns:
                # Don't add unnecessary Projections
                return

            return type(parent)(
                type(self)(
                    self.frame[columns],
                    partial(self.projected_operation, self.func, parent.columns),
                    self.meta[parent.columns],
                    *self.operands[3:],
                ),
                *parent.operands[1:],
            )


def _get_meta_ufunc(dfs, args, func):
    dasks = [arg for arg in args if isinstance(arg, (Expr, Array))]

    if len(dfs) >= 2 and not all(hasattr(d, "npartitions") for d in dasks):
        # should not occur in current funcs
        msg = "elemwise with 2 or more DataFrames and Scalar is not supported"
        raise NotImplementedError(msg)
    # For broadcastable series, use no rows.
    parts = [
        (
            d._meta
            if d.ndim == 0
            else (
                np.empty((), dtype=d.dtype)
                if isinstance(d, Array)
                else meta_nonempty(d._meta)
            )
        )
        for d in dasks
    ]

    other = [
        (i, arg) for i, arg in enumerate(args) if not isinstance(arg, (Expr, Array))
    ]

    return partial_by_order(*parts, function=func, other=other)


class UFuncElemwise(MapPartitions):
    _parameters = [
        "frame",
        "func",
        "meta",
        "transform_divisions",
        "kwargs",
    ]
    enforce_metadata = False

    def __str__(self):
        return f"UFunc({funcname(self.func)})"

    @functools.cached_property
    def args(self):
        return self.operands[len(self._parameters) :]

    @functools.cached_property
    def _dfs(self):
        return [df for df in self.args if isinstance(df, Expr) and df.ndim > 0]

    @functools.cached_property
    def _meta(self):
        if self.operand("meta") is not no_default:
            meta = self.operand("meta")
        else:
            meta = _get_meta_ufunc(self._dfs, self.args, self.func)
        return make_meta(meta)

    def _divisions(self):
        if (
            self.transform_divisions
            and isinstance(self._dfs[0], Index)
            and len(self._dfs) == 1
        ):
            try:
                divisions = self.func(
                    *[
                        pd.Index(arg.divisions) if arg is self._dfs[0] else arg
                        for arg in self.args
                    ],
                    **self.kwargs,
                )
                if isinstance(divisions, pd.Index):
                    divisions = methods.tolist(divisions)
            except Exception:
                pass
            else:
                if not valid_divisions(divisions):
                    divisions = [None] * (self._dfs[0].npartitions + 1)
                return divisions

        return self._dfs[0].divisions


class MapOverlapAlign(Expr):
    _parameters = [
        "frame",
        "func",
        "before",
        "after",
        "meta",
        "enforce_metadata",
        "transform_divisions",
        "clear_divisions",
        "align_dataframes",
        "token",
        "kwargs",
    ]
    _defaults = {
        "meta": None,
        "enfore_metadata": True,
        "transform_divisions": True,
        "kwargs": None,
        "clear_divisions": False,
        "align_dataframes": False,
        "token": None,
    }

    @functools.cached_property
    def _meta(self):
        meta = self.operand("meta")
        args = [self.frame._meta] + [
            arg._meta if isinstance(arg, Expr) else arg
            for arg in self.operands[len(self._parameters) :]
        ]
        return _get_meta_map_partitions(
            args,
            [self.dependencies()[0]],
            self.func,
            self.kwargs,
            meta,
            self.kwargs.pop("parent_meta", None),
        )

    def _divisions(self):
        args = [self.frame] + self.operands[len(self._parameters) :]
        return calc_divisions_for_align(*args, allow_shuffle=False)

    def _lower(self):
        args = [self.frame] + self.operands[len(self._parameters) :]
        args = maybe_align_partitions(*args, divisions=self._divisions())
        return MapOverlap(
            args[0],
            self.func,
            self.before,
            self.after,
            self._meta,
            self.enforce_metadata,
            self.transform_divisions,
            self.clear_divisions,
            self.align_dataframes,
            self.token,
            self.kwargs,
            *args[1:],
        )


class MapOverlap(MapPartitions):
    _parameters = [
        "frame",
        "func",
        "before",
        "after",
        "meta",
        "enforce_metadata",
        "transform_divisions",
        "clear_divisions",
        "align_dataframes",
        "token",
        "kwargs",
    ]
    _defaults: dict = {
        "meta": None,
        "enfore_metadata": True,
        "transform_divisions": True,
        "kwargs": None,
        "clear_divisions": False,
        "align_dataframes": False,
        "token": None,
    }

    @functools.cached_property
    def _kwargs(self) -> dict:
        kwargs = self.kwargs
        if kwargs is None:
            kwargs = {}
        return kwargs

    @property
    def args(self):
        return (
            [self.frame]
            + [self.func, self.before, self.after]
            + self.operands[len(self._parameters) :]
        )

    @functools.cached_property
    def _meta(self):
        meta = self.operand("meta")
        args = [self.frame._meta] + [
            arg._meta if isinstance(arg, Expr) else arg
            for arg in self.operands[len(self._parameters) :]
        ]
        return _get_meta_map_partitions(
            args,
            [self.dependencies()[0]],
            self.func,
            self.kwargs,
            meta,
            self.kwargs.pop("parent_meta", None),
        )

    @functools.cached_property
    def before(self):
        before = self.operand("before")
        if isinstance(before, str):
            return pd.to_timedelta(before)
        return before

    @functools.cached_property
    def after(self):
        after = self.operand("after")
        if isinstance(after, str):
            return pd.to_timedelta(after)
        return after

    def _lower(self):
        overlapped = CreateOverlappingPartitions(self.frame, self.before, self.after)

        return MapPartitions(
            overlapped,
            _overlap_chunk,
            self._meta,
            self.enforce_metadata,
            self.transform_divisions,
            self.clear_divisions,
            self.align_dataframes,
            None,
            None,
            self.token,
            self._kwargs,
            len(self.args[1:]),
            *self.args[1:],
        )


class CreateOverlappingPartitions(Expr):
    _parameters = ["frame", "before", "after"]

    @functools.cached_property
    def _meta(self):
        return self.frame._meta

    def _divisions(self):
        # Keep divisions alive, MapPartitions will handle the actual division logic
        return self.frame.divisions

    def _layer(self) -> dict:
        dsk, prevs, nexts = {}, [], []  # type: ignore

        name_prepend = "overlap-prepend" + self.frame._name
        if self.before:
            prevs.append(None)
            if isinstance(self.before, numbers.Integral):
                before = self.before
                for i in range(self.frame.npartitions - 1):
                    dsk[(name_prepend, i)] = (M.tail, (self.frame._name, i), before)
                    prevs.append((name_prepend, i))
            elif isinstance(self.before, datetime.timedelta):
                # Assumes monotonic (increasing?) index
                divs = pd.Series(self.frame.divisions)
                deltas = divs.diff().iloc[1:-1]

                # In the first case window-size is larger than at least one partition, thus it is
                # necessary to calculate how many partitions must be used for each rolling task.
                # Otherwise, these calculations can be skipped (faster)

                if (self.before > deltas).any():
                    pt_z = divs[0]
                    for i in range(self.frame.npartitions - 1):
                        # Select all indexes of relevant partitions between the current partition and
                        # the partition with the highest division outside the rolling window (before)
                        pt_i = divs[i + 1]

                        # lower-bound the search to the first division
                        lb = max(pt_i - self.before, pt_z)

                        first, j = divs[i], i
                        while first > lb and j > 0:
                            first = first - deltas[j]
                            j = j - 1

                        dsk[(name_prepend, i)] = (  # type: ignore
                            _tail_timedelta,
                            (self.frame._name, i + 1),
                            [(self.frame._name, k) for k in range(j, i + 1)],
                            self.before,
                        )
                        prevs.append((name_prepend, i))
                else:
                    for i in range(self.frame.npartitions - 1):
                        dsk[(name_prepend, i)] = (  # type: ignore
                            _tail_timedelta,
                            (self.frame._name, i + 1),
                            [(self.frame._name, i)],
                            self.before,
                        )
                        prevs.append((name_prepend, i))
        else:
            prevs.extend([None] * self.frame.npartitions)  # type: ignore

        name_append = "overlap-append" + self.frame._name
        if self.after:
            if isinstance(self.after, numbers.Integral):
                after = self.after
                for i in range(1, self.frame.npartitions):
                    dsk[(name_append, i)] = (M.head, (self.frame._name, i), after)
                    nexts.append((name_append, i))
            else:
                # We don't want to look at the divisions, so take twice the step and
                # validate later.
                after = 2 * self.after
                for i in range(1, self.frame.npartitions):
                    dsk[(name_append, i)] = (  # type: ignore
                        _head_timedelta,
                        (self.frame._name, i - 1),
                        (self.frame._name, i),
                        after,
                    )
                    nexts.append((name_append, i))

            nexts.append(None)  # type: ignore

        else:
            nexts.extend([None] * self.frame.npartitions)  # type: ignore

        for i, (prev, next) in enumerate(zip(prevs, nexts)):
            dsk[(self._name, i)] = (  # type: ignore
                _combined_parts,
                prev,
                (self.frame._name, i),
                next,
                self.before,
                self.after,
            )
        return dsk


def _tail_timedelta(current, prev_, before):
    selected = methods.concat(
        [prev[prev.index > (current.index.min() - before)] for prev in prev_]
    )
    return selected


def _overlap_chunk(df, func, before, after, *args, **kwargs):
    return overlap_chunk(func, before, after, df, *args, **kwargs)


def _combined_parts(prev_part, current_part, next_part, before, after):
    msg = (
        "Partition size is less than overlapping "
        "window size. Try using ``df.repartition`` "
        "to increase the partition size."
    )

    if prev_part is not None:
        if isinstance(before, numbers.Integral):
            if prev_part.shape[0] != before:
                raise NotImplementedError(msg)
        else:
            prev_part_input = prev_part
            prev_part = _tail_timedelta(current_part, [prev_part], before)
            if (
                len(prev_part_input) == len(prev_part)
                and len(prev_part_input) > 0
                and not isinstance(before, datetime.timedelta)
            ):
                raise NotImplementedError(msg)

    if next_part is not None:
        if isinstance(after, numbers.Integral):
            if next_part.shape[0] != after:
                raise NotImplementedError(msg)
        else:
            next_part_input = next_part
            next_part = _head_timedelta(current_part, next_part, after)
            if len(next_part_input) == len(next_part) and len(next_part_input) > 0:
                raise NotImplementedError(msg)

    parts = [p for p in (prev_part, current_part, next_part) if p is not None]
    combined = methods.concat(parts)

    return CombinedOutput(
        (
            combined,
            len(prev_part) if prev_part is not None and len(prev_part) > 0 else None,
            len(next_part) if next_part is not None and len(next_part) > 0 else None,
        )
    )


class _Align(Blockwise):
    _parameters = ["frame", "other", "join", "axis", "fill_value"]
    _defaults = {"join": "outer", "fill_value": None, "axis": None}
    _keyword_only = ["join", "fill_value", "axis"]
    operation = M.align
    _preserves_partitioning_information = True

    def _divisions(self):
        # Aligning, so take first frames divisions
        return self.frame._divisions()


class AlignGetitem(Blockwise):
    _parameters = ["frame", "position"]
    operation = operator.getitem
    _preserves_partitioning_information = True

    @functools.cached_property
    def _meta(self):
        return self.frame._meta[self.position]

    def _divisions(self):
        return self.frame.divisions


class ScalarToSeries(Blockwise):
    _parameters = ["frame", "index"]
    _defaults = {"index": 0}

    @staticmethod
    def operation(value, index=0):
        return pd.Series(value, index=[index])


class DropnaSeries(Blockwise):
    _parameters = ["frame"]
    operation = M.dropna
    _preserves_partitioning_information = True


class DropnaFrame(Blockwise):
    _parameters = ["frame", "how", "subset", "thresh"]
    _defaults = {"how": no_default, "subset": None, "thresh": no_default}
    _keyword_only = ["how", "subset", "thresh"]
    operation = M.dropna
    _preserves_partitioning_information = True

    def _simplify_up(self, parent, dependents):
        if isinstance(parent, Projection) and self.subset is not None:
            columns = determine_column_projection(
                self, parent, dependents, additional_columns=self.subset
            )
            columns = [col for col in self.frame.columns if col in columns]

            if columns == self.frame.columns:
                # Don't add unnecessary Projections
                return

            return type(parent)(
                type(self)(self.frame[columns], *self.operands[1:]),
                *parent.operands[1:],
            )


class CombineFirst(Blockwise):
    _parameters = ["frame", "other"]
    operation = M.combine_first

    @functools.cached_property
    def _meta(self):
        return make_meta(
            self.operation(
                meta_nonempty(self.frame._meta),
                meta_nonempty(self.other._meta),
            ),
        )

    def _simplify_up(self, parent, dependents):
        if isinstance(parent, Projection):
            columns = determine_column_projection(self, parent, dependents)
            frame_columns = [col for col in self.frame.columns if col in columns]
            other_columns = [col for col in self.other.columns if col in columns]
            if (
                self.frame.columns == frame_columns
                and self.other.columns == other_columns
            ):
                return

            return type(parent)(
                type(self)(self.frame[frame_columns], self.other[other_columns]),
                *parent.operands[1:],
            )


class Sample(Blockwise):
    _parameters = ["frame", "state_data", "frac", "replace"]
    operation = staticmethod(methods.sample)

    @functools.cached_property
    def _meta(self):
        args = [self.operands[0]._meta] + [self.operands[1][0]] + self.operands[2:]
        return self.operation(*args)

    def _task(self, name: Key, index: int) -> Task:
        args = [self._blockwise_arg(self.frame, index)] + [
            self.state_data[index],
            self.frac,
            self.operand("replace"),
        ]
        return Task(name, self.operation, *args)


class Query(Blockwise):
    _parameters = ["frame", "_expr", "expr_kwargs"]
    _defaults: dict[str, Any] = {"expr_kwargs": {}}  # type: ignore
    _keyword_only = ["expr_kwargs"]
    operation = M.query

    @functools.cached_property
    def _kwargs(self) -> dict:
        return {**self.expr_kwargs}


class MemoryUsagePerPartition(Blockwise):
    _parameters = ["frame", "index", "deep"]
    _defaults = {"index": True, "deep": False}

    @staticmethod
    def operation(*args, **kwargs):
        if is_series_like(args[0]):
            return args[0]._constructor([total_mem_usage(*args, **kwargs)])
        return args[0]._constructor_sliced([total_mem_usage(*args, **kwargs)])

    def _divisions(self):
        return (None,) * (self.frame.npartitions + 1)


class DropDuplicatesBlockwise(Blockwise):
    _parameters = ["frame"]
    operation = M.drop_duplicates
    _preserves_partitioning_information = True


class Elemwise(Blockwise):
    """
    This doesn't really do anything, but we anticipate that future
    optimizations, like `len` will care about which operations preserve length
    """

    _is_length_preserving = True

    def _simplify_up(self, parent, dependents):
        if isinstance(parent, Filter) and self._filter_passthrough_available(
            parent, dependents
        ):
            predicate = None
            if self.frame.ndim == 1 and self.ndim == 2:
                name = self.frame._meta.name
                # Avoid Projection since we are already a Series
                subs = Projection(self, name)
                predicate = parent.predicate.substitute(subs, self.frame)
            return self._filter_simplification(parent, predicate)
        return super()._simplify_up(parent, dependents)


class RenameFrame(Elemwise):
    _parameters = ["frame", "columns"]

    @functools.cached_property
    def unique_partition_mapping_columns_from_shuffle(self):
        result = set()
        columns = self.operand("columns")
        for elem in self.frame.unique_partition_mapping_columns_from_shuffle:
            if isinstance(elem, tuple):
                subset = self.frame._meta[list(elem)].rename(columns=columns)
                result.add(tuple(subset.columns))
            else:
                # scalar
                subset = self.frame._meta[[elem]]
                result.add(subset.columns[0])
        return result

    @staticmethod
    def operation(df, columns):
        return df.rename(columns=columns)

    def _simplify_up(self, parent, dependents):
        if isinstance(parent, Projection) and isinstance(
            self.operand("columns"), Mapping
        ):
            reverse_mapping = {val: key for key, val in self.operand("columns").items()}

            columns = determine_column_projection(self, parent, dependents)
            columns = _convert_to_list(columns)
            frame_columns = set(self.frame.columns)
            columns = [
                (
                    reverse_mapping[col]
                    if col in reverse_mapping and reverse_mapping[col] in frame_columns
                    else col
                )
                for col in columns
            ]
            columns = [col for col in self.frame.columns if col in columns]
            if columns == self.frame.columns:
                return

            return type(parent)(
                type(self)(self.frame[columns], *self.operands[1:]),
                *parent.operands[1:],
            )


class ColumnsSetter(RenameFrame):
    _preserves_partitioning_information = True

    @staticmethod
    def operation(df, columns):
        return _rename(columns, df)


class _DeepCopy(Elemwise):
    _parameters = ["frame"]
    _projection_passthrough = True
    _filter_passthrough = True
    _preserves_partitioning_information = True

    @staticmethod
    def operation(df):
        return df.copy(deep=True)


class ToBackend(Elemwise):
    _parameters = ["frame", "options"]
    _projection_passthrough = True
    _filter_passthrough = True
    _preserves_partitioning_information = True


class RenameSeries(Elemwise):
    _parameters = ["frame", "index", "sorted_index"]
    _defaults = {"sorted_index": False}
    _filter_passthrough = True
    _preserves_partitioning_information = True

    @functools.cached_property
    def _meta(self):
        args = [
            meta_nonempty(op._meta) if isinstance(op, Expr) else op for op in self._args
        ]
        return make_meta(self.operation(*args, **self._kwargs))

    @staticmethod
    def operation(df, index, sorted_index):
        if is_series_like(df):
            return df.rename(index=index)
        return df.rename(name=index)

    def _divisions(self):
        index = self.operand("index")
        if is_scalar(index) and not isinstance(index, Callable):
            return self.frame.divisions
        elif self.sorted_index and self.frame.known_divisions:
            old = pd.Series(1, index=self.frame.divisions)
            new_divisions = old.rename(index).index
            if not new_divisions.is_monotonic_increasing:
                raise ValueError(
                    "The renamer creates an Index with non-monotonic divisions. "
                    "This is not allowed. Please set sorted_index=False."
                )
            return tuple(new_divisions.tolist())
        else:
            return (None,) * (self.frame.npartitions + 1)


class Fillna(Elemwise):
    _projection_passthrough = True
    _parameters = ["frame", "value"]
    _defaults = {"value": None}
    operation = M.fillna


class Replace(Elemwise):
    _projection_passthrough = True
    _parameters = ["frame", "to_replace", "value", "regex"]
    _defaults = {"to_replace": None, "value": no_default, "regex": False}
    _keyword_only = ["value", "regex"]
    operation = M.replace


class Isin(Elemwise):
    _projection_passthrough = True
    _parameters = ["frame", "values"]
    operation = M.isin

    @functools.cached_property
    def _meta(self):
        return make_meta(
            meta_nonempty(self.frame._meta).isin(
                meta_nonempty(self.frame._meta).iloc[[0]]
            )
        )

    def _broadcast_dep(self, dep: Expr):
        return dep.npartitions == 1


class Clip(Elemwise):
    _projection_passthrough = True
    _parameters = ["frame", "lower", "upper", "axis"]
    _defaults = {"lower": None, "upper": None, "axis": None}
    _keyword_only = ["axis"]
    operation = M.clip

    def _simplify_up(self, parent, dependents):
        if isinstance(parent, Projection):
            return plain_column_projection(self, parent, dependents)


class ArrowStringConversion(Elemwise):
    _projection_passthrough = True
    _filter_passthrough = True
    _parameters = ["frame"]
    operation = staticmethod(to_pyarrow_string)
    _preserves_partitioning_information = True


class Between(Elemwise):
    _parameters = ["frame", "left", "right", "inclusive"]
    _defaults = {"inclusive": "both"}
    operation = M.between


class ToTimestamp(Elemwise):
    _projection_passthrough = True
    _parameters = ["frame", "freq", "how"]
    _defaults = {"freq": None, "how": "start"}
    operation = M.to_timestamp
    _filter_passthrough = True
    _preserves_partitioning_information = True

    def _divisions(self):
        return tuple(
            pd.Index(self.frame.divisions).to_timestamp(freq=self.freq, how=self.how)
        )


class CombineSeries(Elemwise):
    _parameters = ["frame", "other", "func", "fill_value"]
    _defaults: dict = {"fill_value": None}
    operation = M.combine

    @functools.cached_property
    def _meta(self):
        return make_meta(
            meta_nonempty(self.frame._meta).combine(
                meta_nonempty(self.other._meta), func=self.func
            )
        )


class CombineFrame(CombineSeries):
    _parameters = CombineSeries._parameters + ["overwrite"]
    _defaults = {"fill_value": None, "overwrite": True}


class ToNumeric(Elemwise):
    _parameters = ["frame", "errors", "downcast", "meta"]
    _defaults = {"errors": "raise", "downcast": None, "meta": None}
    _keyword_only = ["meta"]
    operation = staticmethod(pd.to_numeric)

    @functools.cached_property
    def _kwargs(self):
        kwargs = super()._kwargs
        kwargs.pop("meta", None)
        return kwargs

    @functools.cached_property
    def _meta(self):
        if self.operand("meta") is not None:
            return self.operand("meta")
        return super()._meta


class ToDatetime(Elemwise):
    _parameters = ["frame", "kwargs", "meta"]
    _defaults = {"kwargs": None}
    _keyword_only = ["kwargs", "meta"]

    @functools.cached_property
    def _meta(self):
        return self.operand("meta")

    @staticmethod
    def operation(*args, **kwargs):
        return get_meta_library(args[0]).to_datetime(*args, **kwargs)

    @functools.cached_property
    def _kwargs(self):
        if (kwargs := self.operand("kwargs")) is None:
            return {}
        return kwargs


class ToTimedelta(Elemwise):
    _parameters = ["frame", "unit", "errors"]
    _defaults = {"unit": None, "errors": "raise"}
    operation = staticmethod(pd.to_timedelta)


class AsType(Elemwise):
    """A good example of writing a trivial blockwise operation"""

    _parameters = ["frame", "dtypes"]
    operation = M.astype
    _filter_passthrough = True

    @functools.cached_property
    def _meta(self):
        def _cat_dtype_without_categories(dtype):
            return (
                isinstance(pd.api.types.pandas_dtype(dtype), pd.CategoricalDtype)
                and getattr(dtype, "categories", None) is None
            )

        meta = super()._meta
        dtypes = self.operand("dtypes")
        if hasattr(dtypes, "items"):
            set_unknown = [
                k for k, v in dtypes.items() if _cat_dtype_without_categories(v)
            ]
            meta = clear_known_categories(meta, cols=set_unknown)

        elif _cat_dtype_without_categories(dtypes):
            meta = clear_known_categories(meta)
        return meta

    def _simplify_up(self, parent, dependents):
        if isinstance(parent, Filter) and self._filter_passthrough_available(
            parent, dependents
        ):
            return self._filter_simplification(parent)
        if isinstance(parent, Projection):
            dtypes = self.operand("dtypes")
            columns = determine_column_projection(self, parent, dependents)
            if isinstance(dtypes, dict):
                dtypes = {key: val for key, val in dtypes.items() if key in columns}
                if not dtypes:
                    return type(parent)(self.frame, *parent.operands[1:])
            if isinstance(columns, list):
                columns = [col for col in self.frame.columns if col in columns]
            if self.frame.columns == columns:
                return
            result = type(self)(self.frame[columns], dtypes)
            if not isinstance(columns, list):
                return result
            return type(parent)(result, *parent.operands[1:])


class IsNa(Elemwise):
    _projection_passthrough = True
    _parameters = ["frame"]
    operation = M.isna


class Mask(Elemwise):
    _projection_passthrough = True
    _parameters = ["frame", "cond", "other"]
    _defaults = {"other": np.nan}
    operation = M.mask


class Round(Elemwise):
    _projection_passthrough = True
    _parameters = ["frame", "decimals"]
    operation = M.round


class Where(Elemwise):
    _projection_passthrough = True
    _parameters = ["frame", "cond", "other"]
    _defaults = {"other": np.nan}
    operation = M.where


def _check_divisions(df, i, division_min, division_max, last):
    if not len(df):
        return df
    if is_index_like(df):
        index = df
    else:
        try:
            index = df.index.get_level_values(0)
        except AttributeError:
            index = df.index
    # Check divisions
    real_min = index.min()
    real_max = index.max()
    # Upper division of the last partition is often set to
    # the max value. For all other partitions, the upper
    # division should be greater than the maximum value.
    valid_min = valid_max = True
    if not pd.isna(division_min):
        valid_min = real_min >= division_min
    if not pd.isna(division_max):
        valid_max = (real_max <= division_max) if last else (real_max < division_max)
    if not (valid_min and valid_max):
        raise RuntimeError(
            f"`enforce_runtime_divisions` failed for partition {i}."
            f" Expected a range of [{division_min}, {division_max}), "
            f" but the real range was [{real_min}, {real_max}]."
        )
    return df


class EnforceRuntimeDivisions(Blockwise):
    _parameters = ["frame"]
    operation = staticmethod(_check_divisions)
    _preserves_partitioning_information = True

    @functools.cached_property
    def _meta(self):
        return self.frame._meta

    def _task(self, name: Key, index: int) -> Task:
        args = [self._blockwise_arg(op, index) for op in self._args]
        args = args + [
            index,
            self.divisions[index],
            self.divisions[index + 1],
            index == (self.npartitions - 1),
        ]
        return Task(name, self.operation, *args)


class Abs(Elemwise):
    _projection_passthrough = True
    _parameters = ["frame"]
    operation = M.abs


class RenameAxis(Elemwise):
    _projection_passthrough = True
    _filter_passthrough = True
    _parameters = ["frame", "mapper", "index", "columns", "axis"]
    _defaults = {
        "mapper": no_default,
        "index": no_default,
        "columns": no_default,
        "axis": 0,
    }
    _keyword_only = ["mapper", "index", "columns", "axis"]
    operation = M.rename_axis
    _preserves_partitioning_information = True


class NotNull(Elemwise):
    _parameters = ["frame"]
    operation = M.notnull
    _projection_passthrough = True


class ToFrame(Elemwise):
    _parameters = ["frame", "name"]
    _defaults = {"name": no_default}
    _keyword_only = ["name"]
    operation = M.to_frame
    _filter_passthrough = True

    @functools.cached_property
    def unique_partition_mapping_columns_from_shuffle(self):
        result = set()
        name_mapping = dict(zip(self.frame.columns, self.columns))
        for elem in self.frame.unique_partition_mapping_columns_from_shuffle:
            if isinstance(elem, tuple):
                result.add(tuple(name_mapping.get(v, v) for v in elem))
            else:
                result.add(name_mapping.get(elem, elem))
        return result


class ToFrameIndex(ToFrame):
    _parameters = ["frame", "index", "name"]
    _defaults = {"name": no_default, "index": True}
    _keyword_only = ["name", "index"]
    operation = M.to_frame
    _filter_passthrough = True


class ToSeriesIndex(ToFrameIndex):
    _defaults = {"name": no_default, "index": None}
    operation = M.to_series
    _preserves_partitioning_information = True


def pd_split(df, p, random_state=None, shuffle=False):
    p = list(p)
    if shuffle:
        if not isinstance(random_state, np.random.RandomState):
            random_state = np.random.RandomState(random_state)
        df = df.sample(frac=1.0, random_state=random_state)
    index = pseudorandom(len(df), p, random_state)
    if df.ndim == 1:
        df = df.to_frame()
    return df.assign(_split=index)


class Split(Elemwise):
    _parameters = ["frame", "frac", "random_state", "shuffle"]
    _keyword_only = ["random_state", "shuffle"]
    operation = staticmethod(pd_split)

    @functools.cached_property
    def _kwargs(self) -> dict:
        return {"shuffle": self.shuffle}

    @functools.cached_property
    def random_state_data(self):
        return random_state_data(self.frame.npartitions, self.random_state)

    def _task(self, name: Key, index: int) -> Task:
        args = [self._blockwise_arg(op, index) for op in self._args]
        kwargs = self._kwargs.copy()
        kwargs["random_state"] = self.random_state_data[index]
        return Task(name, self.operation, *args, **kwargs)


def _random_split_take(df, i, ndim):
    df = df[df["_split"] == i].drop(columns="_split")
    if ndim == 1:
        return df[df.columns[0]]
    return df


class SplitTake(Blockwise):
    _parameters = ["frame", "i", "ndim"]
    operation = staticmethod(_random_split_take)


class Apply(Elemwise):
    """A good example of writing a less-trivial blockwise operation"""

    _parameters = ["frame", "function", "args", "meta", "kwargs"]
    _defaults = {"args": (), "kwargs": {}}
    operation = M.apply

    @functools.cached_property
    def _meta(self):
        return make_meta(self.operand("meta"), parent_meta=self.frame._meta)

    def _task(self, name: Key, index: int) -> Task:
        return Task(
            name,
            M.apply,
            TaskRef((self.frame._name, index)),
            self.function,
            *self.args,
            **self.kwargs,
        )


class Map(Elemwise):
    _projection_passthrough = True
    _parameters = ["frame", "arg", "na_action", "meta", "is_monotonic"]
    _defaults = {"na_action": None, "meta": None, "is_monotonic": True}
    _keyword_only = ["meta", "is_monotonic"]
    operation = M.map

    @functools.cached_property
    def _meta(self):
        if self.operand("meta") is None:
            args = [
                meta_nonempty(op._meta) if isinstance(op, Expr) else op
                for op in self._args
            ]
            return make_meta(self.operation(*args, **self._kwargs))
        return make_meta(
            self.operand("meta"),
            parent_meta=self.frame._meta,
            index=getattr(self.frame._meta, "index", None),  # could be an index
        )

    @functools.cached_property
    def _kwargs(self) -> dict:
        return {}

    def _divisions(self):
        if not self.is_monotonic:
            # Implement this consistently with dask.dataframe, e.g. add option to
            # control monotonic map func
            return (None,) * len(self.frame.divisions)
        elif is_index_like(self.frame._meta):
            return tuple(
                pd.Series(self.frame.divisions).map(self.arg, na_action=self.na_action)
            )
        elif isinstance(self.arg, Expr):
            if self.arg.divisions == self.frame.divisions:
                return self.frame.divisions
            else:
                # We only get here when we have only one partition
                return (None,) * (self.frame.npartitions + 1)
        return super()._divisions()


class VarColumns(Elemwise):
    _parameters = ["frame", "skipna", "ddof", "numeric_only"]
    _defaults = {"skipna": True, "ddof": 1, "numeric_only": False}
    _keyword_only = ["skipna", "ddof", "numeric_only"]
    operation = M.var
    _is_length_preserving = True

    @functools.cached_property
    def _kwargs(self) -> dict:
        return {"axis": 1, **super()._kwargs}


class NUniqueColumns(Elemwise):
    _parameters = ["frame", "axis", "dropna"]
    _defaults = {"axis": 1, "dropna": True}
    operation = M.nunique


class Sqrt(Elemwise):
    _parameters = ["frame"]
    operation = np.sqrt


class ExplodeSeries(Blockwise):
    _parameters = ["frame"]
    operation = M.explode


class ExplodeFrame(ExplodeSeries):
    _parameters = ["frame", "column"]

    def _simplify_up(self, parent, dependents):
        if isinstance(parent, Projection):
            return plain_column_projection(self, parent, dependents, [self.column])


class Drop(Elemwise):
    _parameters = ["frame", "columns", "errors"]
    _defaults = {"errors": "raise"}
    operation = staticmethod(drop_by_shallow_copy)
    _preserves_partitioning_information = True

    def _simplify_down(self):
        col_op = self.operand("columns")
        if is_scalar(col_op):
            col_op = [col_op]
        columns = [col for col in self.frame.columns if col not in col_op]
        return Projection(self.frame, columns)


def assign(df, *pairs):
    pairs = dict(partition(2, pairs))
    df = df.copy(deep=False)
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="DataFrame is highly fragmented *",
            category=PerformanceWarning,
        )
        for name, val in pairs.items():
            if isinstance(val, Callable):
                val = val(df)
            df[name] = val
    return df


class Assign(Elemwise):
    """Column Assignment"""

    _parameters = ["frame"]
    operation = staticmethod(assign)

    @functools.cached_property
    def unique_partition_mapping_columns_from_shuffle(self):
        keys = set(self.keys)
        return {
            col
            for col in self.frame.unique_partition_mapping_columns_from_shuffle
            if not isinstance(col, tuple)
            and col not in keys
            or not set(col).intersection(keys)
        }

    @functools.cached_property
    def keys(self):
        return self.operands[1::2]

    @functools.cached_property
    def vals(self):
        return self.operands[2::2]

    @functools.cached_property
    def _meta(self):
        args = [op._meta if isinstance(op, Expr) else op for op in self._args]
        return make_meta(self.operation(*args, **self._kwargs))

    def _tree_repr_argument_construction(self, i, op, header):
        if i == 0:
            return super()._tree_repr_argument_construction(i, op, header)
        if i % 2 == 1:
            sep = "" if i == 1 else ","
            header += f"{sep} {repr(op)[1:-1]}="
        else:
            header += f"{repr(op)}"
        return header

    def _node_label_args(self):
        return self.operands

    def _remove_common_columns(self, other):
        if set(self.keys) & set(other.keys):
            keys = set(self.keys)
            operands = [[k, v] for k, v in zip(other.keys, other.vals) if k not in keys]
            return [other.frame] + list(flatten(operands)) + self.operands[1:]
        else:
            return other.operands + self.operands[1:]

    def _simplify_down(self):
        if isinstance(self.frame, Assign):
            if self._check_for_previously_created_column(self.frame):
                # don't squash if we are using a column that was previously created
                return
            return Assign(*self._remove_common_columns(self.frame))
        elif isinstance(self.frame, Projection) and isinstance(
            self.frame.frame, Assign
        ):
            if self._check_for_previously_created_column(self.frame.frame):
                return
            new_columns = self.frame.operands[1].copy()
            new_columns.extend([k for k in self.keys if k not in new_columns])
            return Projection(
                Assign(*self._remove_common_columns(self.frame.frame)), new_columns
            )

    def _check_for_previously_created_column(self, child):
        input_columns = []
        for v in self.vals:
            if isinstance(v, Expr):
                input_columns.extend(v.columns)
        return bool(set(input_columns) & set(child.keys))

    def _simplify_up(self, parent, dependents):
        if isinstance(parent, Projection):
            columns = determine_column_projection(self, parent, dependents)
            columns = _convert_to_list(columns)

            cols = set(columns) - set(self.keys)
            if cols == set(self.frame.columns):
                # Protect against pushing the same projection twice
                return

            diff = set(self.keys) - set(columns)
            if len(diff) == len(self.keys):
                return type(parent)(self.frame, *parent.operands[1:])
            elif len(diff) > 0:
                new_args = []
                for k, v in zip(self.keys, self.vals):
                    if k in columns:
                        new_args.extend([k, v])
            else:
                new_args = self.operands[1:]

            columns = [col for col in self.frame.columns if col in cols]
            return type(parent)(
                type(self)(self.frame[sorted(columns)], *new_args),
                *parent.operands[1:],
            )


class Eval(Elemwise):
    _parameters = ["frame", "_expr", "expr_kwargs"]
    _defaults = {"expr_kwargs": {}}
    _keyword_only = ["expr_kwargs"]
    operation = M.eval

    @functools.cached_property
    def _kwargs(self) -> dict:
        return {**self.expr_kwargs}


class CaseWhen(Elemwise):
    _parameters = ["frame"]

    @functools.cached_property
    def caselist(self):
        c = self.operands[1:]
        return [(c[i], c[i + 1]) for i in range(0, len(c), 2)]

    @functools.cached_property
    def _meta(self):
        c = self.operands[1:]
        caselist = [
            (
                meta_nonempty(c[i]._meta) if isinstance(c[i], Expr) else c[i],
                (
                    meta_nonempty(c[i + 1]._meta)
                    if isinstance(c[i + 1], Expr)
                    else c[i + 1]
                ),
            )
            for i in range(0, len(c), 2)
        ]
        return make_meta(meta_nonempty(self.frame._meta).case_when(caselist))

    @staticmethod
    def operation(ser, *caselist):
        caselist = [(caselist[i], caselist[i + 1]) for i in range(0, len(caselist), 2)]
        return ser.case_when(list(caselist))


class Filter(Blockwise):
    _projection_passthrough = True
    _filter_passthrough = True
    _parameters = ["frame", "predicate"]
    operation = operator.getitem
    _preserves_partitioning_information = True

    def _simplify_up(self, parent, dependents):
        if isinstance(self.predicate, Or):
            result = rewrite_filters(self.predicate)
            if result._name != self.predicate._name:
                return type(parent)(
                    type(self)(self.frame, result), *parent.operands[1:]
                )

        if isinstance(parent, (FilterAlign, Filter)) and not isinstance(
            self.frame, (FilterAlign, Filter)
        ):
            if not self.frame._filter_passthrough_available(self, dependents):
                # We want to collect filters again when we can't move them
                # anymore. Otherwise, a chain of Filters might block Projections
                # We start with the first Filter, e.g. if my child isn't a filter
                # anymore. Filters above that don't know if the first Filter can
                # still move further
                if is_filter_pushdown_available(
                    self, parent, dependents, allow_reduction=False
                ):
                    # We can only squash 2 filters together if the predicate of parent
                    # does not directly depend on self, e.g. if
                    # sum is in the predicate of parent, then removing self would
                    # alter the condition of parent because the sum changes, this is
                    # only relevant in broadcasting cases
                    return self.frame[
                        self.predicate & parent.predicate.substitute(self, self.frame)
                    ]
        if isinstance(parent, Projection):
            if self.frame._filter_passthrough_available(self, dependents):
                # We can't push Projections through filters if the preceding operation
                # allows us to push filters further down the graph because Projections
                # block filter pushdown
                if not isinstance(self.frame, (FilterAlign, Filter)):
                    return
                elif is_filter_pushdown_available(
                    self.frame, self, dependents, allow_reduction=False
                ):
                    return

            return plain_column_projection(self, parent, dependents)
        if isinstance(parent, Index):
            return self.frame.index[self.predicate]


class Projection(Elemwise):
    """Column Selection"""

    _parameters = ["frame", "columns"]
    operation = operator.getitem

    @functools.cached_property
    def unique_partition_mapping_columns_from_shuffle(self):
        col_op = self.operand("columns")
        columns = set(col_op) if isinstance(col_op, list) else {col_op}
        return {
            c
            for c in self.frame.unique_partition_mapping_columns_from_shuffle
            if c in columns or isinstance(c, tuple) and set(c).issubset(columns)
        }

    @property
    def columns(self):
        cols = self.operand("columns")
        if isinstance(cols, list):
            return cols
        elif isinstance(cols, pd.Index):
            return list(cols)
        else:
            return [cols]

    @functools.cached_property
    def _meta(self):
        if is_dataframe_like(self.frame._meta):
            return super()._meta
        # if we are not a DataFrame and have a scalar, we reduce to a scalar
        if not isinstance(self.operand("columns"), (list, slice)) and not hasattr(
            self.operand("columns"), "dtype"
        ):
            return meta_nonempty(self.frame._meta).iloc[0]
        # Avoid column selection for Series/Index
        return self.frame._meta

    def _node_label_args(self):
        return [self.frame, self.operand("columns")]

    def __str__(self):
        base = str(self.frame)
        if " " in base:
            base = "(" + base + ")"
        return f"{base}[{repr(self.operand('columns'))}]"

    def _divisions(self):
        if self.ndim == 0:
            return (None, None)
        return super()._divisions()

    def _simplify_down(self):
        if (
            str(self.frame.columns) == str(self.columns)
            and self._meta.ndim == self.frame._meta.ndim
        ):
            # TODO: we should get more precise around Expr.columns types
            return self.frame
        if isinstance(self.frame, Projection):
            # df[a][b]
            a = self.frame.operand("columns")
            b = self.operand("columns")

            if not isinstance(a, list):
                # df[scalar][b] -> First selection coerces to Series
                return
            elif isinstance(b, list):
                assert all(bb in a for bb in b)
            else:
                assert b in a

            return self.frame.frame[b]


class Index(Elemwise):
    """Column Selection"""

    _parameters = ["frame"]
    operation = getattr

    @functools.cached_property
    def _meta(self):
        meta = self.frame._meta
        # Handle scalar results
        if is_series_like(meta) or is_dataframe_like(meta):
            return self.frame._meta.index
        return meta

    @property
    def _projection_columns(self):
        return []

    def _task(self, name: Key, index: int) -> Task:
        return Task(
            name,
            getattr,
            TaskRef((self.frame._name, index)),
            "index",
        )

    @functools.cached_property
    def unique_partition_mapping_columns_from_shuffle(self):
        name = self.frame._meta.index.name
        if name in self.frame.unique_partition_mapping_columns_from_shuffle:
            return {name}
        elif (name,) in self.frame.unique_partition_mapping_columns_from_shuffle:
            return {(name,)}
        else:
            return set()


def _return_input(df, divisions=None):
    return df


class ClearDivisions(Elemwise):
    _parameters = ["frame"]
    operation = staticmethod(_return_input)
    _preserves_partitioning_information = True

    def _divisions(self):
        return (None,) * (self.frame.npartitions + 1)


class SetDivisions(Elemwise):
    _parameters = ["frame", "divisions"]
    operation = staticmethod(_return_input)
    _preserves_partitioning_information = True

    def _divisions(self):
        return self.operand("divisions")


class ResolveOverlappingDivisions(Expr):
    _parameters = ["frame", "mins", "maxes", "lens"]

    @functools.cached_property
    def _meta(self):
        return self.frame._meta

    def _divisions(self):
        non_empties = [i for i, length in enumerate(self.lens) if length != 0]
        if len(non_empties) == 0:
            return (None, None)

        return tuple(self.mins) + (self.maxes[-1],)

    def _layer(self):
        non_empties = [i for i, length in enumerate(self.lens) if length != 0]
        # If all empty, collapse into one partition
        if len(non_empties) == 0:
            return {(self._name, 0): (self.frame._name, 0)}

        # drop empty partitions by mapping each partition in a new graph to a particular
        # partition on the old graph.
        dsk = {
            (self._name, i): (self.frame._name, div)
            for i, div in enumerate(non_empties)
        }
        ddf_keys = list(dsk.values())

        overlap = [
            i for i in range(1, len(self.mins)) if self.mins[i] >= self.maxes[i - 1]
        ]
        divisions = self.divisions

        frames = []
        for i in overlap:
            # `frames` is a list of data from previous partitions that we may want to
            # move to partition i.  Here, we add "overlap" from the previous partition
            # (i-1) to this list.
            frames.append((get_overlap, ddf_keys[i - 1], divisions[i]))

            # Make sure that any data added from partition i-1 to `frames` is removed
            # from partition i-1.
            dsk[(self._name, i - 1)] = (
                drop_overlap,
                dsk[(self._name, i - 1)],
                divisions[i],
            )

            # We do not want to move "overlap" from the previous partition (i-1) into
            # this partition (i) if the data from this partition will need to be moved
            # to the next partition (i+1) anyway.  If we concatenate data too early,
            # we may lose rows (https://github.com/dask/dask/issues/6972).
            if divisions[i] == divisions[i + 1] and i + 1 in overlap:
                continue

            frames.append(ddf_keys[i])
            dsk[(self._name, i)] = (methods.concat, frames)
            frames = []
        return dsk


class Lengths(Expr):
    """Returns a tuple of partition lengths"""

    _parameters = ["frame"]

    @functools.cached_property
    def _meta(self):
        return tuple()

    def _divisions(self):
        return (None, None)

    def _simplify_down(self):
        if isinstance(self.frame, Elemwise):
            child = max(self.frame.dependencies(), key=lambda expr: expr.npartitions)
            return Lengths(child)

    def _layer(self):
        name = "part-" + self._name
        dsk = {
            (name, i): (len, (self.frame._name, i))
            for i in range(self.frame.npartitions)
        }
        dsk[(self._name, 0)] = (tuple, list(dsk.keys()))
        return dsk


class ResetIndex(Elemwise):
    """Reset the index of a Series or DataFrame"""

    _parameters = ["frame", "drop", "name"]
    _defaults = {"drop": False, "name": no_default}
    _keyword_only = ["drop", "name"]
    operation = M.reset_index
    _filter_passthrough = True
    _preserves_partitioning_information = True

    @functools.cached_property
    def _kwargs(self) -> dict:
        kwargs = {"drop": self.drop}
        if self.operand("name") is not no_default:
            kwargs.update({"name": self.operand("name")})
        return kwargs

    def _divisions(self):
        return (None,) * (self.frame.npartitions + 1)

    def _simplify_up(self, parent, dependents):
        if isinstance(parent, Filter) and self._filter_passthrough_available(
            parent, dependents
        ):
            parents = [
                p().columns
                for p in dependents[self._name]
                if p() is not None and not isinstance(p(), Filter)
            ]
            predicate = None
            if not set(flatten(parents, list)).issubset(set(self.frame.columns)):
                # one of the filters is the Index
                name = self.operand("name") or self.frame._meta.index.name
                if name is no_default and self.frame._meta.index.name is None:
                    name = "index"
                elif self.frame._meta.index.name is not None:
                    name = self.frame._meta.index.name
                # replace the projection of the former index with the actual index
                subs = Projection(self, name)
                predicate = parent.predicate.substitute(subs, Index(self.frame))
            elif self.frame.ndim == 1 and not self.operand("drop"):
                name = self.frame._meta.name
                # Avoid Projection since we are already a Series
                subs = Projection(self, name)
                predicate = parent.predicate.substitute(subs, self.frame)
            return self._filter_simplification(parent, predicate)

        if isinstance(parent, Projection):
            if self.frame.ndim == 1 and not self.drop:
                if isinstance(parent.operand("columns"), list):
                    # Don't bother, dimensionality changes are tricky here and
                    # potential improvement is tiny
                    return
                col = parent.operand("columns")
                if col in (self.name, "index", self.frame._meta.index.name):
                    return
                if all(
                    isinstance(d(), Projection) and d().operand("columns") == col
                    for d in dependents[self._name]
                ):
                    return type(self)(self.frame, True, self.name)
                return
            result = plain_column_projection(self, parent, dependents)
            if result is not None and set(result.columns) != set(result.frame.columns):
                result = result.substitute_parameters({"drop": True})
            return result


class AddPrefixSeries(Elemwise):
    _parameters = ["frame", "prefix"]
    operation = M.add_prefix
    _filter_passthrough = True
    _preserves_partitioning_information = True

    def _divisions(self):
        return tuple(self.prefix + str(division) for division in self.frame.divisions)


class AddSuffixSeries(AddPrefixSeries):
    _parameters = ["frame", "suffix"]
    operation = M.add_suffix
    _preserves_partitioning_information = True

    def _divisions(self):
        return tuple(str(division) + self.suffix for division in self.frame.divisions)


class AddPrefix(Elemwise):
    _parameters = ["frame", "prefix"]
    operation = M.add_prefix

    @functools.cached_property
    def unique_partition_mapping_columns_from_shuffle(self):
        return {
            (
                f"{self.prefix}{c}"
                if not isinstance(c, tuple)
                else tuple(self.prefix + t for t in c)
            )
            for c in self.frame.unique_partition_mapping_columns_from_shuffle
        }

    def _convert_columns(self, columns):
        len_prefix = len(self.prefix)
        return [col[len_prefix:] for col in columns]

    def _simplify_up(self, parent, dependents):
        if isinstance(parent, Projection):
            columns = determine_column_projection(self, parent, dependents)
            columns = self._convert_columns(_convert_to_list(columns))
            if set(columns) == set(self.frame.columns):
                return

            columns = [col for col in self.frame.columns if col in columns]
            return type(parent)(
                type(self)(self.frame[columns], self.operands[1]),
                parent.operand("columns"),
            )


class AddSuffix(AddPrefix):
    _parameters = ["frame", "suffix"]
    operation = M.add_suffix

    @functools.cached_property
    def unique_partition_mapping_columns_from_shuffle(self):
        return {
            (
                f"{c}{self.suffix}"
                if not isinstance(c, tuple)
                else tuple(t + self.suffix for t in c)
            )
            for c in self.frame.unique_partition_mapping_columns_from_shuffle
        }

    def _convert_columns(self, columns):
        len_suffix = len(self.suffix)
        return [col[:-len_suffix] for col in columns]


class AssignIndex(Elemwise):
    _parameters = ["frame", "value"]
    operation = staticmethod(methods.assign_index)
    _preserves_partitioning_information = True

    def _divisions(self):
        return self.value.divisions


class Head(Expr):
    """Take the first `n` rows of the first partition"""

    _parameters = ["frame", "n", "npartitions"]
    _defaults = {"n": 5, "npartitions": 1}

    @functools.cached_property
    def _meta(self):
        return self.frame._meta

    @functools.cached_property
    def npartitions(self):
        return 1

    def _divisions(self):
        if self.operand("npartitions") <= -1:
            return self.frame.divisions[0], self.frame.divisions[-1]
        return (
            self.frame.divisions[0],
            self.frame.divisions[self.operand("npartitions")],
        )

    def _task(self, name: Key, index: int) -> Task:
        raise NotImplementedError()

    def _simplify_down(self):
        if isinstance(self.frame, Elemwise):
            operands = [
                (
                    Head(op, self.n, self.operand("npartitions"))
                    if isinstance(op, Expr) and not isinstance(op, _DelayedExpr)
                    else op
                )
                for op in self.frame.operands
            ]
            return type(self.frame)(*operands)
        if isinstance(self.frame, Head):
            return Head(
                self.frame.frame, min(self.n, self.frame.n), self.operand("npartitions")
            )

    def _simplify_up(self, parent, dependents):
        from dask.dataframe.dask_expr import Repartition

        if isinstance(parent, Repartition) and parent.new_partitions == 1:
            return self

    def _lower(self):
        if not isinstance(self, BlockwiseHead):
            # Lower to Blockwise
            npartitions = self.operand("npartitions")
            if self.operand("npartitions") > self.frame.npartitions:
                raise ValueError(
                    f"only {self.frame.npartitions} partitions, head received {npartitions}"
                )
            partitions = self._partitions
            if is_index_like(self._meta):
                return BlockwiseHeadIndex(
                    Partitions(self.frame, partitions), self.n, safe=False
                )

            safe = True if npartitions == 1 and self.frame.npartitions != 1 else False
            frame = BlockwiseHead(
                Partitions(self.frame, partitions), self.n, npartitions, safe
            )
            if npartitions != 1:
                from dask.dataframe.dask_expr import Repartition

                safe = npartitions != self.frame.npartitions and npartitions != -1
                frame = BlockwiseHead(
                    Repartition(frame, new_partitions=1), self.n, 1, safe
                )
            return frame

    @property
    def _partitions(self):
        if isinstance(self, PartitionsFiltered):
            partitions = self.frame._partitions
        else:
            partitions = list(range(self.frame.npartitions))
        if self.operand("npartitions") > -1:
            partitions = partitions[: self.operand("npartitions")]
        return partitions


class BlockwiseHead(Head, Blockwise):
    """Take the first `n` rows of every partition

    Typically used after `Partitions(..., [0])` to take
    the first `n` rows of an entire collection.
    """

    _parameters = ["frame", "n", "npartitions", "safe"]
    _preserves_partitioning_information = True

    def _simplify_down(self):
        return

    def _simplify_up(self, parent, dependents):
        return

    @functools.cached_property
    def npartitions(self):
        return len(self._divisions()) - 1

    def _divisions(self):
        return self.frame.divisions[: len(self._partitions) + 1]

    def _task(self, name: Key, index: int) -> Task:
        if self.safe:
            op = safe_head
        else:
            op = M.head
        return Task(name, op, TaskRef((self.frame._name, index)), self.n)


class BlockwiseHeadIndex(BlockwiseHead):
    def _task(self, name: Key, index: int) -> Task:
        return Task(
            name, operator.getitem, TaskRef((self.frame._name, index)), slice(0, self.n)
        )


class Tail(Expr):
    """Take the last `n` rows of the last partition"""

    _parameters = ["frame", "n"]
    _defaults = {"n": 5}

    @functools.cached_property
    def _meta(self):
        return self.frame._meta

    def _divisions(self):
        return self.frame.divisions[-2:]

    def _task(self, name: Key, index: int) -> Task:
        raise NotImplementedError()

    def _simplify_down(self):
        if isinstance(self.frame, Elemwise):
            operands = [
                Tail(op, self.n) if isinstance(op, Expr) else op
                for op in self.frame.operands
            ]
            return type(self.frame)(*operands)
        if isinstance(self.frame, Tail):
            return Tail(self.frame.frame, min(self.n, self.frame.n))

    def _simplify_up(self, parent, dependents):
        from dask.dataframe.dask_expr import Repartition

        if isinstance(parent, Repartition) and parent.new_partitions == 1:
            return self

    def _lower(self):
        if not isinstance(self, BlockwiseTail):
            # Lower to Blockwise
            if is_index_like(self._meta):
                return BlockwiseTailIndex(
                    Partitions(self.frame, [self.frame.npartitions - 1]), self.n
                )
            return BlockwiseTail(
                Partitions(self.frame, [self.frame.npartitions - 1]), self.n
            )


class BlockwiseTail(Tail, Blockwise):
    """Take the last `n` rows of every partition

    Typically used after `Partitions(..., [-1])` to take
    the last `n` rows of an entire collection.
    """

    _preserves_partitioning_information = True

    def _divisions(self):
        return self.frame.divisions

    def _task(self, name: Key, index: int) -> Task:
        return Task(name, M.tail, TaskRef((self.frame._name, index)), self.n)


class BlockwiseTailIndex(BlockwiseTail):
    def _task(self, name: Key, index: int) -> Task:
        return Task(
            name,
            operator.getitem,
            TaskRef((self.frame._name, index)),
            slice(-self.n, None),
        )


class Binop(Elemwise):
    _parameters = ["left", "right"]

    @functools.cached_property
    def _broadcastable(self):
        deps = self.dependencies()
        return (
            1 in {dep.npartitions for dep in deps}
            and len({dep.ndim for dep in deps}) == 2
        )

    def __str__(self):
        return f"{self.left} {self._operator_repr} {self.right}"

    def _simplify_up(self, parent, dependents):
        if isinstance(parent, Projection):
            changed = False
            columns = determine_column_projection(self, parent, dependents)
            columns = _convert_to_list(columns)
            columns = [col for col in self.columns if col in columns]
            if (
                isinstance(self.left, Expr)
                and self.left.ndim > 1
                and self.left.columns != columns
            ):
                left = self.left[columns]  # TODO: filter just the correct columns
                changed = True
            else:
                left = self.left
            if (
                isinstance(self.right, Expr)
                and self.right.ndim > 1
                and self.right.columns != columns
            ):
                right = self.right[columns]  # TODO: filter just the correct columns
                changed = True
            else:
                right = self.right
            if not changed:
                return

            return type(parent)(type(self)(left, right), *parent.operands[1:])

    def _node_label_args(self):
        return [self.left, self.right]

    def _divisions(self):
        if is_index_like(self._meta):
            left_divisions = (
                pd.Series(self.left.divisions)
                if isinstance(self.left, Expr)
                else self.left
            )
            right_divisions = (
                pd.Series(self.right.divisions)
                if isinstance(self.right, Expr)
                else self.right
            )

            return tuple(self.operation(left_divisions, right_divisions))
        elif self._broadcastable and len(self.dependencies()) == 2:
            if self.left.ndim < self.right.ndim:
                return self.right.divisions
            else:
                return self.left.divisions
        else:
            return super()._divisions()


class Add(Binop):
    operation = operator.add
    _operator_repr = "+"


class MethodOperator(Binop):
    _parameters = ["name", "left", "right", "axis", "level", "fill_value"]
    _defaults = {"axis": "columns", "level": None, "fill_value": None}
    _keyword_only = ["axis", "level", "fill_value"]

    @property
    def _operator_repr(self):
        return self.name

    @staticmethod
    def operation(name, left, right, **kwargs):
        return getattr(left, name)(right, **kwargs)


class Sub(Binop):
    operation = operator.sub
    _operator_repr = "-"


class Mul(Binop):
    operation = operator.mul
    _operator_repr = "*"

    def _simplify_down(self):
        if (
            isinstance(self.right, Mul)
            and isinstance(self.left, numbers.Number)
            and isinstance(self.right.left, numbers.Number)
        ):
            return (self.left * self.right.left) * self.right.right


class Pow(Binop):
    operation = operator.pow
    _operator_repr = "**"


class Div(Binop):
    operation = operator.truediv
    _operator_repr = "/"


class LT(Binop):
    operation = operator.lt
    _operator_repr = "<"


class BinOpSeries(Binop):
    _parameters = ["left", "right", "level", "fill_value"]
    _defaults = {"fill_value": None, "level": None}


class BinOpFrame(Binop):
    _parameters = ["left", "right", "axis"]
    _defaults = {"axis": 1}


class LTSeries(BinOpSeries):
    operation = M.lt
    _operator_repr = "<"


class LTFrame(BinOpFrame):
    operation = M.lt
    _operator_repr = "<"


class LESeries(BinOpSeries):
    operation = M.le
    _operator_repr = "<="


class LEFrame(BinOpFrame):
    operation = M.le
    _operator_repr = "<="


class GTSeries(BinOpSeries):
    operation = M.gt
    _operator_repr = ">"


class GTFrame(BinOpFrame):
    operation = M.gt
    _operator_repr = ">"


class GESeries(BinOpSeries):
    operation = M.ge
    _operator_repr = ">="


class GEFrame(BinOpFrame):
    operation = M.ge
    _operator_repr = ">="


class NESeries(BinOpSeries):
    operation = M.ne
    _operator_repr = "!="


class NEFrame(BinOpFrame):
    operation = M.ne
    _operator_repr = "!="


class EQSeries(BinOpSeries):
    operation = M.eq
    _operator_repr = "=="


class EQFrame(BinOpFrame):
    operation = M.eq
    _operator_repr = "=="


class LE(Binop):
    operation = operator.le
    _operator_repr = "<="


class GT(Binop):
    operation = operator.gt
    _operator_repr = ">"


class GE(Binop):
    operation = operator.ge
    _operator_repr = ">="


class EQ(Binop):
    operation = operator.eq
    _operator_repr = "=="


class NE(Binop):
    operation = operator.ne
    _operator_repr = "!="


class And(Binop):
    operation = operator.and_
    _operator_repr = "&"


class Or(Binop):
    operation = operator.or_
    _operator_repr = "|"


class XOr(Binop):
    operation = operator.xor
    _operator_repr = "^"


class Mod(Binop):
    operation = operator.mod
    _operator_repr = "%"


class FloorDiv(Binop):
    operation = operator.floordiv
    _operator_repr = "//"


class Unaryop(Elemwise):
    _parameters = ["frame"]

    def __str__(self):
        return f"{self._operator_repr} {self.frame}"

    def _simplify_up(self, parent, dependents):
        if isinstance(parent, Projection):
            if isinstance(self.frame, Expr):
                return plain_column_projection(self, parent, dependents)
            else:
                frame = self.frame
            return type(self)(frame)

    def _node_label_args(self):
        return [self.frame]

    def _divisions(self):
        if is_index_like(self._meta):
            return (None,) * (self.frame.npartitions + 1)
        else:
            return super()._divisions()


class Invert(Unaryop):
    operation = operator.inv
    _operator_repr = "~"


class Neg(Unaryop):
    operation = operator.neg
    _operator_repr = "-"


class Pos(Unaryop):
    operation = operator.pos
    _operator_repr = "+"


class Partitions(Expr):
    """Select one or more partitions"""

    _parameters = ["frame", "partitions"]

    @functools.cached_property
    def _meta(self):
        return self.frame._meta

    def _divisions(self):
        divisions = []
        for part in self.partitions:
            divisions.append(self.frame.divisions[part])
        divisions.append(self.frame.divisions[part + 1])
        return tuple(divisions)

    def _task(self, name: Key, index: int) -> Task:
        return Alias(name, (self.frame._name, self.partitions[index]))  # type: ignore

    def _simplify_down(self):
        from dask.dataframe.dask_expr import SetIndexBlockwise
        from dask.dataframe.tseries.resample import ResampleAggregation

        if isinstance(self.frame, Blockwise) and not isinstance(
            self.frame, (BlockwiseIO, Fused, SetIndexBlockwise, ResampleAggregation)
        ):
            operands = [
                (
                    Partitions(op, self.partitions)
                    if (isinstance(op, Expr) and not self.frame._broadcast_dep(op))
                    else op
                )
                for op in self.frame.operands
            ]
            return type(self.frame)(*operands)
        elif isinstance(self.frame, PartitionsFiltered):
            if self.frame._partitions:
                partitions = [self.frame._partitions[p] for p in self.partitions]
            else:
                partitions = self.partitions
            # We assume that expressions defining a special "_partitions"
            # parameter can internally capture the same logic as `Partitions`
            return self.frame.substitute_parameters({"_partitions": partitions})

    def _node_label_args(self):
        return [self.frame, self.partitions]


class PartitionsFiltered(Expr):
    """Mixin class for partition filtering

    A ``PartitionsFiltered`` subclass must define a ``_partitions`` parameter. When
    ``_partitions`` is defined, the following expressions must produce the same output
    for :cls:`PartitionsFiltered`:

    - ``cls(expr: Expr, ..., _partitions)``
    - ``Partitions(cls(expr: Expr, ...), _partitions)``

    In order to leverage the default ``Expr._layer`` method, subclasses should define
    ``_filtered_task`` instead of ``_task``.
    """

    @property
    def _filtered(self) -> bool:
        """Whether output partitions have been filtered"""
        return self.operand("_partitions") is not None

    @property
    def _partitions(self) -> list | tuple | range:
        """Selected partition indices"""
        if self._filtered:
            return self.operand("_partitions")
        else:
            return range(self.npartitions)

    @functools.cached_property
    def divisions(self):
        # Common case: Use self._divisions()
        full_divisions = super().divisions
        if not self._filtered:
            return full_divisions

        # Specific case: Specific partitions were selected
        new_divisions = []
        for part in self._partitions:
            new_divisions.append(full_divisions[part])
        new_divisions.append(full_divisions[part + 1])
        return tuple(new_divisions)

    @property
    def npartitions(self):
        if self._filtered:
            return len(self._partitions)
        return super().npartitions

    def _task(self, name: Key, index: int) -> Task:
        return self._filtered_task(name, self._partitions[index])

    def _filtered_task(self, name: Key, index: int) -> Task:
        raise NotImplementedError()


class _DelayedExpr(Expr):
    _parameters = ["obj"]

    def __str__(self):
        return f"{type(self).__name__}({str(self.obj)})"

    @property
    def _name(self):
        return self.obj.key

    def _layer(self) -> dict:
        dc = self.obj.__dask_optimize__(self.obj.dask, self.obj.key).to_dict().copy()
        dc[(self.obj.key, 0)] = dc[self.obj.key]
        dc.pop(self.obj.key)
        return dc

    def _divisions(self):
        return (None, None)

    @property
    def ndim(self):
        return 0


class DelayedsExpr(Expr):
    def __str__(self):
        return f"{type(self).__name__}({str(self.operands[0])})"

    @functools.cached_property
    def _name(self):
        return "delayed-container-" + self.deterministic_token

    def _layer(self) -> dict:
        from dask.delayed import Delayed

        if isinstance(self.operands[0], TaskRef):
            tasks = [
                Alias((self._name, ix), fut.key) for ix, fut in enumerate(self.operands)
            ]
            dsk = {t.key: t for t in tasks}
        elif isinstance(self.operands[0], Delayed):
            expr = collections_to_expr(self.operands).optimize()
            keys = expr.__dask_keys__()
            dsk = expr.__dask_graph__()
            # Many APIs in dask-expr are not honoring __dask_keys__ but are instead
            # assuming they can just construct the keys themselves by walking the
            # partitions. Therefore we'll have to remap the key names and can't just
            # expose __dask_keys__()
            for ix, actual_key in enumerate(keys):
                dsk[(self._name, ix)] = Alias((self._name, ix), actual_key[0])
        else:
            raise TypeError("Expected a Delayed or Future object")

        return dsk

    def _divisions(self):
        return (None,) * (len(self.operands) + 1)

    @property
    def ndim(self):
        return 0


def is_broadcastable(dfs, s):
    """
    This Series is broadcastable against another dataframe in the sequence
    """

    def compare(s, df):
        try:
            return s.divisions == (min(df.columns), max(df.columns))
        except (TypeError, ValueError):
            return False

    return (
        s.ndim == 1
        and s.npartitions == 1
        and s.known_divisions
        and any(compare(s, df) for df in dfs if df.ndim == 2)
        or s.ndim == 0
    )


def are_co_aligned(*exprs):
    """Do inputs come from the same parents, modulo blockwise?"""

    from dask.dataframe.dask_expr._cumulative import CumulativeAggregations
    from dask.dataframe.dask_expr._reductions import Reduction

    seen = set()
    # Scalars can always be broadcasted
    stack = [e for e in exprs if e.ndim > 0]
    ancestors = []
    while stack:
        e = stack.pop()
        if e._name in seen:
            continue
        seen.add(e._name)

        if isinstance(e, IO):
            ancestors.append(e)
        elif e.ndim == 0:
            # Scalars are valid ancestors that are always broadcastable,
            # so don't walk through them
            continue
        elif isinstance(e, (_DelayedExpr, Isin)):
            continue
        elif isinstance(e, (Blockwise, CumulativeAggregations, Reduction)):
            # TODO: Capture this in inheritance logic
            dependencies = e.dependencies()
            stack.extend(dependencies)
        else:
            ancestors.append(e)

    unique_ancestors = {
        # Account for column projection within IO expressions
        _tokenize_partial(item, ["columns", "_series", "_dataset_info_cache"])
        for item in ancestors
    }
    # Don't check divisions or npartitions at all
    return len(unique_ancestors) <= 1


## Utilities for Expr fusion


def is_valid_blockwise_op(expr):
    return isinstance(expr, Blockwise) and not isinstance(
        expr, (FromPandas, FromArray, FromDelayed)
    )


def optimize_blockwise_fusion(expr):
    """Traverse the expression graph and apply fusion"""

    def _fusion_pass(expr):
        # Full pass to find global dependencies
        seen = set()
        stack = [expr]
        dependents = defaultdict(set)
        dependencies = {}
        expr_mapping = {}

        while stack:
            next = stack.pop()

            if next._name in seen:
                continue
            seen.add(next._name)

            if is_valid_blockwise_op(next):
                dependencies[next._name] = set()
                if next._name not in dependents:
                    dependents[next._name] = set()
                    expr_mapping[next._name] = next

            for operand in next.dependencies():
                stack.append(operand)
                if is_valid_blockwise_op(operand):
                    if next._name in dependencies:
                        dependencies[next._name].add(operand._name)
                    dependents[operand._name].add(next._name)
                    expr_mapping[operand._name] = operand
                    expr_mapping[next._name] = next

        # Traverse each "root" until we find a fusable sub-group.
        # Here we use root to refer to a Blockwise Expr node that
        # has no Blockwise dependents
        roots = [
            expr_mapping[k]
            for k, v in dependents.items()
            if v == set()
            or all(not is_valid_blockwise_op(expr_mapping[_expr]) for _expr in v)
        ]
        while roots:
            root = roots.pop()
            seen = set()
            stack = [root]
            group = []
            while stack:
                next = stack.pop()

                if next._name in seen:
                    continue
                seen.add(next._name)

                group.append(next)
                for dep_name in dependencies[next._name]:
                    dep = expr_mapping[dep_name]

                    stack_names = {s._name for s in stack}
                    group_names = {g._name for g in group}
                    if (
                        dep.npartitions == root.npartitions or next._broadcast_dep(dep)
                    ) and not (dependents[dep._name] - stack_names - group_names):
                        # All of deps dependents are contained
                        # in the local group (or the local stack
                        # of expr nodes that we know we will be
                        # adding to the local group).
                        # All nodes must also have the same number
                        # of partitions, since broadcasting within
                        # a group is not allowed.
                        stack.append(dep)
                    elif dependencies[dep._name] and dep._name not in [
                        r._name for r in roots
                    ]:
                        # Couldn't fuse dep, but we may be able to
                        # use it as a new root on the next pass
                        roots.append(dep)

            # Replace fusable sub-group
            if len(group) > 1:
                group_deps = []
                local_names = [_expr._name for _expr in group]
                for _expr in group:
                    group_deps += [
                        operand
                        for operand in _expr.dependencies()
                        if operand._name not in local_names
                    ]
                _ret = expr.substitute(group[0], Fused(group, *group_deps))
                return _ret, not roots

        # Return original expr if no fusable sub-groups were found
        return expr, True

    while True:
        original_name = expr._name
        expr, done = _fusion_pass(expr)
        if done or expr._name == original_name:
            break

    return expr


class Diff(MapOverlap):
    _parameters = ["frame", "periods"]
    _defaults = {"periods": 1}
    func = M.diff
    enforce_metadata = True
    transform_divisions = False
    clear_divisions = False
    align_dataframes = True

    def _divisions(self):
        return self.frame.divisions

    @functools.cached_property
    def _meta(self):
        return make_meta(meta_nonempty(self.frame._meta).diff(**self.kwargs))

    def _simplify_up(self, parent, dependents):
        if isinstance(parent, Projection):
            return plain_column_projection(self, parent, dependents)

    @functools.cached_property
    def kwargs(self):
        return dict(periods=self.periods)

    @property
    def before(self):
        return self.periods if self.periods > 0 else 0

    @property
    def after(self):
        return 0 if self.periods > 0 else -self.periods


class FillnaCheck(Blockwise):
    _parameters = ["frame", "method", "skip_check"]
    operation = staticmethod(methods.fillna_check)
    _projection_passthrough = True

    @functools.cached_property
    def _meta(self):
        return self.frame._meta

    def _task(self, name: Key, index: int) -> Task:
        args = [self._blockwise_arg(op, index) for op in self._args]
        args[-1] = index != self.skip_check(self.frame)
        return Task(name, self.operation, *args)


class FFill(MapOverlap):
    _parameters = ["frame", "limit"]
    _defaults = {"limit": None}
    func = M.ffill
    enforce_metadata = True
    transform_divisions = False
    clear_divisions = False
    align_dataframes = True

    def _divisions(self):
        return self.frame.divisions

    @functools.cached_property
    def _meta(self):
        return self.frame._meta

    def _simplify_up(self, parent, dependents):
        if isinstance(parent, Projection):
            return plain_column_projection(self, parent, dependents)

    @functools.cached_property
    def kwargs(self):
        return dict(limit=self.limit)

    @property
    def before(self):
        return 1 if self.limit is None else self.limit

    @property
    def after(self):
        return 0


class BFill(FFill):
    func = M.bfill

    @property
    def before(self):
        # bfill is the opposite direction of ffill, so
        # we swap before with after of ffill.
        return super().after

    @property
    def after(self):
        # bfill is the opposite direction of ffill, so
        # we swap after with before of ffill.
        return super().before


class Shift(MapOverlap):
    _parameters = ["frame", "periods", "freq"]
    _defaults = {"periods": 1, "freq": None}

    func = M.shift
    enforce_metadata = True
    transform_divisions = False
    align_dataframes = True

    @functools.cached_property
    def clear_divisions(self):
        # TODO We can do better if freq is given, but this needs adjustments in
        #  map_partitions
        return True if self._divisions()[0] is None or self.freq is not None else False

    def _divisions(self):
        if self.freq is None:
            return self.frame.divisions
        divisions = _calc_maybe_new_divisions(self.frame, self.periods, self.freq)
        if divisions is None:
            divisions = (None,) * (self.frame.npartitions + 1)
        return divisions

    @functools.cached_property
    def _meta(self):
        return make_meta(meta_nonempty(self.frame._meta).shift(**self.kwargs))

    @functools.cached_property
    def kwargs(self):
        return dict(periods=self.periods, freq=self.freq)

    def _simplify_up(self, parent, dependents):
        if isinstance(parent, Projection):
            return plain_column_projection(self, parent, dependents)

    @property
    def before(self):
        return self.periods if self.periods > 0 else 0

    @property
    def after(self):
        return 0 if self.periods > 0 else -self.periods


class ShiftIndex(Blockwise):
    _parameters = ["frame", "periods", "freq"]
    _defaults = {"periods": 1, "freq": None}
    _keyword_only = ["freq"]
    operation = M.shift

    def _divisions(self):
        freq = self.freq
        if freq is None:
            freq = self._meta.freq
        divisions = _calc_maybe_new_divisions(self.frame, self.periods, freq)
        if divisions is None:
            divisions = (None,) * (self.frame.npartitions + 1)
        return divisions

    @functools.cached_property
    def _kwargs(self) -> dict:
        return {"freq": self.freq} if self.freq is not None else {}


class MaybeAlignPartitions(Expr):
    _projection_passthrough = False
    _expr_cls: AnyType | None = None

    def _divisions(self):
        if {df.npartitions for df in self.args} == {1}:
            divs = []
            for df in self.args:
                divs.extend(list(df.divisions))
            try:
                return min(divs), max(divs)
            except TypeError:
                # either unknown divisions or int-str mix
                return None, None
        return calc_divisions_for_align(*self.args)

    def _simplify_up(self, parent, dependents):
        if isinstance(parent, Projection) and self._projection_passthrough:
            return plain_column_projection(self, parent, dependents)

    @functools.cached_property
    def args(self):
        dfs = [op for op in self.operands if isinstance(op, Expr)]
        return [op for op in dfs if not is_broadcastable(dfs, op)]

    def _lower(self):
        # This can be expensive when something that has expensive division
        # calculation is in the Expression
        dfs = self.args
        if (
            len(dfs) == 1
            or all(
                dfs[0].divisions == df.divisions and df.known_divisions for df in dfs
            )
            or len(self.divisions) == 2
            and max(map(lambda x: len(x.divisions), dfs)) == 2
        ):
            return self._expr_cls(*self.operands)
        elif self.divisions[0] is None:
            # We have to shuffle
            npartitions = max(df.npartitions for df in dfs)
            dtypes = {df._meta.index.dtype for df in dfs}
            if not _are_dtypes_shuffle_compatible(dtypes):
                raise TypeError(
                    "DataFrames are not aligned. We need to shuffle to align partitions "
                    "with each other. This is not possible because the indexes of the "
                    f"DataFrames have differing dtypes={dtypes}. Please ensure that "
                    "all Indexes have the same dtype or align manually for this to "
                    "work."
                )

            from dask.dataframe.dask_expr._shuffle import RearrangeByColumn

            args = [
                (
                    RearrangeByColumn(df, None, npartitions, index_shuffle=True)
                    if isinstance(df, Expr)
                    else df
                )
                for df in self.operands
            ]
            return self._expr_cls(*args)

        args = maybe_align_partitions(*self.operands, divisions=self.divisions)
        return self._expr_cls(*args)

    @functools.cached_property
    def _meta(self):
        return self._expr_cls(*self.operands)._meta


def _are_dtypes_shuffle_compatible(dtypes):
    if len(dtypes) == 1:
        return True
    return all(pd.api.types.is_numeric_dtype(d) for d in dtypes)


class CombineFirstAlign(MaybeAlignPartitions):
    _parameters = ["frame", "other"]
    _expr_cls = CombineFirst

    def _simplify_up(self, parent, dependents):
        # TODO: de-duplicate
        if isinstance(parent, Projection):
            columns = determine_column_projection(self, parent, dependents)
            frame_columns = [col for col in self.frame.columns if col in columns]
            other_columns = [col for col in self.other.columns if col in columns]
            if (
                self.frame.columns == frame_columns
                and self.other.columns == other_columns
            ):
                return

            return type(parent)(
                type(self)(self.frame[frame_columns], self.other[other_columns]),
                *parent.operands[1:],
            )


class FillnaAlign(MaybeAlignPartitions):
    _projection_passthrough = True
    _parameters = ["frame", "value"]
    _expr_cls = Fillna


class AlignAlignPartitions(MaybeAlignPartitions):
    _parameters = ["frame", "other", "join", "axis", "fill_value"]
    _expr_cls = _Align


class CombineSeriesAlign(MaybeAlignPartitions):
    _parameters = ["frame", "other", "func", "fill_value"]
    _expr_cls = CombineSeries


class CombineFrameAlign(MaybeAlignPartitions):
    _parameters = ["frame", "other", "func", "fill_value", "overwrite"]
    _expr_cls = CombineSeries


class FilterAlign(MaybeAlignPartitions):
    _projection_passthrough = True
    _filter_passthrough = True
    _parameters = ["frame", "predicate"]
    _expr_cls = Filter

    def _simplify_up(self, parent, dependents):
        return Filter._simplify_up(self, parent, dependents)


class AssignAlign(MaybeAlignPartitions):
    _parameters = ["frame", "column", "value"]
    _expr_cls = Assign

    def _simplify_up(self, parent, dependents):
        # TODO: de-duplicate
        if isinstance(parent, Projection):
            columns = determine_column_projection(self, parent, dependents)
            if not isinstance(columns, list):
                columns = [columns]

            cols = set(columns) - {self.column}
            if cols == set(self.frame.columns):
                # Protect against pushing the same projection twice
                return

            diff = {self.column} - set(columns)
            if len(diff) == 1:
                return type(parent)(self.frame, *parent.operands[1:])
            else:
                new_args = self.operands[1:]

            columns = [col for col in self.frame.columns if col in cols]
            return type(parent)(
                type(self)(self.frame[sorted(columns)], *new_args),
                *parent.operands[1:],
            )


class MaskAlign(MaybeAlignPartitions):
    _parameters = ["frame", "cond", "other"]
    _expr_cls: AnyType = Mask


class WhereAlign(MaskAlign):
    _expr_cls = Where


class MapAlign(MaybeAlignPartitions):
    _parameters = ["frame", "other", "op", "na_action", "meta"]
    _projection_passthrough = False
    _expr_cls = Map


class MapIndexAlign(MapAlign):
    _parameters = MaskAlign._parameters + ["is_monotonic"]


class OpAlignPartitions(MaybeAlignPartitions):
    _parameters = ["frame", "other", "op"]
    _projection_passthrough = True

    @functools.cached_property
    def _meta(self):
        return getattr(self.frame._meta, self.op)(self.other._meta)

    def _lower(self):
        # This can be expensive when something that has expensive division
        # calculation is in the Expression
        dfs = self.args
        if (
            len(dfs) == 1
            or all(dfs[0].divisions == df.divisions for df in dfs)
            or len(self.divisions) == 2
            and max(map(lambda x: len(x.divisions), dfs)) == 2
        ):
            return self._op(self.frame, self.op, self.other, *self.operands[3:])

        from dask.dataframe.dask_expr._repartition import RepartitionDivisions

        frame = RepartitionDivisions(
            self.frame, new_divisions=self.divisions, force=True
        )
        other = RepartitionDivisions(
            self.other, new_divisions=self.divisions, force=True
        )
        return self._op(frame, self.op, other, *self.operands[3:])

    @staticmethod
    def _op(frame, op, other, *args, **kwargs):
        return getattr(frame, op)(other)


class MethodOperatorAlign(OpAlignPartitions):
    _parameters = ["frame", "other", "op", "axis", "level", "fill_value"]

    @staticmethod
    def _op(frame, op, other, *args, **kwargs):
        return MethodOperator(op, frame, other, *args, **kwargs)


class UFuncAlign(MaybeAlignPartitions):
    _parameters = ["frame", "func", "meta", "kwargs"]
    enforce_metadata = False

    def __str__(self):
        return f"UFunc({funcname(self.func)})"

    @functools.cached_property
    def args(self):
        return self.operands[len(self._parameters) :]

    @functools.cached_property
    def _dfs(self):
        return [df for df in self.args if isinstance(df, Expr)]

    @functools.cached_property
    def _meta(self):
        if self.operand("meta") is not no_default:
            return self.operand("meta")
        return _get_meta_ufunc(self._dfs, self.args, self.func)

    def _lower(self):
        args = maybe_align_partitions(*self.args, divisions=self._divisions())
        dfs = [x for x in args if isinstance(x, Expr) and x.ndim > 0]
        return UFuncElemwise(dfs[0], self.func, self._meta, False, self.kwargs, *args)


class Fused(Blockwise):
    """Fused ``Blockwise`` expression

    A ``Fused`` corresponds to the fusion of multiple
    ``Blockwise`` expressions into a single ``Expr`` object.
    Before graph-materialization time, the behavior of this
    object should be identical to that of the first element
    of ``Fused.exprs`` (i.e. the top-most expression in
    the fused group).

    Parameters
    ----------
    exprs : List[Expr]
        Group of original ``Expr`` objects being fused together.
    *dependencies:
        List of external ``Expr`` dependencies. External-``Expr``
        dependencies correspond to any ``Expr`` operand that is
        not already included in ``exprs``. Note that these
        dependencies should be defined in the order of the ``Expr``
        objects that require them (in ``exprs``). These
        dependencies do not include literal operands, because those
        arguments should already be captured in the fused subgraphs.
    """

    _parameters = ["exprs"]

    @functools.cached_property
    def _meta(self):
        return self.exprs[0]._meta

    def _tree_repr_lines(self, indent=0, recursive=True):
        header = f"Fused({self._name[-5:]}):"
        if not recursive:
            return [header]

        seen = set()
        lines = []
        stack = [(self.exprs[0], 2)]
        fused_group = [_expr._name for _expr in self.exprs]
        dependencies = {dep._name: dep for dep in self.dependencies()}
        while stack:
            expr, _indent = stack.pop()

            if expr._name in seen:
                continue
            seen.add(expr._name)

            line = expr._tree_repr_lines(_indent, recursive=False)[0]
            lines.append(line.replace(" ", "|", 1))
            for dep in expr.dependencies():
                if dep._name in fused_group:
                    stack.append((dep, _indent + 2))
                elif dep._name in dependencies:
                    dependencies.pop(dep._name)
                    lines.extend(dep._tree_repr_lines(_indent + 2))

        for dep in dependencies.values():
            lines.extend(dep._tree_repr_lines(2))

        lines = [header] + lines
        lines = [" " * indent + line for line in lines]

        return lines

    def __str__(self):
        exprs = sorted(self.exprs, key=M._depth)
        names = [expr._name.split("-")[0] for expr in exprs]
        if len(names) > 4:
            return names[0] + "-fused-" + names[-1]
        else:
            return "-".join(names)

    @functools.cached_property
    def _name(self):
        return f"{str(self)}-{self.deterministic_token}"

    def _divisions(self):
        return self.exprs[0]._divisions()

    def _broadcast_dep(self, dep: Expr):
        # Always broadcast single-partition dependencies in Fused
        return dep.npartitions == 1

    def _task(self, name: Key, index: int) -> Task:
        internal_tasks = []
        for _expr in self.exprs:
            if self._broadcast_dep(_expr):
                subname = (_expr._name, 0)
            else:
                subname = (_expr._name, index)
            t = _expr._task(subname, subname[1])

            assert t.key == subname
            internal_tasks.append(t)
        return Task.fuse(*internal_tasks, key=name)  # type: ignore

    @staticmethod
    def _execute_internal_graph(internal_tasks, dependencies, outkey):
        cache = dict(dependencies)
        res = execute_graph(internal_tasks, cache=cache, keys=[outkey])
        return res[outkey]


# Used for sorting with None
@functools.total_ordering
class MinType:
    def __le__(self, other):
        return True


def determine_column_projection(
    expr: Expr,
    parent: Expr,
    dependents: dict[str, Collection[weakref.ref[BaseExpr]]],
    additional_columns: list | None = None,
) -> object:
    if isinstance(parent, Index):
        column_union = []
    else:
        column_union = parent.columns.copy()
    parents: list[Expr]
    parents = [inst for x in dependents[expr._name] if isinstance((inst := x()), Expr)]

    seen = set()
    for p in parents:
        if p._name in seen:
            continue
        seen.add(p._name)
        column_union.extend(p._projection_columns)

    if additional_columns is not None:
        column_union.extend(flatten(additional_columns, container=list))

    # We can end up with MultiIndex columns from groupby ops, needs to be
    # accounted for in the sort
    flattened_columns = set(column_union)
    try:
        column_union = sorted(flattened_columns)
    except TypeError:
        # mixed type columns
        column_union = _sort_mixed(pd.Index(list(flattened_columns))).tolist()
    if (
        len(column_union) == 1
        and parent.ndim == 1
        and all(p.ndim == 1 for p in parents)
    ):
        return column_union[0]
    return column_union


def _sort_mixed(values):
    """order ints before strings before nulls in 1d arrays"""
    str_pos = np.array([isinstance(x, str) for x in values], dtype=bool)
    tuple_pos = np.array([isinstance(x, tuple) for x in values], dtype=bool)
    null_pos = np.array([pd.isna(x) for x in values], dtype=bool)
    num_pos = ~str_pos & ~null_pos & ~tuple_pos
    str_argsort = np.argsort(values[str_pos])
    tuple_argsort = np.argsort(values[tuple_pos])
    num_argsort = np.argsort(values[num_pos])
    # convert boolean arrays to positional indices, then order by underlying values
    str_locs = str_pos.nonzero()[0].take(str_argsort)
    tuple_locs = tuple_pos.nonzero()[0].take(tuple_argsort)
    num_locs = num_pos.nonzero()[0].take(num_argsort)
    null_locs = null_pos.nonzero()[0]
    locs = np.concatenate([num_locs, str_locs, tuple_locs, null_locs])
    return values.take(locs)


def plain_column_projection(expr, parent, dependents, additional_columns=None):
    column_union = determine_column_projection(
        expr, parent, dependents, additional_columns=additional_columns
    )
    if isinstance(column_union, list):
        column_union = [col for col in expr.frame.columns if col in column_union]
    elif column_union not in expr.frame.columns:
        # we are accessing the index
        column_union = []

    if column_union == expr.frame.columns or not column_union and expr.ndim < 2:
        # this projection is for the index, but the elements are unknown, so
        # don't project
        return
    result = type(expr)(expr.frame[column_union], *expr.operands[1:])
    if column_union == parent.operand("columns"):
        return result
    return type(parent)(result, parent.operand("columns"))


def is_filter_pushdown_available(expr, parent, dependents, allow_reduction=True):
    parents = [x() for x in dependents[expr._name] if x() is not None]
    filters = {e._name for e in parents if isinstance(e, Filter)}
    if len(filters) != 1:
        # Don't push down if not exactly one Filter
        return False
    if len(parents) == 1:
        return True

    # We have to see if the non-filter ops are all exclusively part of the predicates
    others = {e._name for e in parents if not isinstance(e, Filter)}
    return _check_dependents_are_predicates(
        expr, others, parent, dependents, allow_reduction
    )


def rewrite_filters(predicate):
    """Rewriting a filter to decompose OR clauses. If a predicate part is part of
    all OR clauses, we can move it to the front so that we can push it down.
    """
    or_components = _get_predicate_components(predicate, [])
    if len(or_components) == 1:
        return predicate
    result = _replace_common_or_components(or_components[0], or_components[1:])
    if result is None:
        return predicate
    return result


def _get_predicate_components(predicate, components, type_=Or):
    if not isinstance(predicate, type_):
        components.append(predicate)
        return components
    if isinstance(predicate.left, type_):
        components = _get_predicate_components(predicate.left, components, type_)
    else:
        components.append(predicate.left)
    if isinstance(predicate.right, type_):
        components = _get_predicate_components(predicate.right, components, type_)
    else:
        components.append(predicate.right)
    return components


def _convert_mapping(components):
    return dict(zip([e._name for e in components], components))


def _replace_common_or_components(expr, or_components):
    and_component = _get_predicate_components(expr, [], type_=And)
    mapping = _convert_mapping(and_component)
    and_components = [
        _get_predicate_components(c, [], type_=And) for c in or_components
    ]
    and_components = list(map(_convert_mapping, and_components))

    replacements = []
    for c in mapping.keys():
        if all(c in comp for comp in and_components):
            # We can pull this component out if it's part of all or components
            replacements.append(c)
    if len(replacements) == 0:
        # exit if we can't replace anything
        return

    outer_component = mapping[replacements[0]]
    for r in replacements[1:]:
        # construct the outer component
        outer_component = outer_component & mapping[r]

    #
    result_components = []
    for comp in [mapping] + and_components:
        keep_components = [c for c in comp if c not in replacements]
        if len(keep_components) == 0:
            # Just return outer_component if we can replace a whole OR component
            return outer_component
        result_component = comp[keep_components[0]]
        for c in keep_components[1:]:
            result_component = result_component & comp[c]
        result_components.append(result_component)

    or_component = result_components[0]
    for c in result_components[1:]:
        or_component = or_component | c
    return outer_component & or_component


def _check_dependents_are_predicates(
    expr, other_names, parent: Expr, dependents, allow_reduction=True
):
    # singleton approach should make this easier

    # Walk down the predicate side from the filter to see if we can arrive at
    # other_names without hitting an expression that has other dependents that
    # are not part of the predicate, see test_filter_pushdown_unavailable
    allowed_expressions = {parent._name}
    stack = parent.dependencies()
    seen = set()
    all_dependents = set()

    while stack:
        e = stack.pop()
        if expr._name == e._name:
            continue

        if e._name in seen:
            continue
        seen.add(e._name)

        if isinstance(e, _DelayedExpr):
            continue

        all_dependents.update(
            {x()._name for x in dependents[e._name] if x() is not None}
        )

        if not allow_reduction:
            if isinstance(e, (ApplyConcatApply, TreeReduce, ShuffleReduce)):
                return False

        allowed_expressions.add(e._name)
        stack.extend(e.dependencies())

    return all_dependents.issubset(allowed_expressions) and other_names.issubset(
        allowed_expressions
    )


def calc_divisions_for_align(*exprs, allow_shuffle=True):
    dfs = [df for df in exprs if isinstance(df, Expr) and df.ndim > 0]
    if not all(df.known_divisions for df in dfs):
        return (None,) * (max(df.npartitions for df in dfs) + 1)
    if all(dfs[0].divisions == df.divisions for df in dfs):
        return dfs[0].divisions
    divisions = list(unique(merge_sorted(*[df.divisions for df in dfs])))
    if len(divisions) == 1:  # single value for index
        divisions = (divisions[0], divisions[0])
    return divisions


def maybe_align_partitions(*exprs, divisions):
    from dask.dataframe.dask_expr._repartition import Repartition

    return [
        (
            Repartition(df, new_divisions=divisions, force=True)
            if isinstance(df, Expr) and df.ndim > 0
            else df
        )
        for df in exprs
    ]


def _extract_meta(x, nonempty=False):
    """
    Extract internal cache data (``_meta``) from dd.DataFrame / dd.Series
    """
    if isinstance(x, Expr):
        return meta_nonempty(x._meta) if nonempty else x._meta
    elif isinstance(x, list):
        return [_extract_meta(_x, nonempty) for _x in x]
    elif isinstance(x, tuple):
        return tuple(_extract_meta(_x, nonempty) for _x in x)
    elif isinstance(x, dict):
        res = {}
        for k in x:
            res[k] = _extract_meta(x[k], nonempty)
        return res
    elif hasattr(x, "expr"):
        return _extract_meta(x.expr, nonempty)
    else:
        return x


def emulate(func, *args, udf=False, **kwargs):
    """
    Apply a function using args / kwargs. If arguments contain dd.DataFrame /
    dd.Series, using internal cache (``_meta``) for calculation
    """
    with raise_on_meta_error(funcname(func), udf=udf):
        return func(*_extract_meta(args, True), **_extract_meta(kwargs, True))


def _get_meta_map_partitions(args, dfs, func, kwargs, meta, parent_meta):
    """
    Helper to generate metadata for map_partitions and map_overlap output.
    """
    meta_index = getattr(make_meta(dfs[0]), "index", None) if dfs else None
    if parent_meta is None and dfs:
        parent_meta = dfs[0]._meta
    if meta is no_default:
        # Use non-normalized kwargs here, as we want the real values (not
        # delayed values)
        a = [meta_nonempty(arg._meta) if isinstance(arg, Expr) else arg for arg in args]
        meta = emulate(func, *a, udf=True, **kwargs)
        meta_is_emulated = True
    else:
        meta = make_meta(meta, index=meta_index, parent_meta=parent_meta)
        meta_is_emulated = False

    if not (has_parallel_type(meta) or is_arraylike(meta) and meta.shape) and not all(
        isinstance(arg, Expr) and arg.ndim == 0 for arg in args
    ):
        if not meta_is_emulated:
            warnings.warn(
                "Meta is not valid, `map_partitions` and `map_overlap` expects output to be a pandas object. "
                "Try passing a pandas object as meta or a dict or tuple representing the "
                "(name, dtype) of the columns. In the future the meta you passed will not work.",
                FutureWarning,
            )
        # If `meta` is not a pandas object, the concatenated results will be a
        # different type
        meta = make_meta(_concat([meta]), index=meta_index)

    # Ensure meta is empty series
    meta = make_meta(meta, parent_meta=parent_meta)

    return meta


from dask.dataframe.dask_expr._reductions import (
    All,
    Any,
    ApplyConcatApply,
    Count,
    IdxMax,
    IdxMin,
    Max,
    Mean,
    Min,
    Mode,
    NBytes,
    NuniqueApprox,
    Prod,
    ShuffleReduce,
    Size,
    Sum,
    TreeReduce,
    Var,
)
from dask.dataframe.dask_expr.io import IO, BlockwiseIO, FromArray, FromPandas
from dask.dataframe.dask_expr.io._delayed import FromDelayed
