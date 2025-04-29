from __future__ import annotations

import datetime
import functools
import inspect
import warnings
from collections.abc import Callable, Hashable, Iterable, Mapping
from functools import wraps
from numbers import Integral, Number
from typing import Any, ClassVar, Literal

import numpy as np
import pandas as pd
import pyarrow as pa
from fsspec.utils import stringify_path
from packaging.version import parse as parse_version
from pandas import CategoricalDtype
from pandas.api.types import is_bool_dtype, is_datetime64_any_dtype, is_numeric_dtype
from pandas.api.types import is_scalar as pd_is_scalar
from pandas.api.types import is_timedelta64_dtype
from pandas.core.dtypes.common import is_extension_array_dtype
from pyarrow import fs as pa_fs
from tlz import first

import dask.array as da
import dask.dataframe.dask_expr._backends  # noqa: F401
import dask.dataframe.methods as methods
from dask import compute, get_annotations
from dask._collections import new_collection
from dask._expr import OptimizerStage
from dask.array import Array
from dask.base import DaskMethodsMixin, is_dask_collection, named_schedulers
from dask.core import flatten
from dask.dataframe._compat import PANDAS_GE_210, PANDAS_GE_220, PANDAS_VERSION
from dask.dataframe.accessor import CachedAccessor
from dask.dataframe.core import (
    _concat,
    _convert_to_numeric,
    _repr_data_series,
    _sqrt_and_convert_to_timedelta,
    check_divisions,
    has_parallel_type,
    is_dataframe_like,
    is_series_like,
    meta_warning,
)
from dask.dataframe.dask_expr import _expr as expr
from dask.dataframe.dask_expr._backends import dataframe_creation_dispatch
from dask.dataframe.dask_expr._categorical import (
    CategoricalAccessor,
    Categorize,
    GetCategories,
)
from dask.dataframe.dask_expr._concat import Concat
from dask.dataframe.dask_expr._datetime import DatetimeAccessor
from dask.dataframe.dask_expr._describe import DescribeNonNumeric, DescribeNumeric
from dask.dataframe.dask_expr._expr import (
    BFill,
    Diff,
    Eval,
    FFill,
    FillnaCheck,
    Query,
    Shift,
    ToDatetime,
    ToNumeric,
    ToTimedelta,
    no_default,
)
from dask.dataframe.dask_expr._merge import JoinRecursive, Merge
from dask.dataframe.dask_expr._quantile import SeriesQuantile
from dask.dataframe.dask_expr._quantiles import RepartitionQuantiles
from dask.dataframe.dask_expr._reductions import (
    Corr,
    Cov,
    CustomReduction,
    DropDuplicates,
    IndexCount,
    IsMonotonicDecreasing,
    IsMonotonicIncreasing,
    Len,
    MemoryUsageFrame,
    MemoryUsageIndex,
    Moment,
    NLargest,
    NSmallest,
    PivotTable,
    Unique,
    ValueCounts,
)
from dask.dataframe.dask_expr._repartition import Repartition, RepartitionFreq
from dask.dataframe.dask_expr._shuffle import (
    RearrangeByColumn,
    SetIndex,
    SetIndexBlockwise,
    SortValues,
)
from dask.dataframe.dask_expr._str_accessor import StringAccessor
from dask.dataframe.dask_expr._util import (
    PANDAS_GE_300,
    _BackendData,
    _convert_to_list,
    _get_shuffle_preferring_order,
    _is_any_real_numeric_dtype,
    _maybe_from_pandas,
    _raise_if_object_series,
    _validate_axis,
    get_specified_shuffle,
)
from dask.dataframe.dask_expr.io import FromPandasDivisions, FromScalars
from dask.dataframe.dispatch import (
    get_parallel_type,
    is_categorical_dtype,
    make_meta,
    meta_nonempty,
)
from dask.dataframe.multi import warn_dtype_mismatch
from dask.dataframe.utils import (
    AttributeNotImplementedError,
    has_known_categories,
    index_summary,
    insert_meta_param_description,
    is_scalar,
    meta_frame_constructor,
    meta_series_constructor,
    pyarrow_strings_enabled,
)
from dask.delayed import Delayed, delayed
from dask.tokenize import _tokenize_deterministic
from dask.utils import (
    IndexCallable,
    M,
    derived_from,
    get_default_shuffle_method,
    get_meta_library,
    is_arraylike,
    key_split,
    maybe_pluralize,
    memory_repr,
    put_lines,
    random_state_data,
    typename,
)
from dask.widgets import get_template

#
# Utilities to wrap Expr API
# (Helps limit boiler-plate code in collection APIs)
#


def _wrap_expr_api(*args, wrap_api=None, **kwargs):
    # Use Expr API, but convert to/from Expr objects
    assert wrap_api is not None
    result = wrap_api(
        *[arg.expr if isinstance(arg, FrameBase) else arg for arg in args],
        **kwargs,
    )
    if isinstance(result, expr.Expr):
        return new_collection(result)
    return result


def _wrap_expr_op(self, other, op=None):
    # Wrap expr operator
    assert op is not None
    if isinstance(other, FrameBase):
        other = other.expr
    elif isinstance(other, da.Array):
        other = from_dask_array(other, index=self.index, columns=self.columns)
        if self.ndim == 1 and len(self.columns):
            other = other[self.columns[0]]

    if (
        not isinstance(other, expr.Expr)
        and is_dataframe_like(other)
        or is_series_like(other)
    ):
        other = self._create_alignable_frame(other).expr

    if not isinstance(other, expr.Expr):
        return new_collection(getattr(self.expr, op)(other))
    elif expr.are_co_aligned(self.expr, other):
        return new_collection(getattr(self.expr, op)(other))
    else:
        return new_collection(expr.OpAlignPartitions(self, other, op))


def _wrap_expr_method_operator(name, class_):
    """
    Add method operators to Series or DataFrame like DataFrame.add.
    _wrap_expr_method_operator("add", DataFrame)
    """
    if class_ == DataFrame:

        def method(self, other, axis="columns", level=None, fill_value=None):
            if level is not None:
                raise NotImplementedError("level must be None")

            axis = _validate_axis(axis)

            if (
                is_dataframe_like(other)
                or is_series_like(other)
                and axis in (0, "index")
            ) and not is_dask_collection(other):
                other = self._create_alignable_frame(other)

            if axis in (1, "columns"):
                if isinstance(other, Series):
                    msg = f"Unable to {name} dd.Series with axis=1"
                    raise ValueError(msg)

            frame = self
            if isinstance(other, FrameBase) and not expr.are_co_aligned(
                self.expr, other.expr
            ):
                return new_collection(
                    expr.MethodOperatorAlign(
                        op=name,
                        frame=frame,
                        other=other,
                        axis=axis,
                        level=level,
                        fill_value=fill_value,
                    )
                )

            return new_collection(
                expr.MethodOperator(
                    name=name,
                    left=frame,
                    right=other,
                    axis=axis,
                    level=level,
                    fill_value=fill_value,
                )
            )

    elif class_ == Series:

        def method(self, other, level=None, fill_value=None, axis=0):  # type: ignore
            if level is not None:
                raise NotImplementedError("level must be None")

            axis = _validate_axis(axis)

            if is_series_like(other) and not is_dask_collection(other):
                other = self._create_alignable_frame(other)

            frame = self
            if isinstance(other, FrameBase) and not expr.are_co_aligned(
                self.expr, other.expr
            ):
                return new_collection(
                    expr.MethodOperatorAlign(
                        op=name,
                        frame=frame,
                        other=other,
                        axis=axis,
                        level=level,
                        fill_value=fill_value,
                    )
                )

            return new_collection(
                expr.MethodOperator(
                    name=name,
                    left=frame,
                    right=other,
                    axis=axis,
                    fill_value=fill_value,
                    level=level,
                )
            )

    else:
        raise NotImplementedError(f"Cannot create method operator for {class_=}")

    method.__name__ = name
    return method


def _wrap_unary_expr_op(self, op=None):
    # Wrap expr operator
    assert op is not None
    return new_collection(getattr(self.expr, op)())


_WARN_ANNOTATIONS = True
#
# Collection classes
#


class FrameBase(DaskMethodsMixin):
    """Base class for Expr-backed Collections"""

    __dask_scheduler__ = staticmethod(
        named_schedulers.get("threads", named_schedulers["sync"])
    )
    __dask_optimize__ = staticmethod(lambda dsk, keys, **kwargs: dsk)

    def __init__(self, expr):
        global _WARN_ANNOTATIONS
        if _WARN_ANNOTATIONS and (annot := get_annotations()):
            _WARN_ANNOTATIONS = False
            warnings.warn(
                f"Dask annotations {annot} detected. Annotations will be ignored when using query-planning."
            )
        self._expr = expr

    @property
    def expr(self) -> expr.Expr:
        return self._expr

    @property
    def _meta(self):
        return self.expr._meta

    @functools.cached_property
    def _meta_nonempty(self):
        return meta_nonempty(self._meta)

    @property
    def divisions(self):
        """
        Tuple of ``npartitions + 1`` values, in ascending order, marking the
        lower/upper bounds of each partition's index. Divisions allow Dask
        to know which partition will contain a given value, significantly
        speeding up operations like `loc`, `merge`, and `groupby` by not
        having to search the full dataset.

        Example: for ``divisions = (0, 10, 50, 100)``, there are three partitions,
        where the index in each partition contains values [0, 10), [10, 50),
        and [50, 100], respectively. Dask therefore knows ``df.loc[45]``
        will be in the second partition.

        When every item in ``divisions`` is ``None``, the divisions are unknown.
        Most operations can still be performed, but some will be much slower,
        and a few may fail.

        It is not supported to set ``divisions`` directly. Instead, use ``set_index``,
        which sorts and splits the data as needed.
        See https://docs.dask.org/en/latest/dataframe-design.html#partitions.
        """
        return self.expr.divisions

    @property
    def npartitions(self):
        """Return number of partitions"""
        return self.expr.npartitions

    @property
    def dtypes(self):
        """Return data types"""
        return self.expr._meta.dtypes

    @property
    def size(self):
        """Size of the Series or DataFrame as a Delayed object.

        Examples
        --------
        >>> series.size  # doctest: +SKIP
        <dask_expr.expr.Scalar: expr=df.size(), dtype=int64>
        """
        return new_collection(self.expr.size)

    @property
    def columns(self):
        return self._meta.columns

    @columns.setter
    def columns(self, columns):
        if len(columns) != len(self.columns):
            # surface pandas error
            self._expr._meta.columns = columns
        self._expr = expr.ColumnsSetter(self, columns)

    def clear_divisions(self):
        """Forget division information.

        This is useful if the divisions are no longer meaningful.
        """
        return new_collection(expr.ClearDivisions(self))

    def __len__(self):
        return new_collection(Len(self)).compute()

    @property
    def nbytes(self):
        raise NotImplementedError("nbytes is not implemented on DataFrame")

    def __reduce__(self):
        return new_collection, (self._expr,)

    def __getitem__(self, other):
        if isinstance(other, FrameBase):
            if not expr.are_co_aligned(self.expr, other.expr):
                return new_collection(expr.FilterAlign(self, other))
            return new_collection(self.expr.__getitem__(other.expr))
        elif isinstance(other, slice):
            from pandas.api.types import is_float_dtype

            is_integer_slice = any(
                isinstance(i, Integral) for i in (other.start, other.step, other.stop)
            )
            if (
                self.ndim == 2
                and is_integer_slice
                and (not is_float_dtype(self.index.dtype) or PANDAS_GE_300)
            ):
                return self.iloc[other]
            else:
                return self.loc[other]
        if isinstance(other, np.ndarray) or is_series_like(other):
            other = list(other)
        elif isinstance(other, list):
            other = other.copy()
        elif isinstance(other, np.generic):
            other = other.item()
        return new_collection(self.expr.__getitem__(other))

    def __dask_tokenize__(self):
        return type(self).__name__, self._expr._name

    def __repr__(self):
        data = self._repr_data().to_string(max_rows=5)
        _str_fmt = """Dask {klass} Structure:
{data}
Dask Name: {name}, {n_expr}
Expr={expr}"""
        if not isinstance(self, Series) and not len(self.columns):
            data = data.partition("\n")[-1].replace("Index", "Divisions")
            _str_fmt = f"Empty {_str_fmt}"
        n_expr = len({e._name for e in self.expr.walk()})

        return _str_fmt.format(
            klass=self.__class__.__name__,
            data=data,
            name=key_split(self._name),
            n_expr=maybe_pluralize(n_expr, "expression"),
            expr=self.expr,
        )

    def __bool__(self):
        raise ValueError(
            f"The truth value of a {self.__class__.__name__} is ambiguous. "
            "Use a.any() or a.all()."
        )

    def __array__(self, dtype=None, **kwargs):
        return np.array(self.compute())

    def persist(self, fuse=True, **kwargs):
        out = self.optimize(fuse=fuse)
        return DaskMethodsMixin.persist(out, **kwargs)

    def compute(self, fuse=True, concatenate=True, **kwargs):
        """Compute this DataFrame.

        This turns a lazy Dask DataFrame into an in-memory pandas DataFrame.
        The entire dataset must fit into memory before calling this operation.

        The optimizer runs over the DataFrame before triggering the computation.
        The optimizer injects a repartition operation that reduces the partition
        count to 1 to enable better optimization strategies.

        Parameters
        ----------
        fuse : bool, default True
            Whether to fuse the expression tree before computing. Fusing significantly
            reduces the number of tasks and improves performance. It shouldn't be
            disabled unless absolutely necessary.
        concatenate : bool, default True
            Whether to concatenate all partitions into a single one before computing.
            Concatenating enables more powerful optimizations but it also incurs additional
            data transfer cost. Generally, it should be enabled.
        kwargs
            Extra keywords to forward to the base compute function.

        See Also
        --------
        dask.compute
        """
        out = self
        if not isinstance(out, Scalar) and concatenate:
            out = out.repartition(npartitions=1)
        out = out.optimize(fuse=fuse)
        return DaskMethodsMixin.compute(out, **kwargs)

    def analyze(self, filename: str | None = None, format: str | None = None) -> None:
        """Outputs statistics about every node in the expression.

        analyze optimizes the expression and triggers a computation. It records statistics
        like memory usage per partition to analyze how data flow through the graph.

        .. warning::
            analyze adds plugins to the scheduler and the workers that have a non-trivial
            cost. This method should not be used in production workflows.

        Parameters
        ----------
        filename: str, None
            File to store the graph representation.
        format: str, default is png
            File format for the graph representation.

        Returns
        -------
            None, but writes a graph representation of the expression enriched with
            statistics to disk.
        """
        out = self
        if not isinstance(out, Scalar):
            out = out.repartition(npartitions=1)
        return out.expr.analyze(filename=filename, format=format)

    def explain(self, stage: OptimizerStage = "fused", format: str | None = None):
        """Create a graph representation of the Expression.

        explain runs the optimizer and creates a graph of the optimized expression
        with graphviz. No computation is triggered.

        Parameters
        ----------
        stage: {"logical", "simplified-logical", "tuned-logical", "physical", "simplified-physical", "fused"}
            The optimizer stage that is returned. Default is "fused".

            - logical: outputs the expression as is
            - simplified-logical: simplifies the expression which includes predicate
              pushdown and column projection.
            - tuned-logical: applies additional optimizations like partition squashing
            - physical: outputs the physical expression; this expression can actually
              be computed
            - simplified-physical: runs another simplification after the physical
              plan is generated
            - fused: fuses the physical expression to reduce the nodes in thr graph.

            .. warning::
                The optimizer stages are subject to change.
        format: str, default None
            The format of the output. Default is "png".

        Returns
        -------
            None, but opens a new window with the graph visualization and outputs
            a file with the graph representation.
        """
        out = self
        if not isinstance(out, Scalar):
            out = out.repartition(npartitions=1)
        return out.expr.explain(stage, format)

    def pprint(self):
        """Outputs a string representation of the DataFrame.

        The expression is returned as is. Please run optimize manually if necessary.

        Returns
        -------
            None, the representation is put into stdout.
        """
        return self.expr.pprint()

    @property
    def dask(self):
        return self.__dask_graph__()

    def __dask_graph__(self):
        out = self.expr
        out = out.lower_completely()
        return out.__dask_graph__()

    def __dask_keys__(self):
        out = self.expr
        out = out.lower_completely()
        return out.__dask_keys__()

    def simplify(self):
        return new_collection(self.expr.simplify())

    def lower_once(self):
        return new_collection(self.expr.lower_once({}))

    def optimize(self, fuse: bool = True):
        """Optimizes the DataFrame.

        Runs the optimizer with all steps over the DataFrame and wraps the result in a
        new DataFrame collection. Only use this method if you want to analyze the
        optimized expression.

        Parameters
        ----------
        fuse: bool, default True
            Whether to fuse the expression tree after running the optimizer.
            It is often easier to look at the non-fused expression when analyzing
            the result.

        Returns
        -------
            The optimized Dask Dataframe
        """
        return new_collection(self.expr.optimize(fuse=fuse))

    def __dask_postcompute__(self):
        state = new_collection(self.expr.lower_completely())
        if type(self) != type(state):
            return state.__dask_postcompute__()
        return _concat, ()

    def __dask_postpersist__(self):
        state = new_collection(self.expr.lower_completely())
        return from_graph, (
            state._meta,
            state.divisions,
            state.__dask_keys__(),
            key_split(state._name),
        )

    def __getattr__(self, key):
        try:
            # Prioritize `FrameBase` attributes
            return object.__getattribute__(self, key)
        except AttributeError as err:
            try:
                # Fall back to `expr` API
                # (Making sure to convert to/from Expr)
                val = getattr(self.expr, key)
                if callable(val):
                    return functools.partial(_wrap_expr_api, wrap_api=val)
                return val
            except AttributeError:
                # Raise original error
                raise err

    def visualize(self, tasks: bool = False, **kwargs):  # type: ignore
        """Visualize the expression or task graph

        Parameters
        ----------
        tasks:
            Whether to visualize the task graph. By default
            the expression graph will be visualized instead.
        """
        if tasks:
            return super().visualize(**kwargs)
        return self.expr.visualize(**kwargs)

    @property
    def known_divisions(self):
        """Whether the divisions are known.

        This check can be expensive if the division calculation is expensive.
        DataFrame.set_index is a good example where the calculation needs an
        inspection of the data.
        """
        return self.expr.known_divisions

    @property
    def index(self):
        """Return dask Index instance"""
        return new_collection(self.expr.index)

    @index.setter
    def index(self, value):
        assert expr.are_co_aligned(
            self.expr, value.expr
        ), "value needs to be aligned with the index"
        _expr = expr.AssignIndex(self, value)
        self._expr = _expr

    def reset_index(self, drop: bool = False):
        """Reset the index to the default index.

        Note that unlike in ``pandas``, the reset index for a Dask DataFrame will
        not be monotonically increasing from 0. Instead, it will restart at 0
        for each partition (e.g. ``index1 = [0, ..., 10], index2 = [0, ...]``).
        This is due to the inability to statically know the full length of the
        index.

        For DataFrame with multi-level index, returns a new DataFrame with
        labeling information in the columns under the index names, defaulting
        to 'level_0', 'level_1', etc. if any are None. For a standard index,
        the index name will be used (if set), otherwise a default 'index' or
        'level_0' (if 'index' is already taken) will be used.

        Parameters
        ----------
        drop : boolean, default False
            Do not try to insert index into dataframe columns.
        """
        return new_collection(expr.ResetIndex(self, drop))

    def head(self, n: int = 5, npartitions=1, compute: bool = True):
        """First n rows of the dataset

        Parameters
        ----------
        n : int, optional
            The number of rows to return. Default is 5.
        npartitions : int, optional
            Elements are only taken from the first ``npartitions``, with a
            default of 1. If there are fewer than ``n`` rows in the first
            ``npartitions`` a warning will be raised and any found rows
            returned. Pass -1 to use all partitions.
        compute : bool, optional
            Whether to compute the result, default is True.
        """
        out = new_collection(expr.Head(self, n=n, npartitions=npartitions))
        if compute:
            out = out.compute()
        return out

    def tail(self, n: int = 5, compute: bool = True):
        """Last n rows of the dataset

        Caveat, the only checks the last n rows of the last partition.
        """
        out = new_collection(expr.Tail(self, n=n))
        if compute:
            out = out.compute()
        return out

    def copy(self, deep: bool = False):
        """Make a copy of the dataframe

        This is strictly a shallow copy of the underlying computational graph.
        It does not affect the underlying data

        Parameters
        ----------
        deep : boolean, default False
            The deep value must be `False` and it is declared as a parameter just for
            compatibility with third-party libraries like cuDF and pandas
        """
        if deep is not False:
            raise ValueError(
                "The `deep` value must be False. This is strictly a shallow copy "
                "of the underlying computational graph."
            )
        return new_collection(self.expr)

    @derived_from(pd.DataFrame)
    def isin(self, values):
        if isinstance(self, DataFrame):
            # DataFrame.isin does weird alignment stuff
            bad_types = (FrameBase, pd.Series, pd.DataFrame)
        else:
            bad_types = (FrameBase,)
        if isinstance(values, bad_types):
            if (
                isinstance(values, FrameBase)
                and values.ndim == 1
                and values.npartitions == 1
            ):
                # Can broadcast
                return new_collection(expr.Isin(self, values=values))
            raise NotImplementedError("Passing a %r to `isin`" % typename(type(values)))

        # We wrap values in a delayed for two reasons:
        # - avoid serializing data in every task
        # - avoid cost of traversal of large list in optimizations
        if isinstance(values, list):
            # Motivated by https://github.com/dask/dask/issues/9411.  This appears to be
            # caused by https://github.com/dask/distributed/issues/6368, and further
            # exacerbated by the fact that the list contains duplicates.  This is a patch until
            # we can create a better fix for Serialization.
            try:
                values = list(set(values))
            except TypeError:
                pass
            if not any(is_dask_collection(v) for v in values):
                try:
                    values = np.fromiter(values, dtype=object)
                except ValueError:
                    # Numpy 1.23 supports creating arrays of iterables, while lower
                    # version 1.21.x and 1.22.x do not
                    pass

        return new_collection(
            expr.Isin(
                self,
                values=expr._DelayedExpr(
                    delayed(values, name="delayed-" + _tokenize_deterministic(values))
                ),
            )
        )

    def _partitions(self, index):
        # Used by `partitions` for partition-wise slicing

        # Convert index to list
        if isinstance(index, int):
            index = [index]
        index = np.arange(self.npartitions, dtype=object)[index].tolist()

        # Check that selection makes sense
        assert set(index).issubset(range(self.npartitions))

        # Return selected partitions
        return new_collection(expr.Partitions(self, index))

    @property
    def partitions(self):
        """Slice dataframe by partitions

        This allows partitionwise slicing of a Dask Dataframe.  You can perform normal
        Numpy-style slicing, but now rather than slice elements of the array you
        slice along partitions so, for example, ``df.partitions[:5]`` produces a new
        Dask Dataframe of the first five partitions. Valid indexers are integers, sequences
        of integers, slices, or boolean masks.

        Examples
        --------
        >>> df.partitions[0]  # doctest: +SKIP
        >>> df.partitions[:3]  # doctest: +SKIP
        >>> df.partitions[::10]  # doctest: +SKIP

        Returns
        -------
        A Dask DataFrame
        """
        return IndexCallable(self._partitions)

    def get_partition(self, n):
        """
        Get a dask DataFrame/Series representing the `nth` partition.

        Parameters
        ----------
        n : int
            The 0-indexed partition number to select.

        Returns
        -------
        Dask DataFrame or Series
            The same type as the original object.

        See Also
        --------
        DataFrame.partitions
        """
        if not 0 <= n < self.npartitions:
            msg = f"n must be 0 <= n < {self.npartitions}"
            raise ValueError(msg)
        return self.partitions[n]

    def shuffle(
        self,
        on: str | list | no_default = no_default,  # type: ignore
        ignore_index: bool = False,
        npartitions: int | None = None,
        shuffle_method: str | None = None,
        on_index: bool = False,
        **options,
    ):
        """Rearrange DataFrame into new partitions

        Uses hashing of `on` to map rows to output partitions. After this
        operation, rows with the same value of `on` will be in the same
        partition.

        Parameters
        ----------
        on : str, list of str, or Series, Index, or DataFrame
            Column names to shuffle by.
        ignore_index : optional
            Whether to ignore the index. Default is ``False``.
        npartitions : optional
            Number of output partitions. The partition count will
            be preserved by default.
        shuffle_method : optional
            Desired shuffle method. Default chosen at optimization time.
        on_index : bool, default False
            Whether to shuffle on the index. Mutually exclusive with 'on'.
            Set this to ``True`` if 'on' is not provided.
        **options : optional
            Algorithm-specific options.

        Notes
        -----
        This does not preserve a meaningful index/partitioning scheme. This
        is not deterministic if done in parallel.

        Examples
        --------
        >>> df = df.shuffle(df.columns[0])  # doctest: +SKIP
        """
        if on is no_default and not on_index:  # type: ignore
            raise TypeError(
                "Must shuffle on either columns or the index; currently shuffling on "
                "neither. Pass column(s) to 'on' or set 'on_index' to True."
            )
        elif on is not no_default and on_index:
            raise TypeError(
                "Cannot shuffle on both columns and the index. Do not pass column(s) "
                "to 'on' or set 'on_index' to False."
            )

        # Preserve partition count by default
        npartitions = npartitions or self.npartitions

        if isinstance(on, FrameBase):
            if not expr.are_co_aligned(self.expr, on.expr):
                raise TypeError(
                    "index must be aligned with the DataFrame to use as shuffle index."
                )
        else:
            if pd.api.types.is_list_like(on) and not is_dask_collection(on):
                on = list(on)
            elif isinstance(on, str) or isinstance(on, int):
                on = [on]
            elif on_index:
                on = []  # type: ignore
            bad_cols = [
                index_col
                for index_col in on
                if (index_col not in self.columns) and (index_col != self.index.name)
            ]
            if bad_cols:
                raise KeyError(
                    f"Cannot shuffle on {bad_cols}, column(s) not in dataframe to shuffle"
                )

        if (shuffle_method or get_default_shuffle_method()) == "p2p":
            from distributed.shuffle._arrow import check_dtype_support

            check_dtype_support(self._meta)

            if any(not isinstance(c, str) for c in self._meta.columns):
                unsupported = {
                    c: type(c) for c in self._meta.columns if not isinstance(c, str)
                }
                raise TypeError(
                    f"p2p requires all column names to be str, found: {unsupported}",
                )
        # Returned shuffled result
        return new_collection(
            RearrangeByColumn(
                self,
                on,
                npartitions,
                ignore_index,
                get_specified_shuffle(shuffle_method),
                options,
                index_shuffle=on_index,
            )
        )

    @derived_from(pd.DataFrame)
    def resample(self, rule, closed=None, label=None):
        from dask.dataframe.tseries.resample import Resampler

        return Resampler(self, rule, **{"closed": closed, "label": label})

    def rolling(self, window, **kwargs):
        """Provides rolling transformations.

        Parameters
        ----------
        window : int, str, offset
           Size of the moving window. This is the number of observations used
           for calculating the statistic. When not using a ``DatetimeIndex``,
           the window size must not be so large as to span more than one
           adjacent partition. If using an offset or offset alias like '5D',
           the data must have a ``DatetimeIndex``
        min_periods : int, default None
            Minimum number of observations in window required to have a value
            (otherwise result is NA).
        center : boolean, default False
            Set the labels at the center of the window.
        win_type : string, default None
            Provide a window type. The recognized window types are identical
            to pandas.
        axis : int, str, None, default 0
            This parameter is deprecated with ``pandas>=2.1``.

        Returns
        -------
        a Rolling object on which to call a method to compute a statistic
        """
        from dask.dataframe.dask_expr._rolling import Rolling

        return Rolling(self, window, **kwargs)

    @insert_meta_param_description(pad=12)
    def map_partitions(
        self,
        func,
        *args,
        meta=no_default,
        enforce_metadata=True,
        transform_divisions=True,
        clear_divisions=False,
        align_dataframes=False,
        parent_meta=None,
        **kwargs,
    ):
        """Apply a Python function to each partition

        Parameters
        ----------
        func : function
            Function applied to each partition.
        args, kwargs :
            Arguments and keywords to pass to the function. Arguments and
            keywords may contain ``FrameBase`` or regular python objects.
            DataFrame-like args (both dask and pandas) must have the same
            number of partitions as ``self`` or comprise a single partition.
            Key-word arguments, Single-partition arguments, and general
            python-object arguments will be broadcasted to all partitions.
        enforce_metadata : bool, default True
            Whether to enforce at runtime that the structure of the DataFrame
            produced by ``func`` actually matches the structure of ``meta``.
            This will rename and reorder columns for each partition, and will
            raise an error if this doesn't work, but it won't raise if dtypes
            don't match.
        transform_divisions : bool, default True
            Whether to apply the function onto the divisions and apply those
            transformed divisions to the output.
        clear_divisions : bool, default False
            Whether divisions should be cleared. If True, `transform_divisions`
            will be ignored.
        $META

        Examples
        --------
        Given a DataFrame, Series, or Index, such as:

        >>> import pandas as pd
        >>> import dask.dataframe as dd
        >>> df = pd.DataFrame({'x': [1, 2, 3, 4, 5],
        ...                    'y': [1., 2., 3., 4., 5.]})
        >>> ddf = dd.from_pandas(df, npartitions=2)

        One can use ``map_partitions`` to apply a function on each partition.
        Extra arguments and keywords can optionally be provided, and will be
        passed to the function after the partition.

        Here we apply a function with arguments and keywords to a DataFrame,
        resulting in a Series:

        >>> def myadd(df, a, b=1):
        ...     return df.x + df.y + a + b
        >>> res = ddf.map_partitions(myadd, 1, b=2)
        >>> res.dtype
        dtype('float64')

        Here we apply a function to a Series resulting in a Series:

        >>> res = ddf.x.map_partitions(lambda x: len(x)) # ddf.x is a Dask Series Structure
        >>> res.dtype
        dtype('int64')

        By default, dask tries to infer the output metadata by running your
        provided function on some fake data. This works well in many cases, but
        can sometimes be expensive, or even fail. To avoid this, you can
        manually specify the output metadata with the ``meta`` keyword. This
        can be specified in many forms, for more information see
        ``dask.dataframe.utils.make_meta``.

        Here we specify the output is a Series with no name, and dtype
        ``float64``:

        >>> res = ddf.map_partitions(myadd, 1, b=2, meta=(None, 'f8'))

        Here we map a function that takes in a DataFrame, and returns a
        DataFrame with a new column:

        >>> res = ddf.map_partitions(lambda df: df.assign(z=df.x * df.y))
        >>> res.dtypes
        x      int64
        y    float64
        z    float64
        dtype: object

        As before, the output metadata can also be specified manually. This
        time we pass in a ``dict``, as the output is a DataFrame:

        >>> res = ddf.map_partitions(lambda df: df.assign(z=df.x * df.y),
        ...                          meta={'x': 'i8', 'y': 'f8', 'z': 'f8'})

        In the case where the metadata doesn't change, you can also pass in
        the object itself directly:

        >>> res = ddf.map_partitions(lambda df: df.head(), meta=ddf)

        Also note that the index and divisions are assumed to remain unchanged.
        If the function you're mapping changes the index/divisions, you'll need
        to pass ``clear_divisions=True``.

        >>> ddf.map_partitions(func, clear_divisions=True)  # doctest: +SKIP

        Your map function gets information about where it is in the dataframe by
        accepting a special ``partition_info`` keyword argument.

        >>> def func(partition, partition_info=None):
        ...     pass

        This will receive the following information:

        >>> partition_info  # doctest: +SKIP
        {'number': 1, 'division': 3}

        For each argument and keyword arguments that are dask dataframes you will
        receive the number (n) which represents the nth partition of the dataframe
        and the division (the first index value in the partition). If divisions
        are not known (for instance if the index is not sorted) then you will get
        None as the division.
        """
        return map_partitions(
            func,
            self,
            *args,
            meta=meta,
            enforce_metadata=enforce_metadata,
            transform_divisions=transform_divisions,
            clear_divisions=clear_divisions,
            align_dataframes=align_dataframes,
            parent_meta=parent_meta,
            **kwargs,
        )

    @insert_meta_param_description(pad=12)
    def map_overlap(
        self,
        func,
        before,
        after,
        *args,
        meta=no_default,
        enforce_metadata=True,
        transform_divisions=True,
        clear_divisions=False,
        align_dataframes=False,
        **kwargs,
    ):
        """Apply a function to each partition, sharing rows with adjacent partitions.

        This can be useful for implementing windowing functions such as
        ``df.rolling(...).mean()`` or ``df.diff()``.

        Parameters
        ----------
        func : function
            Function applied to each partition.
        before : int, timedelta or string timedelta
            The rows to prepend to partition ``i`` from the end of
            partition ``i - 1``.
        after : int, timedelta or string timedelta
            The rows to append to partition ``i`` from the beginning
            of partition ``i + 1``.
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

        Notes
        -----
        Given positive integers ``before`` and ``after``, and a function
        ``func``, ``map_overlap`` does the following:

        1. Prepend ``before`` rows to each partition ``i`` from the end of
           partition ``i - 1``. The first partition has no rows prepended.

        2. Append ``after`` rows to each partition ``i`` from the beginning of
           partition ``i + 1``. The last partition has no rows appended.

        3. Apply ``func`` to each partition, passing in any extra ``args`` and
           ``kwargs`` if provided.

        4. Trim ``before`` rows from the beginning of all but the first
           partition.

        5. Trim ``after`` rows from the end of all but the last partition.

        Examples
        --------
        Given a DataFrame, Series, or Index, such as:

        >>> import pandas as pd
        >>> import dask.dataframe as dd
        >>> df = pd.DataFrame({'x': [1, 2, 4, 7, 11],
        ...                    'y': [1., 2., 3., 4., 5.]})
        >>> ddf = dd.from_pandas(df, npartitions=2)

        A rolling sum with a trailing moving window of size 2 can be computed by
        overlapping 2 rows before each partition, and then mapping calls to
        ``df.rolling(2).sum()``:

        >>> ddf.compute()
            x    y
        0   1  1.0
        1   2  2.0
        2   4  3.0
        3   7  4.0
        4  11  5.0
        >>> ddf.map_overlap(lambda df: df.rolling(2).sum(), 2, 0).compute()
              x    y
        0   NaN  NaN
        1   3.0  3.0
        2   6.0  5.0
        3  11.0  7.0
        4  18.0  9.0

        The pandas ``diff`` method computes a discrete difference shifted by a
        number of periods (can be positive or negative). This can be
        implemented by mapping calls to ``df.diff`` to each partition after
        prepending/appending that many rows, depending on sign:

        >>> def diff(df, periods=1):
        ...     before, after = (periods, 0) if periods > 0 else (0, -periods)
        ...     return df.map_overlap(lambda df, periods=1: df.diff(periods),
        ...                           periods, 0, periods=periods)
        >>> diff(ddf, 1).compute()
             x    y
        0  NaN  NaN
        1  1.0  1.0
        2  2.0  1.0
        3  3.0  1.0
        4  4.0  1.0

        If you have a ``DatetimeIndex``, you can use a ``pd.Timedelta`` for time-
        based windows or any ``pd.Timedelta`` convertible string:

        >>> ts = pd.Series(range(10), index=pd.date_range('2017', periods=10))
        >>> dts = dd.from_pandas(ts, npartitions=2)
        >>> dts.map_overlap(lambda df: df.rolling('2D').sum(),
        ...                 pd.Timedelta('2D'), 0).compute()
        2017-01-01     0.0
        2017-01-02     1.0
        2017-01-03     3.0
        2017-01-04     5.0
        2017-01-05     7.0
        2017-01-06     9.0
        2017-01-07    11.0
        2017-01-08    13.0
        2017-01-09    15.0
        2017-01-10    17.0
        Freq: D, dtype: float64
        """
        return map_overlap(
            func,
            self,
            before,
            after,
            *args,
            meta=meta,
            enforce_metadata=enforce_metadata,
            transform_divisions=transform_divisions,
            clear_divisions=clear_divisions,
            align_dataframes=align_dataframes,
            **kwargs,
        )

    def repartition(
        self,
        divisions: tuple | None = None,
        npartitions: int | None = None,
        partition_size: str | None = None,
        freq=None,
        force: bool = False,
    ):
        """Repartition a collection

        Exactly one of `divisions`, `npartitions` or `partition_size` should be
        specified. A ``ValueError`` will be raised when that is not the case.

        Parameters
        ----------
        divisions : list, optional
            The "dividing lines" used to split the dataframe into partitions.
            For ``divisions=[0, 10, 50, 100]``, there would be three output partitions,
            where the new index contained [0, 10), [10, 50), and [50, 100), respectively.
            See https://docs.dask.org/en/latest/dataframe-design.html#partitions.
        npartitions : int, Callable, optional
            Approximate number of partitions of output. The number of
            partitions used may be slightly lower than npartitions depending
            on data distribution, but will never be higher.
            The Callable gets the number of partitions of the input as an argument
            and should return an int.
        partition_size : str, optional
            Max number of bytes of memory for each partition. Use numbers or strings
            like 5MB. If specified npartitions and divisions will be ignored. Note that
            the size reflects the number of bytes used as computed by
            pandas.DataFrame.memory_usage, which will not necessarily match the size
            when storing to disk.

            .. warning::

               This keyword argument triggers computation to determine
               the memory size of each partition, which may be expensive.

        force : bool, default False
            Allows the expansion of the existing divisions.
            If False then the new divisions' lower and upper bounds must be
            the same as the old divisions'.
        freq : str, pd.Timedelta
            A period on which to partition timeseries data like ``'7D'`` or
            ``'12h'`` or ``pd.Timedelta(hours=12)``.  Assumes a datetime index.

        Notes
        -----
        Exactly one of `divisions`, `npartitions`, `partition_size`, or `freq`
        should be specified. A ``ValueError`` will be raised when that is
        not the case.

        Also note that ``len(divisons)`` is equal to ``npartitions + 1``. This is because ``divisions``
        represents the upper and lower bounds of each partition. The first item is the
        lower bound of the first partition, the second item is the lower bound of the
        second partition and the upper bound of the first partition, and so on.
        The second-to-last item is the lower bound of the last partition, and the last
        (extra) item is the upper bound of the last partition.

        Examples
        --------
        >>> df = df.repartition(npartitions=10)  # doctest: +SKIP
        >>> df = df.repartition(divisions=[0, 5, 10, 20])  # doctest: +SKIP
        >>> df = df.repartition(freq='7d')  # doctest: +SKIP

        See Also
        --------
        DataFrame.memory_usage_per_partition
        pandas.DataFrame.memory_usage
        """

        if (
            sum(
                [
                    divisions is not None,
                    npartitions is not None,
                    partition_size is not None,
                    freq is not None,
                ]
            )
            != 1
        ):
            raise ValueError(
                "Please provide exactly one of the ``npartitions=`` or "
                "``divisions=`` keyword arguments."
            )
        if divisions is not None:
            check_divisions(divisions)
        if freq is not None:
            if not isinstance(self.divisions[0], pd.Timestamp):
                raise TypeError("Can only repartition on frequency for timeseries")
            return new_collection(RepartitionFreq(self, freq))
        else:
            return new_collection(
                Repartition(self, npartitions, divisions, force, partition_size, freq)
            )

    def to_dask_array(
        self, lengths=None, meta=None, optimize: bool = True, **optimize_kwargs
    ) -> Array:
        """Convert a dask DataFrame to a dask array.

        Parameters
        ----------
        lengths : bool or Sequence of ints, optional
            How to determine the chunks sizes for the output array.
            By default, the output array will have unknown chunk lengths
            along the first axis, which can cause some later operations
            to fail.

            * True : immediately compute the length of each partition
            * Sequence : a sequence of integers to use for the chunk sizes
              on the first axis. These values are *not* validated for
              correctness, beyond ensuring that the number of items
              matches the number of partitions.
        meta : object, optional
            An optional `meta` parameter can be passed for dask to override the
            default metadata on the underlying dask array.
        optimize : bool
            Whether to optimize the expression before converting to an Array.

        Returns
        -------
        A Dask Array
        """
        if lengths is True:
            lengths = tuple(self.map_partitions(len, enforce_metadata=False).compute())

        arr = self.values

        chunks = self._validate_chunks(arr, lengths)
        arr._chunks = chunks

        if meta is not None:
            arr._meta = meta

        return arr

    @property
    def values(self):
        """Return a dask.array of the values of this dataframe

        Warning: This creates a dask.array without precise shape information.
        Operations that depend on shape information, like slicing or reshaping,
        will not work.
        """
        if is_extension_array_dtype(self._meta.values):
            warnings.warn(
                "Dask currently has limited support for converting pandas extension dtypes "
                f"to arrays. Converting {self._meta.values.dtype} to object dtype.",
                UserWarning,
            )
        return self.map_partitions(methods.values)

    def __divmod__(self, other):
        result = self.expr.__divmod__(other)
        return new_collection(result[0]), new_collection(result[1])

    def __rdivmod__(self, other):
        result = self.expr.__rdivmod__(other)
        return new_collection(result[0]), new_collection(result[1])

    def __abs__(self):
        return self.abs()

    @derived_from(pd.DataFrame)
    def sum(
        self,
        axis=0,
        skipna=True,
        numeric_only=False,
        min_count=0,
        split_every=False,
        **kwargs,
    ):
        axis = self._validate_axis(axis)
        if axis == 1:
            return self.map_partitions(
                M.sum,
                skipna=skipna,
                numeric_only=numeric_only,
                axis=axis,
                min_count=min_count,
            )

        result = new_collection(self.expr.sum(skipna, numeric_only, split_every, axis))
        return self._apply_min_count(result, min_count)

    def _apply_min_count(self, result, min_count):
        if min_count:
            cond = self.notnull().sum() >= min_count
            cond_meta = cond._meta
            if not is_series_like(cond_meta):
                result = result.to_series()
                cond = cond.to_series()

            result = result.where(cond, other=np.nan)
            if not is_series_like(cond_meta):
                return result.min()
            else:
                return result
        else:
            return result

    @derived_from(pd.DataFrame)
    def prod(
        self,
        axis=0,
        skipna=True,
        numeric_only=False,
        min_count=0,
        split_every=False,
        **kwargs,
    ):
        axis = self._validate_axis(axis)
        if axis == 1:
            return self.map_partitions(
                M.prod,
                skipna=skipna,
                numeric_only=numeric_only,
                axis=axis,
                min_count=min_count,
            )
        result = new_collection(self.expr.prod(skipna, numeric_only, split_every, axis))
        return self._apply_min_count(result, min_count)

    product = prod

    @derived_from(pd.DataFrame)
    def var(
        self,
        axis=0,
        skipna=True,
        ddof=1,
        numeric_only=False,
        split_every=False,
        **kwargs,
    ):
        _raise_if_object_series(self, "var")
        axis = self._validate_axis(axis)
        self._meta.var(axis=axis, skipna=skipna, numeric_only=numeric_only)
        frame = self
        if is_dataframe_like(self._meta) and numeric_only:
            frame = frame[list(self._meta.var(numeric_only=True).index)]
        return new_collection(
            frame.expr.var(axis, skipna, ddof, numeric_only, split_every=split_every)
        )

    @derived_from(pd.DataFrame)
    def std(
        self,
        axis=0,
        skipna=True,
        ddof=1,
        numeric_only=False,
        split_every=False,
        **kwargs,
    ):
        _raise_if_object_series(self, "std")
        axis = self._validate_axis(axis)
        numeric_dd = self
        meta = meta_nonempty(self._meta).std(
            axis=axis, skipna=skipna, ddof=ddof, numeric_only=numeric_only
        )
        needs_time_conversion, time_cols = False, None
        if is_dataframe_like(self._meta):
            if axis == 0:
                numeric_dd = numeric_dd[list(meta.index)]
            else:
                numeric_dd = numeric_dd.copy()

            if numeric_only is True:
                _meta = numeric_dd._meta.select_dtypes(include=[np.number])
            else:
                _meta = numeric_dd._meta
            time_cols = _meta.select_dtypes(include=["datetime", "timedelta"]).columns
            if len(time_cols) > 0:
                if axis == 1 and len(time_cols) != len(self.columns):
                    numeric_dd = from_pandas(
                        meta_frame_constructor(self)(
                            {"_": meta_series_constructor(self)([np.nan])},
                            index=self.index,
                        ),
                        npartitions=self.npartitions,
                    )
                else:
                    needs_time_conversion = True
                    if axis == 1:
                        numeric_dd = numeric_dd.astype(f"datetime64[{meta.array.unit}]")
                    for col in time_cols:
                        numeric_dd[col] = _convert_to_numeric(numeric_dd[col], skipna)
        else:
            needs_time_conversion = is_datetime64_any_dtype(self._meta)
            if needs_time_conversion:
                numeric_dd = _convert_to_numeric(self, skipna)

        units = None
        if needs_time_conversion and time_cols is not None:
            units = [getattr(self._meta[c].array, "unit", None) for c in time_cols]

        if axis == 1:
            _kwargs = (
                {}
                if not needs_time_conversion
                else {"unit": meta.array.unit, "dtype": meta.dtype}
            )
            return numeric_dd.map_partitions(
                M.std if not needs_time_conversion else _sqrt_and_convert_to_timedelta,
                meta=meta,
                axis=axis,
                skipna=skipna,
                ddof=ddof,
                enforce_metadata=False,
                numeric_only=numeric_only,
                **_kwargs,
            )

        result = numeric_dd.var(
            skipna=skipna, ddof=ddof, numeric_only=numeric_only, split_every=split_every
        )

        if needs_time_conversion:
            sqrt_func_kwargs = {
                "is_df_like": is_dataframe_like(self._meta),
                "time_cols": time_cols,
                "axis": axis,
                "dtype": getattr(meta, "dtype", None),
                "unit": getattr(meta, "unit", None),
                "units": units,
            }
            sqrt_func = _sqrt_and_convert_to_timedelta
        else:
            sqrt_func_kwargs = {}
            sqrt_func = np.sqrt

        result = result.map_partitions(
            sqrt_func,
            meta=meta,
            enforce_metadata=False,
            parent_meta=self._meta,
            **sqrt_func_kwargs,
        )
        return result

    def enforce_runtime_divisions(self):
        """Enforce the current divisions at runtime.

        Injects a layer into the Task Graph that checks that the current divisions
        match the expected divisions at runtime.
        """
        if not self.known_divisions:
            raise ValueError("No known divisions to enforce!")
        return new_collection(expr.EnforceRuntimeDivisions(self))

    @derived_from(pd.DataFrame)
    def skew(
        self,
        axis=0,
        bias=True,
        nan_policy="propagate",
        numeric_only=False,
    ):
        """
        .. note::

           This implementation follows the dask.array.stats implementation
           of skewness and calculates skewness without taking into account
           a bias term for finite sample size, which corresponds to the
           default settings of the scipy.stats skewness calculation. However,
           Pandas corrects for this, so the values differ by a factor of
           (n * (n - 1)) ** 0.5 / (n - 2), where n is the number of samples.

           Further, this method currently does not support filtering out NaN
           values, which is again a difference to Pandas.
        """
        _raise_if_object_series(self, "skew")
        if axis is None:
            raise ValueError("`axis=None` isn't currently supported for `skew`")
        axis = self._validate_axis(axis)

        if is_dataframe_like(self):
            # Let pandas raise errors if necessary
            meta = self._meta_nonempty.skew(axis=axis, numeric_only=numeric_only)
        else:
            meta = self._meta_nonempty.skew()

        if axis == 1:
            return self.map_partitions(
                M.skew,
                meta=meta,
                axis=axis,
                enforce_metadata=False,
            )

        if not bias:
            raise NotImplementedError("bias=False is not implemented.")
        if nan_policy != "propagate":
            raise NotImplementedError(
                "`nan_policy` other than 'propagate' have not been implemented."
            )

        frame = self
        if frame.ndim > 1:
            frame = frame.select_dtypes(
                include=["number", "bool"], exclude=[np.timedelta64]
            )
        m2 = new_collection(Moment(frame, order=2))
        m3 = new_collection(Moment(frame, order=3))
        result = m3 / m2**1.5
        if result.ndim == 1:
            result = result.fillna(0.0)
        return result

    @derived_from(pd.DataFrame)
    def kurtosis(
        self,
        axis=0,
        fisher=True,
        bias=True,
        nan_policy="propagate",
        numeric_only=False,
    ):
        """
        .. note::

           This implementation follows the dask.array.stats implementation
           of kurtosis and calculates kurtosis without taking into account
           a bias term for finite sample size, which corresponds to the
           default settings of the scipy.stats kurtosis calculation. This differs
           from pandas.

           Further, this method currently does not support filtering out NaN
           values, which is again a difference to Pandas.
        """
        _raise_if_object_series(self, "kurtosis")
        if axis is None:
            raise ValueError("`axis=None` isn't currently supported for `skew`")
        axis = self._validate_axis(axis)

        if is_dataframe_like(self):
            # Let pandas raise errors if necessary
            meta = self._meta_nonempty.kurtosis(axis=axis, numeric_only=numeric_only)
        else:
            meta = self._meta_nonempty.kurtosis()

        if axis == 1:
            return map_partitions(
                M.kurtosis,
                self,
                meta=meta,
                token=self._token_prefix + "kurtosis",
                axis=axis,
                enforce_metadata=False,
            )

        if not bias:
            raise NotImplementedError("bias=False is not implemented.")
        if nan_policy != "propagate":
            raise NotImplementedError(
                "`nan_policy` other than 'propagate' have not been implemented."
            )

        frame = self
        if frame.ndim > 1:
            frame = frame.select_dtypes(
                include=["number", "bool"], exclude=[np.timedelta64]
            )
        m2 = new_collection(Moment(frame, order=2))
        m4 = new_collection(Moment(frame, order=4))
        result = m4 / m2**2.0
        if result.ndim == 1:
            result = result.fillna(0.0)
        if fisher:
            return result - 3
        else:
            return result

    kurt = kurtosis

    @derived_from(pd.DataFrame)
    def sem(
        self, axis=None, skipna=True, ddof=1, split_every=False, numeric_only=False
    ):
        axis = self._validate_axis(axis)
        _raise_if_object_series(self, "sem")
        if axis == 1:
            return self.map_partitions(
                M.sem,
                axis=axis,
                skipna=skipna,
                ddof=ddof,
                numeric_only=numeric_only,
            )
        meta = self._meta.sem(skipna=skipna, ddof=ddof, numeric_only=numeric_only)
        frame = self
        if self.ndim == 2:
            frame = self[list(meta.index)]

        v = frame.var(skipna=skipna, ddof=ddof, split_every=split_every)
        n = frame.count(split_every=split_every)
        result = map_partitions(
            np.sqrt,
            v / n,
            meta=meta,
            enforce_metadata=False,
            parent_meta=self._meta,
        )
        return result

    def _prepare_cov_corr(self, min_periods, numeric_only):
        if min_periods is None:
            min_periods = 2
        elif min_periods < 2:
            raise ValueError("min_periods must be >= 2")

        self._meta.cov(numeric_only=numeric_only, min_periods=min_periods)

        frame = self
        if numeric_only:
            numerics = self._meta._get_numeric_data()
            if len(numerics.columns) != len(self.columns):
                frame = frame[list(numerics.columns)]
        return frame, min_periods

    def _cov(
        self, min_periods=None, numeric_only=False, split_every=False, scalar=False
    ):
        frame, min_periods = self._prepare_cov_corr(min_periods, numeric_only)
        return new_collection(Cov(frame, min_periods, split_every, scalar))

    def _corr(
        self,
        method="pearson",
        min_periods=None,
        numeric_only=False,
        split_every=False,
        scalar=False,
    ):
        if method != "pearson":
            raise NotImplementedError("Only Pearson correlation has been implemented")
        frame, min_periods = self._prepare_cov_corr(min_periods, numeric_only)
        return new_collection(Corr(frame, min_periods, split_every, scalar))

    @derived_from(pd.DataFrame)
    def mean(
        self, axis=0, skipna=True, numeric_only=False, split_every=False, **kwargs
    ):
        _raise_if_object_series(self, "mean")
        axis = self._validate_axis(axis)
        if axis == 1:
            return self.map_partitions(
                M.mean, skipna=skipna, numeric_only=numeric_only, axis=axis
            )
        return new_collection(
            self.expr.mean(skipna, numeric_only, split_every=split_every, axis=axis)
        )

    @derived_from(pd.DataFrame)
    def max(self, axis=0, skipna=True, numeric_only=False, split_every=False, **kwargs):
        axis = self._validate_axis(axis)
        if axis == 1:
            return self.map_partitions(
                M.max, skipna=skipna, numeric_only=numeric_only, axis=axis
            )
        return new_collection(self.expr.max(skipna, numeric_only, split_every, axis))

    @derived_from(pd.DataFrame)
    def any(self, axis=0, skipna=True, split_every=False, **kwargs):
        axis = self._validate_axis(axis)
        if axis == 1:
            return self.map_partitions(M.any, skipna=skipna, axis=axis)
        return new_collection(self.expr.any(skipna, split_every))

    @derived_from(pd.DataFrame)
    def all(self, axis=0, skipna=True, split_every=False, **kwargs):
        axis = self._validate_axis(axis)
        if axis == 1:
            return self.map_partitions(M.all, skipna=skipna, axis=axis)
        return new_collection(self.expr.all(skipna, split_every))

    @derived_from(pd.DataFrame)
    def idxmin(self, axis=0, skipna=True, numeric_only=False, split_every=False):
        axis = self._validate_axis(axis)
        if axis == 1:
            return self.map_partitions(
                M.idxmin, skipna=skipna, numeric_only=numeric_only, axis=axis
            )
        return new_collection(self.expr.idxmin(skipna, numeric_only, split_every))

    @derived_from(pd.DataFrame)
    def idxmax(self, axis=0, skipna=True, numeric_only=False, split_every=False):
        axis = self._validate_axis(axis)
        if axis == 1:
            return self.map_partitions(
                M.idxmax, skipna=skipna, numeric_only=numeric_only, axis=axis
            )
        return new_collection(self.expr.idxmax(skipna, numeric_only, split_every))

    @derived_from(pd.DataFrame)
    def min(self, axis=0, skipna=True, numeric_only=False, split_every=False, **kwargs):
        axis = self._validate_axis(axis)
        if axis == 1:
            return self.map_partitions(
                M.min, skipna=skipna, numeric_only=numeric_only, axis=axis
            )
        return new_collection(self.expr.min(skipna, numeric_only, split_every, axis))

    @derived_from(pd.DataFrame)
    def count(self, axis=0, numeric_only=False, split_every=False):
        axis = self._validate_axis(axis)
        if axis == 1:
            return self.map_partitions(M.count, numeric_only=numeric_only, axis=axis)
        return new_collection(self.expr.count(numeric_only, split_every))

    @derived_from(pd.DataFrame)
    def abs(self):
        # Raise pandas errors
        _raise_if_object_series(self, "abs")
        meta_nonempty(self._meta).abs()
        return new_collection(self.expr.abs())

    @derived_from(pd.DataFrame)
    def astype(self, dtypes):
        return new_collection(self.expr.astype(dtypes))

    @derived_from(pd.DataFrame)
    def combine_first(self, other):
        other = self._create_alignable_frame(other, "outer").expr
        return new_collection(self.expr.combine_first(other))

    @derived_from(pd.DataFrame)
    def to_timestamp(self, freq=None, how="start"):
        return new_collection(self.expr.to_timestamp(freq, how))

    @derived_from(pd.DataFrame)
    def isna(self):
        return new_collection(self.expr.isna())

    def random_split(self, frac, random_state=None, shuffle=False):
        """Pseudorandomly split dataframe into different pieces row-wise

        Parameters
        ----------
        frac : list
            List of floats that should sum to one.
        random_state : int or np.random.RandomState
            If int or None create a new RandomState with this as the seed.
            Otherwise draw from the passed RandomState.
        shuffle : bool, default False
            If set to True, the dataframe is shuffled (within partition)
            before the split.

        Examples
        --------

        50/50 split

        >>> a, b = df.random_split([0.5, 0.5])  # doctest: +SKIP

        80/10/10 split, consistent random_state

        >>> a, b, c = df.random_split([0.8, 0.1, 0.1], random_state=123)  # doctest: +SKIP

        See Also
        --------
        dask.DataFrame.sample
        """
        if not np.allclose(sum(frac), 1):
            raise ValueError("frac should sum to 1")
        frame = expr.Split(self, frac, random_state, shuffle)

        out = []
        for i in range(len(frac)):
            out.append(new_collection(expr.SplitTake(frame, i, self.ndim)))
        return out

    @derived_from(pd.DataFrame)
    def round(self, decimals=0):
        return new_collection(self.expr.round(decimals))

    @derived_from(pd.DataFrame)
    def where(self, cond, other=np.nan):
        cond = self._create_alignable_frame(cond)
        other = self._create_alignable_frame(other)
        cond = cond.expr if isinstance(cond, FrameBase) else cond
        other = other.expr if isinstance(other, FrameBase) else other
        return new_collection(self.expr.where(cond, other))

    @derived_from(pd.DataFrame)
    def mask(self, cond, other=np.nan):
        cond = self._create_alignable_frame(cond)
        other = self._create_alignable_frame(other)
        cond = cond.expr if isinstance(cond, FrameBase) else cond
        other = other.expr if isinstance(other, FrameBase) else other
        return new_collection(self.expr.mask(cond, other))

    @derived_from(pd.DataFrame)
    def replace(self, to_replace=None, value=no_default, regex=False):
        return new_collection(self.expr.replace(to_replace, value, regex))

    @derived_from(pd.DataFrame)
    def ffill(self, axis=0, limit=None):
        axis = _validate_axis(axis)
        if axis == 1:
            return self.map_partitions(M.ffill, axis=axis, limit=limit)
        frame = self
        if limit is None:
            frame = FillnaCheck(self, "ffill", lambda x: 0)
        return new_collection(FFill(frame, limit))

    @derived_from(pd.DataFrame)
    def bfill(self, axis=0, limit=None):
        axis = _validate_axis(axis)
        if axis == 1:
            return self.map_partitions(M.bfill, axis=axis, limit=limit)
        frame = self
        if limit is None:
            frame = FillnaCheck(self, "bfill", lambda x: x.npartitions - 1)
        return new_collection(BFill(frame, limit))

    @derived_from(pd.DataFrame)
    def fillna(self, value=None, axis=None):
        axis = self._validate_axis(axis)
        if axis == 1:
            return self.map_partitions(M.fillna, value, axis=axis)
        if isinstance(value, FrameBase):
            value = value.expr
        return new_collection(self.expr.fillna(value))

    @derived_from(pd.DataFrame)
    def shift(self, periods=1, freq=None, axis=0):
        if not isinstance(periods, Integral):
            raise TypeError("periods must be an integer")

        axis = _validate_axis(axis)
        if axis == 0:
            return new_collection(Shift(self, periods, freq))

        return self.map_partitions(
            func=Shift.func,
            enforce_metadata=False,
            transform_divisions=False,
            periods=periods,
            axis=axis,
            freq=freq,
        )

    @derived_from(pd.DataFrame)
    def diff(self, periods=1, axis=0):
        """
        .. note::

           Pandas currently uses an ``object``-dtype column to represent
           boolean data with missing values. This can cause issues for
           boolean-specific operations, like ``|``. To enable boolean-
           specific operations, at the cost of metadata that doesn't match
           pandas, use ``.astype(bool)`` after the ``shift``.
        """
        axis = _validate_axis(axis)
        if axis == 0:
            return new_collection(Diff(self, periods))
        return self.map_partitions(
            func=Diff.func,
            enforce_metadata=False,
            transform_divisions=False,
            clear_divisions=False,
            periods=periods,
            axis=axis,
        )

    @derived_from(pd.DataFrame)
    def rename_axis(
        self, mapper=no_default, index=no_default, columns=no_default, axis=0
    ):
        return new_collection(self.expr.rename_axis(mapper, index, columns, axis))

    def _create_alignable_frame(self, other, join="outer"):
        if not is_dask_collection(other) and (
            is_series_like(other) or is_dataframe_like(other)
        ):
            if join in ("inner", "left"):
                npartitions = 1
            else:
                # We have to trigger alignment, otherwise pandas will add
                # the same values to every partition
                npartitions = 2
            other = from_pandas(other, npartitions=npartitions)
        return other

    @derived_from(pd.DataFrame)
    def align(self, other, join="outer", axis=None, fill_value=None):
        other = self._create_alignable_frame(other, join)
        return self.expr.align(other.expr, join, axis, fill_value)

    def nunique_approx(self, split_every=None):
        """Approximate number of unique rows.

        This method uses the HyperLogLog algorithm for cardinality
        estimation to compute the approximate number of unique rows.
        The approximate error is 0.406%.

        Parameters
        ----------
        split_every : int, optional
            Group partitions into groups of this size while performing a
            tree-reduction. If set to False, no tree-reduction will be used.
            Default is 8.

        Returns
        -------
        a float representing the approximate number of elements
        """
        return new_collection(self.expr.nunique_approx(split_every=split_every))

    @derived_from(pd.DataFrame)
    def cumsum(self, axis=0, skipna=True, **kwargs):
        if axis == 1:
            return self.map_partitions(M.cumsum, axis=axis, skipna=skipna)
        return new_collection(self.expr.cumsum(skipna=skipna))

    @derived_from(pd.DataFrame)
    def cumprod(self, axis=0, skipna=True, **kwargs):
        if axis == 1:
            return self.map_partitions(M.cumprod, axis=axis, skipna=skipna)
        return new_collection(self.expr.cumprod(skipna=skipna))

    @derived_from(pd.DataFrame)
    def cummax(self, axis=0, skipna=True):
        if axis == 1:
            return self.map_partitions(M.cummax, axis=axis, skipna=skipna)
        return new_collection(self.expr.cummax(skipna=skipna))

    @derived_from(pd.DataFrame)
    def cummin(self, axis=0, skipna=True):
        if axis == 1:
            return self.map_partitions(M.cummin, axis=axis, skipna=skipna)
        return new_collection(self.expr.cummin(skipna=skipna))

    def reduction(
        self,
        chunk,
        aggregate=None,
        combine=None,
        meta=no_default,
        token=None,
        split_every=None,
        chunk_kwargs=None,
        aggregate_kwargs=None,
        combine_kwargs=None,
        **kwargs,
    ):
        """Generic row-wise reductions.

        Parameters
        ----------
        chunk : callable
            Function to operate on each partition. Should return a
            ``pandas.DataFrame``, ``pandas.Series``, or a scalar.
        aggregate : callable, optional
            Function to operate on the concatenated result of ``chunk``. If not
            specified, defaults to ``chunk``. Used to do the final aggregation
            in a tree reduction.

            The input to ``aggregate`` depends on the output of ``chunk``.
            If the output of ``chunk`` is a:

            - scalar: Input is a Series, with one row per partition.
            - Series: Input is a DataFrame, with one row per partition. Columns
              are the rows in the output series.
            - DataFrame: Input is a DataFrame, with one row per partition.
              Columns are the columns in the output dataframes.

            Should return a ``pandas.DataFrame``, ``pandas.Series``, or a
            scalar.
        combine : callable, optional
            Function to operate on intermediate concatenated results of
            ``chunk`` in a tree-reduction. If not provided, defaults to
            ``aggregate``. The input/output requirements should match that of
            ``aggregate`` described above.
        $META
        token : str, optional
            The name to use for the output keys.
        split_every : int, optional
            Group partitions into groups of this size while performing a
            tree-reduction. If set to False, no tree-reduction will be used,
            and all intermediates will be concatenated and passed to
            ``aggregate``. Default is 8.
        chunk_kwargs : dict, optional
            Keyword arguments to pass on to ``chunk`` only.
        aggregate_kwargs : dict, optional
            Keyword arguments to pass on to ``aggregate`` only.
        combine_kwargs : dict, optional
            Keyword arguments to pass on to ``combine`` only.
        kwargs :
            All remaining keywords will be passed to ``chunk``, ``combine``,
            and ``aggregate``.

        Examples
        --------
        >>> import pandas as pd
        >>> import dask.dataframe as dd
        >>> df = pd.DataFrame({'x': range(50), 'y': range(50, 100)})
        >>> ddf = dd.from_pandas(df, npartitions=4)

        Count the number of rows in a DataFrame. To do this, count the number
        of rows in each partition, then sum the results:

        >>> res = ddf.reduction(lambda x: x.count(),
        ...                     aggregate=lambda x: x.sum())
        >>> res.compute()
        x    50
        y    50
        dtype: int64

        Count the number of rows in a Series with elements greater than or
        equal to a value (provided via a keyword).

        >>> def count_greater(x, value=0):
        ...     return (x >= value).sum()
        >>> res = ddf.x.reduction(count_greater, aggregate=lambda x: x.sum(),
        ...                       chunk_kwargs={'value': 25})
        >>> res.compute()
        np.int64(25)

        Aggregate both the sum and count of a Series at the same time:

        >>> def sum_and_count(x):
        ...     return pd.Series({'count': x.count(), 'sum': x.sum()},
        ...                      index=['count', 'sum'])
        >>> res = ddf.x.reduction(sum_and_count, aggregate=lambda x: x.sum())
        >>> res.compute()
        count      50
        sum      1225
        dtype: int64

        Doing the same, but for a DataFrame. Here ``chunk`` returns a
        DataFrame, meaning the input to ``aggregate`` is a DataFrame with an
        index with non-unique entries for both 'x' and 'y'. We groupby the
        index, and sum each group to get the final result.

        >>> def sum_and_count(x):
        ...     return pd.DataFrame({'count': x.count(), 'sum': x.sum()},
        ...                         columns=['count', 'sum'])
        >>> res = ddf.reduction(sum_and_count,
        ...                     aggregate=lambda x: x.groupby(level=0).sum())
        >>> res.compute()
           count   sum
        x     50  1225
        y     50  3725
        """
        if split_every is not None and split_every < 2 and split_every is not False:
            raise ValueError("split_every must be at least 2")

        if combine is None:
            if combine_kwargs:
                raise ValueError("`combine_kwargs` provided with no `combine`")

        chunk_kwargs = chunk_kwargs.copy() if chunk_kwargs else {}
        chunk_kwargs.update(kwargs)
        chunk_kwargs["func"] = chunk

        combine_kwargs = combine_kwargs.copy() if combine_kwargs else {}
        combine_kwargs.update(kwargs)
        combine_kwargs["func"] = combine or aggregate or chunk

        aggregate_kwargs = aggregate_kwargs.copy() if aggregate_kwargs else {}
        aggregate_kwargs.update(kwargs)
        aggregate_kwargs["func"] = aggregate or chunk

        return new_collection(
            CustomReduction(
                self,
                meta,
                chunk_kwargs,
                aggregate_kwargs,
                combine_kwargs,
                split_every,
                token,
            )
        )

    def memory_usage_per_partition(self, index: bool = True, deep: bool = False):
        """Return the memory usage of each partition

        Parameters
        ----------
        index : bool, default True
            Specifies whether to include the memory usage of the index in
            returned Series.
        deep : bool, default False
            If True, introspect the data deeply by interrogating
            ``object`` dtypes for system-level memory consumption, and include
            it in the returned values.

        Returns
        -------
        Series
            A Series whose index is the partition number and whose values
            are the memory usage of each partition in bytes.
        """
        return new_collection(self.expr.memory_usage_per_partition(index, deep))

    @property
    def loc(self):
        """Purely label-location based indexer for selection by label.

        >>> df.loc["b"]  # doctest: +SKIP
        >>> df.loc["b":"d"]  # doctest: +SKIP
        """
        from dask.dataframe.dask_expr._indexing import LocIndexer

        return LocIndexer(self)

    @derived_from(pd.DataFrame)
    def notnull(self):
        return new_collection(expr.NotNull(self))

    @derived_from(pd.DataFrame)
    def isnull(self):
        return ~self.notnull()

    def compute_current_divisions(self, col=None, set_divisions: bool = False):
        """Compute the current divisions of the DataFrame.

        This method triggers immediate computation. If you find yourself running this command
        repeatedly for the same dataframe, we recommend storing the result
        so you don't have to rerun it.

        If the column or index values overlap between partitions, raises ``ValueError``.
        To prevent this, make sure the data are sorted by the column or index.

        Parameters
        ----------
        col : string, optional
            Calculate the divisions for a non-index column by passing in the name of the column.
            If col is not specified, the index will be used to calculate divisions.
            In this case, if the divisions are already known, they will be returned
            immediately without computing.
        set_divisions : bool, default False
            Whether to set the computed divisions into the DataFrame. If False, the divisions
            of the DataFrame are unchanged.

        Examples
        --------
        >>> import dask
        >>> ddf = dask.datasets.timeseries(start="2021-01-01", end="2021-01-07", freq="1h").clear_divisions()
        >>> divisions = ddf.compute_current_divisions()
        >>> print(divisions)  # doctest: +NORMALIZE_WHITESPACE
        (Timestamp('2021-01-01 00:00:00'),
         Timestamp('2021-01-02 00:00:00'),
         Timestamp('2021-01-03 00:00:00'),
         Timestamp('2021-01-04 00:00:00'),
         Timestamp('2021-01-05 00:00:00'),
         Timestamp('2021-01-06 00:00:00'),
         Timestamp('2021-01-06 23:00:00'))

        >>> ddf = ddf.reset_index().clear_divisions()
        >>> divisions = ddf.compute_current_divisions("timestamp")
        >>> print(divisions)  # doctest: +NORMALIZE_WHITESPACE
        (Timestamp('2021-01-01 00:00:00'),
         Timestamp('2021-01-02 00:00:00'),
         Timestamp('2021-01-03 00:00:00'),
         Timestamp('2021-01-04 00:00:00'),
         Timestamp('2021-01-05 00:00:00'),
         Timestamp('2021-01-06 00:00:00'),
         Timestamp('2021-01-06 23:00:00'))

        >>> ddf = ddf.set_index("timestamp", divisions=divisions, sorted=True)
        """
        if col is None and self.known_divisions:
            if set_divisions:
                return self
            return self.divisions

        if col is not None and set_divisions:
            raise NotImplementedError(
                "Can't set divisions of non-index, call set_index instead."
            )

        if col is not None:
            frame = self[col]
        else:
            frame = self.index

        mins, maxes, lens = _compute_partition_stats(frame, allow_overlap=set_divisions)
        divisions = tuple(mins) + (maxes[-1],)
        if not set_divisions:
            return divisions
        if len(mins) == len(self.divisions) - 1:
            if not any(mins[i] >= maxes[i - 1] for i in range(1, len(mins))):
                return new_collection(expr.SetDivisions(self, divisions))

        return new_collection(expr.ResolveOverlappingDivisions(self, mins, maxes, lens))

    @classmethod
    def from_dict(
        cls, data, *, npartitions=1, orient="columns", dtype=None, columns=None
    ):
        """
        Construct a Dask DataFrame from a Python Dictionary

        See Also
        --------
        dask.dataframe.from_dict
        """
        return from_dict(data, npartitions, orient, dtype=dtype, columns=columns)

    def to_json(self, filename, *args, **kwargs):
        """See dd.to_json docstring for more information"""
        from dask.dataframe.io import to_json

        return to_json(self, filename, *args, **kwargs)

    def to_sql(
        self,
        name: str,
        uri: str,
        schema=None,
        if_exists: str = "fail",
        index: bool = True,
        index_label=None,
        chunksize=None,
        dtype=None,
        method=None,
        compute=True,
        parallel=False,
        engine_kwargs=None,
    ):
        from dask.dataframe.io.sql import to_sql

        return to_sql(
            self,
            name,
            uri,
            schema=schema,
            if_exists=if_exists,
            index=index,
            index_label=index_label,
            chunksize=chunksize,
            dtype=dtype,
            method=method,
            compute=compute,
            parallel=parallel,
            engine_kwargs=engine_kwargs,
        )

    def to_orc(self, path, *args, **kwargs):
        """See dd.to_orc docstring for more information"""
        from dask.dataframe.io.orc import to_orc

        return to_orc(self, path, *args, **kwargs)

    def to_csv(self, filename, **kwargs):
        """See dd.to_csv docstring for more information"""
        from dask.dataframe.io.csv import to_csv

        return to_csv(self, filename, **kwargs)

    def to_records(self, index=False, lengths=None):
        from dask.dataframe.dask_expr.io.records import to_records

        if lengths is True:
            lengths = tuple(self.map_partitions(len).compute())
        records = to_records(self)

        chunks = self._validate_chunks(records, lengths)
        records._chunks = (chunks[0],)

        return records

    def _validate_chunks(self, arr, lengths):
        from collections.abc import Sequence

        from dask.array.core import normalize_chunks

        if isinstance(lengths, Sequence):
            lengths = tuple(lengths)

            if len(lengths) != self.npartitions:
                raise ValueError(
                    "The number of items in 'lengths' does not match the number of "
                    f"partitions. {len(lengths)} != {self.npartitions}"
                )

            if self.ndim == 1:
                chunks = normalize_chunks((lengths,))
            else:
                chunks = normalize_chunks((lengths, (len(self.columns),)))

            return chunks
        elif lengths is not None:
            raise ValueError(f"Unexpected value for 'lengths': '{lengths}'")

        return arr._chunks

    def to_bag(self, index=False, format="tuple"):
        """Create a Dask Bag from a Series"""
        from dask.dataframe.dask_expr.io.bag import to_bag

        return to_bag(self, index, format=format)

    def to_hdf(self, path_or_buf, key, mode="a", append=False, **kwargs):
        """See dd.to_hdf docstring for more information"""
        from dask.dataframe.io.hdf import to_hdf

        return to_hdf(self, path_or_buf, key, mode, append, **kwargs)

    def to_delayed(self, optimize_graph=True):
        """Convert into a list of ``dask.delayed`` objects, one per partition.

        Parameters
        ----------
        optimize_graph : bool, optional
            If True [default], the graph is optimized before converting into
            ``dask.delayed`` objects.

        Examples
        --------
        >>> partitions = df.to_delayed()  # doctest: +SKIP

        See Also
        --------
        dask_expr.from_delayed
        """
        from dask.highlevelgraph import HighLevelGraph

        if optimize_graph:
            frame = self.optimize()
        else:
            frame = self
        keys = frame.__dask_keys__()
        graph = frame.__dask_graph__()
        layer = "delayed-" + frame._name
        graph = HighLevelGraph.from_collections(layer, graph, dependencies=())
        return [Delayed(k, graph, layer=layer) for k in keys]

    def to_backend(self, backend: str | None = None, **kwargs):
        """Move to a new DataFrame backend

        Parameters
        ----------
        backend : str, Optional
            The name of the new backend to move to. The default
            is the current "dataframe.backend" configuration.

        Returns
        -------
        DataFrame, Series or Index
        """
        from dask.dataframe.dask_expr._backends import dataframe_creation_dispatch

        # Get desired backend
        backend = backend or dataframe_creation_dispatch.backend
        # Check that "backend" has a registered entrypoint
        backend_entrypoint = dataframe_creation_dispatch.dispatch(backend)
        # Call `DataFrameBackendEntrypoint.to_backend`
        return backend_entrypoint.to_backend(self, **kwargs)

    @derived_from(pd.Series)
    def dot(self, other, meta=no_default):
        if not isinstance(other, FrameBase):
            raise TypeError("The second operand must be a dask dataframe")

        if isinstance(other, DataFrame):
            s = self.map_partitions(M.dot, other, meta=meta)
            return s.groupby(by=s.index).apply(
                lambda x: x.sum(skipna=False), meta=s._meta_nonempty
            )

        return self.map_partitions(_dot_series, other, meta=meta).sum(skipna=False)

    @derived_from(pd.DataFrame)
    def pipe(self, func, *args, **kwargs):
        if isinstance(func, tuple):
            func, target = func
            if target in kwargs:
                raise ValueError(
                    "%s is both the pipe target and a keyword argument" % target
                )
            kwargs[target] = self
            return func(*args, **kwargs)
        else:
            return func(self, *args, **kwargs)

    def sample(self, n=None, frac=None, replace=False, random_state=None):
        """Random sample of items

        Parameters
        ----------
        n : int, optional
            Number of items to return is not supported by dask. Use frac
            instead.
        frac : float, optional
            Approximate fraction of items to return. This sampling fraction is
            applied to all partitions equally. Note that this is an
            **approximate fraction**. You should not expect exactly ``len(df) * frac``
            items to be returned, as the exact number of elements selected will
            depend on how your data is partitioned (but should be pretty close
            in practice).
        replace : boolean, optional
            Sample with or without replacement. Default = False.
        random_state : int or ``np.random.RandomState``
            If an int, we create a new RandomState with this as the seed;
            Otherwise we draw from the passed RandomState.

        See Also
        --------
        DataFrame.random_split
        pandas.DataFrame.sample
        """
        if n is not None:
            msg = (
                "sample does not support the number of sampled items "
                "parameter, 'n'. Please use the 'frac' parameter instead."
            )
            if isinstance(n, Number) and 0 <= n <= 1:
                warnings.warn(msg)
                frac = n
            else:
                raise ValueError(msg)

        if frac is None:
            raise ValueError("frac must not be None")

        if random_state is None:
            random_state = np.random.RandomState()

        state_data = random_state_data(self.npartitions, random_state)
        return new_collection(
            expr.Sample(self, state_data=state_data, frac=frac, replace=replace)
        )

    def _repr_data(self):
        raise NotImplementedError

    @property
    def _repr_divisions(self):
        name = f"npartitions={self.npartitions}"
        if self.known_divisions:
            divisions = pd.Index(self.divisions, name=name)
        else:
            # avoid to be converted to NaN
            divisions = pd.Index([""] * (self.npartitions + 1), name=name)
        return divisions


def _dot_series(*args, **kwargs):
    # .sum() is invoked on each partition before being applied to all
    # partitions. The return type is expected to be a series, not a numpy object
    return meta_series_constructor(args[0])(M.dot(*args, **kwargs))


# Add operator attributes
for op in [
    "__add__",
    "__radd__",
    "__sub__",
    "__rsub__",
    "__mul__",
    "__rmul__",
    "__mod__",
    "__rmod__",
    "__floordiv__",
    "__rfloordiv__",
    "__truediv__",
    "__rtruediv__",
    "__pow__",
    "__rpow__",
    "__lt__",
    "__rlt__",
    "__gt__",
    "__rgt__",
    "__le__",
    "__rle__",
    "__ge__",
    "__rge__",
    "__eq__",
    "__ne__",
    "__and__",
    "__rand__",
    "__or__",
    "__ror__",
    "__xor__",
    "__rxor__",
]:
    setattr(FrameBase, op, functools.partialmethod(_wrap_expr_op, op=op))

for op in [
    "__invert__",
    "__neg__",
    "__pos__",
]:
    setattr(FrameBase, op, functools.partialmethod(_wrap_unary_expr_op, op=op))


class DataFrame(FrameBase):
    """DataFrame-like Expr Collection.

    The constructor takes the expression that represents the query as input. The class
    is not meant to be instantiated directly. Instead, use one of the IO connectors from
    Dask.
    """

    _accessors: ClassVar[set[str]] = set()
    _partition_type = pd.DataFrame

    @property
    def shape(self):
        return self.size // max(len(self.columns), 1), len(self.columns)

    @property
    def ndim(self):
        """Return dimensionality"""
        return 2

    @property
    def empty(self):
        # __getattr__ will be called after we raise this, so we'll raise it again from there
        raise AttributeNotImplementedError(
            "Checking whether a Dask DataFrame has any rows may be expensive. "
            "However, checking the number of columns is fast. "
            "Depending on which of these results you need, use either "
            "`len(df.index) == 0` or `len(df.columns) == 0`"
        )

    @derived_from(pd.DataFrame)
    def items(self):
        for i, name in enumerate(self.columns):
            yield (name, self.iloc[:, i])

    @property
    def axes(self):
        return [self.index, self.columns]

    def __contains__(self, key):
        return key in self._meta

    def __iter__(self):
        return iter(self._meta)

    def __dataframe__(self, *args, **kwargs):
        from dask.dataframe.dask_expr._interchange import DaskDataFrameInterchange

        return DaskDataFrameInterchange(self)

    @derived_from(pd.DataFrame)
    def iterrows(self):
        frame = self.optimize()
        for i in range(self.npartitions):
            df = frame.get_partition(i).compute()
            yield from df.iterrows()

    @derived_from(pd.DataFrame)
    def itertuples(self, index=True, name="Pandas"):
        frame = self.optimize()
        for i in range(self.npartitions):
            df = frame.get_partition(i).compute()
            yield from df.itertuples(index=index, name=name)

    @property
    def _elemwise(self):
        return elemwise

    def __array_ufunc__(self, numpy_ufunc, method, *inputs, **kwargs):
        out = kwargs.get("out", ())
        for x in inputs + out:
            # ufuncs work with 0-dimensional NumPy ndarrays
            # so we don't want to raise NotImplemented
            if isinstance(x, np.ndarray) and x.shape == ():
                continue
            elif not isinstance(
                x, (Number, Scalar, FrameBase, Array, pd.DataFrame, pd.Series, pd.Index)
            ):
                return NotImplemented

        if method == "__call__":
            if numpy_ufunc.signature is not None:
                return NotImplemented
            if numpy_ufunc.nout > 1:
                # ufuncs with multiple output values
                # are not yet supported for frames
                return NotImplemented
            else:
                return elemwise(numpy_ufunc, *inputs, **kwargs)
        else:
            # ufunc methods are not yet supported for frames
            return NotImplemented

    def __array_wrap__(self, array, context=None):
        if isinstance(context, tuple) and len(context) > 0:
            if isinstance(context[1][0], np.ndarray) and context[1][0].shape == ():
                index = None
            else:
                index = context[1][0].index
        else:
            try:
                import inspect

                method_name = f"`{inspect.stack()[3][3]}`"
            except IndexError:
                method_name = "This method"

            raise NotImplementedError(
                f"{method_name} is not implemented for `dask.dataframe.DataFrame`."
            )

        return meta_frame_constructor(self)(array, index=index, columns=self.columns)

    def _ipython_key_completions_(self):
        return methods.tolist(self.columns)

    def _repr_html_(self):
        return self.to_html()

    @derived_from(pd.DataFrame)
    def assign(self, **pairs):
        result = self
        args = []
        for k, v in pairs.items():
            v = _maybe_from_pandas([v])[0]
            if not isinstance(k, str):
                raise TypeError(f"Column name cannot be type {type(k)}")

            if callable(v):
                result = new_collection(expr.Assign(result, *args))
                args = []
                result = new_collection(expr.Assign(result, k, v(result)))
                continue

            elif isinstance(v, (Scalar, Series)):
                if isinstance(v, Series):
                    if not expr.are_co_aligned(result.expr, v.expr):
                        if len(args) > 0:
                            result = expr.Assign(result, *args)
                            args = []
                        result = new_collection(expr.AssignAlign(result, k, v.expr))
                        continue

            elif not isinstance(v, FrameBase) and isinstance(v, Hashable):
                pass
            elif isinstance(v, Array):
                if len(v.shape) > 1:
                    raise ValueError("Array assignment only supports 1-D arrays")
                if v.npartitions != result.npartitions:
                    raise ValueError(
                        "Number of partitions do not match "
                        f"({v.npartitions} != {result.npartitions})"
                    )
                v = from_dask_array(v, index=result.index, meta=result._meta)
            else:
                raise TypeError(f"Column assignment doesn't support type {type(v)}")
            args.extend([k, v])

        if len(args) > 0:
            result = new_collection(expr.Assign(result, *args))

        return result

    @derived_from(pd.DataFrame)
    def clip(self, lower=None, upper=None, axis=None, **kwargs):
        axis = self._validate_axis(axis)
        if axis == 1:
            return self.map_partitions(M.clip, lower, upper, axis=axis)
        return new_collection(self.expr.clip(lower, upper, axis))

    def merge(
        self,
        right,
        how="inner",
        on=None,
        left_on=None,
        right_on=None,
        left_index=False,
        right_index=False,
        suffixes=("_x", "_y"),
        indicator=False,
        shuffle_method=None,
        npartitions=None,
        broadcast=None,
    ):
        """Merge the DataFrame with another DataFrame

        This will merge the two datasets, either on the indices, a certain column
        in each dataset or the index in one dataset and the column in another.

        Parameters
        ----------
        right: dask.dataframe.DataFrame
        how : {'left', 'right', 'outer', 'inner', 'leftsemi'}, default: 'inner'
            How to handle the operation of the two objects:

            - left: use calling frame's index (or column if on is specified)
            - right: use other frame's index
            - outer: form union of calling frame's index (or column if on is
              specified) with other frame's index, and sort it
              lexicographically
            - inner: form intersection of calling frame's index (or column if
              on is specified) with other frame's index, preserving the order
              of the calling's one
            - leftsemi: Choose all rows in left where the join keys can be found
              in right. Won't duplicate rows if the keys are duplicated in right.
              Drops all columns from right.

        on : label or list
            Column or index level names to join on. These must be found in both
            DataFrames. If on is None and not merging on indexes then this
            defaults to the intersection of the columns in both DataFrames.
        left_on : label or list, or array-like
            Column to join on in the left DataFrame. Other than in pandas
            arrays and lists are only support if their length is 1.
        right_on : label or list, or array-like
            Column to join on in the right DataFrame. Other than in pandas
            arrays and lists are only support if their length is 1.
        left_index : boolean, default False
            Use the index from the left DataFrame as the join key.
        right_index : boolean, default False
            Use the index from the right DataFrame as the join key.
        suffixes : 2-length sequence (tuple, list, ...)
            Suffix to apply to overlapping column names in the left and
            right side, respectively
        indicator : boolean or string, default False
            If True, adds a column to output DataFrame called "_merge" with
            information on the source of each row. If string, column with
            information on source of each row will be added to output DataFrame,
            and column will be named value of string. Information column is
            Categorical-type and takes on a value of "left_only" for observations
            whose merge key only appears in `left` DataFrame, "right_only" for
            observations whose merge key only appears in `right` DataFrame,
            and "both" if the observation’s merge key is found in both.
        npartitions: int or None, optional
            The ideal number of output partitions. This is only utilised when
            performing a hash_join (merging on columns only). If ``None`` then
            ``npartitions = max(lhs.npartitions, rhs.npartitions)``.
            Default is ``None``.
        shuffle_method: {'disk', 'tasks', 'p2p'}, optional
            Either ``'disk'`` for single-node operation or ``'tasks'`` and
            ``'p2p'``` for distributed operation.  Will be inferred by your
            current scheduler.
        broadcast: boolean or float, optional
            Whether to use a broadcast-based join in lieu of a shuffle-based
            join for supported cases.  By default, a simple heuristic will be
            used to select the underlying algorithm. If a floating-point value
            is specified, that number will be used as the ``broadcast_bias``
            within the simple heuristic (a large number makes Dask more likely
            to choose the ``broacast_join`` code path). See ``broadcast_join``
            for more information.

        Notes
        -----

        There are three ways to join dataframes:

        1. Joining on indices. In this case the divisions are
           aligned using the function ``dask.dataframe.multi.align_partitions``.
           Afterwards, each partition is merged with the pandas merge function.

        2. Joining one on index and one on column. In this case the divisions of
           dataframe merged by index (:math:`d_i`) are used to divide the column
           merged dataframe (:math:`d_c`) one using
           ``dask.dataframe.multi.rearrange_by_divisions``. In this case the
           merged dataframe (:math:`d_m`) has the exact same divisions
           as (:math:`d_i`). This can lead to issues if you merge multiple rows from
           (:math:`d_c`) to one row in (:math:`d_i`).

        3. Joining both on columns. In this case a hash join is performed using
           ``dask.dataframe.multi.hash_join``.

        In some cases, you may see a ``MemoryError`` if the ``merge`` operation requires
        an internal ``shuffle``, because shuffling places all rows that have the same
        index in the same partition. To avoid this error, make sure all rows with the
        same ``on``-column value can fit on a single partition.
        """
        return merge(
            self,
            right,
            how,
            on,
            left_on,
            right_on,
            left_index,
            right_index,
            suffixes,
            indicator,
            shuffle_method,
            npartitions=npartitions,
            broadcast=broadcast,
        )

    @derived_from(pd.DataFrame)
    def join(
        self,
        other,
        on=None,
        how="left",
        lsuffix="",
        rsuffix="",
        shuffle_method=None,
        npartitions=None,
    ):
        if not isinstance(other, list) and not is_dask_collection(other):
            other = from_pandas(other, npartitions=1)
        if (
            not isinstance(other, list)
            and not is_dataframe_like(other._meta)
            and hasattr(other._meta, "name")
        ):
            other = new_collection(expr.ToFrame(other))

        if not isinstance(other, FrameBase):
            if not isinstance(other, list) or not all(
                isinstance(o, FrameBase) for o in other
            ):
                raise ValueError("other must be DataFrame or list of DataFrames")
            if how not in ("outer", "left"):
                raise ValueError("merge_multi only supports left or outer joins")

            other = [
                from_pandas(o, npartitions=1) if not is_dask_collection(o) else o
                for o in other
            ]

            return new_collection(
                JoinRecursive([self.expr] + [o.expr for o in other], how=how)
            )

        return self.merge(
            right=other,
            left_index=on is None,
            right_index=True,
            left_on=on,
            how=how,
            suffixes=(lsuffix, rsuffix),
            shuffle_method=shuffle_method,
            npartitions=npartitions,
        )

    @derived_from(pd.DataFrame)
    def groupby(
        self, by, group_keys=True, sort=None, observed=None, dropna=None, **kwargs
    ):
        from dask.dataframe.dask_expr._groupby import GroupBy

        if isinstance(by, FrameBase) and not isinstance(by, Series):
            raise ValueError(
                f"`by` must be a column name or list of columns, got {by}."
            )

        return GroupBy(
            self,
            by,
            group_keys=group_keys,
            sort=sort,
            observed=observed,
            dropna=dropna,
            **kwargs,
        )

    def __setitem__(self, key, value):
        if isinstance(key, (tuple, list)) and isinstance(value, DataFrame):
            out = self.assign(**{k: value[c] for k, c in zip(key, value.columns)})

        elif isinstance(key, pd.Index) and not isinstance(value, DataFrame):
            out = self.assign(**{k: value for k in list(key)})
        elif (
            is_dataframe_like(key)
            or is_series_like(key)
            or isinstance(key, (DataFrame, Series))
        ):
            out = self.where(~key, value)
        elif not isinstance(key, str):
            raise NotImplementedError(f"Item assignment with {type(key)} not supported")
        else:
            out = self.assign(**{key: value})
        self._expr = out._expr

    def __delitem__(self, key):
        columns = [c for c in self.columns if c != key]
        out = self[columns]
        self._expr = out._expr

    def __setattr__(self, key, value):
        try:
            columns = object.__getattribute__(self, "_expr").columns
        except AttributeError:
            columns = ()

        # exclude protected attributes from setitem
        if key in columns and key not in [
            "divisions",
            "dask",
            "_name",
            "_meta",
            "_expr",
        ]:
            self[key] = value
        else:
            object.__setattr__(self, key, value)

    def __getattr__(self, key):
        try:
            # Prioritize `DataFrame` attributes
            return object.__getattribute__(self, key)
        except AttributeError as err:
            try:
                # Check if key is in columns if key
                # is not a normal attribute
                if key in self.expr._meta.columns:
                    return new_collection(self.expr[key])
                raise err
            except AttributeError:
                # Fall back to `BaseFrame.__getattr__`
                return super().__getattr__(key)

    def __dir__(self):
        o = set(dir(type(self)))
        o.update(self.__dict__)
        o.update(set(dir(expr.Expr)))
        o.update(c for c in self.columns if (isinstance(c, str) and c.isidentifier()))
        return list(o)

    def map(self, func, na_action=None, meta=None):
        if not PANDAS_GE_210:
            raise NotImplementedError(
                f"DataFrame.map requires pandas>=2.1.0, but pandas={PANDAS_VERSION} is "
                "installed."
            )
        if meta is None:
            meta = expr.emulate(M.map, self, func, na_action=na_action, udf=True)
            warnings.warn(meta_warning(meta, method="map"))
        return new_collection(expr.Map(self, arg=func, na_action=na_action, meta=meta))

    @derived_from(pd.DataFrame)
    def nlargest(self, n=5, columns=None, split_every=None):
        return new_collection(
            NLargest(self, n=n, _columns=columns, split_every=split_every)
        )

    @derived_from(pd.DataFrame)
    def nsmallest(self, n=5, columns=None, split_every=None):
        return new_collection(
            NSmallest(self, n=n, _columns=columns, split_every=split_every)
        )

    @derived_from(pd.DataFrame)
    def memory_usage(self, deep=False, index=True):
        return new_collection(MemoryUsageFrame(self, deep=deep, _index=index))

    @derived_from(pd.DataFrame)
    def combine(self, other, func, fill_value=None, overwrite=True):
        other = self._create_alignable_frame(other, "outer")
        if not expr.are_co_aligned(self.expr, other.expr):
            return new_collection(
                expr.CombineFrameAlign(self, other, func, fill_value, overwrite)
            )
        return new_collection(
            expr.CombineFrame(self, other, func, fill_value, overwrite)
        )

    @derived_from(
        pd.DataFrame,
        inconsistencies="keep=False will raise a ``NotImplementedError``",
    )
    def drop_duplicates(
        self,
        subset=None,
        split_every=None,
        split_out=True,
        shuffle_method=None,
        ignore_index=False,
        keep="first",
    ):
        shuffle_method = _get_shuffle_preferring_order(shuffle_method)
        if keep is False:
            raise NotImplementedError("drop_duplicates with keep=False")
        # Fail early if subset is not valid, e.g. missing columns
        subset = _convert_to_list(subset)
        meta_nonempty(self._meta).drop_duplicates(subset=subset, keep=keep)
        return new_collection(
            DropDuplicates(
                self,
                subset=subset,
                ignore_index=ignore_index,
                split_out=split_out,
                split_every=split_every,
                shuffle_method=shuffle_method,
                keep=keep,
            )
        )

    @insert_meta_param_description(pad=12)
    def apply(self, function, *args, meta=no_default, axis=0, **kwargs):
        """Parallel version of pandas.DataFrame.apply

        This mimics the pandas version except for the following:

        1.  Only ``axis=1`` is supported (and must be specified explicitly).
        2.  The user should provide output metadata via the `meta` keyword.

        Parameters
        ----------
        func : function
            Function to apply to each column/row
        axis : {0 or 'index', 1 or 'columns'}, default 0
            - 0 or 'index': apply function to each column (NOT SUPPORTED)
            - 1 or 'columns': apply function to each row
        $META
        args : tuple
            Positional arguments to pass to function in addition to the array/series

        Additional keyword arguments will be passed as keywords to the function

        Returns
        -------
        applied : Series or DataFrame

        Examples
        --------
        >>> import pandas as pd
        >>> import dask.dataframe as dd
        >>> df = pd.DataFrame({'x': [1, 2, 3, 4, 5],
        ...                    'y': [1., 2., 3., 4., 5.]})
        >>> ddf = dd.from_pandas(df, npartitions=2)

        Apply a function to row-wise passing in extra arguments in ``args`` and
        ``kwargs``:

        >>> def myadd(row, a, b=1):
        ...     return row.sum() + a + b
        >>> res = ddf.apply(myadd, axis=1, args=(2,), b=1.5)  # doctest: +SKIP

        By default, dask tries to infer the output metadata by running your
        provided function on some fake data. This works well in many cases, but
        can sometimes be expensive, or even fail. To avoid this, you can
        manually specify the output metadata with the ``meta`` keyword. This
        can be specified in many forms, for more information see
        ``dask.dataframe.utils.make_meta``.

        Here we specify the output is a Series with name ``'x'``, and dtype
        ``float64``:

        >>> res = ddf.apply(myadd, axis=1, args=(2,), b=1.5, meta=('x', 'f8'))

        In the case where the metadata doesn't change, you can also pass in
        the object itself directly:

        >>> res = ddf.apply(lambda row: row + 1, axis=1, meta=ddf)

        See Also
        --------
        DataFrame.map_partitions
        """
        axis = self._validate_axis(axis)
        if axis == 0:
            msg = (
                "Dask DataFrame.apply only supports axis=1\n"
                "  Try: df.apply(func, axis=1)"
            )
            raise NotImplementedError(msg)
        if meta is no_default:
            meta = expr.emulate(
                M.apply, self, function, args=args, udf=True, axis=axis, **kwargs
            )
            warnings.warn(meta_warning(meta))
        return new_collection(
            self.expr.apply(function, *args, meta=meta, axis=axis, **kwargs)
        )

    @derived_from(pd.DataFrame)
    def dropna(self, how=no_default, subset=None, thresh=no_default):
        if how is not no_default and thresh is not no_default:
            raise TypeError(
                "You cannot set both the how and thresh arguments at the same time."
            )
        subset = _convert_to_list(subset)
        return new_collection(
            expr.DropnaFrame(self, how=how, subset=subset, thresh=thresh)
        )

    @classmethod
    def _validate_axis(cls, axis=0, numeric_axis: bool = True) -> None | Literal[0, 1]:
        if axis not in (0, 1, "index", "columns", None):
            raise ValueError(f"No axis named {axis}")
        if numeric_axis:
            num_axis: dict[str | None, Literal[0, 1]] = {"index": 0, "columns": 1}
            return num_axis.get(axis, axis)
        else:
            return axis

    @derived_from(pd.DataFrame, ua_args=["index"])
    def rename(self, index=None, columns=None):
        if index is not None:
            raise ValueError("Cannot rename index.")
        return new_collection(expr.RenameFrame(self, columns=columns))

    @derived_from(pd.DataFrame)
    def squeeze(self, axis=None):
        if axis in [None, 1]:
            if len(self.columns) == 1:
                return self[self.columns[0]]
            else:
                return self

        elif axis == 0:
            raise NotImplementedError(
                f"{type(self)} does not support squeeze along axis 0"
            )

        else:
            raise ValueError(f"No axis {axis} for object type {type(self)}")

    @derived_from(pd.DataFrame)
    def explode(self, column):
        column = _convert_to_list(column)
        return new_collection(expr.ExplodeFrame(self, column=column))

    @derived_from(pd.DataFrame)
    def drop(self, labels=None, axis=0, columns=None, errors="raise"):
        if columns is None and labels is None:
            raise TypeError("must either specify 'columns' or 'labels'")

        axis = _validate_axis(axis)

        if axis == 1:
            columns = labels or columns
        elif axis == 0 and columns is None:
            raise NotImplementedError(
                "Drop currently only works for axis=1 or when columns is not None"
            )
        return new_collection(expr.Drop(self, columns=columns, errors=errors))

    def to_parquet(self, path, **kwargs):
        from dask.dataframe.dask_expr.io.parquet import to_parquet

        return to_parquet(self, path, **kwargs)

    @derived_from(pd.DataFrame)
    def select_dtypes(self, include=None, exclude=None):
        columns = list(
            self._meta.select_dtypes(include=include, exclude=exclude).columns
        )
        return new_collection(self.expr[columns])

    @derived_from(pd.DataFrame)
    def eval(self, expr, **kwargs):
        if "inplace" in kwargs:
            raise NotImplementedError("inplace is not supported for eval")
        return new_collection(Eval(self, _expr=expr, expr_kwargs=kwargs))

    def set_index(
        self,
        other,
        drop=True,
        sorted=False,
        npartitions: int | None = None,
        divisions=None,
        sort: bool = True,
        shuffle_method=None,
        upsample: float = 1.0,
        partition_size: float = 128e6,
        append: bool = False,
        **options,
    ):
        """Set the DataFrame index (row labels) using an existing column.

        If ``sort=False``, this function operates exactly like ``pandas.set_index``
        and sets the index on the DataFrame. If ``sort=True`` (default),
        this function also sorts the DataFrame by the new index. This can have a
        significant impact on performance, because joins, groupbys, lookups, etc.
        are all much faster on that column. However, this performance increase
        comes with a cost, sorting a parallel dataset requires expensive shuffles.
        Often we ``set_index`` once directly after data ingest and filtering and
        then perform many cheap computations off of the sorted dataset.

        With ``sort=True``, this function is much more expensive. Under normal
        operation this function does an initial pass over the index column to
        compute approximate quantiles to serve as future divisions. It then passes
        over the data a second time, splitting up each input partition into several
        pieces and sharing those pieces to all of the output partitions now in
        sorted order.

        In some cases we can alleviate those costs, for example if your dataset is
        sorted already then we can avoid making many small pieces or if you know
        good values to split the new index column then we can avoid the initial
        pass over the data. For example if your new index is a datetime index and
        your data is already sorted by day then this entire operation can be done
        for free. You can control these options with the following parameters.

        Parameters
        ----------
        other: string or Dask Series
            Column to use as index.
        drop: boolean, default True
            Delete column to be used as the new index.
        sorted: bool, optional
            If the index column is already sorted in increasing order.
            Defaults to False
        npartitions: int, None, or 'auto'
            The ideal number of output partitions. If None, use the same as
            the input. If 'auto' then decide by memory use.
            Only used when ``divisions`` is not given. If ``divisions`` is given,
            the number of output partitions will be ``len(divisions) - 1``.
        divisions: list, optional
            The "dividing lines" used to split the new index into partitions.
            For ``divisions=[0, 10, 50, 100]``, there would be three output partitions,
            where the new index contained [0, 10), [10, 50), and [50, 100), respectively.
            See https://docs.dask.org/en/latest/dataframe-design.html#partitions.
            If not given (default), good divisions are calculated by immediately computing
            the data and looking at the distribution of its values. For large datasets,
            this can be expensive.
            Note that if ``sorted=True``, specified divisions are assumed to match
            the existing partitions in the data; if this is untrue you should
            leave divisions empty and call ``repartition`` after ``set_index``.
        inplace: bool, optional
            Modifying the DataFrame in place is not supported by Dask.
            Defaults to False.
        sort: bool, optional
            If ``True``, sort the DataFrame by the new index. Otherwise
            set the index on the individual existing partitions.
            Defaults to ``True``.
        shuffle_method: {'disk', 'tasks', 'p2p'}, optional
            Either ``'disk'`` for single-node operation or ``'tasks'`` and
            ``'p2p'`` for distributed operation.  Will be inferred by your
            current scheduler.
        compute: bool, default False
            Whether or not to trigger an immediate computation. Defaults to False.
            Note, that even if you set ``compute=False``, an immediate computation
            will still be triggered if ``divisions`` is ``None``.
        partition_size: int, optional
            Desired size of each partitions in bytes.
            Only used when ``npartitions='auto'``

        Examples
        --------
        >>> import dask
        >>> ddf = dask.datasets.timeseries(start="2021-01-01", end="2021-01-07", freq="1h").reset_index()
        >>> ddf2 = ddf.set_index("x")
        >>> ddf2 = ddf.set_index(ddf.x)
        >>> ddf2 = ddf.set_index(ddf.timestamp, sorted=True)

        A common case is when we have a datetime column that we know to be
        sorted and is cleanly divided by day.  We can set this index for free
        by specifying both that the column is pre-sorted and the particular
        divisions along which is is separated

        >>> import pandas as pd
        >>> divisions = pd.date_range(start="2021-01-01", end="2021-01-07", freq='1D')
        >>> divisions
        DatetimeIndex(['2021-01-01', '2021-01-02', '2021-01-03', '2021-01-04',
                       '2021-01-05', '2021-01-06', '2021-01-07'],
                      dtype='datetime64[ns]', freq='D')

        Note that ``len(divisons)`` is equal to ``npartitions + 1``. This is because ``divisions``
        represents the upper and lower bounds of each partition. The first item is the
        lower bound of the first partition, the second item is the lower bound of the
        second partition and the upper bound of the first partition, and so on.
        The second-to-last item is the lower bound of the last partition, and the last
        (extra) item is the upper bound of the last partition.

        >>> ddf2 = ddf.set_index("timestamp", sorted=True, divisions=divisions.tolist())

        If you'll be running `set_index` on the same (or similar) datasets repeatedly,
        you could save time by letting Dask calculate good divisions once, then copy-pasting
        them to reuse. This is especially helpful running in a Jupyter notebook:

        >>> ddf2 = ddf.set_index("name")  # slow, calculates data distribution
        >>> ddf2.divisions  # doctest: +SKIP
        ["Alice", "Laura", "Ursula", "Zelda"]
        >>> # ^ Now copy-paste this and edit the line above to:
        >>> # ddf2 = ddf.set_index("name", divisions=["Alice", "Laura", "Ursula", "Zelda"])
        """
        if isinstance(other, list) and len(other) == 1:
            other = other[0]
        if isinstance(other, list):
            if any([isinstance(c, FrameBase) for c in other]):
                raise TypeError("List[FrameBase] not supported by set_index")
            else:
                raise NotImplementedError(
                    "Dask dataframe does not yet support multi-indexes.\n"
                    f"You tried to index with this index: {other}\n"
                    "Indexes must be single columns only."
                )
        if isinstance(other, DataFrame):
            raise NotImplementedError(
                "Dask dataframe does not yet support multi-indexes.\n"
                f"You tried to index with a frame with these columns: {list(other.columns)}\n"
                "Indexes must be single columns only."
            )
        if isinstance(other, Series):
            if other._name == self.index._name:
                return self
        elif other == self.index.name:
            return self

        if divisions is not None:
            check_divisions(divisions)

        if (sorted or not sort) and npartitions is not None:
            raise ValueError(
                "Specifying npartitions with sort=False or sorted=True is not "
                "supported. Call `repartition` afterwards."
            )

        if sorted:
            if divisions is not None and len(divisions) - 1 != self.npartitions:
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
            result = new_collection(
                SetIndexBlockwise(
                    self, other, drop, new_divisions=divisions, append=append
                )
            )
            return result.compute_current_divisions(set_divisions=True)
        elif not sort:
            return new_collection(
                SetIndexBlockwise(self, other, drop, None, append=append)
            )

        return new_collection(
            SetIndex(
                self,
                other,
                drop,
                user_divisions=divisions,
                npartitions=npartitions,
                upsample=upsample,
                partition_size=partition_size,
                shuffle_method=get_specified_shuffle(shuffle_method),
                append=append,
                options=options,
            )
        )

    def sort_values(
        self,
        by: str | list[str],
        npartitions: int | None = None,
        ascending: bool | list[bool] = True,
        na_position: Literal["first"] | Literal["last"] = "last",
        partition_size: float = 128e6,
        sort_function: Callable[[pd.DataFrame], pd.DataFrame] | None = None,
        sort_function_kwargs: Mapping[str, Any] | None = None,
        upsample: float = 1.0,
        ignore_index: bool | None = False,
        shuffle_method: str | None = None,
        **options,
    ):
        """Sort the dataset by a single column.

        Sorting a parallel dataset requires expensive shuffles and is generally
        not recommended. See ``set_index`` for implementation details.

        Parameters
        ----------
        by: str or list[str]
            Column(s) to sort by.
        npartitions: int, None, or 'auto'
            The ideal number of output partitions. If None, use the same as
            the input. If 'auto' then decide by memory use.
        ascending: bool, optional
            Sort ascending vs. descending.
            Defaults to True.
        na_position: {'last', 'first'}, optional
            Puts NaNs at the beginning if 'first', puts NaN at the end if 'last'.
            Defaults to 'last'.
        sort_function: function, optional
            Sorting function to use when sorting underlying partitions.
            If None, defaults to ``M.sort_values`` (the partition library's
            implementation of ``sort_values``).
        sort_function_kwargs: dict, optional
            Additional keyword arguments to pass to the partition sorting function.
            By default, ``by``, ``ascending``, and ``na_position`` are provided.

        Examples
        --------
        >>> df2 = df.sort_values('x')  # doctest: +SKIP
        """
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

        if not isinstance(ascending, bool) and not len(ascending) == len(by):
            raise ValueError(f"Length of {ascending=} != length of {by=}")

        return new_collection(
            SortValues(
                self,
                by,
                ascending,
                na_position,
                npartitions,
                partition_size,
                sort_function,
                sort_function_kwargs,
                upsample,
                ignore_index,
                get_specified_shuffle(shuffle_method),
                options=options,
            )
        )

    def query(self, expr, **kwargs):
        """Filter dataframe with complex expression

        Blocked version of pd.DataFrame.query

        Parameters
        ----------
        expr: str
            The query string to evaluate.
            You can refer to column names that are not valid Python variable names
            by surrounding them in backticks.
            Dask does not fully support referring to variables using the '@' character,
            use f-strings or the ``local_dict`` keyword argument instead.

        Notes
        -----
        This is like the sequential version except that this will also happen
        in many threads.  This may conflict with ``numexpr`` which will use
        multiple threads itself.  We recommend that you set ``numexpr`` to use a
        single thread:

        .. code-block:: python

            import numexpr
            numexpr.set_num_threads(1)

        See also
        --------
        pandas.DataFrame.query
        pandas.eval

        Examples
        --------
        >>> import pandas as pd
        >>> import dask.dataframe as dd
        >>> df = pd.DataFrame({'x': [1, 2, 1, 2],
        ...                    'y': [1, 2, 3, 4],
        ...                    'z z': [4, 3, 2, 1]})
        >>> ddf = dd.from_pandas(df, npartitions=2)

        Refer to column names directly:

        >>> ddf.query('y > x').compute()
           x  y  z z
        2  1  3    2
        3  2  4    1

        Refer to column name using backticks:

        >>> ddf.query('`z z` > x').compute()
           x  y  z z
        0  1  1    4
        1  2  2    3
        2  1  3    2

        Refer to variable name using f-strings:

        >>> value = 1
        >>> ddf.query(f'x == {value}').compute()
           x  y  z z
        0  1  1    4
        2  1  3    2

        Refer to variable name using ``local_dict``:

        >>> ddf.query('x == @value', local_dict={"value": value}).compute()
           x  y  z z
        0  1  1    4
        2  1  3    2
        """
        return new_collection(Query(self, expr, kwargs))

    @derived_from(pd.DataFrame)
    def mode(self, dropna=True, split_every=False, numeric_only=False):
        modes = []
        for _, col in self.items():
            if numeric_only and not pd.api.types.is_numeric_dtype(col.dtype):
                continue
            modes.append(col.mode(dropna=dropna, split_every=split_every))
        return concat(modes, axis=1)

    @derived_from(pd.DataFrame)
    def add_prefix(self, prefix):
        return new_collection(expr.AddPrefix(self, prefix))

    @derived_from(pd.DataFrame)
    def add_suffix(self, suffix):
        return new_collection(expr.AddSuffix(self, suffix))

    def pivot_table(self, index, columns, values, aggfunc="mean"):
        """
        Create a spreadsheet-style pivot table as a DataFrame. Target ``columns``
        must have category dtype to infer result's ``columns``.
        ``index``, ``columns``, ``values`` and ``aggfunc`` must be all scalar.

        Parameters
        ----------
        values : scalar
            column to aggregate
        index : scalar
            column to be index
        columns : scalar
            column to be columns
        aggfunc : {'mean', 'sum', 'count'}, default 'mean'

        Returns
        -------
        table : DataFrame
        """
        return pivot_table(self, index, columns, values, aggfunc)

    @property
    def iloc(self):
        """Purely integer-location based indexing for selection by position.

        Only indexing the column positions is supported. Trying to select
        row positions will raise a ValueError.

        See :ref:`dataframe.indexing` for more.

        Examples
        --------
        >>> df.iloc[:, [2, 0, 1]]  # doctest: +SKIP
        """
        from dask.dataframe.dask_expr._indexing import ILocIndexer

        return ILocIndexer(self)

    def _comparison_op(self, expr_cls, other, level, axis):
        if level is not None:
            raise NotImplementedError("level must be None")
        axis = self._validate_axis(axis)
        return new_collection(expr_cls(self, other, axis))

    def lt(self, other, level=None, axis=0):
        return self._comparison_op(expr.LTFrame, other, level, axis)

    def le(self, other, level=None, axis=0):
        return self._comparison_op(expr.LEFrame, other, level, axis)

    def gt(self, other, level=None, axis=0):
        return self._comparison_op(expr.GTFrame, other, level, axis)

    def ge(self, other, level=None, axis=0):
        return self._comparison_op(expr.GEFrame, other, level, axis)

    def ne(self, other, level=None, axis=0):
        return self._comparison_op(expr.NEFrame, other, level, axis)

    def eq(self, other, level=None, axis=0):
        return self._comparison_op(expr.EQFrame, other, level, axis)

    def categorize(self, columns=None, index=None, split_every=None, **kwargs):
        """Convert columns of the DataFrame to category dtype.

        .. warning:: This method eagerly computes the categories of the chosen columns.

        Parameters
        ----------
        columns : list, optional
            A list of column names to convert to categoricals. By default any
            column with an object dtype is converted to a categorical, and any
            unknown categoricals are made known.
        index : bool, optional
            Whether to categorize the index. By default, object indices are
            converted to categorical, and unknown categorical indices are made
            known. Set True to always categorize the index, False to never.
        split_every : int, optional
            Group partitions into groups of this size while performing a
            tree-reduction. If set to False, no tree-reduction will be used.
        kwargs
            Keyword arguments are passed on to compute.
        """
        df = self
        meta = df._meta
        if columns is None:
            columns = list(meta.select_dtypes(["object", "string", "category"]).columns)
        elif is_scalar(columns):
            columns = [columns]

        # Filter out known categorical columns
        columns = [
            c
            for c in columns
            if not (is_categorical_dtype(meta[c]) and has_known_categories(meta[c]))
        ]

        if index is not False:
            if is_categorical_dtype(meta.index):
                index = not has_known_categories(meta.index)
            elif index is None:
                index = str(meta.index.dtype) in ("object", "string")

        # Nothing to do
        if not len(columns) and index is False:
            return df

        from dask.dataframe.dask_expr._collection import new_collection

        # Eagerly compute the categories
        categories, index = new_collection(
            GetCategories(self, columns=columns, index=index, split_every=split_every)
        ).compute()

        # Some operations like get_dummies() rely on the order of categories
        categories = {k: v.sort_values() for k, v in categories.items()}

        # Categorize each partition
        return new_collection(Categorize(self, categories, index))

    @derived_from(pd.DataFrame)
    def nunique(self, axis=0, dropna=True, split_every=False):
        if axis == 1:
            return new_collection(expr.NUniqueColumns(self, axis=axis, dropna=dropna))
        else:
            return concat(
                [
                    col.nunique(dropna=dropna, split_every=split_every).to_series(name)
                    for name, col in self.items()
                ]
            )

    def quantile(self, q=0.5, axis=0, numeric_only=False, method="default"):
        """Approximate row-wise and precise column-wise quantiles of DataFrame

        Parameters
        ----------
        q : list/array of floats, default 0.5 (50%)
            Iterable of numbers ranging from 0 to 1 for the desired quantiles
        axis : {0, 1, 'index', 'columns'} (default 0)
            0 or 'index' for row-wise, 1 or 'columns' for column-wise
        method : {'default', 'tdigest', 'dask'}, optional
            What method to use. By default will use dask's internal custom
            algorithm (``'dask'``).  If set to ``'tdigest'`` will use tdigest
            for floats and ints and fallback to the ``'dask'`` otherwise.
        """
        allowed_methods = ["default", "dask", "tdigest"]
        if method not in allowed_methods:
            raise ValueError("method can only be 'default', 'dask' or 'tdigest'")
        meta = make_meta(
            meta_nonempty(self._meta).quantile(
                q=q, numeric_only=numeric_only, axis=axis
            )
        )

        if axis == 1:
            if isinstance(q, list):
                # Not supported, the result will have current index as columns
                raise ValueError("'q' must be scalar when axis=1 is specified")

            return self.map_partitions(
                M.quantile,
                q,
                axis,
                enforce_metadata=False,
                meta=meta,
                numeric_only=numeric_only,
            )

        if numeric_only:
            frame = self.select_dtypes(
                "number", exclude=[np.timedelta64, np.datetime64]
            )
        else:
            frame = self

        collections = []
        for _, col in frame.items():
            collections.append(col.quantile(q=q, method=method))

        if len(collections) > 0 and isinstance(collections[0], Scalar):
            return _from_scalars(collections, meta, frame.expr.columns)

        return concat(collections, axis=1)

    @derived_from(pd.DataFrame)
    def median(self, axis=0, numeric_only=False):
        if axis == 1 or self.npartitions == 1:
            return self.median_approximate(axis=axis, numeric_only=numeric_only)
        raise NotImplementedError(
            "Dask doesn't implement an exact median in all cases as this is hard to do in parallel. "
            "See the `median_approximate` method instead, which uses an approximate algorithm."
        )

    def median_approximate(self, axis=0, method="default", numeric_only=False):
        """Return the approximate median of the values over the requested axis.

        Parameters
        ----------
        axis : {0, 1, "index", "columns"} (default 0)
            0 or ``"index"`` for row-wise, 1 or ``"columns"`` for column-wise
        method : {'default', 'tdigest', 'dask'}, optional
            What method to use. By default will use Dask's internal custom
            algorithm (``"dask"``).  If set to ``"tdigest"`` will use tdigest
            for floats and ints and fallback to the ``"dask"`` otherwise.
        """
        return self.quantile(
            axis=axis, method=method, numeric_only=numeric_only
        ).rename(None)

    @derived_from(pd.DataFrame)
    def describe(
        self,
        split_every=False,
        percentiles=None,
        percentiles_method="default",
        include=None,
        exclude=None,
    ):
        # TODO: duplicated columns
        if include is None and exclude is None:
            _include = [np.number, np.timedelta64, np.datetime64]
            columns = self._meta.select_dtypes(include=_include).columns
            if len(columns) == 0:
                columns = self._meta.columns
        elif include == "all":
            if exclude is not None:
                raise ValueError("exclude must be None when include is 'all'")
            columns = self._meta.columns
        else:
            columns = self._meta.select_dtypes(include=include, exclude=exclude).columns

        stats = [
            self[col].describe(
                split_every=split_every,
                percentiles=percentiles,
                percentiles_method=percentiles_method,
            )
            for col in columns
        ]
        return concat(stats, axis=1)

    @derived_from(pd.DataFrame)
    def pop(self, item):
        out = self[item]
        self._expr = expr.Drop(self, columns=[item])
        return out

    def info(self, buf=None, verbose=False, memory_usage=False):
        """
        Concise summary of a Dask DataFrame
        """
        if buf is None:
            import sys

            buf = sys.stdout
        lines = [str(type(self)).replace("._collection", "")]

        if len(self.columns) == 0:
            lines.append(f"{type(self.index._meta).__name__}: 0 entries")
            lines.append(f"Empty {type(self).__name__}")
            put_lines(buf, lines)
            return

        # Group and execute the required computations
        computations = {}
        if verbose:
            computations.update({"index": self.index, "count": self.count()})
        if memory_usage:
            computations["memory_usage"] = self.memory_usage(deep=True, index=True)

        computations = dict(zip(computations.keys(), compute(*computations.values())))

        if verbose:
            import textwrap

            index = computations["index"]
            counts = computations["count"]
            lines.append(index_summary(index))
            lines.append(f"Data columns (total {len(self.columns)} columns):")

            from pandas.io.formats.printing import pprint_thing

            space = max(len(pprint_thing(k)) for k in self.columns) + 1
            column_width = max(space, 7)

            header = (
                textwrap.dedent(
                    """\
             #   {{column:<{column_width}}} Non-Null Count  Dtype
            ---  {{underl:<{column_width}}} --------------  -----"""
                )
                .format(column_width=column_width)
                .format(column="Column", underl="------")
            )
            column_template = textwrap.dedent(
                """\
            {{i:^3}}  {{name:<{column_width}}} {{count}} non-null      {{dtype}}""".format(
                    column_width=column_width
                )
            )
            column_info = [
                column_template.format(
                    i=pprint_thing(i),
                    name=pprint_thing(name),
                    count=pprint_thing(count),
                    dtype=pprint_thing(dtype),
                )
                for i, (name, count, dtype) in enumerate(
                    # NOTE: Use `counts.values` for cudf support
                    zip(self.columns, counts.values, self.dtypes)
                )
            ]
            lines.extend(header.split("\n"))
        else:
            column_info = [index_summary(self.columns, name="Columns")]

        lines.extend(column_info)
        dtype_counts = [
            "%s(%d)" % k for k in sorted(self.dtypes.value_counts().items(), key=str)
        ]
        lines.append("dtypes: {}".format(", ".join(dtype_counts)))

        if memory_usage:
            memory_int = computations["memory_usage"].sum()
            lines.append(f"memory usage: {memory_repr(memory_int)}\n")

        put_lines(buf, lines)

    @derived_from(pd.DataFrame)
    def cov(self, min_periods=None, numeric_only=False, split_every=False):
        return self._cov(min_periods, numeric_only, split_every)

    @derived_from(pd.DataFrame)
    def corr(
        self,
        method="pearson",
        min_periods=None,
        numeric_only=False,
        split_every=False,
    ):
        return self._corr(method, min_periods, numeric_only, split_every)

    @derived_from(pd.DataFrame)
    def to_string(self, max_rows=5):
        # option_context doesn't affect
        return self._repr_data().to_string(max_rows=max_rows, show_dimensions=False)

    @derived_from(pd.DataFrame)
    def to_html(self, max_rows=5):
        # pd.Series doesn't have html repr
        data = self._repr_data().to_html(max_rows=max_rows, show_dimensions=False)
        n_expr = len({e._name for e in self.walk()})
        return get_template("dataframe.html.j2").render(
            data=data,
            name=self._name,
            layers=maybe_pluralize(n_expr, "expression"),
        )

    @derived_from(pd.DataFrame)
    def melt(
        self,
        id_vars=None,
        value_vars=None,
        var_name=None,
        value_name="value",
        col_level=None,
    ):
        return melt(
            self,
            id_vars=id_vars,
            value_vars=value_vars,
            var_name=var_name,
            value_name=value_name,
            col_level=col_level,
        )

    def _repr_data(self):
        meta = self._meta
        index = self._repr_divisions
        cols = meta.columns
        if len(cols) == 0:
            series_df = pd.DataFrame([[]] * len(index), columns=cols, index=index)
        else:
            series_df = pd.concat(
                [_repr_data_series(s, index=index) for _, s in meta.items()], axis=1
            )
        return series_df


class Series(FrameBase):
    """Series-like Expr Collection.

    The constructor takes the expression that represents the query as input. The class
    is not meant to be instantiated directly. Instead, use one of the IO connectors from
    Dask.
    """

    _accessors: ClassVar[set[str]] = set()
    _partition_type = pd.Series

    @property
    def shape(self):
        """
        Return a tuple representing the dimensionality of the DataFrame.

        The number of rows is a Delayed result. The number of columns
        is a concrete integer.
        """
        return (self.size,)

    @property
    def axes(self):
        return [self.index]

    @property
    def ndim(self):
        """Return dimensionality"""
        return 1

    @property
    def _elemwise(self):
        return elemwise

    def __dir__(self):
        o = set(dir(type(self)))
        o.update(self.__dict__)
        o.update(set(dir(expr.Expr)))
        for accessor in ["cat", "str"]:
            if not hasattr(self._meta, accessor):
                o.remove(accessor)
        return list(o)

    def __contains__(self, item):
        raise NotImplementedError(
            "Using 'in' to test for membership is not supported. Use the values instead"
        )

    @derived_from(pd.Series)
    def __iter__(self):
        frame = self.optimize()
        for i in range(self.npartitions):
            s = frame.get_partition(i).compute()
            yield from s

    def __getitem__(self, key):
        if isinstance(key, Series) or self.npartitions == 1:
            return super().__getitem__(key)
        return self.loc[key]

    @property
    def name(self):
        return self.expr._meta.name

    @name.setter
    def name(self, name):
        self._expr = self.rename(index=name)._expr

    @property
    def dtype(self):
        return self.expr._meta.dtype

    @property
    def nbytes(self):
        """Number of bytes"""
        return new_collection(self.expr.nbytes)

    def __array_ufunc__(self, numpy_ufunc, method, *inputs, **kwargs):
        out = kwargs.get("out", ())
        for x in inputs + out:
            # ufuncs work with 0-dimensional NumPy ndarrays
            # so we don't want to raise NotImplemented
            if isinstance(x, np.ndarray) and x.shape == ():
                continue
            elif not isinstance(
                x, (Number, Scalar, FrameBase, Array, pd.DataFrame, pd.Series, pd.Index)
            ):
                return NotImplemented

        if method == "__call__":
            if numpy_ufunc.signature is not None:
                return NotImplemented
            if numpy_ufunc.nout > 1:
                # ufuncs with multiple output values
                # are not yet supported for frames
                return NotImplemented
            else:
                return elemwise(numpy_ufunc, *inputs, **kwargs)
        else:
            # ufunc methods are not yet supported for frames
            return NotImplemented

    def __array_wrap__(self, array, context=None):
        if isinstance(context, tuple) and len(context) > 0:
            if isinstance(context[1][0], np.ndarray) and context[1][0].shape == ():
                index = None
            else:
                index = context[1][0].index
        else:
            try:
                import inspect

                method_name = f"`{inspect.stack()[3][3]}`"
            except IndexError:
                method_name = "This method"
            raise NotImplementedError(
                f"{method_name} is not implemented for `dask.dataframe.Series`."
            )

        return meta_series_constructor(self)(array, index=index, name=self.name)

    @derived_from(pd.Series)
    def map(self, arg, na_action=None, meta=None):
        if isinstance(arg, Series):
            if not expr.are_co_aligned(self.expr, arg.expr):
                if meta is None:
                    warnings.warn(meta_warning(meta, method="map"))
                return new_collection(
                    expr.MapAlign(self, arg, op=None, na_action=na_action, meta=meta)
                )
        if meta is None:
            meta = expr.emulate(M.map, self, arg, na_action=na_action, udf=True)
            warnings.warn(meta_warning(meta, method="map"))
        return new_collection(expr.Map(self, arg=arg, na_action=na_action, meta=meta))

    @derived_from(pd.Series)
    def clip(self, lower=None, upper=None, axis=None, **kwargs):
        axis = self._validate_axis(axis)
        return new_collection(self.expr.clip(lower, upper, axis))

    @derived_from(pd.Series)
    def to_frame(self, name=no_default):
        return new_collection(expr.ToFrame(self, name=name))

    def _comparison_op(self, expr_cls, other, level, fill_value, axis):
        if level is not None:
            raise NotImplementedError("level must be None")
        self._validate_axis(axis)
        return new_collection(expr_cls(self, other, fill_value=fill_value))

    def lt(self, other, level=None, fill_value=None, axis=0):
        return self._comparison_op(expr.LTSeries, other, level, fill_value, axis)

    def le(self, other, level=None, fill_value=None, axis=0):
        return self._comparison_op(expr.LESeries, other, level, fill_value, axis)

    def gt(self, other, level=None, fill_value=None, axis=0):
        return self._comparison_op(expr.GTSeries, other, level, fill_value, axis)

    def ge(self, other, level=None, fill_value=None, axis=0):
        return self._comparison_op(expr.GESeries, other, level, fill_value, axis)

    def ne(self, other, level=None, fill_value=None, axis=0):
        return self._comparison_op(expr.NESeries, other, level, fill_value, axis)

    def eq(self, other, level=None, fill_value=None, axis=0):
        return self._comparison_op(expr.EQSeries, other, level, fill_value, axis)

    @derived_from(pd.Series)
    def value_counts(
        self,
        sort=None,
        ascending=False,
        dropna=True,
        normalize=False,
        split_every=None,
        split_out=no_default,
    ):
        if split_out is no_default:
            if isinstance(self.dtype, CategoricalDtype):
                # unobserved or huge categories will lead to oom errors
                if self.cat.known:
                    split_out = 1 + len(self.dtype.categories) // 100_000
                else:
                    split_out = True
            else:
                split_out = True
        if split_out == 1 and split_out is not True and sort is None:
            sort = True

        length = None
        if (split_out > 1 or split_out is True) and normalize:
            frame = self if not dropna else self.dropna()
            length = Len(frame)

        return new_collection(
            ValueCounts(
                self, sort, ascending, dropna, normalize, split_every, split_out, length
            )
        )

    @derived_from(pd.Series)
    def mode(self, dropna=True, split_every=False):
        return new_collection(self.expr.mode(dropna, split_every))

    @derived_from(pd.Series)
    def nlargest(self, n=5, split_every=None):
        return new_collection(NLargest(self, n=n, split_every=split_every))

    @derived_from(pd.Series)
    def nsmallest(self, n=5, split_every=None):
        return new_collection(NSmallest(self, n=n, split_every=split_every))

    @derived_from(pd.Series)
    def memory_usage(self, deep=False, index=True):
        return new_collection(MemoryUsageFrame(self, deep=deep, _index=index))

    def unique(self, split_every=None, split_out=True, shuffle_method=None):
        """
        Return Series of unique values in the object. Includes NA values.

        Returns
        -------
        uniques : Series
        """
        shuffle_method = _get_shuffle_preferring_order(shuffle_method)
        return new_collection(Unique(self, split_every, split_out, shuffle_method))

    @derived_from(pd.Series)
    def nunique(self, dropna=True, split_every=False, split_out=True):
        uniqs = self.drop_duplicates(split_every=split_every, split_out=split_out)
        if dropna:
            # count mimics pandas behavior and excludes NA values
            if isinstance(uniqs, Index):
                uniqs = uniqs.to_series()
            return uniqs.count()
        else:
            return uniqs.size

    def drop_duplicates(
        self,
        ignore_index=False,
        split_every=None,
        split_out=True,
        shuffle_method=None,
        keep="first",
    ):
        shuffle_method = _get_shuffle_preferring_order(shuffle_method)
        if keep is False:
            raise NotImplementedError("drop_duplicates with keep=False")
        return new_collection(
            DropDuplicates(
                self,
                ignore_index=ignore_index,
                split_out=split_out,
                split_every=split_every,
                shuffle_method=shuffle_method,
                keep=keep,
            )
        )

    @insert_meta_param_description(pad=12)
    def apply(self, function, *args, meta=no_default, axis=0, **kwargs):
        """Parallel version of pandas.Series.apply

        Parameters
        ----------
        func : function
            Function to apply
        $META
        args : tuple
            Positional arguments to pass to function in addition to the value.

        Additional keyword arguments will be passed as keywords to the function.

        Returns
        -------
        applied : Series or DataFrame if func returns a Series.

        Examples
        --------
        >>> import dask.dataframe as dd
        >>> s = pd.Series(range(5), name='x')
        >>> ds = dd.from_pandas(s, npartitions=2)

        Apply a function elementwise across the Series, passing in extra
        arguments in ``args`` and ``kwargs``:

        >>> def myadd(x, a, b=1):
        ...     return x + a + b
        >>> res = ds.apply(myadd, args=(2,), b=1.5)  # doctest: +SKIP

        By default, dask tries to infer the output metadata by running your
        provided function on some fake data. This works well in many cases, but
        can sometimes be expensive, or even fail. To avoid this, you can
        manually specify the output metadata with the ``meta`` keyword. This
        can be specified in many forms, for more information see
        ``dask.dataframe.utils.make_meta``.

        Here we specify the output is a Series with name ``'x'``, and dtype
        ``float64``:

        >>> res = ds.apply(myadd, args=(2,), b=1.5, meta=('x', 'f8'))

        In the case where the metadata doesn't change, you can also pass in
        the object itself directly:

        >>> res = ds.apply(lambda x: x + 1, meta=ds)

        See Also
        --------
        Series.map_partitions
        """
        self._validate_axis(axis)
        if meta is no_default:
            meta = expr.emulate(M.apply, self, function, args=args, udf=True, **kwargs)
            warnings.warn(meta_warning(meta))
        return new_collection(self.expr.apply(function, *args, meta=meta, **kwargs))

    @classmethod
    def _validate_axis(cls, axis=0, numeric_axis: bool = True) -> None | Literal[0, 1]:
        if axis not in (0, "index", None):
            raise ValueError(f"No axis named {axis} for Series")
        if numeric_axis:
            num_axis: dict[str | None, Literal[0, 1]] = {"index": 0}
            return num_axis.get(axis, axis)
        else:
            return axis

    @derived_from(pd.Series)
    def squeeze(self):
        return self

    @derived_from(pd.Series)
    def dropna(self):
        return new_collection(expr.DropnaSeries(self))

    @derived_from(pd.Series)
    def between(self, left, right, inclusive="both"):
        return new_collection(
            expr.Between(self, left=left, right=right, inclusive=inclusive)
        )

    @derived_from(pd.Series)
    def combine(self, other, func, fill_value=None):
        other = self._create_alignable_frame(other, "outer")
        if not expr.are_co_aligned(self.expr, other.expr):
            return new_collection(
                expr.CombineSeriesAlign(self, other, func, fill_value)
            )
        return new_collection(expr.CombineSeries(self, other, func, fill_value))

    @derived_from(pd.Series)
    def explode(self):
        return new_collection(expr.ExplodeSeries(self))

    @derived_from(pd.Series)
    def add_prefix(self, prefix):
        return new_collection(expr.AddPrefixSeries(self, prefix))

    @derived_from(pd.Series)
    def add_suffix(self, suffix):
        return new_collection(expr.AddSuffixSeries(self, suffix))

    cat = CachedAccessor("cat", CategoricalAccessor)
    dt = CachedAccessor("dt", DatetimeAccessor)
    str = CachedAccessor("str", StringAccessor)

    def _repartition_quantiles(self, npartitions, upsample=1.0, random_state=None):
        return new_collection(
            RepartitionQuantiles(self, npartitions, upsample, random_state)
        )

    @derived_from(pd.Series)
    def groupby(self, by, **kwargs):
        from dask.dataframe.dask_expr._groupby import SeriesGroupBy

        return SeriesGroupBy(self, by, **kwargs)

    def rename(self, index, sorted_index=False):
        """Alter Series index labels or name

        Function / dict values must be unique (1-to-1). Labels not contained in
        a dict / Series will be left as-is. Extra labels listed don't throw an
        error.

        Alternatively, change ``Series.name`` with a scalar value.

        Parameters
        ----------
        index : scalar, hashable sequence, dict-like or callable, optional
            If dict-like or callable, the transformation is applied to the
            index. Scalar or hashable sequence-like will alter the
            ``Series.name`` attribute.
        inplace : boolean, default False
            Whether to return a new Series or modify this one inplace.
        sorted_index : bool, default False
            If true, the output ``Series`` will have known divisions inferred
            from the input series and the transformation. Ignored for
            non-callable/dict-like ``index`` or when the input series has
            unknown divisions. Note that this may only be set to ``True`` if
            you know that the transformed index is monotonically increasing. Dask
            will check that transformed divisions are monotonic, but cannot
            check all the values between divisions, so incorrectly setting this
            can result in bugs.

        Returns
        -------
        renamed : Series

        See Also
        --------
        pandas.Series.rename
        """
        return new_collection(expr.RenameSeries(self, index, sorted_index))

    def quantile(self, q=0.5, method="default"):
        """Approximate quantiles of Series

        Parameters
        ----------
        q : list/array of floats, default 0.5 (50%)
            Iterable of numbers ranging from 0 to 1 for the desired quantiles
        method : {'default', 'tdigest', 'dask'}, optional
            What method to use. By default will use dask's internal custom
            algorithm (``'dask'``).  If set to ``'tdigest'`` will use tdigest
            for floats and ints and fallback to the ``'dask'`` otherwise.
        """
        _raise_if_object_series(self, "quantile")
        allowed_methods = ["default", "dask", "tdigest"]
        if method not in allowed_methods:
            raise ValueError("method can only be 'default', 'dask' or 'tdigest'")
        return new_collection(SeriesQuantile(self, q, method))

    @derived_from(pd.Series)
    def median(self):
        if self.npartitions == 1:
            return self.median_approximate()
        raise NotImplementedError(
            "Dask doesn't implement an exact median in all cases as this is hard to do in parallel. "
            "See the `median_approximate` method instead, which uses an approximate algorithm."
        )

    def median_approximate(self, method="default"):
        """Return the approximate median of the values over the requested axis.

        Parameters
        ----------
        method : {'default', 'tdigest', 'dask'}, optional
            What method to use. By default will use Dask's internal custom
            algorithm (``"dask"``).  If set to ``"tdigest"`` will use tdigest
            for floats and ints and fallback to the ``"dask"`` otherwise.
        """
        return self.quantile(method=method)

    @derived_from(pd.Series)
    def cov(self, other, min_periods=None, split_every=False):
        if not isinstance(other, Series):
            raise TypeError("other must be a dask.dataframe.Series")
        df = concat([self, other], axis=1)
        return df._cov(min_periods, True, split_every, scalar=True)

    @derived_from(pd.Series)
    def corr(self, other, method="pearson", min_periods=None, split_every=False):
        if not isinstance(other, Series):
            raise TypeError("other must be a dask.dataframe.Series")
        df = concat([self, other], axis=1)
        return df._corr(method, min_periods, True, split_every, scalar=True)

    @derived_from(pd.Series)
    def autocorr(self, lag=1, split_every=False):
        if not isinstance(lag, Integral):
            raise TypeError("lag must be an integer")
        return self.corr(self if lag == 0 else self.shift(lag), split_every=split_every)

    @derived_from(pd.Series)
    def describe(
        self,
        split_every=False,
        percentiles=None,
        percentiles_method="default",
        include=None,
        exclude=None,
    ):
        if (
            is_numeric_dtype(self.dtype)
            and not is_bool_dtype(self.dtype)
            or is_timedelta64_dtype(self.dtype)
            or is_datetime64_any_dtype(self.dtype)
        ):
            return new_collection(
                DescribeNumeric(self, split_every, percentiles, percentiles_method)
            )
        else:
            return new_collection(
                DescribeNonNumeric(self, split_every, percentiles, percentiles_method)
            )

    @property
    @derived_from(pd.Series)
    def is_monotonic_increasing(self):
        return new_collection(IsMonotonicIncreasing(self))

    @property
    @derived_from(pd.Series)
    def is_monotonic_decreasing(self):
        return new_collection(IsMonotonicDecreasing(self))

    @derived_from(pd.Series)
    def to_string(self, max_rows=5):
        # option_context doesn't affect
        return self._repr_data().to_string(max_rows=max_rows)

    def _repr_data(self):
        return _repr_data_series(self._meta, self._repr_divisions)

    if PANDAS_GE_220:

        @derived_from(pd.Series)
        def case_when(self, caselist):
            if not isinstance(caselist, list):
                raise TypeError("The caselist argument should be a list")
            caselist = list(flatten([[c, v] for c, v in caselist], container=list))
            return new_collection(expr.CaseWhen(self, *caselist))


for name in [
    "add",
    "sub",
    "mul",
    "div",
    "divide",
    "truediv",
    "floordiv",
    "mod",
    "pow",
    "radd",
    "rsub",
    "rmul",
    "rdiv",
    "rtruediv",
    "rfloordiv",
    "rmod",
    "rpow",
]:
    assert not hasattr(DataFrame, name), name
    setattr(DataFrame, name, _wrap_expr_method_operator(name, DataFrame))

    assert not hasattr(Series, name), name
    setattr(Series, name, _wrap_expr_method_operator(name, Series))


class Index(Series):
    """Index-like Expr Collection.

    The constructor takes the expression that represents the query as input. The class
    is not meant to be instantiated directly. Instead, use one of the IO connectors from
    Dask.
    """

    _accessors: ClassVar[set[str]] = set()
    _partition_type = pd.Index

    _dt_attributes = {
        "nanosecond",
        "microsecond",
        "millisecond",
        "dayofyear",
        "minute",
        "hour",
        "day",
        "dayofweek",
        "second",
        "week",
        "weekday",
        "weekofyear",
        "month",
        "quarter",
        "year",
    }

    _cat_attributes = {
        "known",
        "as_known",
        "as_unknown",
        "add_categories",
        "categories",
        "remove_categories",
        "reorder_categories",
        "as_ordered",
        "codes",
        "remove_unused_categories",
        "set_categories",
        "as_unordered",
        "ordered",
        "rename_categories",
    }

    def __getattr__(self, key):
        if (
            isinstance(self._meta.dtype, pd.CategoricalDtype)
            and key in self._cat_attributes
        ):
            return getattr(self.cat, key)
        elif key in self._dt_attributes:
            return getattr(self.dt, key)

        if hasattr(super(), key):  # Doesn't trigger super().__getattr__
            # Not a magic attribute. This is a real method or property of Series that
            # has been overridden by RaiseAttributeError().
            raise AttributeError(
                f"{self.__class__.__name__!r} object has no attribute {key!r}"
            )
        return super().__getattr__(key)

    def __repr__(self):
        return f"<dask_expr.expr.Index: expr={self.expr}>"

    def __array_wrap__(self, array, context=None):
        return pd.Index(array, name=self.name)

    @derived_from(pd.Index)
    def to_series(self, index=None, name=no_default):
        if index is not None:
            raise NotImplementedError
        return new_collection(expr.ToSeriesIndex(self, index=index, name=name))

    @derived_from(pd.Index, ua_args=["index"])
    def to_frame(self, index=True, name=no_default):
        if not index:
            raise NotImplementedError
        return new_collection(expr.ToFrameIndex(self, index=index, name=name))

    @derived_from(pd.Index)
    def memory_usage(self, deep=False):
        return new_collection(MemoryUsageIndex(self, deep=deep))

    def shift(self, periods=1, freq=None):
        return new_collection(expr.ShiftIndex(self, periods, freq))

    @derived_from(pd.Index)
    def map(self, arg, na_action=None, meta=None, is_monotonic=False):
        """
        Note that this method clears any known divisions.

        If your mapping function is monotonically increasing then use `is_monotonic`
        to apply the maping function to the old divisions and assign the new
        divisions to the output.

        """
        if isinstance(arg, Series):
            if not expr.are_co_aligned(self.expr, arg.expr):
                if meta is None:
                    warnings.warn(meta_warning(meta, method="map"))
                return new_collection(
                    expr.MapIndexAlign(self, arg, na_action, meta, is_monotonic)
                )
        if meta is None:
            meta = expr.emulate(M.map, self, arg, na_action=na_action, udf=True)
            warnings.warn(meta_warning(meta, method="map"))
        return new_collection(
            expr.Map(
                self, arg=arg, na_action=na_action, meta=meta, is_monotonic=is_monotonic
            )
        )

    def __dir__(self):
        o = set(dir(type(self)))
        o.update(self.__dict__)
        o.update(set(dir(expr.Expr)))
        o.update(self._dt_attributes)
        if isinstance(self.dtype, pd.CategoricalDtype):
            o.update(self._cat_attributes)
        return list(o)

    # Methods and properties of Series that are not implemented on Index

    def count(self, split_every=False):
        return new_collection(IndexCount(self, split_every))

    @property  # type: ignore
    def index(self):
        raise AttributeError("'Index' object has no attribute 'index'")

    def sum(self, *args, **kwargs):
        raise AttributeError("'Index' object has no attribute 'sum'")

    def prod(self, *args, **kwargs):
        raise AttributeError("'Index' object has no attribute 'prod'")

    def mean(self, *args, **kwargs):
        raise AttributeError("'Index' object has no attribute 'mean'")

    def std(self, *args, **kwargs):
        raise AttributeError("'Index' object has no attribute 'std'")

    def var(self, *args, **kwargs):
        raise AttributeError("'Index' object has no attribute 'var'")

    def idxmax(self, *args, **kwargs):
        raise AttributeError("'Index' object has no attribute 'idxmax'")

    def idxmin(self, *args, **kwargs):
        raise AttributeError("'Index' object has no attribute 'idxmin'")


class Scalar(FrameBase):
    """Scalar Expr Collection"""

    def __repr__(self):
        return f"<dask_expr.expr.Scalar: expr={self.expr}, dtype={self.dtype}>"

    def __bool__(self):
        raise TypeError(
            f"Trying to convert {self} to a boolean value. Because Dask objects are "
            "lazily evaluated, they cannot be converted to a boolean value or used "
            "in boolean conditions like if statements. Try calling .compute() to "
            "force computation prior to converting to a boolean value or using in "
            "a conditional statement."
        )

    def __dask_postcompute__(self):
        return first, ()

    def to_series(self, index=0) -> Series:
        return new_collection(expr.ScalarToSeries(self, index=index))

    def __array__(self):
        # array interface is required to support pandas instance + Scalar
        # Otherwise, above op results in pd.Series of Scalar (object dtype)
        return np.asarray(self.compute())

    @functools.cached_property
    def dtype(self):
        return pd.Series(self._meta).dtype

    def to_delayed(self, optimize_graph=True):
        return super().to_delayed(optimize_graph=optimize_graph)[0]


def optimize(collection, fuse=True):
    return new_collection(expr.optimize(collection.expr, fuse=fuse))


def from_pandas(data, npartitions=None, sort=True, chunksize=None):
    """
    Construct a Dask DataFrame from a Pandas DataFrame

    This splits an in-memory Pandas dataframe into several parts and constructs
    a dask.dataframe from those parts on which Dask.dataframe can operate in
    parallel.  By default, the input dataframe will be sorted by the index to
    produce cleanly-divided partitions (with known divisions).  To preserve the
    input ordering, make sure the input index is monotonically-increasing. The
    ``sort=False`` option will also avoid reordering, but will not result in
    known divisions.

    Parameters
    ----------
    data : pandas.DataFrame or pandas.Series
        The DataFrame/Series with which to construct a Dask DataFrame/Series
    npartitions : int, optional, default 1
        The number of partitions of the index to create. Note that if there
        are duplicate values or insufficient elements in ``data.index``, the
        output may have fewer partitions than requested.
    chunksize : int, optional
        The desired number of rows per index partition to use. Note that
        depending on the size and index of the dataframe, actual partition
        sizes may vary.
    sort: bool, default True
        Sort the input by index first to obtain cleanly divided partitions
        (with known divisions).  If False, the input will not be sorted, and
        all divisions will be set to None. Default is True.

    Returns
    -------
    dask.DataFrame or dask.Series
        A dask DataFrame/Series partitioned along the index

    Examples
    --------
    >>> from dask.dataframe import from_pandas
    >>> df = pd.DataFrame(dict(a=list('aabbcc'), b=list(range(6))),
    ...                   index=pd.date_range(start='20100101', periods=6))
    >>> ddf = from_pandas(df, npartitions=3)
    >>> ddf.divisions  # doctest: +NORMALIZE_WHITESPACE
    (Timestamp('2010-01-01 00:00:00'),
     Timestamp('2010-01-03 00:00:00'),
     Timestamp('2010-01-05 00:00:00'),
     Timestamp('2010-01-06 00:00:00'))
    >>> ddf = from_pandas(df.a, npartitions=3)  # Works with Series too!
    >>> ddf.divisions  # doctest: +NORMALIZE_WHITESPACE
    (Timestamp('2010-01-01 00:00:00'),
     Timestamp('2010-01-03 00:00:00'),
     Timestamp('2010-01-05 00:00:00'),
     Timestamp('2010-01-06 00:00:00'))

    Raises
    ------
    TypeError
        If something other than a ``pandas.DataFrame`` or ``pandas.Series`` is
        passed in.

    See Also
    --------
    from_array : Construct a dask.DataFrame from an array that has record dtype
    read_csv : Construct a dask.DataFrame from a CSV file
    """
    if chunksize is not None and npartitions is not None:
        raise ValueError("Exactly one of npartitions and chunksize must be specified.")
    elif chunksize is None and npartitions is None:
        npartitions = 1

    if not has_parallel_type(data):
        raise TypeError("Input must be a pandas DataFrame or Series.")

    if data.index.isna().any() and not _is_any_real_numeric_dtype(data.index):
        raise NotImplementedError(
            "Index in passed data is non-numeric and contains nulls, which Dask does not entirely support.\n"
            "Consider passing `data.loc[~data.isna()]` instead."
        )

    if npartitions is not None and not isinstance(npartitions, int):
        raise TypeError(
            "Please provide npartitions as an int, or possibly as None if you specify chunksize."
        )
    elif chunksize is not None and not isinstance(chunksize, int):
        raise TypeError(
            "Please provide chunksize as an int, or possibly as None if you specify npartitions."
        )

    from dask.dataframe.dask_expr.io.io import FromPandas

    return new_collection(
        FromPandas(
            _BackendData(data.copy()),
            npartitions=npartitions,
            sort=sort,
            chunksize=chunksize,
            pyarrow_strings_enabled=pyarrow_strings_enabled(),
        )
    )


def from_array(arr, chunksize=50_000, columns=None, meta=None):
    """Read any sliceable array into a Dask Dataframe

    Uses getitem syntax to pull slices out of the array.  The array need not be
    a NumPy array but must support slicing syntax

        x[50000:100000]

    and have 2 dimensions:

        x.ndim == 2

    or have a record dtype:

        x.dtype == [('name', 'O'), ('balance', 'i8')]

    Parameters
    ----------
    x : array_like
    chunksize : int, optional
        The number of rows per partition to use.
    columns : list or string, optional
        list of column names if DataFrame, single string if Series
    meta : object, optional
        An optional `meta` parameter can be passed for dask
        to specify the concrete dataframe type to use for partitions of
        the Dask dataframe. By default, pandas DataFrame is used.

    Returns
    -------
    dask.DataFrame or dask.Series
        A dask DataFrame/Series
    """
    import dask.array as da

    if isinstance(arr, da.Array):
        return from_dask_array(arr, columns=columns, meta=meta)

    from dask.dataframe.dask_expr.io.io import FromArray

    result = FromArray(
        arr,
        chunksize=chunksize,
        original_columns=columns,
        meta=meta,
    )
    if pyarrow_strings_enabled() and arr.dtype.kind in "OU":
        result = expr.ArrowStringConversion(result)
    return new_collection(result)


def from_graph(layer, _meta, divisions, keys, name_prefix):
    from dask.dataframe.dask_expr.io.io import FromGraph

    return new_collection(
        FromGraph(
            layer=layer,
            _meta=_meta,
            divisions=divisions,
            keys=keys,
            name_prefix=name_prefix,
        )
    )


@dataframe_creation_dispatch.register_inplace("pandas")
def from_dict(
    data,
    npartitions,
    orient="columns",
    dtype=None,
    columns=None,
    constructor=pd.DataFrame,
):
    """
    Construct a Dask DataFrame from a Python Dictionary

    Parameters
    ----------
    data : dict
        Of the form {field : array-like} or {field : dict}.
    npartitions : int
        The number of partitions of the index to create. Note that depending on
        the size and index of the dataframe, the output may have fewer
        partitions than requested.
    orient : {'columns', 'index', 'tight'}, default 'columns'
        The "orientation" of the data. If the keys of the passed dict
        should be the columns of the resulting DataFrame, pass 'columns'
        (default). Otherwise if the keys should be rows, pass 'index'.
        If 'tight', assume a dict with keys
        ['index', 'columns', 'data', 'index_names', 'column_names'].
    dtype: bool
        Data type to force, otherwise infer.
    columns: string, optional
        Column labels to use when ``orient='index'``. Raises a ValueError
        if used with ``orient='columns'`` or ``orient='tight'``.
    constructor: class, default pd.DataFrame
        Class with which ``from_dict`` should be called with.

    Examples
    --------
    >>> import dask.dataframe as dd
    >>> ddf = dd.from_dict({"num1": [1, 2, 3, 4], "num2": [7, 8, 9, 10]}, npartitions=2)
    """

    collection_types = {type(v) for v in data.values() if is_dask_collection(v)}
    if collection_types:
        raise NotImplementedError(
            "from_dict doesn't currently support Dask collections as inputs. "
            f"Objects of type {collection_types} were given in the input dict."
        )

    return from_pandas(
        constructor.from_dict(data, orient, dtype, columns),
        npartitions,
    )


def from_dask_array(x, columns=None, index=None, meta=None) -> DataFrame:
    """Create a Dask DataFrame from a Dask Array.

    Converts a 2d array into a DataFrame and a 1d array into a Series.

    Parameters
    ----------
    x : da.Array
    columns : list or string
        list of column names if DataFrame, single string if Series
    index : dask.dataframe.Index, optional
        An optional *dask* Index to use for the output Series or DataFrame.

        The default output index depends on whether `x` has any unknown
        chunks. If there are any unknown chunks, the output has ``None``
        for all the divisions (one per chunk). If all the chunks are known,
        a default index with known divisions is created.

        Specifying `index` can be useful if you're conforming a Dask Array
        to an existing dask Series or DataFrame, and you would like the
        indices to match.
    meta : object, optional
        An optional `meta` parameter can be passed for dask
        to specify the concrete dataframe type to be returned.
        By default, pandas DataFrame is used.

    Examples
    --------
    >>> import dask.array as da
    >>> import dask.dataframe as dd
    >>> x = da.ones((4, 2), chunks=(2, 2))
    >>> df = dd.io.from_dask_array(x, columns=['a', 'b'])
    >>> df.compute()
         a    b
    0  1.0  1.0
    1  1.0  1.0
    2  1.0  1.0
    3  1.0  1.0

    See Also
    --------
    dask.bag.to_dataframe: from dask.bag
    dask.dataframe.DataFrame.values: Reverse conversion
    dask.dataframe.DataFrame.to_records: Reverse conversion
    """
    from dask.dataframe.io import from_dask_array

    if columns is not None and isinstance(columns, list) and not len(columns):
        columns = None
    return from_dask_array(x, columns=columns, index=index, meta=meta)


@dataframe_creation_dispatch.register_inplace("pandas")
def read_parquet(
    path=None,
    columns=None,
    filters=None,
    categories=None,
    index=None,
    storage_options=None,
    dtype_backend=None,
    calculate_divisions=False,
    ignore_metadata_file=False,
    metadata_task_size=None,
    split_row_groups="infer",
    blocksize="default",
    aggregate_files=None,
    parquet_file_extension=(".parq", ".parquet", ".pq"),
    filesystem="fsspec",
    engine=None,
    arrow_to_pandas=None,
    **kwargs,
):
    """
    Read a Parquet file into a Dask DataFrame

    This reads a directory of Parquet data into a Dask.dataframe, one file per
    partition.  It selects the index among the sorted columns if any exist.

    .. note::
        Dask automatically resizes partitions to ensure that each partition is of
        adequate size. The optimizer uses the ratio of selected columns to total
        columns to squash multiple files into one partition.

        Additionally, the Optimizer uses a minimum size per partition (default 75MB)
        to avoid too many small partitions. This configuration can be set with

        >>> dask.config.set({"dataframe.parquet.minimum-partition-size": "100MB"})  # doctest: +SKIP

    .. note::
        Specifying ``filesystem="arrow"`` leverages a complete reimplementation of
        the Parquet reader that is solely based on PyArrow. It is significantly faster
        than the legacy implementation, but doesn't yet support all features.

    Parameters
    ----------
    path : str or list
        Source directory for data, or path(s) to individual parquet files.
        Prefix with a protocol like ``s3://`` to read from alternative
        filesystems. To read from multiple files you can pass a globstring or a
        list of paths, with the caveat that they must all have the same
        protocol.
    columns : str or list, default None
        Field name(s) to read in as columns in the output. By default all
        non-index fields will be read (as determined by the pandas parquet
        metadata, if present). Provide a single field name instead of a list to
        read in the data as a Series.
    filters : Union[List[Tuple[str, str, Any]], List[List[Tuple[str, str, Any]]]], default None
        List of filters to apply, like ``[[('col1', '==', 0), ...], ...]``.
        Using this argument will result in row-wise filtering of the final partitions.

        Predicates can be expressed in disjunctive normal form (DNF). This means that
        the inner-most tuple describes a single column predicate. These inner predicates
        are combined with an AND conjunction into a larger predicate. The outer-most
        list then combines all of the combined filters with an OR disjunction.

        Predicates can also be expressed as a ``List[Tuple]``. These are evaluated
        as an AND conjunction. To express OR in predicates, one must use the
        (preferred for "pyarrow") ``List[List[Tuple]]`` notation.
    index : str, list or False, default None
        Field name(s) to use as the output frame index. By default will be
        inferred from the pandas parquet file metadata, if present. Use ``False``
        to read all fields as columns.
    categories : list or dict, default None
        For any fields listed here, if the parquet encoding is Dictionary,
        the column will be created with dtype category. Use only if it is
        guaranteed that the column is encoded as dictionary in all row-groups.
        If a list, assumes up to 2**16-1 labels; if a dict, specify the number
        of labels expected; if None, will load categories automatically for
        data written by dask, not otherwise.
    storage_options : dict, default None
        Key/value pairs to be passed on to the file-system backend, if any.
        Note that the default file-system backend can be configured with the
        ``filesystem`` argument, described below.
    open_file_options : dict, default None
        Key/value arguments to be passed along to ``AbstractFileSystem.open``
        when each parquet data file is open for reading. Experimental
        (optimized) "precaching" for remote file systems (e.g. S3, GCS) can
        be enabled by adding ``{"method": "parquet"}`` under the
        ``"precache_options"`` key. Also, a custom file-open function can be
        used (instead of ``AbstractFileSystem.open``), by specifying the
        desired function under the ``"open_file_func"`` key.
    dtype_backend : {'numpy_nullable', 'pyarrow'}, defaults to NumPy backed DataFrames
        Which dtype_backend to use, e.g. whether a DataFrame should have NumPy arrays,
        nullable dtypes are used for all dtypes that have a nullable implementation
        when 'numpy_nullable' is set, pyarrow is used for all dtypes if 'pyarrow'
        is set.
        ``dtype_backend="pyarrow"`` requires ``pandas`` 1.5+.
    calculate_divisions : bool, default False
        Whether to use min/max statistics from the footer metadata (or global
        ``_metadata`` file) to calculate divisions for the output DataFrame
        collection. Divisions will not be calculated if statistics are missing.
        This option will be ignored if ``index`` is not specified and there is
        no physical index column specified in the custom "pandas" Parquet
        metadata. Note that ``calculate_divisions=True`` may be extremely slow
        when no global ``_metadata`` file is present, especially when reading
        from remote storage. Set this to ``True`` only when known divisions
        are needed for your workload (see :ref:`dataframe-design-partitions`).
    ignore_metadata_file : bool, default False
        Whether to ignore the global ``_metadata`` file (when one is present).
        If ``True``, or if the global ``_metadata`` file is missing, the parquet
        metadata may be gathered and processed in parallel. Parallel metadata
        processing is currently supported for ``ArrowDatasetEngine`` only.
    metadata_task_size : int, default configurable
        If parquet metadata is processed in parallel (see ``ignore_metadata_file``
        description above), this argument can be used to specify the number of
        dataset files to be processed by each task in the Dask graph.  If this
        argument is set to ``0``, parallel metadata processing will be disabled.
        The default values for local and remote filesystems can be specified
        with the "metadata-task-size-local" and "metadata-task-size-remote"
        config fields, respectively (see "dataframe.parquet").
    split_row_groups : 'infer', 'adaptive', bool, or int, default 'infer'
        If True, then each output dataframe partition will correspond to a single
        parquet-file row-group. If False, each partition will correspond to a
        complete file.  If a positive integer value is given, each dataframe
        partition will correspond to that number of parquet row-groups (or fewer).
        If 'adaptive', the metadata of each file will be used to ensure that every
        partition satisfies ``blocksize``. If 'infer' (the default), the
        uncompressed storage-size metadata in the first file will be used to
        automatically set ``split_row_groups`` to either 'adaptive' or ``False``.
    blocksize : int or str, default 'default'
        The desired size of each output ``DataFrame`` partition in terms of total
        (uncompressed) parquet storage space. This argument is currently used to
        set the default value of ``split_row_groups`` (using row-group metadata
        from a single file), and will be ignored if ``split_row_groups`` is not
        set to 'infer' or 'adaptive'. Default is 256 MiB.
    aggregate_files : bool or str, default None
        WARNING: Passing a string argument to ``aggregate_files`` will result
        in experimental behavior. This behavior may change in the future.

        Whether distinct file paths may be aggregated into the same output
        partition. This parameter is only used when `split_row_groups` is set to
        'infer', 'adaptive' or to an integer >1. A setting of True means that any
        two file paths may be aggregated into the same output partition, while
        False means that inter-file aggregation is prohibited.

        For "hive-partitioned" datasets, a "partition"-column name can also be
        specified. In this case, we allow the aggregation of any two files
        sharing a file path up to, and including, the corresponding directory name.
        For example, if ``aggregate_files`` is set to ``"section"`` for the
        directory structure below, ``03.parquet`` and ``04.parquet`` may be
        aggregated together, but ``01.parquet`` and ``02.parquet`` cannot be.
        If, however, ``aggregate_files`` is set to ``"region"``, ``01.parquet``
        may be aggregated with ``02.parquet``, and ``03.parquet`` may be aggregated
        with ``04.parquet``::

            dataset-path/
            ├── region=1/
            │   ├── section=a/
            │   │   └── 01.parquet
            │   ├── section=b/
            │   └── └── 02.parquet
            └── region=2/
                ├── section=a/
                │   ├── 03.parquet
                └── └── 04.parquet

        Note that the default behavior of ``aggregate_files`` is ``False``.
    parquet_file_extension: str, tuple[str], or None, default (".parq", ".parquet", ".pq")
        A file extension or an iterable of extensions to use when discovering
        parquet files in a directory. Files that don't match these extensions
        will be ignored. This argument only applies when ``paths`` corresponds
        to a directory and no ``_metadata`` file is present (or
        ``ignore_metadata_file=True``). Passing in ``parquet_file_extension=None``
        will treat all files in the directory as parquet files.

        The purpose of this argument is to ensure that the engine will ignore
        unsupported metadata files (like Spark's '_SUCCESS' and 'crc' files).
        It may be necessary to change this argument if the data files in your
        parquet dataset do not end in ".parq", ".parquet", or ".pq".
    filesystem: "fsspec", "arrow", or fsspec.AbstractFileSystem backend to use.
    dataset: dict, default None
        Dictionary of options to use when creating a ``pyarrow.dataset.Dataset`` object.
        These options may include a "filesystem" key to configure the desired
        file-system backend. However, the top-level ``filesystem`` argument will always
        take precedence.

        **Note**: The ``dataset`` options may include a "partitioning" key.
        However, since ``pyarrow.dataset.Partitioning``
        objects cannot be serialized, the value can be a dict of key-word
        arguments for the ``pyarrow.dataset.partitioning`` API
        (e.g. ``dataset={"partitioning": {"flavor": "hive", "schema": ...}}``).
        Note that partitioned columns will not be converted to categorical
        dtypes when a custom partitioning schema is specified in this way.
    read: dict, default None
        Dictionary of options to pass through to ``engine.read_partitions``
        using the ``read`` key-word argument.
    arrow_to_pandas: dict, default None
        Dictionary of options to use when converting from ``pyarrow.Table`` to
        a pandas ``DataFrame`` object. Only used by the "arrow" engine.
    **kwargs: dict (of dicts)
        Options to pass through to ``engine.read_partitions`` as stand-alone
        key-word arguments. Note that these options will be ignored by the
        engines defined in ``dask.dataframe``, but may be used by other custom
        implementations.

    Examples
    --------
    >>> df = dd.read_parquet('s3://bucket/my-parquet-data')  # doctest: +SKIP

    See Also
    --------
    to_parquet
    pyarrow.parquet.ParquetDataset
    """
    from dask.dataframe.dask_expr.io.parquet import (
        ReadParquetFSSpec,
        ReadParquetPyarrowFS,
        _set_parquet_engine,
    )

    if not isinstance(path, str):
        path = stringify_path(path)

    kwargs["dtype_backend"] = dtype_backend
    if arrow_to_pandas:
        kwargs["arrow_to_pandas"] = arrow_to_pandas

    if filters is not None:
        for filter in flatten(filters, container=list):
            col, op, val = filter
            if op == "in" and not isinstance(val, (set, list, tuple)):
                raise TypeError("Value of 'in' filter must be a list, set or tuple.")

    if (
        isinstance(filesystem, pa_fs.FileSystem)
        or isinstance(filesystem, str)
        and filesystem.lower() in ("arrow", "pyarrow")
    ):
        if parse_version(pa.__version__) < parse_version("15.0.0"):
            raise ValueError(
                "pyarrow>=15.0.0 is required to use the pyarrow filesystem."
            )
        if metadata_task_size is not None:
            raise NotImplementedError(
                "metadata_task_size is not supported when using the pyarrow filesystem."
            )
        if split_row_groups != "infer":
            raise NotImplementedError(
                "split_row_groups is not supported when using the pyarrow filesystem."
            )
        if blocksize is not None and blocksize != "default":
            raise NotImplementedError(
                "blocksize is not supported when using the pyarrow filesystem."
            )
        if aggregate_files is not None:
            raise NotImplementedError(
                "aggregate_files is not supported when using the pyarrow filesystem."
            )
        if parquet_file_extension != (".parq", ".parquet", ".pq"):
            raise NotImplementedError(
                "parquet_file_extension is not supported when using the pyarrow filesystem."
            )
        if engine is not None:
            raise NotImplementedError(
                "engine is not supported when using the pyarrow filesystem."
            )

        return new_collection(
            ReadParquetPyarrowFS(
                path,
                columns=_convert_to_list(columns),
                filters=filters,
                categories=categories,
                index=index,
                calculate_divisions=calculate_divisions,
                storage_options=storage_options,
                filesystem=filesystem,
                ignore_metadata_file=ignore_metadata_file,
                arrow_to_pandas=arrow_to_pandas,
                pyarrow_strings_enabled=pyarrow_strings_enabled(),
                kwargs=kwargs,
                _series=isinstance(columns, str),
            )
        )

    return new_collection(
        ReadParquetFSSpec(
            path,
            columns=_convert_to_list(columns),
            filters=filters,
            categories=categories,
            index=index,
            storage_options=storage_options,
            calculate_divisions=calculate_divisions,
            ignore_metadata_file=ignore_metadata_file,
            metadata_task_size=metadata_task_size,
            split_row_groups=split_row_groups,
            blocksize=blocksize,
            aggregate_files=aggregate_files,
            parquet_file_extension=parquet_file_extension,
            filesystem=filesystem,
            engine=_set_parquet_engine(engine),
            kwargs=kwargs,
            _series=isinstance(columns, str),
        )
    )


def concat(
    dfs,
    axis=0,
    join="outer",
    ignore_unknown_divisions=False,
    ignore_order=False,
    interleave_partitions=False,
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
        return dfs[0]

    if join not in ("inner", "outer"):
        raise ValueError("'join' must be 'inner' or 'outer'")

    dfs = [from_pandas(df) if not isinstance(df, FrameBase) else df for df in dfs]

    if axis == 1:
        dfs = [df for df in dfs if len(df.columns) > 0 or isinstance(df, Series)]

    return new_collection(
        Concat(
            join,
            ignore_order,
            axis,
            ignore_unknown_divisions,
            interleave_partitions,
            kwargs,
            *dfs,
        )
    )


def melt(
    frame,
    id_vars=None,
    value_vars=None,
    var_name=None,
    value_name="value",
    col_level=None,
):
    return map_partitions(
        M.melt,
        frame,
        clear_divisions=True,
        id_vars=id_vars,
        value_vars=value_vars,
        var_name=var_name,
        value_name=value_name,
        col_level=col_level,
    )


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
    shuffle_method=None,
    npartitions=None,
    broadcast=None,
):
    for o in [on, left_on, right_on]:
        if isinstance(o, FrameBase):
            raise NotImplementedError()
    if not on and not left_on and not right_on and not left_index and not right_index:
        on = [c for c in left.columns if c in right.columns]
        if not on:
            left_index = right_index = True

    if on and not left_on and not right_on:
        left_on = right_on = on

    if pd.api.types.is_list_like(left_on) and not isinstance(left_on, FrameBase):
        left_on = list(left_on)
    if pd.api.types.is_list_like(right_on) and not isinstance(right_on, FrameBase):
        right_on = list(right_on)

    supported_how = ("left", "right", "outer", "inner", "leftsemi")
    if how not in supported_how:
        raise ValueError(
            f"dask.dataframe.merge does not support how='{how}'."
            f"Options are: {supported_how}."
        )

    if how == "leftsemi":
        if right_index or any(
            o not in right.columns for o in _convert_to_list(right_on)
        ):
            raise NotImplementedError(
                "how='leftsemi' does not support right_index=True or on columns from the index"
            )
        else:
            right = right[_convert_to_list(right_on)].rename(
                columns=dict(zip(_convert_to_list(right_on), _convert_to_list(left_on)))
            )
            right_on = left_on

    # Transform pandas objects into dask.dataframe objects
    if not is_dask_collection(left):
        if right_index and left_on:  # change to join on index
            left = left.set_index(left[left_on])
            left_on = None
            left_index = True
        left = from_pandas(left, npartitions=1)

    if not is_dask_collection(right):
        if left_index and right_on:  # change to join on index
            right = right.set_index(right[right_on])
            right_on = None
            right_index = True
        right = from_pandas(right, npartitions=1)

    assert is_dataframe_like(right._meta)
    if left_on and right_on:
        warn_dtype_mismatch(left, right, left_on, right_on)

    result = new_collection(
        Merge(
            left,
            right,
            how=how,
            left_on=left_on,
            right_on=right_on,
            left_index=left_index,
            right_index=right_index,
            suffixes=suffixes,
            indicator=indicator,
            shuffle_method=get_specified_shuffle(shuffle_method),
            _npartitions=npartitions,
            broadcast=broadcast,
        )
    )
    if left._meta.index.name != right._meta.index.name:
        return result.rename_axis(index=result._meta.index.name)
    else:
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
        kwargs["left_on"] = left_on
        kwargs["right_on"] = right_on
    del kwargs["on"]

    for o in [left_on, right_on]:
        if isinstance(o, FrameBase):
            raise NotImplementedError(
                "Dask collections not currently allowed in merge columns"
            )

    if not is_dask_collection(left):
        left = from_pandas(left, npartitions=1)

    if not is_dask_collection(right):
        right = from_pandas(right, npartitions=1)

    if by is not None:
        if left_by is not None or right_by is not None:
            raise ValueError(
                "Can only pass argument 'by' OR 'left_by' and 'right_by', not a combination of both."
            )
        kwargs["left_by"] = kwargs["right_by"] = by
    del kwargs["by"]

    if left_by is None and right_by is not None:
        raise ValueError("Must specify both left_on and right_on if one is specified.")
    if left_by is not None and right_by is None:
        raise ValueError("Must specify both left_on and right_on if one is specified.")

    from dask.dataframe.dask_expr._merge_asof import MergeAsof

    return new_collection(MergeAsof(left, right, **kwargs))


def from_map(
    func,
    *iterables,
    args=None,
    meta=no_default,
    divisions=None,
    label=None,
    enforce_metadata=False,
    **kwargs,
):
    """Create a DataFrame collection from a custom function map.

    ``from_map`` is the preferred option when reading from data sources
    that are not natively supported by Dask or if the data source
    requires custom handling before handing things of to Dask DataFrames.
    Examples are things like binary files or other unstructured data that
    doesn't have an IO connector.

    ``from_map`` supports column projection by the optimizer. The optimizer
    tries to push column selections into the from_map call if the function
    supports a ``columns`` argument.

    Parameters
    ----------
    func : callable
        Function used to create each partition. Column projection will be
        enabled if the function has a ``columns`` keyword argument.
    *iterables : Iterable objects
        Iterable objects to map to each output partition. All iterables must
        be the same length. This length determines the number of partitions
        in the output collection (only one element of each iterable will
        be passed to ``func`` for each partition).
    args : list or tuple, optional
        Positional arguments to broadcast to each output partition. Note
        that these arguments will always be passed to ``func`` after the
        ``iterables`` positional arguments.
    $META
    divisions : tuple, str, optional
        Partition boundaries along the index.
        For tuple, see https://docs.dask.org/en/latest/dataframe-design.html#partitions
        For string 'sorted' will compute the delayed values to find index
        values.  Assumes that the indexes are mutually sorted.
        If None, then won't use index information
    label : str, optional
        String to use as the function-name label in the output
        collection-key names.
    token : str, optional
        String to use as the "token" in the output collection-key names.
    enforce_metadata : bool, default True
        Whether to enforce at runtime that the structure of the DataFrame
        produced by ``func`` actually matches the structure of ``meta``.
        This will rename and reorder columns for each partition,
        and will raise an error if this doesn't work,
        but it won't raise if dtypes don't match.
    **kwargs:
        Key-word arguments to broadcast to each output partition. These
        same arguments will be passed to ``func`` for every output partition.

    Examples
    --------
    >>> import pandas as pd
    >>> import dask.dataframe as dd
    >>> func = lambda x, size=0: pd.Series([x] * size)
    >>> inputs = ["A", "B"]
    >>> dd.from_map(func, inputs, size=2).compute()
    0    A
    1    A
    0    B
    1    B
    dtype: string

    The optimizer will identify a column selection that happens after from_map
    and push the columns argument into the actual map call to drop unnecessary
    columns as early as possible.

    >>> def map_function(x, columns=None):
    ...     df = pd.DataFrame({"a": [1, 2], "b": x})
    ...     if columns is not None:
    ...         df = df[columns]
    ...     return df
    >>> dd.from_map(map_function, [1, 2])["b"].compute()
    0    1
    1    1
    0    2
    1    2
    Name: b, dtype: int64

    This API can also be used as an alternative to other file-based
    IO functions, like ``read_csv`` (which are already just
    ``from_map`` wrapper functions):

    >>> import pandas as pd
    >>> import dask.dataframe as dd
    >>> paths = ["0.csv", "1.csv", "2.csv"]
    >>> dd.from_map(pd.read_csv, paths).head()  # doctest: +SKIP
                        name
    timestamp
    2000-01-01 00:00:00   Laura
    2000-01-01 00:00:01  Oliver
    2000-01-01 00:00:02   Alice
    2000-01-01 00:00:03  Victor
    2000-01-01 00:00:04     Bob

    Since ``from_map`` allows you to map an arbitrary function
    to any number of iterable objects, it can be a very convenient
    means of implementing functionality that may be missing
    from other DataFrame-creation methods. For example, if you
    happen to have apriori knowledge about the number of rows
    in each of the files in a dataset, you can generate a
    DataFrame collection with a global RangeIndex:

    >>> import pandas as pd
    >>> import numpy as np
    >>> import dask.dataframe as dd
    >>> paths = ["0.csv", "1.csv", "2.csv"]
    >>> file_sizes = [86400, 86400, 86400]
    >>> def func(path, row_offset):
    ...     # Read parquet file and set RangeIndex offset
    ...     df = pd.read_csv(path)
    ...     return df.set_index(
    ...         pd.RangeIndex(row_offset, row_offset+len(df))
    ...     )
    >>> def get_ddf(paths, file_sizes):
    ...     offsets = [0] + list(np.cumsum(file_sizes))
    ...     return dd.from_map(
    ...         func, paths, offsets[:-1], divisions=offsets
    ...     )
    >>> ddf = get_ddf(paths, file_sizes)  # doctest: +SKIP
    >>> ddf.index  # doctest: +SKIP
    Dask Index Structure:
    npartitions=3
    0         int64
    86400       ...
    172800      ...
    259200      ...
    dtype: int64
    Dask Name: myfunc, 6 tasks
    """
    from dask.dataframe.dask_expr.io import FromMap, FromMapProjectable
    from dask.dataframe.io.utils import DataFrameIOFunction

    if "token" in kwargs:
        # This option doesn't really make sense in dask-expr
        raise NotImplementedError("dask_expr does not support a token argument.")

    lengths = set()
    iterables = list(iterables)
    for i, iterable in enumerate(iterables):
        if not isinstance(iterable, Iterable):
            raise ValueError(
                f"All elements of `iterables` must be Iterable, got {type(iterable)}"
            )
        try:
            lengths.add(len(iterable))
        except (AttributeError, TypeError):
            iterables[i] = list(iterable)
            lengths.add(len(iterables[i]))
    if len(lengths) == 0:
        raise ValueError("`from_map` requires at least one Iterable input")
    elif len(lengths) > 1:
        raise ValueError("All `iterables` must have the same length")
    if lengths == {0}:
        raise ValueError("All `iterables` must have a non-zero length")

    # Check if `func` supports column projection
    allow_projection = False
    columns_arg_required = False
    if param := inspect.signature(func).parameters.get("columns", None):
        allow_projection = True
        columns_arg_required = param.default is param.empty
        if meta is no_default and columns_arg_required:
            raise TypeError(
                "Argument `func` of `from_map` has a required `columns` "
                " parameter and not `meta` provided."
                "Either provide `meta` yourself or make `columns` an optional argument."
            )
    elif isinstance(func, DataFrameIOFunction):
        warnings.warn(
            "dask_expr does not support the DataFrameIOFunction "
            "protocol for column projection. To enable column "
            "projection, please ensure that the signature of `func` "
            "includes a `columns=` keyword argument instead."
        )
    else:
        allow_projection = False

    args = [] if args is None else args
    kwargs = {} if kwargs is None else kwargs
    if allow_projection:
        columns = kwargs.pop("columns", None)
        result = new_collection(
            FromMapProjectable(
                func,
                iterables,
                columns,
                args,
                kwargs,
                columns_arg_required,
                meta,
                enforce_metadata,
                divisions,
                label,
            )
        )
    else:
        result = new_collection(
            FromMap(
                func,
                iterables,
                args,
                kwargs,
                meta,
                enforce_metadata,
                divisions,
                label,
            )
        )
    if pyarrow_strings_enabled():
        return new_collection(expr.ArrowStringConversion(result))
    return result


def repartition(df, divisions, force=False):
    """Repartition dataframe along new divisions

    Dask.DataFrame objects are partitioned along their index.  Often when
    multiple dataframes interact we need to align these partitionings.  The
    ``repartition`` function constructs a new DataFrame object holding the same
    data but partitioned on different values.  It does this by performing a
    sequence of ``loc`` and ``concat`` calls to split and merge the previous
    generation of partitions.

    Parameters
    ----------

    divisions : list
        List of partitions to be used
    force : bool, default False
        Allows the expansion of the existing divisions.
        If False then the new divisions lower and upper bounds must be
        the same as the old divisions.

    Examples
    --------

    >>> df = df.repartition([0, 5, 10, 20])  # doctest: +SKIP

    Also works on Pandas objects

    >>> ddf = dd.repartition(df, [0, 5, 10, 20])  # doctest: +SKIP
    """
    if isinstance(df, FrameBase):
        return df.repartition(divisions=divisions, force=force)
    elif is_dataframe_like(df) or is_series_like(df):
        return new_collection(
            FromPandasDivisions(
                _BackendData(df),
                divisions=divisions,
                pyarrow_strings_enabled=pyarrow_strings_enabled(),
            )
        )
    else:
        raise NotImplementedError(f"repartition is not implemented for {type(df)}.")


def pivot_table(df, index, columns, values, aggfunc="mean"):
    """
    Create a spreadsheet-style pivot table as a DataFrame. Target ``columns``
    must have category dtype to infer result's ``columns``.
    ``index``, ``columns``, and ``aggfunc`` must be all scalar.
    ``values`` can be scalar or list-like.

    Parameters
    ----------
    df : DataFrame
    index : scalar
        column to be index
    columns : scalar
        column to be columns
    values : scalar or list(scalar)
        column(s) to aggregate
    aggfunc : {'mean', 'sum', 'count', 'first', 'last'}, default 'mean'

    Returns
    -------
    table : DataFrame

    See Also
    --------
    pandas.DataFrame.pivot_table
    """
    if not is_scalar(index) or index not in df._meta.columns:
        raise ValueError("'index' must be the name of an existing column")
    if not is_scalar(columns) or columns not in df._meta.columns:
        raise ValueError("'columns' must be the name of an existing column")
    if not methods.is_categorical_dtype(df._meta[columns]):
        raise ValueError("'columns' must be category dtype")
    if not has_known_categories(df._meta[columns]):
        raise ValueError("'columns' must have known categories")

    if not (
        is_scalar(values)
        and values in df._meta.columns
        or not is_scalar(values)
        and all(is_scalar(x) and x in df._meta.columns for x in values)
    ):
        raise ValueError("'values' must refer to an existing column or columns")

    available_aggfuncs = ["mean", "sum", "count", "first", "last"]

    if not is_scalar(aggfunc) or aggfunc not in available_aggfuncs:
        raise ValueError(
            "aggfunc must be either " + ", ".join(f"'{x}'" for x in available_aggfuncs)
        )

    return new_collection(
        PivotTable(df, index=index, columns=columns, values=values, aggfunc=aggfunc)
    )


@derived_from(pd, ua_args=["downcast"])
def to_numeric(arg, errors="raise", downcast=None, meta=None):
    """
    Return type depends on input. Delayed if scalar, otherwise same as input.
    For errors, only "raise" and "coerce" are allowed.
    """
    if errors not in ("raise", "coerce"):
        raise ValueError("invalid error value specified")

    if pd_is_scalar(arg):
        if meta is not None:
            raise KeyError("``meta`` is not allowed when input is a scalar.")
        return delayed(pd.to_numeric, pure=True)(arg, errors=errors, downcast=downcast)

    if is_arraylike(arg):
        return new_collection(
            ToNumeric(
                from_array(arg).astype(arg.dtype), errors=errors, downcast=downcast
            )
        ).to_dask_array(meta=meta)

    if is_series_like(arg):
        return new_collection(
            ToNumeric(frame=arg, errors=errors, downcast=downcast, meta=meta)
        )

    raise TypeError(
        "arg must be a list, tuple, dask.array.Array, or dask.dataframe.Series"
    )


@wraps(pd.to_datetime)
def to_datetime(arg, meta=None, **kwargs):
    tz_kwarg = {"tz": "utc"} if kwargs.get("utc") else {}

    (arg,) = _maybe_from_pandas([arg])

    if meta is None:
        if isinstance(arg, Index):
            meta = get_meta_library(arg).DatetimeIndex([], **tz_kwarg)
            meta.name = arg.name
            if PANDAS_GE_300:
                meta = meta.as_unit(
                    get_meta_library(arg).to_datetime(meta_nonempty(arg._meta)).unit
                )
        elif not (is_dataframe_like(arg) or is_series_like(arg)):
            raise NotImplementedError(
                "dask.dataframe.to_datetime does not support "
                "non-index-able arguments (like scalars)"
            )
        else:
            meta = meta_series_constructor(arg)([pd.Timestamp("2000", **tz_kwarg)])
            meta.index = meta.index.astype(arg.index.dtype)
            meta.index.name = arg.index.name
    else:
        meta = make_meta(meta)

    kwargs.pop("infer_datetime_format", None)

    return new_collection(ToDatetime(frame=arg, kwargs=kwargs, meta=meta))


@wraps(pd.to_timedelta)
def to_timedelta(arg, unit=None, errors="raise"):
    if not isinstance(arg, Series):
        raise TypeError("arg must be a Series")
    return new_collection(ToTimedelta(frame=arg, unit=unit, errors=errors))


def _from_scalars(scalars, meta, names):
    return new_collection(FromScalars(meta, names, *scalars))


@insert_meta_param_description
def map_partitions(
    func,
    *args,
    meta=no_default,
    enforce_metadata=True,
    transform_divisions=True,
    clear_divisions=False,
    align_dataframes=False,
    parent_meta=None,
    **kwargs,
):
    """Apply Python function on each DataFrame partition.

    Parameters
    ----------
    func : function
        Function applied to each partition.
    args, kwargs :
        Arguments and keywords to pass to the function.  At least one of the
        args should be a Dask.dataframe. Arguments and keywords may contain
        ``Scalar``, ``Delayed`` or regular python objects. DataFrame-like args
        (both dask and pandas) will be repartitioned to align (if necessary)
        before applying the function (see ``align_dataframes`` to control).
    enforce_metadata : bool, default True
        Whether to enforce at runtime that the structure of the DataFrame
        produced by ``func`` actually matches the structure of ``meta``.
        This will rename and reorder columns for each partition,
        and will raise an error if this doesn't work,
        but it won't raise if dtypes don't match.
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
    """
    if align_dataframes:
        # TODO: Handle alignment?
        # Perhaps we only handle the case that all `Expr` operands
        # have the same number of partitions or can be broadcasted
        # within `MapPartitions`. If so, the `map_partitions` API
        # will need to call `Repartition` on operands that are not
        # aligned with `self.expr`.
        raise NotImplementedError()
    new_expr = expr.MapPartitions(
        args[0],
        func,
        meta,
        enforce_metadata,
        transform_divisions,
        clear_divisions,
        align_dataframes,
        parent_meta,
        kwargs.pop("token", None),
        kwargs,
        *args[1:],
    )
    return new_collection(new_expr)


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
    clear_divisions=False,
    align_dataframes=False,
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
    if isinstance(before, str):
        before = pd.to_timedelta(before)
    if isinstance(after, str):
        after = pd.to_timedelta(after)

    if isinstance(before, datetime.timedelta) or isinstance(after, datetime.timedelta):
        if isinstance(df, FrameBase):
            inferred_type = df.index._meta_nonempty.inferred_type
        else:
            inferred_type = df.index.inferred_type

        if not is_datetime64_any_dtype(inferred_type):
            raise TypeError(
                "Must have a `DatetimeIndex` when using string offset "
                "for `before` and `after`"
            )

    elif not (
        isinstance(before, Integral)
        and before >= 0
        and isinstance(after, Integral)
        and after >= 0
    ):
        raise ValueError("before and after must be positive integers")

    df = _maybe_from_pandas([df])[0]
    args = _maybe_from_pandas(args)

    if align_dataframes:
        dfs = [df] + args
        dfs = [df for df in dfs if isinstance(df, FrameBase)]
        if len(dfs) > 1 and not expr.are_co_aligned(*dfs):
            return new_collection(
                expr.MapOverlapAlign(
                    df,
                    func,
                    before,
                    after,
                    meta,
                    enforce_metadata,
                    transform_divisions,
                    clear_divisions,
                    align_dataframes,
                    kwargs.pop("token", None),
                    kwargs,
                    *args,
                )
            )

    new_expr = expr.MapOverlap(
        df,
        func,
        before,
        after,
        meta,
        enforce_metadata,
        transform_divisions,
        clear_divisions,
        align_dataframes,
        kwargs.pop("token", None),
        kwargs,
        *args,
    )
    return new_collection(new_expr)


def isna(arg):
    if isinstance(arg, FrameBase):
        return arg.isna()
    else:
        return from_pandas(arg).isna()


def elemwise(op, *args, meta=no_default, out=None, transform_divisions=True, **kwargs):
    """Elementwise operation for Dask dataframes

    Parameters
    ----------
    op: callable
        Function to apply across input dataframes
    *args: DataFrames, Series, Scalars, Arrays,
        The arguments of the operation
    meta: pd.DataFrame, pd.Series (optional)
        Valid metadata for the operation.  Will evaluate on a small piece of
        data if not provided.
    transform_divisions: boolean
        If the input is a ``dask.dataframe.Index`` we normally will also apply
        the function onto the divisions and apply those transformed divisions
        to the output.  You can pass ``transform_divisions=False`` to override
        this behavior
    out : ``dask.array`` or ``None``
        If out is a dask.DataFrame, dask.Series or dask.Scalar then
        this overwrites the contents of it with the result
    **kwargs: scalars

    Examples
    --------
    >>> elemwise(operator.add, df.x, df.y)  # doctest: +SKIP
    """

    args = _maybe_from_pandas(args)
    dfs = [df for df in args if isinstance(df, FrameBase)]
    if len(dfs) <= 1 or expr.are_co_aligned(*dfs):
        result = new_collection(
            expr.UFuncElemwise(dfs[0], op, meta, transform_divisions, kwargs, *args)
        )
    else:
        result = new_collection(expr.UFuncAlign(dfs[0], op, meta, kwargs, *args))
    return handle_out(out, result)


def handle_out(out, result):
    """Handle out parameters

    If out is a dask.DataFrame, dask.Series or dask.Scalar then
    this overwrites the contents of it with the result. The method
    replaces the expression of the out parameter with the result
    from this operation to perform something akin to an inplace
    modification.
    """
    if isinstance(out, tuple):
        if len(out) == 1:
            out = out[0]
        elif len(out) > 1:
            raise NotImplementedError(
                "The `out` parameter with length > 1 is not supported"
            )
        else:
            out = None

    if out is not None and out.__class__ != result.__class__:
        raise TypeError(
            "Mismatched types between result and out parameter. "
            "out=%s, result=%s" % (str(type(out)), str(type(result)))
        )

    if isinstance(out, DataFrame):
        if len(out.columns) != len(result.columns):
            raise ValueError(
                "Mismatched columns count between result and out parameter. "
                "out=%s, result=%s" % (str(len(out.columns)), str(len(result.columns)))
            )

    if isinstance(out, (Series, DataFrame, Scalar)):
        out._expr = result._expr
    elif out is not None:
        msg = (
            "The out parameter is not fully supported."
            " Received type %s, expected %s "
            % (
                typename(type(out)),
                typename(type(result)),
            )
        )
        raise NotImplementedError(msg)
    else:
        return result


def _compute_partition_stats(
    column: Series, allow_overlap: bool = False
) -> tuple[list, list, list[int]]:
    """For a given column, compute the min, max, and len of each partition.

    And make sure that the partitions are sorted relative to each other.
    NOTE: this does not guarantee that every partition is internally sorted.
    """
    mins = column.map_partitions(M.min, meta=column)
    maxes = column.map_partitions(M.max, meta=column)
    lens = column.map_partitions(len, meta=column)
    mins, maxes, lens = compute(mins, maxes, lens)
    mins = mins.bfill().tolist()
    maxes = maxes.bfill().tolist()
    non_empty_mins = [m for m, length in zip(mins, lens) if length != 0]
    non_empty_maxes = [m for m, length in zip(maxes, lens) if length != 0]
    if (
        sorted(non_empty_mins) != non_empty_mins
        or sorted(non_empty_maxes) != non_empty_maxes
    ):
        raise ValueError(
            f"Partitions are not sorted ascending by {column.name or 'the index'}. ",
            f"In your dataset the (min, max, len) values of {column.name or 'the index'} "
            f"for each partition are: {list(zip(mins, maxes, lens))}",
        )
    if not allow_overlap and any(
        a <= b for a, b in zip(non_empty_mins[1:], non_empty_maxes[:-1])
    ):
        warnings.warn(
            "Partitions have overlapping values, so divisions are non-unique. "
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


@get_parallel_type.register(FrameBase)
def get_parallel_type_frame(o):
    return get_parallel_type(o._meta)
