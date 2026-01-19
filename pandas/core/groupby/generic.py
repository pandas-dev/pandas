"""
Define the SeriesGroupBy and DataFrameGroupBy
classes that hold the groupby interfaces (and some implementations).

These are user facing as the result of the ``df.groupby(...)`` operations,
which here returns a DataFrameGroupBy object.
"""

from __future__ import annotations

from collections import abc
from collections.abc import Callable
import dataclasses
from functools import partial
from textwrap import dedent
from typing import (
    TYPE_CHECKING,
    Any,
    Literal,
    TypeAlias,
    TypeVar,
    cast,
)
import warnings

import numpy as np

from pandas._libs import Interval
from pandas._libs.hashtable import duplicated
from pandas.errors import (
    Pandas4Warning,
    SpecificationError,
)
from pandas.util._decorators import set_module
from pandas.util._exceptions import find_stack_level

from pandas.core.dtypes.common import (
    ensure_int64,
    is_bool,
    is_dict_like,
    is_integer_dtype,
    is_list_like,
    is_numeric_dtype,
    is_scalar,
)
from pandas.core.dtypes.dtypes import (
    CategoricalDtype,
    IntervalDtype,
)
from pandas.core.dtypes.inference import is_hashable
from pandas.core.dtypes.missing import (
    isna,
    notna,
)

from pandas.core import algorithms
from pandas.core.apply import (
    GroupByApply,
    maybe_mangle_lambdas,
    reconstruct_func,
    validate_func_kwargs,
)
import pandas.core.common as com
from pandas.core.frame import DataFrame
from pandas.core.groupby import base
from pandas.core.groupby.groupby import (
    GroupBy,
    GroupByPlot,
)
from pandas.core.indexes.api import (
    Index,
    MultiIndex,
    all_indexes_same,
    default_index,
)
from pandas.core.series import Series
from pandas.core.sorting import get_group_index
from pandas.core.util.numba_ import maybe_use_numba

from pandas.plotting import boxplot_frame_groupby

if TYPE_CHECKING:
    from collections.abc import (
        Hashable,
        Sequence,
    )

    from pandas._typing import (
        ArrayLike,
        BlockManager,
        CorrelationMethod,
        IndexLabel,
        Manager,
        SingleBlockManager,
        TakeIndexer,
    )

    from pandas import Categorical
    from pandas.core.generic import NDFrame

# TODO(typing) the return value on this callable should be any *scalar*.
AggScalar: TypeAlias = str | Callable[..., Any]
# TODO: validate types on ScalarResult and move to _typing
# Blocked from using by https://github.com/python/mypy/issues/1484
# See note at _mangle_lambda_list
ScalarResult = TypeVar("ScalarResult")


@set_module("pandas")
@dataclasses.dataclass
class NamedAgg:
    """
    Helper for column specific aggregation with control over output column names.

    Parameters
    ----------
    column : Hashable
        Column label in the DataFrame to apply aggfunc.
    aggfunc : function or str
        Function to apply to the provided column. If string, the name of a built-in
        pandas function.
    *args, **kwargs : Any
        Optional positional and keyword arguments passed to ``aggfunc``.

    See Also
    --------
    DataFrame.groupby : Group DataFrame using a mapper or by a Series of columns.

    Examples
    --------
    >>> df = pd.DataFrame({"key": [1, 1, 2], "a": [-1, 0, 1], 1: [10, 11, 12]})
    >>> agg_a = pd.NamedAgg(column="a", aggfunc="min")
    >>> agg_1 = pd.NamedAgg(column=1, aggfunc=lambda x: np.mean(x))
    >>> df.groupby("key").agg(result_a=agg_a, result_1=agg_1)
        result_a  result_1
    key
    1          -1      10.5
    2           1      12.0

    >>> def n_between(ser, low, high, **kwargs):
    ...     return ser.between(low, high, **kwargs).sum()

    >>> agg_between = pd.NamedAgg("a", n_between, 0, 1)
    >>> df.groupby("key").agg(count_between=agg_between)
        count_between
    key
    1               1
    2               1

    >>> agg_between_kw = pd.NamedAgg("a", n_between, 0, 1, inclusive="both")
    >>> df.groupby("key").agg(count_between_kw=agg_between_kw)
        count_between_kw
    key
    1                   1
    2                   1
    """

    column: Hashable
    aggfunc: AggScalar
    args: tuple[Any, ...] = ()
    kwargs: dict[str, Any] = dataclasses.field(default_factory=dict)

    def __init__(
        self,
        column: Hashable,
        aggfunc: Callable[..., Any] | str,
        *args: Any,
        **kwargs: Any,
    ) -> None:
        self.column = column
        self.aggfunc = aggfunc
        self.args = args
        self.kwargs = kwargs

    def __getitem__(self, key: int) -> Any:
        """Provide backward-compatible tuple-style access."""
        if key == 0:
            return self.column
        elif key == 1:
            return self.aggfunc
        elif key == 2:
            return self.args
        elif key == 3:
            return self.kwargs
        raise IndexError("index out of range")


@set_module("pandas.api.typing")
class SeriesGroupBy(GroupBy[Series]):
    def _wrap_agged_manager(self, mgr: Manager) -> Series:
        out = self.obj._constructor_from_mgr(mgr, axes=mgr.axes)
        out._name = self.obj.name
        return out

    def _get_data_to_aggregate(
        self, *, numeric_only: bool = False, name: str | None = None
    ) -> SingleBlockManager:
        ser = self._obj_with_exclusions
        single = ser._mgr
        if numeric_only and not is_numeric_dtype(ser.dtype):
            # GH#41291 match Series behavior
            kwd_name = "numeric_only"
            raise TypeError(
                f"Cannot use {kwd_name}=True with "
                f"{type(self).__name__}.{name} and non-numeric dtypes."
            )
        return single

    def apply(self, func, *args, **kwargs) -> Series:
        """
        Apply function ``func`` group-wise and combine the results together.

        The function passed to ``apply`` must take a series as its first
        argument and return a DataFrame, Series or scalar. ``apply`` will
        then take care of combining the results back together into a single
        dataframe or series. ``apply`` is therefore a highly flexible
        grouping method.

        While ``apply`` is a very flexible method, its downside is that
        using it can be quite a bit slower than using more specific methods
        like ``agg`` or ``transform``. Pandas offers a wide range of method that will
        be much faster than using ``apply`` for their specific purposes, so try to
        use them before reaching for ``apply``.

        Parameters
        ----------
        func : callable
            A callable that takes a series as its first argument, and
            returns a dataframe, a series or a scalar. In addition the
            callable may take positional and keyword arguments.

        *args : tuple
            Optional positional arguments to pass to ``func``.

        **kwargs : dict
            Optional keyword arguments to pass to ``func``.

        Returns
        -------
        Series or DataFrame
            A pandas object with the result of applying ``func`` to each group.

        See Also
        --------
        pipe : Apply function to the full GroupBy object instead of to each
            group.
        aggregate : Apply aggregate function to the GroupBy object.
        transform : Apply function column-by-column to the GroupBy object.
        Series.apply : Apply a function to a Series.
        DataFrame.apply : Apply a function to each row or column of a DataFrame.

        Notes
        -----
        The resulting dtype will reflect the return value of the passed ``func``,
        see the examples below.

        Functions that mutate the passed object can produce unexpected
        behavior or errors and are not supported. See :ref:`gotchas.udf-mutation`
        for more details.

        Examples
        --------
        >>> s = pd.Series([0, 1, 2], index="a a b".split())
        >>> g1 = s.groupby(s.index, group_keys=False)
        >>> g2 = s.groupby(s.index, group_keys=True)

        From ``s`` above we can see that ``g`` has two groups, ``a`` and ``b``.
        Notice that ``g1`` have ``g2`` have two groups, ``a`` and ``b``, and only
        differ in their ``group_keys`` argument. Calling `apply` in various ways,
        we can get different grouping results:

        Example 1: The function passed to `apply` takes a Series as
        its argument and returns a Series.  `apply` combines the result for
        each group together into a new Series.

        The resulting dtype will reflect the return value of the passed ``func``.

        >>> g1.apply(lambda x: x * 2 if x.name == "a" else x / 2)
        a    0.0
        a    2.0
        b    1.0
        dtype: float64

        In the above, the groups are not part of the index. We can have them included
        by using ``g2`` where ``group_keys=True``:

        >>> g2.apply(lambda x: x * 2 if x.name == "a" else x / 2)
        a  a    0.0
           a    2.0
        b  b    1.0
        dtype: float64

        Example 2: The function passed to `apply` takes a Series as
        its argument and returns a scalar. `apply` combines the result for
        each group together into a Series, including setting the index as
        appropriate:

        >>> g1.apply(lambda x: x.max() - x.min())
        a    1
        b    0
        dtype: int64

        The ``group_keys`` argument has no effect here because the result is not
        like-indexed (i.e. :ref:`a transform <groupby.transform>`) when compared
        to the input.

        >>> g2.apply(lambda x: x.max() - x.min())
        a    1
        b    0
        dtype: int64
        """
        return super().apply(func, *args, **kwargs)

    def aggregate(self, func=None, *args, engine=None, engine_kwargs=None, **kwargs):
        """
        Aggregate using one or more operations.

        The ``aggregate`` method enables flexible and efficient aggregation of grouped
        data using a variety of functions, including built-in, user-defined, and
        optimized JIT-compiled functions.

        Parameters
        ----------
        func : function, str, list or None
            Function to use for aggregating the data. If a function, must either
            work when passed a Series or when passed to Series.apply.

            Accepted combinations are:

            - function
            - string function name
            - list of functions and/or function names, e.g. ``[np.sum, 'mean']``
            - None, in which case ``**kwargs`` are used with Named Aggregation. Here
              the output has one column for each element in ``**kwargs``. The name of
              the column is keyword, whereas the value determines the aggregation
              used to compute the values in the column.

              Can also accept a Numba JIT function with
              ``engine='numba'`` specified. Only passing a single function is supported
              with this engine.

              If the ``'numba'`` engine is chosen, the function must be
              a user defined function with ``values`` and ``index`` as the
              first and second arguments respectively in the function signature.
              Each group's index will be passed to the user defined function
              and optionally available for use.

        *args
            Positional arguments to pass to func.
        engine : str, default None
            * ``'cython'`` : Runs the function through C-extensions from cython.
            * ``'numba'`` : Runs the function through JIT compiled code from numba.
            * ``None`` : Defaults to ``'cython'`` or globally setting
                ``compute.use_numba``

        engine_kwargs : dict, default None
            * For ``'cython'`` engine, there are no accepted ``engine_kwargs``
            * For ``'numba'`` engine, the engine can accept ``nopython``, ``nogil``
              and ``parallel`` dictionary keys. The values must either be ``True`` or
              ``False``. The default ``engine_kwargs`` for the ``'numba'`` engine is
              ``{'nopython': True, 'nogil': False, 'parallel': False}`` and will be
              applied to the function

        **kwargs
            * If ``func`` is None, ``**kwargs`` are used to define the output names and
              aggregations via Named Aggregation. See ``func`` entry.
            * Otherwise, keyword arguments to be passed into func.

        Returns
        -------
        Series
            Aggregated Series based on the grouping and the applied aggregation
            functions.

        See Also
        --------
        SeriesGroupBy.apply : Apply function func group-wise
            and combine the results together.
        SeriesGroupBy.transform : Transforms the Series on each group
            based on the given function.
        Series.aggregate : Aggregate using one or more operations.

        Notes
        -----
        When using ``engine='numba'``, there will be no "fall back" behavior internally.
        The group data and group index will be passed as numpy arrays to the JITed
        user defined function, and no alternative execution attempts will be tried.

        Functions that mutate the passed object can produce unexpected
        behavior or errors and are not supported. See :ref:`gotchas.udf-mutation`
        for more details.

        The resulting dtype will reflect the return value of the passed ``func``,
        see the examples below.

        Examples
        --------
        >>> s = pd.Series([1, 2, 3, 4])

        >>> s
        0    1
        1    2
        2    3
        3    4
        dtype: int64

        >>> s.groupby([1, 1, 2, 2]).min()
        1    1
        2    3
        dtype: int64

        >>> s.groupby([1, 1, 2, 2]).agg("min")
        1    1
        2    3
        dtype: int64

        >>> s.groupby([1, 1, 2, 2]).agg(["min", "max"])
           min  max
        1    1    2
        2    3    4

        The output column names can be controlled by passing
        the desired column names and aggregations as keyword arguments.

        >>> s.groupby([1, 1, 2, 2]).agg(
        ...     minimum="min",
        ...     maximum="max",
        ... )
           minimum  maximum
        1        1        2
        2        3        4

        The resulting dtype will reflect the return value of the aggregating
        function.

        >>> s.groupby([1, 1, 2, 2]).agg(lambda x: x.astype(float).min())
        1    1.0
        2    3.0
        dtype: float64
        """
        relabeling = func is None
        columns = None
        if relabeling:
            columns, func = validate_func_kwargs(kwargs)
            kwargs = {}

        if isinstance(func, str):
            if maybe_use_numba(engine) and engine is not None:
                # Not all agg functions support numba, only propagate numba kwargs
                # if user asks for numba, and engine is not None
                # (if engine is None, the called function will handle the case where
                # numba is requested via the global option)
                kwargs["engine"] = engine
            if engine_kwargs is not None:
                kwargs["engine_kwargs"] = engine_kwargs
            return getattr(self, func)(*args, **kwargs)

        elif isinstance(func, abc.Iterable):
            # Catch instances of lists / tuples
            # but not the class list / tuple itself.
            func = maybe_mangle_lambdas(func)
            kwargs["engine"] = engine
            kwargs["engine_kwargs"] = engine_kwargs
            ret = self._aggregate_multiple_funcs(func, *args, **kwargs)
            if relabeling:
                # columns is not narrowed by mypy from relabeling flag
                assert columns is not None  # for mypy
                ret.columns = columns
            if not self.as_index:
                ret = ret.reset_index()
            return ret

        else:
            if maybe_use_numba(engine):
                return self._aggregate_with_numba(
                    func, *args, engine_kwargs=engine_kwargs, **kwargs
                )

            if self.ngroups == 0:
                # e.g. test_evaluate_with_empty_groups without any groups to
                #  iterate over, we have no output on which to do dtype
                #  inference. We default to using the existing dtype.
                #  xref GH#51445
                obj = self._obj_with_exclusions
                return self._wrap_aggregated_output(
                    self.obj._constructor(
                        [],
                        name=self.obj.name,
                        index=self._grouper.result_index,
                        dtype=obj.dtype,
                    )
                )
            return self._python_agg_general(func, *args, **kwargs)

    agg = aggregate

    def _python_agg_general(self, func, *args, **kwargs):
        f = lambda x: func(x, *args, **kwargs)

        obj = self._obj_with_exclusions
        result = self._grouper.agg_series(obj, f)
        res = obj._constructor(result, name=obj.name)
        return self._wrap_aggregated_output(res)

    def _aggregate_multiple_funcs(self, arg, *args, **kwargs) -> DataFrame:
        if isinstance(arg, dict):
            raise SpecificationError("nested renamer is not supported")

        if any(isinstance(x, (tuple, list)) for x in arg):
            arg = ((x, x) if not isinstance(x, (tuple, list)) else x for x in arg)
        else:
            # list of functions / function names
            columns = (com.get_callable_name(f) or f for f in arg)
            arg = zip(columns, arg, strict=True)

        results: dict[base.OutputKey, DataFrame | Series] = {}
        with com.temp_setattr(self, "as_index", True):
            # Combine results using the index, need to adjust index after
            # if as_index=False (GH#50724)
            for idx, (name, func) in enumerate(arg):
                key = base.OutputKey(label=name, position=idx)
                results[key] = self.aggregate(func, *args, **kwargs)

        if any(isinstance(x, DataFrame) for x in results.values()):
            from pandas import concat

            res_df = concat(
                results.values(), axis=1, keys=[key.label for key in results]
            )
            return res_df

        indexed_output = {key.position: val for key, val in results.items()}
        output = self.obj._constructor_expanddim(indexed_output, index=None)
        output.columns = Index(key.label for key in results)

        return output

    def _wrap_applied_output(
        self,
        data: Series,
        values: list[Any],
        not_indexed_same: bool = False,
        is_transform: bool = False,
    ) -> DataFrame | Series:
        """
        Wrap the output of SeriesGroupBy.apply into the expected result.

        Parameters
        ----------
        data : Series
            Input data for groupby operation.
        values : List[Any]
            Applied output for each group.
        not_indexed_same : bool, default False
            Whether the applied outputs are not indexed the same as the group axes.

        Returns
        -------
        DataFrame or Series
        """
        if len(values) == 0:
            # GH #6265
            if is_transform:
                # GH#47787 see test_group_on_empty_multiindex
                res_index = data.index
            elif not self.group_keys:
                res_index = None
            else:
                res_index = self._grouper.result_index

            return self.obj._constructor(
                [],
                name=self.obj.name,
                index=res_index,
                dtype=data.dtype,
            )
        assert values is not None

        if isinstance(values[0], dict):
            # GH #823 #24880
            index = self._grouper.result_index
            res_df = self.obj._constructor_expanddim(values, index=index)
            # if self.observed is False,
            # keep all-NaN rows created while re-indexing
            res_ser = res_df.stack()
            res_ser.name = self.obj.name
            return res_ser
        elif isinstance(values[0], (Series, DataFrame)):
            result = self._concat_objects(
                values,
                not_indexed_same=not_indexed_same,
                is_transform=is_transform,
            )
            if isinstance(result, Series):
                result.name = self.obj.name
            if not self.as_index and not_indexed_same:
                result = self._insert_inaxis_grouper(result)
                result.index = default_index(len(result))
            return result.__finalize__(self.obj, method="groupby")
        else:
            # GH #6265 #24880
            result = self.obj._constructor(
                data=values, index=self._grouper.result_index, name=self.obj.name
            )
            if not self.as_index:
                result = self._insert_inaxis_grouper(result)
                result.index = default_index(len(result))
            return result.__finalize__(self.obj, method="groupby")

    __examples_series_doc = dedent(
        """
    >>> ser = pd.Series([390.0, 350.0, 30.0, 20.0],
    ...                 index=["Falcon", "Falcon", "Parrot", "Parrot"],
    ...                 name="Max Speed")
    >>> grouped = ser.groupby([1, 1, 2, 2])
    >>> grouped.transform(lambda x: (x - x.mean()) / x.std())
        Falcon    0.707107
        Falcon   -0.707107
        Parrot    0.707107
        Parrot   -0.707107
        Name: Max Speed, dtype: float64

    Broadcast result of the transformation

    >>> grouped.transform(lambda x: x.max() - x.min())
    Falcon    40.0
    Falcon    40.0
    Parrot    10.0
    Parrot    10.0
    Name: Max Speed, dtype: float64

    >>> grouped.transform("mean")
    Falcon    370.0
    Falcon    370.0
    Parrot     25.0
    Parrot     25.0
    Name: Max Speed, dtype: float64

    The resulting dtype will reflect the return value of the passed ``func``,
    for example:

    >>> grouped.transform(lambda x: x.astype(int).max())
    Falcon    390
    Falcon    390
    Parrot     30
    Parrot     30
    Name: Max Speed, dtype: int64
    """
    )

    def transform(self, func, *args, engine=None, engine_kwargs=None, **kwargs):
        """
        Call function producing a same-indexed Series on each group.

        Returns a Series having the same indexes as the original object
        filled with the transformed values.

        Parameters
        ----------
        func : function, str
            Function to apply to each group.
            See the Notes section below for requirements.

            Accepted inputs are:

            - String
            - Python function
            - Numba JIT function with ``engine='numba'`` specified.

            Only passing a single function is supported with this engine.
            If the ``'numba'`` engine is chosen, the function must be
            a user defined function with ``values`` and ``index`` as the
            first and second arguments respectively in the function signature.
            Each group's index will be passed to the user defined function
            and optionally available for use.

            If a string is chosen, then it needs to be the name
            of the groupby method you want to use.
        *args
            Positional arguments to pass to func.
        engine : str, default None
            * ``'cython'`` : Runs the function through C-extensions from cython.
            * ``'numba'`` : Runs the function through JIT compiled code from numba.
            * ``None`` : Defaults to ``'cython'``
              or the global setting ``compute.use_numba``

        engine_kwargs : dict, default None
            * For ``'cython'`` engine, there are no accepted ``engine_kwargs``
            * For ``'numba'`` engine, the engine can accept ``nopython``, ``nogil``
              and ``parallel`` dictionary keys. The values must either be ``True`` or
              ``False``. The default ``engine_kwargs`` for the ``'numba'`` engine is
              ``{'nopython': True, 'nogil': False, 'parallel': False}`` and will be
              applied to the function

        **kwargs
            Keyword arguments to be passed into func.

        Returns
        -------
        Series
            Series with the same indexes as the original object filled
            with transformed values.

        See Also
        --------
        Series.groupby.apply : Apply function ``func`` group-wise and combine
            the results together.
        Series.groupby.aggregate : Aggregate using one or more operations.
        Series.transform : Call ``func`` on self producing a Series with the
            same axis shape as self.

        Notes
        -----
        Each group is endowed the attribute 'name' in case you need to know
        which group you are working on.

        The current implementation imposes three requirements on f:

        * f must return a value that either has the same shape as the input
          subframe or can be broadcast to the shape of the input subframe.
          For example, if `f` returns a scalar it will be broadcast to have the
          same shape as the input subframe.
        * if this is a DataFrame, f must support application column-by-column
          in the subframe. If f also supports application to the entire subframe,
          then a fast path is used starting from the second chunk.
        * f must not mutate groups. Mutation is not supported and may
          produce unexpected results. See :ref:`gotchas.udf-mutation` for more details.

        When using ``engine='numba'``, there will be no "fall back" behavior internally.
        The group data and group index will be passed as numpy arrays to the JITed
        user defined function, and no alternative execution attempts will be tried.

        The resulting dtype will reflect the return value of the passed ``func``,
        see the examples below.

        .. versionchanged:: 2.0.0

            When using ``.transform`` on a grouped DataFrame and
            the transformation function returns a DataFrame,
            pandas now aligns the result's index with the input's index.
            You can call ``.to_numpy()`` on the result of
            the transformation function to avoid alignment.

        Examples
        --------

        >>> ser = pd.Series(
        ...     [390.0, 350.0, 30.0, 20.0],
        ...     index=["Falcon", "Falcon", "Parrot", "Parrot"],
        ...     name="Max Speed",
        ... )
        >>> grouped = ser.groupby([1, 1, 2, 2])
        >>> grouped.transform(lambda x: (x - x.mean()) / x.std())
            Falcon    0.707107
            Falcon   -0.707107
            Parrot    0.707107
            Parrot   -0.707107
            Name: Max Speed, dtype: float64

        Broadcast result of the transformation

        >>> grouped.transform(lambda x: x.max() - x.min())
        Falcon    40.0
        Falcon    40.0
        Parrot    10.0
        Parrot    10.0
        Name: Max Speed, dtype: float64

        >>> grouped.transform("mean")
        Falcon    370.0
        Falcon    370.0
        Parrot     25.0
        Parrot     25.0
        Name: Max Speed, dtype: float64

        The resulting dtype will reflect the return value of the passed ``func``,
        for example:

        >>> grouped.transform(lambda x: x.astype(int).max())
        Falcon    390
        Falcon    390
        Parrot     30
        Parrot     30
        Name: Max Speed, dtype: int64
        """
        return self._transform(
            func, *args, engine=engine, engine_kwargs=engine_kwargs, **kwargs
        )

    def _cython_transform(self, how: str, numeric_only: bool = False, **kwargs):
        obj = self._obj_with_exclusions

        try:
            result = self._grouper._cython_operation(
                "transform", obj._values, how, 0, **kwargs
            )
        except NotImplementedError as err:
            # e.g. test_groupby_raises_string
            raise TypeError(f"{how} is not supported for {obj.dtype} dtype") from err

        return obj._constructor(result, index=self.obj.index, name=obj.name)

    def _transform_general(
        self, func: Callable, engine, engine_kwargs, *args, **kwargs
    ) -> Series:
        """
        Transform with a callable `func`.
        """
        if maybe_use_numba(engine):
            return self._transform_with_numba(
                func, *args, engine_kwargs=engine_kwargs, **kwargs
            )
        assert callable(func)
        klass = type(self.obj)

        results = []
        for name, group in self._grouper.get_iterator(
            self._obj_with_exclusions,
        ):
            # this setattr is needed for test_transform_lambda_with_datetimetz
            object.__setattr__(group, "name", name)
            res = func(group, *args, **kwargs)

            results.append(klass(res, index=group.index))

        # check for empty "results" to avoid concat ValueError
        if results:
            from pandas.core.reshape.concat import concat

            concatenated = concat(results, ignore_index=True)
            result = self._set_result_index_ordered(concatenated)
        else:
            result = self.obj._constructor(dtype=np.float64)

        result.name = self.obj.name
        return result

    def filter(self, func, dropna: bool = True, *args, **kwargs):
        """
        Filter elements from groups that don't satisfy a criterion.

        Elements from groups are filtered if they do not satisfy the
        boolean criterion specified by func.

        Parameters
        ----------
        func : function
            Criterion to apply to each group. Should return True or False.
        dropna : bool, optional
            Drop groups that do not pass the filter. True by default; if False,
            groups that evaluate False are filled with NaNs.
        *args : tuple
            Optional positional arguments to pass to `func`.
        **kwargs : dict
            Optional keyword arguments to pass to `func`.

        Returns
        -------
        Series
            The filtered subset of the original Series.

        See Also
        --------
        Series.filter: Filter elements of ungrouped Series.
        DataFrameGroupBy.filter : Filter elements from groups base on criterion.

        Notes
        -----
        Functions that mutate the passed object can produce unexpected
        behavior or errors and are not supported. See :ref:`gotchas.udf-mutation`
        for more details.

        Examples
        --------
        >>> df = pd.DataFrame(
        ...     {
        ...         "A": ["foo", "bar", "foo", "bar", "foo", "bar"],
        ...         "B": [1, 2, 3, 4, 5, 6],
        ...         "C": [2.0, 5.0, 8.0, 1.0, 2.0, 9.0],
        ...     }
        ... )
        >>> grouped = df.groupby("A")
        >>> df.groupby("A").B.filter(lambda x: x.mean() > 3.0)
        1    2
        3    4
        5    6
        Name: B, dtype: int64
        """
        if isinstance(func, str):
            wrapper = lambda x: getattr(x, func)(*args, **kwargs)
        else:
            wrapper = lambda x: func(x, *args, **kwargs)

        # Interpret np.nan as False.
        def true_and_notna(x) -> bool:
            b = wrapper(x)
            return notna(b) and b

        try:
            indices = [
                self._get_index(name)
                for name, group in self._grouper.get_iterator(self._obj_with_exclusions)
                if true_and_notna(group)
            ]
        except (ValueError, TypeError) as err:
            raise TypeError("the filter must return a boolean result") from err

        filtered = self._apply_filter(indices, dropna)
        return filtered

    def nunique(self, dropna: bool = True) -> Series | DataFrame:
        """
        Return number of unique elements in the group.

        Parameters
        ----------
        dropna : bool, default True
            Don't include NaN in the counts.

        Returns
        -------
        Series
            Number of unique values within each group.

        See Also
        --------
        core.resample.Resampler.nunique : Method nunique for Resampler.

        Examples
        --------
        >>> lst = ["a", "a", "b", "b"]
        >>> ser = pd.Series([1, 2, 3, 3], index=lst)
        >>> ser
        a    1
        a    2
        b    3
        b    3
        dtype: int64
        >>> ser.groupby(level=0).nunique()
        a    2
        b    1
        dtype: int64
        """
        ids = self._grouper.ids
        ngroups = self._grouper.ngroups
        val = self.obj._values
        codes, uniques = algorithms.factorize(val, use_na_sentinel=dropna, sort=False)

        if self._grouper.has_dropped_na:
            mask = ids >= 0
            ids = ids[mask]
            codes = codes[mask]

        group_index = get_group_index(
            labels=[ids, codes],
            shape=(ngroups, len(uniques)),
            sort=False,
            xnull=dropna,
        )

        if dropna:
            mask = group_index >= 0
            if (~mask).any():
                ids = ids[mask]
                group_index = group_index[mask]

        mask = duplicated(group_index, "first")
        res = np.bincount(ids[~mask], minlength=ngroups)
        res = ensure_int64(res)

        ri = self._grouper.result_index
        result: Series | DataFrame = self.obj._constructor(
            res, index=ri, name=self.obj.name
        )
        if not self.as_index:
            result = self._insert_inaxis_grouper(result)
            result.index = default_index(len(result))
        return result

    def describe(self, percentiles=None, include=None, exclude=None) -> Series:
        """
        Generate descriptive statistics.

        Descriptive statistics include those that summarize the central
        tendency, dispersion and shape of a
        dataset's distribution, excluding ``NaN`` values.

        Analyzes both numeric and object series, as well
        as ``DataFrame`` column sets of mixed data types. The output
        will vary depending on what is provided. Refer to the notes
        below for more detail.

        Parameters
        ----------
        percentiles : list-like of numbers, optional
            The percentiles to include in the output. All should
            fall between 0 and 1. The default, ``None``, will automatically
            return the 25th, 50th, and 75th percentiles.
        include : 'all', list-like of dtypes or None (default), optional
            A white list of data types to include in the result. Ignored
            for ``Series``. Here are the options:

            - 'all' : All columns of the input will be included in the output.
            - A list-like of dtypes : Limits the results to the
              provided data types.
              To limit the result to numeric types submit
              ``numpy.number``. To limit it instead to object columns submit
              the ``numpy.object`` data type. Strings
              can also be used in the style of
              ``select_dtypes`` (e.g. ``df.describe(include=['O'])``). To
              select pandas categorical columns, use ``'category'``
            - None (default) : The result will include all numeric columns.
        exclude : list-like of dtypes or None (default), optional,
            A black list of data types to omit from the result. Ignored
            for ``Series``. Here are the options:

            - A list-like of dtypes : Excludes the provided data types
              from the result. To exclude numeric types submit
              ``numpy.number``. To exclude object columns submit the data
              type ``numpy.object``. Strings can also be used in the style of
              ``select_dtypes`` (e.g. ``df.describe(exclude=['O'])``). To
              exclude pandas categorical columns, use ``'category'``
            - None (default) : The result will exclude nothing.

        Returns
        -------
        Series or DataFrame
            Summary statistics of the Series or Dataframe provided.

        See Also
        --------
        DataFrame.count: Count number of non-NA/null observations.
        DataFrame.max: Maximum of the values in the object.
        DataFrame.min: Minimum of the values in the object.
        DataFrame.mean: Mean of the values.
        DataFrame.std: Standard deviation of the observations.
        DataFrame.select_dtypes: Subset of a DataFrame including/excluding
            columns based on their dtype.

        Notes
        -----
        For numeric data, the result's index will include ``count``,
        ``mean``, ``std``, ``min``, ``max`` as well as lower, ``50`` and
        upper percentiles. By default the lower percentile is ``25`` and the
        upper percentile is ``75``. The ``50`` percentile is the
        same as the median.

        For object data (e.g. strings), the result's index
        will include ``count``, ``unique``, ``top``, and ``freq``. The ``top``
        is the most common value. The ``freq`` is the most common value's
        frequency.

        If multiple object values have the highest count, then the
        ``count`` and ``top`` results will be arbitrarily chosen from
        among those with the highest count.

        For mixed data types provided via a ``DataFrame``, the default is to
        return only an analysis of numeric columns. If the DataFrame consists
        only of object and categorical data without any numeric columns, the
        default is to return an analysis of both the object and categorical
        columns. If ``include='all'`` is provided as an option, the result
        will include a union of attributes of each type.

        The `include` and `exclude` parameters can be used to limit
        which columns in a ``DataFrame`` are analyzed for the output.
        The parameters are ignored when analyzing a ``Series``.

        Examples
        --------
        Describing a numeric ``Series``.

        >>> s = pd.Series([1, 2, 3, 4])

        >>> s
        0    1
        1    2
        2    3
        3    4
        dtype: int64

        >>> s.groupby([1, 1, 2, 2]).describe()
           count  mean       std  min   25%  50%   75%  max
        1    2.0   1.5  0.707107  1.0  1.25  1.5  1.75  2.0
        2    2.0   3.5  0.707107  3.0  3.25  3.5  3.75  4.0
        """
        return super().describe(
            percentiles=percentiles, include=include, exclude=exclude
        )

    def value_counts(
        self,
        normalize: bool = False,
        sort: bool = True,
        ascending: bool = False,
        bins=None,
        dropna: bool = True,
    ) -> Series | DataFrame:
        """
        Return a Series or DataFrame containing counts of unique rows.

        Parameters
        ----------
        normalize : bool, default False
            Return proportions rather than frequencies.
        sort : bool, default True
            Sort by frequencies.
        ascending : bool, default False
            Sort in ascending order.
        bins : int or list of ints, optional
            Rather than count values, group them into half-open bins,
            a convenience for pd.cut, only works with numeric data.
        dropna : bool, default True
            Don't include counts of rows that contain NA values.

        Returns
        -------
        Series or DataFrame
            Series if the groupby ``as_index`` is True, otherwise DataFrame.

        See Also
        --------
        Series.value_counts: Equivalent method on Series.
        DataFrame.value_counts: Equivalent method on DataFrame.
        DataFrameGroupBy.value_counts: Equivalent method on DataFrameGroupBy.

        Notes
        -----
        - If the groupby ``as_index`` is True then the returned Series will have a
          MultiIndex with one level per input column.
        - If the groupby ``as_index`` is False then the returned DataFrame will have an
          additional column with the value_counts. The column is labelled 'count' or
          'proportion', depending on the ``normalize`` parameter.

        By default, rows that contain any NA values are omitted from
        the result.

        By default, the result will be in descending order so that the
        first element of each group is the most frequently-occurring row.

        Examples
        --------
        >>> s = pd.Series(
        ...     [1, 1, 2, 3, 2, 3, 3, 1, 1, 3, 3, 3],
        ...     index=["A", "A", "A", "A", "A", "A", "B", "B", "B", "B", "B", "B"],
        ... )
        >>> s
        A    1
        A    1
        A    2
        A    3
        A    2
        A    3
        B    3
        B    1
        B    1
        B    3
        B    3
        B    3
        dtype: int64
        >>> g1 = s.groupby(s.index)
        >>> g1.value_counts(bins=2)
        A  (0.997, 2.0]    4
           (2.0, 3.0]      2
        B  (2.0, 3.0]      4
           (0.997, 2.0]    2
        Name: count, dtype: int64
        >>> g1.value_counts(normalize=True)
        A  1    0.333333
           2    0.333333
           3    0.333333
        B  3    0.666667
           1    0.333333
        Name: proportion, dtype: float64
        """
        name = "proportion" if normalize else "count"

        if bins is None:
            result = self._value_counts(
                normalize=normalize, sort=sort, ascending=ascending, dropna=dropna
            )
            result.name = name
            return result

        from pandas.core.reshape.merge import get_join_indexers
        from pandas.core.reshape.tile import cut

        ids = self._grouper.ids
        val = self.obj._values

        index_names = [*self._grouper.names, self.obj.name]

        if isinstance(val.dtype, CategoricalDtype) or (
            bins is not None and not np.iterable(bins)
        ):
            # scalar bins cannot be done at top level
            # in a backward compatible way
            # GH38672 relates to categorical dtype
            ser = self.apply(
                Series.value_counts,
                normalize=normalize,
                sort=sort,
                ascending=ascending,
                bins=bins,
            )
            ser.name = name
            ser.index.names = index_names
            return ser

        # groupby removes null keys from groupings
        mask = ids != -1
        ids, val = ids[mask], val[mask]

        lab: Index | np.ndarray
        if bins is None:
            lab, lev = algorithms.factorize(val, sort=True)
            llab = lambda lab, inc: lab[inc]
        else:
            # lab is a Categorical with categories an IntervalIndex
            cat_ser = cut(Series(val, copy=False), bins, include_lowest=True)
            cat_obj = cast("Categorical", cat_ser._values)
            lev = cat_obj.categories
            lab = lev.take(
                cat_obj.codes,
                allow_fill=True,
                fill_value=lev._na_value,
            )
            llab = lambda lab, inc: lab[inc]._multiindex.codes[-1]

        if isinstance(lab.dtype, IntervalDtype):
            # TODO: should we do this inside II?
            lab_interval = cast(Interval, lab)

            sorter = np.lexsort((lab_interval.left, lab_interval.right, ids))
        else:
            sorter = np.lexsort((lab, ids))

        ids, lab = ids[sorter], lab[sorter]

        # group boundaries are where group ids change
        idchanges = 1 + np.nonzero(ids[1:] != ids[:-1])[0]
        idx = np.r_[0, idchanges]
        if not len(ids):
            idx = idchanges

        # new values are where sorted labels change
        lchanges = llab(lab, slice(1, None)) != llab(lab, slice(None, -1))
        inc = np.r_[True, lchanges]
        if not len(val):
            inc = lchanges
        inc[idx] = True  # group boundaries are also new values
        out = np.diff(np.nonzero(np.r_[inc, True])[0])  # value counts

        # num. of times each group should be repeated
        rep = partial(np.repeat, repeats=np.add.reduceat(inc, idx))

        # multi-index components
        if isinstance(self._grouper.result_index, MultiIndex):
            codes = list(self._grouper.result_index.codes)
        else:
            codes = [
                algorithms.factorize(
                    self._grouper.result_index,
                    sort=self._grouper._sort,
                    use_na_sentinel=self._grouper.dropna,
                )[0]
            ]
        codes = [rep(level_codes) for level_codes in codes] + [llab(lab, inc)]
        levels = [*self._grouper.levels, lev]

        if dropna:
            mask = codes[-1] != -1
            if mask.all():
                dropna = False
            else:
                out, codes = out[mask], [level_codes[mask] for level_codes in codes]

        if normalize:
            out = out.astype("float")
            d = np.diff(np.r_[idx, len(ids)])
            if dropna:
                m = ids[lab == -1]
                np.add.at(d, m, -1)
                acc = rep(d)[mask]
            else:
                acc = rep(d)
            out /= acc

        if sort and bins is None:
            cat = ids[inc][mask] if dropna else ids[inc]
            sorter = np.lexsort((out if ascending else -out, cat))
            out, codes[-1] = out[sorter], codes[-1][sorter]

        if bins is not None:
            # for compat. with libgroupby.value_counts need to ensure every
            # bin is present at every index level, null filled with zeros
            diff = np.zeros(len(out), dtype="bool")
            for level_codes in codes[:-1]:
                diff |= np.r_[True, level_codes[1:] != level_codes[:-1]]

            ncat, nbin = diff.sum(), len(levels[-1])

            left = [np.repeat(np.arange(ncat), nbin), np.tile(np.arange(nbin), ncat)]

            right = [diff.cumsum() - 1, codes[-1]]

            # error: Argument 1 to "get_join_indexers" has incompatible type
            # "List[ndarray[Any, Any]]"; expected "List[Union[Union[ExtensionArray,
            # ndarray[Any, Any]], Index, Series]]
            _, idx = get_join_indexers(
                left,  # type: ignore[arg-type]
                right,
                sort=False,
                how="left",
            )
            if idx is not None:
                out = np.where(idx != -1, out[idx], 0)

            if sort:
                sorter = np.lexsort((out if ascending else -out, left[0]))
                out, left[-1] = out[sorter], left[-1][sorter]

            # build the multi-index w/ full levels
            def build_codes(lev_codes: np.ndarray) -> np.ndarray:
                return np.repeat(lev_codes[diff], nbin)

            codes = [build_codes(lev_codes) for lev_codes in codes[:-1]]
            codes.append(left[-1])

        mi = MultiIndex(
            levels=levels, codes=codes, names=index_names, verify_integrity=False
        )

        if is_integer_dtype(out.dtype):
            out = ensure_int64(out)
        result = self.obj._constructor(out, index=mi, name=name)
        if not self.as_index:
            result = result.reset_index()
        return result

    def take(
        self,
        indices: TakeIndexer,
        **kwargs,
    ) -> Series:
        """
        Return the elements in the given *positional* indices in each group.

        This means that we are not indexing according to actual values in
        the index attribute of the object. We are indexing according to the
        actual position of the element in the object.

        If a requested index does not exist for some group, this method will raise.
        To get similar behavior that ignores indices that don't exist, see
        :meth:`.SeriesGroupBy.nth`.

        Parameters
        ----------
        indices : array-like
            An array of ints indicating which positions to take in each group.

        **kwargs
            For compatibility with :meth:`numpy.take`. Has no effect on the
            output.

        Returns
        -------
        Series
            A Series containing the elements taken from each group.

        See Also
        --------
        Series.take : Take elements from a Series along an axis.
        Series.loc : Select a subset of a DataFrame by labels.
        Series.iloc : Select a subset of a DataFrame by positions.
        numpy.take : Take elements from an array along an axis.
        SeriesGroupBy.nth : Similar to take, won't raise if indices don't exist.

        Examples
        --------
        >>> df = pd.DataFrame(
        ...     [
        ...         ("falcon", "bird", 389.0),
        ...         ("parrot", "bird", 24.0),
        ...         ("lion", "mammal", 80.5),
        ...         ("monkey", "mammal", np.nan),
        ...         ("rabbit", "mammal", 15.0),
        ...     ],
        ...     columns=["name", "class", "max_speed"],
        ...     index=[4, 3, 2, 1, 0],
        ... )
        >>> df
             name   class  max_speed
        4  falcon    bird      389.0
        3  parrot    bird       24.0
        2    lion  mammal       80.5
        1  monkey  mammal        NaN
        0  rabbit  mammal       15.0
        >>> gb = df["name"].groupby([1, 1, 2, 2, 2])

        Take elements at rows 0 and 1 in each group.

        >>> gb.take([0, 1])
        1  4    falcon
           3    parrot
        2  2      lion
           1    monkey
        Name: name, dtype: str

        We may take elements using negative integers for positive indices,
        starting from the end of the object, just like with Python lists.

        >>> gb.take([-1, -2])
        1  3    parrot
           4    falcon
        2  0    rabbit
           1    monkey
        Name: name, dtype: str
        """
        result = self._op_via_apply("take", indices=indices, **kwargs)
        return result

    def skew(
        self,
        skipna: bool = True,
        numeric_only: bool = False,
        **kwargs,
    ) -> Series:
        """
        Return unbiased skew within groups.

        Normalized by N-1.

        Parameters
        ----------
        skipna : bool, default True
            Exclude NA/null values when computing the result.

        numeric_only : bool, default False
            Include only float, int, boolean columns. Not implemented for Series.

        **kwargs
            Additional keyword arguments to be passed to the function.

        Returns
        -------
        Series
            Unbiased skew within groups.

        See Also
        --------
        Series.skew : Return unbiased skew over requested axis.

        Examples
        --------
        >>> ser = pd.Series(
        ...     [390.0, 350.0, 357.0, np.nan, 22.0, 20.0, 30.0],
        ...     index=[
        ...         "Falcon",
        ...         "Falcon",
        ...         "Falcon",
        ...         "Falcon",
        ...         "Parrot",
        ...         "Parrot",
        ...         "Parrot",
        ...     ],
        ...     name="Max Speed",
        ... )
        >>> ser
        Falcon    390.0
        Falcon    350.0
        Falcon    357.0
        Falcon      NaN
        Parrot     22.0
        Parrot     20.0
        Parrot     30.0
        Name: Max Speed, dtype: float64
        >>> ser.groupby(level=0).skew()
        Falcon    1.525174
        Parrot    1.457863
        Name: Max Speed, dtype: float64
        >>> ser.groupby(level=0).skew(skipna=False)
        Falcon         NaN
        Parrot    1.457863
        Name: Max Speed, dtype: float64
        """

        return self._cython_agg_general(
            "skew", alt=None, skipna=skipna, numeric_only=numeric_only, **kwargs
        )

    def kurt(
        self,
        skipna: bool = True,
        numeric_only: bool = False,
        **kwargs,
    ) -> Series:
        """
        Return unbiased kurtosis within groups.

        Parameters
        ----------
        skipna : bool, default True
            Exclude NA/null values when computing the result.

        numeric_only : bool, default False
            Include only float, int, boolean columns. Not implemented for Series.

        **kwargs
            Additional keyword arguments to be passed to the function.

        Returns
        -------
        Series
            Unbiased kurtosis within groups.

        See Also
        --------
        Series.kurt : Return unbiased kurtosis over requested axis.

        Examples
        --------
        >>> ser = pd.Series(
        ...     [390.0, 350.0, 357.0, 333.0, np.nan, 22.0, 20.0, 30.0, 40.0, 41.0],
        ...     index=[
        ...         "Falcon",
        ...         "Falcon",
        ...         "Falcon",
        ...         "Falcon",
        ...         "Falcon",
        ...         "Parrot",
        ...         "Parrot",
        ...         "Parrot",
        ...         "Parrot",
        ...         "Parrot",
        ...     ],
        ...     name="Max Speed",
        ... )
        >>> ser
        Falcon    390.0
        Falcon    350.0
        Falcon    357.0
        Falcon    333.0
        Falcon      NaN
        Parrot     22.0
        Parrot     20.0
        Parrot     30.0
        Parrot     40.0
        Parrot     41.0
        Name: Max Speed, dtype: float64
        >>> ser.groupby(level=0).kurt()
        Falcon    1.622109
        Parrot   -2.878714
        Name: Max Speed, dtype: float64
        >>> ser.groupby(level=0).kurt(skipna=False)
        Falcon         NaN
        Parrot   -2.878714
        Name: Max Speed, dtype: float64
        """

        def alt(obj):
            # This should not be reached since the cython path should raise
            #  TypeError and not NotImplementedError.
            raise TypeError(f"'kurt' is not supported for dtype={obj.dtype}")

        return self._cython_agg_general(
            "kurt", alt=alt, skipna=skipna, numeric_only=numeric_only, **kwargs
        )

    @property
    def plot(self) -> GroupByPlot:
        """
        Make plots of groups from a Series.

        Uses the backend specified by the option ``plotting.backend``.
        By default, matplotlib is used.

        Returns
        -------
        GroupByPlot
            A plotting object that can be used to create plots for each group.

        See Also
        --------
        Series.plot : Make plots of Series.

        Examples
        --------
        >>> ser = pd.Series([1, 2, 3, 4, 5], index=["a", "a", "b", "b", "c"])
        >>> g = ser.groupby(level=0)
        >>> g.plot()  # doctest: +SKIP
        """
        result = GroupByPlot(self)
        return result

    def nlargest(
        self, n: int = 5, keep: Literal["first", "last", "all"] = "first"
    ) -> Series:
        """
        Return the largest `n` elements.

        Parameters
        ----------
        n : int, default 5
            Return this many descending sorted values.
        keep : {'first', 'last', 'all'}, default 'first'
            When there are duplicate values that cannot all fit in a
            Series of `n` elements:

            - ``first`` : return the first `n` occurrences in order
              of appearance.
            - ``last`` : return the last `n` occurrences in reverse
              order of appearance.
            - ``all`` : keep all occurrences. This can result in a Series of
              size larger than `n`.

        Returns
        -------
        Series
            The `n` largest values in the Series, sorted in decreasing order.

        See Also
        --------
        Series.nsmallest: Get the `n` smallest elements.
        Series.sort_values: Sort Series by values.
        Series.head: Return the first `n` rows.

        Notes
        -----
        Faster than ``.sort_values(ascending=False).head(n)`` for small `n`
        relative to the size of the ``Series`` object.

        Examples
        --------
        >>> s = pd.Series([1, 2, 3, 4, 5, 6])

        >>> s
        0    1
        1    2
        2    3
        3    4
        4    5
        5    6
        dtype: int64

        >>> s.groupby([1, 1, 1, 2, 2, 2]).nlargest(n=2)
        1  2    3
           1    2
        2  5    6
           4    5
        dtype: int64
        """
        f = partial(Series.nlargest, n=n, keep=keep)
        data = self._obj_with_exclusions
        # Don't change behavior if result index happens to be the same, i.e.
        # already ordered and n >= all group sizes.
        result = self._python_apply_general(f, data, not_indexed_same=True)
        return result

    def nsmallest(
        self, n: int = 5, keep: Literal["first", "last", "all"] = "first"
    ) -> Series:
        """
        Return the smallest `n` elements.

        Parameters
        ----------
        n : int, default 5
            Return this many ascending sorted values.
        keep : {'first', 'last', 'all'}, default 'first'
            When there are duplicate values that cannot all fit in a
            Series of `n` elements:

            - ``first`` : return the first `n` occurrences in order
              of appearance.
            - ``last`` : return the last `n` occurrences in reverse
              order of appearance.
            - ``all`` : keep all occurrences. This can result in a Series of
              size larger than `n`.

        Returns
        -------
        Series
            The `n` smallest values in the Series, sorted in increasing order.

        See Also
        --------
        Series.nlargest: Get the `n` largest elements.
        Series.sort_values: Sort Series by values.
        Series.head: Return the first `n` rows.

        Notes
        -----
        Faster than ``.sort_values().head(n)`` for small `n` relative to
        the size of the ``Series`` object.

        Examples
        --------
        >>> s = pd.Series([1, 2, 3, 4, 5, 6])

        >>> s
        0    1
        1    2
        2    3
        3    4
        4    5
        5    6
        dtype: int64

        >>> s.groupby([1, 1, 1, 2, 2, 2]).nsmallest(n=2)
        1  0    1
           1    2
        2  3    4
           4    5
        dtype: int64
        """
        f = partial(Series.nsmallest, n=n, keep=keep)
        data = self._obj_with_exclusions
        # Don't change behavior if result index happens to be the same, i.e.
        # already ordered and n >= all group sizes.
        result = self._python_apply_general(f, data, not_indexed_same=True)
        return result

    def idxmin(self, skipna: bool = True) -> Series:
        """
        Return the row label of the minimum value.

        If multiple values equal the minimum, the first row label with that
        value is returned.

        Parameters
        ----------
        skipna : bool, default True
            Exclude NA values.

        Returns
        -------
        Series
            Indexes of minima in each group.

        Raises
        ------
        ValueError
            When there are no valid values for a group. Then can happen if:

                * There is an unobserved group and ``observed=False``.
                * All values for a group are NA.
                * Some values for a group are NA and ``skipna=False``.

        .. versionchanged:: 3.0.0
            Previously if all values for a group are NA or some values for a group are
            NA and ``skipna=False``, this method would return NA. Now it raises instead.

        See Also
        --------
        numpy.argmin : Return indices of the minimum values
            along the given axis.
        DataFrame.idxmin : Return index of first occurrence of minimum
            over requested axis.
        Series.idxmax : Return index *label* of the first occurrence
            of maximum of values.

        Examples
        --------
        >>> ser = pd.Series(
        ...     [1, 2, 3, 4],
        ...     index=pd.DatetimeIndex(
        ...         ["2023-01-01", "2023-01-15", "2023-02-01", "2023-02-15"]
        ...     ),
        ... )
        >>> ser
        2023-01-01    1
        2023-01-15    2
        2023-02-01    3
        2023-02-15    4
        dtype: int64

        >>> ser.groupby(["a", "a", "b", "b"]).idxmin()
        a   2023-01-01
        b   2023-02-01
        dtype: datetime64[us]
        """
        return self._idxmax_idxmin("idxmin", skipna=skipna)

    def idxmax(self, skipna: bool = True) -> Series:
        """
        Return the row label of the maximum value.

        If multiple values equal the maximum, the first row label with that
        value is returned.

        Parameters
        ----------
        skipna : bool, default True
            Exclude NA values.

        Returns
        -------
        Series
            Indexes of maxima in each group.

        Raises
        ------
        ValueError
            When there are no valid values for a group. Then can happen if:

                * There is an unobserved group and ``observed=False``.
                * All values for a group are NA.
                * Some values for a group are NA and ``skipna=False``.

        .. versionchanged:: 3.0.0
            Previously if all values for a group are NA or some values for a group are
            NA and ``skipna=False``, this method would return NA. Now it raises instead.

        See Also
        --------
        numpy.argmax : Return indices of the maximum values
            along the given axis.
        DataFrame.idxmax : Return index of first occurrence of maximum
            over requested axis.
        Series.idxmin : Return index *label* of the first occurrence
            of minimum of values.

        Examples
        --------
        >>> ser = pd.Series(
        ...     [1, 2, 3, 4],
        ...     index=pd.DatetimeIndex(
        ...         ["2023-01-01", "2023-01-15", "2023-02-01", "2023-02-15"]
        ...     ),
        ... )
        >>> ser
        2023-01-01    1
        2023-01-15    2
        2023-02-01    3
        2023-02-15    4
        dtype: int64

        >>> ser.groupby(["a", "a", "b", "b"]).idxmax()
        a   2023-01-15
        b   2023-02-15
        dtype: datetime64[us]
        """
        return self._idxmax_idxmin("idxmax", skipna=skipna)

    def corr(
        self,
        other: Series,
        method: CorrelationMethod = "pearson",
        min_periods: int | None = None,
    ) -> Series:
        """
        Compute correlation between each group and another Series.

        Parameters
        ----------
        other : Series
            Series to compute correlation with.
        method : {'pearson', 'kendall', 'spearman'}, default 'pearson'
            Method of correlation to use.
        min_periods : int, optional
            Minimum number of observations required per pair of columns to
            have a valid result.

        Returns
        -------
        Series
            Correlation value for each group.

        See Also
        --------
        Series.corr : Equivalent method on ``Series``.

        Examples
        --------
        >>> s = pd.Series([1, 2, 3, 4], index=[0, 0, 1, 1])
        >>> g = s.groupby([0, 0, 1, 1])
        >>> g.corr()  # doctest: +SKIP
        """
        result = self._op_via_apply(
            "corr", other=other, method=method, min_periods=min_periods
        )
        return result

    def cov(
        self, other: Series, min_periods: int | None = None, ddof: int | None = 1
    ) -> Series:
        """
        Compute covariance between each group and another Series.

        Parameters
        ----------
        other : Series
            Series to compute covariance with.
        min_periods : int, optional
            Minimum number of observations required per pair of columns to
            have a valid result.
        ddof : int, optional
            Delta degrees of freedom for variance calculation.

        Returns
        -------
        Series
            Covariance value for each group.

        See Also
        --------
        Series.cov : Equivalent method on ``Series``.

        Examples
        --------
        >>> s = pd.Series([1, 2, 3, 4], index=[0, 0, 1, 1])
        >>> g = s.groupby([0, 0, 1, 1])
        >>> g.cov()  # doctest: +SKIP
        """
        result = self._op_via_apply(
            "cov", other=other, min_periods=min_periods, ddof=ddof
        )
        return result

    @property
    def is_monotonic_increasing(self) -> Series:
        """
        Return whether each group's values are monotonically increasing.

        Returns
        -------
        Series

        See Also
        --------
        SeriesGroupBy.is_monotonic_decreasing : Return whether each group's values
            are monotonically decreasing.

        Examples
        --------
        >>> s = pd.Series([2, 1, 3, 4], index=["Falcon", "Falcon", "Parrot", "Parrot"])
        >>> s.groupby(level=0).is_monotonic_increasing
        Falcon    False
        Parrot     True
        dtype: bool
        """
        return self.apply(lambda ser: ser.is_monotonic_increasing)

    @property
    def is_monotonic_decreasing(self) -> Series:
        """
        Return whether each group's values are monotonically decreasing.

        Returns
        -------
        Series

        See Also
        --------
        SeriesGroupBy.is_monotonic_increasing : Return whether each group's values
            are monotonically increasing.

        Examples
        --------
        >>> s = pd.Series([2, 1, 3, 4], index=["Falcon", "Falcon", "Parrot", "Parrot"])
        >>> s.groupby(level=0).is_monotonic_decreasing
        Falcon     True
        Parrot    False
        dtype: bool
        """
        return self.apply(lambda ser: ser.is_monotonic_decreasing)

    def hist(
        self,
        by=None,
        ax=None,
        grid: bool = True,
        xlabelsize: int | None = None,
        xrot: float | None = None,
        ylabelsize: int | None = None,
        yrot: float | None = None,
        figsize: tuple[float, float] | None = None,
        bins: int | Sequence[int] = 10,
        backend: str | None = None,
        legend: bool = False,
        **kwargs,
    ):
        """
        Draw histogram for each group's values using :meth:`Series.hist` API.

        Parameters
        ----------
        by : object, optional
            Grouping key.
        ax : matplotlib.axes.Axes, optional
            Axis to draw the histogram on.
        grid : bool, default True
            Show axis grid lines.
        xlabelsize : int, default None
            X axis label size.
        xrot : float, default None
            Rotation for x ticks.
        ylabelsize : int, default None
            Y axis label size.
        yrot : float, default None
            Rotation for y ticks.
        figsize : tuple, optional
            Figure size in inches.
        bins : int or sequence, default 10
            Number of histogram bins or bin edges.
        backend : str or callable or None, optional
            Plotting backend to use (e.g. 'matplotlib'). If None, use the default
            plotting backend.
        legend : bool, default False
            Whether to draw the legend.
        **kwargs
            Additional keyword arguments passed to :meth:`Series.hist`.

        Returns
        -------
        matplotlib.axes.Axes or ndarray of Axes
            The returned matplotlib axes or array of axes depending on input.

        See Also
        --------
        Series.hist : Equivalent histogram plotting method on Series.

        Examples
        --------
        >>> df = pd.DataFrame({"val": [1, 2, 2, 3, 3, 3]}, index=[0, 0, 1, 1, 2, 2])
        >>> g = df["val"].groupby([0, 0, 1, 1, 2, 2])
        >>> g.hist()  # doctest: +SKIP
        """
        result = self._op_via_apply(
            "hist",
            by=by,
            ax=ax,
            grid=grid,
            xlabelsize=xlabelsize,
            xrot=xrot,
            ylabelsize=ylabelsize,
            yrot=yrot,
            figsize=figsize,
            bins=bins,
            backend=backend,
            legend=legend,
            **kwargs,
        )
        return result

    @property
    def dtype(self) -> Series:
        """
        Return the dtype object of the underlying data for each group.

        Mirrors :meth:`Series.dtype` applied group-wise.

        Returns
        -------
        Series
            Dtype of each group's values.
        """
        return self.apply(lambda ser: ser.dtype)

    def unique(self) -> Series:
        """
        Return unique values for each group.

        It returns unique values for each of the grouped values. Returned in
        order of appearance. Hash table-based unique, therefore does NOT sort.

        Returns
        -------
        Series
            Unique values for each of the grouped values.

        See Also
        --------
        Series.unique : Return unique values of Series object.

        Examples
        --------
        >>> df = pd.DataFrame(
        ...     [
        ...         ("Chihuahua", "dog", 6.1),
        ...         ("Beagle", "dog", 15.2),
        ...         ("Chihuahua", "dog", 6.9),
        ...         ("Persian", "cat", 9.2),
        ...         ("Chihuahua", "dog", 7),
        ...         ("Persian", "cat", 8.8),
        ...     ],
        ...     columns=["breed", "animal", "height_in"],
        ... )
        >>> df
               breed     animal   height_in
        0  Chihuahua        dog         6.1
        1     Beagle        dog        15.2
        2  Chihuahua        dog         6.9
        3    Persian        cat         9.2
        4  Chihuahua        dog         7.0
        5    Persian        cat         8.8
        >>> ser = df.groupby("animal")["breed"].unique()
        >>> ser
        animal
        cat              [Persian]
        dog    [Chihuahua, Beagle]
        Name: breed, dtype: object
        """
        result = self._op_via_apply("unique")
        return result


@set_module("pandas.api.typing")
class DataFrameGroupBy(GroupBy[DataFrame]):
    def aggregate(self, func=None, *args, engine=None, engine_kwargs=None, **kwargs):
        """
        Aggregate using one or more operations.

        The ``aggregate`` function allows the application of one or more aggregation
        operations on groups of data within a DataFrameGroupBy object. It supports
        various aggregation methods, including user-defined functions and predefined
        functions such as 'sum', 'mean', etc.

        Parameters
        ----------
        func : function, str, list, dict or None
            Function to use for aggregating the data. If a function, must either
            work when passed a DataFrame or when passed to DataFrame.apply.

            Accepted combinations are:

            - function
            - string function name
            - list of functions and/or function names, e.g. ``[np.sum, 'mean']``
            - dict of index labels -> functions, function names or list of such.
            - None, in which case ``**kwargs`` are used with Named Aggregation. Here the
              output has one column for each element in ``**kwargs``. The name of the
              column is keyword, whereas the value determines the aggregation used to
              compute the values in the column.

              Can also accept a Numba JIT function with
              ``engine='numba'`` specified. Only passing a single function is supported
              with this engine.

              If the ``'numba'`` engine is chosen, the function must be
              a user defined function with ``values`` and ``index`` as the
              first and second arguments respectively in the function signature.
              Each group's index will be passed to the user defined function
              and optionally available for use.

        *args
            Positional arguments to pass to func.
        engine : str, default None
            * ``'cython'`` : Runs the function through C-extensions from cython.
            * ``'numba'`` : Runs the function through JIT compiled code from numba.
            * ``None`` : Defaults to ``'cython'`` or globally setting
                ``compute.use_numba``

        engine_kwargs : dict, default None
            * For ``'cython'`` engine, there are no accepted ``engine_kwargs``
            * For ``'numba'`` engine, the engine can accept ``nopython``, ``nogil``
              and ``parallel`` dictionary keys. The values must either be ``True`` or
              ``False``. The default ``engine_kwargs`` for the ``'numba'`` engine is
              ``{'nopython': True, 'nogil': False, 'parallel': False}`` and will be
              applied to the function

        **kwargs
            * If ``func`` is None, ``**kwargs`` are used to define the output names and
              aggregations via Named Aggregation. See ``func`` entry.
            * Otherwise, keyword arguments to be passed into func.

        Returns
        -------
        DataFrame
            Aggregated DataFrame based on the grouping and the applied aggregation
            functions.

        See Also
        --------
        DataFrame.groupby.apply : Apply function func group-wise
            and combine the results together.
        DataFrame.groupby.transform : Transforms the Series on each group
            based on the given function.
        DataFrame.aggregate : Aggregate using one or more operations.

        Notes
        -----
        When using ``engine='numba'``, there will be no "fall back" behavior internally.
        The group data and group index will be passed as numpy arrays to the JITed
        user defined function, and no alternative execution attempts will be tried.

        Functions that mutate the passed object can produce unexpected
        behavior or errors and are not supported. See :ref:`gotchas.udf-mutation`
        for more details.

        The resulting dtype will reflect the return value of the passed ``func``,
        see the examples below.

        Examples
        --------
        >>> data = {
        ...     "A": [1, 1, 2, 2],
        ...     "B": [1, 2, 3, 4],
        ...     "C": [0.362838, 0.227877, 1.267767, -0.562860],
        ... }
        >>> df = pd.DataFrame(data)
        >>> df
           A  B         C
        0  1  1  0.362838
        1  1  2  0.227877
        2  2  3  1.267767
        3  2  4 -0.562860

        The aggregation is for each column.

        >>> df.groupby("A").agg("min")
           B         C
        A
        1  1  0.227877
        2  3 -0.562860

        Multiple aggregations

        >>> df.groupby("A").agg(["min", "max"])
            B             C
          min max       min       max
        A
        1   1   2  0.227877  0.362838
        2   3   4 -0.562860  1.267767

        Select a column for aggregation

        >>> df.groupby("A").B.agg(["min", "max"])
           min  max
        A
        1    1    2
        2    3    4

        User-defined function for aggregation

        >>> df.groupby("A").agg(lambda x: sum(x) + 2)
            B          C
        A
        1       5       2.590715
        2       9       2.704907

        Different aggregations per column

        >>> df.groupby("A").agg({"B": ["min", "max"], "C": "sum"})
            B             C
          min max       sum
        A
        1   1   2  0.590715
        2   3   4  0.704907

        To control the output names with different aggregations per column,
        pandas supports "named aggregation"

        >>> df.groupby("A").agg(
        ...     b_min=pd.NamedAgg(column="B", aggfunc="min"),
        ...     c_sum=pd.NamedAgg(column="C", aggfunc="sum"),
        ... )
           b_min     c_sum
        A
        1      1  0.590715
        2      3  0.704907

        - The keywords are the *output* column names
        - The values are tuples whose first element is the column to select
          and the second element is the aggregation to apply to that column.
          Pandas provides the ``pandas.NamedAgg`` namedtuple with the fields
          ``['column', 'aggfunc']`` to make it clearer what the arguments are.
          As usual, the aggregation can be a callable or a string alias.

        See :ref:`groupby.aggregate.named` for more.

        The resulting dtype will reflect the return value of the aggregating
        function.

        >>> df.groupby("A")[["B"]].agg(lambda x: x.astype(float).min())
              B
        A
        1   1.0
        2   3.0
        """
        relabeling, func, columns, order = reconstruct_func(func, **kwargs)
        func = maybe_mangle_lambdas(func)

        if maybe_use_numba(engine):
            # Not all agg functions support numba, only propagate numba kwargs
            # if user asks for numba
            kwargs["engine"] = engine
            kwargs["engine_kwargs"] = engine_kwargs

        op = GroupByApply(self, func, args=args, kwargs=kwargs)
        result = op.agg()
        if not is_dict_like(func) and result is not None:
            # GH #52849
            if not self.as_index and is_list_like(func):
                return result.reset_index()
            else:
                return result
        elif relabeling:
            # this should be the only (non-raising) case with relabeling
            # used reordered index of columns
            result = cast(DataFrame, result)
            result = result.iloc[:, order]
            result = cast(DataFrame, result)
            # error: Incompatible types in assignment (expression has type
            # "Optional[List[str]]", variable has type
            # "Union[Union[Union[ExtensionArray, ndarray[Any, Any]],
            # Index, Series], Sequence[Any]]")
            result.columns = columns  # type: ignore[assignment]

        if result is None:
            # Remove the kwargs we inserted
            # (already stored in engine, engine_kwargs arguments)
            if "engine" in kwargs:
                del kwargs["engine"]
                del kwargs["engine_kwargs"]
            # at this point func is not a str, list-like, dict-like,
            # or a known callable(e.g. sum)
            if maybe_use_numba(engine):
                return self._aggregate_with_numba(
                    func, *args, engine_kwargs=engine_kwargs, **kwargs
                )
            # grouper specific aggregations
            if self._grouper.nkeys > 1:
                # test_groupby_as_index_series_scalar gets here with 'not self.as_index'
                return self._python_agg_general(func, *args, **kwargs)
            elif args or kwargs:
                # test_pass_args_kwargs gets here (with and without as_index)
                # can't return early
                result = self._aggregate_frame(func, *args, **kwargs)

            else:
                # try to treat as if we are passing a list
                gba = GroupByApply(self, [func], args=(), kwargs={})
                try:
                    result = gba.agg()

                except ValueError as err:
                    if "No objects to concatenate" not in str(err):
                        raise
                    # _aggregate_frame can fail with e.g. func=Series.mode,
                    # where it expects 1D values but would be getting 2D values
                    # In other tests, using aggregate_frame instead of GroupByApply
                    #  would give correct values but incorrect dtypes
                    #  object vs float64 in test_cython_agg_empty_buckets
                    #  float64 vs int64 in test_category_order_apply
                    result = self._aggregate_frame(func)

                else:
                    # GH#32040, GH#35246
                    # e.g. test_groupby_as_index_select_column_sum_empty_df
                    result = cast(DataFrame, result)
                    result.columns = self._obj_with_exclusions.columns.copy()

        if not self.as_index:
            result = self._insert_inaxis_grouper(result)
            result.index = default_index(len(result))

        return result

    agg = aggregate

    def _python_agg_general(self, func, *args, **kwargs):
        f = lambda x: func(x, *args, **kwargs)

        if self.ngroups == 0:
            # e.g. test_evaluate_with_empty_groups different path gets different
            #  result dtype in empty case.
            return self._python_apply_general(f, self._selected_obj, is_agg=True)

        obj = self._obj_with_exclusions

        if not len(obj.columns):
            # e.g. test_margins_no_values_no_cols
            return self._python_apply_general(f, self._selected_obj)

        output: dict[int, ArrayLike] = {}
        for idx, (name, ser) in enumerate(obj.items()):
            result = self._grouper.agg_series(ser, f)
            output[idx] = result

        res = self.obj._constructor(output)
        res.columns = obj.columns.copy(deep=False)
        return self._wrap_aggregated_output(res)

    def _aggregate_frame(self, func, *args, **kwargs) -> DataFrame:
        if self._grouper.nkeys != 1:
            raise AssertionError("Number of keys must be 1")

        obj = self._obj_with_exclusions

        result: dict[Hashable, NDFrame | np.ndarray] = {}
        for name, grp_df in self._grouper.get_iterator(obj):
            fres = func(grp_df, *args, **kwargs)
            result[name] = fres

        result_index = self._grouper.result_index
        out = self.obj._constructor(result, index=obj.columns, columns=result_index)
        out = out.T

        return out

    def _wrap_applied_output(
        self,
        data: DataFrame,
        values: list,
        not_indexed_same: bool = False,
        is_transform: bool = False,
    ):
        if len(values) == 0:
            if is_transform:
                # GH#47787 see test_group_on_empty_multiindex
                res_index = data.index
            elif not self.group_keys:
                res_index = None
            else:
                res_index = self._grouper.result_index

            result = self.obj._constructor(index=res_index, columns=data.columns)
            result = result.astype(data.dtypes)
            return result

        # GH12824
        # using values[0] here breaks test_groupby_apply_none_first
        first_not_none = next(com.not_none(*values), None)

        if first_not_none is None:
            # GH9684 - All values are None, return an empty frame
            # GH57775 - Ensure that columns and dtypes from original frame are kept.
            result = self.obj._constructor(columns=data.columns)
            result = result.astype(data.dtypes)
            return result
        elif isinstance(first_not_none, DataFrame):
            return self._concat_objects(
                values,
                not_indexed_same=not_indexed_same,
                is_transform=is_transform,
            )

        key_index = self._grouper.result_index if self.as_index else None

        if isinstance(first_not_none, (np.ndarray, Index)):
            # GH#1738: values is list of arrays of unequal lengths
            #  fall through to the outer else clause
            # TODO: sure this is right?  we used to do this
            #  after raising AttributeError above
            # GH 18930
            if not is_hashable(self._selection):
                # error: Need type annotation for "name"
                name = tuple(self._selection)  # type: ignore[var-annotated, arg-type]
            else:
                # error: Incompatible types in assignment
                # (expression has type "Hashable", variable
                # has type "Tuple[Any, ...]")
                name = self._selection  # type: ignore[assignment]
            return self.obj._constructor_sliced(values, index=key_index, name=name)
        elif not isinstance(first_not_none, Series):
            # values are not series or array-like but scalars
            # self._selection not passed through to Series as the
            # result should not take the name of original selection
            # of columns
            if self.as_index:
                return self.obj._constructor_sliced(values, index=key_index)
            else:
                result = self.obj._constructor(values, columns=[self._selection])
                result = self._insert_inaxis_grouper(result)
                return result
        else:
            # values are Series
            return self._wrap_applied_output_series(
                values,
                not_indexed_same,
                first_not_none,
                key_index,
                is_transform,
            )

    def _wrap_applied_output_series(
        self,
        values: list[Series],
        not_indexed_same: bool,
        first_not_none,
        key_index: Index | None,
        is_transform: bool,
    ) -> DataFrame | Series:
        kwargs = first_not_none._construct_axes_dict()
        backup = Series(**kwargs)
        values = [x if (x is not None) else backup for x in values]

        all_indexed_same = all_indexes_same(x.index for x in values)

        if not all_indexed_same:
            # GH 8467
            return self._concat_objects(
                values,
                not_indexed_same=True,
                is_transform=is_transform,
            )

        # Combine values
        # vstack+constructor is faster than concat and handles MI-columns
        stacked_values = np.vstack([np.asarray(v) for v in values])

        index = key_index
        columns = first_not_none.index.copy()
        if columns.name is None:
            # GH6124 - propagate name of Series when it's consistent
            names = {v.name for v in values}
            if len(names) == 1:
                columns.name = next(iter(names))

        if stacked_values.dtype == object:
            # We'll have the DataFrame constructor do inference
            stacked_values = stacked_values.tolist()
        result = self.obj._constructor(stacked_values, index=index, columns=columns)

        if not self.as_index:
            result = self._insert_inaxis_grouper(result)

        return result.__finalize__(self.obj, method="groupby")

    def _cython_transform(
        self,
        how: str,
        numeric_only: bool = False,
        **kwargs,
    ) -> DataFrame:
        # We have multi-block tests
        #  e.g. test_rank_min_int, test_cython_transform_frame
        #  test_transform_numeric_ret
        mgr: BlockManager = self._get_data_to_aggregate(
            numeric_only=numeric_only, name=how
        )

        def arr_func(bvalues: ArrayLike) -> ArrayLike:
            return self._grouper._cython_operation(
                "transform", bvalues, how, 1, **kwargs
            )

        res_mgr = mgr.apply(arr_func)

        res_df = self.obj._constructor_from_mgr(res_mgr, axes=res_mgr.axes)
        return res_df

    def _transform_general(self, func, engine, engine_kwargs, *args, **kwargs):
        if maybe_use_numba(engine):
            return self._transform_with_numba(
                func, *args, engine_kwargs=engine_kwargs, **kwargs
            )
        from pandas.core.reshape.concat import concat

        applied = []
        obj = self._obj_with_exclusions
        gen = self._grouper.get_iterator(obj)
        fast_path, slow_path = self._define_paths(func, *args, **kwargs)

        # Determine whether to use slow or fast path by evaluating on the first group.
        # Need to handle the case of an empty generator and process the result so that
        # it does not need to be computed again.
        try:
            name, group = next(gen)
        except StopIteration:
            pass
        else:
            # 2023-02-27 No tests broken by disabling this pinning
            object.__setattr__(group, "name", name)
            try:
                path, res = self._choose_path(fast_path, slow_path, group)
            except ValueError as err:
                # e.g. test_transform_with_non_scalar_group
                msg = "transform must return a scalar value for each group"
                raise ValueError(msg) from err
            if group.size > 0:
                res = _wrap_transform_general_frame(self.obj, group, res)
                applied.append(res)

        # Compute and process with the remaining groups
        for name, group in gen:
            if group.size == 0:
                continue
            # 2023-02-27 No tests broken by disabling this pinning
            object.__setattr__(group, "name", name)
            res = path(group)

            res = _wrap_transform_general_frame(self.obj, group, res)
            applied.append(res)

        concat_index = obj.columns
        concatenated = concat(
            applied, axis=0, verify_integrity=False, ignore_index=True
        )
        concatenated = concatenated.reindex(concat_index, axis=1)
        return self._set_result_index_ordered(concatenated)

    __examples_dataframe_doc = dedent(
        """
    >>> df = pd.DataFrame({'A' : ['foo', 'bar', 'foo', 'bar',
    ...                           'foo', 'bar'],
    ...                    'B' : ['one', 'one', 'two', 'three',
    ...                           'two', 'two'],
    ...                    'C' : [1, 5, 5, 2, 5, 5],
    ...                    'D' : [2.0, 5., 8., 1., 2., 9.]})
    >>> grouped = df.groupby('A')[['C', 'D']]
    >>> grouped.transform(lambda x: (x - x.mean()) / x.std())
            C         D
    0 -1.154701 -0.577350
    1  0.577350  0.000000
    2  0.577350  1.154701
    3 -1.154701 -1.000000
    4  0.577350 -0.577350
    5  0.577350  1.000000

    Broadcast result of the transformation

    >>> grouped.transform(lambda x: x.max() - x.min())
        C    D
    0  4.0  6.0
    1  3.0  8.0
    2  4.0  6.0
    3  3.0  8.0
    4  4.0  6.0
    5  3.0  8.0

    >>> grouped.transform("mean")
        C    D
    0  3.666667  4.0
    1  4.000000  5.0
    2  3.666667  4.0
    3  4.000000  5.0
    4  3.666667  4.0
    5  4.000000  5.0

    The resulting dtype will reflect the return value of the passed ``func``,
    for example:

    >>> grouped.transform(lambda x: x.astype(int).max())
    C  D
    0  5  8
    1  5  9
    2  5  8
    3  5  9
    4  5  8
    5  5  9
    """
    )

    def transform(self, func, *args, engine=None, engine_kwargs=None, **kwargs):
        """
        Call function producing a same-indexed DataFrame on each group.

        Returns a DataFrame having the same indexes as the original object
        filled with the transformed values.

        Parameters
        ----------
        func : function, str
            Function to apply to each group.
            See the Notes section below for requirements.

            Accepted inputs are:

            - String
            - Python function
            - Numba JIT function with ``engine='numba'`` specified.

            Only passing a single function is supported with this engine.
            If the ``'numba'`` engine is chosen, the function must be
            a user defined function with ``values`` and ``index`` as the
            first and second arguments respectively in the function signature.
            Each group's index will be passed to the user defined function
            and optionally available for use.

            If a string is chosen, then it needs to be the name
            of the groupby method you want to use.
        *args
            Positional arguments to pass to func.
        engine : str, default None
            * ``'cython'`` : Runs the function through C-extensions from cython.
            * ``'numba'`` : Runs the function through JIT compiled code from numba.
            * ``None`` : Defaults to ``'cython'``
              or the global setting ``compute.use_numba``

        engine_kwargs : dict, default None
            * For ``'cython'`` engine, there are no accepted ``engine_kwargs``
            * For ``'numba'`` engine, the engine can accept ``nopython``, ``nogil``
              and ``parallel`` dictionary keys. The values must either be ``True`` or
              ``False``. The default ``engine_kwargs`` for the ``'numba'`` engine is
              ``{'nopython': True, 'nogil': False, 'parallel': False}`` and will be
              applied to the function

        **kwargs
            Keyword arguments to be passed into func.

        Returns
        -------
        DataFrame
            DataFrame with the same indexes as the original object filled
            with transformed values.

        See Also
        --------
        DataFrame.groupby.apply : Apply function ``func`` group-wise and combine
            the results together.
        DataFrame.groupby.aggregate : Aggregate using one or more operations.
        DataFrame.transform : Call ``func`` on self producing a DataFrame with the
            same axis shape as self.

        Notes
        -----
        Each group is endowed the attribute 'name' in case you need to know
        which group you are working on.

        The current implementation imposes three requirements on f:

        * f must return a value that either has the same shape as the input
          subframe or can be broadcast to the shape of the input subframe.
          For example, if `f` returns a scalar it will be broadcast to have the
          same shape as the input subframe.
        * if this is a DataFrame, f must support application column-by-column
          in the subframe. If f also supports application to the entire subframe,
          then a fast path is used starting from the second chunk.
        * f must not mutate groups. Mutation is not supported and may
          produce unexpected results. See :ref:`gotchas.udf-mutation` for more details.

        When using ``engine='numba'``, there will be no "fall back" behavior internally.
        The group data and group index will be passed as numpy arrays to the JITed
        user defined function, and no alternative execution attempts will be tried.

        The resulting dtype will reflect the return value of the passed ``func``,
        see the examples below.

        .. versionchanged:: 2.0.0

            When using ``.transform`` on a grouped DataFrame
            and the transformation function returns a DataFrame,
            pandas now aligns the result's index with the input's index.
            You can call ``.to_numpy()`` on the result of the
            transformation function to avoid alignment.

        Examples
        --------

        >>> df = pd.DataFrame(
        ...     {
        ...         "A": ["foo", "bar", "foo", "bar", "foo", "bar"],
        ...         "B": ["one", "one", "two", "three", "two", "two"],
        ...         "C": [1, 5, 5, 2, 5, 5],
        ...         "D": [2.0, 5.0, 8.0, 1.0, 2.0, 9.0],
        ...     }
        ... )
        >>> grouped = df.groupby("A")[["C", "D"]]
        >>> grouped.transform(lambda x: (x - x.mean()) / x.std())
                  C         D
        0 -1.154701 -0.577350
        1  0.577350  0.000000
        2  0.577350  1.154701
        3 -1.154701 -1.000000
        4  0.577350 -0.577350
        5  0.577350  1.000000

        Broadcast result of the transformation

        >>> grouped.transform(lambda x: x.max() - x.min())
             C    D
        0  4.0  6.0
        1  3.0  8.0
        2  4.0  6.0
        3  3.0  8.0
        4  4.0  6.0
        5  3.0  8.0

        >>> grouped.transform("mean")
                  C    D
        0  3.666667  4.0
        1  4.000000  5.0
        2  3.666667  4.0
        3  4.000000  5.0
        4  3.666667  4.0
        5  4.000000  5.0

        The resulting dtype will reflect the return value of the passed ``func``,
        for example:

        >>> grouped.transform(lambda x: x.astype(int).max())
           C  D
        0  5  8
        1  5  9
        2  5  8
        3  5  9
        4  5  8
        5  5  9
        """
        return self._transform(
            func, *args, engine=engine, engine_kwargs=engine_kwargs, **kwargs
        )

    def _define_paths(self, func, *args, **kwargs):
        if isinstance(func, str):
            fast_path = lambda group: getattr(group, func)(*args, **kwargs)
            slow_path = lambda group: group.apply(
                lambda x: getattr(x, func)(*args, **kwargs), axis=0
            )
        else:
            fast_path = lambda group: func(group, *args, **kwargs)
            slow_path = lambda group: group.apply(
                lambda x: func(x, *args, **kwargs), axis=0
            )
        return fast_path, slow_path

    def _choose_path(self, fast_path: Callable, slow_path: Callable, group: DataFrame):
        path = slow_path
        res = slow_path(group)

        if self.ngroups == 1:
            # no need to evaluate multiple paths when only
            # a single group exists
            return path, res

        # if we make it here, test if we can use the fast path
        try:
            res_fast = fast_path(group)
        except AssertionError:
            raise  # pragma: no cover
        except Exception:
            # GH#29631 For user-defined function, we can't predict what may be
            #  raised; see test_transform.test_transform_fastpath_raises
            return path, res

        # verify fast path returns either:
        # a DataFrame with columns equal to group.columns
        # OR a Series with index equal to group.columns
        if isinstance(res_fast, DataFrame):
            if not res_fast.columns.equals(group.columns):
                return path, res
        elif isinstance(res_fast, Series):
            if not res_fast.index.equals(group.columns):
                return path, res
        else:
            return path, res

        if res_fast.equals(res):
            path = fast_path

        return path, res

    def filter(self, func, dropna: bool = True, *args, **kwargs) -> DataFrame:
        """
        Filter elements from groups that don't satisfy a criterion.

        Elements from groups are filtered if they do not satisfy the
        boolean criterion specified by func.

        Parameters
        ----------
        func : function
            Criterion to apply to each group. Should return True or False.
        dropna : bool
            Drop groups that do not pass the filter. True by default; if False,
            groups that evaluate False are filled with NaNs.
        *args : tuple
            Additional positional arguments to pass to `func`.
        **kwargs : dict
            Additional keyword arguments to pass to `func`.

        Returns
        -------
        DataFrame
            The filtered subset of the original DataFrame.

        See Also
        --------
        DataFrame.filter: Filter elements of ungrouped DataFrame.
        SeriesGroupBy.filter : Filter elements from groups base on criterion.

        Notes
        -----
        Each subframe is endowed the attribute 'name' in case you need to know
        which group you are working on.

        Functions that mutate the passed object can produce unexpected
        behavior or errors and are not supported. See :ref:`gotchas.udf-mutation`
        for more details.

        Examples
        --------
        >>> df = pd.DataFrame(
        ...     {
        ...         "A": ["foo", "bar", "foo", "bar", "foo", "bar"],
        ...         "B": [1, 2, 3, 4, 5, 6],
        ...         "C": [2.0, 5.0, 8.0, 1.0, 2.0, 9.0],
        ...     }
        ... )
        >>> grouped = df.groupby("A")
        >>> grouped.filter(lambda x: x["B"].mean() > 3.0)
             A  B    C
        1  bar  2  5.0
        3  bar  4  1.0
        5  bar  6  9.0
        """
        indices = []

        obj = self._selected_obj
        gen = self._grouper.get_iterator(obj)

        for name, group in gen:
            # 2023-02-27 no tests are broken this pinning, but it is documented in the
            #  docstring above.
            object.__setattr__(group, "name", name)

            res = func(group, *args, **kwargs)

            try:
                res = res.squeeze()
            except AttributeError:  # allow e.g., scalars and frames to pass
                pass

            # interpret the result of the filter
            if is_bool(res) or (is_scalar(res) and isna(res)):
                if notna(res) and res:
                    indices.append(self._get_index(name))
            else:
                # non scalars aren't allowed
                raise TypeError(
                    f"filter function returned a {type(res).__name__}, "
                    "but expected a scalar bool"
                )

        return self._apply_filter(indices, dropna)

    def __getitem__(self, key) -> DataFrameGroupBy | SeriesGroupBy:
        # per GH 23566
        if isinstance(key, tuple) and len(key) > 1:
            # if len == 1, then it becomes a SeriesGroupBy and this is actually
            # valid syntax, so don't raise
            raise ValueError(
                "Cannot subset columns with a tuple with more than one element. "
                "Use a list instead."
            )
        return super().__getitem__(key)

    def _gotitem(self, key, ndim: int, subset=None):
        """
        sub-classes to define
        return a sliced object

        Parameters
        ----------
        key : string / list of selections
        ndim : {1, 2}
            requested ndim of result
        subset : object, default None
            subset to act on
        """
        if ndim == 2:
            if subset is None:
                subset = self.obj
            return DataFrameGroupBy(
                subset,
                self.keys,
                level=self.level,
                grouper=self._grouper,
                exclusions=self.exclusions,
                selection=key,
                as_index=self.as_index,
                sort=self.sort,
                group_keys=self.group_keys,
                observed=self.observed,
                dropna=self.dropna,
            )
        elif ndim == 1:
            if subset is None:
                subset = self.obj[key]
            return SeriesGroupBy(
                subset,
                self.keys,
                level=self.level,
                grouper=self._grouper,
                exclusions=self.exclusions,
                selection=key,
                as_index=self.as_index,
                sort=self.sort,
                group_keys=self.group_keys,
                observed=self.observed,
                dropna=self.dropna,
            )

        raise AssertionError("invalid ndim for _gotitem")

    def _get_data_to_aggregate(
        self, *, numeric_only: bool = False, name: str | None = None
    ) -> BlockManager:
        obj = self._obj_with_exclusions
        mgr = obj._mgr
        if numeric_only:
            mgr = mgr.get_numeric_data()
        return mgr

    def _wrap_agged_manager(self, mgr: BlockManager) -> DataFrame:
        return self.obj._constructor_from_mgr(mgr, axes=mgr.axes)

    def _apply_to_column_groupbys(self, func) -> DataFrame:
        from pandas.core.reshape.concat import concat

        obj = self._obj_with_exclusions
        columns = obj.columns
        sgbs = (
            SeriesGroupBy(
                obj.iloc[:, i],
                selection=colname,
                grouper=self._grouper,
                exclusions=self.exclusions,
                observed=self.observed,
            )
            for i, colname in enumerate(obj.columns)
        )
        results = [func(sgb) for sgb in sgbs]

        if not results:
            # concat would raise
            res_df = DataFrame([], columns=columns, index=self._grouper.result_index)
        else:
            res_df = concat(results, keys=columns, axis=1)

        if not self.as_index:
            res_df.index = default_index(len(res_df))
            res_df = self._insert_inaxis_grouper(res_df)
        return res_df

    def nunique(self, dropna: bool = True) -> DataFrame:
        """
        Return DataFrame with counts of unique elements in each position.

        Parameters
        ----------
        dropna : bool, default True
            Don't include NaN in the counts.

        Returns
        -------
        nunique: DataFrame
            Counts of unique elements in each position.

        See Also
        --------
        DataFrame.nunique : Count number of distinct elements in specified axis.

        Examples
        --------
        >>> df = pd.DataFrame(
        ...     {
        ...         "id": ["spam", "egg", "egg", "spam", "ham", "ham"],
        ...         "value1": [1, 5, 5, 2, 5, 5],
        ...         "value2": list("abbaxy"),
        ...     }
        ... )
        >>> df
             id  value1 value2
        0  spam       1      a
        1   egg       5      b
        2   egg       5      b
        3  spam       2      a
        4   ham       5      x
        5   ham       5      y

        >>> df.groupby("id").nunique()
              value1  value2
        id
        egg        1       1
        ham        1       2
        spam       2       1

        Check for rows with the same id but conflicting values:

        >>> df.groupby("id").filter(lambda g: (g.nunique() > 1).any())
             id  value1 value2
        0  spam       1      a
        3  spam       2      a
        4   ham       5      x
        5   ham       5      y
        """
        return self._apply_to_column_groupbys(lambda sgb: sgb.nunique(dropna))

    def idxmax(
        self,
        skipna: bool = True,
        numeric_only: bool = False,
    ) -> DataFrame:
        """
        Return index of first occurrence of maximum in each group.

        Parameters
        ----------
        skipna : bool, default True
            Exclude NA values.
        numeric_only : bool, default False
            Include only `float`, `int` or `boolean` data.

        Returns
        -------
        DataFrame
            Indexes of maxima in each column according to the group.

        Raises
        ------
        ValueError
            When there are no valid values for a group. Then can happen if:

                * There is an unobserved group and ``observed=False``.
                * All values for a group are NA.
                * Some values for a group are NA and ``skipna=False``.

        .. versionchanged:: 3.0.0
            Previously if all values for a group are NA or some values for a group are
            NA and ``skipna=False``, this method would return NA. Now it raises instead.

        See Also
        --------
        Series.idxmax : Return index of the maximum element.
        DataFrame.idxmax : Indexes of maxima along the specified axis.

        Notes
        -----
        This method is the DataFrame version of ``ndarray.argmax``.

        Examples
        --------
        Consider a dataset containing food consumption in Argentina.

        >>> df = pd.DataFrame(
        ...     {
        ...         "consumption": [10.51, 103.11, 55.48],
        ...         "co2_emissions": [37.2, 19.66, 1712],
        ...         "food_type": ["meat", "plant", "meat"],
        ...     },
        ...     index=["Pork", "Wheat Products", "Beef"],
        ... )

        >>> df
                        consumption  co2_emissions food_type
        Pork                  10.51         37.20       meat
        Wheat Products       103.11         19.66      plant
        Beef                  55.48       1712.00       meat

        By default, it returns the index for the maximum value in each column
        according to the group.

        >>> df.groupby("food_type").idxmax()
                      consumption   co2_emissions
        food_type
        meat                 Beef            Beef
        plant      Wheat Products  Wheat Products
        """
        return self._idxmax_idxmin("idxmax", numeric_only=numeric_only, skipna=skipna)

    def idxmin(
        self,
        skipna: bool = True,
        numeric_only: bool = False,
    ) -> DataFrame:
        """
        Return index of first occurrence of minimum in each group.

        Parameters
        ----------
        skipna : bool, default True
            Exclude NA values.
        numeric_only : bool, default False
            Include only `float`, `int` or `boolean` data.

        Returns
        -------
        DataFrame
            Indexes of minima in each column according to the group.

        Raises
        ------
        ValueError
            When there are no valid values for a group. Then can happen if:

                * There is an unobserved group and ``observed=False``.
                * All values for a group are NA.
                * Some values for a group are NA and ``skipna=False``.

        .. versionchanged:: 3.0.0
            Previously if all values for a group are NA or some values for a group are
            NA and ``skipna=False``, this method would return NA. Now it raises instead.

        See Also
        --------
        Series.idxmin : Return index of the minimum element.
        DataFrame.idxmin : Indexes of minima along the specified axis.

        Notes
        -----
        This method is the DataFrame version of ``ndarray.argmin``.

        Examples
        --------
        Consider a dataset containing food consumption in Argentina.

        >>> df = pd.DataFrame(
        ...     {
        ...         "consumption": [10.51, 103.11, 55.48],
        ...         "co2_emissions": [37.2, 19.66, 1712],
        ...         "food_type": ["meat", "plant", "meat"],
        ...     },
        ...     index=["Pork", "Wheat Products", "Beef"],
        ... )

        >>> df
                        consumption  co2_emissions food_type
        Pork                  10.51         37.20       meat
        Wheat Products       103.11         19.66      plant
        Beef                  55.48       1712.00       meat

        By default, it returns the index for the minimum value in each column
        according to the group.

        >>> df.groupby("food_type").idxmin()
                      consumption   co2_emissions
        food_type
        meat                 Pork            Pork
        plant      Wheat Products  Wheat Products
        """
        return self._idxmax_idxmin("idxmin", numeric_only=numeric_only, skipna=skipna)

    boxplot = boxplot_frame_groupby

    def value_counts(
        self,
        subset: Sequence[Hashable] | None = None,
        normalize: bool = False,
        sort: bool = True,
        ascending: bool = False,
        dropna: bool = True,
    ) -> DataFrame | Series:
        """
        Return a Series or DataFrame containing counts of unique rows.

        Parameters
        ----------
        subset : list-like, optional
            Columns to use when counting unique combinations.
        normalize : bool, default False
            Return proportions rather than frequencies.
        sort : bool, default True
            Stable sort by frequencies when True. When False, non-grouping
            columns will appear in the order they occur in within groups.

            .. versionchanged:: 3.0.0

                In prior versions, ``sort=False`` would sort the non-grouping columns
                by label.
        ascending : bool, default False
            Sort in ascending order.
        dropna : bool, default True
            Don't include counts of rows that contain NA values.

        Returns
        -------
        Series or DataFrame
            Series if the groupby ``as_index`` is True, otherwise DataFrame.

        See Also
        --------
        Series.value_counts: Equivalent method on Series.
        DataFrame.value_counts: Equivalent method on DataFrame.
        SeriesGroupBy.value_counts: Equivalent method on SeriesGroupBy.

        Notes
        -----
        - If the groupby ``as_index`` is True then the returned Series will have a
          MultiIndex with one level per input column.
        - If the groupby ``as_index`` is False then the returned DataFrame will have an
          additional column with the value_counts. The column is labelled 'count' or
          'proportion', depending on the ``normalize`` parameter.

        By default, rows that contain any NA values are omitted from
        the result.

        By default, the result will be in descending order so that the
        first element of each group is the most frequently-occurring row.

        Examples
        --------
        >>> df = pd.DataFrame(
        ...     {
        ...         "gender": ["male", "male", "female", "male", "female", "male"],
        ...         "education": ["low", "medium", "high", "low", "high", "low"],
        ...         "country": ["US", "FR", "US", "FR", "FR", "FR"],
        ...     }
        ... )

        >>> df
                gender  education   country
        0       male    low         US
        1       male    medium      FR
        2       female  high        US
        3       male    low         FR
        4       female  high        FR
        5       male    low         FR

        >>> df.groupby("gender").value_counts()
        gender  education  country
        female  high       US         1
                           FR         1
        male    low        FR         2
                           US         1
                medium     FR         1
        Name: count, dtype: int64

        >>> df.groupby("gender").value_counts(ascending=True)
        gender  education  country
        female  high       US         1
                           FR         1
        male    low        US         1
                medium     FR         1
                low        FR         2
        Name: count, dtype: int64

        >>> df.groupby("gender").value_counts(normalize=True)
        gender  education  country
        female  high       US         0.50
                           FR         0.50
        male    low        FR         0.50
                           US         0.25
                medium     FR         0.25
        Name: proportion, dtype: float64

        >>> df.groupby("gender", as_index=False).value_counts()
           gender education country  count
        0  female      high      US      1
        1  female      high      FR      1
        2    male       low      FR      2
        3    male       low      US      1
        4    male    medium      FR      1

        >>> df.groupby("gender", as_index=False).value_counts(normalize=True)
           gender education country  proportion
        0  female      high      US        0.50
        1  female      high      FR        0.50
        2    male       low      FR        0.50
        3    male       low      US        0.25
        4    male    medium      FR        0.25
        """
        return self._value_counts(subset, normalize, sort, ascending, dropna)

    def take(
        self,
        indices: TakeIndexer,
        **kwargs,
    ) -> DataFrame:
        """
        Return the elements in the given *positional* indices in each group.

        This means that we are not indexing according to actual values in
        the index attribute of the object. We are indexing according to the
        actual position of the element in the object.

        If a requested index does not exist for some group, this method will raise.
        To get similar behavior that ignores indices that don't exist, see
        :meth:`.DataFrameGroupBy.nth`.

        Parameters
        ----------
        indices : array-like
            An array of ints indicating which positions to take.

        **kwargs
            For compatibility with :meth:`numpy.take`. Has no effect on the
            output.

        Returns
        -------
        DataFrame
            A DataFrame containing the elements taken from each group.

        See Also
        --------
        DataFrame.take : Take elements from a Series along an axis.
        DataFrame.loc : Select a subset of a DataFrame by labels.
        DataFrame.iloc : Select a subset of a DataFrame by positions.
        numpy.take : Take elements from an array along an axis.

        Examples
        --------
        >>> df = pd.DataFrame(
        ...     [
        ...         ("falcon", "bird", 389.0),
        ...         ("parrot", "bird", 24.0),
        ...         ("lion", "mammal", 80.5),
        ...         ("monkey", "mammal", np.nan),
        ...         ("rabbit", "mammal", 15.0),
        ...     ],
        ...     columns=["name", "class", "max_speed"],
        ...     index=[4, 3, 2, 1, 0],
        ... )
        >>> df
             name   class  max_speed
        4  falcon    bird      389.0
        3  parrot    bird       24.0
        2    lion  mammal       80.5
        1  monkey  mammal        NaN
        0  rabbit  mammal       15.0
        >>> gb = df.groupby([1, 1, 2, 2, 2])

        Take elements at rows 0 and 1.

        Note how the indices selected in the result do not correspond to
        our input indices 0 and 1. That's because we are selecting the 0th
        and 1st rows, not rows whose indices equal 0 and 1.

        >>> gb.take([0, 1])
               name   class  max_speed
        1 4  falcon    bird      389.0
          3  parrot    bird       24.0
        2 2    lion  mammal       80.5
          1  monkey  mammal        NaN

        The order of the specified indices influences the order in the result.
        Here, the order is swapped from the previous example.

        >>> gb.take([1, 0])
               name   class  max_speed
        1 3  parrot    bird       24.0
          4  falcon    bird      389.0
        2 1  monkey  mammal        NaN
          2    lion  mammal       80.5

        We may take elements using negative integers for positive indices,
        starting from the end of the object, just like with Python lists.

        >>> gb.take([-1, -2])
               name   class  max_speed
        1 3  parrot    bird       24.0
          4  falcon    bird      389.0
        2 0  rabbit  mammal       15.0
          1  monkey  mammal        NaN
        """
        result = self._op_via_apply("take", indices=indices, **kwargs)
        return result

    def skew(
        self,
        skipna: bool = True,
        numeric_only: bool = False,
        **kwargs,
    ) -> DataFrame:
        """
        Return unbiased skew within groups.

        Normalized by N-1.

        Parameters
        ----------
        skipna : bool, default True
            Exclude NA/null values when computing the result.

        numeric_only : bool, default False
            Include only float, int, boolean columns.

        **kwargs
            Additional keyword arguments to be passed to the function.

        Returns
        -------
        DataFrame
            Unbiased skew within groups.

        See Also
        --------
        DataFrame.skew : Return unbiased skew over requested axis.

        Examples
        --------
        >>> arrays = [
        ...     ["falcon", "parrot", "cockatoo", "kiwi", "lion", "monkey", "rabbit"],
        ...     ["bird", "bird", "bird", "bird", "mammal", "mammal", "mammal"],
        ... ]
        >>> index = pd.MultiIndex.from_arrays(arrays, names=("name", "class"))
        >>> df = pd.DataFrame(
        ...     {"max_speed": [389.0, 24.0, 70.0, np.nan, 80.5, 21.5, 15.0]},
        ...     index=index,
        ... )
        >>> df
                        max_speed
        name     class
        falcon   bird        389.0
        parrot   bird         24.0
        cockatoo bird         70.0
        kiwi     bird          NaN
        lion     mammal       80.5
        monkey   mammal       21.5
        rabbit   mammal       15.0
        >>> gb = df.groupby(["class"])
        >>> gb.skew()
                max_speed
        class
        bird     1.628296
        mammal   1.669046
        >>> gb.skew(skipna=False)
                max_speed
        class
        bird          NaN
        mammal   1.669046
        """

        def alt(obj):
            # This should not be reached since the cython path should raise
            #  TypeError and not NotImplementedError.
            raise TypeError(f"'skew' is not supported for dtype={obj.dtype}")

        return self._cython_agg_general(
            "skew", alt=alt, skipna=skipna, numeric_only=numeric_only, **kwargs
        )

    def kurt(
        self,
        skipna: bool = True,
        numeric_only: bool = False,
        **kwargs,
    ) -> DataFrame:
        """
        Return unbiased kurtosis within groups.

        Parameters
        ----------
        skipna : bool, default True
            Exclude NA/null values when computing the result.

        numeric_only : bool, default False
            Include only float, int, boolean columns.

        **kwargs
            Additional keyword arguments to be passed to the function.

        Returns
        -------
        DataFrame
            Unbiased kurtosis within groups.

        See Also
        --------
        DataFrame.kurt : Return unbiased kurtosis over requested axis.

        Examples
        --------
        >>> arrays = [
        ...     [
        ...         "falcon",
        ...         "parrot",
        ...         "cockatoo",
        ...         "kiwi",
        ...         "eagle",
        ...         "lion",
        ...         "monkey",
        ...         "rabbit",
        ...         "dog",
        ...         "wolf",
        ...     ],
        ...     [
        ...         "bird",
        ...         "bird",
        ...         "bird",
        ...         "bird",
        ...         "bird",
        ...         "mammal",
        ...         "mammal",
        ...         "mammal",
        ...         "mammal",
        ...         "mammal",
        ...     ],
        ... ]
        >>> index = pd.MultiIndex.from_arrays(arrays, names=("name", "class"))
        >>> df = pd.DataFrame(
        ...     {
        ...         "max_speed": [
        ...             389.0,
        ...             24.0,
        ...             70.0,
        ...             np.nan,
        ...             350.0,
        ...             80.5,
        ...             21.5,
        ...             15.0,
        ...             40.0,
        ...             50.0,
        ...         ]
        ...     },
        ...     index=index,
        ... )
        >>> df
                         max_speed
        name     class
        falcon   bird        389.0
        parrot   bird         24.0
        cockatoo bird         70.0
        kiwi     bird          NaN
        eagle    bird        350.0
        lion     mammal       80.5
        monkey   mammal       21.5
        rabbit   mammal       15.0
        dog      mammal       40.0
        wolf     mammal       50.0
        >>> gb = df.groupby(["class"])
        >>> gb.kurt()
                max_speed
        class
        bird    -5.493277
        mammal   0.204125
        >>> gb.kurt(skipna=False)
                max_speed
        class
        bird          NaN
        mammal   0.204125
        """

        return self._cython_agg_general(
            "kurt", alt=None, skipna=skipna, numeric_only=numeric_only, **kwargs
        )

    @property
    def plot(self) -> GroupByPlot:
        """
        Make plots of groups from a DataFrame.

        Uses the backend specified by the option ``plotting.backend``.
        By default, matplotlib is used.

        Returns
        -------
        GroupByPlot
            A plotting object that can be used to create plots for each group.

        See Also
        --------
        DataFrame.plot : Make plots of DataFrame.

        Examples
        --------
        >>> df = pd.DataFrame(
        ...     {"A": [1, 2, 3, 4], "B": [5, 6, 7, 8]}, index=["a", "a", "b", "b"]
        ... )
        >>> g = df.groupby(level=0)
        >>> g.plot()  # doctest: +SKIP
        """
        result = GroupByPlot(self)
        return result

    def corr(
        self,
        method: str | Callable[[np.ndarray, np.ndarray], float] = "pearson",
        min_periods: int = 1,
        numeric_only: bool = False,
    ) -> DataFrame:
        """
        Compute pairwise correlation of columns, excluding NA/null values.

        Parameters
        ----------
        method : {'pearson', 'kendall', 'spearman'} or callable
            Method of correlation:

            * pearson : standard correlation coefficient
            * kendall : Kendall Tau correlation coefficient
            * spearman : Spearman rank correlation
            * callable: callable with input two 1d ndarrays
                and returning a float. Note that the returned matrix from corr
                will have 1 along the diagonals and will be symmetric
                regardless of the callable's behavior.
        min_periods : int, optional
            Minimum number of observations required per pair of columns
            to have a valid result. Currently only available for Pearson
            and Spearman correlation.
        numeric_only : bool, default False
            Include only `float`, `int` or `boolean` data.

            .. versionchanged:: 2.0.0
                The default value of ``numeric_only`` is now ``False``.

        Returns
        -------
        DataFrame
            Correlation matrix.

        See Also
        --------
        DataFrame.corrwith : Compute pairwise correlation with another
            DataFrame or Series.
        Series.corr : Compute the correlation between two Series.

        Notes
        -----
        Pearson, Kendall and Spearman correlation are currently computed using
        pairwise complete observations.

        * `Pearson correlation coefficient <https://en.wikipedia.org/wiki/Pearson_correlation_coefficient>`_
        * `Kendall rank correlation coefficient <https://en.wikipedia.org/wiki/Kendall_rank_correlation_coefficient>`_
        * `Spearman's rank correlation coefficient <https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient>`_

        Examples
        --------
        >>> df = pd.DataFrame(
        ...     {
        ...         "age": [2, 3, 4, 6, 6, 1, 2, 1],
        ...         "weight": [2.1, 3.2, 4.1, 6.5, 3.3, 2.1, 4.1, 1.9],
        ...         "pet": ["dog", "cat", "dog", "cat", "dog", "cat", "dog", "cat"],
        ...     }
        ... )
        >>> df
           age  weight  pet
        0    2     2.1  dog
        1    3     3.2  cat
        2    4     4.1  dog
        3    6     6.5  cat
        4    6     3.3  dog
        5    1     2.1  cat
        6    2     4.1  dog
        7    1     1.9  cat
        >>> df.groupby("pet").corr()
                         age    weight
        pet
        cat age     1.000000  0.989321
            weight  0.989321  1.000000
        dog age     1.000000  0.184177
            weight  0.184177  1.000000
        """
        result = self._op_via_apply(
            "corr", method=method, min_periods=min_periods, numeric_only=numeric_only
        )
        return result

    def cov(
        self,
        min_periods: int | None = None,
        ddof: int | None = 1,
        numeric_only: bool = False,
    ) -> DataFrame:
        """
        Compute pairwise covariance of columns, excluding NA/null values.

        Compute the pairwise covariance among the series of a DataFrame.
        The returned data frame is the `covariance matrix
        <https://en.wikipedia.org/wiki/Covariance_matrix>`__ of the columns
        of the DataFrame.

        Both NA and null values are automatically excluded from the
        calculation. (See the note below about bias from missing values.)
        A threshold can be set for the minimum number of
        observations for each value created. Comparisons with observations
        below this threshold will be returned as ``NaN``.

        This method is generally used for the analysis of time series data to
        understand the relationship between different measures
        across time.

        Parameters
        ----------
        min_periods : int, optional
            Minimum number of observations required per pair of columns
            to have a valid result.

        ddof : int, default 1
            Delta degrees of freedom.  The divisor used in calculations
            is ``N - ddof``, where ``N`` represents the number of elements.
            This argument is applicable only when no ``nan`` is in the dataframe.

        numeric_only : bool, default False
            Include only `float`, `int` or `boolean` data.

            .. versionchanged:: 2.0.0
                The default value of ``numeric_only`` is now ``False``.

        Returns
        -------
        DataFrame
            The covariance matrix of the series of the DataFrame.

        See Also
        --------
        Series.cov : Compute covariance with another Series.
        core.window.ewm.ExponentialMovingWindow.cov : Exponential weighted sample
            covariance.
        core.window.expanding.Expanding.cov : Expanding sample covariance.
        core.window.rolling.Rolling.cov : Rolling sample covariance.

        Notes
        -----
        Returns the covariance matrix of the DataFrame's time series.
        The covariance is normalized by N-ddof.

        For DataFrames that have Series that are missing data (assuming that
        data is `missing at random
        <https://en.wikipedia.org/wiki/Missing_data#Missing_at_random>`__)
        the returned covariance matrix will be an unbiased estimate
        of the variance and covariance between the member Series.

        However, for many applications this estimate may not be acceptable
        because the estimate covariance matrix is not guaranteed to be positive
        semi-definite. This could lead to estimate correlations having
        absolute values which are greater than one, and/or a non-invertible
        covariance matrix. See `Estimation of covariance matrices
        <https://en.wikipedia.org/w/index.php?title=Estimation_of_covariance_
        matrices>`__ for more details.

        Examples
        --------
        >>> df = pd.DataFrame(
        ...     {
        ...         "age": [2, 3, 4, 6, 6, 1, 2, 1],
        ...         "weight": [2.1, 3.2, 4.1, 6.5, 3.3, 2.1, 4.1, 1.9],
        ...         "pet": ["dog", "cat", "dog", "cat", "dog", "cat", "dog", "cat"],
        ...     }
        ... )
        >>> df
           age  weight  pet
        0    2     2.1  dog
        1    3     3.2  cat
        2    4     4.1  dog
        3    6     6.5  cat
        4    6     3.3  dog
        5    1     2.1  cat
        6    2     4.1  dog
        7    1     1.9  cat
        >>> df.groupby("pet").cov()
                         age    weight
        pet
        cat age     5.583333  4.975000
            weight  4.975000  4.529167
        dog age     3.666667  0.333333
            weight  0.333333  0.893333
        """
        result = self._op_via_apply(
            "cov", min_periods=min_periods, ddof=ddof, numeric_only=numeric_only
        )
        return result

    def hist(
        self,
        column: IndexLabel | None = None,
        by=None,
        grid: bool = True,
        xlabelsize: int | None = None,
        xrot: float | None = None,
        ylabelsize: int | None = None,
        yrot: float | None = None,
        ax=None,
        sharex: bool = False,
        sharey: bool = False,
        figsize: tuple[float, float] | None = None,
        layout: tuple[int, int] | None = None,
        bins: int | Sequence[int] = 10,
        backend: str | None = None,
        legend: bool = False,
        **kwargs,
    ):
        """
        Make a histogram of the DataFrame's columns.

        A `histogram`_ is a representation of the distribution of data.
        This function calls :meth:`matplotlib.pyplot.hist`, on each series in
        the DataFrame, resulting in one histogram per column.

        .. _histogram: https://en.wikipedia.org/wiki/Histogram

        Parameters
        ----------
        column : str or sequence, optional
            If passed, will be used to limit data to a subset of columns.
        by : object, optional
            If passed, then used to form histograms for separate groups.
        grid : bool, default True
            Whether to show axis grid lines.
        xlabelsize : int, default None
            If specified changes the x-axis label size.
        xrot : float, default None
            Rotation of x axis labels. For example, a value of 90 displays the
            x labels rotated 90 degrees clockwise.
        ylabelsize : int, default None
            If specified changes the y-axis label size.
        yrot : float, default None
            Rotation of y axis labels. For example, a value of 90 displays the
            y labels rotated 90 degrees clockwise.
        ax : Matplotlib axes object, default None
            The axes to plot the histogram on.
        sharex : bool, default True if ax is None else False
            In case subplots=True, share x axis and set some x axis labels to
            invisible; defaults to True if ax is None otherwise False if an ax
            is passed in.
            Note that passing in both an ax and sharex=True will alter all x axis
            labels for all subplots in a figure.
        sharey : bool, default False
            In case subplots=True, share y axis and set some y axis labels to
            invisible.
        figsize : tuple, optional
            The size in inches of the figure to create. Uses the value in
            `matplotlib.rcParams` by default.
        layout : tuple, optional
            Tuple of (rows, columns) for the layout of the histograms.
        bins : int or sequence, default 10
            Number of histogram bins to be used. If an integer is given, bins + 1
            bin edges are calculated and returned. If bins is a sequence, gives
            bin edges, including left edge of first bin and right edge of last
            bin. In this case, bins is returned unmodified.

        backend : str, default None
            Backend to use instead of the backend specified in the option
            ``plotting.backend``. For instance, 'matplotlib'. Alternatively, to
            specify the ``plotting.backend`` for the whole session, set
            ``pd.options.plotting.backend``.

        legend : bool, default False
            Whether to show the legend.

        **kwargs
            All other plotting keyword arguments to be passed to
            :meth:`matplotlib.pyplot.hist`.

        Returns
        -------
        matplotlib.Axes or numpy.ndarray
            A ``matplotlib.Axes`` object or an array of ``Axes`` objects, depending on
            the layout and grouping.

        See Also
        --------
        matplotlib.pyplot.hist : Plot a histogram using matplotlib.

        Examples
        --------
        This example draws a histogram based on the length and width of
        some animals, displayed in three bins

        .. plot::
            :context: close-figs

            >>> data = {
            ...     "length": [1.5, 0.5, 1.2, 0.9, 3],
            ...     "width": [0.7, 0.2, 0.15, 0.2, 1.1],
            ... }
            >>> index = ["pig", "rabbit", "duck", "chicken", "horse"]
            >>> df = pd.DataFrame(data, index=index)
            >>> hist = df.groupby("length").hist(bins=3)
        """
        result = self._op_via_apply(
            "hist",
            column=column,
            by=by,
            grid=grid,
            xlabelsize=xlabelsize,
            xrot=xrot,
            ylabelsize=ylabelsize,
            yrot=yrot,
            ax=ax,
            sharex=sharex,
            sharey=sharey,
            figsize=figsize,
            layout=layout,
            bins=bins,
            backend=backend,
            legend=legend,
            **kwargs,
        )
        return result

    def corrwith(
        self,
        other: DataFrame | Series,
        drop: bool = False,
        method: CorrelationMethod = "pearson",
        numeric_only: bool = False,
    ) -> DataFrame:
        """
        Compute pairwise correlation.

        .. deprecated:: 3.0.0

        Pairwise correlation is computed between rows or columns of
        DataFrame with rows or columns of Series or DataFrame. DataFrames
        are first aligned along both axes before computing the
        correlations.

        Parameters
        ----------
        other : DataFrame, Series
            Object with which to compute correlations.
        drop : bool, default False
            Drop missing indices from result.
        method : {'pearson', 'kendall', 'spearman'} or callable
            Method of correlation:

            * pearson : standard correlation coefficient
            * kendall : Kendall Tau correlation coefficient
            * spearman : Spearman rank correlation
            * callable: callable with input two 1d ndarrays
                and returning a float.

        numeric_only : bool, default False
            Include only `float`, `int` or `boolean` data.

            .. versionchanged:: 2.0.0
                The default value of ``numeric_only`` is now ``False``.

        Returns
        -------
        Series
            Pairwise correlations.

        See Also
        --------
        DataFrame.corr : Compute pairwise correlation of columns.

        Examples
        --------
        >>> df1 = pd.DataFrame(
        ...     {
        ...         "Day": [1, 1, 1, 2, 2, 2, 3, 3, 3],
        ...         "Data": [6, 6, 8, 5, 4, 2, 7, 3, 9],
        ...     }
        ... )
        >>> df2 = pd.DataFrame(
        ...     {
        ...         "Day": [1, 1, 1, 2, 2, 2, 3, 3, 3],
        ...         "Data": [5, 3, 8, 3, 1, 1, 2, 3, 6],
        ...     }
        ... )

        >>> df1.groupby("Day").corrwith(df2)
                 Data  Day
        Day
        1    0.917663  NaN
        2    0.755929  NaN
        3    0.576557  NaN
        """
        warnings.warn(
            "DataFrameGroupBy.corrwith is deprecated",
            Pandas4Warning,
            stacklevel=find_stack_level(),
        )
        result = self._op_via_apply(
            "corrwith",
            other=other,
            drop=drop,
            method=method,
            numeric_only=numeric_only,
        )
        return result


def _wrap_transform_general_frame(
    obj: DataFrame, group: DataFrame, res: DataFrame | Series
) -> DataFrame:
    from pandas import concat

    if isinstance(res, Series):
        # we need to broadcast across the
        # other dimension; this will preserve dtypes
        # GH14457
        if res.index.is_(obj.index):
            res_frame = concat([res] * len(group.columns), axis=1, ignore_index=True)
            res_frame.columns = group.columns
            res_frame.index = group.index
        else:
            res_frame = obj._constructor(
                np.tile(res.values, (len(group.index), 1)),
                columns=group.columns,
                index=group.index,
            )
        assert isinstance(res_frame, DataFrame)
        return res_frame
    elif isinstance(res, DataFrame) and not res.index.is_(group.index):
        return res._align_frame(group)[0]
    else:
        return res
