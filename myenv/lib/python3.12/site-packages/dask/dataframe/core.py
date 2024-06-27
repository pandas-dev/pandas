from __future__ import annotations

import operator
import warnings
from collections.abc import Callable, Iterator, Mapping, Sequence
from functools import partial, wraps
from numbers import Integral, Number
from operator import getitem
from pprint import pformat
from typing import Any, ClassVar, Literal, cast

import numpy as np
import pandas as pd
from pandas.api.types import (
    is_bool_dtype,
    is_datetime64_any_dtype,
    is_extension_array_dtype,
    is_numeric_dtype,
    is_timedelta64_dtype,
)
from tlz import first, merge, partition_all, remove, unique

import dask.array as da
from dask import config, core
from dask.array.core import Array, normalize_arg
from dask.bag import map_partitions as map_bag_partitions
from dask.base import (
    DaskMethodsMixin,
    dont_optimize,
    is_dask_collection,
    named_schedulers,
    tokenize,
)
from dask.blockwise import Blockwise, BlockwiseDep, BlockwiseDepDict, blockwise
from dask.context import globalmethod
from dask.dataframe import methods
from dask.dataframe._compat import (
    PANDAS_GE_140,
    PANDAS_GE_150,
    PANDAS_GE_200,
    PANDAS_GE_210,
    PANDAS_GE_300,
    PANDAS_VERSION,
    check_convert_dtype_deprecation,
    check_nuisance_columns_warning,
    check_numeric_only_deprecation,
    check_reductions_runtime_warning,
    is_string_dtype,
)
from dask.dataframe.accessor import CachedAccessor, DatetimeAccessor, StringAccessor
from dask.dataframe.categorical import CategoricalAccessor, categorize
from dask.dataframe.dispatch import (
    get_parallel_type,
    group_split_dispatch,
    hash_object_dispatch,
    meta_nonempty,
)
from dask.dataframe.optimize import optimize
from dask.dataframe.utils import (
    AttributeNotImplementedError,
    check_matching_columns,
    check_numeric_only_valid,
    clear_known_categories,
    drop_by_shallow_copy,
    get_numeric_only_kwargs,
    has_known_categories,
    index_summary,
    insert_meta_param_description,
    is_dataframe_like,
    is_index_like,
    is_series_like,
    make_meta,
    meta_frame_constructor,
    meta_series_constructor,
    pyarrow_strings_enabled,
    raise_on_meta_error,
    valid_divisions,
)
from dask.delayed import Delayed, delayed, unpack_collections
from dask.highlevelgraph import HighLevelGraph
from dask.layers import DataFrameTreeReduction
from dask.typing import Graph, NestedKeys, no_default
from dask.utils import (
    F,
    IndexCallable,
    M,
    OperatorMethodMixin,
    _deprecated,
    _deprecated_kwarg,
    apply,
    derived_from,
    funcname,
    get_meta_library,
    has_keyword,
    is_arraylike,
    iter_chunks,
    key_split,
    maybe_pluralize,
    memory_repr,
    parse_bytes,
    partial_by_order,
    pseudorandom,
    put_lines,
    random_state_data,
    typename,
)
from dask.widgets import get_template

DEFAULT_GET = named_schedulers.get("threads", named_schedulers["sync"])

GROUP_KEYS_DEFAULT: bool | None = True
if PANDAS_GE_150 and not PANDAS_GE_200:
    GROUP_KEYS_DEFAULT = None

pd.set_option("compute.use_numexpr", False)


def _numeric_only(func):
    """Decorator for methods that accept a numeric_only kwarg"""

    @wraps(func)
    def wrapper(self, *args, **kwargs):
        # numeric_only is None by default - in that case self = self.
        if kwargs.get("numeric_only") is False:
            raise NotImplementedError(
                "'numeric_only=False' is not implemented in Dask."
            )
        elif kwargs.get("numeric_only") is True:
            self = self._get_numeric_data()
        return func(self, *args, **kwargs)

    return wrapper


def _numeric_data(func):
    """Modified version of the above decorator, right now only used with std. We don't
    need raising NotImplementedError there, because it's handled by
    _numeric_only_maybe_warn instead. This is a temporary solution that needs
    more time to be generalized."""

    @wraps(func)
    def wrapper(self, *args, **kwargs):
        if kwargs.get("numeric_only") is True:
            self = self._get_numeric_data()
        return func(self, *args, **kwargs)

    return wrapper


def _numeric_only_maybe_warn(df, numeric_only, default=None):
    """Update numeric_only to get rid of no_default, and possibly warn about default value.
    TODO: should move to numeric_only decorator. See https://github.com/dask/dask/pull/9952
    """
    if is_dataframe_like(df):
        warn_numeric_only = False
        if numeric_only is no_default:
            if PANDAS_GE_200:
                numeric_only = False
            else:
                warn_numeric_only = True

        numerics = df._meta._get_numeric_data()
        has_non_numerics = len(numerics.columns) < len(df._meta.columns)
        if has_non_numerics:
            if numeric_only is False:
                raise NotImplementedError(
                    "'numeric_only=False' is not implemented in Dask."
                )
            elif warn_numeric_only:
                warnings.warn(
                    "The default value of numeric_only in dask will be changed to False in "
                    "the future when using dask with pandas 2.0",
                    FutureWarning,
                )
    if numeric_only is no_default and default is not None:
        numeric_only = default
    return {} if numeric_only is no_default else {"numeric_only": numeric_only}


def _concat(args, ignore_index=False):
    if not args:
        return args
    if isinstance(first(core.flatten(args)), np.ndarray):
        return da.core.concatenate3(args)
    if not has_parallel_type(args[0]):
        try:
            return pd.Series(args)
        except Exception:
            return args
    # We filter out empty partitions here because pandas frequently has
    # inconsistent dtypes in results between empty and non-empty frames.
    # Ideally this would be handled locally for each operation, but in practice
    # this seems easier. TODO: don't do this.
    args2 = [i for i in args if len(i)]
    return (
        args[0]
        if not args2
        else methods.concat(args2, uniform=True, ignore_index=ignore_index)
    )


def _determine_split_out_shuffle(shuffle_method, split_out):
    """Determine the default shuffle behavior based on split_out"""
    if shuffle_method is None:
        if split_out > 1:
            # FIXME: This is using a different default but it is not fully
            # understood why this is a better choice.
            # For more context, see
            # https://github.com/dask/dask/pull/9826/files#r1072395307
            # https://github.com/dask/distributed/issues/5502
            return config.get("dataframe.shuffle.method", None) or "tasks"
        else:
            return False
    if shuffle_method is True:
        return config.get("dataframe.shuffle.method", None) or "tasks"
    return shuffle_method


def _dummy_numpy_dispatcher(
    *arg_names: Literal["dtype", "out"], deprecated: bool = False
) -> Callable[[F], F]:
    """Decorator to handle the out= and dtype= keyword arguments.

    These parameters are deprecated in all dask.dataframe reduction methods
    and will be soon completely disallowed.
    However, these methods must continue accepting 'out=None' and/or 'dtype=None'
    indefinitely in order to support numpy dispatchers. For example,
    ``np.mean(df)`` calls ``df.mean(out=None, dtype=None)``.

    Parameters
    ----------
    deprecated: bool
        If True, warn if not None and then pass the parameter to the wrapped function
        If False, raise error if not None; do not pass the parameter down.

    See Also
    --------
    _deprecated_kwarg
    """

    def decorator(func: F) -> F:
        @wraps(func)
        def wrapper(*args, **kwargs):
            for name in arg_names:
                if deprecated:
                    if kwargs.get(name, None) is not None:
                        warnings.warn(
                            f"the '{name}' keyword is deprecated and "
                            "will be removed in a future version.",
                            FutureWarning,
                            stacklevel=2,
                        )
                else:
                    if kwargs.pop(name, None) is not None:
                        raise ValueError(f"the '{name}' keyword is not supported")

            return func(*args, **kwargs)

        return cast(F, wrapper)

    return decorator


def finalize(results):
    return _concat(results)


class Scalar(DaskMethodsMixin, OperatorMethodMixin):
    """A Dask object to represent a pandas scalar"""

    def __init__(self, dsk, name, meta, divisions=None):
        # divisions is ignored, only present to be compatible with other
        # objects.
        if not isinstance(dsk, HighLevelGraph):
            dsk = HighLevelGraph.from_collections(name, dsk, dependencies=[])
        self.dask = dsk
        self._name = name
        self._parent_meta = pd.Series(dtype="float64")

        meta = make_meta(meta, parent_meta=self._parent_meta)
        if is_dataframe_like(meta) or is_series_like(meta) or is_index_like(meta):
            raise TypeError(
                f"Expected meta to specify scalar, got {typename(type(meta))}"
            )
        self._meta = meta

    def __dask_graph__(self) -> Graph:
        return self.dask

    def __dask_keys__(self) -> NestedKeys:
        return [self.key]

    def __dask_tokenize__(self):
        return self._name

    def __dask_layers__(self) -> Sequence[str]:
        return (self._name,)

    __dask_optimize__ = globalmethod(
        optimize, key="dataframe_optimize", falsey=dont_optimize
    )
    __dask_scheduler__ = staticmethod(DEFAULT_GET)

    def __dask_postcompute__(self):
        return first, ()

    def __dask_postpersist__(self):
        return self._rebuild, ()

    def _rebuild(self, dsk, *, rename=None):
        name = self._name
        if rename:
            name = rename.get(name, name)
        return Scalar(dsk, name, self._meta, self.divisions)

    @property
    def _meta_nonempty(self):
        return self._meta

    @property
    def dtype(self):
        return self._meta.dtype

    def __dir__(self):
        o = set(dir(type(self)))
        o.update(self.__dict__)
        if not hasattr(self._meta, "dtype"):
            o.remove("dtype")  # dtype only in `dir` if available
        return list(o)

    @property
    def divisions(self):
        """Dummy divisions to be compat with Series and DataFrame"""
        return (None, None)

    def __repr__(self):
        name = self._name if len(self._name) < 10 else self._name[:7] + "..."
        if hasattr(self._meta, "dtype"):
            extra = ", dtype=%s" % self._meta.dtype
        else:
            extra = ", type=%s" % type(self._meta).__name__
        return f"dd.Scalar<{name}{extra}>"

    def __array__(self):
        # array interface is required to support pandas instance + Scalar
        # Otherwise, above op results in pd.Series of Scalar (object dtype)
        return np.asarray(self.compute())

    @property
    def _args(self):
        return (self.dask, self._name, self._meta)

    def __getstate__(self):
        return self._args

    def __setstate__(self, state):
        self.dask, self._name, self._meta = state

    def __bool__(self):
        raise TypeError(
            f"Trying to convert {self} to a boolean value. Because Dask objects are "
            "lazily evaluated, they cannot be converted to a boolean value or used "
            "in boolean conditions like if statements. Try calling .compute() to "
            "force computation prior to converting to a boolean value or using in "
            "a conditional statement."
        )

    @property
    def key(self):
        return (self._name, 0)

    @classmethod
    def _get_unary_operator(cls, op):
        def f(self):
            name = funcname(op) + "-" + tokenize(self)
            dsk = {(name, 0): (op, (self._name, 0))}
            meta = op(self._meta_nonempty)
            graph = HighLevelGraph.from_collections(name, dsk, dependencies=[self])
            return Scalar(graph, name, meta)

        return f

    @classmethod
    def _get_binary_operator(cls, op, inv=False):
        return lambda self, other: _scalar_binary(op, self, other, inv=inv)

    def to_delayed(self, optimize_graph=True):
        """Convert into a ``dask.delayed`` object.

        Parameters
        ----------
        optimize_graph : bool, optional
            If True [default], the graph is optimized before converting into
            ``dask.delayed`` objects.
        """
        dsk = self.__dask_graph__()
        layer = self.__dask_layers__()[0]
        if optimize_graph:
            dsk = self.__dask_optimize__(dsk, self.__dask_keys__())
            layer = "delayed-" + self._name
            dsk = HighLevelGraph.from_collections(layer, dsk, dependencies=())
        return Delayed(self.key, dsk, layer=layer)


def _scalar_binary(op, self, other, inv=False):
    name = f"{funcname(op)}-{tokenize(self, other)}"
    dependencies = [self]

    dsk = {}
    return_type = get_parallel_type(other)

    if isinstance(other, Scalar):
        dependencies.append(other)
        other_key = (other._name, 0)
    elif is_dask_collection(other):
        return NotImplemented
    else:
        other_key = other

    dsk[(name, 0)] = (
        (op, other_key, (self._name, 0)) if inv else (op, (self._name, 0), other_key)
    )

    other_meta = make_meta(other, parent_meta=self._parent_meta)
    other_meta_nonempty = meta_nonempty(other_meta)
    if inv:
        meta = op(other_meta_nonempty, self._meta_nonempty)
    else:
        meta = op(self._meta_nonempty, other_meta_nonempty)

    graph = HighLevelGraph.from_collections(name, dsk, dependencies=dependencies)
    if return_type is not Scalar:
        return return_type(graph, name, meta, [other.index.min(), other.index.max()])
    else:
        return Scalar(graph, name, meta)


class _Frame(DaskMethodsMixin, OperatorMethodMixin):
    """Superclass for DataFrame and Series

    Parameters
    ----------
    dsk: dict
        The dask graph to compute this DataFrame
    name: str
        The key prefix that specifies which keys in the dask comprise this
        particular DataFrame / Series
    meta: pandas.DataFrame, pandas.Series, or pandas.Index
        An empty pandas object with names, dtypes, and indices matching the
        expected output.
    divisions: tuple of index values
        Values along which we partition our blocks on the index
    """

    def __init__(self, dsk, name, meta, divisions):
        if not isinstance(dsk, HighLevelGraph):
            dsk = HighLevelGraph.from_collections(name, dsk, dependencies=[])
        self.dask = dsk
        self._name = name
        meta = make_meta(meta)
        if not self._is_partition_type(meta):
            raise TypeError(
                f"Expected meta to specify type {type(self).__name__}, got type "
                f"{typename(type(meta))}"
            )
        self._meta = meta
        self.divisions = tuple(divisions)

        # Optionally cast object dtypes to `pyarrow` strings.
        # By default, if `pyarrow` and `pandas>=2` are installed,
        # we convert to pyarrow strings.
        # Disable for read_parquet since the reader takes care of this
        # conversion
        if pyarrow_strings_enabled() and "read-parquet" not in name:
            from dask.dataframe._pyarrow import check_pyarrow_string_supported

            check_pyarrow_string_supported()

            from dask.dataframe._pyarrow import (
                is_object_string_dataframe,
                is_object_string_dtype,
                is_object_string_index,
                is_object_string_series,
                to_pyarrow_string,
            )

            if (
                is_object_string_dataframe(meta)
                or is_object_string_series(meta)
                or is_object_string_index(meta)
            ):
                # Prior to pandas=1.4, `pd.Index` couldn't contain extension dtypes.
                # Here we don't cast objects to pyarrow strings where only the index
                # contains non-pyarrow string data.
                if not PANDAS_GE_140 and (
                    (
                        is_object_string_dataframe(meta)
                        and not any(is_object_string_dtype(d) for d in meta.dtypes)
                    )
                    or (
                        is_object_string_series(meta)
                        and not is_object_string_dtype(meta.dtype)
                    )
                    or is_object_string_index(meta)
                ):
                    return

                # this is an internal call, and if we enforce metadata,
                # it may interfere when reading csv with enforce=False
                result = self.map_partitions(
                    to_pyarrow_string, enforce_metadata=False, token="to_pyarrow_string"
                )
                self.dask = result.dask
                self._name = result._name
                self._meta = result._meta
                self.divisions = result.divisions

    def __dask_graph__(self) -> Graph:
        return self.dask

    def __dask_keys__(self) -> NestedKeys:
        return [(self._name, i) for i in range(self.npartitions)]

    def __dask_layers__(self) -> Sequence[str]:
        return (self._name,)

    def __dask_tokenize__(self):
        return self._name

    __dask_optimize__ = globalmethod(
        optimize, key="dataframe_optimize", falsey=dont_optimize
    )
    __dask_scheduler__ = staticmethod(DEFAULT_GET)

    def __dask_postcompute__(self):
        return finalize, ()

    def __dask_postpersist__(self):
        return self._rebuild, ()

    def _rebuild(self, dsk, *, rename=None):
        name = self._name
        if rename:
            name = rename.get(name, name)
        return type(self)(dsk, name, self._meta, self.divisions)

    @property
    def _constructor(self):
        return new_dd_object

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

        It is uncommon to set ``divisions`` directly. Instead, use ``set_index``,
        which sorts and splits the data as needed.
        See https://docs.dask.org/en/latest/dataframe-design.html#partitions.
        """
        return self._divisions

    @divisions.setter
    def divisions(self, value):
        if not isinstance(value, tuple):
            raise TypeError("divisions must be a tuple")

        if hasattr(self, "_divisions") and len(value) != len(self._divisions):
            n = len(self._divisions)
            raise ValueError(
                f"This dataframe has npartitions={n - 1}, divisions should be a "
                f"tuple of length={n}, got {len(value)}"
            )

        if None in value:
            if any(v is not None for v in value):
                raise ValueError(
                    "divisions may not contain a mix of None and non-None values"
                )
        else:
            # Known divisions, check monotonically increasing

            # XXX: if the index dtype is an ordered categorical dtype, then we skip the
            # sortedness check, since the order is dtype dependent
            index_dtype = getattr(self._meta, "index", self._meta).dtype
            if not (
                isinstance(index_dtype, pd.CategoricalDtype) and index_dtype.ordered
            ):
                if value != tuple(sorted(value)):
                    raise ValueError("divisions must be sorted")

        self._divisions = value

    @property
    def npartitions(self) -> int:
        """Return number of partitions"""
        return len(self.divisions) - 1

    @property
    @derived_from(pd.DataFrame)
    def attrs(self):
        return self._meta.attrs

    @attrs.setter
    def attrs(self, value):
        self._meta.attrs = dict(value)

    @property
    def size(self):
        """Size of the Series or DataFrame as a Delayed object.

        Examples
        --------
        >>> series.size  # doctest: +SKIP
        dd.Scalar<size-ag..., dtype=int64>
        """
        return self.reduction(
            methods.size, np.sum, token="size", meta=int, split_every=False
        )

    @property
    def _meta_nonempty(self):
        """A non-empty version of `_meta` with fake data."""
        return meta_nonempty(self._meta)

    @property
    def _args(self):
        return (self.dask, self._name, self._meta, self.divisions)

    def __getstate__(self):
        return self._args

    def __setstate__(self, state):
        self.dask, self._name, self._meta, self._divisions = state

    def copy(self, deep=False):
        """Make a copy of the dataframe

        This is strictly a shallow copy of the underlying computational graph.
        It does not affect the underlying data

        Parameters
        ----------
        deep : boolean, default False
            The deep value must be `False` and it is declared as a parameter just for
            compatibility with third-party libraries like cuDF
        """
        if deep is not False:
            raise ValueError(
                "The `deep` value must be False. This is strictly a shallow copy "
                "of the underlying computational graph."
            )
        return new_dd_object(self.dask, self._name, self._meta, self.divisions)

    def __array__(self, dtype=None, **kwargs):
        self._computed = self.compute()
        x = np.array(self._computed)
        return x

    def __array_wrap__(self, array, context=None):
        raise NotImplementedError

    def __array_ufunc__(self, numpy_ufunc, method, *inputs, **kwargs):
        out = kwargs.get("out", ())
        for x in inputs + out:
            # ufuncs work with 0-dimensional NumPy ndarrays
            # so we don't want to raise NotImplemented
            if isinstance(x, np.ndarray) and x.shape == ():
                continue
            elif not isinstance(
                x, (Number, Scalar, _Frame, Array, pd.DataFrame, pd.Series, pd.Index)
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

    @property
    def _elemwise(self):
        return elemwise

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

    def __repr__(self):
        data = self._repr_data().to_string(max_rows=5, show_dimensions=False)
        _str_fmt = """Dask {klass} Structure:
{data}
Dask Name: {name}, {layers}"""
        if len(self.columns) == 0:
            data = data.partition("\n")[-1].replace("Index", "Divisions")
            _str_fmt = f"Empty {_str_fmt}"
        return _str_fmt.format(
            klass=self.__class__.__name__,
            data=data,
            name=key_split(self._name),
            layers=maybe_pluralize(len(self.dask.layers), "graph layer"),
        )

    @property
    def index(self):
        """Return dask Index instance"""
        return self.map_partitions(
            getattr,
            "index",
            token=key_split(self._name) + "-index",
            meta=self._meta.index,
            enforce_metadata=False,
        )

    @index.setter
    def index(self, value):
        self.divisions = value.divisions
        result = map_partitions(
            methods.assign_index, self, value, enforce_metadata=False
        )
        self.dask = result.dask
        self._name = result._name
        self._meta = result._meta

    def reset_index(self, drop=False):
        """Reset the index to the default index.

        Note that unlike in ``pandas``, the reset ``dask.dataframe`` index will
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
        return self.map_partitions(
            M.reset_index, drop=drop, enforce_metadata=False
        ).clear_divisions()

    @property
    def known_divisions(self):
        """Whether divisions are already known"""
        return len(self.divisions) > 0 and self.divisions[0] is not None

    def clear_divisions(self):
        """Forget division information"""
        divisions = (None,) * (self.npartitions + 1)
        return type(self)(self.dask, self._name, self._meta, divisions)

    def enforce_runtime_divisions(self):
        """Enforce the current divisions at runtime"""
        if not self.known_divisions:
            raise ValueError("No known divisions to enforce!")

        def _check_divisions(df, expect):
            # Check divisions
            id, expect_min, expect_max, last = expect
            real_min = df.index.min()
            real_max = df.index.max()
            # Upper division of the last partition is often set to
            # the max value. For all other partitions, the upper
            # division should be greater than the maximum value.
            valid_min = real_min >= expect_min
            valid_max = (real_max <= expect_max) if last else (real_max < expect_max)
            if not (valid_min and valid_max):
                raise RuntimeError(
                    f"`enforce_runtime_divisions` failed for partition {id}."
                    f" Expected a range of [{expect_min}, {expect_max}), "
                    f" but the real range was [{real_min}, {real_max}]."
                )
            return df

        return self.map_partitions(
            _check_divisions,
            BlockwiseDepDict(
                {
                    (i,): (i, dmin, dmax, i == (self.npartitions - 1))
                    for i, (dmin, dmax) in enumerate(
                        zip(self.divisions[:-1], self.divisions[1:])
                    )
                }
            ),
            meta=self._meta,
            enforce_metadata=False,
        )

    def compute_current_divisions(self, col=None):
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

        >>> ddf.divisions = divisions
        >>> ddf.known_divisions
        True

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
            return self.divisions

        from dask.dataframe.shuffle import compute_divisions

        return compute_divisions(self, col=col)

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

        Examples
        --------
        >>> import dask
        >>> ddf = dask.datasets.timeseries(start="2021-01-01", end="2021-01-07", freq="1h")
        >>> ddf.get_partition(0)  # doctest: +NORMALIZE_WHITESPACE
        Dask DataFrame Structure:
                         name     id        x        y
        npartitions=1
        2021-01-01     string  int64  float64  float64
        2021-01-02        ...    ...      ...      ...
        Dask Name: get-partition, 3 graph layers

        See Also
        --------
        DataFrame.partitions
        """
        if 0 <= n < self.npartitions:
            name = f"get-partition-{str(n)}-{self._name}"
            divisions = self.divisions[n : n + 2]
            layer = {(name, 0): (self._name, n)}
            graph = HighLevelGraph.from_collections(name, layer, dependencies=[self])
            return new_dd_object(graph, name, self._meta, divisions)
        else:
            msg = f"n must be 0 <= n < {self.npartitions}"
            raise ValueError(msg)

    def _drop_duplicates_shuffle(
        self, split_out, split_every, shuffle_method, ignore_index, **kwargs
    ):
        # Private method that drops duplicate rows using a
        # shuffle-based algorithm.
        # Used by `_Frame.drop_duplicates`.

        # Make sure we have a DataFrame to shuffle
        if isinstance(self, Index):
            df = self.to_frame(name=self.name or "__index__")
        elif isinstance(self, Series):
            df = self.to_frame(name=self.name or "__series__")
        else:
            df = self

        # Choose appropriate shuffle partitioning
        split_every = 8 if split_every is None else split_every
        shuffle_npartitions = max(
            df.npartitions // (split_every or df.npartitions),
            split_out,
        )

        # Deduplicate, then shuffle, then deduplicate again
        chunk = M.drop_duplicates
        deduplicated = (
            df.map_partitions(
                chunk,
                token="drop-duplicates-chunk",
                meta=df._meta,
                ignore_index=ignore_index,
                enforce_metadata=False,
                transform_divisions=False,
                **kwargs,
            )
            .shuffle(
                kwargs.get("subset", None) or list(df.columns),
                ignore_index=ignore_index,
                npartitions=shuffle_npartitions,
                shuffle_method=shuffle_method,
            )
            .map_partitions(
                chunk,
                meta=df._meta,
                ignore_index=ignore_index,
                token="drop-duplicates-agg",
                transform_divisions=False,
                **kwargs,
            )
        )

        # Convert back to Series/Index if necessary
        if isinstance(self, Index):
            deduplicated = deduplicated.set_index(
                self.name or "__index__", sort=False
            ).index
            if deduplicated.name == "__index__":
                deduplicated.name = None
        elif isinstance(self, Series):
            deduplicated = deduplicated[self.name or "__series__"]
            if deduplicated.name == "__series__":
                deduplicated.name = None

        # Return `split_out` partitions
        return deduplicated.repartition(npartitions=split_out)

    @_deprecated_kwarg("shuffle", "shuffle_method")
    @derived_from(
        pd.DataFrame,
        inconsistencies="keep=False will raise a ``NotImplementedError``",
    )
    def drop_duplicates(
        self,
        subset=None,
        split_every=None,
        split_out=1,
        shuffle_method=None,
        ignore_index=False,
        **kwargs,
    ):
        if subset is not None:
            # Let pandas error on bad inputs
            self._meta_nonempty.drop_duplicates(subset=subset, **kwargs)
            kwargs["subset"] = subset
            split_out_setup = split_out_on_cols
            split_out_setup_kwargs = {"cols": subset}
        else:
            self._meta_nonempty.drop_duplicates(**kwargs)
            split_out_setup = split_out_setup_kwargs = None

        if kwargs.get("keep", True) is False:
            raise NotImplementedError("drop_duplicates with keep=False")

        # Check if we should use a shuffle-based algorithm,
        # which is typically faster when we are not reducing
        # to a small number of partitions
        shuffle_method = _determine_split_out_shuffle(shuffle_method, split_out)
        if shuffle_method:
            return self._drop_duplicates_shuffle(
                split_out,
                split_every,
                shuffle_method,
                ignore_index,
                **kwargs,
            )

        # Use general ACA reduction
        # (Usually best when split_out == 1)
        chunk = M.drop_duplicates
        return aca(
            self,
            chunk=chunk,
            aggregate=chunk,
            meta=self._meta,
            token="drop-duplicates",
            split_every=split_every,
            split_out=split_out,
            split_out_setup=split_out_setup,
            split_out_setup_kwargs=split_out_setup_kwargs,
            ignore_index=ignore_index,
            **kwargs,
        )

    def __len__(self):
        return self.reduction(
            len, np.sum, token="len", meta=int, split_every=False
        ).compute()

    def __bool__(self):
        raise ValueError(
            f"The truth value of a {self.__class__.__name__} is ambiguous. "
            "Use a.any() or a.all()."
        )

    __nonzero__ = __bool__  # python 2

    def _scalarfunc(self, cast_type):
        def wrapper():
            raise TypeError(f"cannot convert the series to {cast_type}")

        return wrapper

    def __float__(self):
        return self._scalarfunc(float)

    def __int__(self):
        return self._scalarfunc(int)

    __long__ = __int__  # python 2

    def __complex__(self):
        return self._scalarfunc(complex)

    @insert_meta_param_description(pad=12)
    def map_partitions(self, func, *args, **kwargs):
        """Apply Python function on each DataFrame partition.

        Note that the index and divisions are assumed to remain unchanged.

        Parameters
        ----------
        func : function
            The function applied to each partition. If this function accepts
            the special ``partition_info`` keyword argument, it will receive
            information on the partition's relative location within the
            dataframe.
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
        to clear them afterwards:

        >>> ddf.map_partitions(func).clear_divisions()  # doctest: +SKIP

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
        return map_partitions(func, self, *args, **kwargs)

    @insert_meta_param_description(pad=12)
    def map_overlap(self, func, before, after, *args, **kwargs):
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
        from dask.dataframe.rolling import map_overlap

        return map_overlap(func, self, before, after, *args, **kwargs)

    def memory_usage_per_partition(self, index=True, deep=False):
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
        return self.map_partitions(
            total_mem_usage, index=index, deep=deep
        ).clear_divisions()

    @insert_meta_param_description(pad=12)
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
        25

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
        if aggregate is None:
            aggregate = chunk

        if combine is None:
            if combine_kwargs:
                raise ValueError("`combine_kwargs` provided with no `combine`")
            combine = aggregate
            combine_kwargs = aggregate_kwargs

        chunk_kwargs = chunk_kwargs.copy() if chunk_kwargs else {}
        chunk_kwargs["aca_chunk"] = chunk

        combine_kwargs = combine_kwargs.copy() if combine_kwargs else {}
        combine_kwargs["aca_combine"] = combine

        aggregate_kwargs = aggregate_kwargs.copy() if aggregate_kwargs else {}
        aggregate_kwargs["aca_aggregate"] = aggregate

        return aca(
            self,
            chunk=_reduction_chunk,
            aggregate=_reduction_aggregate,
            combine=_reduction_combine,
            meta=meta,
            token=token,
            split_every=split_every,
            chunk_kwargs=chunk_kwargs,
            aggregate_kwargs=aggregate_kwargs,
            combine_kwargs=combine_kwargs,
            **kwargs,
        )

    @derived_from(pd.DataFrame)
    def pipe(self, func, *args, **kwargs):
        # Taken from pandas:
        # https://github.com/pydata/pandas/blob/master/pandas/core/generic.py#L2698-L2707
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

    def random_split(self, frac, random_state=None, shuffle=False):
        """Pseudorandomly split dataframe into different pieces row-wise

        Parameters
        ----------
        frac : list
            List of floats that should sum to one.
        random_state : int or np.random.RandomState
            If int create a new RandomState with this as the seed.
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
        state_data = random_state_data(self.npartitions, random_state)
        token = tokenize(self, frac, random_state)
        name = "split-" + token
        layer = {
            (name, i): (pd_split, (self._name, i), frac, state, shuffle)
            for i, state in enumerate(state_data)
        }

        out = []
        for i in range(len(frac)):
            name2 = "split-%d-%s" % (i, token)
            dsk2 = {
                (name2, j): (getitem, (name, j), i) for j in range(self.npartitions)
            }
            graph = HighLevelGraph.from_collections(
                name2, merge(dsk2, layer), dependencies=[self]
            )
            out_df = type(self)(graph, name2, self._meta, self.divisions)
            out.append(out_df)
        return out

    def head(self, n=5, npartitions=1, compute=True):
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
        if npartitions <= -1:
            npartitions = self.npartitions
        # No need to warn if we're already looking at all partitions
        safe = npartitions != self.npartitions
        return self._head(n=n, npartitions=npartitions, compute=compute, safe=safe)

    def _head(self, n, npartitions, compute, safe):
        if npartitions <= -1:
            npartitions = self.npartitions
        if npartitions > self.npartitions:
            raise ValueError(
                f"only {self.npartitions} partitions, head received {npartitions}"
            )

        name = f"head-{npartitions}-{n}-{self._name}"
        if safe:
            head = safe_head
        else:
            head = M.head

        if npartitions > 1:
            name_p = f"head-partial-{n}-{self._name}"

            dsk = {}
            for i in range(npartitions):
                dsk[(name_p, i)] = (M.head, (self._name, i), n)

            concat = (_concat, [(name_p, i) for i in range(npartitions)])
            dsk[(name, 0)] = (head, concat, n)
        else:
            dsk = {(name, 0): (head, (self._name, 0), n)}

        graph = HighLevelGraph.from_collections(name, dsk, dependencies=[self])
        result = new_dd_object(
            graph, name, self._meta, [self.divisions[0], self.divisions[npartitions]]
        )

        if compute:
            result = result.compute()
        return result

    def tail(self, n=5, compute=True):
        """Last n rows of the dataset

        Caveat, the only checks the last n rows of the last partition.
        """
        name = "tail-%d-%s" % (n, self._name)
        dsk = {(name, 0): (M.tail, (self._name, self.npartitions - 1), n)}

        graph = HighLevelGraph.from_collections(name, dsk, dependencies=[self])
        result = new_dd_object(graph, name, self._meta, self.divisions[-2:])

        if compute:
            result = result.compute()
        return result

    @property
    def loc(self):
        """Purely label-location based indexer for selection by label.

        >>> df.loc["b"]  # doctest: +SKIP
        >>> df.loc["b":"d"]  # doctest: +SKIP
        """
        from dask.dataframe.indexing import _LocIndexer

        return _LocIndexer(self)

    def _partitions(self, index):
        if not isinstance(index, tuple):
            index = (index,)
        from dask.array.slicing import normalize_index

        index = normalize_index(index, (self.npartitions,))
        index = tuple(slice(k, k + 1) if isinstance(k, Number) else k for k in index)
        name = "blocks-" + tokenize(self, index)
        new_keys = np.array(self.__dask_keys__(), dtype=object)[index].tolist()

        divisions = [self.divisions[i] for _, i in new_keys] + [
            self.divisions[new_keys[-1][1] + 1]
        ]
        dsk = {(name, i): tuple(key) for i, key in enumerate(new_keys)}

        graph = HighLevelGraph.from_collections(name, dsk, dependencies=[self])
        return new_dd_object(graph, name, self._meta, divisions)

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

    # Note: iloc is implemented only on DataFrame

    def repartition(
        self,
        divisions=None,
        npartitions=None,
        partition_size=None,
        freq=None,
        force=False,
    ):
        """Repartition dataframe along new divisions

        Parameters
        ----------
        divisions : list, optional
            The "dividing lines" used to split the dataframe into partitions.
            For ``divisions=[0, 10, 50, 100]``, there would be three output partitions,
            where the new index contained [0, 10), [10, 50), and [50, 100), respectively.
            See https://docs.dask.org/en/latest/dataframe-design.html#partitions.
            Only used if npartitions and partition_size isn't specified.
            For convenience if given an integer this will defer to npartitions
            and if given a string it will defer to partition_size (see below)
        npartitions : int, optional
            Approximate number of partitions of output. Only used if partition_size
            isn't specified. The number of partitions used may be slightly
            lower than npartitions depending on data distribution, but will never be
            higher.
        partition_size: int or string, optional
            Max number of bytes of memory for each partition. Use numbers or
            strings like 5MB. If specified npartitions and divisions will be
            ignored. Note that the size reflects the number of bytes used as
            computed by ``pandas.DataFrame.memory_usage``, which will not
            necessarily match the size when storing to disk.

            .. warning::

               This keyword argument triggers computation to determine
               the memory size of each partition, which may be expensive.

        freq : str, pd.Timedelta
            A period on which to partition timeseries data like ``'7D'`` or
            ``'12h'`` or ``pd.Timedelta(hours=12)``.  Assumes a datetime index.
        force : bool, default False
            Allows the expansion of the existing divisions.
            If False then the new divisions' lower and upper bounds must be
            the same as the old divisions'.

        Notes
        -----
        Exactly one of `divisions`, `npartitions`, `partition_size`, or `freq`
        should be specified. A ``ValueError`` will be raised when that is
        not the case.

        Also note that ``len(divisions)`` is equal to ``npartitions + 1``. This is because ``divisions``
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
        if isinstance(divisions, int):
            warnings.warn(
                "divisions is an integer and will be inferred as npartitions instead. "
                "This automatic inference is deprecated and will change in the future. "
                f"Please set npartitions={divisions} instead.",
                FutureWarning,
                stacklevel=2,
            )
            npartitions = divisions
            divisions = None
        if isinstance(divisions, str):
            warnings.warn(
                "divisions is a string and will be inferred as partition_size instead. "
                "This automatic inference is deprecated and will change in the future. "
                f"Please set partition_size={divisions} instead.",
                FutureWarning,
                stacklevel=2,
            )
            partition_size = divisions
            divisions = None
        if (
            sum(
                [
                    partition_size is not None,
                    divisions is not None,
                    npartitions is not None,
                    freq is not None,
                ]
            )
            != 1
        ):
            raise ValueError(
                "Please provide exactly one of ``npartitions=``, ``freq=``, "
                "``divisions=``, ``partition_size=`` keyword arguments"
            )

        if partition_size is not None:
            return repartition_size(self, partition_size)
        elif npartitions is not None:
            return repartition_npartitions(self, npartitions)
        elif divisions is not None:
            return repartition(self, divisions, force=force)
        elif freq is not None:
            return repartition_freq(self, freq=freq)

    @_deprecated_kwarg("shuffle", "shuffle_method")
    def shuffle(
        self,
        on,
        npartitions=None,
        max_branch=None,
        shuffle_method=None,
        ignore_index=False,
        compute=None,
    ):
        """Rearrange DataFrame into new partitions

        Uses hashing of `on` to map rows to output partitions. After this
        operation, rows with the same value of `on` will be in the same
        partition.

        Parameters
        ----------
        on : str, list of str, or Series, Index, or DataFrame
            Column(s) or index to be used to map rows to output partitions
        npartitions : int, optional
            Number of partitions of output. Partition count will not be
            changed by default.
        max_branch: int, optional
            The maximum number of splits per input partition. Used within
            the staged shuffling algorithm.
        shuffle_method: {'disk', 'tasks', 'p2p'}, optional
            Either ``'disk'`` for single-node operation or ``'tasks'`` and
            ``'p2p'`` for distributed operation.  Will be inferred by your
            current scheduler.
        ignore_index: bool, default False
            Ignore index during shuffle.  If ``True``, performance may improve,
            but index values will not be preserved.
        compute: bool
            Whether or not to trigger an immediate computation. Defaults to False.

        Notes
        -----
        This does not preserve a meaningful index/partitioning scheme. This
        is not deterministic if done in parallel.

        Examples
        --------
        >>> df = df.shuffle(df.columns[0])  # doctest: +SKIP
        """
        from dask.dataframe.shuffle import shuffle as dd_shuffle

        return dd_shuffle(
            self,
            on,
            npartitions=npartitions,
            max_branch=max_branch,
            shuffle_method=shuffle_method,
            ignore_index=ignore_index,
            compute=compute,
        )

    def _limit_fillna(self, method=None, *, limit=None, skip_check=None, meta=None):
        if limit is None:
            name = "fillna-chunk-" + tokenize(self, method)
            dsk = {
                (name, i): (
                    methods.fillna_check,
                    (self._name, i),
                    method,
                    i != skip_check,
                )
                for i in range(self.npartitions)
            }
            graph = HighLevelGraph.from_collections(name, dsk, dependencies=[self])
            return new_dd_object(graph, name, meta, self.divisions)
        else:
            return self

    @_deprecated_kwarg("method", None, comment="Use ffill or bfill instead.")
    @derived_from(pd.DataFrame)
    def fillna(self, value=None, method=None, limit=None, axis=None):
        if method is None and limit is not None:
            raise NotImplementedError("fillna with set limit and method=None")

        axis = self._validate_axis(axis)
        test_value = (
            value._meta_nonempty if isinstance(value, (_Frame, Scalar)) else value
        )

        # let it raise a FutureWarning if `method` is not None
        meta = self._meta_nonempty.fillna(
            value=test_value, method=method, limit=limit, axis=axis
        )

        if method is None:
            # Control whether or not dask's partition alignment happens.
            # We don't want for a pandas Series.
            # We do want it for a dask Series
            if is_series_like(value) and not is_dask_collection(value):
                args = ()
                kwargs = {"value": value}
            else:
                args = (value,)
                kwargs = {}

            return self.map_partitions(
                M.fillna,
                *args,
                limit=limit,
                axis=axis,
                meta=meta,
                enforce_metadata=False,
                **kwargs,
            )
        elif method in ("pad", "ffill"):
            return self.ffill(limit=limit, axis=axis)
        else:
            return self.bfill(limit=limit, axis=axis)

    @derived_from(pd.DataFrame)
    def ffill(self, axis=None, limit=None):
        axis = self._validate_axis(axis)
        meta = self._meta_nonempty.ffill(limit=limit, axis=axis)

        if axis == 1:
            return self.map_partitions(
                M.ffill, limit=limit, axis=axis, meta=meta, enforce_metadata=False
            )

        before, after = 1 if limit is None else limit, 0
        parts = self._limit_fillna("ffill", limit=limit, skip_check=0, meta=meta)
        return parts.map_overlap(M.ffill, before, after, limit=limit, meta=meta)

    @derived_from(pd.DataFrame)
    def bfill(self, axis=None, limit=None):
        axis = self._validate_axis(axis)
        meta = self._meta_nonempty.bfill(limit=limit, axis=axis)

        if axis == 1:
            return self.map_partitions(
                M.bfill, limit=limit, axis=axis, meta=meta, enforce_metadata=False
            )

        before, after = 0, 1 if limit is None else limit
        parts = self._limit_fillna(
            "bfill", limit=limit, skip_check=self.npartitions - 1, meta=meta
        )
        return parts.map_overlap(M.bfill, before, after, limit=limit, meta=meta)

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

        name = "sample-" + tokenize(self, frac, replace, random_state)

        state_data = random_state_data(self.npartitions, random_state)
        dsk = {
            (name, i): (methods.sample, (self._name, i), state, frac, replace)
            for i, state in enumerate(state_data)
        }

        graph = HighLevelGraph.from_collections(name, dsk, dependencies=[self])
        return new_dd_object(graph, name, self._meta, self.divisions)

    @derived_from(pd.DataFrame)
    def replace(self, to_replace=None, value=None, regex=False):
        # In PANDAS_GE_140 pandas starts using no_default instead of None
        value_kwarg = {"value": value} if value is not None else {}
        return self.map_partitions(
            M.replace,
            to_replace=to_replace,
            **value_kwarg,
            regex=regex,
            enforce_metadata=False,
        )

    def to_dask_array(self, lengths=None, meta=None):
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

        Returns
        -------
        """
        if lengths is True:
            lengths = tuple(self.map_partitions(len, enforce_metadata=False).compute())

        arr = self.values

        chunks = self._validate_chunks(arr, lengths)
        arr._chunks = chunks

        if meta is not None:
            arr._meta = meta

        return arr

    def to_hdf(self, path_or_buf, key, mode="a", append=False, **kwargs):
        """See dd.to_hdf docstring for more information"""
        from dask.dataframe.io import to_hdf

        return to_hdf(self, path_or_buf, key, mode, append, **kwargs)

    def to_csv(self, filename, **kwargs):
        """See dd.to_csv docstring for more information"""
        from dask.dataframe.io import to_csv

        return to_csv(self, filename, **kwargs)

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
        """See dd.to_sql docstring for more information"""
        from dask.dataframe.io import to_sql

        return to_sql(
            self,
            name=name,
            uri=uri,
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

    def to_json(self, filename, *args, **kwargs):
        """See dd.to_json docstring for more information"""
        from dask.dataframe.io import to_json

        return to_json(self, filename, *args, **kwargs)

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
        dask.dataframe.from_delayed
        """
        keys = self.__dask_keys__()
        graph = self.__dask_graph__()
        layer = self.__dask_layers__()[0]
        if optimize_graph:
            graph = self.__dask_optimize__(graph, self.__dask_keys__())
            layer = "delayed-" + self._name
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
        from dask.dataframe.io import to_backend

        return to_backend(self, backend=backend, **kwargs)

    @classmethod
    def _get_unary_operator(cls, op):
        return lambda self: elemwise(op, self)

    @classmethod
    def _get_binary_operator(cls, op, inv=False):
        if inv:
            return lambda self, other: elemwise(op, other, self)
        else:
            return lambda self, other: elemwise(op, self, other)

    @_deprecated_kwarg("axis", None)
    def rolling(
        self, window, min_periods=None, center=False, win_type=None, axis=no_default
    ):
        """Provides rolling transformations.

        Parameters
        ----------
        window : int, str, offset
           Size of the moving window. This is the number of observations used
           for calculating the statistic. When not using a ``DatetimeIndex``,
           the window size must not be so large as to span more than one
           adjacent partition. If using an offset or offset alias like '5D',
           the data must have a ``DatetimeIndex``

           .. versionchanged:: 0.15.0

              Now accepts offsets and string offset aliases

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
        from dask.dataframe.rolling import Rolling

        if isinstance(window, Integral):
            if window < 0:
                raise ValueError("window must be >= 0")

        if min_periods is not None:
            if not isinstance(min_periods, Integral):
                raise ValueError("min_periods must be an integer")
            if min_periods < 0:
                raise ValueError("min_periods must be >= 0")

        return Rolling(
            self,
            window=window,
            min_periods=min_periods,
            center=center,
            win_type=win_type,
            axis=axis,
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
        axis = self._validate_axis(axis)
        if not isinstance(periods, Integral):
            raise TypeError("periods must be an integer")

        if axis == 1:
            return self.map_partitions(
                M.diff, token="diff", periods=periods, axis=1, enforce_metadata=False
            )

        before, after = (periods, 0) if periods > 0 else (0, -periods)
        return self.map_overlap(M.diff, before, after, token="diff", periods=periods)

    @derived_from(pd.DataFrame)
    def shift(self, periods=1, freq=None, axis=0):
        axis = self._validate_axis(axis)
        if not isinstance(periods, Integral):
            raise TypeError("periods must be an integer")

        if axis == 1:
            return self.map_partitions(
                M.shift,
                token="shift",
                periods=periods,
                freq=freq,
                axis=1,
                enforce_metadata=False,
            )

        if freq is None:
            before, after = (periods, 0) if periods > 0 else (0, -periods)
            return self.map_overlap(
                M.shift, before, after, token="shift", periods=periods
            )

        # Let pandas error on invalid arguments
        meta = self._meta_nonempty.shift(periods, freq=freq)
        out = self.map_partitions(
            M.shift,
            token="shift",
            periods=periods,
            freq=freq,
            meta=meta,
            enforce_metadata=False,
            transform_divisions=False,
        )
        return maybe_shift_divisions(out, periods, freq=freq)

    def _reduction_agg(
        self,
        name,
        axis=None,
        skipna=True,
        split_every=False,
        out=None,
        numeric_only=None,
        none_is_zero=True,
    ):
        axis = self._validate_axis(axis, none_is_zero=none_is_zero)

        if has_keyword(getattr(self._meta_nonempty, name), "numeric_only"):
            numeric_only_kwargs = {"numeric_only": numeric_only}
        else:
            numeric_only_kwargs = {}

        with check_numeric_only_deprecation(name, True):
            meta = getattr(self._meta_nonempty, name)(
                axis=axis, skipna=skipna, **numeric_only_kwargs
            )

        token = self._token_prefix + name
        if axis == 1:
            result = self.map_partitions(
                _getattr_numeric_only,
                meta=meta,
                token=token,
                skipna=skipna,
                axis=axis,
                _dask_method_name=name,
                **numeric_only_kwargs,
            )
            return handle_out(out, result)
        else:
            result = self.reduction(
                _getattr_numeric_only,
                meta=meta,
                token=token,
                skipna=skipna,
                axis=axis,
                split_every=split_every,
                _dask_method_name=name,
                **numeric_only_kwargs,
            )
            if isinstance(self, DataFrame) and isinstance(result, Series):
                result.divisions = (self.columns.min(), self.columns.max())
            return handle_out(out, result)

    @derived_from(pd.DataFrame)
    def add_prefix(self, prefix):
        res = self.map_partitions(M.add_prefix, prefix)
        if self.known_divisions and is_series_like(self):
            res.divisions = tuple(prefix + str(division) for division in self.divisions)
        return res

    @derived_from(pd.DataFrame)
    def add_suffix(self, suffix):
        res = self.map_partitions(M.add_suffix, suffix)
        if self.known_divisions and is_series_like(self):
            res.divisions = tuple(str(division) + suffix for division in self.divisions)
        return res

    @derived_from(pd.DataFrame)
    def abs(self):
        _raise_if_object_series(self, "abs")
        meta = self._meta_nonempty.abs()
        return self.map_partitions(M.abs, meta=meta, enforce_metadata=False)

    @_dummy_numpy_dispatcher("out", deprecated=True)
    @derived_from(pd.DataFrame)
    def all(self, axis=None, skipna=True, split_every=False, out=None):
        return self._reduction_agg(
            "all", axis=axis, skipna=skipna, split_every=split_every, out=out
        )

    @_dummy_numpy_dispatcher("out", deprecated=True)
    @derived_from(pd.DataFrame)
    def any(self, axis=None, skipna=True, split_every=False, out=None):
        return self._reduction_agg(
            "any", axis=axis, skipna=skipna, split_every=split_every, out=out
        )

    @_dummy_numpy_dispatcher("dtype", "out", deprecated=True)
    @derived_from(pd.DataFrame)
    def sum(
        self,
        axis=None,
        skipna=True,
        split_every=False,
        dtype=None,
        out=None,
        min_count=None,
        numeric_only=None,
    ):
        result = self._reduction_agg(
            "sum",
            axis=axis,
            skipna=skipna,
            split_every=split_every,
            out=out,
            numeric_only=numeric_only,
        )
        if min_count:
            cond = self.notnull().sum(axis=axis) >= min_count
            if is_series_like(cond):
                return result.where(cond, other=np.nan)
            else:
                return _scalar_binary(
                    lambda x, y: result if x is y else np.nan, cond, True
                )
        else:
            return result

    @_dummy_numpy_dispatcher("dtype", "out", deprecated=True)
    @derived_from(pd.DataFrame)
    def prod(
        self,
        axis=None,
        skipna=True,
        split_every=False,
        dtype=None,
        out=None,
        min_count=None,
        numeric_only=None,
    ):
        result = self._reduction_agg(
            "prod",
            axis=axis,
            skipna=skipna,
            split_every=split_every,
            out=out,
            numeric_only=numeric_only,
        )
        if min_count:
            cond = self.notnull().sum(axis=axis) >= min_count
            if is_series_like(cond):
                return result.where(cond, other=np.nan)
            else:
                return _scalar_binary(
                    lambda x, y: result if x is y else np.nan, cond, True
                )
        else:
            return result

    product = prod  # aliased dd.product

    @_dummy_numpy_dispatcher("out", deprecated=True)
    @derived_from(pd.DataFrame)
    def max(self, axis=0, skipna=True, split_every=False, out=None, numeric_only=None):
        if (
            PANDAS_GE_140
            and not PANDAS_GE_200
            and axis is None
            and isinstance(self, DataFrame)
        ):
            warnings.warn(
                "In a future version, DataFrame.max(axis=None) will return a scalar max over the entire DataFrame. "
                "To retain the old behavior, use 'frame.max(axis=0)' or just 'frame.max()'",
                FutureWarning,
            )
            axis = 0

        return self._reduction_agg(
            "max",
            axis=axis,
            skipna=skipna,
            split_every=split_every,
            out=out,
            # Starting in pandas 2.0, `axis=None` does a full aggregation across both axes
            none_is_zero=not PANDAS_GE_200,
            numeric_only=numeric_only,
        )

    @_dummy_numpy_dispatcher("out", deprecated=True)
    @derived_from(pd.DataFrame)
    def min(self, axis=0, skipna=True, split_every=False, out=None, numeric_only=None):
        if (
            PANDAS_GE_140
            and not PANDAS_GE_200
            and axis is None
            and isinstance(self, DataFrame)
        ):
            warnings.warn(
                "In a future version, DataFrame.min(axis=None) will return a scalar min over the entire DataFrame. "
                "To retain the old behavior, use 'frame.min(axis=0)' or just 'frame.min()'",
                FutureWarning,
            )
            axis = 0

        return self._reduction_agg(
            "min",
            axis=axis,
            skipna=skipna,
            split_every=split_every,
            out=out,
            # Starting in pandas 2.0, `axis=None` does a full aggregation across both axes
            none_is_zero=not PANDAS_GE_200,
            numeric_only=numeric_only,
        )

    @derived_from(pd.DataFrame)
    def idxmax(
        self, axis=None, skipna=True, split_every=False, numeric_only=no_default
    ):
        fn = "idxmax"
        axis = self._validate_axis(axis)
        numeric_only_kwargs = check_numeric_only_valid(numeric_only, "idxmax")
        meta = self._meta_nonempty.idxmax(
            axis=axis, skipna=skipna, **numeric_only_kwargs
        )
        if axis == 1:
            return map_partitions(
                M.idxmax,
                self,
                meta=meta,
                token=self._token_prefix + fn,
                skipna=skipna,
                axis=axis,
                enforce_metadata=False,
                **numeric_only_kwargs,
            )
        else:
            scalar = not is_series_like(meta)
            result = aca(
                [self],
                chunk=idxmaxmin_chunk,
                aggregate=idxmaxmin_agg,
                combine=idxmaxmin_combine,
                meta=meta,
                aggregate_kwargs={"scalar": scalar},
                token=self._token_prefix + fn,
                split_every=split_every,
                skipna=skipna,
                fn=fn,
                **numeric_only_kwargs,
            )
            if isinstance(self, DataFrame):
                result.divisions = (min(self.columns), max(self.columns))
            return result

    @derived_from(pd.DataFrame)
    def idxmin(
        self, axis=None, skipna=True, split_every=False, numeric_only=no_default
    ):
        fn = "idxmin"
        axis = self._validate_axis(axis)
        numeric_only_kwargs = check_numeric_only_valid(numeric_only, "idxmax")
        meta = self._meta_nonempty.idxmax(
            axis=axis, skipna=skipna, **numeric_only_kwargs
        )
        if axis == 1:
            return map_partitions(
                M.idxmin,
                self,
                meta=meta,
                token=self._token_prefix + fn,
                skipna=skipna,
                axis=axis,
                enforce_metadata=False,
                **numeric_only_kwargs,
            )
        else:
            scalar = not is_series_like(meta)
            result = aca(
                [self],
                chunk=idxmaxmin_chunk,
                aggregate=idxmaxmin_agg,
                combine=idxmaxmin_combine,
                meta=meta,
                aggregate_kwargs={"scalar": scalar},
                token=self._token_prefix + fn,
                split_every=split_every,
                skipna=skipna,
                fn=fn,
                **numeric_only_kwargs,
            )
            if isinstance(self, DataFrame):
                result.divisions = (min(self.columns), max(self.columns))
            return result

    @derived_from(pd.DataFrame)
    def count(self, axis=None, split_every=False, numeric_only=False):
        # This method is shared by DataFrame / Series, but only DataFrame
        # supports `numeric_only=`. Handle accordingly here.
        numeric_only_kwargs = {}
        if is_dataframe_like(self):
            numeric_only_kwargs = get_numeric_only_kwargs(numeric_only)
        axis = self._validate_axis(axis)
        token = self._token_prefix + "count"
        if axis == 1:
            meta = self._meta_nonempty.count(axis=axis, **numeric_only_kwargs)
            return self.map_partitions(
                M.count,
                meta=meta,
                token=token,
                axis=axis,
                enforce_metadata=False,
                **numeric_only_kwargs,
            )
        else:
            meta = self._meta_nonempty.count()

            # Need the astype(int) for empty dataframes, which sum to float dtype
            result = self.reduction(
                M.count,
                aggregate=_count_aggregate,
                meta=meta,
                token=token,
                split_every=split_every,
                chunk_kwargs=numeric_only_kwargs,
            )
            if isinstance(self, DataFrame):
                result.divisions = (self.columns.min(), self.columns.max())
            return result

    @derived_from(pd.DataFrame)
    def mode(self, dropna=True, split_every=False):
        mode_series = self.reduction(
            chunk=M.value_counts,
            combine=M.sum,
            aggregate=_mode_aggregate,
            split_every=split_every,
            chunk_kwargs={"dropna": dropna},
            aggregate_kwargs={"dropna": dropna},
        )
        mode_series.name = self.name
        return mode_series

    @_dummy_numpy_dispatcher("dtype", "out", deprecated=True)
    @_numeric_only
    @derived_from(pd.DataFrame)
    def mean(
        self,
        axis=0,
        skipna=True,
        split_every=False,
        dtype=None,
        out=None,
        numeric_only=None,
    ):
        if (
            PANDAS_GE_140
            and not PANDAS_GE_200
            and axis is None
            and isinstance(self, DataFrame)
        ):
            warnings.warn(
                "In a future version, DataFrame.mean(axis=None) will return a scalar mean over the entire DataFrame. "
                "To retain the old behavior, use 'frame.mean(axis=0)' or just 'frame.mean()'",
                FutureWarning,
            )
        axis = self._validate_axis(axis, none_is_zero=not PANDAS_GE_200)
        _raise_if_object_series(self, "mean")
        # NOTE: Do we want to warn here?
        with check_numeric_only_deprecation(), check_nuisance_columns_warning():
            meta = self._meta_nonempty.mean(
                axis=axis, skipna=skipna, numeric_only=numeric_only
            )
        if axis == 1:
            result = map_partitions(
                M.mean,
                self,
                meta=meta,
                token=self._token_prefix + "mean",
                axis=axis,
                skipna=skipna,
                enforce_metadata=False,
                numeric_only=numeric_only,
            )
            return handle_out(out, result)
        else:
            num = self._get_numeric_data()
            s = num.sum(skipna=skipna, split_every=split_every)
            n = num.count(split_every=split_every)
            # Starting in pandas 2.0, `axis=None` does a full aggregation across both axes
            if PANDAS_GE_200 and axis is None and isinstance(self, DataFrame):
                result = s.sum() / n.sum()
            else:
                name = self._token_prefix + "mean-%s" % tokenize(self, axis, skipna)
                result = map_partitions(
                    methods.mean_aggregate,
                    s,
                    n,
                    token=name,
                    meta=meta,
                    enforce_metadata=False,
                    parent_meta=self._meta,
                )
                if isinstance(self, DataFrame):
                    result.divisions = (self.columns.min(), self.columns.max())
            return handle_out(out, result)

    def median_approximate(
        self,
        axis=None,
        method="default",
    ):
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
        return self.quantile(q=0.5, axis=axis, method=method).rename(None)

    @derived_from(pd.DataFrame)
    def median(self, axis=None, method="default"):
        if axis in (1, "columns") or self.npartitions == 1:
            # Can provide an exact median in these cases
            return self.median_approximate(axis=axis, method=method)
        raise NotImplementedError(
            "Dask doesn't implement an exact median in all cases as this is hard to do in parallel. "
            "See the `median_approximate` method instead, which uses an approximate algorithm."
        )

    @_dummy_numpy_dispatcher("dtype", "out", deprecated=True)
    @derived_from(pd.DataFrame)
    def var(
        self,
        axis=None,
        skipna=True,
        ddof=1,
        split_every=False,
        dtype=None,
        out=None,
        numeric_only=no_default,
    ):
        axis = self._validate_axis(axis)
        _raise_if_object_series(self, "var")
        numeric_only_kwargs = get_numeric_only_kwargs(numeric_only)
        with check_numeric_only_deprecation(
            "var", True
        ), check_reductions_runtime_warning():
            meta = self._meta_nonempty.var(
                axis=axis,
                skipna=skipna,
                **numeric_only_kwargs,
            )
        if axis == 1:
            result = map_partitions(
                M.var,
                self,
                meta=meta,
                token=self._token_prefix + "var",
                axis=axis,
                skipna=skipna,
                ddof=ddof,
                enforce_metadata=False,
                **numeric_only_kwargs,
            )
            return handle_out(out, result)
        else:
            if self.ndim == 1:
                result = self._var_1d(self, skipna, ddof, split_every)
                return handle_out(out, result)

            # pandas 1.0+ does not implement var on timedelta
            result = self._var_numeric(skipna, ddof, split_every)

            if isinstance(self, DataFrame):
                result.divisions = (self.columns.min(), self.columns.max())
            return handle_out(out, result)

    def _var_numeric(self, skipna=True, ddof=1, split_every=False):
        num = self.select_dtypes(include=["number", "bool"], exclude=[np.timedelta64])

        values_dtype = num.values.dtype
        array_values = num.values

        if not np.issubdtype(values_dtype, np.number):
            array_values = num.values.astype("f8")

        var = da.nanvar if skipna or skipna is None else da.var
        array_var = var(array_values, axis=0, ddof=ddof, split_every=split_every)

        name = self._token_prefix + "var-numeric" + tokenize(num, split_every)
        cols = num._meta.columns if is_dataframe_like(num) else None

        with check_reductions_runtime_warning():
            meta_computation = num._meta_nonempty.var(axis=0)
        var_shape = meta_computation.shape
        array_var_name = (array_var._name,) + (0,) * len(var_shape)

        layer = {(name, 0): (methods.wrap_var_reduction, array_var_name, cols)}
        graph = HighLevelGraph.from_collections(name, layer, dependencies=[array_var])

        return new_dd_object(graph, name, meta_computation, divisions=[None, None])

    def _var_1d(self, column, skipna=True, ddof=1, split_every=False):
        is_timedelta = is_timedelta64_dtype(column._meta)

        if is_timedelta:
            if not skipna:
                is_nan = column.isna()
                column = column.astype("i8")
                column = column.mask(is_nan)
            else:
                column = column.dropna().astype("i8")

        if pd.api.types.is_extension_array_dtype(column._meta_nonempty):
            # Don't have to worry about non-numeric, this raises earlier
            column = column.astype("f8")

        elif not np.issubdtype(column.dtype, np.number):
            column = column.astype("f8")

        name = self._token_prefix + "var-1d-" + tokenize(column, split_every)

        var = da.nanvar if skipna or skipna is None else da.var
        array_var = var(column.values, axis=0, ddof=ddof, split_every=split_every)

        layer = {(name, 0): (methods.wrap_var_reduction, (array_var._name,), None)}
        graph = HighLevelGraph.from_collections(name, layer, dependencies=[array_var])

        return new_dd_object(
            graph, name, column._meta_nonempty.var(), divisions=[None, None]
        )

    @_dummy_numpy_dispatcher("dtype", "out", deprecated=True)
    @_numeric_data
    @derived_from(pd.DataFrame)
    def std(
        self,
        axis=None,
        skipna=True,
        ddof=1,
        split_every=False,
        dtype=None,
        out=None,
        numeric_only=no_default,
    ):
        axis = self._validate_axis(axis)
        _raise_if_object_series(self, "std")
        _raise_if_not_series_or_dataframe(self, "std")
        numeric_kwargs = get_numeric_only_kwargs(numeric_only)

        with check_numeric_only_deprecation(
            "std", True
        ), check_reductions_runtime_warning():
            meta = self._meta_nonempty.std(axis=axis, skipna=skipna, **numeric_kwargs)
        is_df_like = is_dataframe_like(self._meta)
        needs_time_conversion = False
        numeric_dd = self

        if is_df_like:
            time_cols = self._meta.select_dtypes(
                include=["datetime", "timedelta"]
            ).columns
            if len(time_cols) > 0:
                (
                    numeric_dd,
                    needs_time_conversion,
                ) = self._convert_time_cols_to_numeric(time_cols, axis, meta, skipna)
        else:
            needs_time_conversion = is_datetime64_any_dtype(self._meta)
            if needs_time_conversion:
                numeric_dd = _convert_to_numeric(self, skipna)

        if axis == 1:
            result = map_partitions(
                M.std if not needs_time_conversion else _sqrt_and_convert_to_timedelta,
                numeric_dd,
                meta=meta,
                token=self._token_prefix + "std",
                axis=axis,
                skipna=skipna,
                ddof=ddof,
                enforce_metadata=False,
                **numeric_kwargs,
                parent_meta=self._meta,
            )
            return handle_out(out, result)

        # Case where axis=0 or axis=None
        v = numeric_dd.var(skipna=skipna, ddof=ddof, split_every=split_every)
        name = self._token_prefix + "std"

        if needs_time_conversion:
            sqrt_func_kwargs = {
                "is_df_like": is_df_like,
                "time_cols": time_cols if is_df_like else None,
                "axis": axis,
                "dtype": getattr(meta, "dtype", None),
            }
            sqrt_func = _sqrt_and_convert_to_timedelta
        else:
            sqrt_func_kwargs = {}
            sqrt_func = np.sqrt

        result = map_partitions(
            sqrt_func,
            v,
            meta=meta,
            token=name,
            enforce_metadata=False,
            parent_meta=self._meta,
            **sqrt_func_kwargs,
        )

        # Try to match the Pandas result dtype
        if is_df_like and hasattr(meta, "dtype"):
            result = result.astype(meta.dtype)

        return handle_out(out, result)

    def _convert_time_cols_to_numeric(self, time_cols, axis, meta, skipna):
        from dask.dataframe.io import from_pandas

        needs_time_conversion = True

        # Ensure all columns are correct type. Need to shallow copy since cols will be modified
        if axis == 0:
            numeric_dd = self[meta.index].copy()
        else:
            numeric_dd = self.copy()

        # Mix of datetimes with other numeric types produces NaNs for each value in std() series
        if axis == 1 and len(time_cols) != len(self.columns):
            # This is faster than converting each column to numeric when it's not necessary
            # since each standard deviation will just be NaN
            needs_time_conversion = False
            numeric_dd = from_pandas(
                meta_frame_constructor(self)(
                    {"_": meta_series_constructor(self)([np.nan])},
                    index=self.index,
                ),
                npartitions=self.npartitions,
            )
        else:
            # Convert timedelta and datetime columns to integer types so we can use var
            for col in time_cols:
                numeric_dd[col] = _convert_to_numeric(numeric_dd[col], skipna)

        return numeric_dd, needs_time_conversion

    @_dummy_numpy_dispatcher("out", deprecated=True)
    @derived_from(pd.DataFrame)
    def skew(
        self,
        axis=0,
        bias=True,
        nan_policy="propagate",
        out=None,
        numeric_only=no_default,
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
        if PANDAS_GE_200 and axis is None:
            raise ValueError(
                "`axis=None` isn't currently supported for `skew` when using `pandas >=2` "
                f"(pandas={str(PANDAS_VERSION)} is installed)."
            )
        axis = self._validate_axis(axis)
        _raise_if_object_series(self, "skew")
        numeric_only_kwargs = get_numeric_only_kwargs(numeric_only)

        if is_dataframe_like(self):
            # Let pandas raise errors if necessary
            meta = self._meta_nonempty.skew(axis=axis, **numeric_only_kwargs)
        else:
            meta = self._meta_nonempty.skew()

        if axis == 1:
            result = map_partitions(
                M.skew,
                self,
                meta=meta,
                token=self._token_prefix + "skew",
                axis=axis,
                enforce_metadata=False,
            )
            return handle_out(out, result)
        else:
            if self.ndim == 1:
                result = self._skew_1d(self, bias=bias, nan_policy=nan_policy)
                return handle_out(out, result)
            else:
                result = self._skew_numeric(bias=bias, nan_policy=nan_policy)

            if isinstance(self, DataFrame):
                result.divisions = (self.columns.min(), self.columns.max())

            return handle_out(out, result)

    def _skew_1d(self, column, bias=True, nan_policy="propagate"):
        """1D version of the skew calculation.

        Uses the array version from da.stats in case we are passing in a single series
        """
        # import depends on scipy, not installed by default
        from dask.array import stats as da_stats

        if pd.api.types.is_extension_array_dtype(column._meta_nonempty):
            column = column.astype("f8")

        elif not np.issubdtype(column.dtype, np.number):
            column = column.astype("f8")

        name = self._token_prefix + "skew-1d-" + tokenize(column)

        array_skew = da_stats.skew(
            column.values, axis=0, bias=bias, nan_policy=nan_policy
        )

        layer = {(name, 0): (methods.wrap_skew_reduction, (array_skew._name,), None)}
        graph = HighLevelGraph.from_collections(name, layer, dependencies=[array_skew])

        return new_dd_object(
            graph, name, column._meta_nonempty.skew(), divisions=[None, None]
        )

    def _skew_numeric(self, bias=True, nan_policy="propagate"):
        """Method for dataframes with numeric columns.

        Maps the array version from da.stats onto the numeric array of columns.
        """
        # import depends on scipy, not installed by default
        from dask.array import stats as da_stats

        num = self.select_dtypes(include=["number", "bool"], exclude=[np.timedelta64])

        values_dtype = num.values.dtype
        array_values = num.values

        if not np.issubdtype(values_dtype, np.number):
            array_values = num.values.astype("f8")

        array_skew = da_stats.skew(
            array_values, axis=0, bias=bias, nan_policy=nan_policy
        )

        name = self._token_prefix + "var-numeric" + tokenize(num)
        cols = num._meta.columns if is_dataframe_like(num) else None

        with check_reductions_runtime_warning():
            skew_shape = num._meta_nonempty.var(axis=0).shape
        array_skew_name = (array_skew._name,) + (0,) * len(skew_shape)

        layer = {(name, 0): (methods.wrap_skew_reduction, array_skew_name, cols)}
        graph = HighLevelGraph.from_collections(name, layer, dependencies=[array_skew])

        return new_dd_object(
            graph, name, num._meta_nonempty.skew(), divisions=[None, None]
        )

    @_dummy_numpy_dispatcher("out", deprecated=True)
    @derived_from(pd.DataFrame)
    def kurtosis(
        self,
        axis=0,
        fisher=True,
        bias=True,
        nan_policy="propagate",
        out=None,
        numeric_only=no_default,
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
        if PANDAS_GE_200 and axis is None:
            raise ValueError(
                "`axis=None` isn't currently supported for `kurtosis` when using `pandas >=2` "
                f"(pandas={str(PANDAS_VERSION)} is installed)."
            )
        axis = self._validate_axis(axis)
        _raise_if_object_series(self, "kurtosis")
        numeric_only_kwargs = get_numeric_only_kwargs(numeric_only)

        if is_dataframe_like(self):
            # Let pandas raise errors if necessary
            meta = self._meta_nonempty.kurtosis(axis=axis, **numeric_only_kwargs)
        else:
            meta = self._meta_nonempty.kurtosis()

        if axis == 1:
            result = map_partitions(
                M.kurtosis,
                self,
                meta=meta,
                token=self._token_prefix + "kurtosis",
                axis=axis,
                enforce_metadata=False,
            )
            return handle_out(out, result)
        else:
            if self.ndim == 1:
                result = self._kurtosis_1d(
                    self, fisher=fisher, bias=bias, nan_policy=nan_policy
                )
                return handle_out(out, result)
            else:
                result = self._kurtosis_numeric(
                    fisher=fisher, bias=bias, nan_policy=nan_policy
                )

            if isinstance(self, DataFrame):
                result.divisions = (self.columns.min(), self.columns.max())

            return handle_out(out, result)

    def _kurtosis_1d(self, column, fisher=True, bias=True, nan_policy="propagate"):
        """1D version of the kurtosis calculation.

        Uses the array version from da.stats in case we are passing in a single series
        """
        # import depends on scipy, not installed by default
        from dask.array import stats as da_stats

        if pd.api.types.is_extension_array_dtype(column._meta_nonempty):
            column = column.astype("f8")

        elif not np.issubdtype(column.dtype, np.number):
            column = column.astype("f8")

        name = self._token_prefix + "kurtosis-1d-" + tokenize(column)

        array_kurtosis = da_stats.kurtosis(
            column.values, axis=0, fisher=fisher, bias=bias, nan_policy=nan_policy
        )

        layer = {
            (name, 0): (methods.wrap_kurtosis_reduction, (array_kurtosis._name,), None)
        }
        graph = HighLevelGraph.from_collections(
            name, layer, dependencies=[array_kurtosis]
        )

        return new_dd_object(
            graph, name, column._meta_nonempty.kurtosis(), divisions=[None, None]
        )

    def _kurtosis_numeric(self, fisher=True, bias=True, nan_policy="propagate"):
        """Method for dataframes with numeric columns.

        Maps the array version from da.stats onto the numeric array of columns.
        """
        # import depends on scipy, not installed by default
        from dask.array import stats as da_stats

        num = self.select_dtypes(include=["number", "bool"], exclude=[np.timedelta64])

        values_dtype = num.values.dtype
        array_values = num.values

        if not np.issubdtype(values_dtype, np.number):
            array_values = num.values.astype("f8")

        array_kurtosis = da_stats.kurtosis(
            array_values, axis=0, fisher=fisher, bias=bias, nan_policy=nan_policy
        )

        name = self._token_prefix + "kurtosis-numeric" + tokenize(num)
        cols = num._meta.columns if is_dataframe_like(num) else None

        with check_reductions_runtime_warning():
            kurtosis_shape = num._meta_nonempty.var(axis=0).shape
        array_kurtosis_name = (array_kurtosis._name,) + (0,) * len(kurtosis_shape)

        layer = {
            (name, 0): (methods.wrap_kurtosis_reduction, array_kurtosis_name, cols)
        }
        graph = HighLevelGraph.from_collections(
            name, layer, dependencies=[array_kurtosis]
        )

        return new_dd_object(
            graph, name, num._meta_nonempty.kurtosis(), divisions=[None, None]
        )

    @_numeric_only
    @derived_from(pd.DataFrame)
    def sem(self, axis=None, skipna=True, ddof=1, split_every=False, numeric_only=None):
        axis = self._validate_axis(axis)
        _raise_if_object_series(self, "sem")
        with check_numeric_only_deprecation():
            meta = self._meta_nonempty.sem(
                axis=axis, skipna=skipna, ddof=ddof, numeric_only=numeric_only
            )
        if axis == 1:
            return map_partitions(
                M.sem,
                self,
                meta=meta,
                token=self._token_prefix + "sem",
                axis=axis,
                skipna=skipna,
                ddof=ddof,
                parent_meta=self._meta,
                numeric_only=numeric_only,
            )
        else:
            num = self._get_numeric_data()
            v = num.var(skipna=skipna, ddof=ddof, split_every=split_every)
            n = num.count(split_every=split_every)
            name = self._token_prefix + "sem"
            result = map_partitions(
                np.sqrt,
                v / n,
                meta=meta,
                token=name,
                enforce_metadata=False,
                parent_meta=self._meta,
            )

            if isinstance(self, DataFrame):
                result.divisions = (self.columns.min(), self.columns.max())
            return result

    def quantile(self, q=0.5, axis=0, numeric_only=no_default, method="default"):
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
        axis = self._validate_axis(axis)
        keyname = "quantiles-concat--" + tokenize(self, q, axis)
        numeric_kwargs = get_numeric_only_kwargs(numeric_only)

        with check_numeric_only_deprecation("quantile", True):
            meta = self._meta.quantile(q, axis=axis, **numeric_kwargs)

        if axis == 1:
            if isinstance(q, list):
                # Not supported, the result will have current index as columns
                raise ValueError("'q' must be scalar when axis=1 is specified")

            return map_partitions(
                _getattr_numeric_only,
                self,
                q,
                axis,
                token=keyname,
                enforce_metadata=False,
                **numeric_kwargs,
                meta=(q, "f8"),
                parent_meta=self._meta,
                _dask_method_name="quantile",
            )
        else:
            _raise_if_object_series(self, "quantile")
            num = (
                self._get_numeric_data()
                if numeric_only is True
                or (not PANDAS_GE_200 and numeric_only is no_default)
                else self
            )
            quantiles = tuple(
                quantile(num.iloc[:, i], q, method) for i in range(len(num.columns))
            )

            qnames = [(_q._name, 0) for _q in quantiles]

            if isinstance(quantiles[0], Scalar):
                layer = {
                    (keyname, 0): (type(meta), qnames, num.columns, None, meta.name)
                }
                graph = HighLevelGraph.from_collections(
                    keyname, layer, dependencies=quantiles
                )
                divisions = (min(num.columns), max(num.columns))
                return Series(graph, keyname, meta, divisions)
            else:
                layer = {(keyname, 0): (methods.concat, qnames, 1)}
                graph = HighLevelGraph.from_collections(
                    keyname, layer, dependencies=quantiles
                )
                return new_dd_object(graph, keyname, meta, quantiles[0].divisions)

    @derived_from(pd.DataFrame)
    def describe(
        self,
        split_every=False,
        percentiles=None,
        percentiles_method="default",
        include=None,
        exclude=None,
        datetime_is_numeric=no_default,
    ):
        if PANDAS_GE_200:
            if datetime_is_numeric is no_default:
                datetime_is_numeric = True
                datetime_is_numeric_kwarg = {}
            else:
                raise TypeError(
                    "datetime_is_numeric is removed in pandas>=2.0.0, datetime data will always be "
                    "summarized as numeric"
                )
        else:
            datetime_is_numeric = (
                False if datetime_is_numeric is no_default else datetime_is_numeric
            )
            datetime_is_numeric_kwarg = {"datetime_is_numeric": datetime_is_numeric}

        if self._meta.ndim == 1:
            meta = self._meta_nonempty.describe(
                percentiles=percentiles,
                include=include,
                exclude=exclude,
                **datetime_is_numeric_kwarg,
            )
            output = self._describe_1d(
                self, split_every, percentiles, percentiles_method, datetime_is_numeric
            )
            output._meta = meta
            return output
        elif (include is None) and (exclude is None):
            _include = [np.number, np.timedelta64]
            if datetime_is_numeric:
                _include.append(np.datetime64)
            data = self._meta.select_dtypes(include=_include)

            # when some numerics/timedeltas are found, by default keep them
            if len(data.columns) == 0:
                chosen_columns_indexes = list(range(len(self._meta.columns)))
            else:
                # check if there are timedelta, boolean, or datetime columns
                _include = [np.timedelta64, bool]
                if datetime_is_numeric:
                    _include.append(np.datetime64)
                bools_and_times = self._meta.select_dtypes(include=_include)
                if len(bools_and_times.columns) == 0:
                    return self._describe_numeric(
                        self,
                        split_every,
                        percentiles,
                        percentiles_method,
                    )
                else:
                    chosen_columns_indexes = self._get_columns_indexes_based_on_dtypes(
                        data
                    )
        elif include == "all":
            if exclude is not None:
                msg = "exclude must be None when include is 'all'"
                raise ValueError(msg)
            chosen_columns_indexes = list(range(len(self._meta.columns)))
        else:
            data = self._meta.select_dtypes(include=include, exclude=exclude)
            chosen_columns_indexes = self._get_columns_indexes_based_on_dtypes(data)

        stats = [
            self._describe_1d(
                self.iloc[:, col_idx],
                split_every,
                percentiles,
                percentiles_method,
                datetime_is_numeric,
            )
            for col_idx in chosen_columns_indexes
        ]
        stats_names = [(s._name, 0) for s in stats]

        name = "describe--" + tokenize(self, split_every)
        layer = {(name, 0): (methods.describe_aggregate, stats_names)}
        graph = HighLevelGraph.from_collections(name, layer, dependencies=stats)
        meta = self._meta_nonempty.describe(
            include=include, exclude=exclude, **datetime_is_numeric_kwarg
        )
        return new_dd_object(graph, name, meta, divisions=[None, None])

    def _get_columns_indexes_based_on_dtypes(self, subset):
        meta = self._meta.dtypes.reset_index()
        meta.index.name = "indexer"
        return (
            meta.reset_index()
            .merge(subset.dtypes.reset_index(), how="inner")["indexer"]
            .values
        )

    def _describe_1d(
        self,
        data,
        split_every=False,
        percentiles=None,
        percentiles_method="default",
        datetime_is_numeric=False,
    ):
        if is_bool_dtype(data._meta):
            return self._describe_nonnumeric_1d(
                data, split_every=split_every, datetime_is_numeric=datetime_is_numeric
            )
        elif is_numeric_dtype(data._meta):
            return self._describe_numeric(
                data,
                split_every=split_every,
                percentiles=percentiles,
                percentiles_method=percentiles_method,
            )
        elif is_timedelta64_dtype(data._meta):
            return self._describe_numeric(
                data.dropna(),
                split_every=split_every,
                percentiles=percentiles,
                percentiles_method=percentiles_method,
                is_timedelta_column=True,
            )
        elif is_datetime64_any_dtype(data._meta) and datetime_is_numeric:
            return self._describe_numeric(
                data.dropna(),
                split_every=split_every,
                percentiles=percentiles,
                percentiles_method=percentiles_method,
                is_datetime_column=True,
            )
        else:
            return self._describe_nonnumeric_1d(
                data, split_every=split_every, datetime_is_numeric=datetime_is_numeric
            )

    def _describe_numeric(
        self,
        data,
        split_every=False,
        percentiles=None,
        percentiles_method="default",
        is_timedelta_column=False,
        is_datetime_column=False,
    ):
        from dask.dataframe.numeric import to_numeric

        if is_timedelta_column or is_datetime_column:
            num = to_numeric(data)
        else:
            num = data._get_numeric_data()

        if data.ndim == 2 and len(num.columns) == 0:
            raise ValueError("DataFrame contains only non-numeric data.")
        elif data.ndim == 1 and data.dtype == "object":
            raise ValueError("Cannot compute ``describe`` on object dtype.")
        if percentiles is None:
            percentiles = [0.25, 0.5, 0.75]
        else:
            # always include the the 50%tle to calculate the median
            # unique removes duplicates and sorts quantiles
            percentiles = np.array(percentiles)
            percentiles = np.append(percentiles, 0.5)
            percentiles = np.unique(percentiles)
            percentiles = list(percentiles)
        stats = [
            num.count(split_every=split_every),
            num.mean(split_every=split_every),
            num.std(split_every=split_every),
            num.min(split_every=split_every),
            num.quantile(percentiles, method=percentiles_method),
            num.max(split_every=split_every),
        ]
        stats_names = [(s._name, 0) for s in stats]

        colname = data._meta.name if is_series_like(data._meta) else None

        name = "describe-numeric--" + tokenize(num, split_every)
        layer = {
            (name, 0): (
                methods.describe_numeric_aggregate,
                stats_names,
                colname,
                is_timedelta_column,
                is_datetime_column,
            )
        }
        graph = HighLevelGraph.from_collections(name, layer, dependencies=stats)
        meta = num._meta_nonempty.describe()
        return new_dd_object(graph, name, meta, divisions=[None, None])

    def _describe_nonnumeric_1d(
        self, data, split_every=False, datetime_is_numeric=False
    ):
        from dask.dataframe.numeric import to_numeric

        vcounts = data.value_counts(split_every=split_every)
        count_nonzero = vcounts[vcounts != 0]
        count_unique = count_nonzero.size

        stats = [
            # nunique
            count_unique,
            # count
            data.count(split_every=split_every),
            # most common value
            vcounts._head(1, npartitions=1, compute=False, safe=False),
        ]

        if is_datetime64_any_dtype(data._meta) and not datetime_is_numeric:
            min_ts = to_numeric(data.dropna()).min(split_every=split_every)
            max_ts = to_numeric(data.dropna()).max(split_every=split_every)
            stats.extend([min_ts, max_ts])

        stats_names = [(s._name, 0) for s in stats]
        colname = data._meta.name

        name = "describe-nonnumeric-1d--" + tokenize(data, split_every)
        layer = {
            (name, 0): (methods.describe_nonnumeric_aggregate, stats_names, colname)
        }
        graph = HighLevelGraph.from_collections(name, layer, dependencies=stats)

        if not PANDAS_GE_200:
            datetime_is_numeric_kwarg = {"datetime_is_numeric": datetime_is_numeric}
        else:
            datetime_is_numeric_kwarg = {}
        meta = data._meta_nonempty.describe(**datetime_is_numeric_kwarg)
        return new_dd_object(graph, name, meta, divisions=[None, None])

    def _cum_agg(
        self,
        op_name,
        chunk,
        aggregate,
        axis,
        skipna=True,
        chunk_kwargs=None,
        out=None,  # Deprecated
    ):
        """Wrapper for cumulative operation"""

        axis = self._validate_axis(axis)

        if axis == 1:
            name = f"{self._token_prefix}{op_name}(axis=1)"
            result = self.map_partitions(chunk, token=name, **chunk_kwargs)
            return handle_out(out, result)
        else:
            # cumulate each partitions
            name1 = f"{self._token_prefix}{op_name}-map"
            cumpart = map_partitions(
                chunk, self, token=name1, meta=self, **chunk_kwargs
            )

            name2 = f"{self._token_prefix}{op_name}-take-last"
            cumlast = map_partitions(
                _take_last,
                cumpart,
                skipna,
                meta=meta_series_constructor(self)([], dtype="float"),
                token=name2,
            )

            suffix = tokenize(self)
            name = f"{self._token_prefix}{op_name}-{suffix}"
            cname = f"{self._token_prefix}{op_name}-cum-last-{suffix}"

            # aggregate cumulated partitions and its previous last element
            layer = {}
            layer[(name, 0)] = (cumpart._name, 0)

            for i in range(1, self.npartitions):
                # store each cumulative step to graph to reduce computation
                if i == 1:
                    layer[(cname, i)] = (cumlast._name, i - 1)
                else:
                    # aggregate with previous cumulation results
                    layer[(cname, i)] = (
                        methods._cum_aggregate_apply,
                        aggregate,
                        (cname, i - 1),
                        (cumlast._name, i - 1),
                    )
                layer[(name, i)] = (aggregate, (cumpart._name, i), (cname, i))
            graph = HighLevelGraph.from_collections(
                name, layer, dependencies=[cumpart, cumlast]
            )
            result = new_dd_object(graph, name, chunk(self._meta), self.divisions)
            return handle_out(out, result)

    @_dummy_numpy_dispatcher("dtype", "out", deprecated=True)
    @derived_from(pd.DataFrame)
    def cumsum(self, axis=None, skipna=True, dtype=None, out=None):
        return self._cum_agg(
            "cumsum",
            chunk=M.cumsum,
            aggregate=methods.cumsum_aggregate,
            axis=axis,
            skipna=skipna,
            chunk_kwargs=dict(axis=axis, skipna=skipna),
            out=out,
        )

    @_dummy_numpy_dispatcher("dtype", "out", deprecated=True)
    @derived_from(pd.DataFrame)
    def cumprod(self, axis=None, skipna=True, dtype=None, out=None):
        return self._cum_agg(
            "cumprod",
            chunk=M.cumprod,
            aggregate=methods.cumprod_aggregate,
            axis=axis,
            skipna=skipna,
            chunk_kwargs=dict(axis=axis, skipna=skipna),
            out=out,
        )

    @_dummy_numpy_dispatcher("out", deprecated=True)
    @derived_from(pd.DataFrame)
    def cummax(self, axis=None, skipna=True, out=None):
        return self._cum_agg(
            "cummax",
            chunk=M.cummax,
            aggregate=methods.cummax_aggregate,
            axis=axis,
            skipna=skipna,
            chunk_kwargs=dict(axis=axis, skipna=skipna),
            out=out,
        )

    @_dummy_numpy_dispatcher("out", deprecated=True)
    @derived_from(pd.DataFrame)
    def cummin(self, axis=None, skipna=True, out=None):
        return self._cum_agg(
            "cummin",
            chunk=M.cummin,
            aggregate=methods.cummin_aggregate,
            axis=axis,
            skipna=skipna,
            chunk_kwargs=dict(axis=axis, skipna=skipna),
            out=out,
        )

    def _validate_condition(self, cond):
        cond_res = cond(self._meta) if callable(cond) else cond
        if not (
            is_dask_collection(cond_res)
            or is_dataframe_like(cond_res)
            or is_series_like(cond_res)
            or is_index_like(cond_res)
        ):
            raise ValueError(
                f"Condition should be an object that can be aligned with {self.__class__}, "
                f" which includes Dask or pandas collections, DataFrames or Series, or a Callable."
            )

    @derived_from(pd.DataFrame)
    def where(self, cond, other=np.nan):
        # cond and other may be dask instance,
        # passing map_partitions via keyword will not be aligned
        self._validate_condition(cond)
        return map_partitions(M.where, self, cond, other, enforce_metadata=False)

    @derived_from(pd.DataFrame)
    def mask(self, cond, other=np.nan):
        self._validate_condition(cond)
        return map_partitions(M.mask, self, cond, other, enforce_metadata=False)

    @derived_from(pd.DataFrame)
    def notnull(self):
        return self.map_partitions(M.notnull, enforce_metadata=False)

    @derived_from(pd.DataFrame)
    def isnull(self):
        return self.map_partitions(M.isnull, enforce_metadata=False)

    @derived_from(pd.DataFrame)
    def isna(self):
        return self.map_partitions(M.isna, enforce_metadata=False)

    @derived_from(pd.DataFrame)
    def isin(self, values):
        if is_dataframe_like(self._meta):
            # DataFrame.isin does weird alignment stuff
            bad_types = (_Frame, pd.Series, pd.DataFrame)
        else:
            bad_types = (_Frame,)
        if isinstance(values, bad_types):
            raise NotImplementedError("Passing a %r to `isin`" % typename(type(values)))
        meta = self._meta_nonempty.isin(values)
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
        return self.map_partitions(
            M.isin, delayed(values), meta=meta, enforce_metadata=False
        )

    @derived_from(pd.DataFrame)
    def astype(self, dtype):
        meta = self._meta.astype(dtype)
        if hasattr(dtype, "items"):
            set_unknown = [
                k
                for k, v in dtype.items()
                if (isinstance(pd.api.types.pandas_dtype(v), pd.CategoricalDtype))
                and getattr(v, "categories", None) is None
            ]
            meta = clear_known_categories(meta, cols=set_unknown)
        elif (
            isinstance(pd.api.types.pandas_dtype(dtype), pd.CategoricalDtype)
        ) and getattr(dtype, "categories", None) is None:
            meta = clear_known_categories(meta)
        return self.map_partitions(
            M.astype, dtype=dtype, meta=meta, enforce_metadata=False
        )

    if not PANDAS_GE_200:

        @derived_from(pd.Series)
        def append(self, other, interleave_partitions=False):
            if PANDAS_GE_140:
                warnings.warn(
                    "The frame.append method is deprecated and will be removed from"
                    "dask in a future version. Use dask.dataframe.concat instead.",
                    FutureWarning,
                )
            # because DataFrame.append will override the method,
            # wrap by pd.Series.append docstring
            from dask.dataframe.multi import concat

            if isinstance(other, (list, dict)):
                msg = "append doesn't support list or dict input"
                raise NotImplementedError(msg)

            return concat(
                [self, other], join="outer", interleave_partitions=interleave_partitions
            )

    @derived_from(pd.Series)
    def dot(self, other, meta=no_default):
        if not isinstance(other, _Frame):
            raise TypeError("The second operand must be a dask array or dask dataframe")

        if isinstance(other, DataFrame):
            s = self.map_partitions(M.dot, other, token="dot", meta=meta)
            return s.groupby(by=s.index).apply(
                lambda x: x.sum(skipna=False), meta=s._meta_nonempty
            )

        def _dot_series(*args, **kwargs):
            # .sum() is invoked on each partition before being applied to all
            # partitions. The return type is expected to be a series, not a numpy object
            return meta_series_constructor(self)(M.dot(*args, **kwargs))

        return self.map_partitions(_dot_series, other, token="dot", meta=meta).sum(
            skipna=False
        )

    @derived_from(pd.DataFrame)
    def align(self, other, join="outer", axis=None, fill_value=None):
        meta1, meta2 = _emulate(
            M.align, self, other, join, axis=axis, fill_value=fill_value
        )
        aligned = self.map_partitions(
            M.align,
            other,
            join=join,
            axis=axis,
            fill_value=fill_value,
            enforce_metadata=False,
        )

        token = tokenize(self, other, join, axis, fill_value)

        name1 = "align1-" + token
        dsk1 = {
            (name1, i): (getitem, key, 0)
            for i, key in enumerate(aligned.__dask_keys__())
        }
        dsk1.update(aligned.dask)
        result1 = new_dd_object(dsk1, name1, meta1, aligned.divisions)

        name2 = "align2-" + token
        dsk2 = {
            (name2, i): (getitem, key, 1)
            for i, key in enumerate(aligned.__dask_keys__())
        }
        dsk2.update(aligned.dask)
        result2 = new_dd_object(dsk2, name2, meta2, aligned.divisions)

        return result1, result2

    @derived_from(pd.DataFrame)
    def combine(self, other, func, fill_value=None, overwrite=True):
        return self.map_partitions(
            M.combine, other, func, fill_value=fill_value, overwrite=overwrite
        )

    @derived_from(pd.DataFrame)
    def combine_first(self, other):
        return self.map_partitions(M.combine_first, other)

    @classmethod
    def _bind_operator_method(cls, name, op, original=pd.DataFrame):
        """bind operator method like DataFrame.add to this class"""
        raise NotImplementedError

    @derived_from(pd.DataFrame)
    def resample(self, rule, closed=None, label=None):
        from dask.dataframe.tseries.resample import Resampler

        return Resampler(self, rule, closed=closed, label=label)

    @_deprecated(
        message=(
            "Will be removed in a future version. "
            "Please create a mask and filter using .loc instead"
        )
    )
    @derived_from(pd.DataFrame)
    def first(self, offset):
        # Let pandas error on bad args
        self._meta_nonempty.first(offset)

        if not self.known_divisions:
            raise ValueError("`first` is not implemented for unknown divisions")

        offset = pd.tseries.frequencies.to_offset(offset)
        date = self.divisions[0] + offset
        end = self.loc._get_partitions(date)

        is_anchored = offset.is_anchored()

        include_right = is_anchored or not hasattr(offset, "delta")

        if end == self.npartitions - 1:
            divs = self.divisions
        else:
            divs = self.divisions[: end + 1] + (date,)

        name = "first-" + tokenize(self, offset)
        dsk = {(name, i): (self._name, i) for i in range(end)}
        dsk[(name, end)] = (
            methods.boundary_slice,
            (self._name, end),
            None,
            date,
            include_right,
            True,
        )
        graph = HighLevelGraph.from_collections(name, dsk, dependencies=[self])
        return new_dd_object(graph, name, self, divs)

    @_deprecated(
        message=(
            "Will be removed in a future version. "
            "Please create a mask and filter using .loc instead"
        )
    )
    @derived_from(pd.DataFrame)
    def last(self, offset):
        # Let pandas error on bad args
        self._meta_nonempty.last(offset)

        if not self.known_divisions:
            raise ValueError("`last` is not implemented for unknown divisions")

        offset = pd.tseries.frequencies.to_offset(offset)
        date = self.divisions[-1] - offset
        start = self.loc._get_partitions(date)

        if start == 0:
            divs = self.divisions
        else:
            divs = (date,) + self.divisions[start + 1 :]

        name = "last-" + tokenize(self, offset)
        dsk = {
            (name, i + 1): (self._name, j + 1)
            for i, j in enumerate(range(start, self.npartitions))
        }
        dsk[(name, 0)] = (
            methods.boundary_slice,
            (self._name, start),
            date,
            None,
            True,
            False,
        )
        graph = HighLevelGraph.from_collections(name, dsk, dependencies=[self])
        return new_dd_object(graph, name, self, divs)

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
        from dask.dataframe import hyperloglog  # here to avoid circular import issues

        return aca(
            [self],
            chunk=hyperloglog.compute_hll_array,
            combine=hyperloglog.reduce_state,
            aggregate=hyperloglog.estimate_count,
            split_every=split_every,
            b=16,
            meta=float,
        )

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

    def _validate_chunks(self, arr, lengths):
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

    def _is_index_level_reference(self, key):
        """
        Test whether a key is an index level reference

        To be considered an index level reference, `key` must match the index name
        and must NOT match the name of any column (if a dataframe).
        """
        return (
            self.index.name is not None
            and not is_dask_collection(key)
            and (np.isscalar(key) or isinstance(key, tuple))
            and key == self.index.name
            and key not in getattr(self, "columns", ())
        )

    def _contains_index_name(self, columns_or_index):
        """
        Test whether the input contains a reference to the index of the DataFrame/Series
        """
        if isinstance(columns_or_index, list):
            return any(self._is_index_level_reference(n) for n in columns_or_index)
        else:
            return self._is_index_level_reference(columns_or_index)


def _raise_if_object_series(x, funcname):
    """
    Utility function to raise an error if an object column does not support
    a certain operation like `mean`.
    """
    if isinstance(x, Series) and hasattr(x, "dtype"):
        if x.dtype == object:
            raise ValueError("`%s` not supported with object series" % funcname)
        elif is_string_dtype(x):
            raise ValueError("`%s` not supported with string series" % funcname)


class Series(_Frame):
    """Parallel Pandas Series

    Do not use this class directly.  Instead use functions like
    ``dd.read_csv``, ``dd.read_parquet``, or ``dd.from_pandas``.

    Parameters
    ----------

    dsk: dict
        The dask graph to compute this Series
    _name: str
        The key prefix that specifies which keys in the dask comprise this
        particular Series
    meta: pandas.Series
        An empty ``pandas.Series`` with names, dtypes, and index matching the
        expected output.
    divisions: tuple of index values
        Values along which we partition our blocks on the index

    See Also
    --------
    dask.dataframe.DataFrame
    """

    _partition_type = pd.Series
    _is_partition_type = staticmethod(is_series_like)
    _token_prefix = "series-"
    _accessors: ClassVar[set[str]] = set()

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

    @property
    def axes(self):
        return [self.index]

    @property
    def name(self):
        return self._meta.name

    @name.setter
    def name(self, name):
        self._meta.name = name
        renamed = _rename_dask(self, name)
        # update myself
        self.dask = renamed.dask
        self._name = renamed._name

    @property
    def ndim(self):
        """Return dimensionality"""
        return 1

    @property
    def shape(self):
        """
        Return a tuple representing the dimensionality of a Series.

        The single element of the tuple is a Delayed result.

        Examples
        --------
        >>> series.shape  # doctest: +SKIP
        (dd.Scalar<size-ag..., dtype=int64>,)
        """
        return (self.size,)

    @property
    def dtype(self):
        """Return data type"""
        return self._meta.dtype

    dt = CachedAccessor("dt", DatetimeAccessor)

    cat = CachedAccessor("cat", CategoricalAccessor)

    str = CachedAccessor("str", StringAccessor)

    def __dir__(self):
        o = set(dir(type(self)))
        o.update(self.__dict__)
        # Remove the `cat` and `str` accessors if not available. We can't
        # decide this statically for the `dt` accessor, as it works on
        # datetime-like things as well.
        for accessor in ["cat", "str"]:
            if not hasattr(self._meta, accessor):
                o.remove(accessor)
        return list(o)

    @property
    def nbytes(self):
        """Number of bytes"""
        return self.reduction(
            methods.nbytes, np.sum, token="nbytes", meta=int, split_every=False
        )

    def _repr_data(self):
        return _repr_data_series(self._meta, self._repr_divisions)

    def __repr__(self):
        """have to overwrite footer"""
        if self.name is not None:
            footer = f"Name: {self.name}, dtype: {self.dtype}"
        else:
            footer = f"dtype: {self.dtype}"

        return """Dask {klass} Structure:
{data}
{footer}
Dask Name: {name}, {layers}""".format(
            klass=self.__class__.__name__,
            data=self.to_string(),
            footer=footer,
            name=key_split(self._name),
            layers=maybe_pluralize(len(self.dask.layers), "graph layer"),
        )

    def rename(self, index=None, inplace=no_default, sorted_index=False):
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
        from pandas.api.types import is_dict_like, is_list_like, is_scalar

        import dask.dataframe as dd

        if inplace is not no_default:
            warnings.warn(
                "'inplace' argument for dask series will be removed in future versions",
                FutureWarning,
            )
        else:
            inplace = False

        if is_scalar(index) or (
            is_list_like(index)
            and not is_dict_like(index)
            and not isinstance(index, dd.Series)
        ):
            res = self if inplace else self.copy()
            res.name = index
        else:
            res = self.map_partitions(M.rename, index, enforce_metadata=False)
            if self.known_divisions:
                if sorted_index and (callable(index) or is_dict_like(index)):
                    old = pd.Series(range(self.npartitions + 1), index=self.divisions)
                    new = old.rename(index).index
                    if not new.is_monotonic_increasing:
                        msg = (
                            "sorted_index=True, but the transformed index "
                            "isn't monotonic_increasing"
                        )
                        raise ValueError(msg)
                    res._divisions = tuple(methods.tolist(new))
                else:
                    res = res.clear_divisions()
            if inplace:
                self.dask = res.dask
                self._name = res._name
                self._divisions = res.divisions
                self._meta = res._meta
                res = self
        return res

    @derived_from(pd.Series)
    def round(self, decimals=0):
        return elemwise(M.round, self, decimals)

    @derived_from(pd.Series)
    def to_timestamp(self, freq=None, how="start", axis=0):
        df = elemwise(M.to_timestamp, self, freq, how, axis)
        df.divisions = tuple(pd.Index(self.divisions).to_timestamp(freq=freq, how=how))
        return df

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
        return quantile(self, q, method=method)

    def median_approximate(self, method="default"):
        """Return the approximate median of the values over the requested axis.

        Parameters
        ----------
        method : {'default', 'tdigest', 'dask'}, optional
            What method to use. By default will use Dask's internal custom
            algorithm (``"dask"``).  If set to ``"tdigest"`` will use tdigest
            for floats and ints and fallback to the ``"dask"`` otherwise.
        """
        return self.quantile(q=0.5, method=method)

    @derived_from(pd.Series)
    def median(self, method="default"):
        if self.npartitions == 1:
            # Can provide an exact median in these cases
            return self.median_approximate(method=method)
        raise NotImplementedError(
            "Dask doesn't implement an exact median in all cases as this is hard to do in parallel. "
            "See the `median_approximate` method instead, which uses an approximate algorithm."
        )

    def _repartition_quantiles(self, npartitions, upsample=1.0):
        """Approximate quantiles of Series used for repartitioning"""
        from dask.dataframe.partitionquantiles import partition_quantiles

        return partition_quantiles(self, npartitions, upsample=upsample)

    def __getitem__(self, key):
        if isinstance(key, Series) and self.divisions == key.divisions:
            name = "index-%s" % tokenize(self, key)
            dsk = partitionwise_graph(operator.getitem, name, self, key)
            graph = HighLevelGraph.from_collections(name, dsk, dependencies=[self, key])
            return Series(graph, name, self._meta, self.divisions)
        return self.loc[key]

    @derived_from(pd.DataFrame)
    def _get_numeric_data(self, how="any", subset=None):
        return self

    if not PANDAS_GE_200:

        @derived_from(pd.Series)
        def iteritems(self):
            if PANDAS_GE_150:
                warnings.warn(
                    "iteritems is deprecated and will be removed in a future version. "
                    "Use .items instead.",
                    FutureWarning,
                )
            # We use the `_` generator below to ensure the deprecation warning above
            # is raised when `.iteritems()` is called, not when the first `next(<generator>)`
            # iteration happens

            def _(self):
                for i in range(self.npartitions):
                    s = self.get_partition(i).compute()
                    yield from s.items()

            return _(self)

    @derived_from(pd.Series)
    def __iter__(self):
        for i in range(self.npartitions):
            s = self.get_partition(i).compute()
            yield from s

    @_deprecated(
        message=(
            "Using the ``in`` operator to test for membership in Series is "
            "deprecated. To test for membership in the index use "
            "``(s.index == key).any()``. Similarly to test for membership in "
            "the values use ``(s == key).any()``"
        )
    )
    def __contains__(self, key):
        return (self == key).any().compute()

    @classmethod
    def _validate_axis(cls, axis=0, none_is_zero: bool = True) -> None | Literal[0, 1]:
        if axis not in (0, "index", None):
            raise ValueError(f"No axis named {axis}")
        # convert to numeric axis
        numeric_axis: dict[str | None, Literal[0, 1]] = {"index": 0}
        if none_is_zero:
            numeric_axis[None] = 0
        return numeric_axis.get(axis, axis)

    @derived_from(pd.Series)
    def groupby(
        self,
        by=None,
        group_keys=GROUP_KEYS_DEFAULT,
        sort=None,
        observed=None,
        dropna=None,
        **kwargs,
    ):
        from dask.dataframe.groupby import SeriesGroupBy

        return SeriesGroupBy(
            self,
            by=by,
            group_keys=group_keys,
            sort=sort,
            observed=observed,
            dropna=dropna,
            **kwargs,
        )

    @derived_from(pd.Series)
    def count(self, split_every=False):
        return super().count(split_every=split_every)

    @derived_from(pd.Series)
    def mode(self, dropna=True, split_every=False):
        return super().mode(dropna=dropna, split_every=split_every)

    @derived_from(pd.Series)
    def explode(self):
        meta = self._meta.explode()
        return self.map_partitions(M.explode, meta=meta, enforce_metadata=False)

    def unique(self, split_every=None, split_out=1):
        """
        Return Series of unique values in the object. Includes NA values.

        Returns
        -------
        uniques : Series
        """
        return aca(
            self,
            chunk=methods.unique,
            aggregate=methods.unique,
            meta=self._meta,
            token="unique",
            split_every=split_every,
            series_name=self.name,
            split_out=split_out,
        )

    @derived_from(pd.Series)
    def nunique(self, split_every=None, dropna=True):
        uniqs = self.drop_duplicates(split_every=split_every)
        if dropna:
            # count mimics pandas behavior and excludes NA values
            return uniqs.count()
        else:
            return uniqs.size

    @derived_from(pd.Series)
    def value_counts(
        self,
        sort=None,
        ascending=False,
        dropna=True,
        normalize=False,
        split_every=None,
        split_out=1,
    ):
        kwargs = {"sort": sort, "ascending": ascending, "dropna": dropna}
        aggregate_kwargs = {"normalize": normalize}
        if split_out > 1:
            aggregate_kwargs["total_length"] = (
                len(self) if dropna is False else len(self.dropna())
            )

        return aca(
            self,
            chunk=M.value_counts,
            aggregate=methods.value_counts_aggregate,
            combine=methods.value_counts_combine,
            meta=self._meta.value_counts(normalize=normalize),
            token="value-counts",
            split_every=split_every,
            split_out=split_out,
            split_out_setup=split_out_on_index,
            aggregate_kwargs=aggregate_kwargs,
            **kwargs,
        )

    @derived_from(pd.Series)
    def nlargest(self, n=5, split_every=None):
        return aca(
            self,
            chunk=M.nlargest,
            aggregate=M.nlargest,
            meta=self._meta,
            token="series-nlargest",
            split_every=split_every,
            n=n,
        )

    @derived_from(pd.Series)
    def nsmallest(self, n=5, split_every=None):
        return aca(
            self,
            chunk=M.nsmallest,
            aggregate=M.nsmallest,
            meta=self._meta,
            token="series-nsmallest",
            split_every=split_every,
            n=n,
        )

    @derived_from(pd.Series)
    def isin(self, values):
        # Added just to get the different docstring for Series
        return super().isin(values)

    @insert_meta_param_description(pad=12)
    @derived_from(pd.Series)
    def map(self, arg, na_action=None, meta=no_default):
        if is_series_like(arg) and is_dask_collection(arg):
            return series_map(self, arg)
        if not (
            isinstance(arg, dict)
            or callable(arg)
            or is_series_like(arg)
            and not is_dask_collection(arg)
        ):
            raise TypeError(
                f"arg must be pandas.Series, dict or callable. Got {type(arg)}"
            )
        name = "map-" + tokenize(self, arg, na_action)
        dsk = {
            (name, i): (M.map, k, arg, na_action)
            for i, k in enumerate(self.__dask_keys__())
        }
        graph = HighLevelGraph.from_collections(name, dsk, dependencies=[self])
        if meta is no_default:
            meta = _emulate(M.map, self, arg, na_action=na_action, udf=True)
        else:
            meta = make_meta(
                meta,
                index=getattr(make_meta(self), "index", None),
                parent_meta=self._meta,
            )

        return type(self)(graph, name, meta, self.divisions)

    @derived_from(pd.Series)
    def dropna(self):
        return self.map_partitions(M.dropna, enforce_metadata=False)

    @derived_from(pd.Series)
    def between(self, left, right, inclusive="both"):
        return self.map_partitions(
            M.between, left=left, right=right, inclusive=inclusive
        )

    @_dummy_numpy_dispatcher("out")
    @derived_from(pd.Series)
    def clip(self, lower=None, upper=None, axis=None):
        if axis not in (None, 0):
            raise ValueError(f"Series.clip does not support axis={axis}")
        # np.clip may pass out
        return self.map_partitions(
            M.clip,
            lower=lower,
            upper=upper,
            enforce_metadata=False,
            # See: https://github.com/rapidsai/cudf/issues/14369
            **({} if axis is None else {"axis": axis}),
        )

    @derived_from(pd.Series)
    def align(self, other, join="outer", axis=None, fill_value=None):
        return super().align(other, join=join, axis=axis, fill_value=fill_value)

    @derived_from(pd.Series)
    def combine(self, other, func, fill_value=None):
        return self.map_partitions(M.combine, other, func, fill_value=fill_value)

    @derived_from(pd.Series)
    def squeeze(self):
        return self

    @derived_from(pd.Series)
    def combine_first(self, other):
        return self.map_partitions(M.combine_first, other)

    def to_bag(self, index=False, format="tuple"):
        """Create a Dask Bag from a Series"""
        from dask.dataframe.io import to_bag

        return to_bag(self, index, format=format)

    @derived_from(pd.Series)
    def to_frame(self, name=None):
        args = [] if name is None else [name]
        return self.map_partitions(M.to_frame, *args, meta=self._meta.to_frame(*args))

    @derived_from(pd.Series)
    def to_string(self, max_rows=5):
        # option_context doesn't affect
        return self._repr_data().to_string(max_rows=max_rows)

    @classmethod
    def _bind_operator_method(cls, name, op, original=pd.Series):
        """bind operator method like Series.add to this class"""

        def meth(self, other, level=None, fill_value=None, axis=0):
            if level is not None:
                raise NotImplementedError("level must be None")
            axis = self._validate_axis(axis)
            meta = _emulate(op, self, other, axis=axis, fill_value=fill_value)
            return map_partitions(
                op, self, other, meta=meta, axis=axis, fill_value=fill_value
            )

        meth.__name__ = name
        setattr(cls, name, derived_from(original)(meth))

    @classmethod
    def _bind_comparison_method(cls, name, comparison, original=pd.Series):
        """bind comparison method like Series.eq to this class"""

        def meth(self, other, level=None, fill_value=None, axis=0):
            if level is not None:
                raise NotImplementedError("level must be None")
            axis = self._validate_axis(axis)
            if fill_value is None:
                return elemwise(comparison, self, other, axis=axis)
            else:
                op = partial(comparison, fill_value=fill_value)
                return elemwise(op, self, other, axis=axis)

        meth.__name__ = name
        setattr(cls, name, derived_from(original)(meth))

    @insert_meta_param_description(pad=12)
    def apply(self, func, convert_dtype=no_default, meta=no_default, args=(), **kwds):
        """Parallel version of pandas.Series.apply

        Parameters
        ----------
        func : function
            Function to apply
        convert_dtype : boolean, default True
            Try to find better dtype for elementwise function results.
            If False, leave as dtype=object.
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
        dask.Series.map_partitions
        """
        if convert_dtype is not no_default:
            kwds["convert_dtype"] = convert_dtype

        # let pandas trigger any warnings, such as convert_dtype warning
        self._meta.apply(func, args=args, **kwds)

        if meta is no_default:
            with check_convert_dtype_deprecation():
                meta = _emulate(
                    M.apply, self._meta_nonempty, func, args=args, udf=True, **kwds
                )
            warnings.warn(meta_warning(meta))

        return map_partitions(methods.apply, self, func, args=args, meta=meta, **kwds)

    @derived_from(pd.Series)
    def cov(self, other, min_periods=None, split_every=False):
        from dask.dataframe.multi import concat

        if not isinstance(other, Series):
            raise TypeError("other must be a dask.dataframe.Series")
        df = concat([self, other], axis=1)
        return _cov_corr(
            df, min_periods, scalar=True, numeric_only=False, split_every=split_every
        )

    @derived_from(pd.Series)
    def corr(self, other, method="pearson", min_periods=None, split_every=False):
        from dask.dataframe.multi import concat

        if not isinstance(other, Series):
            raise TypeError("other must be a dask.dataframe.Series")
        if method != "pearson":
            raise NotImplementedError("Only Pearson correlation has been implemented")
        df = concat([self, other], axis=1)
        return _cov_corr(
            df,
            min_periods,
            corr=True,
            scalar=True,
            numeric_only=False,
            split_every=split_every,
        )

    @derived_from(pd.Series)
    def autocorr(self, lag=1, split_every=False):
        if not isinstance(lag, Integral):
            raise TypeError("lag must be an integer")
        return self.corr(self if lag == 0 else self.shift(lag), split_every=split_every)

    @derived_from(pd.Series)
    def memory_usage(self, index=True, deep=False):
        return self.reduction(
            M.memory_usage,
            M.sum,
            chunk_kwargs={"index": index, "deep": deep},
            split_every=False,
            token=self._token_prefix + "memory-usage",
        )

    def __divmod__(self, other):
        res1 = self // other
        res2 = self % other
        return res1, res2

    def __rdivmod__(self, other):
        res1 = other // self
        res2 = other % self
        return res1, res2

    if not PANDAS_GE_200:

        @property
        @derived_from(pd.Series)
        def is_monotonic(self):
            if PANDAS_GE_150:
                warnings.warn(
                    "is_monotonic is deprecated and will be removed in a future version. "
                    "Use is_monotonic_increasing instead.",
                    FutureWarning,
                )
            return self.is_monotonic_increasing

    @property
    @derived_from(pd.Series)
    def is_monotonic_increasing(self):
        return aca(
            self,
            chunk=methods.monotonic_increasing_chunk,
            combine=methods.monotonic_increasing_combine,
            aggregate=methods.monotonic_increasing_aggregate,
            meta=bool,
            token="monotonic_increasing",
        )

    @property
    @derived_from(pd.Series)
    def is_monotonic_decreasing(self):
        return aca(
            self,
            chunk=methods.monotonic_decreasing_chunk,
            combine=methods.monotonic_decreasing_combine,
            aggregate=methods.monotonic_decreasing_aggregate,
            meta=bool,
            token="monotonic_decreasing",
        )

    @_deprecated(
        message=(
            "Will be removed in a future version. "
            "Use `Series.astype()` as an alternative to change the dtype."
        )
    )
    @derived_from(pd.Series)
    def view(self, dtype):
        return self._view(dtype)

    def _view(self, dtype):
        meta = self._meta.view(dtype)
        return self.map_partitions(M.view, dtype, meta=meta)


class Index(Series):
    _partition_type = pd.Index
    _is_partition_type = staticmethod(is_index_like)
    _token_prefix = "index-"
    _accessors: ClassVar[set[str]] = set()

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

    _monotonic_attributes = {
        "is_monotonic",
        "is_monotonic_increasing",
        "is_monotonic_decreasing",
    }

    def __getattr__(self, key):
        if (
            isinstance(self._meta.dtype, pd.CategoricalDtype)
            and key in self._cat_attributes
        ):
            return getattr(self.cat, key)
        elif key in self._dt_attributes:
            return getattr(self.dt, key)
        elif key in self._monotonic_attributes:
            return getattr(self, key)
        raise AttributeError("'Index' object has no attribute %r" % key)

    def __dir__(self):
        out = super().__dir__()
        out.extend(self._dt_attributes)
        if isinstance(self.dtype, pd.CategoricalDtype):
            out.extend(self._cat_attributes)
        return out

    # Typing: https://github.com/python/mypy/issues/4125
    @property  # type: ignore
    def index(self):
        raise AttributeError(
            f"{self.__class__.__name__!r} object has no attribute 'index'"
        )

    def __array_wrap__(self, array, context=None):
        return pd.Index(array, name=self.name)

    def head(self, n=5, compute=True):
        """First n items of the Index.

        Caveat, this only checks the first partition.
        """
        name = "head-%d-%s" % (n, self._name)
        dsk = {(name, 0): (operator.getitem, (self._name, 0), slice(0, n))}
        graph = HighLevelGraph.from_collections(name, dsk, dependencies=[self])

        result = new_dd_object(graph, name, self._meta, self.divisions[:2])

        if compute:
            result = result.compute()
        return result

    @derived_from(pd.Index)
    def max(self, split_every=False):
        return self.reduction(
            M.max,
            meta=self._meta_nonempty.max(),
            token=self._token_prefix + "max",
            split_every=split_every,
        )

    @derived_from(pd.Index)
    def min(self, split_every=False):
        return self.reduction(
            M.min,
            meta=self._meta_nonempty.min(),
            token=self._token_prefix + "min",
            split_every=split_every,
        )

    def count(self, split_every=False):
        return self.reduction(
            methods.index_count,
            np.sum,
            token="index-count",
            meta=int,
            split_every=split_every,
        )

    @derived_from(pd.Index)
    def shift(self, periods=1, freq=None):
        if isinstance(self._meta, pd.PeriodIndex):
            if freq is not None:
                raise ValueError("PeriodIndex doesn't accept `freq` argument")
            meta = self._meta_nonempty.shift(periods)
            out = self.map_partitions(
                M.shift, periods, meta=meta, token="shift", transform_divisions=False
            )
        else:
            # Pandas will raise for other index types that don't implement shift
            meta = self._meta_nonempty.shift(periods, freq=freq)
            out = self.map_partitions(
                M.shift,
                periods,
                token="shift",
                meta=meta,
                freq=freq,
                transform_divisions=False,
            )
        if freq is None:
            freq = meta.freq
        return maybe_shift_divisions(out, periods, freq=freq)

    @derived_from(pd.Index)
    def to_series(self):
        return self.map_partitions(
            M.to_series,
            meta=self._meta.to_series(),
            transform_divisions=False,
        )

    @derived_from(pd.Index, ua_args=["index"])
    def to_frame(self, index=True, name=None):
        if not index:
            raise NotImplementedError()
        args = [index] if name is None else [index, name]

        return self.map_partitions(
            M.to_frame,
            *args,
            meta=self._meta.to_frame(*args),
            transform_divisions=False,
        )

    @insert_meta_param_description(pad=12)
    @derived_from(pd.Index)
    def map(self, arg, na_action=None, meta=no_default, is_monotonic=False):
        """
        Note that this method clears any known divisions.

        If your mapping function is monotonically increasing then use `is_monotonic`
        to apply the mapping function to the old divisions and assign the new
        divisions to the output.

        """
        applied = super().map(arg, na_action=na_action, meta=meta)
        if is_monotonic and self.known_divisions:
            applied.divisions = tuple(
                pd.Series(self.divisions).map(arg, na_action=na_action)
            )
        else:
            applied = applied.clear_divisions()
        return applied

    if not PANDAS_GE_200:

        @property
        @derived_from(pd.Index)
        def is_monotonic(self):
            if PANDAS_GE_150:
                warnings.warn(
                    "is_monotonic is deprecated and will be removed in a future version. "
                    "Use is_monotonic_increasing instead.",
                    FutureWarning,
                )
            return super().is_monotonic_increasing

    @property
    @derived_from(pd.Index)
    def is_monotonic_increasing(self):
        return super().is_monotonic_increasing

    @property
    @derived_from(pd.Index)
    def is_monotonic_decreasing(self):
        return super().is_monotonic_decreasing

    @derived_from(pd.Index)
    def memory_usage(self, deep=False):
        return self.reduction(
            M.memory_usage,
            M.sum,
            chunk_kwargs={"deep": deep},
            split_every=False,
            token=self._token_prefix + "memory-usage",
        )

    @_deprecated_kwarg("shuffle", "shuffle_method")
    @derived_from(
        pd.Index,
        inconsistencies="keep=False will raise a ``NotImplementedError``",
    )
    def drop_duplicates(
        self,
        split_every=None,
        split_out=1,
        shuffle_method=None,
        **kwargs,
    ):
        if not self.known_divisions:
            # Use base class if we have unknown divisions
            return super().drop_duplicates(
                split_every=split_every,
                split_out=split_out,
                shuffle_method=shuffle_method,
                **kwargs,
            )

        # Let pandas error on bad inputs
        self._meta_nonempty.drop_duplicates(**kwargs)

        # Raise error for unsupported `keep`
        if kwargs.get("keep", True) is False:
            raise NotImplementedError("drop_duplicates with keep=False")

        # Simple `drop_duplicates` case that we are acting on
        # an Index with known divisions
        chunk = M.drop_duplicates
        repartition_npartitions = max(
            self.npartitions // (split_every or self.npartitions),
            split_out,
        )
        assert self.known_divisions, "Requires known divisions"
        return (
            self.map_partitions(
                chunk,
                token="drop-duplicates-chunk",
                meta=self._meta,
                transform_divisions=False,
                **kwargs,
            )
            .repartition(npartitions=repartition_npartitions)
            .map_partitions(
                chunk,
                token="drop-duplicates-agg",
                meta=self._meta,
                transform_divisions=False,
                **kwargs,
            )
            .repartition(npartitions=split_out)
        )


class DataFrame(_Frame):
    """
    Parallel Pandas DataFrame

    Do not use this class directly.  Instead use functions like
    ``dd.read_csv``, ``dd.read_parquet``, or ``dd.from_pandas``.

    Parameters
    ----------
    dsk: dict
        The dask graph to compute this DataFrame
    name: str
        The key prefix that specifies which keys in the dask comprise this
        particular DataFrame
    meta: pandas.DataFrame
        An empty ``pandas.DataFrame`` with names, dtypes, and index matching
        the expected output.
    divisions: tuple of index values
        Values along which we partition our blocks on the index
    """

    _partition_type = pd.DataFrame
    _is_partition_type = staticmethod(is_dataframe_like)
    _token_prefix = "dataframe-"
    _accessors: ClassVar[set[str]] = set()

    def __init__(self, dsk, name, meta, divisions):
        super().__init__(dsk, name, meta, divisions)
        if self.dask.layers[name].collection_annotations is None:
            self.dask.layers[name].collection_annotations = {
                "npartitions": self.npartitions,
                "columns": [col for col in self.columns],
                "type": typename(type(self)),
                "dataframe_type": typename(type(self._meta)),
                "series_dtypes": {
                    col: self._meta[col].dtype
                    if hasattr(self._meta[col], "dtype")
                    else None
                    for col in self._meta.columns
                },
            }
        else:
            self.dask.layers[name].collection_annotations.update(
                {
                    "npartitions": self.npartitions,
                    "columns": [col for col in self.columns],
                    "type": typename(type(self)),
                    "dataframe_type": typename(type(self._meta)),
                    "series_dtypes": {
                        col: self._meta[col].dtype
                        if hasattr(self._meta[col], "dtype")
                        else None
                        for col in self._meta.columns
                    },
                }
            )

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

    @property
    def axes(self):
        return [self.index, self.columns]

    @property
    def columns(self):
        return self._meta.columns

    @columns.setter
    def columns(self, columns):
        renamed = _rename_dask(self, columns)
        self._meta = renamed._meta
        self._name = renamed._name
        self.dask = renamed.dask

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
        from dask.dataframe.indexing import _iLocIndexer

        # For dataframes with unique column names, this will be transformed into a __getitem__ call
        return _iLocIndexer(self)

    def __len__(self):
        try:
            s = self.iloc[:, 0]
        except IndexError:
            return super().__len__()
        else:
            return len(s)

    def __contains__(self, key):
        return key in self._meta

    @property
    def empty(self):
        # __getattr__ will be called after we raise this, so we'll raise it again from there
        raise AttributeNotImplementedError(
            "Checking whether a Dask DataFrame has any rows may be expensive. "
            "However, checking the number of columns is fast. "
            "Depending on which of these results you need, use either "
            "`len(df.index) == 0` or `len(df.columns) == 0`"
        )

    def __getitem__(self, key):
        name = "getitem-%s" % tokenize(self, key)
        if np.isscalar(key) or isinstance(key, (tuple, str)):
            if isinstance(self._meta.index, (pd.DatetimeIndex, pd.PeriodIndex)):
                if key not in self._meta.columns:
                    warnings.warn(
                        "Indexing a DataFrame with a datetimelike index using a single "
                        "string to slice the rows, like `frame[string]`, is deprecated "
                        "and will be removed in a future version. Use `frame.loc[string]` "
                        "instead.",
                        FutureWarning,
                    )
                    return self.loc[key]

            # error is raised from pandas
            meta = self._meta[_extract_meta(key)]
            dsk = partitionwise_graph(operator.getitem, name, self, key)
            graph = HighLevelGraph.from_collections(name, dsk, dependencies=[self])
            return new_dd_object(graph, name, meta, self.divisions)
        elif isinstance(key, slice):
            from pandas.api.types import is_float_dtype

            is_integer_slice = any(
                isinstance(i, Integral) for i in (key.start, key.step, key.stop)
            )
            # Slicing with integer labels is always iloc based except for a
            # float indexer for some reason
            if is_integer_slice and not is_float_dtype(self.index.dtype):
                # NOTE: this always fails currently, as iloc is mostly
                # unsupported, but we call it anyway here for future-proofing
                # and error-attribution purposes
                return self.iloc[key]
            else:
                return self.loc[key]

        if isinstance(key, (np.ndarray, list)) or (
            not is_dask_collection(key) and (is_series_like(key) or is_index_like(key))
        ):
            # error is raised from pandas
            meta = self._meta[_extract_meta(key)]

            dsk = partitionwise_graph(operator.getitem, name, self, key)
            graph = HighLevelGraph.from_collections(name, dsk, dependencies=[self])
            return new_dd_object(graph, name, meta, self.divisions)
        if isinstance(key, Series):
            # do not perform dummy calculation, as columns will not be changed.
            if self.divisions != key.divisions:
                from dask.dataframe.multi import _maybe_align_partitions

                self, key = _maybe_align_partitions([self, key])
            dsk = partitionwise_graph(operator.getitem, name, self, key)
            graph = HighLevelGraph.from_collections(name, dsk, dependencies=[self, key])
            return new_dd_object(graph, name, self, self.divisions)
        if isinstance(key, DataFrame):
            return self.where(key, np.nan)

        raise NotImplementedError(key)

    def __setitem__(self, key, value):
        if isinstance(key, (tuple, list)) and isinstance(value, DataFrame):
            df = self.assign(**{k: value[c] for k, c in zip(key, value.columns)})

        elif isinstance(key, pd.Index) and not isinstance(value, DataFrame):
            key = list(key)
            df = self.assign(**{k: value for k in key})
        elif (
            is_dataframe_like(key)
            or is_series_like(key)
            or isinstance(key, (DataFrame, Series))
        ):
            df = self.where(~key, value)
        elif not isinstance(key, str):
            raise NotImplementedError(f"Item assignment with {type(key)} not supported")
        else:
            df = self.assign(**{key: value})

        self.dask = df.dask
        self._name = df._name
        self._meta = df._meta
        self._divisions = df.divisions

    def __delitem__(self, key):
        result = self.drop([key], axis=1)
        self.dask = result.dask
        self._name = result._name
        self._meta = result._meta

    def __setattr__(self, key, value):
        try:
            columns = object.__getattribute__(self, "_meta").columns
        except AttributeError:
            columns = ()

        # exclude protected attributes from setitem
        if key in columns and key not in ["divisions", "dask", "_name", "_meta"]:
            self[key] = value
        else:
            object.__setattr__(self, key, value)

    def __getattr__(self, key):
        if key in self.columns:
            return self[key]
        elif key == "empty":
            # self.empty raises AttributeNotImplementedError, which includes
            # AttributeError, which means we end up here in self.__getattr__,
            # because DataFrame.__getattribute__ doesn't think the attribute exists
            # and uses __getattr__ as the fallback. So, get `self.empty` more
            # forcefully via object.__getattribute__ to raise informative error.
            object.__getattribute__(self, key)
        else:
            raise AttributeError("'DataFrame' object has no attribute %r" % key)

    def __dir__(self):
        o = set(dir(type(self)))
        o.update(self.__dict__)
        o.update(c for c in self.columns if (isinstance(c, str) and c.isidentifier()))
        return list(o)

    def __iter__(self):
        return iter(self._meta)

    def _ipython_key_completions_(self):
        return methods.tolist(self.columns)

    @property
    def ndim(self):
        """Return dimensionality"""
        return 2

    @property
    def shape(self):
        """
        Return a tuple representing the dimensionality of the DataFrame.

        The number of rows is a Delayed result. The number of columns
        is a concrete integer.

        Examples
        --------
        >>> df.size  # doctest: +SKIP
        (Delayed('int-07f06075-5ecc-4d77-817e-63c69a9188a8'), 2)
        """
        col_size = len(self.columns)
        if col_size == 0:
            return (self.index.shape[0], 0)
        row_size = delayed(int)(self.size / col_size)
        return (row_size, col_size)

    @property
    def dtypes(self):
        """Return data types"""
        return self._meta.dtypes

    @derived_from(pd.DataFrame)
    def select_dtypes(self, include=None, exclude=None):
        cs = self._meta.select_dtypes(include=include, exclude=exclude)
        indexer = self._get_columns_indexes_based_on_dtypes(cs)
        return self.iloc[:, indexer]

    def sort_values(
        self,
        by: str | list[str],
        npartitions: int | Literal["auto"] | None = None,
        ascending: bool = True,
        na_position: Literal["first"] | Literal["last"] = "last",
        sort_function: Callable[[pd.DataFrame], pd.DataFrame] | None = None,
        sort_function_kwargs: Mapping[str, Any] | None = None,
        **kwargs,
    ) -> DataFrame:
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
        from dask.dataframe.shuffle import sort_values

        return sort_values(
            self,
            by,
            ascending=ascending,
            npartitions=npartitions,
            na_position=na_position,
            sort_function=sort_function,
            sort_function_kwargs=sort_function_kwargs,
            **kwargs,
        )

    def set_index(
        self,
        other: str | Series,
        drop: bool = True,
        sorted: bool = False,
        npartitions: int | Literal["auto"] | None = None,
        divisions: Sequence | None = None,
        inplace: bool = False,
        sort: bool = True,
        **kwargs,
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

        Note that ``len(divisions)`` is equal to ``npartitions + 1``. This is because ``divisions``
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

        if inplace:
            raise NotImplementedError("The inplace= keyword is not supported")
        pre_sorted = sorted
        del sorted

        # Check other can be translated to column name or column object, possibly flattening it
        if not isinstance(other, str):
            # It may refer to several columns
            if isinstance(other, Sequence):  # type: ignore[unreachable]
                # Accept ["a"], but not [["a"]]
                if len(other) == 1 and (  # type: ignore[unreachable]
                    isinstance(other[0], str) or not isinstance(other[0], Sequence)
                ):
                    other = other[0]
                else:
                    raise NotImplementedError(
                        "Dask dataframe does not yet support multi-indexes.\n"
                        f"You tried to index with this index: {other}\n"
                        "Indexes must be single columns only."
                    )

            # Or be a frame directly
            elif isinstance(other, DataFrame):  # type: ignore[unreachable]
                raise NotImplementedError(
                    "Dask dataframe does not yet support multi-indexes.\n"
                    f"You tried to index with a frame with these columns: {list(other.columns)}\n"
                    "Indexes must be single columns only."
                )

        # If already a series
        if isinstance(other, Series):
            # If it's already the index, there's nothing to do
            if other._name == self.index._name:
                warnings.warn(
                    "New index has same name as existing, this is a no-op.", UserWarning
                )
                return self

        # If the name of a column/index
        else:
            # With the same name as the index, there's nothing to do either
            if other == self.index.name:
                warnings.warn(
                    "New index has same name as existing, this is a no-op.", UserWarning
                )
                return self

            # If a missing column, KeyError
            if other not in self.columns:
                raise KeyError(
                    f"Data has no column '{other}': use any column of {list(self.columns)}"
                )

        # Check divisions
        if divisions is not None:
            check_divisions(divisions)
        elif (
            isinstance(other, Index)
            and other.known_divisions
            and other.npartitions == self.npartitions
        ):
            # If the index has the same number of partitions and known
            # divisions, then we can treat it as pre-sorted with known
            # divisions
            pre_sorted = True
            divisions = other.divisions

        # If index is already sorted, take advantage of that with set_sorted_index
        if pre_sorted:
            from dask.dataframe.shuffle import set_sorted_index

            return set_sorted_index(
                self, other, drop=drop, divisions=divisions, **kwargs
            )
        else:
            from dask.dataframe.shuffle import set_index

            return set_index(
                self,
                other,
                drop=drop,
                npartitions=npartitions,
                divisions=divisions,
                sort=sort,
                **kwargs,
            )

    @derived_from(pd.DataFrame)
    def pop(self, item):
        out = self[item]
        del self[item]
        return out

    @derived_from(pd.DataFrame)
    def nlargest(self, n=5, columns=no_default, split_every=None):
        if columns is no_default:
            raise TypeError(
                "DataFrame.nlargest() missing required positional argument: 'columns'"
            )
        token = "dataframe-nlargest"
        return aca(
            self,
            chunk=M.nlargest,
            aggregate=M.nlargest,
            meta=self._meta,
            token=token,
            split_every=split_every,
            n=n,
            columns=columns,
        )

    @derived_from(pd.DataFrame)
    def nsmallest(self, n=5, columns=no_default, split_every=None):
        if columns is no_default:
            raise TypeError(
                "DataFrame.nsmallest() missing required positional argument: 'columns'"
            )
        token = "dataframe-nsmallest"
        return aca(
            self,
            chunk=M.nsmallest,
            aggregate=M.nsmallest,
            meta=self._meta,
            token=token,
            split_every=split_every,
            n=n,
            columns=columns,
        )

    @derived_from(pd.DataFrame)
    def groupby(
        self,
        by=None,
        group_keys=GROUP_KEYS_DEFAULT,
        sort=None,
        observed=None,
        dropna=None,
        **kwargs,
    ):
        from dask.dataframe.groupby import DataFrameGroupBy

        return DataFrameGroupBy(
            self,
            by=by,
            group_keys=group_keys,
            sort=sort,
            observed=observed,
            dropna=dropna,
            **kwargs,
        )

    @wraps(categorize)
    def categorize(self, columns=None, index=None, split_every=None, **kwargs):
        return categorize(
            self, columns=columns, index=index, split_every=split_every, **kwargs
        )

    @derived_from(pd.DataFrame)
    def assign(self, **kwargs):
        data = self.copy()
        for k, v in kwargs.items():
            if not (
                isinstance(v, Scalar)
                or is_series_like(v)
                or callable(v)
                or pd.api.types.is_scalar(v)
                or is_index_like(v)
                or isinstance(v, Array)
            ):
                raise TypeError(
                    f"Column assignment doesn't support type {typename(type(v))}"
                )
            if callable(v):
                kwargs[k] = v(data)
            if isinstance(v, Array):
                from dask.dataframe.io import from_dask_array

                if len(v.shape) > 1:
                    raise ValueError("Array assignment only supports 1-D arrays")
                if v.npartitions != data.npartitions:
                    raise ValueError(
                        "Number of partitions do not match "
                        f"({v.npartitions} != {data.npartitions})"
                    )
                kwargs[k] = from_dask_array(v, index=data.index, meta=data._meta)

            pairs = [k, kwargs[k]]

            # Figure out columns of the output
            df2 = data._meta_nonempty.assign(
                **_extract_meta({k: kwargs[k]}, nonempty=True)
            )
            data = elemwise(methods.assign, data, *pairs, meta=df2)

        return data

    @derived_from(pd.DataFrame, ua_args=["index"])
    def rename(self, index=None, columns=None):
        if index is not None:
            raise ValueError("Cannot rename index.")

        # *args here is index, columns but columns arg is already used
        return self.map_partitions(M.rename, None, columns=columns)

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
        return self.map_partitions(M.query, expr, **kwargs)

    @derived_from(pd.DataFrame)
    def eval(self, expr, inplace=None, **kwargs):
        if inplace is None:
            inplace = False
        else:
            warnings.warn(
                "`inplace` is deprecated and will be removed in a futuere version.",
                FutureWarning,
                2,
            )

        if "=" in expr and inplace in (True, None):
            raise NotImplementedError(
                "Inplace eval not supported. Please use inplace=False"
            )
        meta = self._meta.eval(expr, inplace=inplace, **kwargs)
        return self.map_partitions(M.eval, expr, meta=meta, inplace=inplace, **kwargs)

    @derived_from(pd.DataFrame)
    def dropna(self, how=no_default, subset=None, thresh=no_default):
        # These keywords are incompatible with each other.
        # Don't allow them both to be set.
        if how is not no_default and thresh is not no_default:
            raise TypeError(
                "You cannot set both the how and thresh arguments at the same time."
            )

        # Only specify `how` or `thresh` keyword if specified by the user
        # so we utilize other `dropna` keyword defaults appropriately
        kwargs = {"subset": subset}
        if how is not no_default:
            kwargs["how"] = how
        elif thresh is not no_default:
            kwargs["thresh"] = thresh

        return self.map_partitions(M.dropna, **kwargs, enforce_metadata=False)

    @_dummy_numpy_dispatcher("out")
    @derived_from(pd.DataFrame)
    def clip(self, lower=None, upper=None, axis=None):
        return self.map_partitions(
            M.clip,
            lower=lower,
            upper=upper,
            enforce_metadata=False,
            # See: https://github.com/rapidsai/cudf/issues/14369
            **({} if axis is None else {"axis": axis}),
        )

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

        elif axis not in [0, 1, None]:
            raise ValueError(f"No axis {axis} for object type {type(self)}")

    @derived_from(pd.DataFrame)
    def to_timestamp(self, freq=None, how="start", axis=0):
        df = elemwise(M.to_timestamp, self, freq, how, axis)
        df.divisions = tuple(pd.Index(self.divisions).to_timestamp(how=how, freq=freq))
        return df

    @derived_from(pd.DataFrame)
    def explode(self, column):
        meta = self._meta.explode(column)
        return self.map_partitions(M.explode, column, meta=meta, enforce_metadata=False)

    def to_bag(self, index=False, format="tuple"):
        """Convert to a dask Bag of tuples of each row.

        Parameters
        ----------
        index : bool, optional
            If True, the index is included as the first element of each tuple.
            Default is False.
        format : {"tuple", "dict", "frame"}, optional
            Whether to return a bag of tuples, dictionaries, or
            dataframe-like objects. Default is "tuple". If "frame",
            the original partitions of ``df`` will not be transformed
            in any way.
        """
        from dask.dataframe.io import to_bag

        return to_bag(self, index, format)

    def to_parquet(self, path, *args, **kwargs):
        """See dd.to_parquet docstring for more information"""
        from dask.dataframe.io import to_parquet

        return to_parquet(self, path, *args, **kwargs)

    def to_orc(self, path, *args, **kwargs):
        """See dd.to_orc docstring for more information"""
        from dask.dataframe.io import to_orc

        return to_orc(self, path, *args, **kwargs)

    @derived_from(pd.DataFrame)
    def to_string(self, max_rows=5):
        # option_context doesn't affect
        return self._repr_data().to_string(max_rows=max_rows, show_dimensions=False)

    def _get_numeric_data(self, how="any", subset=None):
        # calculate columns to avoid unnecessary calculation
        numerics = self._meta._get_numeric_data()

        if len(numerics.columns) < len(self.columns):
            name = self._token_prefix + "-get_numeric_data"
            return self.map_partitions(M._get_numeric_data, meta=numerics, token=name)
        else:
            # use myself if all numerics
            return self

    @classmethod
    def _validate_axis(cls, axis=0, none_is_zero: bool = True) -> None | Literal[0, 1]:
        if axis not in (0, 1, "index", "columns", None):
            raise ValueError(f"No axis named {axis}")
        # convert to numeric axis
        numeric_axis: dict[str | None, Literal[0, 1]] = {"index": 0, "columns": 1}
        if none_is_zero:
            numeric_axis[None] = 0

        return numeric_axis.get(axis, axis)

    @derived_from(pd.DataFrame)
    def drop(self, labels=None, axis=0, columns=None, errors="raise"):
        axis = self._validate_axis(axis)
        if axis == 0 and columns is not None:
            # Columns must be specified if axis==0
            return self.map_partitions(
                drop_by_shallow_copy, columns, errors=errors, enforce_metadata=False
            )
        elif axis == 1:
            return self.map_partitions(
                drop_by_shallow_copy, labels, errors=errors, enforce_metadata=False
            )
        raise NotImplementedError(
            "Drop currently only works for axis=1 or when columns is not None"
        )

    @_deprecated_kwarg("shuffle", "shuffle_method")
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
        npartitions=None,
        shuffle_method=None,
        broadcast=None,
    ):
        """Merge the DataFrame with another DataFrame

        This will merge the two datasets, either on the indices, a certain column
        in each dataset or the index in one dataset and the column in another.

        Parameters
        ----------
        right: dask.dataframe.DataFrame
        how : {'left', 'right', 'outer', 'inner'}, default: 'inner'
            How to handle the operation of the two objects:

            - left: use calling frame's index (or column if on is specified)
            - right: use other frame's index
            - outer: form union of calling frame's index (or column if on is
              specified) with other frame's index, and sort it
              lexicographically
            - inner: form intersection of calling frame's index (or column if
              on is specified) with other frame's index, preserving the order
              of the calling's one

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

        if not is_dataframe_like(right):
            raise ValueError("right must be DataFrame")

        from dask.dataframe.multi import merge

        return merge(
            self,
            right,
            how=how,
            on=on,
            left_on=left_on,
            right_on=right_on,
            left_index=left_index,
            right_index=right_index,
            suffixes=suffixes,
            npartitions=npartitions,
            indicator=indicator,
            shuffle_method=shuffle_method,
            broadcast=broadcast,
        )

    @_deprecated_kwarg("shuffle", "shuffle_method")
    @derived_from(pd.DataFrame)
    def join(
        self,
        other,
        on=None,
        how="left",
        lsuffix="",
        rsuffix="",
        npartitions=None,
        shuffle_method=None,
    ):
        if is_series_like(other) and hasattr(other, "name"):
            other = other.to_frame()

        if not is_dataframe_like(other):
            if not isinstance(other, list) or not all(
                [is_dataframe_like(o) for o in other]
            ):
                raise ValueError("other must be DataFrame or list of DataFrames")
            if how not in ["outer", "left"]:
                raise ValueError("merge_multi only supports left or outer joins")

            from dask.dataframe.multi import _recursive_pairwise_outer_join

            # If its an outer join we can use the full recursive pairwise join.
            if how == "outer":
                full = [self] + other

                return _recursive_pairwise_outer_join(
                    full,
                    on=on,
                    lsuffix=lsuffix,
                    rsuffix=rsuffix,
                    npartitions=npartitions,
                    shuffle_method=shuffle_method,
                )
            else:
                # Do recursive pairwise join on everything _except_ the last join
                # where we need to do a left join.
                other = _recursive_pairwise_outer_join(
                    other,
                    on=on,
                    lsuffix=lsuffix,
                    rsuffix=rsuffix,
                    npartitions=npartitions,
                    shuffle_method=shuffle_method,
                )

        from dask.dataframe.multi import merge

        return merge(
            self,
            other,
            how=how,
            left_index=on is None,
            right_index=True,
            left_on=on,
            suffixes=(lsuffix, rsuffix),
            npartitions=npartitions,
            shuffle_method=shuffle_method,
        )

    if not PANDAS_GE_200:

        @derived_from(pd.DataFrame)
        def append(self, other, interleave_partitions=False):
            if isinstance(other, Series):
                msg = (
                    "Unable to appending dd.Series to dd.DataFrame."
                    "Use pd.Series to append as row."
                )
                raise ValueError(msg)
            elif is_series_like(other):
                other = other.to_frame().T
            return super().append(other, interleave_partitions=interleave_partitions)

    @derived_from(pd.DataFrame)
    def iterrows(self):
        for i in range(self.npartitions):
            df = self.get_partition(i).compute()
            yield from df.iterrows()

    @derived_from(pd.DataFrame)
    def itertuples(self, index=True, name="Pandas"):
        for i in range(self.npartitions):
            df = self.get_partition(i).compute()
            yield from df.itertuples(index=index, name=name)

    @derived_from(pd.DataFrame)
    def items(self):
        for col_idx, label in enumerate(self.columns):
            yield label, self.iloc[:, col_idx]

    @classmethod
    def _bind_operator_method(cls, name, op, original=pd.DataFrame):
        """bind operator method like DataFrame.add to this class"""

        # name must be explicitly passed for div method whose name is truediv

        def meth(self, other, axis="columns", level=None, fill_value=None):
            if level is not None:
                raise NotImplementedError("level must be None")

            axis = self._validate_axis(axis)

            if axis in (1, "columns"):
                # When axis=1 and other is a series, `other` is transposed
                # and the operator is applied broadcast across rows. This
                # isn't supported with dd.Series.
                if isinstance(other, Series):
                    msg = f"Unable to {name} dd.Series with axis=1"
                    raise ValueError(msg)
                elif is_series_like(other):
                    # Special case for pd.Series to avoid unwanted partitioning
                    # of other. We pass it in as a kwarg to prevent this.
                    meta = _emulate(
                        op, self, other=other, axis=axis, fill_value=fill_value
                    )
                    return map_partitions(
                        op,
                        self,
                        other=other,
                        meta=meta,
                        axis=axis,
                        fill_value=fill_value,
                        enforce_metadata=False,
                    )

            meta = _emulate(op, self, other, axis=axis, fill_value=fill_value)
            return map_partitions(
                op,
                self,
                other,
                meta=meta,
                axis=axis,
                fill_value=fill_value,
                enforce_metadata=False,
            )

        meth.__name__ = name
        setattr(cls, name, derived_from(original)(meth))

    @classmethod
    def _bind_comparison_method(cls, name, comparison, original=pd.DataFrame):
        """bind comparison method like DataFrame.eq to this class"""

        def meth(self, other, axis="columns", level=None):
            if level is not None:
                raise NotImplementedError("level must be None")
            axis = self._validate_axis(axis)
            return elemwise(comparison, self, other, axis=axis)

        meth.__name__ = name
        setattr(cls, name, derived_from(original)(meth))

    @insert_meta_param_description(pad=12)
    def apply(
        self,
        func,
        axis=0,
        broadcast=None,
        raw=False,
        reduce=None,
        args=(),
        meta=no_default,
        result_type=None,
        **kwds,
    ):
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
        dask.DataFrame.map_partitions
        """

        if broadcast is not None:
            warnings.warn(
                "The `broadcast` argument is no longer used/supported. "
                "It will be dropped in a future release.",
                category=FutureWarning,
            )

        axis = self._validate_axis(axis)
        pandas_kwargs = {"axis": axis, "raw": raw, "result_type": result_type}

        kwds.update(pandas_kwargs)

        if axis == 0:
            msg = (
                "dd.DataFrame.apply only supports axis=1\n"
                "  Try: df.apply(func, axis=1)"
            )
            raise NotImplementedError(msg)

        if meta is no_default:
            meta = _emulate(
                M.apply, self._meta_nonempty, func, args=args, udf=True, **kwds
            )
            warnings.warn(meta_warning(meta))
        kwds.update({"parent_meta": self._meta})
        return map_partitions(M.apply, self, func, args=args, meta=meta, **kwds)

    @derived_from(pd.DataFrame)
    def applymap(self, func, meta=no_default):
        # Let pandas raise deprecation warnings
        self._meta.applymap(func)
        return elemwise(methods.applymap, self, func, meta=meta)

    def map(self, func, meta=no_default, na_action=None):
        if not PANDAS_GE_210:
            raise NotImplementedError(
                f"DataFrame.map requires pandas>=2.1.0, but pandas={PANDAS_VERSION} is "
                "installed."
            )
        return elemwise(M.map, self, func, meta=meta, na_action=na_action)

    @derived_from(pd.DataFrame)
    def round(self, decimals=0):
        return elemwise(M.round, self, decimals)

    @derived_from(pd.DataFrame)
    def nunique(self, split_every=False, dropna=True, axis=0):
        if axis == 1:
            # split_every not used for axis=1
            meta = self._meta_nonempty.nunique(axis=axis)
            return self.map_partitions(
                M.nunique,
                meta=meta,
                token="series-nunique",
                axis=axis,
                dropna=dropna,
                enforce_metadata=False,
            )
        else:
            nunique_list = [
                self.iloc[:, i].nunique(split_every=split_every, dropna=dropna)
                for i in range(len(self.columns))
            ]
            name = "series-" + tokenize(*nunique_list)
            dsk = {
                (name, 0): (
                    apply,
                    pd.Series,
                    [[(s._name, 0) for s in nunique_list]],
                    {"index": self.columns},
                )
            }
            graph = HighLevelGraph.from_collections(
                name, dsk, dependencies=nunique_list
            )
            return Series(graph, name, self._meta.nunique(), (None, None))

    @derived_from(pd.DataFrame)
    def mode(self, dropna=True, split_every=False, numeric_only=False):
        mode_series_list = []
        for col_index in range(len(self.columns)):
            col_series = self.iloc[:, col_index]
            if numeric_only and not pd.api.types.is_numeric_dtype(col_series.dtype):
                continue
            mode_series = Series.mode(
                col_series, dropna=dropna, split_every=split_every
            )
            mode_series_list.append(mode_series)

        name = "concat-" + tokenize(*mode_series_list)

        dsk = {
            (name, 0): (
                apply,
                methods.concat,
                [[(df._name, 0) for df in mode_series_list]],
                {"axis": 1},
            )
        }

        meta = methods.concat([df._meta for df in mode_series_list], axis=1)
        graph = HighLevelGraph.from_collections(
            name, dsk, dependencies=mode_series_list
        )
        ddf = new_dd_object(graph, name, meta, divisions=(None, None))

        return ddf

    @derived_from(pd.DataFrame)
    def cov(self, min_periods=None, numeric_only=no_default, split_every=False):
        return _cov_corr(
            self, min_periods, numeric_only=numeric_only, split_every=split_every
        )

    @derived_from(pd.DataFrame)
    def corr(
        self,
        method="pearson",
        min_periods=None,
        numeric_only=no_default,
        split_every=False,
    ):
        if method != "pearson":
            raise NotImplementedError("Only Pearson correlation has been implemented")

        return _cov_corr(
            self, min_periods, True, numeric_only=numeric_only, split_every=split_every
        )

    def info(self, buf=None, verbose=False, memory_usage=False):
        """
        Concise summary of a Dask DataFrame.
        """

        if buf is None:
            import sys

            buf = sys.stdout

        lines = [str(type(self))]

        if len(self.columns) == 0:
            lines.append(f"{type(self.index._meta).__name__}: 0 entries")
            lines.append(f"Empty {type(self).__name__}")
            if PANDAS_GE_150:
                # pandas dataframe started adding a newline when info is called.
                lines.append("")
            put_lines(buf, lines)
            return

        # Group and execute the required computations
        computations = {}
        if verbose:
            memory_usage = True
            computations.update({"index": self.index, "count": self.count()})
        if memory_usage:
            computations.update(
                {"memory_usage": self.map_partitions(M.memory_usage, index=True)}
            )
        computations = dict(
            zip(computations.keys(), da.compute(*computations.values()))
        )

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
                    zip(self.columns, counts, self.dtypes)
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
    def memory_usage(self, index=True, deep=False):
        result = self.map_partitions(M.memory_usage, index=index, deep=deep)
        result = result.groupby(result.index).sum()
        return result

    def pivot_table(self, index=None, columns=None, values=None, aggfunc="mean"):
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
        from dask.dataframe.reshape import pivot_table

        return pivot_table(
            self, index=index, columns=columns, values=values, aggfunc=aggfunc
        )

    def melt(
        self,
        id_vars=None,
        value_vars=None,
        var_name=None,
        value_name="value",
        col_level=None,
    ):
        """
        Unpivots a DataFrame from wide format to long format,
        optionally leaving identifier variables set.

        This function is useful to massage a DataFrame into a format where
        one or more columns are identifier variables (``id_vars``), while
        all other columns, considered measured variables (``value_vars``),
        are "unpivoted" to the row axis, leaving just two non-identifier
        columns, 'variable' and 'value'.

        Parameters
        ----------
        frame : DataFrame
        id_vars : tuple, list, or ndarray, optional
            Column(s) to use as identifier variables.
        value_vars : tuple, list, or ndarray, optional
            Column(s) to unpivot. If not specified, uses all columns that
            are not set as `id_vars`.
        var_name : scalar
            Name to use for the 'variable' column. If None it uses
            ``frame.columns.name`` or 'variable'.
        value_name : scalar, default 'value'
            Name to use for the 'value' column.
        col_level : int or string, optional
            If columns are a MultiIndex then use this level to melt.

        Returns
        -------
        DataFrame
            Unpivoted DataFrame.

        See Also
        --------
        pandas.DataFrame.melt
        """
        from dask.dataframe.reshape import melt

        return melt(
            self,
            id_vars=id_vars,
            value_vars=value_vars,
            var_name=var_name,
            value_name=value_name,
            col_level=col_level,
        )

    def to_records(self, index=False, lengths=None):
        from dask.dataframe.io import to_records

        if lengths is True:
            lengths = tuple(self.map_partitions(len).compute())

        records = to_records(self)

        chunks = self._validate_chunks(records, lengths)
        records._chunks = (chunks[0],)

        return records

    @derived_from(pd.DataFrame)
    def to_html(self, max_rows=5):
        # pd.Series doesn't have html repr
        data = self._repr_data().to_html(max_rows=max_rows, show_dimensions=False)
        return get_template("dataframe.html.j2").render(
            data=data,
            name=self._name,
            layers=maybe_pluralize(len(self.dask.layers), "graph layer"),
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

    def _repr_html_(self):
        data = self._repr_data().to_html(
            max_rows=5, show_dimensions=False, notebook=True
        )
        return get_template("dataframe.html.j2").render(
            data=data,
            name=self._name,
            layers=maybe_pluralize(len(self.dask.layers), "graph layer"),
        )

    def _select_columns_or_index(self, columns_or_index):
        """
        Parameters
        ----------
        columns_or_index
            Column or index name, or a list of these

        Returns
        -------
        dd.DataFrame
            Dask DataFrame with columns corresponding to each column or
            index level in columns_or_index.  If included, the column
            corresponding to the index level is named _index
        """

        # Ensure columns_or_index is a list
        columns_or_index = (
            columns_or_index
            if isinstance(columns_or_index, list)
            else [columns_or_index]
        )

        column_names = [
            n for n in columns_or_index if self._is_column_label_reference(n)
        ]

        selected_df = self[column_names]
        if self._contains_index_name(columns_or_index):
            # Index name was included
            selected_df = selected_df.assign(_index=self.index)

        return selected_df

    def _is_column_label_reference(self, key):
        """
        Test whether a key is a column label reference

        To be considered a column label reference, `key` must match the name of at
        least one column.
        """
        return (
            not is_dask_collection(key)
            and (np.isscalar(key) or isinstance(key, tuple))
            and key in self.columns
        )

    @classmethod
    def from_dict(
        cls, data, *, npartitions, orient="columns", dtype=None, columns=None
    ):
        """
        Construct a Dask DataFrame from a Python Dictionary

        See Also
        --------
        dask.dataframe.from_dict
        """
        from dask.dataframe.io import from_dict

        return from_dict(
            data,
            npartitions,
            orient=orient,
            dtype=dtype,
            columns=columns,
            constructor=cls._partition_type,
        )


# bind operators
# TODO: dynamically bound operators are defeating type annotations
for op in [
    operator.abs,
    operator.add,
    operator.and_,
    operator.eq,
    operator.gt,
    operator.ge,
    operator.inv,
    operator.lt,
    operator.le,
    operator.mod,
    operator.mul,
    operator.ne,
    operator.neg,
    operator.or_,
    operator.pow,
    operator.sub,
    operator.truediv,
    operator.floordiv,
    operator.xor,
]:
    _Frame._bind_operator(op)
    Scalar._bind_operator(op)

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
    meth = getattr(pd.DataFrame, name)
    DataFrame._bind_operator_method(name, meth)

    meth = getattr(pd.Series, name)
    Series._bind_operator_method(name, meth)

for name in ["lt", "gt", "le", "ge", "ne", "eq"]:
    meth = getattr(pd.DataFrame, name)
    DataFrame._bind_comparison_method(name, meth)

    meth = getattr(pd.Series, name)
    Series._bind_comparison_method(name, meth)


def is_broadcastable(dfs, s):
    """
    This Series is broadcastable against another dataframe in the sequence
    """

    def compare(s, df):
        try:
            return s.divisions == (df.columns.min(), df.columns.max())
        except TypeError:
            return False

    return (
        isinstance(s, Series)
        and s.npartitions == 1
        and s.known_divisions
        and any(compare(s, df) for df in dfs if isinstance(df, DataFrame))
    )


def elemwise(op, *args, meta=no_default, out=None, transform_divisions=True, **kwargs):
    """Elementwise operation for Dask dataframes

    Parameters
    ----------
    op: callable
        Function to apply across input dataframes
    *args: DataFrames, Series, Scalars, Arrays, etc.
        The arguments of the operation
    meta: pd.DataFrame, pd.Series (optional)
        Valid metadata for the operation.  Will evaluate on a small piece of
        data if not provided.
    transform_divisions: boolean
        If the input is a ``dask.dataframe.Index`` we normally will also apply
        the function onto the divisions and apply those transformed divisions
        to the output.  You can pass ``transform_divisions=False`` to override
        this behavior
    out : dask.DataFrame, dask.Series, dask.Scalar, or None
        If out is a dask.DataFrame, dask.Series or dask.Scalar then
        this overwrites the contents of it with the result
    **kwargs: scalars

    Examples
    --------
    >>> elemwise(operator.add, df.x, df.y)  # doctest: +SKIP
    """
    _name = funcname(op) + "-" + tokenize(op, *args, **kwargs)

    args = _maybe_from_pandas(args)

    from dask.dataframe.multi import _maybe_align_partitions

    args = _maybe_align_partitions(args)
    dasks = [arg for arg in args if isinstance(arg, (_Frame, Scalar, Array))]
    dfs = [df for df in dasks if isinstance(df, _Frame)]

    # Clean up dask arrays if present
    deps = dasks.copy()
    for i, a in enumerate(dasks):
        if not isinstance(a, Array):
            continue
        # Ensure that they have similar-ish chunk structure
        if not all(not a.chunks or len(a.chunks[0]) == df.npartitions for df in dfs):
            msg = (
                "When combining dask arrays with dataframes they must "
                "match chunking exactly.  Operation: %s" % funcname(op)
            )
            raise ValueError(msg)
        # Rechunk to have a single chunk along all other axes
        if a.ndim > 1:
            a = a.rechunk({i + 1: d for i, d in enumerate(a.shape[1:])})
            dasks[i] = a

    divisions = dfs[0].divisions
    if transform_divisions and isinstance(dfs[0], Index) and len(dfs) == 1:
        try:
            divisions = op(
                *[pd.Index(arg.divisions) if arg is dfs[0] else arg for arg in args],
                **kwargs,
            )
            if isinstance(divisions, pd.Index):
                divisions = methods.tolist(divisions)
        except Exception:
            pass
        else:
            if not valid_divisions(divisions):
                divisions = [None] * (dfs[0].npartitions + 1)

    _is_broadcastable = partial(is_broadcastable, dfs)
    dfs = list(remove(_is_broadcastable, dfs))

    other = [
        (i, arg)
        for i, arg in enumerate(args)
        if not isinstance(arg, (_Frame, Scalar, Array))
    ]

    # adjust the key length of Scalar
    dsk = partitionwise_graph(op, _name, *args, **kwargs)

    graph = HighLevelGraph.from_collections(_name, dsk, dependencies=deps)

    if meta is no_default:
        if len(dfs) >= 2 and not all(hasattr(d, "npartitions") for d in dasks):
            # should not occur in current funcs
            msg = "elemwise with 2 or more DataFrames and Scalar is not supported"
            raise NotImplementedError(msg)
        # For broadcastable series, use no rows.
        parts = [
            d._meta
            if _is_broadcastable(d)
            else np.empty((), dtype=d.dtype)
            if isinstance(d, Array)
            else d._meta_nonempty
            for d in dasks
        ]
        with raise_on_meta_error(funcname(op)):
            meta = partial_by_order(*parts, function=op, other=other)

    result = new_dd_object(graph, _name, meta, divisions)
    return handle_out(out, result)


def handle_out(out, result):
    """Handle out parameters

    If out is a dask.DataFrame, dask.Series or dask.Scalar then
    this overwrites the contents of it with the result
    """
    if isinstance(out, tuple):
        if len(out) == 1:
            out = out[0]
        elif len(out) > 1:
            raise NotImplementedError("The out parameter is not fully supported")
        else:
            out = None

    # Notice, we use .__class__ as opposed to type() in order to support
    # object proxies see <https://github.com/dask/dask/pull/6981>
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
        out._meta = result._meta
        out._name = result._name
        out.dask = result.dask

        if not isinstance(out, Scalar):
            out._divisions = result.divisions
        return result
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


def _maybe_from_pandas(dfs):
    from dask.dataframe.io import from_pandas

    dfs = [
        from_pandas(df, 1)
        if (is_series_like(df) or is_dataframe_like(df)) and not is_dask_collection(df)
        else df
        for df in dfs
    ]
    return dfs


def hash_shard(
    df, nparts, split_out_setup=None, split_out_setup_kwargs=None, ignore_index=False
):
    if split_out_setup:
        h = split_out_setup(df, **(split_out_setup_kwargs or {}))
    else:
        h = df

    h = hash_object_dispatch(h, index=False)
    return group_split_dispatch(df, h % nparts, nparts, ignore_index=ignore_index)


def split_evenly(df, k):
    """Split dataframe into k roughly equal parts"""
    divisions = np.linspace(0, len(df), k + 1).astype(int)
    return {i: df.iloc[divisions[i] : divisions[i + 1]] for i in range(k)}


def split_out_on_index(df):
    h = df.index
    if isinstance(h, pd.MultiIndex):
        h = meta_frame_constructor(df)([], index=h).reset_index()
    return h


def split_out_on_cols(df, cols=None):
    return df[cols]


@insert_meta_param_description
def apply_concat_apply(
    args,
    chunk=None,
    aggregate=None,
    combine=None,
    meta=no_default,
    token=None,
    chunk_kwargs=None,
    aggregate_kwargs=None,
    combine_kwargs=None,
    split_every=None,
    split_out=None,
    split_out_setup=None,
    split_out_setup_kwargs=None,
    sort=None,
    ignore_index=False,
    **kwargs,
):
    """Apply a function to blocks, then concat, then apply again

    Parameters
    ----------
    args :
        Positional arguments for the `chunk` function. All `dask.dataframe`
        objects should be partitioned and indexed equivalently.
    chunk : function [block-per-arg] -> block
        Function to operate on each block of data
    aggregate : function concatenated-block -> block
        Function to operate on the concatenated result of chunk
    combine : function concatenated-block -> block, optional
        Function to operate on intermediate concatenated results of chunk
        in a tree-reduction. If not provided, defaults to aggregate.
    $META
    token : str, optional
        The name to use for the output keys.
    chunk_kwargs : dict, optional
        Keywords for the chunk function only.
    aggregate_kwargs : dict, optional
        Keywords for the aggregate function only.
    combine_kwargs : dict, optional
        Keywords for the combine function only.
    split_every : int, optional
        Group partitions into groups of this size while performing a
        tree-reduction. If set to False, no tree-reduction will be used,
        and all intermediates will be concatenated and passed to ``aggregate``.
        Default is 8.
    split_out : int, optional
        Number of output partitions. Split occurs after first chunk reduction.
    split_out_setup : callable, optional
        If provided, this function is called on each chunk before performing
        the hash-split. It should return a pandas object, where each row
        (excluding the index) is hashed. If not provided, the chunk is hashed
        as is.
    split_out_setup_kwargs : dict, optional
        Keywords for the `split_out_setup` function only.
    sort : bool, default None
        If allowed, sort the keys of the output aggregation.
    ignore_index : bool, default False
        If True, do not preserve index values throughout ACA operations.
    kwargs :
        All remaining keywords will be passed to ``chunk``, ``aggregate``, and
        ``combine``.

    Examples
    --------
    >>> def chunk(a_block, b_block):
    ...     pass

    >>> def agg(df):
    ...     pass

    >>> apply_concat_apply([a, b], chunk=chunk, aggregate=agg)  # doctest: +SKIP
    """
    if split_out is None:
        split_out = 1
    if chunk_kwargs is None:
        chunk_kwargs = dict()
    if aggregate_kwargs is None:
        aggregate_kwargs = dict()
    chunk_kwargs.update(kwargs)
    aggregate_kwargs.update(kwargs)

    if combine is None:
        if combine_kwargs:
            raise ValueError("`combine_kwargs` provided with no `combine`")
        combine = aggregate
        combine_kwargs = aggregate_kwargs
    else:
        if combine_kwargs is None:
            combine_kwargs = dict()
        combine_kwargs.update(kwargs)

    if not isinstance(args, (tuple, list)):
        args = [args]

    dfs = [arg for arg in args if isinstance(arg, _Frame)]

    npartitions = {arg.npartitions for arg in dfs}
    if len(npartitions) > 1:
        raise ValueError("All arguments must have same number of partitions")
    npartitions = npartitions.pop()

    if split_every is None:
        split_every = 8
    elif split_every is False:
        split_every = npartitions
    elif split_every < 2 or not isinstance(split_every, Integral):
        raise ValueError("split_every must be an integer >= 2")

    token_key = tokenize(
        token or (chunk, aggregate),
        meta,
        args,
        chunk_kwargs,
        aggregate_kwargs,
        combine_kwargs,
        split_every,
        split_out,
        split_out_setup,
        split_out_setup_kwargs,
    )

    # Blockwise Chunk Layer
    chunk_name = f"{token or funcname(chunk)}-chunk-{token_key}"
    chunked = map_bag_partitions(
        chunk,
        # Convert _Frame collections to Bag
        *[
            arg.to_bag(format="frame") if isinstance(arg, _Frame) else arg
            for arg in args
        ],
        token=chunk_name,
        **chunk_kwargs,
    )

    # NOTE: `chunked` is now a Bag collection.
    # We don't use a DataFrame collection, because
    # the partitions may not contain dataframe-like
    # objects anymore.

    # Blockwise Split Layer
    if split_out and split_out > 1:
        chunked = chunked.map_partitions(
            hash_shard,
            split_out,
            split_out_setup,
            split_out_setup_kwargs,
            ignore_index,
            token="split-%s" % token_key,
        )

    # Handle sort behavior
    if sort is not None:
        if sort and split_out > 1:
            raise NotImplementedError(
                "Cannot guarantee sorted keys for `split_out>1`."
                " Try using split_out=1, or grouping with sort=False."
            )
        aggregate_kwargs = aggregate_kwargs or {}
        aggregate_kwargs["sort"] = sort

    # Tree-Reduction Layer
    final_name = f"{token or funcname(aggregate)}-agg-{token_key}"
    layer = DataFrameTreeReduction(
        final_name,
        chunked.name,
        npartitions,
        partial(_concat, ignore_index=ignore_index),
        partial(combine, **combine_kwargs) if combine_kwargs else combine,
        finalize_func=partial(aggregate, **aggregate_kwargs)
        if aggregate_kwargs
        else aggregate,
        split_every=split_every,
        split_out=split_out if (split_out and split_out > 1) else None,
        tree_node_name=f"{token or funcname(combine)}-combine-{token_key}",
    )

    if meta is no_default:
        meta_chunk = _emulate(chunk, *args, udf=True, **chunk_kwargs)
        meta = _emulate(
            aggregate, _concat([meta_chunk], ignore_index), udf=True, **aggregate_kwargs
        )
    meta = make_meta(
        meta,
        index=(getattr(make_meta(dfs[0]), "index", None) if dfs else None),
        parent_meta=dfs[0]._meta,
    )

    graph = HighLevelGraph.from_collections(final_name, layer, dependencies=(chunked,))
    divisions = [None] * ((split_out or 1) + 1)
    return new_dd_object(graph, final_name, meta, divisions, parent_meta=dfs[0]._meta)


aca = apply_concat_apply


def _extract_meta(x, nonempty=False):
    """
    Extract internal cache data (``_meta``) from dd.DataFrame / dd.Series
    """
    if isinstance(x, (Scalar, _Frame)):
        return x._meta_nonempty if nonempty else x._meta
    elif isinstance(x, list):
        return [_extract_meta(_x, nonempty) for _x in x]
    elif isinstance(x, tuple):
        return tuple(_extract_meta(_x, nonempty) for _x in x)
    elif isinstance(x, dict):
        res = {}
        for k in x:
            res[k] = _extract_meta(x[k], nonempty)
        return res
    elif isinstance(x, Delayed):
        raise ValueError(
            "Cannot infer dataframe metadata with a `dask.delayed` argument"
        )
    else:
        return x


def _emulate(func, *args, udf=False, **kwargs):
    """
    Apply a function using args / kwargs. If arguments contain dd.DataFrame /
    dd.Series, using internal cache (``_meta``) for calculation
    """
    with raise_on_meta_error(funcname(func), udf=udf), check_numeric_only_deprecation():
        return func(*_extract_meta(args, True), **_extract_meta(kwargs, True))


@insert_meta_param_description
def map_partitions(
    func,
    *args,
    meta=no_default,
    enforce_metadata=True,
    transform_divisions=True,
    align_dataframes=True,
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
    name = kwargs.pop("token", None)
    parent_meta = kwargs.pop("parent_meta", None)

    assert callable(func)
    if name is not None:
        token = tokenize(meta, *args, **kwargs)
    else:
        name = funcname(func)
        token = tokenize(func, meta, *args, **kwargs)
    name = f"{name}-{token}"

    from dask.dataframe.multi import _maybe_align_partitions

    if align_dataframes:
        args = _maybe_from_pandas(args)
        try:
            args = _maybe_align_partitions(args)
        except ValueError as e:
            raise ValueError(
                f"{e}. If you don't want the partitions to be aligned, and are "
                "calling `map_partitions` directly, pass `align_dataframes=False`."
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
    for arg in args:
        if isinstance(arg, _Frame):
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

    divisions = _get_divisions_map_partitions(
        align_dataframes, transform_divisions, dfs, func, args, kwargs
    )

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
            *args2,
            dependencies=dependencies,
            _func=func,
            _meta=meta,
            **kwargs3,
        )
    else:
        kwargs4 = kwargs if simple else kwargs3
        dsk = partitionwise_graph(
            func, name, *args2, **kwargs4, dependencies=dependencies
        )

    graph = HighLevelGraph.from_collections(name, dsk, dependencies=dependencies)
    return new_dd_object(graph, name, meta, divisions)


def _get_divisions_map_partitions(
    align_dataframes, transform_divisions, dfs, func, args, kwargs
):
    """
    Helper to get divisions for map_partitions and map_overlap output.
    """
    if align_dataframes:
        divisions = dfs[0].divisions
    else:
        # Unaligned, dfs is a mix of 1 partition and 1+ partition dataframes,
        # use longest divisions found
        divisions = max((d.divisions for d in dfs), key=len)
    if transform_divisions and isinstance(dfs[0], Index) and len(dfs) == 1:
        try:
            divisions = func(
                *[pd.Index(a.divisions) if a is dfs[0] else a for a in args], **kwargs
            )
            if isinstance(divisions, pd.Index):
                divisions = methods.tolist(divisions)
        except Exception:
            pass
        else:
            if not valid_divisions(divisions):
                divisions = [None] * (dfs[0].npartitions + 1)
    return divisions


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
        meta = _emulate(func, *args, udf=True, **kwargs)
        meta_is_emulated = True
    else:
        meta = make_meta(meta, index=meta_index, parent_meta=parent_meta)
        meta_is_emulated = False

    if not (has_parallel_type(meta) or is_arraylike(meta) and meta.shape) and not all(
        isinstance(arg, Scalar) for arg in args
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


def apply_and_enforce(*args, **kwargs):
    """Apply a function, and enforce the output to match meta

    Ensures the output has the same columns, even if empty."""
    func = kwargs.pop("_func")
    meta = kwargs.pop("_meta")
    df = func(*args, **kwargs)
    if is_dataframe_like(df) or is_series_like(df) or is_index_like(df):
        if not len(df):
            return meta
        if is_dataframe_like(df):
            check_matching_columns(meta, df)
            c = meta.columns
        else:
            c = meta.name
        return _rename(c, df)
    return df


def _rename(columns, df):
    """
    Rename columns of pd.DataFrame or name of pd.Series.
    Not for dd.DataFrame or dd.Series.

    Parameters
    ----------
    columns : tuple, string, pd.DataFrame or pd.Series
        Column names, Series name or pandas instance which has the
        target column names / name.
    df : pd.DataFrame or pd.Series
        target DataFrame / Series to be renamed
    """
    assert not isinstance(df, _Frame)

    if columns is no_default:
        return df

    if isinstance(columns, Iterator):
        columns = list(columns)

    if is_dataframe_like(df):
        if is_dataframe_like(columns):
            columns = columns.columns
        if not isinstance(columns, pd.Index):
            columns = pd.Index(columns)
        if (
            len(columns) == len(df.columns)
            and type(columns) is type(df.columns)
            and columns.dtype == df.columns.dtype
            and columns.equals(df.columns)
        ):
            # if target is identical, rename is not necessary
            return df
        # deep=False doesn't doesn't copy any data/indices, so this is cheap
        df = df.copy(deep=False)
        df.columns = columns
        return df
    elif is_series_like(df) or is_index_like(df):
        if is_series_like(columns) or is_index_like(columns):
            columns = columns.name
        if df.name == columns:
            return df
        return df.rename(columns)
    # map_partition may pass other types
    return df


def _rename_dask(df, names):
    """
    Destructively rename columns of dd.DataFrame or name of dd.Series.
    Not for pd.DataFrame or pd.Series.

    Internally used to overwrite dd.DataFrame.columns and dd.Series.name
    We can't use map_partition because it applies function then rename

    Parameters
    ----------
    df : dd.DataFrame or dd.Series
        target DataFrame / Series to be renamed
    names : tuple, string
        Column names/Series name
    """

    assert isinstance(df, _Frame)
    metadata = _rename(names, df._meta)
    name = f"rename-{tokenize(df, metadata)}"

    dsk = partitionwise_graph(_rename, name, metadata, df)
    graph = HighLevelGraph.from_collections(name, dsk, dependencies=[df])
    return new_dd_object(graph, name, metadata, df.divisions)


def quantile(df, q, method="default"):
    """Approximate quantiles of Series.

    Parameters
    ----------
    q : list/array of floats
        Iterable of numbers ranging from 0 to 100 for the desired quantiles
    method : {'default', 'tdigest', 'dask'}, optional
        What method to use. By default will use dask's internal custom
        algorithm (``'dask'``).  If set to ``'tdigest'`` will use tdigest for
        floats and ints and fallback to the ``'dask'`` otherwise.
    """
    # current implementation needs q to be sorted so
    # sort if array-like, otherwise leave it alone
    q_ndarray = np.array(q)
    if q_ndarray.ndim > 0:
        q_ndarray.sort(kind="mergesort")
        q = q_ndarray

    assert isinstance(df, Series)

    allowed_methods = ["default", "dask", "tdigest"]
    if method not in allowed_methods:
        raise ValueError("method can only be 'default', 'dask' or 'tdigest'")

    if method == "default":
        internal_method = "dask"
    else:
        internal_method = method

    # currently, only Series has quantile method
    if isinstance(df, Index):
        series_typ = df._meta.to_series()._constructor
        meta = df._meta_nonempty.to_series().quantile(q)
    else:
        if is_series_like(df._meta):
            series_typ = df._meta._constructor
        else:
            series_typ = df._meta._constructor_sliced
        meta = df._meta_nonempty.quantile(q)

    if is_series_like(meta):
        # Index.quantile(list-like) must be pd.Series, not pd.Index
        df_name = df.name
        finalize_tsk = lambda tsk: (series_typ, tsk, q, None, df_name)
        return_type = Series
    else:
        # repack as Series object to convert scalars from NumPy back before extracting
        finalize_tsk = lambda tsk: (
            lambda *args, **kwargs: series_typ(*args, **kwargs)[0],
            tsk,
            [0],
        )
        return_type = Scalar
        q = [q]

    # pandas uses quantile in [0, 1]
    # numpy / everyone else uses [0, 100]
    qs = np.asarray(q) * 100
    token = tokenize(df, qs)

    if len(qs) == 0:
        name = "quantiles-" + token
        empty_index = pd.Index([], dtype=float)

        return Series(
            {(name, 0): series_typ([], name=df.name, index=empty_index, dtype="float")},
            name,
            df._meta,
            [None, None],
        )
    else:
        new_divisions = [np.min(q), np.max(q)]

    df = df.dropna()

    if internal_method == "tdigest" and (
        np.issubdtype(df.dtype, np.floating) or np.issubdtype(df.dtype, np.integer)
    ):
        from dask.utils import import_required

        import_required(
            "crick", "crick is a required dependency for using the t-digest method."
        )

        from dask.array.percentile import _percentiles_from_tdigest, _tdigest_chunk

        name = "quantiles_tdigest-1-" + token
        val_dsk = {
            (name, i): (_tdigest_chunk, (getattr, key, "values"))
            for i, key in enumerate(df.__dask_keys__())
        }

        name2 = "quantiles_tdigest-2-" + token
        merge_dsk = {
            (name2, 0): finalize_tsk((_percentiles_from_tdigest, qs, sorted(val_dsk)))
        }
    else:
        from dask.array.dispatch import percentile_lookup as _percentile
        from dask.array.percentile import merge_percentiles

        # Add 0 and 100 during calculation for more robust behavior (hopefully)
        calc_qs = np.pad(qs, 1, mode="constant")
        calc_qs[-1] = 100
        name = "quantiles-1-" + token
        val_dsk = {
            (name, i): (_percentile, key, calc_qs)
            for i, key in enumerate(df.__dask_keys__())
        }

        name2 = "quantiles-2-" + token
        merge_dsk = {
            (name2, 0): finalize_tsk(
                (
                    merge_percentiles,
                    qs,
                    [calc_qs] * df.npartitions,
                    sorted(val_dsk),
                    "lower",
                    None,
                    False,
                )
            )
        }
    dsk = merge(val_dsk, merge_dsk)
    graph = HighLevelGraph.from_collections(name2, dsk, dependencies=[df])
    return return_type(graph, name2, meta, new_divisions)


def _cov_corr(
    df,
    min_periods=None,
    corr=False,
    scalar=False,
    numeric_only=no_default,
    split_every=False,
):
    """DataFrame covariance and pearson correlation.

    Computes pairwise covariance or correlation of columns, excluding NA/null
    values.

    Parameters
    ----------
    df : DataFrame
    min_periods : int, optional
        Minimum number of observations required per pair of columns
        to have a valid result.
    corr : bool, optional
        If True, compute the Pearson correlation. If False [default], compute
        the covariance.
    scalar : bool, optional
        If True, compute covariance between two variables as a scalar. Only
        valid if `df` has 2 columns.  If False [default], compute the entire
        covariance/correlation matrix.
    split_every : int, optional
        Group partitions into groups of this size while performing a
        tree-reduction. If set to False, no tree-reduction will be used.
        Default is False.
    """
    if min_periods is None:
        min_periods = 2
    elif min_periods < 2:
        raise ValueError("min_periods must be >= 2")

    if split_every is False:
        split_every = df.npartitions
    elif split_every < 2 or not isinstance(split_every, Integral):
        raise ValueError("split_every must be an integer >= 2")

    # Handle selecting numeric data and associated deprecation warning
    maybe_warn = False
    if numeric_only is no_default:
        if PANDAS_GE_200:
            numeric_only = False
        elif PANDAS_GE_150:
            maybe_warn = True
            numeric_only = True
        else:
            numeric_only = True

    all_numeric = df._get_numeric_data()._name == df._name
    if maybe_warn and not all_numeric:
        warnings.warn(
            "The default value of numeric_only will be `False` "
            "in a future version of Dask.",
            FutureWarning,
        )

    if numeric_only and not all_numeric:
        df = df._get_numeric_data()

    if scalar and len(df.columns) != 2:
        raise ValueError("scalar only valid for 2 column dataframe")

    token = tokenize(df, min_periods, scalar, split_every)

    funcname = "corr" if corr else "cov"
    a = f"{funcname}-chunk-{df._name}"
    dsk = {
        (a, i): (_cov_corr_chunk, f, corr) for (i, f) in enumerate(df.__dask_keys__())
    }

    prefix = f"{funcname}-combine-{df._name}-"
    k = df.npartitions
    b = a
    depth = 0
    while k > split_every:
        b = prefix + str(depth)
        for part_i, inds in enumerate(partition_all(split_every, range(k))):
            dsk[(b, part_i)] = (_cov_corr_combine, [(a, i) for i in inds], corr)
        k = part_i + 1
        a = b
        depth += 1

    name = f"{funcname}-{token}"
    dsk[(name, 0)] = (
        _cov_corr_agg,
        [(a, i) for i in range(k)],
        df.columns,
        min_periods,
        corr,
        scalar,
        df._meta,
    )
    graph = HighLevelGraph.from_collections(name, dsk, dependencies=[df])
    if scalar:
        return Scalar(graph, name, "f8")
    meta = make_meta(
        [(c, "f8") for c in df.columns],
        index=meta_series_constructor(df)(df.columns),
        parent_meta=df._meta,
    )
    return new_dd_object(graph, name, meta, (df.columns.min(), df.columns.max()))


def _cov_corr_chunk(df, corr=False):
    """Chunk part of a covariance or correlation computation"""
    shape = (df.shape[1], df.shape[1])
    kwargs = {} if PANDAS_GE_300 else {"copy": False}
    df = df.astype("float64", **kwargs)
    sums = np.zeros_like(df.values, shape=shape)
    counts = np.zeros_like(df.values, shape=shape)
    for idx in range(len(df.columns)):
        mask = df.iloc[:, idx].notnull()
        sums[idx] = df[mask].sum().values
        counts[idx] = df[mask].count().values
    # Special case single-row DataFrame cov to avoid warnings from pandas.
    if df.shape[0] == 1:
        cov = np.full_like(sums, np.nan)  # always an all nan result
    else:
        cov = df.cov().values
    dtype = [("sum", sums.dtype), ("count", counts.dtype), ("cov", cov.dtype)]
    if corr:
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            mu = (sums / counts).T
        m = np.zeros_like(df.values, shape=shape)
        mask = df.isnull().values
        for idx in range(len(df.columns)):
            # Avoid using ufunc.outer (not supported by cupy)
            mu_discrepancy = (
                np.subtract(df.iloc[:, idx].values[:, None], mu[idx][None, :]) ** 2
            )
            mu_discrepancy[mask] = np.nan
            m[idx] = np.nansum(mu_discrepancy, axis=0)
        m = m.T
        dtype.append(("m", m.dtype))

    out = {"sum": sums, "count": counts, "cov": cov * (counts - 1)}
    if corr:
        out["m"] = m
    return out


def _cov_corr_combine(data_in, corr=False):
    data = {"sum": None, "count": None, "cov": None}
    if corr:
        data["m"] = None

    for k in data.keys():
        data[k] = [d[k] for d in data_in]
        data[k] = np.concatenate(data[k]).reshape((len(data[k]),) + data[k][0].shape)

    sums = np.nan_to_num(data["sum"])
    counts = data["count"]

    cum_sums = np.cumsum(sums, 0)
    cum_counts = np.cumsum(counts, 0)

    s1 = cum_sums[:-1]
    s2 = sums[1:]
    n1 = cum_counts[:-1]
    n2 = counts[1:]
    with np.errstate(invalid="ignore"):
        d = (s2 / n2) - (s1 / n1)
        C = np.nansum(
            (n1 * n2) / (n1 + n2) * (d * d.transpose((0, 2, 1))), 0
        ) + np.nansum(data["cov"], 0)

    out = {"sum": cum_sums[-1], "count": cum_counts[-1], "cov": C}

    if corr:
        nobs = np.where(cum_counts[-1], cum_counts[-1], np.nan)
        mu = cum_sums[-1] / nobs
        counts_na = np.where(counts, counts, np.nan)
        m = np.nansum(data["m"] + counts * (sums / counts_na - mu) ** 2, axis=0)
        out["m"] = m
    return out


def _cov_corr_agg(data, cols, min_periods=2, corr=False, scalar=False, like_df=None):
    out = _cov_corr_combine(data, corr)
    counts = out["count"]
    C = out["cov"]
    C[counts < min_periods] = np.nan
    if corr:
        m2 = out["m"]
        den = np.sqrt(m2 * m2.T)
    else:
        den = np.where(counts, counts, np.nan) - 1
    with np.errstate(invalid="ignore", divide="ignore"):
        mat = C / den
    if scalar:
        return float(mat[0, 1])
    return (pd.DataFrame if like_df is None else meta_frame_constructor(like_df))(
        mat, columns=cols, index=cols
    )


def pd_split(df, p, random_state=None, shuffle=False):
    """Split DataFrame into multiple pieces pseudorandomly

    >>> df = pd.DataFrame({'a': [1, 2, 3, 4, 5, 6],
    ...                    'b': [2, 3, 4, 5, 6, 7]})

    >>> a, b = pd_split(
    ...     df, [0.5, 0.5], random_state=123, shuffle=True
    ... )  # roughly 50/50 split
    >>> a
       a  b
    3  4  5
    0  1  2
    5  6  7
    >>> b
       a  b
    1  2  3
    4  5  6
    2  3  4
    """
    p = list(p)
    if shuffle:
        if not isinstance(random_state, np.random.RandomState):
            random_state = np.random.RandomState(random_state)
        df = df.sample(frac=1.0, random_state=random_state)
    index = pseudorandom(len(df), p, random_state)
    return [df.iloc[index == i] for i in range(len(p))]


def _take_last(a, skipna=True):
    """
    take last row (Series) of DataFrame / last value of Series
    considering NaN.

    Parameters
    ----------
    a : pd.DataFrame or pd.Series
    skipna : bool, default True
        Whether to exclude NaN

    """

    def _last_valid(s):
        for i in range(1, min(10, len(s) + 1)):
            val = s.iloc[-i]
            if not pd.isnull(val):
                return val
        else:
            nonnull = s[s.notna()]
            if not nonnull.empty:
                return nonnull.iloc[-1]
        return None

    if skipna is False:
        return a.iloc[-1]
    else:
        # take last valid value excluding NaN, NaN location may be different
        # in each column
        if is_dataframe_like(a):
            # create Series from appropriate backend dataframe library
            series_typ = type(a.iloc[0:1, 0])
            if a.empty:
                return series_typ([], dtype="float")
            return series_typ(
                [_last_valid(a.iloc[:, i]) for i in range(len(a.columns))],
                index=a.columns,
            )
        else:
            return _last_valid(a)


def check_divisions(divisions):
    if not isinstance(divisions, (list, tuple)):
        raise ValueError("New division must be list or tuple")
    divisions = list(divisions)
    if len(divisions) == 0:
        raise ValueError("New division must not be empty")
    if divisions != sorted(divisions):
        raise ValueError("New division must be sorted")
    if len(divisions[:-1]) != len(list(unique(divisions[:-1]))):
        msg = "New division must be unique, except for the last element"
        raise ValueError(msg)


def repartition_divisions(a, b, name, out1, out2, force=False):
    """dask graph to repartition dataframe by new divisions

    Parameters
    ----------
    a : tuple
        old divisions
    b : tuple, listmypy
    out2 : str
        name of new dataframe
    force : bool, default False
        Allows the expansion of the existing divisions.
        If False then the new divisions lower and upper bounds must be
        the same as the old divisions.

    Examples
    --------
    >>> from pprint import pprint
    >>> pprint(repartition_divisions([1, 3, 7], [1, 4, 6, 7], 'a', 'b', 'c'))  # doctest: +ELLIPSIS
    {('b', 0): (<function boundary_slice at ...>, ('a', 0), 1, 3, False),
     ('b', 1): (<function boundary_slice at ...>, ('a', 1), 3, 4, False),
     ('b', 2): (<function boundary_slice at ...>, ('a', 1), 4, 6, False),
     ('b', 3): (<function boundary_slice at ...>, ('a', 1), 6, 7, True),
     ('c', 0): (<function concat at ...>, [('b', 0), ('b', 1)]),
     ('c', 1): ('b', 2),
     ('c', 2): ('b', 3)}
    """
    check_divisions(b)

    if len(b) < 2:
        # minimum division is 2 elements, like [0, 0]
        raise ValueError("New division must be longer than 2 elements")

    if force:
        if a[0] < b[0]:
            msg = (
                "left side of the new division must be equal or smaller "
                "than old division"
            )
            raise ValueError(msg)
        if a[-1] > b[-1]:
            msg = (
                "right side of the new division must be equal or larger "
                "than old division"
            )
            raise ValueError(msg)
    else:
        if a[0] != b[0]:
            msg = "left side of old and new divisions are different"
            raise ValueError(msg)
        if a[-1] != b[-1]:
            msg = "right side of old and new divisions are different"
            raise ValueError(msg)

    def _is_single_last_div(x):
        """Whether last division only contains single label"""
        return len(x) >= 2 and x[-1] == x[-2]

    c = [a[0]]
    d = dict()
    low = a[0]

    i, j = 1, 1  # indices for old/new divisions
    k = 0  # index for temp divisions

    last_elem = _is_single_last_div(a)

    # process through old division
    # left part of new division can be processed in this loop
    while i < len(a) and j < len(b):
        if a[i] < b[j]:
            # tuple is something like:
            # (methods.boundary_slice, ('from_pandas-#', 0), 3, 4, False))
            d[(out1, k)] = (methods.boundary_slice, (name, i - 1), low, a[i], False)
            low = a[i]
            i += 1
        elif a[i] > b[j]:
            d[(out1, k)] = (methods.boundary_slice, (name, i - 1), low, b[j], False)
            low = b[j]
            j += 1
        else:
            d[(out1, k)] = (methods.boundary_slice, (name, i - 1), low, b[j], False)
            low = b[j]
            if len(a) == i + 1 or a[i] < a[i + 1]:
                j += 1
            i += 1
        c.append(low)
        k += 1

    # right part of new division can remain
    if a[-1] < b[-1] or b[-1] == b[-2]:
        for _j in range(j, len(b)):
            # always use right-most of old division
            # because it may contain last element
            m = len(a) - 2
            d[(out1, k)] = (methods.boundary_slice, (name, m), low, b[_j], False)
            low = b[_j]
            c.append(low)
            k += 1
    else:
        # even if new division is processed through,
        # right-most element of old division can remain
        if last_elem and i < len(a):
            d[(out1, k)] = (methods.boundary_slice, (name, i - 1), a[i], a[i], False)
            k += 1
        c.append(a[-1])

    # replace last element of tuple with True
    d[(out1, k - 1)] = d[(out1, k - 1)][:-1] + (True,)

    i, j = 0, 1

    last_elem = _is_single_last_div(c)

    while j < len(b):
        tmp = []
        while c[i] < b[j]:
            tmp.append((out1, i))
            i += 1
        while (
            last_elem
            and c[i] == b[-1]
            and (b[-1] != b[-2] or j == len(b) - 1)
            and i < k
        ):
            # append if last split is not included
            tmp.append((out1, i))
            i += 1
        if len(tmp) == 0:
            # dummy slice to return empty DataFrame or Series,
            # which retain original data attributes (columns / name)
            d[(out2, j - 1)] = (methods.boundary_slice, (name, 0), a[0], a[0], False)
        elif len(tmp) == 1:
            d[(out2, j - 1)] = tmp[0]
        else:
            if not tmp:
                raise ValueError(
                    "check for duplicate partitions\nold:\n%s\n\n"
                    "new:\n%s\n\ncombined:\n%s" % (pformat(a), pformat(b), pformat(c))
                )
            d[(out2, j - 1)] = (methods.concat, tmp)
        j += 1
    return d


def repartition_freq(df, freq=None):
    """Repartition a timeseries dataframe by a new frequency"""
    if not isinstance(df.divisions[0], pd.Timestamp):
        raise TypeError("Can only repartition on frequency for timeseries")

    freq = _map_freq_to_period_start(freq)

    try:
        start = df.divisions[0].ceil(freq)
    except ValueError:
        start = df.divisions[0]
    divisions = methods.tolist(
        pd.date_range(start=start, end=df.divisions[-1], freq=freq)
    )
    if not len(divisions):
        divisions = [df.divisions[0], df.divisions[-1]]
    else:
        divisions.append(df.divisions[-1])
        if divisions[0] != df.divisions[0]:
            divisions = [df.divisions[0]] + divisions

    return df.repartition(divisions=divisions)


def _map_freq_to_period_start(freq):
    """Ensure that the frequency pertains to the **start** of a period.

    If e.g. `freq='M'`, then the divisions are:
        - 2021-31-1 00:00:00 (start of February partition)
        - 2021-2-28 00:00:00 (start of March partition)
        - ...

    but this **should** be:
        - 2021-2-1 00:00:00 (start of February partition)
        - 2021-3-1 00:00:00 (start of March partition)
        - ...

    Therefore, we map `freq='M'` to `freq='MS'` (same for quarter and year).
    """

    if not isinstance(freq, str):
        return freq

    offset = pd.tseries.frequencies.to_offset(freq)
    offset_type_name = type(offset).__name__

    if not offset_type_name.endswith("End"):
        return freq

    new_offset = offset_type_name[: -len("End")] + "Begin"
    try:
        new_offset_type = getattr(pd.tseries.offsets, new_offset)
        if "-" in freq:
            _, anchor = freq.split("-")
            anchor = "-" + anchor
        else:
            anchor = ""
        n = str(offset.n) if offset.n != 1 else ""
        return f"{n}{new_offset_type._prefix}{anchor}"
    except AttributeError:
        return freq


def repartition_size(df, size):
    """
    Repartition dataframe so that new partitions have approximately `size` memory usage each
    """
    if isinstance(size, str):
        size = parse_bytes(size)
    size = int(size)

    mem_usages = df.map_partitions(total_mem_usage, deep=True).compute()

    # 1. split each partition that is larger than partition_size
    nsplits = 1 + mem_usages // size
    if np.any(nsplits > 1):
        split_name = f"repartition-split-{size}-{tokenize(df)}"
        df = _split_partitions(df, nsplits, split_name)
        # update mem_usages to account for the split partitions
        split_mem_usages = []
        for n, usage in zip(nsplits, mem_usages):
            split_mem_usages.extend([usage / n] * n)
        mem_usages = pd.Series(split_mem_usages)

    # 2. now that all partitions are less than size, concat them up to size
    assert np.all(mem_usages <= size)
    new_npartitions = list(map(len, iter_chunks(mem_usages, size)))
    new_partitions_boundaries = np.cumsum(new_npartitions)
    new_name = f"repartition-{size}-{tokenize(df)}"
    return _repartition_from_boundaries(df, new_partitions_boundaries, new_name)


def total_mem_usage(df, index=True, deep=False):
    mem_usage = df.memory_usage(index=index, deep=deep)
    if is_series_like(mem_usage):
        mem_usage = mem_usage.sum()
    return mem_usage


def repartition_npartitions(df, npartitions):
    """Repartition dataframe to a smaller number of partitions"""
    new_name = "repartition-%d-%s" % (npartitions, tokenize(df))
    if df.npartitions == npartitions:
        return df
    elif df.npartitions > npartitions:
        npartitions_ratio = df.npartitions / npartitions
        new_partitions_boundaries = [
            int(new_partition_index * npartitions_ratio)
            for new_partition_index in range(npartitions + 1)
        ]
        return _repartition_from_boundaries(df, new_partitions_boundaries, new_name)
    else:
        # Drop duplicates in case last partition has same
        # value for min and max division
        original_divisions = divisions = pd.Series(df.divisions).drop_duplicates()
        if df.known_divisions and (
            is_datetime64_any_dtype(divisions.dtype)
            or is_numeric_dtype(divisions.dtype)
        ):
            if is_datetime64_any_dtype(divisions.dtype):
                divisions = divisions.values.astype("float64")

            if is_series_like(divisions):
                divisions = divisions.values

            n = len(divisions)
            divisions = np.interp(
                x=np.linspace(0, n, npartitions + 1),
                xp=np.linspace(0, n, n),
                fp=divisions,
            )
            if is_datetime64_any_dtype(original_divisions.dtype):
                divisions = methods.tolist(
                    pd.Series(divisions).astype(original_divisions.dtype)
                )
            elif np.issubdtype(original_divisions.dtype, np.integer):
                divisions = divisions.astype(original_divisions.dtype)

            if isinstance(divisions, np.ndarray):
                divisions = divisions.tolist()

            divisions = list(divisions)
            divisions[0] = df.divisions[0]
            divisions[-1] = df.divisions[-1]

            # Ensure the computed divisions are unique
            divisions = list(unique(divisions[:-1])) + [divisions[-1]]
            return df.repartition(divisions=divisions)
        else:
            div, mod = divmod(npartitions, df.npartitions)
            nsplits = [div] * df.npartitions
            nsplits[-1] += mod
            return _split_partitions(df, nsplits, new_name)


def _repartition_from_boundaries(df, new_partitions_boundaries, new_name):
    if not isinstance(new_partitions_boundaries, list):
        new_partitions_boundaries = list(new_partitions_boundaries)
    if new_partitions_boundaries[0] > 0:
        new_partitions_boundaries.insert(0, 0)
    if new_partitions_boundaries[-1] < df.npartitions:
        new_partitions_boundaries.append(df.npartitions)
    dsk = {}
    for i, (start, end) in enumerate(
        zip(new_partitions_boundaries, new_partitions_boundaries[1:])
    ):
        dsk[new_name, i] = (methods.concat, [(df._name, j) for j in range(start, end)])
    divisions = [df.divisions[i] for i in new_partitions_boundaries]
    graph = HighLevelGraph.from_collections(new_name, dsk, dependencies=[df])
    return new_dd_object(graph, new_name, df._meta, divisions)


def _split_partitions(df, nsplits, new_name):
    """Split a Dask dataframe into new partitions

    Parameters
    ----------
    df: DataFrame or Series
    nsplits: List[int]
        Number of target dataframes for each partition
        The length of nsplits should be the same as df.npartitions
    new_name: str

    See Also
    --------
    repartition_npartitions
    repartition_size
    """
    if len(nsplits) != df.npartitions:
        raise ValueError(f"nsplits should have len={df.npartitions}")

    dsk = {}
    split_name = f"split-{tokenize(df, nsplits)}"
    j = 0
    for i, k in enumerate(nsplits):
        if k == 1:
            dsk[new_name, j] = (df._name, i)
            j += 1
        else:
            dsk[split_name, i] = (split_evenly, (df._name, i), k)
            for jj in range(k):
                dsk[new_name, j] = (getitem, (split_name, i), jj)
                j += 1

    divisions = [None] * (1 + sum(nsplits))
    graph = HighLevelGraph.from_collections(new_name, dsk, dependencies=[df])
    return new_dd_object(graph, new_name, df._meta, divisions)


def repartition(df, divisions=None, force=False):
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

    # no-op fastpath for when we already have matching divisions
    if is_dask_collection(df) and df.divisions == tuple(divisions):
        return df

    token = tokenize(df, divisions)
    if isinstance(df, _Frame):
        tmp = "repartition-split-" + token
        out = "repartition-merge-" + token
        dsk = repartition_divisions(
            df.divisions, divisions, df._name, tmp, out, force=force
        )
        graph = HighLevelGraph.from_collections(out, dsk, dependencies=[df])
        return new_dd_object(graph, out, df._meta, divisions)
    elif is_dataframe_like(df) or is_series_like(df):
        name = "repartition-dataframe-" + token
        from dask.dataframe.utils import shard_df_on_index

        dfs = shard_df_on_index(df, divisions[1:-1])
        dsk = {(name, i): df for i, df in enumerate(dfs)}
        return new_dd_object(dsk, name, df, divisions)
    raise ValueError("Data must be DataFrame or Series")


def _reduction_chunk(x, aca_chunk=None, **kwargs):
    o = aca_chunk(x, **kwargs)
    # Return a dataframe so that the concatenated version is also a dataframe
    if not is_series_like(o):
        return o
    result = o.to_frame().T
    if o.dtype.kind == "O":
        # Was coerced to object, so cast back
        return result.infer_objects()
    return result


def _reduction_combine(x, aca_combine=None, **kwargs):
    if isinstance(x, list):
        x = pd.Series(x)
    o = aca_combine(x, **kwargs)
    # Return a dataframe so that the concatenated version is also a dataframe
    return o.to_frame().T if is_series_like(o) else o


def _reduction_aggregate(x, aca_aggregate=None, **kwargs):
    if isinstance(x, list):
        x = pd.Series(x)
    return aca_aggregate(x, **kwargs)


def idxmaxmin_chunk(x, fn=None, skipna=True, numeric_only=False):
    numeric_only_kwargs = (
        {} if not PANDAS_GE_150 or is_series_like(x) else {"numeric_only": numeric_only}
    )
    minmax = "max" if fn == "idxmax" else "min"
    if len(x) > 0:
        idx = getattr(x, fn)(skipna=skipna, **numeric_only_kwargs)
        value = getattr(x, minmax)(skipna=skipna, **numeric_only_kwargs)
    else:
        idx = value = meta_series_constructor(x)([], dtype="i8")
    if is_series_like(idx):
        return meta_frame_constructor(x)({"idx": idx, "value": value})
    return meta_frame_constructor(x)({"idx": [idx], "value": [value]})


def idxmaxmin_row(x, fn=None, skipna=True):
    minmax = "max" if fn == "idxmax" else "min"
    if len(x) > 0:
        x = x.set_index("idx")
        # potentially coerced to object, so cast back
        value = x.value.infer_objects()
        idx = [getattr(value, fn)(skipna=skipna)]
        value = [getattr(value, minmax)(skipna=skipna)]
    else:
        idx = value = meta_series_constructor(x)([], dtype="i8")
    return meta_frame_constructor(x)(
        {
            "idx": meta_series_constructor(x)(idx, dtype=x.index.dtype),
            "value": meta_series_constructor(x)(value, dtype=x.dtypes.iloc[0]),
        }
    )


def idxmaxmin_combine(x, fn=None, skipna=True):
    if len(x) <= 1:
        return x
    return (
        x.groupby(level=0)
        .apply(idxmaxmin_row, fn=fn, skipna=skipna)
        .reset_index(level=1, drop=True)
    )


def idxmaxmin_agg(x, fn=None, skipna=True, scalar=False, numeric_only=no_default):
    res = idxmaxmin_combine(x, fn, skipna=skipna)["idx"]
    if len(res) == 0:
        raise ValueError("attempt to get argmax of an empty sequence")
    if scalar:
        return res[0]
    res.name = None
    return res


def _mode_aggregate(df, dropna):
    value_count_series = df.sum()
    max_val = value_count_series.max(skipna=dropna)
    mode_series = (
        value_count_series[value_count_series == max_val]
        .index.to_series()
        .sort_values()
        .reset_index(drop=True)
    )
    return mode_series


def _count_aggregate(x):
    return x.sum().astype("int64")


def safe_head(df, n):
    r = M.head(df, n)
    if len(r) != n:
        warnings.warn(
            f"Insufficient elements for `head`. {n} elements requested, only {len(r)} "
            "elements available. Try passing larger `npartitions` to `head`."
        )
    return r


def maybe_shift_divisions(df, periods, freq):
    """Maybe shift divisions by periods of size freq

    Used to shift the divisions for the `shift` method. If freq isn't a fixed
    size (not anchored or relative), then the divisions are shifted
    appropriately. Otherwise the divisions are cleared.

    Parameters
    ----------
    df : dd.DataFrame, dd.Series, or dd.Index
    periods : int
        The number of periods to shift.
    freq : DateOffset, timedelta, or time rule string
        The frequency to shift by.
    """
    if isinstance(freq, str):
        freq = pd.tseries.frequencies.to_offset(freq)

    is_offset = isinstance(freq, pd.DateOffset)
    if is_offset:
        if not isinstance(freq, pd.offsets.Tick):
            # Can't infer divisions on relative or anchored offsets, as
            # divisions may now split identical index value.
            # (e.g. index_partitions = [[1, 2, 3], [3, 4, 5]])
            return df.clear_divisions()
    if df.known_divisions:
        divs = pd.Series(range(len(df.divisions)), index=df.divisions)
        divisions = divs.shift(periods, freq=freq).index
        return df.__class__(df.dask, df._name, df._meta, divisions)
    return df


@wraps(pd.to_datetime)
def to_datetime(arg, meta=None, **kwargs):
    tz_kwarg = {"tz": "utc"} if kwargs.get("utc") else {}
    if meta is None:
        if isinstance(arg, Index):
            meta = get_meta_library(arg).DatetimeIndex([], **tz_kwarg)
            meta.name = arg.name
        elif not (is_dataframe_like(arg) or is_series_like(arg)):
            raise NotImplementedError(
                "dask.dataframe.to_datetime does not support "
                "non-index-able arguments (like scalars)"
            )
        else:
            meta = meta_series_constructor(arg)([pd.Timestamp("2000", **tz_kwarg)])
            meta.index = meta.index.astype(arg.index.dtype)
            meta.index.name = arg.index.name
    if PANDAS_GE_200 and "infer_datetime_format" in kwargs:
        warnings.warn(
            "The argument 'infer_datetime_format' is deprecated and will be removed in a future version. "
            "A strict version of it is now the default, see "
            "https://pandas.pydata.org/pdeps/0004-consistent-to-datetime-parsing.html. "
            "You can safely remove this argument.",
            UserWarning,
        )
        kwargs.pop("infer_datetime_format")

    return map_partitions(get_meta_library(arg).to_datetime, arg, meta=meta, **kwargs)


@wraps(pd.to_timedelta)
def to_timedelta(arg, unit=None, errors="raise"):
    meta = meta_series_constructor(arg)([pd.Timedelta(1, unit=unit)])
    return map_partitions(pd.to_timedelta, arg, unit=unit, errors=errors, meta=meta)


if hasattr(pd, "isna"):

    @wraps(pd.isna)
    def isna(arg):
        return map_partitions(pd.isna, arg)


def _repr_data_series(s, index):
    """A helper for creating the ``_repr_data`` property"""
    npartitions = len(index) - 1
    if isinstance(s.dtype, pd.CategoricalDtype):
        if has_known_categories(s):
            dtype = "category[known]"
        else:
            dtype = "category[unknown]"
    else:
        dtype = str(s.dtype)
    return pd.Series([dtype] + ["..."] * npartitions, index=index, name=s.name)


def has_parallel_type(x):
    """Does this object have a dask dataframe equivalent?"""
    return get_parallel_type(x) is not Scalar


def new_dd_object(dsk, name, meta, divisions, parent_meta=None):
    """Generic constructor for dask.dataframe objects.

    Decides the appropriate output class based on the type of `meta` provided.
    """
    if has_parallel_type(meta):
        return get_parallel_type(meta)(dsk, name, meta, divisions)
    elif is_arraylike(meta) and meta.shape:
        import dask.array as da

        chunks = ((np.nan,) * (len(divisions) - 1),) + tuple(
            (d,) for d in meta.shape[1:]
        )
        if len(chunks) > 1:
            if isinstance(dsk, HighLevelGraph):
                layer = dsk.layers[name]
            else:
                # dask-expr provides a dict only
                layer = dsk
            if isinstance(layer, Blockwise):
                layer.new_axes["j"] = chunks[1][0]
                layer.output_indices = layer.output_indices + ("j",)
            else:
                suffix = (0,) * (len(chunks) - 1)
                for i in range(len(chunks[0])):
                    layer[(name, i) + suffix] = layer.pop((name, i))
        return da.Array(dsk, name=name, chunks=chunks, dtype=meta.dtype)
    else:
        return get_parallel_type(meta)(dsk, name, meta, divisions)


def partitionwise_graph(func, layer_name, *args, **kwargs):
    """
    Apply a function partition-wise across arguments to create layer of a graph

    This applies a function, ``func``, in an embarrassingly parallel fashion
    across partitions/chunks in the provided arguments.  It handles Dataframes,
    Arrays, and scalars smoothly, and relies on the ``blockwise`` machinery
    to provide a nicely symbolic graph.

    It is most commonly used in other graph-building functions to create the
    appropriate layer of the resulting dataframe.

    Parameters
    ----------
    func: callable
    layer_name: str
        Descriptive name for the operation. Used as the output name
        in the resulting ``Blockwise`` graph layer.
    *args:
    **kwargs:

    Returns
    -------
    out: Blockwise graph

    Examples
    --------
    >>> subgraph = partitionwise_graph(function, x, y, z=123)  # doctest: +SKIP
    >>> layer = partitionwise_graph(function, df, x, z=123)  # doctest: +SKIP
    >>> graph = HighLevelGraph.from_collections(name, layer, dependencies=[df, x])  # doctest: +SKIP
    >>> result = new_dd_object(graph, name, metadata, df.divisions)  # doctest: +SKIP

    See Also
    --------
    map_partitions
    """
    pairs = []
    numblocks = {}
    for arg in args:
        if isinstance(arg, _Frame):
            pairs.extend([arg._name, "i"])
            numblocks[arg._name] = (arg.npartitions,)
        elif isinstance(arg, Scalar):
            pairs.extend([arg._name, "i"])
            numblocks[arg._name] = (1,)
        elif isinstance(arg, Array):
            if arg.ndim == 1:
                pairs.extend([arg.name, "i"])
            elif arg.ndim == 0:
                pairs.extend([arg.name, ""])
            elif arg.ndim == 2:
                pairs.extend([arg.name, "ij"])
            else:
                raise ValueError("Can't add multi-dimensional array to dataframes")
            numblocks[arg._name] = arg.numblocks
        elif isinstance(arg, BlockwiseDep):
            if len(arg.numblocks) == 1:
                pairs.extend([arg, "i"])
            elif len(arg.numblocks) == 2:
                pairs.extend([arg, "ij"])
            else:
                raise ValueError(
                    f"BlockwiseDep arg {arg!r} has {len(arg.numblocks)} dimensions; only 1 or 2 are supported."
                )
        else:
            pairs.extend([arg, None])
    return blockwise(
        func, layer_name, "i", *pairs, numblocks=numblocks, concatenate=True, **kwargs
    )


def meta_warning(df):
    """
    Provide an informative message when the user is asked to provide metadata
    """
    if is_dataframe_like(df):
        meta_str = {k: str(v) for k, v in df.dtypes.to_dict().items()}
    elif is_series_like(df):
        meta_str = (df.name, str(df.dtype))
    else:
        meta_str = None
    msg = (
        "\nYou did not provide metadata, so Dask is running your "
        "function on a small dataset to guess output types. "
        "It is possible that Dask will guess incorrectly.\n"
        "To provide an explicit output types or to silence this message, "
        "please provide the `meta=` keyword, as described in the map or "
        "apply function that you are using."
    )
    if meta_str:
        msg += (
            "\n"
            "  Before: .apply(func)\n"
            "  After:  .apply(func, meta=%s)\n" % str(meta_str)
        )
    return msg


def prefix_reduction(f, ddf, identity, **kwargs):
    """Computes the prefix sums of f on df

    If df has partitions [P1, P2, ..., Pn], then returns the DataFrame with
    partitions [f(identity, P1),
                f(f(identity, P1), P2),
                f(f(f(identity, P1), P2), P3),
                ...]

    Parameters
    ----------
    f : callable
        an associative function f
    ddf : dd.DataFrame
    identity : pd.DataFrame
        an identity element of f, that is f(identity, df) = f(df, identity) = df
    """
    dsk = dict()
    name = "prefix_reduction-" + tokenize(f, ddf, identity, **kwargs)
    meta = ddf._meta
    n = len(ddf.divisions) - 1
    divisions = [None] * (n + 1)

    N = 1
    while N < n:
        N *= 2
    for i in range(n):
        dsk[(name, i, 1, 0)] = (apply, f, [(ddf._name, i), identity], kwargs)
    for i in range(n, N):
        dsk[(name, i, 1, 0)] = identity

    d = 1
    while d < N:
        for i in range(0, N, 2 * d):
            dsk[(name, i + 2 * d - 1, 2 * d, 0)] = (
                apply,
                f,
                [(name, i + d - 1, d, 0), (name, i + 2 * d - 1, d, 0)],
                kwargs,
            )
        d *= 2

    dsk[(name, N - 1, N, 1)] = identity

    while d > 1:
        d //= 2
        for i in range(0, N, 2 * d):
            dsk[(name, i + d - 1, d, 1)] = (name, i + 2 * d - 1, 2 * d, 1)
            dsk[(name, i + 2 * d - 1, d, 1)] = (
                apply,
                f,
                [(name, i + 2 * d - 1, 2 * d, 1), (name, i + d - 1, d, 0)],
                kwargs,
            )

    for i in range(n):
        dsk[(name, i)] = (apply, f, [(name, i, 1, 1), identity], kwargs)

    graph = HighLevelGraph.from_collections(name, dsk, dependencies=[ddf])
    return new_dd_object(graph, name, meta, divisions)


def suffix_reduction(f, ddf, identity, **kwargs):
    """Computes the suffix sums of f on df

    If df has partitions [P1, P2, ..., Pn], then returns the DataFrame with
    partitions [f(P1, f(P2, ...f(Pn, identity)...)),
                f(P2, ...f(Pn, identity)...),
                ...f(Pn, identity)...,
                ...]

    Parameters
    ----------
    f : callable
        an associative function f
    ddf : dd.DataFrame
    identity : pd.DataFrame
        an identity element of f, that is f(identity, df) = f(df, identity) = df
    kwargs : ??
        keyword arguments of f ??
    """
    dsk = dict()
    name = "suffix_reduction-" + tokenize(f, ddf, identity, **kwargs)
    meta = ddf._meta
    n = len(ddf.divisions) - 1
    divisions = [None] * (n + 1)

    N = 1
    while N < n:
        N *= 2
    for i in range(n):
        dsk[(name, i, 1, 0)] = (apply, f, [(ddf._name, n - 1 - i), identity], kwargs)
    for i in range(n, N):
        dsk[(name, i, 1, 0)] = identity

    d = 1
    while d < N:
        for i in range(0, N, 2 * d):
            dsk[(name, i + 2 * d - 1, 2 * d, 0)] = (
                apply,
                f,
                [(name, i + 2 * d - 1, d, 0), (name, i + d - 1, d, 0)],
                kwargs,
            )
        d *= 2

    dsk[(name, N - 1, N, 1)] = identity

    while d > 1:
        d //= 2
        for i in range(0, N, 2 * d):
            dsk[(name, i + d - 1, d, 1)] = (name, i + 2 * d - 1, 2 * d, 1)
            dsk[(name, i + 2 * d - 1, d, 1)] = (
                apply,
                f,
                [(name, i + d - 1, d, 0), (name, i + 2 * d - 1, 2 * d, 1)],
                kwargs,
            )

    for i in range(n):
        dsk[(name, i)] = (apply, f, [(name, n - 1 - i, 1, 1), identity], kwargs)

    graph = HighLevelGraph.from_collections(name, dsk, dependencies=[ddf])
    return new_dd_object(graph, name, meta, divisions)


def mapseries(base_chunk, concat_map):
    return base_chunk.map(concat_map)


def mapseries_combine(index, concat_result):
    final_series = concat_result.sort_index()
    final_series = index.to_series().map(final_series)
    return final_series


def series_map(base_series, map_series):
    npartitions = base_series.npartitions
    split_out = map_series.npartitions

    dsk = {}

    base_token_key = tokenize(base_series, split_out)
    base_split_prefix = f"base-split-{base_token_key}"
    base_shard_prefix = f"base-shard-{base_token_key}"
    for i, key in enumerate(base_series.__dask_keys__()):
        dsk[(base_split_prefix, i)] = (hash_shard, key, split_out)
        for j in range(split_out):
            dsk[(base_shard_prefix, 0, i, j)] = (getitem, (base_split_prefix, i), j)

    map_token_key = tokenize(map_series)
    map_split_prefix = f"map-split-{map_token_key}"
    map_shard_prefix = f"map-shard-{map_token_key}"
    for i, key in enumerate(map_series.__dask_keys__()):
        dsk[(map_split_prefix, i)] = (
            hash_shard,
            key,
            split_out,
            split_out_on_index,
            None,
        )
        for j in range(split_out):
            dsk[(map_shard_prefix, 0, i, j)] = (getitem, (map_split_prefix, i), j)

    token_key = tokenize(base_series, map_series)
    map_prefix = f"map-series-{token_key}"
    for i in range(npartitions):
        for j in range(split_out):
            dsk[(map_prefix, i, j)] = (
                mapseries,
                (base_shard_prefix, 0, i, j),
                (_concat, [(map_shard_prefix, 0, k, j) for k in range(split_out)]),
            )

    final_prefix = f"map-series-combine-{token_key}"
    for i, key in enumerate(base_series.index.__dask_keys__()):
        dsk[(final_prefix, i)] = (
            mapseries_combine,
            key,
            (_concat, [(map_prefix, i, j) for j in range(split_out)]),
        )

    meta = map_series._meta.copy()
    meta.index = base_series._meta.index
    meta = make_meta(meta)

    dependencies = [base_series, map_series, base_series.index]
    graph = HighLevelGraph.from_collections(
        final_prefix, dsk, dependencies=dependencies
    )
    divisions = list(base_series.divisions)

    return new_dd_object(graph, final_prefix, meta, divisions)


def _convert_to_numeric(series, skipna):
    if PANDAS_GE_200:
        if skipna:
            return series.dropna().astype("i8")

        # series.view("i8") with pd.NaT produces -9223372036854775808 is why we need to do this
        return series.astype("i8").mask(series.isnull(), np.nan)
    else:
        view = "_view" if isinstance(series, Series) else "view"
        if skipna:
            return getattr(series.dropna(), view)("i8")

        # series.view("i8") with pd.NaT produces -9223372036854775808 is why we need to do this
        return getattr(series, view)("i8").mask(series.isnull(), np.nan)


def _sqrt_and_convert_to_timedelta(partition, axis, dtype=None, *args, **kwargs):
    if axis == 1:
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                category=RuntimeWarning,
                message="invalid value encountered in cast",
            )
            return pd.to_timedelta(M.std(partition, *args, axis=axis, **kwargs))

    is_df_like, time_cols = kwargs["is_df_like"], kwargs["time_cols"]

    sqrt = np.sqrt(partition)

    if not is_df_like:
        return pd.to_timedelta(sqrt)

    time_col_mask = sqrt.index.isin(time_cols)
    matching_vals = sqrt[time_col_mask]
    if len(time_cols) > 0:
        sqrt = sqrt.astype(object)
    for time_col, matching_val in zip(time_cols, matching_vals):
        sqrt[time_col] = pd.to_timedelta(matching_val)

    if dtype is not None:
        sqrt = sqrt.astype(dtype)
    return sqrt


def _raise_if_not_series_or_dataframe(x, funcname):
    """
    Utility function to raise an error if an object is not a Series or DataFrame
    """
    if not is_series_like(x) and not is_dataframe_like(x):
        raise NotImplementedError(
            "`%s` is only supported with objects that are Dataframes or Series"
            % funcname
        )


def _getattr_numeric_only(*args, _dask_method_name, **kwargs):
    with check_numeric_only_deprecation():
        return getattr(M, _dask_method_name)(*args, **kwargs)
