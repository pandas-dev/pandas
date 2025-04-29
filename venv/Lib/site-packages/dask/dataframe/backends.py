from __future__ import annotations

import warnings
from collections.abc import Iterable

import numpy as np
import pandas as pd
from pandas.api.types import is_scalar, union_categoricals

from dask.array import Array
from dask.array.dispatch import percentile_lookup
from dask.array.percentile import _percentile
from dask.backends import CreationDispatch, DaskBackendEntrypoint
from dask.dataframe._compat import PANDAS_GE_220, is_any_real_numeric_dtype
from dask.dataframe.dispatch import (
    categorical_dtype_dispatch,
    concat,
    concat_dispatch,
    from_pyarrow_table_dispatch,
    get_parallel_type,
    group_split_dispatch,
    grouper_dispatch,
    hash_object_dispatch,
    is_categorical_dtype_dispatch,
    make_meta_dispatch,
    make_meta_obj,
    meta_lib_from_array,
    meta_nonempty,
    partd_encode_dispatch,
    pyarrow_schema_dispatch,
    to_pandas_dispatch,
    to_pyarrow_table_dispatch,
    tolist_dispatch,
    union_categoricals_dispatch,
)
from dask.dataframe.extensions import make_array_nonempty, make_scalar
from dask.dataframe.utils import (
    _empty_series,
    _nonempty_scalar,
    _scalar_from_dtype,
    is_float_na_dtype,
    is_integer_na_dtype,
)
from dask.sizeof import SimpleSizeof, sizeof
from dask.utils import is_arraylike, is_series_like, typename


class DataFrameBackendEntrypoint(DaskBackendEntrypoint):
    """Dask-DataFrame version of ``DaskBackendEntrypoint``

    See Also
    --------
    PandasBackendEntrypoint
    """

    @staticmethod
    def from_dict(data: dict, *, npartitions: int, **kwargs):
        """Create a DataFrame collection from a dictionary

        Parameters
        ----------
        data : dict
            Of the form {field : array-like} or {field : dict}.
        npartitions : int
            The desired number of output partitions.
        **kwargs :
            Optional backend kwargs.

        See Also
        --------
        dask.dataframe.io.io.from_dict
        """
        raise NotImplementedError

    @staticmethod
    def read_parquet(path: str | list, **kwargs):
        """Read Parquet files into a DataFrame collection

        Parameters
        ----------
        path : str or list
            Source path(s).
        **kwargs :
            Optional backend kwargs.

        See Also
        --------
        dask.dataframe.io.parquet.core.read_parquet
        """
        raise NotImplementedError

    @staticmethod
    def read_json(url_path: str | list, **kwargs):
        """Read json files into a DataFrame collection

        Parameters
        ----------
        url_path : str or list
            Source path(s).
        **kwargs :
            Optional backend kwargs.

        See Also
        --------
        dask.dataframe.io.json.read_json
        """
        raise NotImplementedError

    @staticmethod
    def read_orc(path: str | list, **kwargs):
        """Read ORC files into a DataFrame collection

        Parameters
        ----------
        path : str or list
            Source path(s).
        **kwargs :
            Optional backend kwargs.

        See Also
        --------
        dask.dataframe.io.orc.core.read_orc
        """
        raise NotImplementedError

    @staticmethod
    def read_csv(urlpath: str | list, **kwargs):
        """Read CSV files into a DataFrame collection

        Parameters
        ----------
        urlpath : str or list
            Source path(s).
        **kwargs :
            Optional backend kwargs.

        See Also
        --------
        dask.dataframe.io.csv.read_csv
        """
        raise NotImplementedError

    @staticmethod
    def read_hdf(pattern: str | list, key: str, **kwargs):
        """Read HDF5 files into a DataFrame collection

        Parameters
        ----------
        pattern : str or list
            Source path(s).
        key : str
            Group identifier in the store.
        **kwargs :
            Optional backend kwargs.

        See Also
        --------
        dask.dataframe.io.hdf.read_hdf
        """
        raise NotImplementedError


dataframe_creation_dispatch = CreationDispatch(
    module_name="dataframe",
    default="pandas",
    entrypoint_class=DataFrameBackendEntrypoint,
    name="dataframe_creation_dispatch",
)


##########
# Pandas #
##########


@make_scalar.register(np.dtype)
def _(dtype):
    return _scalar_from_dtype(dtype)


@make_scalar.register(pd.Timestamp)
@make_scalar.register(pd.Timedelta)
@make_scalar.register(pd.Period)
@make_scalar.register(pd.Interval)
def _(x):
    return x


@make_meta_dispatch.register((pd.Series, pd.DataFrame))
def _(x, index=None):
    out = x.iloc[:0].copy(deep=True)
    # index isn't copied by default in pandas, even if deep=true
    out.index = out.index.copy(deep=True)
    return out


@make_meta_dispatch.register(pd.Index)
def _(x, index=None):
    return x[0:0].copy(deep=True)


meta_object_types: tuple[type, ...] = (pd.Series, pd.DataFrame, pd.Index, pd.MultiIndex)
try:
    import scipy.sparse as sp

    meta_object_types += (sp.spmatrix,)
except ImportError:
    pass


@pyarrow_schema_dispatch.register((pd.DataFrame,))
def get_pyarrow_schema_pandas(obj, preserve_index=None):
    import pyarrow as pa

    return pa.Schema.from_pandas(obj, preserve_index=preserve_index)


@to_pyarrow_table_dispatch.register((pd.DataFrame,))
def get_pyarrow_table_from_pandas(obj, **kwargs):
    # `kwargs` must be supported by `pyarrow.Table.to_pandas`
    import pyarrow as pa

    return pa.Table.from_pandas(obj, **kwargs)


@from_pyarrow_table_dispatch.register((pd.DataFrame,))
def get_pandas_dataframe_from_pyarrow(meta, table, **kwargs):
    # `kwargs` must be supported by `pyarrow.Table.to_pandas`
    import pyarrow as pa

    def default_types_mapper(pyarrow_dtype: pa.DataType) -> object:
        # Avoid converting strings from `string[pyarrow]` to
        # `string[python]` if we have *any* `string[pyarrow]`
        if (
            pyarrow_dtype in {pa.large_string(), pa.string()}
            and pd.StringDtype("pyarrow") in meta.dtypes.values
        ):
            return pd.StringDtype("pyarrow")
        return None

    types_mapper = kwargs.pop("types_mapper", default_types_mapper)
    return table.to_pandas(types_mapper=types_mapper, **kwargs)


@partd_encode_dispatch.register(pd.DataFrame)
def partd_pandas_blocks(_):
    from partd import PandasBlocks

    return PandasBlocks


@meta_nonempty.register(pd.DatetimeTZDtype)
@make_meta_dispatch.register(pd.DatetimeTZDtype)
def make_meta_pandas_datetime_tz(x, index=None):
    return _nonempty_scalar(x)


@make_meta_obj.register(meta_object_types)
def make_meta_object(x, index=None):
    """Create an empty pandas object containing the desired metadata.

    Parameters
    ----------
    x : dict, tuple, list, pd.Series, pd.DataFrame, pd.Index, dtype, scalar
        To create a DataFrame, provide a `dict` mapping of `{name: dtype}`, or
        an iterable of `(name, dtype)` tuples. To create a `Series`, provide a
        tuple of `(name, dtype)`. If a pandas object, names, dtypes, and index
        should match the desired output. If a dtype or scalar, a scalar of the
        same dtype is returned.
    index :  pd.Index, optional
        Any pandas index to use in the metadata. If none provided, a
        `RangeIndex` will be used.

    Examples
    --------

    >>> make_meta_object([('a', 'i8'), ('b', 'O')])
    Empty DataFrame
    Columns: [a, b]
    Index: []
    >>> make_meta_object(('a', 'f8'))
    Series([], Name: a, dtype: float64)
    >>> make_meta_object('i8')
    np.int64(1)
    """

    if is_arraylike(x) and x.shape:
        return x[:0]

    if index is not None:
        index = make_meta_dispatch(index)

    if isinstance(x, dict):
        return pd.DataFrame(
            {c: _empty_series(c, d, index=index) for (c, d) in x.items()}, index=index
        )
    if isinstance(x, tuple) and len(x) == 2:
        return _empty_series(x[0], x[1], index=index)
    elif isinstance(x, Iterable) and not isinstance(x, str):
        if not all(isinstance(i, tuple) and len(i) == 2 for i in x):
            raise ValueError(f"Expected iterable of tuples of (name, dtype), got {x}")
        return pd.DataFrame(
            {c: _empty_series(c, d, index=index) for (c, d) in x},
            columns=[c for c, d in x],
            index=index,
        )
    elif not hasattr(x, "dtype") and x is not None:
        # could be a string, a dtype object, or a python type. Skip `None`,
        # because it is implicitly converted to `dtype('f8')`, which we don't
        # want here.
        try:
            dtype = np.dtype(x)
            return _scalar_from_dtype(dtype)
        except Exception:
            # Continue on to next check
            pass

    if is_scalar(x):
        return _nonempty_scalar(x)

    raise TypeError(f"Don't know how to create metadata from {x}")


@meta_nonempty.register(object)
def meta_nonempty_object(x):
    """Create a nonempty pandas object from the given metadata.

    Returns a pandas DataFrame, Series, or Index that contains two rows
    of fake data.
    """
    if is_scalar(x):
        return _nonempty_scalar(x)
    else:
        raise TypeError(
            "Expected Pandas-like Index, Series, DataFrame, or scalar, "
            f"got {typename(type(x))}"
        )


@meta_nonempty.register(pd.DataFrame)
def meta_nonempty_dataframe(x):
    idx = meta_nonempty(x.index)
    dt_s_dict = dict()
    data = dict()
    for i in range(len(x.columns)):
        series = x.iloc[:, i]
        dt = series.dtype
        if dt not in dt_s_dict:
            dt_s_dict[dt] = _nonempty_series(x.iloc[:, i], idx=idx)
        data[i] = dt_s_dict[dt]
    res = pd.DataFrame(data, index=idx, columns=np.arange(len(x.columns)))
    res.columns = x.columns
    res.attrs = x.attrs
    return res


@meta_nonempty.register(pd.Index)
def _nonempty_index(idx):
    typ = type(idx)
    if typ is pd.RangeIndex:
        return pd.RangeIndex(2, name=idx.name, dtype=idx.dtype)
    elif is_any_real_numeric_dtype(idx):
        return typ([1, 2], name=idx.name, dtype=idx.dtype)
    elif typ is pd.DatetimeIndex:
        start = "1970-01-01"
        # Need a non-monotonic decreasing index to avoid issues with
        # partial string indexing see https://github.com/dask/dask/issues/2389
        # and https://github.com/pandas-dev/pandas/issues/16515
        # This doesn't mean `_meta_nonempty` should ever rely on
        # `self.monotonic_increasing` or `self.monotonic_decreasing`
        try:
            return pd.date_range(
                start=start,
                periods=2,
                freq=idx.freq,
                tz=idx.tz,
                name=idx.name,
                unit=idx.unit,
            )
        except ValueError:  # older pandas versions
            data = [start, "1970-01-02"] if idx.freq is None else None
            return pd.DatetimeIndex(
                data, start=start, periods=2, freq=idx.freq, tz=idx.tz, name=idx.name
            )
    elif typ is pd.PeriodIndex:
        return pd.period_range(
            start="1970-01-01", periods=2, freq=idx.freq, name=idx.name
        )
    elif typ is pd.TimedeltaIndex:
        start = np.timedelta64(1, "D")
        try:
            return pd.timedelta_range(
                start=start, periods=2, freq=idx.freq, name=idx.name
            )
        except ValueError:  # older pandas versions
            start = np.timedelta64(1, "D")
            data = [start, start + 1] if idx.freq is None else None
            return pd.TimedeltaIndex(
                data, start=start, periods=2, freq=idx.freq, name=idx.name
            )
    elif typ is pd.CategoricalIndex:
        if len(idx.categories) == 0:
            data = pd.Categorical(_nonempty_index(idx.categories), ordered=idx.ordered)
        else:
            data = pd.Categorical.from_codes(
                [-1, 0], categories=idx.categories, ordered=idx.ordered
            )
        return pd.CategoricalIndex(data, name=idx.name)
    elif typ is pd.MultiIndex:
        levels = [_nonempty_index(l) for l in idx.levels]
        codes = [[0, 0] for i in idx.levels]
        try:
            return pd.MultiIndex(levels=levels, codes=codes, names=idx.names)
        except TypeError:  # older pandas versions
            return pd.MultiIndex(levels=levels, labels=codes, names=idx.names)
    elif typ is pd.Index:
        if type(idx.dtype) in make_array_nonempty._lookup:
            return pd.Index(
                make_array_nonempty(idx.dtype), dtype=idx.dtype, name=idx.name
            )
        elif idx.dtype == bool:
            # pd 1.5 introduce bool dtypes and respect non-uniqueness
            return pd.Index([True, False], name=idx.name)
        else:
            # for pd 1.5 in the case of bool index this would be cast as [True, True]
            # breaking uniqueness
            return pd.Index(["a", "b"], name=idx.name, dtype=idx.dtype)

    raise TypeError(f"Don't know how to handle index of type {typename(type(idx))}")


@meta_nonempty.register(pd.Series)
def _nonempty_series(s, idx=None):
    # TODO: Use register dtypes with make_array_nonempty
    if idx is None:
        idx = _nonempty_index(s.index)
    dtype = s.dtype
    if len(s) > 0:
        # use value from meta if provided
        data = [s.iloc[0]] * 2
    elif isinstance(dtype, pd.DatetimeTZDtype):
        entry = pd.Timestamp("1970-01-01", tz=dtype.tz)
        data = pd.array([entry, entry], dtype=dtype)
    elif isinstance(dtype, pd.CategoricalDtype):
        if len(s.cat.categories):
            data = [s.cat.categories[0]] * 2
            cats = s.cat.categories
        else:
            data = _nonempty_index(s.cat.categories)
            cats = s.cat.categories[:0]
        data = pd.Categorical(data, categories=cats, ordered=s.cat.ordered)
    elif is_integer_na_dtype(dtype):
        data = pd.array([1, None], dtype=dtype)
    elif is_float_na_dtype(dtype):
        data = pd.array([1.0, None], dtype=dtype)
    elif isinstance(dtype, pd.PeriodDtype):
        # pandas 0.24.0+ should infer this to be Series[Period[freq]]
        freq = dtype.freq
        data = [pd.Period("2000", freq), pd.Period("2001", freq)]
    elif isinstance(dtype, pd.SparseDtype):
        entry = _scalar_from_dtype(dtype.subtype)
        data = pd.array([entry, entry], dtype=dtype)
    elif isinstance(dtype, pd.IntervalDtype):
        entry = _scalar_from_dtype(dtype.subtype)
        data = pd.array([entry, entry], dtype=dtype)
    elif type(dtype) in make_array_nonempty._lookup:
        data = make_array_nonempty(dtype)
    else:
        entry = _scalar_from_dtype(dtype)
        data = np.array([entry, entry], dtype=dtype)

    out = pd.Series(data, name=s.name, index=idx)
    out.attrs = s.attrs
    return out


@meta_lib_from_array.register(Array)
def _meta_lib_from_array_da(x):
    # Use x._meta for dask arrays
    return meta_lib_from_array(x._meta)


@meta_lib_from_array.register(np.ndarray)
def _meta_lib_from_array_numpy(x):
    # numpy -> pandas
    return pd


@union_categoricals_dispatch.register(
    (pd.DataFrame, pd.Series, pd.Index, pd.Categorical)
)
def union_categoricals_pandas(to_union, sort_categories=False, ignore_order=False):
    return pd.api.types.union_categoricals(
        to_union, sort_categories=sort_categories, ignore_order=ignore_order
    )


@get_parallel_type.register(pd.Series)
def get_parallel_type_series(_):
    from dask.dataframe.dask_expr._collection import Series

    return Series


@get_parallel_type.register(pd.DataFrame)
def get_parallel_type_dataframe(_):
    from dask.dataframe.dask_expr._collection import DataFrame

    return DataFrame


@get_parallel_type.register(pd.Index)
def get_parallel_type_index(_):
    from dask.dataframe.dask_expr._collection import Index

    return Index


@get_parallel_type.register(object)
def get_parallel_type_object(_):
    from dask.dataframe.dask_expr._collection import Scalar

    return Scalar


@hash_object_dispatch.register((pd.DataFrame, pd.Series, pd.Index))
def hash_object_pandas(
    obj, index=True, encoding="utf8", hash_key=None, categorize=True
):
    return pd.util.hash_pandas_object(
        obj, index=index, encoding=encoding, hash_key=hash_key, categorize=categorize
    )


class ShuffleGroupResult(SimpleSizeof, dict):
    def __sizeof__(self) -> int:
        """
        The result of the shuffle split are typically small dictionaries
        (#keys << 100; typically <= 32) The splits are often non-uniformly
        distributed. Some of the splits may even be empty. Sampling the
        dictionary for size estimation can cause severe errors.

        See also https://github.com/dask/distributed/issues/4962
        """
        total_size = super().__sizeof__()
        for k, df in self.items():
            total_size += sizeof(k)
            total_size += sizeof(df)
        return total_size


@group_split_dispatch.register((pd.DataFrame, pd.Series, pd.Index))
def group_split_pandas(df, c, k, ignore_index=False):
    if is_series_like(c):
        c = c.values
    indexer, locations = pd._libs.algos.groupsort_indexer(
        c.astype(np.intp, copy=False), k
    )
    df2 = df.take(indexer)
    locations = locations.cumsum()
    parts = [
        df2.iloc[a:b].reset_index(drop=True) if ignore_index else df2.iloc[a:b]
        for a, b in zip(locations[:-1], locations[1:])
    ]
    return ShuffleGroupResult(zip(range(k), parts))


@concat_dispatch.register((pd.DataFrame, pd.Series, pd.Index))
def concat_pandas(
    dfs,
    axis=0,
    join="outer",
    uniform=False,
    filter_warning=True,
    ignore_index=False,
    **kwargs,
):
    ignore_order = kwargs.pop("ignore_order", False)

    if axis == 1:
        return pd.concat(dfs, axis=axis, join=join, **kwargs)

    # Support concatenating indices along axis 0
    if isinstance(dfs[0], pd.Index):
        if isinstance(dfs[0], pd.CategoricalIndex):
            for i in range(1, len(dfs)):
                if not isinstance(dfs[i], pd.CategoricalIndex):
                    dfs[i] = dfs[i].astype("category")
            return pd.CategoricalIndex(
                union_categoricals(dfs, ignore_order=ignore_order), name=dfs[0].name
            )
        elif isinstance(dfs[0], pd.MultiIndex):
            first, rest = dfs[0], dfs[1:]
            if all(
                (isinstance(o, pd.MultiIndex) and o.nlevels >= first.nlevels)
                for o in rest
            ):
                arrays = [
                    concat([i._get_level_values(n) for i in dfs])
                    for n in range(first.nlevels)
                ]
                return pd.MultiIndex.from_arrays(arrays, names=first.names)

            to_concat = (first.values,) + tuple(k._values for k in rest)
            new_tuples = np.concatenate(to_concat)
            try:
                return pd.MultiIndex.from_tuples(new_tuples, names=first.names)
            except Exception:
                return pd.Index(new_tuples)
        return dfs[0].append(dfs[1:])

    # Handle categorical index separately
    dfs0_index = dfs[0].index

    has_categoricalindex = isinstance(dfs0_index, pd.CategoricalIndex) or (
        isinstance(dfs0_index, pd.MultiIndex)
        and any(isinstance(i, pd.CategoricalIndex) for i in dfs0_index.levels)
    )

    if has_categoricalindex:
        dfs2 = [df.reset_index(drop=True) for df in dfs]
        ind = concat([df.index for df in dfs])
    else:
        dfs2 = dfs
        ind = None

    # Concatenate the partitions together, handling categories as needed
    if (
        isinstance(dfs2[0], pd.DataFrame)
        if uniform
        else any(isinstance(df, pd.DataFrame) for df in dfs2)
    ):
        if uniform or PANDAS_GE_220:
            dfs3 = dfs2
            cat_mask = dfs2[0].dtypes == "category"
        else:
            # When concatenating mixed dataframes and series on axis 1, Pandas <2.2
            # converts series to dataframes with a single column named 0, then
            # concatenates.
            dfs3 = [
                (
                    df
                    if isinstance(df, pd.DataFrame)
                    else df.to_frame().rename(columns={df.name: 0})
                )
                for df in dfs2
            ]
            # pandas may raise a RuntimeWarning for comparing ints and strs
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                if filter_warning:
                    warnings.simplefilter("ignore", FutureWarning)
                cat_mask = pd.concat(
                    [(df.dtypes == "category").to_frame().T for df in dfs3],
                    join=join,
                    **kwargs,
                ).any()

        if isinstance(cat_mask, pd.Series) and cat_mask.any():
            not_cat = cat_mask[~cat_mask].index
            # this should be aligned, so no need to filter warning
            out = pd.concat(
                [df[df.columns.intersection(not_cat)] for df in dfs3],
                join=join,
                **kwargs,
            )
            temp_ind = out.index
            for col in cat_mask.index.difference(not_cat):
                # Find an example of categoricals in this column
                for df in dfs3:
                    sample = df.get(col)
                    if sample is not None:
                        break
                # Extract partitions, subbing in missing if needed
                parts = []
                for df in dfs3:
                    if col in df.columns:
                        parts.append(df[col])
                    else:
                        codes = np.full(len(df), -1, dtype="i8")
                        data = pd.Categorical.from_codes(
                            codes, sample.cat.categories, sample.cat.ordered
                        )
                        parts.append(data)
                out[col] = union_categoricals(parts, ignore_order=ignore_order)
                # Pandas resets index type on assignment if frame is empty
                # https://github.com/pandas-dev/pandas/issues/17101
                if not len(temp_ind):
                    out.index = temp_ind
            out = out.reindex(columns=cat_mask.index)
        else:
            # pandas may raise a RuntimeWarning for comparing ints and strs
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                if filter_warning:
                    warnings.simplefilter("ignore", FutureWarning)
                out = pd.concat(dfs3, join=join, sort=False)
    else:
        if isinstance(dfs2[0].dtype, pd.CategoricalDtype):
            if ind is None:
                ind = concat([df.index for df in dfs2])
            return pd.Series(
                union_categoricals(dfs2, ignore_order=ignore_order),
                index=ind,
                name=dfs2[0].name,
            )
        with warnings.catch_warnings():
            if filter_warning:
                warnings.simplefilter("ignore", FutureWarning)

            out = pd.concat(dfs2, join=join, **kwargs)
    # Re-add the index if needed
    if ind is not None:
        out.index = ind
    return out


@categorical_dtype_dispatch.register((pd.DataFrame, pd.Series, pd.Index))
def categorical_dtype_pandas(categories=None, ordered=False):
    return pd.api.types.CategoricalDtype(categories=categories, ordered=ordered)


@tolist_dispatch.register((np.ndarray, pd.Series, pd.Index, pd.Categorical))
def tolist_numpy_or_pandas(obj):
    return obj.tolist()


@is_categorical_dtype_dispatch.register(
    (pd.Series, pd.Index, pd.api.extensions.ExtensionDtype, np.dtype)
)
def is_categorical_dtype_pandas(obj):
    if hasattr(obj, "dtype"):
        dtype = obj.dtype
    else:
        dtype = obj
    return isinstance(dtype, pd.CategoricalDtype)


@grouper_dispatch.register((pd.DataFrame, pd.Series))
def get_grouper_pandas(obj):
    return pd.core.groupby.Grouper


@percentile_lookup.register((pd.Series, pd.Index))
def percentile(a, q, interpolation="linear"):
    if isinstance(a.dtype, pd.ArrowDtype):
        a = a.to_numpy()
    return _percentile(a, q, interpolation)


@to_pandas_dispatch.register((pd.DataFrame, pd.Series, pd.Index))
def to_pandas_dispatch_from_pandas(data, **kwargs):
    return data


class PandasBackendEntrypoint(DataFrameBackendEntrypoint):
    """Pandas-Backend Entrypoint Class for Dask-DataFrame

    Note that all DataFrame-creation functions are defined
    and registered 'in-place' within the ``dask.dataframe``
    ``io`` module.
    """

    @classmethod
    def to_backend_dispatch(cls):
        return to_pandas_dispatch

    @classmethod
    def to_backend(cls, data, **kwargs):
        if isinstance(data._meta, (pd.DataFrame, pd.Series, pd.Index)):
            # Already a pandas-backed collection
            return data
        return data.map_partitions(cls.to_backend_dispatch(), **kwargs)


dataframe_creation_dispatch.register_backend("pandas", PandasBackendEntrypoint())


######################################
# cuDF: Pandas Dataframes on the GPU #
######################################


@concat_dispatch.register_lazy("cudf")
@from_pyarrow_table_dispatch.register_lazy("cudf")
@group_split_dispatch.register_lazy("cudf")
@get_parallel_type.register_lazy("cudf")
@hash_object_dispatch.register_lazy("cudf")
@meta_nonempty.register_lazy("cudf")
@make_meta_dispatch.register_lazy("cudf")
@make_meta_obj.register_lazy("cudf")
@percentile_lookup.register_lazy("cudf")
@to_pandas_dispatch.register_lazy("cudf")
@to_pyarrow_table_dispatch.register_lazy("cudf")
@tolist_dispatch.register_lazy("cudf")
def _register_cudf():
    import dask_cudf  # noqa: F401


@meta_lib_from_array.register_lazy("cupy")
@tolist_dispatch.register_lazy("cupy")
def _register_cupy_to_cudf():
    # Handle cupy.ndarray -> cudf.DataFrame dispatching
    try:
        import cudf
        import cupy

        @meta_lib_from_array.register(cupy.ndarray)
        def meta_lib_from_array_cupy(x):
            # cupy -> cudf
            return cudf

        @tolist_dispatch.register(cupy.ndarray)
        def tolist_cupy(x):
            return x.tolist()

    except ImportError:
        pass
