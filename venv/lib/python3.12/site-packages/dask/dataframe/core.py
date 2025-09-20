from __future__ import annotations

import warnings
from collections.abc import Iterator

import numpy as np
import pandas as pd
from tlz import first, unique

import dask.array as da
from dask import core
from dask.base import named_schedulers
from dask.dataframe import methods
from dask.dataframe._compat import PANDAS_GE_300
from dask.dataframe.dispatch import get_parallel_type
from dask.dataframe.utils import (
    check_matching_columns,
    has_known_categories,
    is_dataframe_like,
    is_index_like,
    is_scalar,
    is_series_like,
    meta_frame_constructor,
    meta_series_constructor,
    valid_divisions,
)
from dask.typing import no_default
from dask.utils import M

DEFAULT_GET = named_schedulers.get("threads", named_schedulers["sync"])

pd.set_option("compute.use_numexpr", False)


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


def split_evenly(df, k):
    """Split dataframe into k roughly equal parts"""
    divisions = np.linspace(0, len(df), k + 1).astype(int)
    return {i: df.iloc[divisions[i] : divisions[i + 1]] for i in range(k)}


def _get_divisions_map_partitions(
    align_dataframes, transform_divisions, dfs, func, args, kwargs
):
    """
    Helper to get divisions for map_partitions and map_overlap output.
    """
    from dask.dataframe import Index

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


def apply_and_enforce(*args, **kwargs):
    """Apply a function, and enforce the output to match meta

    Ensures the output has the same columns, even if empty."""
    func = kwargs.pop("_func")
    meta = kwargs.pop("_meta")
    df = func(*args, **kwargs)

    if any(
        bool(is_dataframe_like(obj) or is_series_like(obj) or is_index_like(obj))
        for obj in [df, meta]
    ):
        if is_scalar(df):
            df = pd.Series(df)
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


def total_mem_usage(df, index=True, deep=False):
    mem_usage = df.memory_usage(index=index, deep=deep)
    if is_series_like(mem_usage):
        mem_usage = mem_usage.sum()
    return mem_usage


def idxmaxmin_chunk(x, fn=None, skipna=True, numeric_only=False):
    numeric_only_kwargs = {} if is_series_like(x) else {"numeric_only": numeric_only}
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


def safe_head(df, n):
    r = M.head(df, n)
    if len(r) != n:
        warnings.warn(
            f"Insufficient elements for `head`. {n} elements requested, only {len(r)} "
            "elements available. Try passing larger `npartitions` to `head`."
        )
    return r


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
    from dask.dataframe.dask_expr._collection import Scalar

    return get_parallel_type(x) is not Scalar


def meta_warning(df, method="apply"):
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
        f"please provide the `meta=` keyword, as described in the {method} function "
        f"that you are using."
    )
    if meta_str:
        msg += (
            "\n"
            f"  Before: .{method}(func)\n"
            f"  After:  .{method}(func, meta=%s)\n" % str(meta_str)
        )
    return msg


def _convert_to_numeric(series, skipna):
    if skipna:
        return series.dropna().astype("i8")

    # series.view("i8") with pd.NaT produces -9223372036854775808 is why we need to do this
    return series.astype("i8").mask(series.isnull(), np.nan)


def _sqrt_and_convert_to_timedelta(partition, axis, dtype=None, *args, **kwargs):
    if axis == 1:
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                category=RuntimeWarning,
                message="invalid value encountered in cast",
            )
            unit = kwargs.pop("unit", None)
            result = pd.to_timedelta(
                M.std(partition, *args, axis=axis, **kwargs), unit=unit
            )
            if unit is not None and dtype is not None:
                result = result.astype(dtype)
            return result

    is_df_like, time_cols = kwargs["is_df_like"], kwargs["time_cols"]

    sqrt = np.sqrt(partition)

    if not is_df_like:
        result = pd.to_timedelta(sqrt, unit=kwargs.get("unit"))
        if kwargs.get("unit") is not None:
            result = result.as_unit(kwargs["unit"])
        return result

    time_col_mask = sqrt.index.isin(time_cols)
    matching_vals = sqrt[time_col_mask]
    if len(time_cols) > 0:
        sqrt = sqrt.astype(object)

    units = kwargs.get("units")
    if units is None:
        units = [None] * len(time_cols)
    for time_col, matching_val, unit in zip(time_cols, matching_vals, units):
        result = pd.to_timedelta(matching_val, unit=unit)
        if unit is not None:
            result = result.as_unit(unit)
        sqrt[time_col] = result

    if dtype is not None:
        sqrt = sqrt.astype(dtype)
    return sqrt
