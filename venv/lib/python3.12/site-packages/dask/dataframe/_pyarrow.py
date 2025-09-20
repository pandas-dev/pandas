from __future__ import annotations

from functools import partial

import pandas as pd

from dask._compatibility import import_optional_dependency
from dask.dataframe.utils import is_dataframe_like, is_index_like, is_series_like

pa = import_optional_dependency("pyarrow")


def is_pyarrow_string_dtype(dtype):
    """Is the input dtype a pyarrow string?"""
    if pa is None:
        return False
    return dtype in (pd.StringDtype("pyarrow"), pd.ArrowDtype(pa.string()))


def is_object_string_dtype(dtype):
    """Determine if input is a non-pyarrow string dtype"""
    # in pandas < 2.0, is_string_dtype(DecimalDtype()) returns True
    return (
        pd.api.types.is_string_dtype(dtype)
        and not is_pyarrow_string_dtype(dtype)
        and not pd.api.types.is_dtype_equal(dtype, "decimal")
    )


def is_pyarrow_string_index(x):
    if isinstance(x, pd.MultiIndex):
        return any(is_pyarrow_string_index(level) for level in x.levels)
    return isinstance(x, pd.Index) and is_pyarrow_string_dtype(x.dtype)


def is_object_string_index(x):
    if isinstance(x, pd.MultiIndex):
        return any(is_object_string_index(level) for level in x.levels)
    return isinstance(x, pd.Index) and is_object_string_dtype(x.dtype)


def is_object_string_series(x):
    return isinstance(x, pd.Series) and (
        is_object_string_dtype(x.dtype) or is_object_string_index(x.index)
    )


def is_object_string_dataframe(x):
    return isinstance(x, pd.DataFrame) and (
        any(is_object_string_series(s) for _, s in x.items())
        or is_object_string_index(x.index)
    )


def _to_string_dtype(df, dtype_check, index_check, string_dtype):
    if not (is_dataframe_like(df) or is_series_like(df) or is_index_like(df)):
        return df

    # Guards against importing `pyarrow` at the module level (where it may not be installed)
    if string_dtype == "pyarrow":
        string_dtype = pd.StringDtype("pyarrow")

    # Possibly convert DataFrame/Series/Index to `string[pyarrow]`
    if is_dataframe_like(df):
        dtypes = {
            col: string_dtype for col, dtype in df.dtypes.items() if dtype_check(dtype)
        }
        if dtypes:
            df = df.astype(dtypes)
    elif dtype_check(df.dtype):
        dtypes = string_dtype
        df = df.copy().astype(dtypes)

    # Convert DataFrame/Series index too
    if (is_dataframe_like(df) or is_series_like(df)) and index_check(df.index):
        if isinstance(df.index, pd.MultiIndex):
            levels = {
                i: level.astype(string_dtype)
                for i, level in enumerate(df.index.levels)
                if dtype_check(level.dtype)
            }
            # set verify_integrity=False to preserve index codes
            df.index = df.index.set_levels(
                levels.values(), level=levels.keys(), verify_integrity=False
            )
        else:
            df.index = df.index.astype(string_dtype)
    return df


to_pyarrow_string = partial(
    _to_string_dtype,
    dtype_check=is_object_string_dtype,
    index_check=is_object_string_index,
    string_dtype="pyarrow",
)
to_object_string = partial(
    _to_string_dtype,
    dtype_check=is_pyarrow_string_dtype,
    index_check=is_pyarrow_string_index,
    string_dtype=object,
)
