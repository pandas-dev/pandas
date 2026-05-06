from __future__ import annotations

import datetime as dt
from typing import (
    TYPE_CHECKING,
    Literal,
    cast,
)
import zoneinfo

import numpy as np

from pandas._config import using_string_dtype

from pandas._libs import lib
from pandas._libs.tslibs import timezones
from pandas.compat import (
    pa_version_under18p0,
    pa_version_under19p0,
)
from pandas.compat._optional import import_optional_dependency

from pandas.core.dtypes.common import pandas_dtype

import pandas as pd

if TYPE_CHECKING:
    from collections.abc import (
        Callable,
        Hashable,
        Sequence,
    )

    import pyarrow

    from pandas._typing import (
        DtypeArg,
        DtypeBackend,
    )


pytz = import_optional_dependency("pytz", errors="ignore")


def _arrow_dtype_mapping() -> dict:
    pa = import_optional_dependency("pyarrow")
    return {
        pa.int8(): pd.Int8Dtype(),
        pa.int16(): pd.Int16Dtype(),
        pa.int32(): pd.Int32Dtype(),
        pa.int64(): pd.Int64Dtype(),
        pa.uint8(): pd.UInt8Dtype(),
        pa.uint16(): pd.UInt16Dtype(),
        pa.uint32(): pd.UInt32Dtype(),
        pa.uint64(): pd.UInt64Dtype(),
        pa.bool_(): pd.BooleanDtype(),
        pa.string(): pd.StringDtype(),
        pa.float32(): pd.Float32Dtype(),
        pa.float64(): pd.Float64Dtype(),
        pa.string(): pd.StringDtype(),
        pa.large_string(): pd.StringDtype(),
    }


def _arrow_string_types_mapper() -> Callable:
    pa = import_optional_dependency("pyarrow")

    mapping = {
        pa.string(): pd.StringDtype(na_value=np.nan),
        pa.large_string(): pd.StringDtype(na_value=np.nan),
    }
    if not pa_version_under18p0:
        mapping[pa.string_view()] = pd.StringDtype(na_value=np.nan)

    return mapping.get


def arrow_table_to_pandas(
    table: pyarrow.Table,
    dtype_backend: DtypeBackend | Literal["numpy"] | lib.NoDefault = lib.no_default,
    null_to_int64: bool = False,
    to_pandas_kwargs: dict | None = None,
    dtype: DtypeArg | None = None,
    names: Sequence[Hashable] | None = None,
) -> pd.DataFrame:
    pa = import_optional_dependency("pyarrow")

    to_pandas_kwargs = {} if to_pandas_kwargs is None else to_pandas_kwargs

    types_mapper: type[pd.ArrowDtype] | None | Callable
    if dtype_backend == "numpy_nullable":
        mapping = _arrow_dtype_mapping()
        if null_to_int64:
            # Modify the default mapping to also map null to Int64
            # (to match other engines - only for CSV parser)
            mapping[pa.null()] = pd.Int64Dtype()
        types_mapper = mapping.get
    elif dtype_backend == "pyarrow":
        types_mapper = pd.ArrowDtype
    elif using_string_dtype():
        if pa_version_under19p0:
            types_mapper = _arrow_string_types_mapper()
        elif dtype is not None:
            # GH#56136 Avoid lossy conversion to float64
            # We'll convert to numpy below if
            types_mapper = {
                pa.int8(): pd.Int8Dtype(),
                pa.int16(): pd.Int16Dtype(),
                pa.int32(): pd.Int32Dtype(),
                pa.int64(): pd.Int64Dtype(),
            }.get
        else:
            types_mapper = None
    elif dtype_backend is lib.no_default or dtype_backend == "numpy":
        if dtype is not None:
            # GH#56136 Avoid lossy conversion to float64
            # We'll convert to numpy below if
            types_mapper = {
                pa.int8(): pd.Int8Dtype(),
                pa.int16(): pd.Int16Dtype(),
                pa.int32(): pd.Int32Dtype(),
                pa.int64(): pd.Int64Dtype(),
            }.get
        else:
            types_mapper = None
    else:
        raise NotImplementedError

    df = table.to_pandas(types_mapper=types_mapper, **to_pandas_kwargs)
    df = _post_convert_dtypes(df, dtype_backend, dtype, names)
    df = _normalize_timezone_dtypes(df)
    return df


def _post_convert_dtypes(
    df: pd.DataFrame,
    dtype_backend: DtypeBackend | Literal["numpy"] | lib.NoDefault,
    dtype: DtypeArg | None,
    names: Sequence[Hashable] | None,
) -> pd.DataFrame:
    if dtype is not None and (
        dtype_backend is lib.no_default or dtype_backend == "numpy"
    ):
        # GH#56136 apply any user-provided dtype, and convert any IntegerDtype
        #  columns the user didn't explicitly ask for.
        if isinstance(dtype, dict):
            if names is not None:
                df.columns = names

            cmp_dtypes = {
                pd.Int8Dtype(),
                pd.Int16Dtype(),
                pd.Int32Dtype(),
                pd.Int64Dtype(),
            }
            for col in df.columns:
                if col not in dtype and df[col].dtype in cmp_dtypes:
                    # Any key that the user didn't explicitly specify
                    #  that got converted to IntegerDtype now gets converted
                    #  to numpy dtype.
                    dtype[col] = df[col].dtype.numpy_dtype

            # Ignore non-existent columns from dtype mapping
            # like other parsers do
            dtype = {
                key: pandas_dtype(dtype[key]) for key in dtype if key in df.columns
            }

        else:
            dtype = pandas_dtype(dtype)

        try:
            df = df.astype(dtype)
        except TypeError as err:
            # GH#44901 reraise to keep api consistent
            raise ValueError(str(err)) from err

    if (
        not using_string_dtype()
        and dtype != "str"
        and (dtype_backend is lib.no_default or dtype_backend == "numpy")
    ):
        # Convert any StringDtype columns back to object dtype (pyarrow always
        # uses string dtype even when the infer_string option is False)
        for i in range(len(df.columns)):
            new_col = _maybe_convert_string_to_object(df.iloc[:, i])
            if new_col is not None:
                df.isetitem(i, new_col)

        new_idx = _maybe_convert_string_index_to_object(df.index)
        if new_idx is not None:
            df.index = new_idx
        new_cols = _maybe_convert_string_index_to_object(df.columns)
        if new_cols is not None:
            df.columns = new_cols

    return df


def _maybe_convert_string_to_object(
    data: pd.Series | pd.Index,
) -> pd.Series | pd.Index | None:
    if isinstance(data.dtype, pd.StringDtype) and data.dtype.na_value is np.nan:
        return data.astype("object").fillna(None)
    elif isinstance(data.dtype, pd.CategoricalDtype):
        cat_dtype = data.dtype.categories.dtype
        if isinstance(cat_dtype, pd.StringDtype) and cat_dtype.na_value is np.nan:
            cat_dtype = pd.CategoricalDtype(
                categories=data.dtype.categories.astype("object"),
                ordered=data.dtype.ordered,
            )
            return data.astype(cat_dtype)

    # no conversion needed
    return None


def _maybe_convert_string_index_to_object(index: pd.Index) -> pd.Index | None:
    if isinstance(index, pd.MultiIndex):
        if any(
            isinstance(level.dtype, pd.StringDtype) and level.dtype.na_value is np.nan
            for level in index.levels
        ):
            new_levels = []
            for level in index.levels:
                new_level = _maybe_convert_string_to_object(level)
                if new_level is not None:
                    new_levels.append(new_level)
                else:
                    new_levels.append(level)
            return index.set_levels(new_levels)
        return None

    else:
        return cast("pd.Index | None", _maybe_convert_string_to_object(index))


def _normalize_pytz_timezone(tz: dt.tzinfo) -> dt.tzinfo:
    """
    If the input tz is a pytz timezone, attempt to convert it to "default"
    tzinfo object (zoneinfo or datetime.timezone).
    """
    if not type(tz).__module__.startswith("pytz"):
        # isinstance(col.dtype.tz, pytz.BaseTzInfo) does not included
        # fixed offsets
        return tz

    if timezones.is_utc(tz):
        return dt.timezone.utc

    if tz.zone is not None:  # type: ignore[attr-defined]
        try:
            return zoneinfo.ZoneInfo(tz.zone)  # type: ignore[attr-defined]
        except Exception:
            # some pytz timezones might not be available for zoneinfo
            pass

    if timezones.is_fixed_offset(tz):
        # Convert pytz fixed offset to datetime.timezone
        try:
            offset = tz.utcoffset(None)
            if offset is not None:
                return dt.timezone(offset)
        except Exception:
            pass

    return tz


def _normalize_timezone_index(index: pd.Index) -> pd.Index:
    if isinstance(index, pd.MultiIndex):
        if any(isinstance(level.dtype, pd.DatetimeTZDtype) for level in index.levels):
            levels = [_normalize_timezone_index(level) for level in index.levels]
            return index.set_levels(levels)

        return index

    if isinstance(index.dtype, pd.DatetimeTZDtype):
        normalized_tz = _normalize_pytz_timezone(index.dtype.tz)
        if normalized_tz is not index.dtype.tz:
            return index.tz_convert(normalized_tz)  # type: ignore[attr-defined]

    return index


def _normalize_timezone_dtypes(df: pd.DataFrame) -> pd.DataFrame:
    """
    PyArrow uses pytz by default for timezones, but pandas uses
    zoneinfo / datetime.timezone since pandas 3.0.

    TODO: Starting with pyarrow 25, it will use zoneinfo by default, and then
    this normalization can be skipped (https://github.com/apache/arrow/pull/49694).
    """
    if pytz is not None:
        # Convert any pytz timezones to zoneinfo / fixed offset timezones
        if any(
            isinstance(dtype, pd.DatetimeTZDtype)
            for dtype in df._mgr.get_unique_dtypes()
        ):
            col_indices = df._select_dtypes_indices(pd.DatetimeTZDtype)
            for i in col_indices:
                col = df.iloc[:, i]
                normalized_tz = _normalize_pytz_timezone(col.dtype.tz)
                if normalized_tz is not col.dtype.tz:
                    df.isetitem(i, col.dt.tz_convert(normalized_tz))

        df.index = _normalize_timezone_index(df.index)
        df.columns = _normalize_timezone_index(df.columns)

    return df
