from __future__ import annotations

import contextlib
import datetime as dt
from typing import (
    TYPE_CHECKING,
    Any,
    cast,
)
import warnings
import zoneinfo

import numpy as np

from pandas._config import using_string_dtype

from pandas._libs import lib
from pandas._libs.tslibs import timezones
from pandas.compat import (
    pa_version_under18p0,
    pa_version_under19p0,
    pa_version_under25p0,
)
from pandas.compat._optional import import_optional_dependency
from pandas.errors import Pandas4Warning

from pandas.core.dtypes.common import pandas_dtype
from pandas.core.dtypes.dtypes import (
    ArrowDtype,
    BaseMaskedDtype,
    CategoricalDtype,
    DatetimeTZDtype,
    SparseDtype,
)
from pandas.core.dtypes.missing import isna

from pandas.core.arrays.arrow.util import _arrow_dtype_mapping
from pandas.core.arrays.integer import (
    Int8Dtype,
    Int16Dtype,
    Int32Dtype,
    Int64Dtype,
    IntegerDtype,
)
from pandas.core.arrays.string_ import StringDtype
from pandas.core.indexes.api import (
    CategoricalIndex,
    DatetimeIndex,
    Index,
    MultiIndex,
)
from pandas.core.tools.datetimes import to_datetime

if TYPE_CHECKING:
    from collections.abc import (
        Callable,
        Generator,
        Hashable,
        Sequence,
    )

    import pyarrow

    from pandas._typing import (
        DtypeArg,
        DtypeBackend,
        DtypeObj,
    )

    from pandas.core.dtypes.base import ExtensionDtype

    from pandas import (
        DataFrame,
        Series,
    )
    from pandas.core.arrays import Categorical


pytz = import_optional_dependency("pytz", errors="ignore")


@contextlib.contextmanager
def suppress_pyarrow_values_warning() -> Generator[None]:
    """
    Suppress the deprecation warning pyarrow triggers by calling ``.values``
    on timezone-aware data when converting pandas objects.

    Remove once the minimum pyarrow (26.0) no longer does this, see GH#68426 and
    apache/arrow#51302.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            r"(Series|DatetimeIndex)\.values returning an ndarray that drops "
            "timezone information",
            Pandas4Warning,
        )
        yield


def _arrow_string_types_mapper() -> Callable[[pyarrow.DataType], ExtensionDtype | None]:
    pa = import_optional_dependency("pyarrow")

    mapping = {
        pa.string(): StringDtype(na_value=np.nan),
        pa.large_string(): StringDtype(na_value=np.nan),
    }
    if not pa_version_under18p0:
        mapping[pa.string_view()] = StringDtype(na_value=np.nan)

    return mapping.get


def arrow_table_to_pandas(
    table: pyarrow.Table,
    dtype_backend: DtypeBackend | lib.NoDefault = lib.no_default,
    null_to_int64: bool = False,
    to_pandas_kwargs: dict[str, Any] | None = None,
    dtype: DtypeArg | None = None,
    names: Sequence[Hashable] | None = None,
) -> DataFrame:
    pa = import_optional_dependency("pyarrow")

    to_pandas_kwargs = {} if to_pandas_kwargs is None else to_pandas_kwargs

    types_mapper: (
        type[ArrowDtype] | Callable[[pyarrow.DataType], ExtensionDtype | None] | None
    )
    if dtype_backend == "numpy_nullable":
        mapping = _arrow_dtype_mapping()
        if null_to_int64:
            # Modify the default mapping to also map null to Int64
            # (to match other engines - only for CSV parser)
            mapping[pa.null()] = Int64Dtype()
        types_mapper = mapping.get
    elif dtype_backend == "pyarrow":
        types_mapper = ArrowDtype
    elif using_string_dtype():
        if pa_version_under19p0:
            types_mapper = _arrow_string_types_mapper()
        elif dtype is not None:
            # GH#56136 Avoid lossy conversion to float64
            # We'll convert to numpy below if
            types_mapper = {
                pa.int8(): Int8Dtype(),
                pa.int16(): Int16Dtype(),
                pa.int32(): Int32Dtype(),
                pa.int64(): Int64Dtype(),
            }.get
        else:
            types_mapper = None
    elif dtype_backend is lib.no_default:
        if dtype is not None:
            # GH#56136 Avoid lossy conversion to float64
            # We'll convert to numpy below if
            types_mapper = {
                pa.int8(): Int8Dtype(),
                pa.int16(): Int16Dtype(),
                pa.int32(): Int32Dtype(),
                pa.int64(): Int64Dtype(),
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
    df: DataFrame,
    dtype_backend: DtypeBackend | lib.NoDefault,
    dtype: DtypeArg | None,
    names: Sequence[Hashable] | None,
) -> DataFrame:
    if dtype is not None and dtype_backend is lib.no_default:
        # GH#56136 apply any user-provided dtype, and convert any IntegerDtype
        #  columns the user didn't explicitly ask for.
        if isinstance(dtype, dict):
            if names is not None:
                df.columns = names

            cmp_dtypes = {
                Int8Dtype(),
                Int16Dtype(),
                Int32Dtype(),
                Int64Dtype(),
            }
            for col in df.columns:
                if col not in dtype and df[col].dtype in cmp_dtypes:
                    # Any key that the user didn't explicitly specify
                    #  that got converted to IntegerDtype now gets converted
                    #  to numpy dtype: float64 if it has missing values, as
                    #  pyarrow gives without the IntegerDtype mapper.
                    if df[col].hasnans:
                        dtype[col] = np.dtype("float64")
                    else:
                        dtype[col] = df[col].dtype.numpy_dtype

            # Ignore non-existent columns from dtype mapping
            # like other parsers do
            dtype = {
                key: pandas_dtype(dtype[key]) for key in dtype if key in df.columns
            }

            # GH#65056 pyarrow returns datetime.date objects which trigger
            #  a deprecation warning when compared against DatetimeIndex
            #  categories. Convert to datetime64 first since the user can't
            #  control pyarrow's behavior.
            for col, col_dtype in dtype.items():
                if (
                    isinstance(col_dtype, CategoricalDtype)
                    and col_dtype.categories is not None
                    and isinstance(col_dtype.categories, DatetimeIndex)
                    and df[col].dtype == np.dtype("object")
                    and lib.infer_dtype(df[col]) == "date"
                ):
                    df[col] = to_datetime(df[col])

        else:
            dtype = pandas_dtype(dtype)

        try:
            casted = df.astype(dtype)
        except TypeError as err:
            # GH#44901 reraise to keep api consistent
            raise ValueError(str(err)) from err
        validate_integer_casts(df, casted)
        df = casted

        # GH#56136 IntegerDtype was used above to avoid lossy float64
        #  conversion in pyarrow; convert back to numpy now that the data
        #  is categorical
        col_dtypes = df.dtypes
        for i in range(len(df.columns)):
            col_dtype = col_dtypes.iloc[i]
            if not isinstance(col_dtype, CategoricalDtype):
                continue
            cat_arr_dtype = col_dtype.categories.dtype
            if not isinstance(cat_arr_dtype, IntegerDtype):
                continue
            if isinstance(dtype, dict):
                requested = dtype.get(df.columns[i])
            else:
                requested = dtype
            if (
                isinstance(requested, CategoricalDtype)
                and requested.categories is not None
            ):
                # the user explicitly asked for these categories
                continue
            new_cat_dtype = CategoricalDtype(
                categories=col_dtype.categories.astype(cat_arr_dtype.numpy_dtype),
                ordered=col_dtype.ordered,
            )
            df.isetitem(i, df.iloc[:, i].astype(new_cat_dtype))

    if not using_string_dtype() and dtype != "str" and dtype_backend is lib.no_default:
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


def _numpy_dtype(dtype: DtypeObj) -> np.dtype | None:
    # the numpy dtype underlying a column's dtype, if any
    if isinstance(dtype, np.dtype):
        return dtype
    if isinstance(dtype, SparseDtype):
        return dtype.subtype
    if isinstance(dtype, (BaseMaskedDtype, ArrowDtype)):
        return dtype.numpy_dtype
    return None


def _can_wrap(old: DtypeObj, new: DtypeObj) -> bool:
    """
    Whether casting from ``old`` to ``new`` can wrap around an integer value.
    """
    old_np = _numpy_dtype(old)
    new_np = _numpy_dtype(new)
    if old_np is None or new_np is None:
        return False
    if new_np.kind not in "iu" or old_np.kind not in "iuf":
        return False
    # NA cannot wrap, so e.g. Int64 -> int64 is as safe as int64 -> int64
    return not np.can_cast(old_np, new_np)


def validate_integer_casts(original: DataFrame, casted: DataFrame) -> None:
    """
    Apply validate_integer_cast to each column of a DataFrame cast.
    """
    for i, (old, new) in enumerate(zip(original.dtypes, casted.dtypes, strict=True)):
        # check dtypes first; iloc on every column is slow for wide frames
        if _can_wrap(old, new):
            validate_integer_cast(original.iloc[:, i], new)


def validate_integer_cast(original: Series, dtype: DtypeObj) -> None:
    """
    Raise if casting a column to an integer ``dtype`` wrapped around (GH#38013).

    Discarding a float's fractional part is allowed, as in the python engine.
    """
    if not _can_wrap(original.dtype, dtype):
        return
    lo, hi = original.min(), original.max()
    if isna(lo):
        return
    info = np.iinfo(_numpy_dtype(dtype))
    # int() truncates exactly, so e.g. 2.0**63 does not equal int64 max
    if int(lo) < info.min or int(hi) > info.max:
        raise ValueError(
            f"cannot safely convert passed user dtype of {dtype} for "
            f"{_numpy_dtype(original.dtype)} dtyped data in column {original.name}"
        )


def _maybe_convert_string_to_object(
    data: Series | Index,
) -> Series | Index | None:
    if isinstance(data.dtype, StringDtype) and data.dtype.na_value is np.nan:
        return data.astype("object").fillna(None)
    elif isinstance(data.dtype, CategoricalDtype):
        cat_dtype = data.dtype.categories.dtype
        if isinstance(cat_dtype, StringDtype) and cat_dtype.na_value is np.nan:
            # not astype: for ordered categoricals CategoricalDtype.__eq__
            #  ignores the categories' dtype, so astype would no-op
            cat_arr = cast("Categorical", data._values)
            new_arr = cat_arr.set_categories(data.dtype.categories.astype("object"))
            if isinstance(data, Index):
                return CategoricalIndex(new_arr, name=data.name)
            from pandas.core.series import Series

            return Series(new_arr, index=data.index, name=data.name, copy=False)

    # no conversion needed
    return None


def _maybe_convert_string_index_to_object(index: Index) -> Index | None:
    if isinstance(index, MultiIndex):
        if any(
            isinstance(level.dtype, StringDtype) and level.dtype.na_value is np.nan
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
        return cast("Index | None", _maybe_convert_string_to_object(index))


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
        return dt.UTC

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


def _normalize_timezone_index(index: Index) -> Index:
    if isinstance(index, MultiIndex):
        if any(isinstance(level.dtype, DatetimeTZDtype) for level in index.levels):
            levels = [_normalize_timezone_index(level) for level in index.levels]
            return index.set_levels(levels)

        return index

    if isinstance(index.dtype, DatetimeTZDtype):
        normalized_tz = _normalize_pytz_timezone(index.dtype.tz)
        if normalized_tz is not index.dtype.tz:
            return index.tz_convert(normalized_tz)  # type: ignore[attr-defined]

    return index


def _normalize_timezone_dtypes(df: DataFrame) -> DataFrame:
    """
    PyArrow below version 25 uses pytz by default for timezones, but pandas uses
    zoneinfo / datetime.timezone since pandas 3.0.

    Can be dropped once the minimum supported pyarrow is 25, which uses zoneinfo
    itself (https://github.com/apache/arrow/pull/49694).
    """
    if pytz is not None and pa_version_under25p0:
        # Convert any pytz timezones to zoneinfo / fixed offset timezones
        if any(
            isinstance(dtype, DatetimeTZDtype) for dtype in df._mgr.get_unique_dtypes()
        ):
            col_indices = df._select_dtypes_indices(DatetimeTZDtype)
            for i in col_indices:
                col = df.iloc[:, i]
                normalized_tz = _normalize_pytz_timezone(col.dtype.tz)
                if normalized_tz is not col.dtype.tz:
                    df.isetitem(i, col.dt.tz_convert(normalized_tz))

        df.index = _normalize_timezone_index(df.index)
        df.columns = _normalize_timezone_index(df.columns)

    return df
