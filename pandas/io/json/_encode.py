"""
Convert pandas objects to JSON-native Python objects for ``to_json``.

This module owns the pandas-specific serialization semantics: orients,
date formatting, missing values, label stringification and nested objects.
The result holds only dict, list, str, int, float, bool and None, so that
any JSON encoder can write it. When the engine supports it (see
``EncodeOptions``), float columns may additionally be handed over as numpy
arrays holding NaN.
"""

from __future__ import annotations

from dataclasses import (
    dataclass,
    field,
)
import datetime as dt
from decimal import Decimal
import math
from typing import (
    TYPE_CHECKING,
    Any,
)

import numpy as np

from pandas._libs import (
    lib,
    writers,
)
from pandas._libs.missing import NA
from pandas._libs.tslibs import (
    NaT,
    Timedelta,
    Timestamp,
    iNaT,
)
from pandas._libs.tslibs.np_datetime import astype_overflowsafe

from pandas.core.dtypes.dtypes import DatetimeTZDtype

from pandas.core.arrays import TimedeltaArray
from pandas.core.frame import DataFrame
from pandas.core.indexes.base import Index
from pandas.core.indexes.multi import MultiIndex
from pandas.core.series import Series

from pandas.io.json import _engine

if TYPE_CHECKING:
    from collections.abc import Callable

    from pandas._typing import ArrayLike


@dataclass(frozen=True)
class EncodeOptions:
    iso_dates: bool
    date_unit: str
    default_handler: Callable[[Any], Any] | None = None
    index: bool = True
    # replace NaN / infinity by None; not needed for an engine that writes
    # them as null itself
    nan_to_none: bool = field(default_factory=lambda: not _engine.native_nan())
    # hand numeric columns over as numpy arrays instead of lists
    keep_numpy: bool = field(default_factory=_engine.native_numpy)
    # splice label -> value objects together as pre-encoded text instead of
    # building dicts; only for compact output, as fragments are not indented
    fragments: bool = field(default_factory=_engine.native_fragments)


def encode(obj: Any, orient: str, options: EncodeOptions) -> Any:
    """
    Convert ``obj`` to JSON-native Python objects.

    Parameters
    ----------
    obj : DataFrame, Series, Index, or any JSON-serializable object
    orient : {"split", "records", "index", "columns", "values"}
    options : EncodeOptions
    """
    if isinstance(obj, DataFrame):
        return _encode_frame(obj, orient, options)
    if isinstance(obj, Series):
        return _encode_series(obj, orient, options)
    if isinstance(obj, Index):
        return _encode_index(obj, orient, options)
    return _encode_scalar(obj, orient, options)


def _nested_orient(orient: str) -> str:
    # objects nested inside a split-oriented object are written as plain
    # values; any other orient carries over to nested objects
    return "values" if orient == "split" else orient


def _array_values(arr: ArrayLike) -> tuple[np.ndarray, bool]:
    """
    The numpy array to serialize for a column, and whether its datetime64
    values are UTC (tz-aware data, written with a trailing "Z").
    """
    if isinstance(arr, np.ndarray):
        return arr, False
    if isinstance(arr.dtype, DatetimeTZDtype):
        return arr._ndarray, True
    return arr._values_for_json(), False


# --------------------------------------------------------------------------
# pandas containers


def _encode_frame(df: DataFrame, orient: str, options: EncodeOptions) -> Any:
    nested = _nested_orient(orient)
    columns = []
    for arr in df._iter_column_arrays():
        values, is_utc = _array_values(arr)
        columns.append(_encode_values(values, is_utc, nested, options))
    nrows = len(df)

    if orient == "split":
        result: dict[str, Any] = {
            "columns": _encode_index(df.columns, "values", options)
        }
        if options.index:
            result["index"] = _encode_index(df.index, "values", options)
        result["data"] = _rows(columns, nrows, options)
        return result
    if orient == "values":
        return _rows(columns, nrows, options)
    if orient == "records":
        if not columns:
            return []
        if options.fragments:
            keys = _encode_labels_bytes(df.columns, options)
            return _engine.fragment(_transpose(columns, nrows, options, keys))
        keys = _encode_labels(df.columns, options)
        return writers.build_json_rows(_as_lists(columns), nrows, keys)
    if orient == "index":
        if not columns:
            return {}
        if options.fragments:
            keys = _encode_labels_bytes(df.columns, options)
            row_labels = _encode_labels_bytes(df.index, options)
            return _engine.fragment(
                _transpose(columns, nrows, options, keys, row_labels)
            )
        keys = _encode_labels(df.columns, options)
        records = writers.build_json_rows(_as_lists(columns), nrows, keys)
        return dict(zip(_encode_labels(df.index, options), records, strict=True))
    if orient == "columns":
        if options.fragments:
            row_labels = _encode_labels_bytes(df.index, options)
            encoded = [
                writers.json_zip_object(row_labels, _engine.dumps_bytes(col))
                for col in columns
            ]
            keys = _encode_labels_bytes(df.columns, options)
            return _engine.fragment(writers.json_object_of_fragments(keys, encoded))
        row_labels = _encode_labels(df.index, options)
        return {
            label: dict(zip(row_labels, _as_list(col), strict=True))
            for label, col in zip(
                _encode_labels(df.columns, options), columns, strict=True
            )
        }
    raise ValueError(f"Invalid value '{orient}' for option 'orient'")


def _encode_series(ser: Series, orient: str, options: EncodeOptions) -> Any:
    values, is_utc = _array_values(ser.array)
    data = _encode_values(values, is_utc, _nested_orient(orient), options)

    if orient == "split":
        result: dict[str, Any] = {"name": _encode_scalar(ser.name, "values", options)}
        if options.index:
            result["index"] = _encode_index(ser.index, "values", options)
        result["data"] = _maybe_list(data, options)
        return result
    if orient in ("index", "columns"):
        if options.fragments:
            return _engine.fragment(
                writers.json_zip_object(
                    _encode_labels_bytes(ser.index, options),
                    _engine.dumps_bytes(data),
                )
            )
        return dict(
            zip(_encode_labels(ser.index, options), _as_list(data), strict=True)
        )
    if orient in ("records", "values"):
        return _maybe_list(data, options)
    raise ValueError(f"Invalid value '{orient}' for option 'orient'")


def _encode_index(index: Index, orient: str, options: EncodeOptions) -> Any:
    if orient == "split":
        return {
            "name": _encode_scalar(index.name, "values", options),
            "data": _encode_index(index, "values", options),
        }
    if isinstance(index, MultiIndex):
        values, is_utc = index.values, False
    else:
        values, is_utc = _array_values(index.array)
    return _maybe_list(_encode_values(values, is_utc, "values", options), options)


def _homogeneous_2d(columns: list, options: EncodeOptions) -> np.ndarray | None:
    """Columns as one 2-D numpy array when they share a numpy dtype."""
    if (
        options.keep_numpy
        and all(isinstance(col, np.ndarray) for col in columns)
        and len({col.dtype for col in columns}) == 1
    ):
        return np.ascontiguousarray(np.column_stack(columns))
    return None


def _rows(columns: list, nrows: int, options: EncodeOptions) -> Any:
    """Row-major data for orient="values" / "split"."""
    if not columns:
        return []
    values = _homogeneous_2d(columns, options)
    if values is not None:
        return values
    if options.fragments:
        return _engine.fragment(_transpose(columns, nrows, options))
    return writers.build_json_rows(_as_lists(columns), nrows)


def _transpose(
    columns: list,
    nrows: int,
    options: EncodeOptions,
    keys: bytes | None = None,
    labels: bytes | None = None,
) -> bytes:
    """Row-oriented JSON text built from the JSON text of each column."""
    values = _homogeneous_2d(columns, options)
    if values is not None and keys is not None:
        rows = _engine.dumps_bytes(values)
        if labels is None:
            return writers.json_zip_records(keys, rows)
        return writers.json_zip_records(keys, rows, labels)
    encoded = [_engine.dumps_bytes(col) for col in columns]
    return writers.json_transpose(encoded, keys, labels)


def _as_list(values: np.ndarray | list) -> list:
    if isinstance(values, np.ndarray):
        return values.tolist()
    return values


def _as_lists(columns: list) -> list[list]:
    return [_as_list(col) for col in columns]


def _maybe_list(values: np.ndarray | list, options: EncodeOptions) -> Any:
    if options.keep_numpy:
        return values
    return _as_list(values)


# --------------------------------------------------------------------------
# labels


def _encode_labels(index: Index, options: EncodeOptions) -> list[str]:
    """JSON object keys for the labels of ``index``."""
    if isinstance(index, MultiIndex):
        return [str(label) for label in index]
    values, is_utc = _array_values(index.array)
    kind = values.dtype.kind
    if kind in "Mm":
        encoded = _encode_datetimelike(values, is_utc, options)
        return ["null" if label is None else str(label) for label in encoded]
    if kind == "O":
        return [_encode_object_label(label, options) for label in values]
    return lib.ensure_string_array(values, skipna=False).tolist()


def _encode_labels_bytes(index: Index, options: EncodeOptions) -> bytes:
    """The labels of ``index`` as a JSON array of strings."""
    if not isinstance(index, MultiIndex):
        values, is_utc = _array_values(index.array)
        if values.dtype.kind in "iu":
            # the JSON text of an integer is its str()
            return writers.json_quote_elements(_engine.dumps_bytes(values))
        if values.dtype.kind in "Mm":
            if options.iso_dates:
                return _encode_datetimelike_bytes(
                    values, is_utc, options, nat_as_string=True
                )
            # quoting the null element gives the "null" label
            return writers.json_quote_elements(
                _encode_datetimelike_bytes(values, is_utc, options)
            )
    return _engine.dumps_bytes(_encode_labels(index, options))


def _encode_object_label(label: Any, options: EncodeOptions) -> str:
    if label is NaT or isinstance(
        label, (dt.date, dt.timedelta, np.datetime64, np.timedelta64)
    ):
        encoded = _encode_scalar(label, "values", options)
        return "null" if encoded is None else str(encoded)
    return str(label)


# --------------------------------------------------------------------------
# values


def _encode_values(
    values: np.ndarray, is_utc: bool, orient: str, options: EncodeOptions
) -> np.ndarray | list:
    """
    Encode a 1-D array. Returns a list, or a numpy array of int, uint, bool
    or float64 (the latter possibly holding NaN when the engine writes NaN
    as null itself).
    """
    kind = values.dtype.kind
    if kind in "Mm":
        if options.fragments:
            if not options.iso_dates:
                cast = _cast_datetimelike(values, options).view("i8")
                if not (cast == iNaT).any():
                    return cast
            return _engine.fragment(_encode_datetimelike_bytes(values, is_utc, options))
        return _encode_datetimelike(values, is_utc, options)
    if kind == "f":
        if values.dtype != np.float64:
            values = values.astype(np.float64)
        if options.nan_to_none:
            mask = ~np.isfinite(values)
            if mask.any():
                result = values.astype(object)
                result[mask] = None
                return result.tolist()
        return values
    if kind in "iub":
        return values
    if kind == "O":
        return writers.normalize_json_objects(
            values,
            lambda value: _encode_scalar(value, orient, options),
            options.nan_to_none,
            options.iso_dates,
            options.date_unit,
        )
    # anything else (bytes, str, complex, void) goes through the scalar path
    return _encode_values(values.astype(object), is_utc, orient, options)


def _cast_datetimelike(values: np.ndarray, options: EncodeOptions) -> np.ndarray:
    """
    Cast to ``options.date_unit`` where the output depends on it: always for
    epoch integers, and for ISO strings only when the unit is coarser (a finer
    unit only adds trailing zeros, so the values keep their own unit).
    """
    kind = values.dtype.kind
    if options.iso_dates and kind == "m":
        # durations are written at their own resolution
        return values
    target = np.dtype(f"{kind}8[{options.date_unit}]")
    if values.dtype != target and (
        not options.iso_dates or _is_coarser(target, values.dtype)
    ):
        values = astype_overflowsafe(values, target, copy=False)
    return values


def _encode_datetimelike_bytes(
    values: np.ndarray, is_utc: bool, options: EncodeOptions, nat_as_string=False
) -> bytes:
    """
    JSON array text for datetime64 / timedelta64 values: ISO 8601 strings or
    epoch integers in ``options.date_unit``; NaT becomes null.
    """
    values = _cast_datetimelike(values, options)
    unit = np.datetime_data(values.dtype)[0]
    if unit not in _UNIT_ORDER:
        # e.g. "D" from a numpy scalar; bring it into the supported range
        values = astype_overflowsafe(
            values, np.dtype(f"{values.dtype.kind}8[s]"), copy=False
        )
        unit = "s"
    i8 = values.view("i8")
    if not options.iso_dates:
        return writers.int64_to_json(i8, nat_as_string)
    if values.dtype.kind == "m":
        return writers.timedelta64_to_json(i8, unit, nat_as_string)
    return writers.datetime64_to_json(
        i8, unit, options.date_unit, is_utc, nat_as_string
    )


def _iso_duration(td: Timedelta) -> str:
    """
    ISO 8601 duration with the fractional seconds in groups of three digits
    down to the smallest non-zero unit (the format of
    ``writers.timedelta64_to_json``).
    """
    c = td.components
    if c.nanoseconds:
        frac = f".{c.milliseconds:03d}{c.microseconds:03d}{c.nanoseconds:03d}"
    elif c.microseconds:
        frac = f".{c.milliseconds:03d}{c.microseconds:03d}"
    elif c.milliseconds:
        frac = f".{c.milliseconds:03d}"
    else:
        frac = ""
    return f"P{c.days}DT{c.hours}H{c.minutes}M{c.seconds}{frac}S"


def _encode_datetimelike(
    values: np.ndarray, is_utc: bool, options: EncodeOptions
) -> list:
    """
    Encode datetime64 / timedelta64 values as ISO 8601 strings or epoch
    integers in ``options.date_unit``; NaT becomes None.
    """
    kind = values.dtype.kind
    if np.datetime_data(values.dtype)[0] not in _UNIT_ORDER:
        # e.g. "D" from a numpy scalar; bring it into the supported range
        values = astype_overflowsafe(values, np.dtype(f"{kind}8[s]"), copy=False)
    if options.iso_dates and kind == "m":
        tda = TimedeltaArray._simple_new(values, dtype=values.dtype)
        return [None if td is NaT else _iso_duration(td) for td in tda]

    values = _cast_datetimelike(values, options)
    mask = values.view("i8") == iNaT
    if options.iso_dates:
        timezone = "UTC" if is_utc else "naive"
        result = np.datetime_as_string(
            values, unit=options.date_unit, timezone=timezone
        ).astype(object)
    else:
        result = values.view("i8").astype(object)
    if mask.any():
        result[mask] = None
    return result.tolist()


_UNIT_ORDER = ("s", "ms", "us", "ns")


def _is_coarser(target: np.dtype, source: np.dtype) -> bool:
    target_unit = np.datetime_data(target)[0]
    source_unit = np.datetime_data(source)[0]
    if source_unit not in _UNIT_ORDER:
        # e.g. "D" or "W" from a numpy scalar; always cast those
        return True
    return _UNIT_ORDER.index(target_unit) < _UNIT_ORDER.index(source_unit)


def _encode_key(key: Any) -> str:
    if isinstance(key, str):
        return key
    if isinstance(key, bytes):
        return key.decode("utf-8")
    return str(key)


def _encode_scalar(value: Any, orient: str, options: EncodeOptions) -> Any:
    """
    Encode a single Python object, recursing into containers.
    """
    if value is None or value is NaT or value is NA:
        return None
    if isinstance(value, (bool, str, int)):
        return value
    if isinstance(value, float):
        if options.nan_to_none and not math.isfinite(value):
            return None
        return value
    if isinstance(value, dict):
        return {
            _encode_key(key): _encode_scalar(item, orient, options)
            for key, item in value.items()
        }
    if isinstance(value, (list, tuple, set, frozenset)):
        return [_encode_scalar(item, orient, options) for item in value]
    if isinstance(value, (DataFrame, Series, Index)):
        return encode(value, orient, options)
    if isinstance(value, np.ndarray):
        if value.ndim == 0:
            return _encode_scalar(value[()], orient, options)
        if value.ndim == 1:
            return _maybe_list(_encode_values(value, False, orient, options), options)
        return [_encode_scalar(row, orient, options) for row in value]
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.generic):
        if isinstance(value, (np.datetime64, np.timedelta64)):
            if np.isnat(value):
                return None
            return _encode_datetimelike(np.array([value]), False, options)[0]
        return _encode_scalar(value.item(), orient, options)
    if isinstance(value, dt.date):
        # datetime.datetime is a subclass of datetime.date; tz-aware values
        # are written in UTC
        ts = Timestamp(value)
        if ts is NaT:
            return None
        arr = np.array([ts._value], dtype=f"M8[{ts.unit}]")
        return _encode_datetimelike(arr, ts.tz is not None, options)[0]
    if isinstance(value, dt.time):
        return value.isoformat()
    if isinstance(value, dt.timedelta):
        td = Timedelta(value)
        if td is NaT:
            return None
        arr = np.array([td._value], dtype=f"m8[{td.unit}]")
        return _encode_datetimelike(arr, False, options)[0]
    if isinstance(value, Decimal):
        if value.is_nan():
            return None
        return format(value, "f")
    if options.default_handler is not None:
        return _encode_scalar(options.default_handler(value), orient, options)
    return _encode_attributes(value, orient, options)


def _encode_attributes(value: Any, orient: str, options: EncodeOptions) -> dict:
    """
    Encode an object of unknown type as the mapping of its public,
    non-callable attributes.
    """
    result = {}
    for name in dir(value):
        if name.startswith("_"):
            continue
        try:
            attr = getattr(value, name)
        except Exception:
            continue
        if callable(attr):
            continue
        result[name] = _encode_scalar(attr, orient, options)
    return result
