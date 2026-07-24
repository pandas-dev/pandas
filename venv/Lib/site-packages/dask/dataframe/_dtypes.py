from __future__ import annotations

from datetime import date, time
from decimal import Decimal

import pandas as pd

from dask.dataframe.extensions import make_array_nonempty, make_scalar


@make_array_nonempty.register(pd.DatetimeTZDtype)
def _(dtype):
    return pd.array([pd.Timestamp(1), pd.NaT], dtype=dtype)


@make_scalar.register(pd.DatetimeTZDtype)
def _(x):
    return pd.Timestamp(1, tz=x.tz, unit=x.unit)


@make_array_nonempty.register(pd.StringDtype)
def _(dtype):
    return pd.array(["a", pd.NA], dtype=dtype)


@make_array_nonempty.register(pd.ArrowDtype)
def _make_array_nonempty_pyarrow_dtype(dtype):
    import pyarrow as pa

    if pa.types.is_integer(dtype.pyarrow_dtype):
        data = [1, 2]
    elif pa.types.is_floating(dtype.pyarrow_dtype):
        data = [1.5, 2.5]
    elif pa.types.is_boolean(dtype.pyarrow_dtype):
        data = [True, False]
    elif pa.types.is_string(dtype.pyarrow_dtype) or pa.types.is_large_string(
        dtype.pyarrow_dtype
    ):
        data = ["a", "b"]
    elif pa.types.is_timestamp(dtype.pyarrow_dtype):
        data = [pd.Timestamp("1970-01-01"), pd.Timestamp("1970-01-02")]
    elif pa.types.is_date(dtype.pyarrow_dtype):
        data = [date(1970, 1, 1), date(1970, 1, 2)]
    elif pa.types.is_binary(dtype.pyarrow_dtype) or pa.types.is_large_binary(
        dtype.pyarrow_dtype
    ):
        data = [b"a", b"b"]
    elif pa.types.is_decimal(dtype.pyarrow_dtype):
        data = [Decimal("1"), Decimal("0.0")]
    elif pa.types.is_duration(dtype.pyarrow_dtype):
        data = [pd.Timedelta("1 day"), pd.Timedelta("2 days")]
    elif pa.types.is_time(dtype.pyarrow_dtype):
        data = [time(12, 0), time(0, 12)]
    else:
        data = dtype.empty(2)
    return pd.array(data, dtype=dtype)


@make_scalar.register(str)
def _(x):
    return "s"


@make_array_nonempty.register(pd.BooleanDtype)
def _(dtype):
    return pd.array([True, pd.NA], dtype=dtype)


@make_scalar.register(bool)
def _(x):
    return True
