from __future__ import annotations

from datetime import (
    date,
    datetime,
    time,
    timedelta,
)
from decimal import Decimal

import pytest

from pandas.core.dtypes.dtypes import (
    ArrowDtype,
)

import pandas as pd
import pandas._testing as tm

pa = pytest.importorskip("pyarrow")


@pytest.fixture(params=tm.ALL_PYARROW_DTYPES, ids=str)
def dtype(request):
    return ArrowDtype(pyarrow_dtype=request.param)


@pytest.fixture
def data(dtype):
    pa_dtype = dtype.pyarrow_dtype
    if pa.types.is_boolean(pa_dtype):
        data = [True, False] * 2 + [None] + [True, False] + [None] + [True, False]
    elif pa.types.is_floating(pa_dtype):
        data = [1.0, 0.0] * 2 + [None] + [-2.0, -1.0] + [None] + [0.5, 99.5]
    elif pa.types.is_signed_integer(pa_dtype):
        data = [1, 0] * 2 + [None] + [-2, -1] + [None] + [1, 99]
    elif pa.types.is_unsigned_integer(pa_dtype):
        data = [1, 0] * 2 + [None] + [2, 1] + [None] + [1, 99]
    elif pa.types.is_decimal(pa_dtype):
        data = (
            [Decimal("1"), Decimal("0.0")] * 2
            + [None]
            + [Decimal("-2.0"), Decimal("-1.0")]
            + [None]
            + [Decimal("0.5"), Decimal("33.123")]
        )
    elif pa.types.is_date(pa_dtype):
        data = (
            [date(2022, 1, 1), date(1999, 12, 31)] * 2
            + [None]
            + [date(2022, 1, 1), date(2022, 1, 1)]
            + [None]
            + [date(1999, 12, 31), date(1999, 12, 31)]
        )
    elif pa.types.is_timestamp(pa_dtype):
        data = (
            [datetime(2020, 1, 1, 1, 1, 1, 1), datetime(1999, 1, 1, 1, 1, 1, 1)] * 2
            + [None]
            + [datetime(2020, 1, 1, 1), datetime(1999, 1, 1, 1)]
            + [None]
            + [datetime(2020, 1, 1), datetime(1999, 1, 1)]
        )
    elif pa.types.is_duration(pa_dtype):
        data = (
            [timedelta(1), timedelta(1, 1)] * 2
            + [None]
            + [timedelta(-1), timedelta(0)]
            + [None]
            + [timedelta(-10), timedelta(10)]
        )
    elif pa.types.is_time(pa_dtype):
        data = (
            [time(12, 0), time(0, 12)] * 2
            + [None]
            + [time(0, 0), time(1, 1)]
            + [None]
            + [time(0, 5), time(5, 0)]
        )
    elif pa.types.is_string(pa_dtype):
        data = ["a", "b"] * 2 + [None] + ["1", "2"] + [None] + ["!", ">"]
    elif pa.types.is_binary(pa_dtype):
        data = [b"a", b"b"] * 2 + [None] + [b"1", b"2"] + [None] + [b"!", b">"]
    else:
        raise NotImplementedError
    return pd.array(data, dtype=dtype)


@pytest.fixture
def data_missing(data):
    """Length-2 array with [NA, Valid]"""
    return type(data)._from_sequence([None, data[0]], dtype=data.dtype)


@pytest.fixture
def data_for_grouping(dtype):
    """
    Data for factorization, grouping, and unique tests.

    Expected to be like [B, B, NA, NA, A, A, B, C]

    Where A < B < C and NA is missing
    """
    pa_dtype = dtype.pyarrow_dtype
    if pa.types.is_boolean(pa_dtype):
        A = False
        B = True
        C = True
    elif pa.types.is_floating(pa_dtype):
        A = -1.1
        B = 0.0
        C = 1.1
    elif pa.types.is_signed_integer(pa_dtype):
        A = -1
        B = 0
        C = 1
    elif pa.types.is_unsigned_integer(pa_dtype):
        A = 0
        B = 1
        C = 10
    elif pa.types.is_date(pa_dtype):
        A = date(1999, 12, 31)
        B = date(2010, 1, 1)
        C = date(2022, 1, 1)
    elif pa.types.is_timestamp(pa_dtype):
        A = datetime(1999, 1, 1, 1, 1, 1, 1)
        B = datetime(2020, 1, 1)
        C = datetime(2020, 1, 1, 1)
    elif pa.types.is_duration(pa_dtype):
        A = timedelta(-1)
        B = timedelta(0)
        C = timedelta(1, 4)
    elif pa.types.is_time(pa_dtype):
        A = time(0, 0)
        B = time(0, 12)
        C = time(12, 12)
    elif pa.types.is_string(pa_dtype):
        A = "a"
        B = "b"
        C = "c"
    elif pa.types.is_binary(pa_dtype):
        A = b"a"
        B = b"b"
        C = b"c"
    elif pa.types.is_decimal(pa_dtype):
        A = Decimal("-1.1")
        B = Decimal("0.0")
        C = Decimal("1.1")
    else:
        raise NotImplementedError
    return pd.array([B, B, None, None, A, A, B, C], dtype=dtype)


@pytest.fixture
def data_for_sorting(data_for_grouping):
    """
    Length-3 array with a known sort order.

    This should be three items [B, C, A] with
    A < B < C
    """
    return type(data_for_grouping)._from_sequence(
        [data_for_grouping[0], data_for_grouping[7], data_for_grouping[4]],
        dtype=data_for_grouping.dtype,
    )


@pytest.fixture(params=[True, False])
def as_series(request):
    """
    Boolean fixture to support arr and Series(arr) comparison testing.
    """
    return request.param
