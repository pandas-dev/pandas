from __future__ import annotations

import pickle
from datetime import date, datetime, time, timedelta
from decimal import Decimal

import numpy as np
import pandas as pd
import pandas._testing as tm
import pytest

from dask.dataframe.utils import get_string_dtype

pa = pytest.importorskip("pyarrow")
import dask.dataframe as dd

# Tests are from https://github.com/pandas-dev/pandas/pull/49078


@pytest.fixture
def data(dtype):
    pa_dtype = dtype.pyarrow_dtype
    if pa.types.is_boolean(pa_dtype):
        data = [True, False] * 4 + [None] + [True, False] * 44 + [None] + [True, False]
    elif pa.types.is_floating(pa_dtype):
        data = [1.0, 0.0] * 4 + [None] + [-2.0, -1.0] * 44 + [None] + [0.5, 99.5]
    elif pa.types.is_signed_integer(pa_dtype):
        data = [1, 0] * 4 + [None] + [-2, -1] * 44 + [None] + [1, 99]
    elif pa.types.is_unsigned_integer(pa_dtype):
        data = [1, 0] * 4 + [None] + [2, 1] * 44 + [None] + [1, 99]
    elif pa.types.is_decimal(pa_dtype):
        data = (
            [Decimal("1"), Decimal("0.0")] * 4
            + [None]
            + [Decimal("-2.0"), Decimal("-1.0")] * 44
            + [None]
            + [Decimal("0.5"), Decimal("33.123")]
        )
    elif pa.types.is_date(pa_dtype):
        data = (
            [date(2022, 1, 1), date(1999, 12, 31)] * 4
            + [None]
            + [date(2022, 1, 1), date(2022, 1, 1)] * 44
            + [None]
            + [date(1999, 12, 31), date(1999, 12, 31)]
        )
    elif pa.types.is_timestamp(pa_dtype):
        data = (
            [datetime(2020, 1, 1, 1, 1, 1, 1), datetime(1999, 1, 1, 1, 1, 1, 1)] * 4
            + [None]
            + [datetime(2020, 1, 1, 1), datetime(1999, 1, 1, 1)] * 44
            + [None]
            + [datetime(2020, 1, 1), datetime(1999, 1, 1)]
        )
    elif pa.types.is_duration(pa_dtype):
        data = (
            [timedelta(1), timedelta(1, 1)] * 4
            + [None]
            + [timedelta(-1), timedelta(0)] * 44
            + [None]
            + [timedelta(-10), timedelta(10)]
        )
    elif pa.types.is_time(pa_dtype):
        data = (
            [time(12, 0), time(0, 12)] * 4
            + [None]
            + [time(0, 0), time(1, 1)] * 44
            + [None]
            + [time(0, 5), time(5, 0)]
        )
    elif pa.types.is_string(pa_dtype):
        data = ["a", "b"] * 4 + [None] + ["1", "2"] * 44 + [None] + ["!", ">"]
    elif pa.types.is_binary(pa_dtype):
        data = [b"a", b"b"] * 4 + [None] + [b"1", b"2"] * 44 + [None] + [b"!", b">"]
    else:
        raise NotImplementedError
    return pd.array(data * 100, dtype=dtype)


PYARROW_TYPES = tm.ALL_PYARROW_DTYPES


@pytest.fixture(params=PYARROW_TYPES, ids=str)
def dtype(request):
    return pd.ArrowDtype(pyarrow_dtype=request.param)


def test_pickle_roundtrip(data):
    expected = pd.Series(data)
    expected_sliced = expected.head(2)
    full_pickled = pickle.dumps(expected)
    sliced_pickled = pickle.dumps(expected_sliced)

    # Make sure slicing gives a large reduction in serialized bytes
    assert len(full_pickled) > len(sliced_pickled) * 3

    result = pickle.loads(full_pickled)
    tm.assert_series_equal(result, expected)

    result_sliced = pickle.loads(sliced_pickled)
    tm.assert_series_equal(result_sliced, expected_sliced)


@pytest.mark.parametrize(
    "string_dtype",
    [
        "stringdtype",
        "arrowdtype",
    ],
)
def test_pickle_roundtrip_pyarrow_string_implementations(string_dtype):
    # There are two pyarrow string implementations in pandas.
    # This tests that both implementations have similar serialization performance.
    if string_dtype == "stringdtype":
        string_dtype = pd.StringDtype("pyarrow")
    else:
        string_dtype = pd.ArrowDtype(pa.string())
    expected = pd.Series(map(str, range(1_000)), dtype=string_dtype)
    expected_sliced = expected.head(2)
    full_pickled = pickle.dumps(expected)
    sliced_pickled = pickle.dumps(expected_sliced)

    # Make sure slicing gives a large reduction in serialized bytes
    assert len(full_pickled) > len(sliced_pickled) * 3

    result = pickle.loads(full_pickled)
    tm.assert_series_equal(result, expected)

    result_sliced = pickle.loads(sliced_pickled)
    tm.assert_series_equal(result_sliced, expected_sliced)


def test_inplace_modification_read_only():
    arr = np.array([(1, 2), None, 1], dtype="object")
    base = pd.Series(arr, copy=False, dtype=object, name="a")
    base_copy = pickle.loads(pickle.dumps(base))
    base_copy.values.flags.writeable = False
    dtype = get_string_dtype()
    tm.assert_series_equal(
        dd.from_array(base_copy.values, columns="a").compute(),
        base.astype(dtype),
    )
