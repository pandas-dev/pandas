from __future__ import annotations

import numpy as np
import pytest

from pandas.core.dtypes.cast import find_common_type
from pandas.core.dtypes.dtypes import (
    ArrowDtype,
)

import pandas as pd
import pandas._testing as tm

pa = pytest.importorskip("pyarrow")


def test_concat_timestamp_tz_aware_and_naive():
    # GH#69029 the common dtype went through numpy_dtype, which drops the tz, so
    #  the tz-aware values came back shifted by their UTC offset
    aware = pd.Series(
        pd.array(
            [pd.Timestamp("2016-01-05", tz="US/Eastern")],
            dtype=ArrowDtype(pa.timestamp("ns", tz="US/Eastern")),
        )
    )
    naive = pd.Series(
        pd.array([pd.Timestamp("2016-01-06")], dtype=ArrowDtype(pa.timestamp("ns")))
    )

    result = pd.concat([aware, naive], ignore_index=True)

    assert result.dtype == object
    assert result[0] == pd.Timestamp("2016-01-05", tz="US/Eastern")
    assert result[1] == pd.Timestamp("2016-01-06")


@pytest.mark.parametrize(
    "other, expected",
    [
        (ArrowDtype(pa.timestamp("ns")), None),
        (ArrowDtype(pa.timestamp("ns", tz="UTC")), None),
        (np.dtype("M8[ns]"), None),
        (
            ArrowDtype(pa.timestamp("us", tz="US/Eastern")),
            pa.timestamp("ns", tz="US/Eastern"),
        ),
        (
            pd.DatetimeTZDtype("ns", "US/Eastern"),
            pa.timestamp("ns", tz="US/Eastern"),
        ),
    ],
)
def test_get_common_dtype_timestamp_tz(other, expected):
    # GH#69029 disagreeing tz-awareness has no common ArrowDtype, matching
    #  datetime64; a shared tz survives a unit difference or a numpy spelling
    dtype = ArrowDtype(pa.timestamp("ns", tz="US/Eastern"))
    result = dtype._get_common_dtype([dtype, other])
    assert result == (None if expected is None else ArrowDtype(expected))


def test_concat_null_and_numpy_tz_aware():
    # GH#69029 a null[pyarrow] column names no type, so the tz-aware column decides
    #  the result; it arrow-ifies, as a numpy int64 column beside a null[pyarrow]
    #  one already does
    nullcol = pd.Series(pd.array([None], dtype=ArrowDtype(pa.null())))
    numpy_tz = pd.Series(
        pd.DatetimeIndex(["2016-01-06"]).tz_localize("US/Eastern").as_unit("us")
    )

    result = pd.concat([numpy_tz, nullcol], ignore_index=True)

    assert result.dtype == ArrowDtype(pa.timestamp("us", tz="US/Eastern"))
    assert result[0] == pd.Timestamp("2016-01-06", tz="US/Eastern")
    assert result[1] is pd.NA


def test_concat_timestamp_tz_arrow_and_numpy():
    # GH#69029 the tz was dropped on the way through numpy_dtype, so two equivalent
    #  tz-aware dtypes had no common dtype and the result fell back to object
    dtype = ArrowDtype(pa.timestamp("ns", tz="US/Eastern"))
    arrow = pd.Series(
        pd.array([pd.Timestamp("2016-01-05", tz="US/Eastern")], dtype=dtype)
    )
    numpy = pd.Series(
        pd.DatetimeIndex(["2016-01-06"]).tz_localize("US/Eastern").as_unit("ns")
    )

    result = pd.concat([arrow, numpy], ignore_index=True)

    assert result.dtype == dtype
    assert result[0] == pd.Timestamp("2016-01-05", tz="US/Eastern")
    assert result[1] == pd.Timestamp("2016-01-06", tz="US/Eastern")


def test_get_common_dtype_unresolvable_tz():
    # GH#69029 pa.timestamp does not validate its tz, so it can carry a label
    #  pandas cannot resolve; resolving a common dtype must not raise
    dtype = ArrowDtype(pa.timestamp("s", tz="Z"))
    other = ArrowDtype(pa.timestamp("ns", tz="Z"))

    assert dtype._get_common_dtype([dtype, other]) is None
    assert find_common_type([dtype, other]) == object


def test_concat_empty_arrow_backed_series(dtype):
    # GH#51734
    ser = pd.Series([], dtype=dtype)
    expected = ser.copy()
    result = pd.concat([ser[np.array([], dtype=np.bool_)]])
    tm.assert_series_equal(result, expected)


def test_concat_null_array():
    df = pd.DataFrame({"a": [None, None]}, dtype=ArrowDtype(pa.null()))
    df2 = pd.DataFrame({"a": [0, 1]}, dtype="int64[pyarrow]")

    result = pd.concat([df, df2], ignore_index=True)
    expected = pd.DataFrame({"a": [None, None, 0, 1]}, dtype="int64[pyarrow]")
    tm.assert_frame_equal(result, expected)


@pytest.mark.parametrize(
    "pa_type",
    [
        pa.date32(),
        pa.date64(),
        pa.time64("us"),
        pa.decimal128(7, 3),
        pa.binary(),
        pa.large_string(),
        pa.timestamp("us", "US/Pacific"),
        pa.list_(pa.int64()),
    ],
)
def test_concat_null_array_preserves_dtype(pa_type):
    # GH#62343 the null dtype should not affect the resulting dtype
    dtype = ArrowDtype(pa_type)
    ser = pd.Series([None], dtype=dtype)
    null_ser = pd.Series([None], dtype=ArrowDtype(pa.null()))

    result = pd.concat([ser, null_ser], ignore_index=True)
    expected = pd.Series([None, None], dtype=dtype)
    tm.assert_series_equal(result, expected)


def test_get_common_dtype_all_null():
    # GH#62343
    dtype = ArrowDtype(pa.null())
    assert dtype._get_common_dtype([dtype, dtype]) == dtype
