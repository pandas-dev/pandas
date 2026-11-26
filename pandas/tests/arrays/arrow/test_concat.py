from __future__ import annotations

from datetime import UTC
import zoneinfo

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


def test_get_common_dtype_timestamp_tz_equivalent_different_units():
    arrow_us = ArrowDtype(pa.timestamp("us", tz="UTC"))
    arrow_ns = ArrowDtype(pa.timestamp("ns", tz="UTC"))
    zi_ns = pd.DatetimeTZDtype("ns", zoneinfo.ZoneInfo("UTC"))
    std_ns = pd.DatetimeTZDtype("ns", UTC)
    dateutil_utc_ns = pd.DatetimeTZDtype("ns", "dateutil/UTC")

    expected = arrow_ns
    for numpy_tz in [zi_ns, std_ns, dateutil_utc_ns]:
        assert find_common_type([arrow_us, numpy_tz]) == expected
        assert find_common_type([numpy_tz, arrow_us]) == expected


def test_concat_timestamp_tz_arrow_and_equivalent_numpy_different_units():
    arrow_ns = ArrowDtype(pa.timestamp("ns", tz="UTC"))
    s_pa = pd.Series(
        pd.array(
            [pd.Timestamp("2024-01-02", tz="UTC")],
            dtype=arrow_ns,
        )
    )
    expected = pd.Series(
        pd.array(
            [
                pd.Timestamp("2024-01-01", tz="UTC"),
                pd.Timestamp("2024-01-02", tz="UTC"),
            ],
            dtype=arrow_ns,
        )
    )
    expected_rev = pd.Series(
        pd.array(
            [
                pd.Timestamp("2024-01-02", tz="UTC"),
                pd.Timestamp("2024-01-01", tz="UTC"),
            ],
            dtype=arrow_ns,
        )
    )

    for tz in [zoneinfo.ZoneInfo("UTC"), "dateutil/UTC"]:
        s_np = pd.Series(pd.DatetimeIndex(["2024-01-01"], tz=tz).as_unit("us"))
        res1 = pd.concat([s_np, s_pa], ignore_index=True)
        res2 = pd.concat([s_pa, s_np], ignore_index=True)
        tm.assert_series_equal(res1, expected)
        tm.assert_series_equal(res2, expected_rev)


def test_get_common_dtype_non_utc_dateutil_falls_back():
    dateutil_east = pd.DatetimeTZDtype("us", "dateutil/US/Eastern")
    null_arrow = ArrowDtype(pa.null())
    assert null_arrow._get_common_dtype([null_arrow, dateutil_east]) is None
    assert find_common_type([null_arrow, dateutil_east]) == object


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
