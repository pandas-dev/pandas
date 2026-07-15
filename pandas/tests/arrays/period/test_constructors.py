import numpy as np
import pytest

from pandas._libs.tslibs import iNaT
from pandas._libs.tslibs.offsets import MonthEnd
from pandas._libs.tslibs.period import IncompatibleFrequency
from pandas.errors import Pandas4Warning

import pandas as pd
import pandas._testing as tm
from pandas.core.arrays import (
    PeriodArray,
    period_array,
)


@pytest.mark.parametrize(
    "data, freq, expected",
    [
        ([pd.Period("2017", "D")], None, [17167]),
        ([pd.Period("2017", "D")], "D", [17167]),
        ([2017], "D", [17167]),
        (["2017"], "D", [17167]),
        ([pd.Period("2017", "D")], pd.tseries.offsets.Day(), [17167]),
        ([pd.Period("2017", "D"), None], None, [17167, iNaT]),
        (pd.date_range("2017", periods=3), None, [17167, 17168, 17169]),
        (pd.period_range("2017", periods=4, freq="Q"), None, [188, 189, 190, 191]),
    ],
)
def test_period_array_ok(data, freq, expected):
    result = pd.PeriodIndex(data, freq=freq).array.asi8
    expected = np.asarray(expected, dtype=np.int64)
    tm.assert_numpy_array_equal(result, expected)

    dtype = pd.PeriodDtype(freq) if freq is not None else None
    result = PeriodArray._from_sequence(data, dtype=dtype)
    tm.assert_numpy_array_equal(result.asi8, expected)

    result = pd.PeriodIndex(data, dtype=dtype)
    tm.assert_numpy_array_equal(result.asi8, expected)

    result = pd.PeriodIndex(data, freq=freq)
    tm.assert_numpy_array_equal(result.asi8, expected)


def test_period_array_from_datetime_series_inferred_freq():
    # GH#64241
    data = pd.Series(pd.date_range("2017", periods=3))
    expected = np.asarray([17167, 17168, 17169], dtype=np.int64)

    msg = "Constructing PeriodArray from a Series of datetime64 data"
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        result = pd.PeriodIndex(data, freq=None).array.asi8
    tm.assert_numpy_array_equal(result, expected)

    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        result = PeriodArray._from_sequence(data, dtype=None)
    tm.assert_numpy_array_equal(result.asi8, expected)


def test_period_array_readonly_object():
    # https://github.com/pandas-dev/pandas/issues/25403
    pa = pd.PeriodIndex([pd.Period("2019-01-01")]).array
    arr = np.asarray(pa, dtype="object")
    arr.setflags(write=False)

    result = pd.PeriodIndex(arr).array
    tm.assert_period_array_equal(result, pa)

    result = pd.Series(arr)
    tm.assert_series_equal(result, pd.Series(pa))

    result = pd.DataFrame({"A": arr})
    tm.assert_frame_equal(result, pd.DataFrame({"A": pa}))


def test_from_datetime64_freq_changes():
    # https://github.com/pandas-dev/pandas/issues/23438
    arr = pd.date_range("2017", periods=3, freq="D")
    result = PeriodArray._from_datetime64(arr, freq="M")
    expected = pd.PeriodIndex(
        ["2017-01-01", "2017-01-01", "2017-01-01"], freq="M"
    ).array
    tm.assert_period_array_equal(result, expected)


@pytest.mark.parametrize("freq", ["2M", MonthEnd(2)])
def test_from_datetime64_freq_2M(freq):
    arr = np.array(
        ["2020-01-01T00:00:00", "2020-01-02T00:00:00"], dtype="datetime64[ns]"
    )
    result = PeriodArray._from_datetime64(arr, freq)
    expected = pd.PeriodIndex(["2020-01", "2020-01"], freq=freq).array
    tm.assert_period_array_equal(result, expected)


@pytest.mark.parametrize(
    "data, freq, msg",
    [
        (
            [pd.Period("2017", "D"), pd.Period("2017", "Y")],
            None,
            "Input has different freq",
        ),
        ([pd.Period("2017", "D")], "Y", "Input has different freq"),
    ],
)
def test_period_array_raises(data, freq, msg):
    with pytest.raises(IncompatibleFrequency, match=msg):
        pd.PeriodIndex(data, freq)


def test_period_array_non_period_series_raies():
    ser = pd.Series([1, 2, 3])
    with pytest.raises(TypeError, match="dtype"):
        PeriodArray(ser, dtype="period[D]")


def test_period_array_freq_mismatch():
    arr = pd.PeriodIndex(["2000", "2001"], freq="D").array
    with pytest.raises(IncompatibleFrequency, match="freq"):
        PeriodArray(arr, dtype="period[M]")

    dtype = pd.PeriodDtype(pd.tseries.offsets.MonthEnd())
    with pytest.raises(IncompatibleFrequency, match="freq"):
        PeriodArray(arr, dtype=dtype)


def test_from_sequence_allows_i8():
    # GH#64227 this used to be allowed for PeriodIndex and period_array
    # but not PeriodArray._from_sequence
    arr = period_array(["1975", "1976"], dtype="period[D]")

    expected = pd.PeriodIndex([pd.Period(x, freq=arr.freq) for x in arr.asi8]).array

    result1 = pd.PeriodIndex(arr.asi8, dtype=arr.dtype).array
    result2 = period_array(arr.asi8, dtype=arr.dtype)
    result3 = PeriodArray._from_sequence(arr.asi8, dtype=arr.dtype)
    result4 = PeriodArray._from_sequence(arr.asi8.astype(object), dtype=arr.dtype)
    result5 = PeriodArray._from_sequence(list(arr.asi8), dtype=arr.dtype)

    tm.assert_period_array_equal(result1, expected)
    tm.assert_period_array_equal(result2, expected)
    tm.assert_period_array_equal(result3, expected)
    tm.assert_period_array_equal(result4, expected)
    tm.assert_period_array_equal(result5, expected)


def test_from_sequence_integers_with_na_consistent():
    # GH#64227 an int is interpreted as a calendar year regardless of whether
    #  an NA is present; previously a None flipped ints to raw-ordinal values
    expected = PeriodArray._from_sequence([2000, 2001], dtype="period[D]")
    assert expected.tolist() == [pd.Period("2000", "D"), pd.Period("2001", "D")]

    result = PeriodArray._from_sequence([2000, None], dtype="period[D]")
    tm.assert_period_array_equal(result[:1], expected[:1])
    assert result[1] is pd.NaT

    # the integer NaT sentinel is still respected in the object path
    result = PeriodArray._from_sequence([iNaT, 2001], dtype="period[D]")
    assert result[0] is pd.NaT
    tm.assert_period_array_equal(result[1:], expected[1:])


@pytest.mark.parametrize("dtype", ["Int64", "UInt32", "int64[pyarrow]"])
def test_from_sequence_masked_arrow_integers_with_na(dtype):
    # GH#64227 masked/arrow integer arrays with NA used to raise a misleading
    #  "does not allow floating point" TypeError from the np.asarray-float path
    if "pyarrow" in dtype:
        pytest.importorskip("pyarrow")
    values = pd.array([2000, None], dtype=dtype)
    result = PeriodArray._from_sequence(values, dtype="period[D]")

    expected = PeriodArray._from_sequence(
        [pd.Period("2000", "D"), None], dtype="period[D]"
    )
    tm.assert_period_array_equal(result, expected)


def test_from_td64nat_sequence_raises():
    # GH#44507
    td = pd.NaT.to_numpy("m8[ns]")

    dtype = pd.period_range("2005-01-01", periods=3, freq="D").dtype

    arr = np.array([None], dtype=object)
    arr[0] = td

    msg = "Value must be Period, string, integer, or datetime"
    with pytest.raises(ValueError, match=msg):
        PeriodArray._from_sequence(arr, dtype=dtype)

    with pytest.raises(ValueError, match=msg):
        pd.PeriodIndex(arr, dtype=dtype)
    with pytest.raises(ValueError, match=msg):
        pd.Index(arr, dtype=dtype)
    with pytest.raises(ValueError, match=msg):
        pd.array(arr, dtype=dtype)
    with pytest.raises(ValueError, match=msg):
        pd.Series(arr, dtype=dtype)
    with pytest.raises(ValueError, match=msg):
        pd.DataFrame(arr, dtype=dtype)


def test_period_array_from_datetime64():
    arr = np.array(
        ["2020-01-01T00:00:00", "2020-02-02T00:00:00"], dtype="datetime64[ns]"
    )
    result = PeriodArray._from_datetime64(arr, freq=MonthEnd(2))

    expected = pd.PeriodIndex(["2020-01-01", "2020-02-01"], freq=MonthEnd(2)).array
    tm.assert_period_array_equal(result, expected)
