import copy
import datetime
import zoneinfo

import numpy as np
import pytest

from pandas._libs.missing import is_matching_na

from pandas.core.dtypes.common import is_float

import pandas as pd


@pytest.mark.parametrize(
    "arr, idx",
    [
        ([1, 2, 3, 4], [0, 2, 1, 3]),
        ([1, np.nan, 3, np.nan], [0, 2, 1, 3]),
        (
            [1, np.nan, 3, np.nan],
            pd.MultiIndex.from_tuples([(0, "a"), (1, "b"), (2, "c"), (3, "c")]),
        ),
    ],
)
def test_equals(arr, idx):
    s1 = pd.Series(arr, index=idx)
    s2 = s1.copy()
    assert s1.equals(s2)

    s1[1] = 9
    assert not s1.equals(s2)


@pytest.mark.parametrize(
    "val", [1, 1.1, 1 + 1j, True, "abc", [1, 2], (1, 2), {1, 2}, {"a": 1}, None]
)
def test_equals_list_array(val):
    # GH20676 Verify equals operator for list of Numpy arrays
    arr = np.array([1, 2])
    s1 = pd.Series([arr, arr])
    s2 = s1.copy()
    assert s1.equals(s2)

    s1[1] = val
    assert not s1.equals(s2)


def test_equals_false_negative():
    # GH8437 Verify false negative behavior of equals function for dtype object
    arr = [False, np.nan]
    s1 = pd.Series(arr)
    s2 = s1.copy()
    s3 = pd.Series(index=range(2), dtype=object)
    s4 = s3.copy()
    s5 = s3.copy()
    s6 = s3.copy()

    s3[:-1] = s4[:-1] = s5[0] = s6[0] = False
    assert s1.equals(s1)
    assert s1.equals(s2)
    assert s1.equals(s3)
    assert s1.equals(s4)
    assert s1.equals(s5)
    assert s5.equals(s6)


def test_equals_matching_nas():
    # matching but not identical NAs
    left = pd.Series([np.datetime64("NaT", "ns")], dtype=object)
    right = pd.Series([np.datetime64("NaT", "ns")], dtype=object)
    assert left.equals(right)
    assert pd.Index(left).equals(pd.Index(right))
    assert left.array.equals(right.array)

    left = pd.Series([np.timedelta64("NaT", "ns")], dtype=object)
    right = pd.Series([np.timedelta64("NaT", "ns")], dtype=object)
    assert left.equals(right)
    assert pd.Index(left).equals(pd.Index(right))
    assert left.array.equals(right.array)

    left = pd.Series([np.float64("NaN")], dtype=object)
    right = pd.Series([np.float64("NaN")], dtype=object)
    assert left.equals(right)
    assert pd.Index(left, dtype=left.dtype).equals(pd.Index(right, dtype=right.dtype))
    assert left.array.equals(right.array)


def test_equals_mismatched_nas(nulls_fixture, nulls_fixture2):
    # GH#39650
    left = nulls_fixture
    right = nulls_fixture2
    if hasattr(right, "copy"):
        right = right.copy()
    else:
        right = copy.copy(right)

    ser = pd.Series([left], dtype=object)
    ser2 = pd.Series([right], dtype=object)

    if is_matching_na(left, right):
        assert ser.equals(ser2)
    elif (left is None and is_float(right)) or (right is None and is_float(left)):
        assert ser.equals(ser2)
    else:
        assert not ser.equals(ser2)


def test_equals_none_vs_nan():
    # GH#39650
    ser = pd.Series([1, None], dtype=object)
    ser2 = pd.Series([1, np.nan], dtype=object)

    assert ser.equals(ser2)
    assert pd.Index(ser, dtype=ser.dtype).equals(pd.Index(ser2, dtype=ser2.dtype))
    assert ser.array.equals(ser2.array)


@pytest.mark.parametrize("unit", ["s", "ms", "us"])
@pytest.mark.parametrize(
    "values",
    [
        pd.DatetimeIndex(["2016-01-01", "NaT", "2016-01-03"]),
        pd.DatetimeIndex(["2016-01-01", "NaT", "2016-01-03"], tz="US/Pacific"),
        pd.TimedeltaIndex(["1 day", "NaT", "3 days"]),
    ],
)
def test_equals_mismatched_datetimelike_resolution(values, unit):
    # GH#55694
    left = pd.Series(values.as_unit("ns"))
    right = pd.Series(values.as_unit(unit))
    assert left.dtype != right.dtype
    assert left.equals(right)
    assert right.equals(left)
    assert left.array.equals(right.array)

    assert not left.equals(right.iloc[::-1].reset_index(drop=True))


def test_equals_mismatched_resolution_equivalent_tz():
    # GH#55694 distinct but equivalent tz objects, as in DatetimeTZDtype.__eq__
    dti = pd.date_range("2016-01-01", periods=3)
    left = pd.Series(dti.tz_localize(datetime.UTC))
    right = pd.Series(dti.tz_localize(zoneinfo.ZoneInfo("UTC")).as_unit("s"))
    assert left.dt.tz is not right.dt.tz
    assert left.equals(right)
    assert right.equals(left)
    assert pd.Index(left).equals(pd.Index(right))


def test_equals_mismatched_datetimelike_resolution_not_equal():
    # GH#55694 a differing tz, or values not representable at the finer
    #  resolution, are never equal
    ser = pd.Series(pd.date_range("2016-01-01", periods=3, tz="UTC", unit="ns"))
    assert not ser.equals(ser.dt.tz_convert("US/Pacific").dt.as_unit("s"))
    assert not ser.equals(ser.dt.tz_localize(None).dt.as_unit("s"))

    far = pd.Series(pd.DatetimeIndex(["3000-01-01"]).as_unit("s"))
    assert not far.equals(pd.Series(pd.DatetimeIndex(["2000-01-01"]).as_unit("ns")))
    far_td = pd.Series(np.array([2**62], dtype="m8[s]"))
    assert not far_td.equals(pd.Series(np.array([1], dtype="m8[ns]")))


def test_equals_None_vs_float():
    # GH#44190
    left = pd.Series([-np.inf, np.nan, -1.0, 0.0, 1.0, 10 / 3, np.inf], dtype=object)
    right = pd.Series([None] * len(left))

    # these series were found to be equal due to a bug, check that they are correctly
    # found to not equal
    assert not left.equals(right)
    assert not right.equals(left)
    assert not left.to_frame().equals(right.to_frame())
    assert not right.to_frame().equals(left.to_frame())
    assert not pd.Index(left, dtype="object").equals(pd.Index(right, dtype="object"))
    assert not pd.Index(right, dtype="object").equals(pd.Index(left, dtype="object"))
