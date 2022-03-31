import numpy as np
import pytest

from pandas._libs.tslibs.np_datetime import (
    OutOfBoundsDatetime,
    astype_overflowsafe,
    py_get_unit_from_dtype,
    py_td64_to_tdstruct,
)

import pandas._testing as tm


def test_get_unit_from_dtype():
    # datetime64
    assert py_get_unit_from_dtype(np.dtype("M8[Y]")) == 0
    assert py_get_unit_from_dtype(np.dtype("M8[M]")) == 1
    assert py_get_unit_from_dtype(np.dtype("M8[W]")) == 2
    # B has been deprecated and removed -> no 3
    assert py_get_unit_from_dtype(np.dtype("M8[D]")) == 4
    assert py_get_unit_from_dtype(np.dtype("M8[h]")) == 5
    assert py_get_unit_from_dtype(np.dtype("M8[m]")) == 6
    assert py_get_unit_from_dtype(np.dtype("M8[s]")) == 7
    assert py_get_unit_from_dtype(np.dtype("M8[ms]")) == 8
    assert py_get_unit_from_dtype(np.dtype("M8[us]")) == 9
    assert py_get_unit_from_dtype(np.dtype("M8[ns]")) == 10
    assert py_get_unit_from_dtype(np.dtype("M8[ps]")) == 11
    assert py_get_unit_from_dtype(np.dtype("M8[fs]")) == 12
    assert py_get_unit_from_dtype(np.dtype("M8[as]")) == 13

    # timedelta64
    assert py_get_unit_from_dtype(np.dtype("m8[Y]")) == 0
    assert py_get_unit_from_dtype(np.dtype("m8[M]")) == 1
    assert py_get_unit_from_dtype(np.dtype("m8[W]")) == 2
    # B has been deprecated and removed -> no 3
    assert py_get_unit_from_dtype(np.dtype("m8[D]")) == 4
    assert py_get_unit_from_dtype(np.dtype("m8[h]")) == 5
    assert py_get_unit_from_dtype(np.dtype("m8[m]")) == 6
    assert py_get_unit_from_dtype(np.dtype("m8[s]")) == 7
    assert py_get_unit_from_dtype(np.dtype("m8[ms]")) == 8
    assert py_get_unit_from_dtype(np.dtype("m8[us]")) == 9
    assert py_get_unit_from_dtype(np.dtype("m8[ns]")) == 10
    assert py_get_unit_from_dtype(np.dtype("m8[ps]")) == 11
    assert py_get_unit_from_dtype(np.dtype("m8[fs]")) == 12
    assert py_get_unit_from_dtype(np.dtype("m8[as]")) == 13


def test_td64_to_tdstruct():
    val = 12454636234  # arbitrary value

    res1 = py_td64_to_tdstruct(val, 10)  # ns
    exp1 = {
        "days": 0,
        "hrs": 0,
        "min": 0,
        "sec": 12,
        "ms": 454,
        "us": 636,
        "ns": 234,
        "seconds": 12,
        "microseconds": 454636,
        "nanoseconds": 234,
    }
    assert res1 == exp1

    res2 = py_td64_to_tdstruct(val, 9)  # us
    exp2 = {
        "days": 0,
        "hrs": 3,
        "min": 27,
        "sec": 34,
        "ms": 636,
        "us": 234,
        "ns": 0,
        "seconds": 12454,
        "microseconds": 636234,
        "nanoseconds": 0,
    }
    assert res2 == exp2

    res3 = py_td64_to_tdstruct(val, 8)  # ms
    exp3 = {
        "days": 144,
        "hrs": 3,
        "min": 37,
        "sec": 16,
        "ms": 234,
        "us": 0,
        "ns": 0,
        "seconds": 13036,
        "microseconds": 234000,
        "nanoseconds": 0,
    }
    assert res3 == exp3

    # Note this out of bounds for nanosecond Timedelta
    res4 = py_td64_to_tdstruct(val, 7)  # s
    exp4 = {
        "days": 144150,
        "hrs": 21,
        "min": 10,
        "sec": 34,
        "ms": 0,
        "us": 0,
        "ns": 0,
        "seconds": 76234,
        "microseconds": 0,
        "nanoseconds": 0,
    }
    assert res4 == exp4


class TestAstypeOverflowSafe:
    def test_pass_non_dt64_array(self):
        # check that we raise, not segfault
        arr = np.arange(5)
        dtype = np.dtype("M8[ns]")

        msg = "astype_overflowsafe values must have datetime64 dtype"
        with pytest.raises(TypeError, match=msg):
            astype_overflowsafe(arr, dtype, copy=True)

        with pytest.raises(TypeError, match=msg):
            astype_overflowsafe(arr, dtype, copy=False)

    def test_pass_non_dt64_dtype(self):
        # check that we raise, not segfault
        arr = np.arange(5, dtype="i8").view("M8[D]")
        dtype = np.dtype("m8[ns]")

        msg = "astype_overflowsafe dtype must be datetime64"
        with pytest.raises(TypeError, match=msg):
            astype_overflowsafe(arr, dtype, copy=True)

        with pytest.raises(TypeError, match=msg):
            astype_overflowsafe(arr, dtype, copy=False)

    def test_astype_overflowsafe(self):
        dtype = np.dtype("M8[ns]")

        dt = np.datetime64("2262-04-05", "D")
        arr = dt + np.arange(10, dtype="m8[D]")

        # arr.astype silently overflows, so this
        wrong = arr.astype(dtype)
        roundtrip = wrong.astype(arr.dtype)
        assert not (wrong == roundtrip).all()

        msg = "Out of bounds nanosecond timestamp"
        with pytest.raises(OutOfBoundsDatetime, match=msg):
            astype_overflowsafe(arr, dtype)

        # But converting to microseconds is fine, and we match numpy's results.
        dtype2 = np.dtype("M8[us]")
        result = astype_overflowsafe(arr, dtype2)
        expected = arr.astype(dtype2)
        tm.assert_numpy_array_equal(result, expected)
