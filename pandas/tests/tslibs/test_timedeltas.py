import re

import numpy as np
import pytest

from pandas._libs.tslibs.timedeltas import (
    _python_get_unit_from_dtype,
    array_to_timedelta64,
    delta_to_nanoseconds,
    ints_to_pytimedelta,
)

from pandas import (
    Timedelta,
    offsets,
)
import pandas._testing as tm


@pytest.mark.parametrize(
    "obj,expected",
    [
        (np.timedelta64(14, "D"), 14 * 24 * 3600 * 1e9),
        (Timedelta(minutes=-7), -7 * 60 * 1e9),
        (Timedelta(minutes=-7).to_pytimedelta(), -7 * 60 * 1e9),
        (Timedelta(seconds=1234e-9), 1234),  # GH43764, GH40946
        (
            Timedelta(seconds=1e-9, milliseconds=1e-5, microseconds=1e-1),
            111,
        ),  # GH43764
        (
            Timedelta(days=1, seconds=1e-9, milliseconds=1e-5, microseconds=1e-1),
            24 * 3600e9 + 111,
        ),  # GH43764
        (offsets.Nano(125), 125),
    ],
)
def test_delta_to_nanoseconds(obj, expected):
    result = delta_to_nanoseconds(obj)
    assert result == expected


def test_delta_to_nanoseconds_error():
    obj = np.array([123456789], dtype="m8[ns]")

    with pytest.raises(TypeError, match="<class 'numpy.ndarray'>"):
        delta_to_nanoseconds(obj)

    with pytest.raises(TypeError, match="float"):
        delta_to_nanoseconds(1.5)
    with pytest.raises(TypeError, match="int"):
        delta_to_nanoseconds(1)
    with pytest.raises(TypeError, match="int"):
        delta_to_nanoseconds(np.int64(2))
    with pytest.raises(TypeError, match="int"):
        delta_to_nanoseconds(np.int32(3))


def test_delta_to_nanoseconds_td64_MY_raises():
    msg = (
        "delta_to_nanoseconds does not support Y or M units, "
        "as their duration in nanoseconds is ambiguous"
    )

    td = np.timedelta64(1234, "Y")

    with pytest.raises(ValueError, match=msg):
        delta_to_nanoseconds(td)

    td = np.timedelta64(1234, "M")

    with pytest.raises(ValueError, match=msg):
        delta_to_nanoseconds(td)


@pytest.mark.parametrize("unit", ["Y", "M"])
def test_unsupported_td64_unit_raises(unit):
    # GH 52806
    with pytest.raises(
        ValueError,
        match=f"Unit {unit} is not supported. "
        "Only unambiguous timedelta values durations are supported. "
        "Allowed units are 'W', 'D', 'h', 'm', 's', 'ms', 'us', 'ns'",
    ):
        Timedelta(np.timedelta64(1, unit))


def test_huge_nanoseconds_overflow():
    # GH 32402
    assert delta_to_nanoseconds(Timedelta(1e10)) == 1e10
    assert delta_to_nanoseconds(Timedelta(nanoseconds=1e10)) == 1e10


@pytest.mark.parametrize(
    "kwargs", [{"Seconds": 1}, {"seconds": 1, "Nanoseconds": 1}, {"Foo": 2}]
)
def test_kwarg_assertion(kwargs):
    err_message = (
        "cannot construct a Timedelta from the passed arguments, "
        "allowed keywords are "
        "[weeks, days, hours, minutes, seconds, "
        "milliseconds, microseconds, nanoseconds]"
    )

    with pytest.raises(ValueError, match=re.escape(err_message)):
        Timedelta(**kwargs)


class TestArrayToTimedelta64:
    def test_array_to_timedelta64_string_with_unit_2d_raises(self):
        # check the 'unit is not None and errors != "coerce"' path
        #  in array_to_timedelta64 raises correctly with 2D values
        values = np.array([["1", 2], [3, "4"]], dtype=object)
        with pytest.raises(ValueError, match="unit must not be specified"):
            array_to_timedelta64(values, unit="s")

    def test_array_to_timedelta64_non_object_raises(self):
        # check we raise, not segfault
        values = np.arange(5)

        msg = "'values' must have object dtype"
        with pytest.raises(TypeError, match=msg):
            array_to_timedelta64(values)


@pytest.mark.parametrize("unit", ["s", "ms", "us"])
def test_ints_to_pytimedelta(unit):
    # tests for non-nanosecond cases
    arr = np.arange(6, dtype=np.int64).view(f"m8[{unit}]")

    res = ints_to_pytimedelta(arr, box=False)
    # For non-nanosecond, .astype(object) gives pytimedelta objects
    #  instead of integers
    expected = arr.astype(object)
    tm.assert_numpy_array_equal(res, expected)

    res = ints_to_pytimedelta(arr, box=True)
    expected = np.array([Timedelta(x) for x in arr], dtype=object)
    tm.assert_numpy_array_equal(res, expected)


@pytest.mark.parametrize("unit", ["Y", "M", "ps", "fs", "as"])
def test_ints_to_pytimedelta_unsupported(unit):
    arr = np.arange(6, dtype=np.int64).view(f"m8[{unit}]")

    with pytest.raises(NotImplementedError, match=r"\d{1,2}"):
        ints_to_pytimedelta(arr, box=False)
    msg = "Only resolutions 's', 'ms', 'us', 'ns' are supported"
    with pytest.raises(NotImplementedError, match=msg):
        ints_to_pytimedelta(arr, box=True)


@pytest.mark.parametrize(
    "original_unit,value,new_unit,expected_pytime_days,expected_pytime_seconds,expected_pytime_microseconds",
    [
        ["days", 1, "s", 1, 0, 0],
        ["seconds", 1, "s", 0, 1, 0],
        ["milliseconds", 1, "ms", 0, 0, 1000],
        ["microseconds", 1, "us", 0, 0, 1],
        ["days", 1, "ms", 1, 0, 0],
        ["days", 1, "us", 1, 0, 0],
    ],
)
def test_timedelta_as_unit_conversion(
    original_unit,
    value,
    new_unit,
    expected_pytime_days,
    expected_pytime_seconds,
    expected_pytime_microseconds,
):
    kwargs = {original_unit: value}
    orig_timedelta = Timedelta(**kwargs)

    converted_timedelta = orig_timedelta.as_unit(new_unit)

    assert orig_timedelta._get_pytimedelta_days() == expected_pytime_days
    assert converted_timedelta._get_pytimedelta_days() == expected_pytime_days
    assert orig_timedelta._get_pytimedelta_seconds() == expected_pytime_seconds
    assert converted_timedelta._get_pytimedelta_seconds() == expected_pytime_seconds
    assert (
        orig_timedelta._get_pytimedelta_microseconds() == expected_pytime_microseconds
    )
    assert (
        converted_timedelta._get_pytimedelta_microseconds()
        == expected_pytime_microseconds
    )


@pytest.mark.parametrize(
    "unit,value,expected_pytime_days,expected_pytime_seconds,expected_pytime_microseconds",
    [
        ["s", 1, 0, 1, 0],
        ["s", 86400, 1, 0, 0],
        ["s", 86401, 1, 1, 0],
        ["ms", 1, 0, 0, 1],
        ["ms", 1000, 0, 1, 0],
        ["ms", 86400000, 1, 0, 0],
        ["ms", 86401001, 1, 1, 1],
    ],
)
def test_non_nano_c_api(
    unit,
    value,
    expected_pytime_days,
    expected_pytime_seconds,
    expected_pytime_microseconds,
):
    tmp_delta = Timedelta(days=1)  # just for exposing the function under test
    dtype = np.dtype(f"m8[{unit}]")
    reso = _python_get_unit_from_dtype(dtype)

    uut_delta = tmp_delta._from_value_and_reso(value, reso)

    assert uut_delta._get_pytimedelta_days() == expected_pytime_days
    assert uut_delta._get_pytimedelta_seconds() == expected_pytime_seconds
    assert uut_delta._get_pytimedelta_microseconds() == expected_pytime_microseconds


@pytest.mark.parametrize(
    "unit,value",
    [
        ["s", 86400000000000],  # (86400 s/day)*(1 billion days)
        ["ms", 86400000000000000],  # (86400000 ms/day)*(1 billion days)
        ["s", -86400000000000],
        ["ms", -86400000000000000],
    ],
)
def test_non_nano_c_api_pytimedelta_overflow(unit, value):
    tmp_delta = Timedelta(days=1)  # just for exposing the function under test
    dtype = np.dtype(f"m8[{unit}]")
    reso = _python_get_unit_from_dtype(dtype)

    uut_delta = tmp_delta._from_value_and_reso(value, reso)

    assert uut_delta._get_pytimedelta_days() == 0
    assert uut_delta._get_pytimedelta_seconds() == 0
    assert uut_delta._get_pytimedelta_microseconds() == 0
