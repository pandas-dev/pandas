from datetime import (
    datetime,
    timedelta,
)

import numpy as np
import pytest

from pandas.core.dtypes.cast import (
    maybe_box_native,
    maybe_unbox_numpy_scalar,
)

from pandas import (
    Interval,
    Period,
    Timedelta,
    Timestamp,
)


@pytest.mark.parametrize(
    "obj,expected_dtype",
    [
        (b"\x00\x10", bytes),
        (4, int),
        (np.uint(4), int),
        (np.int32(-4), int),
        (np.uint8(4), int),
        (float(454.98), float),
        (np.float16(0.4), float),
        (np.float64(1.4), float),
        (np.bool_(False), bool),
        (datetime(2005, 2, 25), datetime),
        (np.datetime64("2005-02-25"), Timestamp),
        (Timestamp("2005-02-25"), Timestamp),
        (np.timedelta64(1, "D"), Timedelta),
        (Timedelta(1, "D"), Timedelta),
        (Interval(0, 1), Interval),
        (Period("4Q2005"), Period),
    ],
)
def test_maybe_box_native(obj, expected_dtype):
    boxed_obj = maybe_box_native(obj)
    result_dtype = type(boxed_obj)
    assert result_dtype is expected_dtype


@pytest.mark.parametrize("typecode", np.typecodes["All"])
def test_maybe_unbox_numpy_scalar(typecode, using_python_scalars):
    # https://github.com/pandas-dev/pandas/pull/63016
    if typecode == "?":
        scalar = False
        expected = bool
    elif typecode in "bhilqnpBHILQNP":
        scalar = 0
        expected = int
    elif typecode in "efdg":
        scalar = 0.0
        expected = float
    elif typecode in "FDG":
        scalar = 0.0 + 0.0j
        expected = complex
    elif typecode in "SV":
        scalar = b""
        expected = bytes
    elif typecode == "U":
        scalar = ""
        expected = str
    elif typecode == "O":
        scalar = 0
        expected = int
    elif typecode == "M":
        scalar = datetime(2025, 1, 1)
        expected = Timestamp
    elif typecode == "m":
        scalar = timedelta(seconds=3)
        expected = Timedelta
    else:
        raise ValueError(f"typecode {typecode} not recognized")
    value = np.array([scalar], dtype=typecode)[0]
    result = maybe_unbox_numpy_scalar(value)
    if using_python_scalars:
        assert type(result) == expected
    else:
        assert result is value


def test_maybe_unbox_numpy_scalar_timestamp(unit, using_python_scalars):
    # https://github.com/pandas-dev/pandas/pull/63016
    value = np.datetime64(1, unit)
    expected = Timestamp(1, unit=unit) if using_python_scalars else value
    result = maybe_unbox_numpy_scalar(value)
    assert result == expected
    assert type(result) == type(expected)


def test_maybe_unbox_numpy_scalar_datetime(unit, using_python_scalars):
    # https://github.com/pandas-dev/pandas/pull/63016
    value = np.timedelta64(1, unit)
    expected = Timedelta(1, unit=unit) if using_python_scalars else value
    result = maybe_unbox_numpy_scalar(value)
    assert result == expected
    assert type(result) == type(expected)
