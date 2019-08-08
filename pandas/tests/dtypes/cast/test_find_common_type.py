import numpy as np
import pytest

from pandas.core.dtypes.cast import find_common_type
from pandas.core.dtypes.dtypes import CategoricalDtype, DatetimeTZDtype, PeriodDtype


@pytest.mark.parametrize(
    "source_dtypes,expected_common_dtype",
    [
        ((np.int64,), np.int64),
        ((np.uint64,), np.uint64),
        ((np.float32,), np.float32),
        ((np.object,), np.object),
        # Into ints.
        ((np.int16, np.int64), np.int64),
        ((np.int32, np.uint32), np.int64),
        ((np.uint16, np.uint64), np.uint64),
        # Into floats.
        ((np.float16, np.float32), np.float32),
        ((np.float16, np.int16), np.float32),
        ((np.float32, np.int16), np.float32),
        ((np.uint64, np.int64), np.float64),
        ((np.int16, np.float64), np.float64),
        ((np.float16, np.int64), np.float64),
        # Into others.
        ((np.complex128, np.int32), np.complex128),
        ((np.object, np.float32), np.object),
        ((np.object, np.int16), np.object),
        # Bool with int.
        ((np.dtype("bool"), np.int64), np.object),
        ((np.dtype("bool"), np.int32), np.object),
        ((np.dtype("bool"), np.int16), np.object),
        ((np.dtype("bool"), np.int8), np.object),
        ((np.dtype("bool"), np.uint64), np.object),
        ((np.dtype("bool"), np.uint32), np.object),
        ((np.dtype("bool"), np.uint16), np.object),
        ((np.dtype("bool"), np.uint8), np.object),
        # Bool with float.
        ((np.dtype("bool"), np.float64), np.object),
        ((np.dtype("bool"), np.float32), np.object),
        (
            (np.dtype("datetime64[ns]"), np.dtype("datetime64[ns]")),
            np.dtype("datetime64[ns]"),
        ),
        (
            (np.dtype("timedelta64[ns]"), np.dtype("timedelta64[ns]")),
            np.dtype("timedelta64[ns]"),
        ),
        (
            (np.dtype("datetime64[ns]"), np.dtype("datetime64[ms]")),
            np.dtype("datetime64[ns]"),
        ),
        (
            (np.dtype("timedelta64[ms]"), np.dtype("timedelta64[ns]")),
            np.dtype("timedelta64[ns]"),
        ),
        ((np.dtype("datetime64[ns]"), np.dtype("timedelta64[ns]")), np.object),
        ((np.dtype("datetime64[ns]"), np.int64), np.object),
    ],
)
def test_numpy_dtypes(source_dtypes, expected_common_dtype):
    assert find_common_type(source_dtypes) == expected_common_dtype


def test_raises_empty_input():
    with pytest.raises(ValueError, match="no types given"):
        find_common_type([])


@pytest.mark.parametrize(
    "dtypes,exp_type",
    [
        ([CategoricalDtype()], "category"),
        ([np.object, CategoricalDtype()], np.object),
        ([CategoricalDtype(), CategoricalDtype()], "category"),
    ],
)
def test_categorical_dtype(dtypes, exp_type):
    assert find_common_type(dtypes) == exp_type


def test_datetimetz_dtype_match():
    dtype = DatetimeTZDtype(unit="ns", tz="US/Eastern")
    assert find_common_type([dtype, dtype]) == "datetime64[ns, US/Eastern]"


@pytest.mark.parametrize(
    "dtype2",
    [
        DatetimeTZDtype(unit="ns", tz="Asia/Tokyo"),
        np.dtype("datetime64[ns]"),
        np.object,
        np.int64,
    ],
)
def test_datetimetz_dtype_mismatch(dtype2):
    dtype = DatetimeTZDtype(unit="ns", tz="US/Eastern")
    assert find_common_type([dtype, dtype2]) == np.object
    assert find_common_type([dtype2, dtype]) == np.object


def test_period_dtype_match():
    dtype = PeriodDtype(freq="D")
    assert find_common_type([dtype, dtype]) == "period[D]"


@pytest.mark.parametrize(
    "dtype2",
    [
        DatetimeTZDtype(unit="ns", tz="Asia/Tokyo"),
        PeriodDtype(freq="2D"),
        PeriodDtype(freq="H"),
        np.dtype("datetime64[ns]"),
        np.object,
        np.int64,
    ],
)
def test_period_dtype_mismatch(dtype2):
    dtype = PeriodDtype(freq="D")
    assert find_common_type([dtype, dtype2]) == np.object
    assert find_common_type([dtype2, dtype]) == np.object
