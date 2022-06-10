"""
This file contains a minimal set of tests for compliance with the extension
array interface test suite, and should contain no other tests.
The test suite for the full functionality of the array is located in
`pandas/tests/arrays/`.
The tests in this file are inherited from the BaseExtensionTests, and only
minimal tweaks should be applied to get the tests passing (by overwriting a
parent method).
Additional tests should either be added to one of the BaseExtensionTests
classes (if they are relevant for the extension interface for all dtypes), or
be added to the array-specific tests in `pandas/tests/arrays/`.
"""

from datetime import (
    date,
    datetime,
    time,
    timedelta,
)

import pytest

from pandas.compat import (
    pa_version_under2p0,
    pa_version_under3p0,
)

import pandas as pd
import pandas._testing as tm
from pandas.tests.extension import base

pa = pytest.importorskip("pyarrow", minversion="1.0.1")

from pandas.core.arrays.arrow.dtype import ArrowDtype  # isort:skip


@pytest.fixture(params=tm.ALL_PYARROW_DTYPES)
def dtype(request):
    return ArrowDtype(pyarrow_dtype=request.param)


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
    else:
        raise NotImplementedError
    return pd.array(data, dtype=dtype)


@pytest.fixture
def data_missing(data):
    """Length-2 array with [NA, Valid]"""
    return type(data)._from_sequence([None, data[0]])


@pytest.fixture
def na_value():
    """The scalar missing value for this type. Default 'None'"""
    return pd.NA


class TestBaseCasting(base.BaseCastingTests):
    pass


class TestConstructors(base.BaseConstructorsTests):
    @pytest.mark.xfail(
        reason=(
            "str(dtype) constructs "
            "e.g. in64[pyarrow] like int64 (numpy) "
            "due to StorageExtensionDtype.__str__"
        )
    )
    def test_from_dtype(self, data):
        super().test_from_dtype(data)


@pytest.mark.xfail(
    raises=NotImplementedError, reason="pyarrow.ChunkedArray backing is 1D."
)
class TestDim2Compat(base.Dim2CompatTests):
    pass


@pytest.mark.xfail(
    raises=NotImplementedError, reason="pyarrow.ChunkedArray backing is 1D."
)
class TestNDArrayBacked2D(base.NDArrayBacked2DTests):
    pass


class TestGetitemTests(base.BaseGetitemTests):
    @pytest.mark.xfail(
        reason=(
            "data.dtype.type return pyarrow.DataType "
            "but this (intentionally) returns "
            "Python scalars or pd.Na"
        )
    )
    def test_getitem_scalar(self, data):
        super().test_getitem_scalar(data)

    def test_take_series(self, request, data):
        tz = getattr(data.dtype.pyarrow_dtype, "tz", None)
        unit = getattr(data.dtype.pyarrow_dtype, "unit", None)
        bad_units = ["ns"]
        if pa_version_under2p0:
            bad_units.extend(["s", "ms", "us"])
        if pa_version_under3p0 and tz not in (None, "UTC") and unit in bad_units:
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=(
                        f"Not supported by pyarrow < 3.0 "
                        f"with timestamp type {tz} and {unit}"
                    )
                )
            )
        super().test_take_series(data)

    def test_reindex(self, request, data, na_value):
        tz = getattr(data.dtype.pyarrow_dtype, "tz", None)
        unit = getattr(data.dtype.pyarrow_dtype, "unit", None)
        bad_units = ["ns"]
        if pa_version_under2p0:
            bad_units.extend(["s", "ms", "us"])
        if pa_version_under3p0 and tz not in (None, "UTC") and unit in bad_units:
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=(
                        f"Not supported by pyarrow < 3.0 "
                        f"with timestamp type {tz} and {unit}"
                    )
                )
            )
        super().test_reindex(data, na_value)

    def test_loc_iloc_frame_single_dtype(self, request, using_array_manager, data):
        tz = getattr(data.dtype.pyarrow_dtype, "tz", None)
        unit = getattr(data.dtype.pyarrow_dtype, "unit", None)
        bad_units = ["ns"]
        if pa_version_under2p0:
            bad_units.extend(["s", "ms", "us"])
        if (
            pa_version_under3p0
            and not using_array_manager
            and tz not in (None, "UTC")
            and unit in bad_units
        ):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=(
                        f"Not supported by pyarrow < 3.0 "
                        f"with timestamp type {tz} and {unit}"
                    )
                )
            )
        super().test_loc_iloc_frame_single_dtype(data)


class TestBaseIndex(base.BaseIndexTests):
    pass


def test_arrowdtype_construct_from_string_type_with_parameters():
    with pytest.raises(NotImplementedError, match="Passing pyarrow type"):
        ArrowDtype.construct_from_string("timestamp[s][pyarrow]")
