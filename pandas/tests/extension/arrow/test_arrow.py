from datetime import (
    date,
    datetime,
    time,
    timedelta,
)

import numpy as np
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
        data = [True, None, False, None, False, None]
    elif pa.types.is_floating(pa_dtype):
        data = [1.0, None, 0.0, None, -2.0, None, 0.5, None, 99.9, None]
    elif pa.types.is_signed_integer(pa_dtype):
        data = [1, None, 0, None, -2, None, 10]
    elif pa.types.is_unsigned_integer(pa_dtype):
        data = [1, None, 0, None, 2, None, 10]
    elif pa.types.is_date(pa_dtype):
        data = [
            date(2022, 1, 1),
            None,
            date(1999, 12, 31),
            None,
            date(2000, 1, 1),
            None,
        ]
    elif pa.types.is_timestamp(pa_dtype):
        data = [
            datetime(2020, 1, 1, 1, 1, 1, 1),
            None,
            datetime(1999, 1, 1, 1, 1, 1, 1),
            None,
            datetime(2000, 1, 1, 1, 1, 1, 1),
            None,
        ]
    elif pa.types.is_duration(pa_dtype):
        data = [timedelta(1), None, timedelta(1, 1), None, timedelta(-1), None]
    elif pa.types.is_time(pa_dtype):
        data = [time(12, 0), None, time(0, 12), None, time(0, 0), None]
    else:
        data = []
    return pd.array(data, dtype=dtype)


@pytest.fixture
def data_not_missing(data):
    data = data.take(
        indices=np.full(len(data), -1), allow_fill=True, fill_value=data[0]
    )
    return data


@pytest.fixture
def data_missing(data):
    """Length-2 array with [NA, Valid]"""
    return type(data)._from_sequence([data[1], data[0]])


@pytest.fixture
def na_value():
    """The scalar missing value for this type. Default 'None'"""
    return pd.NA


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

    def test_get(self, data_not_missing):
        super().test_get(data_not_missing)

    def test_take_sequence(self, data_not_missing):
        super().test_take_sequence(data_not_missing)

    def test_take(self, data_not_missing, na_value, na_cmp):
        super().test_take(data_not_missing, na_value, na_cmp)

    def test_take_non_na_fill_value(self, data_missing):
        super().test_take_non_na_fill_value(data_missing)

    def test_reindex_non_na_fill_value(self, data_missing):
        super().test_reindex_non_na_fill_value(data_missing)

    def test_take_series(self, request, data):
        tz = getattr(data._dtype.pyarrow_dtype, "tz", None)
        unit = getattr(data._dtype.pyarrow_dtype, "unit", None)
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
        tz = getattr(data._dtype.pyarrow_dtype, "tz", None)
        unit = getattr(data._dtype.pyarrow_dtype, "unit", None)
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
        tz = getattr(data._dtype.pyarrow_dtype, "tz", None)
        unit = getattr(data._dtype.pyarrow_dtype, "unit", None)
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
