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


@pytest.fixture(params=tm.ALL_PYARROW_DTYPES, ids=str)
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


@pytest.fixture(params=["data", "data_missing"])
def all_data(request, data, data_missing):
    """Parametrized fixture returning 'data' or 'data_missing' integer arrays.

    Used to test dtype conversion with and without missing values.
    """
    if request.param == "data":
        return data
    elif request.param == "data_missing":
        return data_missing


@pytest.fixture
def na_value():
    """The scalar missing value for this type. Default 'None'"""
    return pd.NA


class TestBaseCasting(base.BaseCastingTests):
    pass


class TestConstructors(base.BaseConstructorsTests):
    def test_from_dtype(self, data, request):
        pa_dtype = data.dtype.pyarrow_dtype
        if pa.types.is_timestamp(pa_dtype) and pa_dtype.tz:
            if pa_version_under2p0:
                request.node.add_marker(
                    pytest.mark.xfail(
                        reason=f"timestamp data with tz={pa_dtype.tz} "
                        "converted to integer when pyarrow < 2.0",
                    )
                )
            else:
                request.node.add_marker(
                    pytest.mark.xfail(
                        raises=NotImplementedError,
                        reason=f"pyarrow.type_for_alias cannot infer {pa_dtype}",
                    )
                )
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


class TestBaseDtype(base.BaseDtypeTests):
    def test_construct_from_string_own_name(self, dtype, request):
        pa_dtype = dtype.pyarrow_dtype
        if pa.types.is_timestamp(pa_dtype) and pa_dtype.tz is not None:
            request.node.add_marker(
                pytest.mark.xfail(
                    raises=NotImplementedError,
                    reason=f"pyarrow.type_for_alias cannot infer {pa_dtype}",
                )
            )
        super().test_construct_from_string_own_name(dtype)

    def test_is_dtype_from_name(self, dtype, request):
        pa_dtype = dtype.pyarrow_dtype
        if pa.types.is_timestamp(pa_dtype) and pa_dtype.tz is not None:
            request.node.add_marker(
                pytest.mark.xfail(
                    raises=NotImplementedError,
                    reason=f"pyarrow.type_for_alias cannot infer {pa_dtype}",
                )
            )
        super().test_is_dtype_from_name(dtype)

    def test_construct_from_string(self, dtype, request):
        pa_dtype = dtype.pyarrow_dtype
        if pa.types.is_timestamp(pa_dtype) and pa_dtype.tz is not None:
            request.node.add_marker(
                pytest.mark.xfail(
                    raises=NotImplementedError,
                    reason=f"pyarrow.type_for_alias cannot infer {pa_dtype}",
                )
            )
        super().test_construct_from_string(dtype)

    def test_construct_from_string_another_type_raises(self, dtype):
        msg = r"'another_type' must end with '\[pyarrow\]'"
        with pytest.raises(TypeError, match=msg):
            type(dtype).construct_from_string("another_type")

    def test_get_common_dtype(self, dtype, request):
        pa_dtype = dtype.pyarrow_dtype
        if (
            pa.types.is_date(pa_dtype)
            or pa.types.is_time(pa_dtype)
            or (
                pa.types.is_timestamp(pa_dtype)
                and (pa_dtype.unit != "ns" or pa_dtype.tz is not None)
            )
            or (pa.types.is_duration(pa_dtype) and pa_dtype.unit != "ns")
        ):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=(
                        f"{pa_dtype} does not have associated numpy "
                        f"dtype findable by find_common_type"
                    )
                )
            )
        super().test_get_common_dtype(dtype)


class TestBaseIndex(base.BaseIndexTests):
    pass


class TestBaseInterface(base.BaseInterfaceTests):
    def test_contains(self, data, data_missing, request):
        tz = getattr(data.dtype.pyarrow_dtype, "tz", None)
        unit = getattr(data.dtype.pyarrow_dtype, "unit", None)
        if pa_version_under2p0 and tz not in (None, "UTC") and unit == "us":
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=(
                        f"Not supported by pyarrow < 2.0 "
                        f"with timestamp type {tz} and {unit}"
                    )
                )
            )
        super().test_contains(data, data_missing)

    @pytest.mark.xfail(reason="pyarrow.ChunkedArray does not support views.")
    def test_view(self, data):
        super().test_view(data)


class TestBaseMissing(base.BaseMissingTests):
    def test_fillna_limit_pad(self, data_missing, using_array_manager, request):
        if using_array_manager and pa.types.is_duration(
            data_missing.dtype.pyarrow_dtype
        ):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason="Checking ndim when using arraymanager with duration type"
                )
            )
        super().test_fillna_limit_pad(data_missing)

    def test_fillna_limit_backfill(self, data_missing, using_array_manager, request):
        if using_array_manager and pa.types.is_duration(
            data_missing.dtype.pyarrow_dtype
        ):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason="Checking ndim when using arraymanager with duration type"
                )
            )
        super().test_fillna_limit_backfill(data_missing)

    def test_fillna_series(self, data_missing, using_array_manager, request):
        if using_array_manager and pa.types.is_duration(
            data_missing.dtype.pyarrow_dtype
        ):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason="Checking ndim when using arraymanager with duration type"
                )
            )
        super().test_fillna_series(data_missing)

    def test_fillna_series_method(
        self, data_missing, fillna_method, using_array_manager, request
    ):
        if using_array_manager and pa.types.is_duration(
            data_missing.dtype.pyarrow_dtype
        ):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason="Checking ndim when using arraymanager with duration type"
                )
            )
        super().test_fillna_series_method(data_missing, fillna_method)

    def test_fillna_frame(self, data_missing, using_array_manager, request):
        if using_array_manager and pa.types.is_duration(
            data_missing.dtype.pyarrow_dtype
        ):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason="Checking ndim when using arraymanager with duration type"
                )
            )
        super().test_fillna_frame(data_missing)


class TestBaseSetitem(base.BaseSetitemTests):
    def test_setitem_scalar_series(self, data, box_in_series, request):
        tz = getattr(data.dtype.pyarrow_dtype, "tz", None)
        if pa_version_under2p0 and tz not in (None, "UTC"):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=(f"Not supported by pyarrow < 2.0 with timestamp type {tz}")
                )
            )
        super().test_setitem_scalar_series(data, box_in_series)

    def test_setitem_sequence(self, data, box_in_series, using_array_manager, request):
        tz = getattr(data.dtype.pyarrow_dtype, "tz", None)
        if pa_version_under2p0 and tz not in (None, "UTC"):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=(f"Not supported by pyarrow < 2.0 with timestamp type {tz}")
                )
            )
        elif (
            using_array_manager
            and pa.types.is_duration(data.dtype.pyarrow_dtype)
            and box_in_series
        ):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason="Checking ndim when using arraymanager with duration type"
                )
            )
        super().test_setitem_sequence(data, box_in_series)

    def test_setitem_sequence_mismatched_length_raises(
        self, data, as_array, using_array_manager, request
    ):
        if using_array_manager and pa.types.is_duration(data.dtype.pyarrow_dtype):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason="Checking ndim when using arraymanager with duration type"
                )
            )
        super().test_setitem_sequence_mismatched_length_raises(data, as_array)

    def test_setitem_empty_indexer(
        self, data, box_in_series, using_array_manager, request
    ):
        if (
            using_array_manager
            and pa.types.is_duration(data.dtype.pyarrow_dtype)
            and box_in_series
        ):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason="Checking ndim when using arraymanager with duration type"
                )
            )
        super().test_setitem_empty_indexer(data, box_in_series)

    def test_setitem_sequence_broadcasts(
        self, data, box_in_series, using_array_manager, request
    ):
        tz = getattr(data.dtype.pyarrow_dtype, "tz", None)
        if pa_version_under2p0 and tz not in (None, "UTC"):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=(f"Not supported by pyarrow < 2.0 with timestamp type {tz}")
                )
            )
        elif (
            using_array_manager
            and pa.types.is_duration(data.dtype.pyarrow_dtype)
            and box_in_series
        ):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason="Checking ndim when using arraymanager with duration type"
                )
            )
        super().test_setitem_sequence_broadcasts(data, box_in_series)

    @pytest.mark.parametrize("setter", ["loc", "iloc"])
    def test_setitem_scalar(self, data, setter, using_array_manager, request):
        tz = getattr(data.dtype.pyarrow_dtype, "tz", None)
        if pa_version_under2p0 and tz not in (None, "UTC"):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=(f"Not supported by pyarrow < 2.0 with timestamp type {tz}")
                )
            )
        elif using_array_manager and pa.types.is_duration(data.dtype.pyarrow_dtype):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason="Checking ndim when using arraymanager with duration type"
                )
            )
        super().test_setitem_scalar(data, setter)

    def test_setitem_loc_scalar_mixed(self, data, using_array_manager, request):
        tz = getattr(data.dtype.pyarrow_dtype, "tz", None)
        if pa_version_under2p0 and tz not in (None, "UTC"):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=(f"Not supported by pyarrow < 2.0 with timestamp type {tz}")
                )
            )
        elif using_array_manager and pa.types.is_duration(data.dtype.pyarrow_dtype):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason="Checking ndim when using arraymanager with duration type"
                )
            )
        super().test_setitem_loc_scalar_mixed(data)

    def test_setitem_loc_scalar_single(self, data, using_array_manager, request):
        tz = getattr(data.dtype.pyarrow_dtype, "tz", None)
        if pa_version_under2p0 and tz not in (None, "UTC"):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=(f"Not supported by pyarrow < 2.0 with timestamp type {tz}")
                )
            )
        elif using_array_manager and pa.types.is_duration(data.dtype.pyarrow_dtype):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason="Checking ndim when using arraymanager with duration type"
                )
            )
        super().test_setitem_loc_scalar_single(data)

    def test_setitem_loc_scalar_multiple_homogoneous(
        self, data, using_array_manager, request
    ):
        tz = getattr(data.dtype.pyarrow_dtype, "tz", None)
        if pa_version_under2p0 and tz not in (None, "UTC"):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=(f"Not supported by pyarrow < 2.0 with timestamp type {tz}")
                )
            )
        elif using_array_manager and pa.types.is_duration(data.dtype.pyarrow_dtype):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason="Checking ndim when using arraymanager with duration type"
                )
            )
        super().test_setitem_loc_scalar_multiple_homogoneous(data)

    def test_setitem_iloc_scalar_mixed(self, data, using_array_manager, request):
        tz = getattr(data.dtype.pyarrow_dtype, "tz", None)
        if pa_version_under2p0 and tz not in (None, "UTC"):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=(f"Not supported by pyarrow < 2.0 with timestamp type {tz}")
                )
            )
        elif using_array_manager and pa.types.is_duration(data.dtype.pyarrow_dtype):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason="Checking ndim when using arraymanager with duration type"
                )
            )
        super().test_setitem_iloc_scalar_mixed(data)

    def test_setitem_iloc_scalar_single(self, data, using_array_manager, request):
        tz = getattr(data.dtype.pyarrow_dtype, "tz", None)
        if pa_version_under2p0 and tz not in (None, "UTC"):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=(f"Not supported by pyarrow < 2.0 with timestamp type {tz}")
                )
            )
        elif using_array_manager and pa.types.is_duration(data.dtype.pyarrow_dtype):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason="Checking ndim when using arraymanager with duration type"
                )
            )
        super().test_setitem_iloc_scalar_single(data)

    def test_setitem_iloc_scalar_multiple_homogoneous(
        self, data, using_array_manager, request
    ):
        tz = getattr(data.dtype.pyarrow_dtype, "tz", None)
        if pa_version_under2p0 and tz not in (None, "UTC"):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=(f"Not supported by pyarrow < 2.0 with timestamp type {tz}")
                )
            )
        elif using_array_manager and pa.types.is_duration(data.dtype.pyarrow_dtype):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason="Checking ndim when using arraymanager with duration type"
                )
            )
        super().test_setitem_iloc_scalar_multiple_homogoneous(data)

    @pytest.mark.parametrize(
        "mask",
        [
            np.array([True, True, True, False, False]),
            pd.array([True, True, True, False, False], dtype="boolean"),
            pd.array([True, True, True, pd.NA, pd.NA], dtype="boolean"),
        ],
        ids=["numpy-array", "boolean-array", "boolean-array-na"],
    )
    def test_setitem_mask(
        self, data, mask, box_in_series, using_array_manager, request
    ):
        tz = getattr(data.dtype.pyarrow_dtype, "tz", None)
        if pa_version_under2p0 and tz not in (None, "UTC"):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=(f"Not supported by pyarrow < 2.0 with timestamp type {tz}")
                )
            )
        elif (
            using_array_manager
            and pa.types.is_duration(data.dtype.pyarrow_dtype)
            and box_in_series
        ):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason="Checking ndim when using arraymanager with duration type"
                )
            )
        super().test_setitem_mask(data, mask, box_in_series)

    def test_setitem_mask_boolean_array_with_na(
        self, data, box_in_series, using_array_manager, request
    ):
        tz = getattr(data.dtype.pyarrow_dtype, "tz", None)
        unit = getattr(data.dtype.pyarrow_dtype, "unit", None)
        if pa_version_under2p0 and tz not in (None, "UTC") and unit == "us":
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=(f"Not supported by pyarrow < 2.0 with timestamp type {tz}")
                )
            )
        elif (
            using_array_manager
            and pa.types.is_duration(data.dtype.pyarrow_dtype)
            and box_in_series
        ):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason="Checking ndim when using arraymanager with duration type"
                )
            )
        super().test_setitem_mask_boolean_array_with_na(data, box_in_series)

    @pytest.mark.parametrize(
        "idx",
        [[0, 1, 2], pd.array([0, 1, 2], dtype="Int64"), np.array([0, 1, 2])],
        ids=["list", "integer-array", "numpy-array"],
    )
    def test_setitem_integer_array(
        self, data, idx, box_in_series, using_array_manager, request
    ):
        tz = getattr(data.dtype.pyarrow_dtype, "tz", None)
        if pa_version_under2p0 and tz not in (None, "UTC"):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=(f"Not supported by pyarrow < 2.0 with timestamp type {tz}")
                )
            )
        elif (
            using_array_manager
            and pa.types.is_duration(data.dtype.pyarrow_dtype)
            and box_in_series
        ):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason="Checking ndim when using arraymanager with duration type"
                )
            )
        super().test_setitem_integer_array(data, idx, box_in_series)

    @pytest.mark.parametrize("as_callable", [True, False])
    @pytest.mark.parametrize("setter", ["loc", None])
    def test_setitem_mask_aligned(
        self, data, as_callable, setter, using_array_manager, request
    ):
        tz = getattr(data.dtype.pyarrow_dtype, "tz", None)
        if pa_version_under2p0 and tz not in (None, "UTC"):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=(f"Not supported by pyarrow < 2.0 with timestamp type {tz}")
                )
            )
        elif using_array_manager and pa.types.is_duration(data.dtype.pyarrow_dtype):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason="Checking ndim when using arraymanager with duration type"
                )
            )
        super().test_setitem_mask_aligned(data, as_callable, setter)

    @pytest.mark.parametrize("setter", ["loc", None])
    def test_setitem_mask_broadcast(self, data, setter, using_array_manager, request):
        tz = getattr(data.dtype.pyarrow_dtype, "tz", None)
        if pa_version_under2p0 and tz not in (None, "UTC"):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=(f"Not supported by pyarrow < 2.0 with timestamp type {tz}")
                )
            )
        elif using_array_manager and pa.types.is_duration(data.dtype.pyarrow_dtype):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason="Checking ndim when using arraymanager with duration type"
                )
            )
        super().test_setitem_mask_broadcast(data, setter)

    def test_setitem_tuple_index(self, data, request):
        tz = getattr(data.dtype.pyarrow_dtype, "tz", None)
        if pa_version_under2p0 and tz not in (None, "UTC"):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=(f"Not supported by pyarrow < 2.0 with timestamp type {tz}")
                )
            )
        super().test_setitem_tuple_index(data)

    def test_setitem_slice(self, data, box_in_series, using_array_manager, request):
        tz = getattr(data.dtype.pyarrow_dtype, "tz", None)
        if pa_version_under2p0 and tz not in (None, "UTC"):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=(f"Not supported by pyarrow < 2.0 with timestamp type {tz}")
                )
            )
        elif (
            using_array_manager
            and pa.types.is_duration(data.dtype.pyarrow_dtype)
            and box_in_series
        ):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason="Checking ndim when using arraymanager with duration type"
                )
            )
        super().test_setitem_slice(data, box_in_series)

    def test_setitem_loc_iloc_slice(self, data, using_array_manager, request):
        tz = getattr(data.dtype.pyarrow_dtype, "tz", None)
        if pa_version_under2p0 and tz not in (None, "UTC"):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=(f"Not supported by pyarrow < 2.0 with timestamp type {tz}")
                )
            )
        elif using_array_manager and pa.types.is_duration(data.dtype.pyarrow_dtype):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason="Checking ndim when using arraymanager with duration type"
                )
            )
        super().test_setitem_loc_iloc_slice(data)

    def test_setitem_slice_array(self, data, request):
        tz = getattr(data.dtype.pyarrow_dtype, "tz", None)
        if pa_version_under2p0 and tz not in (None, "UTC"):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=(f"Not supported by pyarrow < 2.0 with timestamp type {tz}")
                )
            )
        super().test_setitem_slice_array(data)

    def test_setitem_with_expansion_dataframe_column(
        self, data, full_indexer, using_array_manager, request
    ):
        # Is there a way to get the full_indexer id "null_slice"?
        is_null_slice = full_indexer(pd.Series(dtype=object)) == slice(None)
        tz = getattr(data.dtype.pyarrow_dtype, "tz", None)
        if pa_version_under2p0 and tz not in (None, "UTC") and not is_null_slice:
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=(f"Not supported by pyarrow < 2.0 with timestamp type {tz}")
                )
            )
        elif (
            using_array_manager
            and pa.types.is_duration(data.dtype.pyarrow_dtype)
            and not is_null_slice
        ):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason="Checking ndim when using arraymanager with duration type"
                )
            )
        super().test_setitem_with_expansion_dataframe_column(data, full_indexer)

    def test_setitem_with_expansion_row(
        self, data, na_value, using_array_manager, request
    ):
        tz = getattr(data.dtype.pyarrow_dtype, "tz", None)
        if pa_version_under2p0 and tz not in (None, "UTC"):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=(f"Not supported by pyarrow < 2.0 with timestamp type {tz}")
                )
            )
        elif using_array_manager and pa.types.is_duration(data.dtype.pyarrow_dtype):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason="Checking ndim when using arraymanager with duration type"
                )
            )
        super().test_setitem_with_expansion_row(data, na_value)

    def test_setitem_frame_2d_values(self, data, using_array_manager, request):
        tz = getattr(data.dtype.pyarrow_dtype, "tz", None)
        if pa_version_under2p0 and tz not in (None, "UTC"):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=(f"Not supported by pyarrow < 2.0 with timestamp type {tz}")
                )
            )
        elif using_array_manager and pa.types.is_duration(data.dtype.pyarrow_dtype):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason="Checking ndim when using arraymanager with duration type"
                )
            )
        super().test_setitem_frame_2d_values(data)

    @pytest.mark.xfail(reason="GH 45419: pyarrow.ChunkedArray does not support views")
    def test_setitem_preserves_views(self, data):
        super().test_setitem_preserves_views(data)


def test_arrowdtype_construct_from_string_type_with_unsupported_parameters():
    with pytest.raises(NotImplementedError, match="Passing pyarrow type"):
        ArrowDtype.construct_from_string("timestamp[s, tz=UTC][pyarrow]")
