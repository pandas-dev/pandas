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
from io import (
    BytesIO,
    StringIO,
)
import pickle

import numpy as np
import pytest

from pandas.compat import (
    PY311,
    is_ci_environment,
    is_platform_windows,
    pa_version_under6p0,
    pa_version_under7p0,
    pa_version_under8p0,
    pa_version_under9p0,
)
from pandas.errors import PerformanceWarning

from pandas.core.dtypes.common import is_any_int_dtype

import pandas as pd
import pandas._testing as tm
from pandas.api.types import (
    is_bool_dtype,
    is_float_dtype,
    is_integer_dtype,
    is_numeric_dtype,
    is_signed_integer_dtype,
    is_unsigned_integer_dtype,
)
from pandas.tests.extension import base

pa = pytest.importorskip("pyarrow", minversion="6.0.0")

from pandas.core.arrays.arrow.array import ArrowExtensionArray

from pandas.core.arrays.arrow.dtype import ArrowDtype  # isort:skip

pytestmark = pytest.mark.filterwarnings(
    "ignore:.* may decrease performance. Upgrade to pyarrow >=7 to possibly"
)


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
    elif pa.types.is_string(pa_dtype):
        data = ["a", "b"] * 4 + [None] + ["1", "2"] * 44 + [None] + ["!", ">"]
    elif pa.types.is_binary(pa_dtype):
        data = [b"a", b"b"] * 4 + [None] + [b"1", b"2"] * 44 + [None] + [b"!", b">"]
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
def data_for_grouping(dtype):
    """
    Data for factorization, grouping, and unique tests.

    Expected to be like [B, B, NA, NA, A, A, B, C]

    Where A < B < C and NA is missing
    """
    pa_dtype = dtype.pyarrow_dtype
    if pa.types.is_boolean(pa_dtype):
        A = False
        B = True
        C = True
    elif pa.types.is_floating(pa_dtype):
        A = -1.1
        B = 0.0
        C = 1.1
    elif pa.types.is_signed_integer(pa_dtype):
        A = -1
        B = 0
        C = 1
    elif pa.types.is_unsigned_integer(pa_dtype):
        A = 0
        B = 1
        C = 10
    elif pa.types.is_date(pa_dtype):
        A = date(1999, 12, 31)
        B = date(2010, 1, 1)
        C = date(2022, 1, 1)
    elif pa.types.is_timestamp(pa_dtype):
        A = datetime(1999, 1, 1, 1, 1, 1, 1)
        B = datetime(2020, 1, 1)
        C = datetime(2020, 1, 1, 1)
    elif pa.types.is_duration(pa_dtype):
        A = timedelta(-1)
        B = timedelta(0)
        C = timedelta(1, 4)
    elif pa.types.is_time(pa_dtype):
        A = time(0, 0)
        B = time(0, 12)
        C = time(12, 12)
    elif pa.types.is_string(pa_dtype):
        A = "a"
        B = "b"
        C = "c"
    elif pa.types.is_binary(pa_dtype):
        A = b"a"
        B = b"b"
        C = b"c"
    else:
        raise NotImplementedError
    return pd.array([B, B, None, None, A, A, B, C], dtype=dtype)


@pytest.fixture
def data_for_sorting(data_for_grouping):
    """
    Length-3 array with a known sort order.

    This should be three items [B, C, A] with
    A < B < C
    """
    return type(data_for_grouping)._from_sequence(
        [data_for_grouping[0], data_for_grouping[7], data_for_grouping[4]]
    )


@pytest.fixture
def data_missing_for_sorting(data_for_grouping):
    """
    Length-3 array with a known sort order.

    This should be three items [B, NA, A] with
    A < B and NA missing.
    """
    return type(data_for_grouping)._from_sequence(
        [data_for_grouping[0], data_for_grouping[2], data_for_grouping[4]]
    )


@pytest.fixture
def data_for_twos(data):
    """Length-100 array in which all the elements are two."""
    pa_dtype = data.dtype.pyarrow_dtype
    if pa.types.is_integer(pa_dtype) or pa.types.is_floating(pa_dtype):
        return pd.array([2] * 100, dtype=data.dtype)
    # tests will be xfailed where 2 is not a valid scalar for pa_dtype
    return data


@pytest.fixture
def na_value():
    """The scalar missing value for this type. Default 'None'"""
    return pd.NA


class TestBaseCasting(base.BaseCastingTests):
    def test_astype_str(self, data, request):
        pa_dtype = data.dtype.pyarrow_dtype
        if pa.types.is_binary(pa_dtype):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=f"For {pa_dtype} .astype(str) decodes.",
                )
            )
        super().test_astype_str(data)


class TestConstructors(base.BaseConstructorsTests):
    def test_from_dtype(self, data, request):
        pa_dtype = data.dtype.pyarrow_dtype
        if (pa.types.is_timestamp(pa_dtype) and pa_dtype.tz) or pa.types.is_string(
            pa_dtype
        ):
            if pa.types.is_string(pa_dtype):
                reason = "ArrowDtype(pa.string()) != StringDtype('pyarrow')"
            else:
                reason = f"pyarrow.type_for_alias cannot infer {pa_dtype}"
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=reason,
                )
            )
        super().test_from_dtype(data)

    def test_from_sequence_pa_array(self, data, request):
        # https://github.com/pandas-dev/pandas/pull/47034#discussion_r955500784
        # data._data = pa.ChunkedArray
        result = type(data)._from_sequence(data._data)
        tm.assert_extension_array_equal(result, data)
        assert isinstance(result._data, pa.ChunkedArray)

        result = type(data)._from_sequence(data._data.combine_chunks())
        tm.assert_extension_array_equal(result, data)
        assert isinstance(result._data, pa.ChunkedArray)

    def test_from_sequence_pa_array_notimplemented(self, request):
        if pa_version_under6p0:
            request.node.add_marker(
                pytest.mark.xfail(
                    raises=AttributeError,
                    reason="month_day_nano_interval not implemented by pyarrow.",
                )
            )
        with pytest.raises(NotImplementedError, match="Converting strings to"):
            ArrowExtensionArray._from_sequence_of_strings(
                ["12-1"], dtype=pa.month_day_nano_interval()
            )

    def test_from_sequence_of_strings_pa_array(self, data, request):
        pa_dtype = data.dtype.pyarrow_dtype
        if pa.types.is_time64(pa_dtype) and pa_dtype.equals("time64[ns]") and not PY311:
            request.node.add_marker(
                pytest.mark.xfail(
                    reason="Nanosecond time parsing not supported.",
                )
            )
        elif pa.types.is_duration(pa_dtype):
            request.node.add_marker(
                pytest.mark.xfail(
                    raises=pa.ArrowNotImplementedError,
                    reason=f"pyarrow doesn't support parsing {pa_dtype}",
                )
            )
        elif pa.types.is_timestamp(pa_dtype) and pa_dtype.tz is not None:
            if pa_version_under7p0:
                request.node.add_marker(
                    pytest.mark.xfail(
                        raises=pa.ArrowNotImplementedError,
                        reason=f"pyarrow doesn't support string cast from {pa_dtype}",
                    )
                )
            elif is_platform_windows() and is_ci_environment():
                request.node.add_marker(
                    pytest.mark.xfail(
                        raises=pa.ArrowInvalid,
                        reason=(
                            "TODO: Set ARROW_TIMEZONE_DATABASE environment variable "
                            "on CI to path to the tzdata for pyarrow."
                        ),
                    )
                )
        elif pa_version_under6p0 and pa.types.is_temporal(pa_dtype):
            request.node.add_marker(
                pytest.mark.xfail(
                    raises=pa.ArrowNotImplementedError,
                    reason=f"pyarrow doesn't support string cast from {pa_dtype}",
                )
            )
        pa_array = data._data.cast(pa.string())
        result = type(data)._from_sequence_of_strings(pa_array, dtype=data.dtype)
        tm.assert_extension_array_equal(result, data)

        pa_array = pa_array.combine_chunks()
        result = type(data)._from_sequence_of_strings(pa_array, dtype=data.dtype)
        tm.assert_extension_array_equal(result, data)


class TestGetitemTests(base.BaseGetitemTests):
    def test_getitem_scalar(self, data):
        # In the base class we expect data.dtype.type; but this (intentionally)
        #  returns Python scalars or pd.NA
        pa_type = data._data.type
        if pa.types.is_integer(pa_type):
            exp_type = int
        elif pa.types.is_floating(pa_type):
            exp_type = float
        elif pa.types.is_string(pa_type):
            exp_type = str
        elif pa.types.is_binary(pa_type):
            exp_type = bytes
        elif pa.types.is_boolean(pa_type):
            exp_type = bool
        elif pa.types.is_duration(pa_type):
            exp_type = timedelta
        elif pa.types.is_timestamp(pa_type):
            if pa_type.unit == "ns":
                exp_type = pd.Timestamp
            else:
                exp_type = datetime
        elif pa.types.is_date(pa_type):
            exp_type = date
        elif pa.types.is_time(pa_type):
            exp_type = time
        else:
            raise NotImplementedError(data.dtype)

        result = data[0]
        assert isinstance(result, exp_type), type(result)

        result = pd.Series(data)[0]
        assert isinstance(result, exp_type), type(result)


class TestBaseAccumulateTests(base.BaseAccumulateTests):
    def check_accumulate(self, s, op_name, skipna):
        result = getattr(s, op_name)(skipna=skipna).astype("Float64")
        expected = getattr(s.astype("Float64"), op_name)(skipna=skipna)
        self.assert_series_equal(result, expected, check_dtype=False)

    @pytest.mark.parametrize("skipna", [True, False])
    def test_accumulate_series_raises(self, data, all_numeric_accumulations, skipna):
        pa_type = data.dtype.pyarrow_dtype
        if (
            (pa.types.is_integer(pa_type) or pa.types.is_floating(pa_type))
            and all_numeric_accumulations == "cumsum"
            and not pa_version_under9p0
        ):
            pytest.skip("These work, are tested by test_accumulate_series.")

        op_name = all_numeric_accumulations
        ser = pd.Series(data)

        with pytest.raises(NotImplementedError):
            getattr(ser, op_name)(skipna=skipna)

    @pytest.mark.parametrize("skipna", [True, False])
    def test_accumulate_series(self, data, all_numeric_accumulations, skipna, request):
        pa_type = data.dtype.pyarrow_dtype
        op_name = all_numeric_accumulations
        ser = pd.Series(data)

        do_skip = False
        if pa.types.is_string(pa_type) or pa.types.is_binary(pa_type):
            if op_name in ["cumsum", "cumprod"]:
                do_skip = True
        elif pa.types.is_temporal(pa_type) and not pa.types.is_duration(pa_type):
            if op_name in ["cumsum", "cumprod"]:
                do_skip = True
        elif pa.types.is_duration(pa_type):
            if op_name == "cumprod":
                do_skip = True

        if do_skip:
            pytest.skip(
                "These should *not* work, we test in test_accumulate_series_raises "
                "that these correctly raise."
            )

        if all_numeric_accumulations != "cumsum" or pa_version_under9p0:
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=f"{all_numeric_accumulations} not implemented",
                    raises=NotImplementedError,
                )
            )
        elif all_numeric_accumulations == "cumsum" and (
            pa.types.is_duration(pa_type) or pa.types.is_boolean(pa_type)
        ):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=f"{all_numeric_accumulations} not implemented for {pa_type}",
                    raises=NotImplementedError,
                )
            )

        self.check_accumulate(ser, op_name, skipna)


class TestBaseNumericReduce(base.BaseNumericReduceTests):
    def check_reduce(self, ser, op_name, skipna):
        pa_dtype = ser.dtype.pyarrow_dtype
        if op_name == "count":
            result = getattr(ser, op_name)()
        else:
            result = getattr(ser, op_name)(skipna=skipna)
        if pa.types.is_boolean(pa_dtype):
            # Can't convert if ser contains NA
            pytest.skip(
                "pandas boolean data with NA does not fully support all reductions"
            )
        elif pa.types.is_integer(pa_dtype) or pa.types.is_floating(pa_dtype):
            ser = ser.astype("Float64")
        if op_name == "count":
            expected = getattr(ser, op_name)()
        else:
            expected = getattr(ser, op_name)(skipna=skipna)
        tm.assert_almost_equal(result, expected)

    @pytest.mark.parametrize("skipna", [True, False])
    def test_reduce_series(self, data, all_numeric_reductions, skipna, request):
        pa_dtype = data.dtype.pyarrow_dtype
        opname = all_numeric_reductions

        ser = pd.Series(data)

        should_work = True
        if pa.types.is_temporal(pa_dtype) and opname in [
            "sum",
            "var",
            "skew",
            "kurt",
            "prod",
        ]:
            if pa.types.is_duration(pa_dtype) and opname in ["sum"]:
                # summing timedeltas is one case that *is* well-defined
                pass
            else:
                should_work = False
        elif (
            pa.types.is_string(pa_dtype) or pa.types.is_binary(pa_dtype)
        ) and opname in [
            "sum",
            "mean",
            "median",
            "prod",
            "std",
            "sem",
            "var",
            "skew",
            "kurt",
        ]:
            should_work = False

        if not should_work:
            # matching the non-pyarrow versions, these operations *should* not
            #  work for these dtypes
            msg = f"does not support reduction '{opname}'"
            with pytest.raises(TypeError, match=msg):
                getattr(ser, opname)(skipna=skipna)

            return

        xfail_mark = pytest.mark.xfail(
            raises=TypeError,
            reason=(
                f"{all_numeric_reductions} is not implemented in "
                f"pyarrow={pa.__version__} for {pa_dtype}"
            ),
        )
        if all_numeric_reductions in {"skew", "kurt"}:
            request.node.add_marker(xfail_mark)
        elif (
            all_numeric_reductions in {"median", "var", "std", "prod", "max", "min"}
            and pa_version_under6p0
        ):
            request.node.add_marker(xfail_mark)
        elif all_numeric_reductions == "sem" and pa_version_under8p0:
            request.node.add_marker(xfail_mark)
        elif (
            all_numeric_reductions in {"sum", "mean"}
            and skipna is False
            and pa_version_under6p0
            and (pa.types.is_integer(pa_dtype) or pa.types.is_floating(pa_dtype))
        ):
            request.node.add_marker(
                pytest.mark.xfail(
                    raises=AssertionError,
                    reason=(
                        f"{all_numeric_reductions} with skip_nulls={skipna} did not "
                        f"return NA for {pa_dtype} with pyarrow={pa.__version__}"
                    ),
                )
            )

        elif all_numeric_reductions in [
            "mean",
            "median",
            "std",
            "sem",
        ] and pa.types.is_temporal(pa_dtype):
            request.node.add_marker(xfail_mark)
        elif all_numeric_reductions in ["sum", "min", "max"] and pa.types.is_duration(
            pa_dtype
        ):
            request.node.add_marker(xfail_mark)
        elif pa.types.is_boolean(pa_dtype) and all_numeric_reductions in {
            "sem",
            "std",
            "var",
            "median",
        }:
            request.node.add_marker(xfail_mark)
        super().test_reduce_series(data, all_numeric_reductions, skipna)


class TestBaseBooleanReduce(base.BaseBooleanReduceTests):
    @pytest.mark.parametrize("skipna", [True, False])
    def test_reduce_series(
        self, data, all_boolean_reductions, skipna, na_value, request
    ):
        pa_dtype = data.dtype.pyarrow_dtype
        xfail_mark = pytest.mark.xfail(
            raises=TypeError,
            reason=(
                f"{all_boolean_reductions} is not implemented in "
                f"pyarrow={pa.__version__} for {pa_dtype}"
            ),
        )
        if not pa.types.is_boolean(pa_dtype):
            request.node.add_marker(xfail_mark)
        op_name = all_boolean_reductions
        ser = pd.Series(data)
        result = getattr(ser, op_name)(skipna=skipna)
        assert result is (op_name == "any")


class TestBaseGroupby(base.BaseGroupbyTests):
    def test_groupby_extension_no_sort(self, data_for_grouping, request):
        pa_dtype = data_for_grouping.dtype.pyarrow_dtype
        if pa.types.is_boolean(pa_dtype):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=f"{pa_dtype} only has 2 unique possible values",
                )
            )
        super().test_groupby_extension_no_sort(data_for_grouping)

    def test_groupby_extension_transform(self, data_for_grouping, request):
        pa_dtype = data_for_grouping.dtype.pyarrow_dtype
        if pa.types.is_boolean(pa_dtype):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=f"{pa_dtype} only has 2 unique possible values",
                )
            )
        with tm.maybe_produces_warning(
            PerformanceWarning,
            pa_version_under7p0 and not pa.types.is_duration(pa_dtype),
            check_stacklevel=False,
        ):
            super().test_groupby_extension_transform(data_for_grouping)

    def test_groupby_extension_apply(
        self, data_for_grouping, groupby_apply_op, request
    ):
        pa_dtype = data_for_grouping.dtype.pyarrow_dtype
        with tm.maybe_produces_warning(
            PerformanceWarning,
            pa_version_under7p0 and not pa.types.is_duration(pa_dtype),
            check_stacklevel=False,
        ):
            super().test_groupby_extension_apply(data_for_grouping, groupby_apply_op)

    @pytest.mark.parametrize("as_index", [True, False])
    def test_groupby_extension_agg(self, as_index, data_for_grouping, request):
        pa_dtype = data_for_grouping.dtype.pyarrow_dtype
        if pa.types.is_boolean(pa_dtype):
            request.node.add_marker(
                pytest.mark.xfail(
                    raises=ValueError,
                    reason=f"{pa_dtype} only has 2 unique possible values",
                )
            )
        with tm.maybe_produces_warning(
            PerformanceWarning,
            pa_version_under7p0 and not pa.types.is_duration(pa_dtype),
            check_stacklevel=False,
        ):
            super().test_groupby_extension_agg(as_index, data_for_grouping)


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
        elif pa.types.is_string(pa_dtype):
            request.node.add_marker(
                pytest.mark.xfail(
                    raises=TypeError,
                    reason=(
                        "Still support StringDtype('pyarrow') "
                        "over ArrowDtype(pa.string())"
                    ),
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
        elif pa.types.is_string(pa_dtype):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=(
                        "Still support StringDtype('pyarrow') "
                        "over ArrowDtype(pa.string())"
                    ),
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
        elif pa.types.is_string(pa_dtype):
            request.node.add_marker(
                pytest.mark.xfail(
                    raises=TypeError,
                    reason=(
                        "Still support StringDtype('pyarrow') "
                        "over ArrowDtype(pa.string())"
                    ),
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
            or pa.types.is_string(pa_dtype)
            or pa.types.is_binary(pa_dtype)
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
    @pytest.mark.xfail(reason="GH 45419: pyarrow.ChunkedArray does not support views.")
    def test_view(self, data):
        super().test_view(data)


class TestBaseMissing(base.BaseMissingTests):
    def test_dropna_array(self, data_missing):
        with tm.maybe_produces_warning(
            PerformanceWarning, pa_version_under6p0, check_stacklevel=False
        ):
            super().test_dropna_array(data_missing)

    def test_fillna_no_op_returns_copy(self, data):
        with tm.maybe_produces_warning(
            PerformanceWarning, pa_version_under7p0, check_stacklevel=False
        ):
            super().test_fillna_no_op_returns_copy(data)

    def test_fillna_series_method(self, data_missing, fillna_method):
        with tm.maybe_produces_warning(
            PerformanceWarning, pa_version_under7p0, check_stacklevel=False
        ):
            super().test_fillna_series_method(data_missing, fillna_method)


class TestBasePrinting(base.BasePrintingTests):
    pass


class TestBaseReshaping(base.BaseReshapingTests):
    @pytest.mark.xfail(reason="GH 45419: pyarrow.ChunkedArray does not support views")
    def test_transpose(self, data):
        super().test_transpose(data)


class TestBaseSetitem(base.BaseSetitemTests):
    @pytest.mark.xfail(reason="GH 45419: pyarrow.ChunkedArray does not support views")
    def test_setitem_preserves_views(self, data):
        super().test_setitem_preserves_views(data)


class TestBaseParsing(base.BaseParsingTests):
    @pytest.mark.parametrize("engine", ["c", "python"])
    def test_EA_types(self, engine, data, request):
        pa_dtype = data.dtype.pyarrow_dtype
        if pa.types.is_boolean(pa_dtype):
            request.node.add_marker(
                pytest.mark.xfail(raises=TypeError, reason="GH 47534")
            )
        elif pa.types.is_timestamp(pa_dtype) and pa_dtype.tz is not None:
            request.node.add_marker(
                pytest.mark.xfail(
                    raises=NotImplementedError,
                    reason=f"Parameterized types with tz={pa_dtype.tz} not supported.",
                )
            )
        elif pa.types.is_timestamp(pa_dtype) and pa_dtype.unit in ("us", "ns"):
            request.node.add_marker(
                pytest.mark.xfail(
                    raises=ValueError,
                    reason="https://github.com/pandas-dev/pandas/issues/49767",
                )
            )
        elif pa.types.is_binary(pa_dtype):
            request.node.add_marker(
                pytest.mark.xfail(reason="CSV parsers don't correctly handle binary")
            )
        df = pd.DataFrame({"with_dtype": pd.Series(data, dtype=str(data.dtype))})
        csv_output = df.to_csv(index=False, na_rep=np.nan)
        if pa.types.is_binary(pa_dtype):
            csv_output = BytesIO(csv_output)
        else:
            csv_output = StringIO(csv_output)
        result = pd.read_csv(
            csv_output, dtype={"with_dtype": str(data.dtype)}, engine=engine
        )
        expected = df
        self.assert_frame_equal(result, expected)


class TestBaseUnaryOps(base.BaseUnaryOpsTests):
    def test_invert(self, data, request):
        pa_dtype = data.dtype.pyarrow_dtype
        if not pa.types.is_boolean(pa_dtype):
            request.node.add_marker(
                pytest.mark.xfail(
                    raises=pa.ArrowNotImplementedError,
                    reason=f"pyarrow.compute.invert does support {pa_dtype}",
                )
            )
        super().test_invert(data)


class TestBaseMethods(base.BaseMethodsTests):
    def test_argsort_missing_array(self, data_missing_for_sorting):
        with tm.maybe_produces_warning(
            PerformanceWarning, pa_version_under7p0, check_stacklevel=False
        ):
            super().test_argsort_missing_array(data_missing_for_sorting)

    @pytest.mark.parametrize("periods", [1, -2])
    def test_diff(self, data, periods, request):
        pa_dtype = data.dtype.pyarrow_dtype
        if pa.types.is_unsigned_integer(pa_dtype) and periods == 1:
            request.node.add_marker(
                pytest.mark.xfail(
                    raises=pa.ArrowInvalid,
                    reason=(
                        f"diff with {pa_dtype} and periods={periods} will overflow"
                    ),
                )
            )
        super().test_diff(data, periods)

    @pytest.mark.filterwarnings("ignore:Falling back:pandas.errors.PerformanceWarning")
    @pytest.mark.parametrize("dropna", [True, False])
    def test_value_counts(self, all_data, dropna, request):
        super().test_value_counts(all_data, dropna)

    def test_value_counts_with_normalize(self, data, request):
        pa_dtype = data.dtype.pyarrow_dtype
        with tm.maybe_produces_warning(
            PerformanceWarning,
            pa_version_under7p0 and not pa.types.is_duration(pa_dtype),
            check_stacklevel=False,
        ):
            super().test_value_counts_with_normalize(data)

    @pytest.mark.xfail(
        pa_version_under6p0,
        raises=NotImplementedError,
        reason="argmin/max only implemented for pyarrow version >= 6.0",
    )
    def test_argmin_argmax(
        self, data_for_sorting, data_missing_for_sorting, na_value, request
    ):
        pa_dtype = data_for_sorting.dtype.pyarrow_dtype
        if pa.types.is_boolean(pa_dtype):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=f"{pa_dtype} only has 2 unique possible values",
                )
            )
        super().test_argmin_argmax(data_for_sorting, data_missing_for_sorting, na_value)

    @pytest.mark.parametrize(
        "op_name, skipna, expected",
        [
            ("idxmax", True, 0),
            ("idxmin", True, 2),
            ("argmax", True, 0),
            ("argmin", True, 2),
            ("idxmax", False, np.nan),
            ("idxmin", False, np.nan),
            ("argmax", False, -1),
            ("argmin", False, -1),
        ],
    )
    def test_argreduce_series(
        self, data_missing_for_sorting, op_name, skipna, expected, request
    ):
        if pa_version_under6p0 and skipna:
            request.node.add_marker(
                pytest.mark.xfail(
                    raises=NotImplementedError,
                    reason="min_max not supported in pyarrow",
                )
            )
        super().test_argreduce_series(
            data_missing_for_sorting, op_name, skipna, expected
        )

    @pytest.mark.parametrize(
        "na_position, expected",
        [
            ("last", np.array([2, 0, 1], dtype=np.dtype("intp"))),
            ("first", np.array([1, 2, 0], dtype=np.dtype("intp"))),
        ],
    )
    def test_nargsort(self, data_missing_for_sorting, na_position, expected):
        with tm.maybe_produces_warning(
            PerformanceWarning, pa_version_under7p0, check_stacklevel=False
        ):
            super().test_nargsort(data_missing_for_sorting, na_position, expected)

    @pytest.mark.parametrize("ascending", [True, False])
    def test_sort_values(self, data_for_sorting, ascending, sort_by_key, request):
        with tm.maybe_produces_warning(
            PerformanceWarning, pa_version_under7p0, check_stacklevel=False
        ):
            super().test_sort_values(data_for_sorting, ascending, sort_by_key)

    @pytest.mark.parametrize("ascending", [True, False])
    def test_sort_values_missing(
        self, data_missing_for_sorting, ascending, sort_by_key
    ):
        with tm.maybe_produces_warning(
            PerformanceWarning, pa_version_under7p0, check_stacklevel=False
        ):
            super().test_sort_values_missing(
                data_missing_for_sorting, ascending, sort_by_key
            )

    @pytest.mark.parametrize("ascending", [True, False])
    def test_sort_values_frame(self, data_for_sorting, ascending, request):
        pa_dtype = data_for_sorting.dtype.pyarrow_dtype
        with tm.maybe_produces_warning(
            PerformanceWarning,
            pa_version_under7p0 and not pa.types.is_duration(pa_dtype),
            check_stacklevel=False,
        ):
            super().test_sort_values_frame(data_for_sorting, ascending)

    def test_factorize(self, data_for_grouping, request):
        pa_dtype = data_for_grouping.dtype.pyarrow_dtype
        if pa.types.is_boolean(pa_dtype):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=f"{pa_dtype} only has 2 unique possible values",
                )
            )
        super().test_factorize(data_for_grouping)

    @pytest.mark.xfail(
        reason="result dtype pyarrow[bool] better than expected dtype object"
    )
    def test_combine_le(self, data_repeated):
        super().test_combine_le(data_repeated)

    def test_combine_add(self, data_repeated, request):
        pa_dtype = next(data_repeated(1)).dtype.pyarrow_dtype
        if pa.types.is_duration(pa_dtype):
            # TODO: this fails on the scalar addition constructing 'expected'
            #  but not in the actual 'combine' call, so may be salvage-able
            mark = pytest.mark.xfail(
                raises=TypeError,
                reason=f"{pa_dtype} cannot be added to {pa_dtype}",
            )
            request.node.add_marker(mark)
            super().test_combine_add(data_repeated)

        elif pa.types.is_temporal(pa_dtype):
            # analogous to datetime64, these cannot be added
            orig_data1, orig_data2 = data_repeated(2)
            s1 = pd.Series(orig_data1)
            s2 = pd.Series(orig_data2)
            with pytest.raises(TypeError):
                s1.combine(s2, lambda x1, x2: x1 + x2)

        else:
            super().test_combine_add(data_repeated)

    def test_searchsorted(self, data_for_sorting, as_series, request):
        pa_dtype = data_for_sorting.dtype.pyarrow_dtype
        if pa.types.is_boolean(pa_dtype):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=f"{pa_dtype} only has 2 unique possible values",
                )
            )
        super().test_searchsorted(data_for_sorting, as_series)

    def test_basic_equals(self, data):
        # https://github.com/pandas-dev/pandas/issues/34660
        assert pd.Series(data).equals(pd.Series(data))


class TestBaseArithmeticOps(base.BaseArithmeticOpsTests):

    divmod_exc = NotImplementedError

    def _patch_combine(self, obj, other, op):
        # BaseOpsUtil._combine can upcast expected dtype
        # (because it generates expected on python scalars)
        # while ArrowExtensionArray maintains original type
        expected = base.BaseArithmeticOpsTests._combine(self, obj, other, op)
        was_frame = False
        if isinstance(expected, pd.DataFrame):
            was_frame = True
            expected_data = expected.iloc[:, 0]
            original_dtype = obj.iloc[:, 0].dtype
        else:
            expected_data = expected
            original_dtype = obj.dtype
        pa_array = pa.array(expected_data._values).cast(original_dtype.pyarrow_dtype)
        pd_array = type(expected_data._values)(pa_array)
        if was_frame:
            expected = pd.DataFrame(
                pd_array, index=expected.index, columns=expected.columns
            )
        else:
            expected = pd.Series(pd_array)
        return expected

    def _is_temporal_supported(self, opname, pa_dtype):
        return not pa_version_under8p0 and (
            opname in ("__add__", "__radd__")
            and pa.types.is_duration(pa_dtype)
            or opname in ("__sub__", "__rsub__")
            and pa.types.is_temporal(pa_dtype)
        )

    def _get_scalar_exception(self, opname, pa_dtype):
        arrow_temporal_supported = self._is_temporal_supported(opname, pa_dtype)
        if opname in {
            "__mod__",
            "__rmod__",
        }:
            exc = NotImplementedError
        elif arrow_temporal_supported:
            exc = None
        elif not (pa.types.is_floating(pa_dtype) or pa.types.is_integer(pa_dtype)):
            exc = pa.ArrowNotImplementedError
        else:
            exc = None
        return exc

    def _get_arith_xfail_marker(self, opname, pa_dtype):
        mark = None

        arrow_temporal_supported = self._is_temporal_supported(opname, pa_dtype)

        if (
            opname == "__rpow__"
            and (pa.types.is_floating(pa_dtype) or pa.types.is_integer(pa_dtype))
            and not pa_version_under6p0
        ):
            mark = pytest.mark.xfail(
                reason=(
                    f"GH#29997: 1**pandas.NA == 1 while 1**pyarrow.NA == NULL "
                    f"for {pa_dtype}"
                )
            )
        elif arrow_temporal_supported:
            mark = pytest.mark.xfail(
                raises=TypeError,
                reason=(
                    f"{opname} not supported between"
                    f"pd.NA and {pa_dtype} Python scalar"
                ),
            )
        elif (
            opname in {"__rtruediv__", "__rfloordiv__"}
            and (pa.types.is_floating(pa_dtype) or pa.types.is_integer(pa_dtype))
            and not pa_version_under6p0
        ):
            mark = pytest.mark.xfail(
                raises=pa.ArrowInvalid,
                reason="divide by 0",
            )

        return mark

    def test_arith_series_with_scalar(
        self, data, all_arithmetic_operators, request, monkeypatch
    ):
        pa_dtype = data.dtype.pyarrow_dtype

        if all_arithmetic_operators == "__rmod__" and (
            pa.types.is_string(pa_dtype) or pa.types.is_binary(pa_dtype)
        ):
            pytest.skip("Skip testing Python string formatting")

        self.series_scalar_exc = self._get_scalar_exception(
            all_arithmetic_operators, pa_dtype
        )

        mark = self._get_arith_xfail_marker(all_arithmetic_operators, pa_dtype)
        if mark is not None:
            request.node.add_marker(mark)

        if all_arithmetic_operators == "__floordiv__" and pa.types.is_integer(pa_dtype):
            # BaseOpsUtil._combine always returns int64, while ArrowExtensionArray does
            # not upcast
            monkeypatch.setattr(TestBaseArithmeticOps, "_combine", self._patch_combine)
        super().test_arith_series_with_scalar(data, all_arithmetic_operators)

    def test_arith_frame_with_scalar(
        self, data, all_arithmetic_operators, request, monkeypatch
    ):
        pa_dtype = data.dtype.pyarrow_dtype

        if all_arithmetic_operators == "__rmod__" and (
            pa.types.is_string(pa_dtype) or pa.types.is_binary(pa_dtype)
        ):
            pytest.skip("Skip testing Python string formatting")

        self.frame_scalar_exc = self._get_scalar_exception(
            all_arithmetic_operators, pa_dtype
        )

        mark = self._get_arith_xfail_marker(all_arithmetic_operators, pa_dtype)
        if mark is not None:
            request.node.add_marker(mark)

        if all_arithmetic_operators == "__floordiv__" and pa.types.is_integer(pa_dtype):
            # BaseOpsUtil._combine always returns int64, while ArrowExtensionArray does
            # not upcast
            monkeypatch.setattr(TestBaseArithmeticOps, "_combine", self._patch_combine)
        super().test_arith_frame_with_scalar(data, all_arithmetic_operators)

    def test_arith_series_with_array(
        self, data, all_arithmetic_operators, request, monkeypatch
    ):
        pa_dtype = data.dtype.pyarrow_dtype

        self.series_array_exc = self._get_scalar_exception(
            all_arithmetic_operators, pa_dtype
        )

        if (
            all_arithmetic_operators
            in (
                "__sub__",
                "__rsub__",
            )
            and pa.types.is_unsigned_integer(pa_dtype)
            and not pa_version_under6p0
        ):
            request.node.add_marker(
                pytest.mark.xfail(
                    raises=pa.ArrowInvalid,
                    reason=(
                        f"Implemented pyarrow.compute.subtract_checked "
                        f"which raises on overflow for {pa_dtype}"
                    ),
                )
            )

        mark = self._get_arith_xfail_marker(all_arithmetic_operators, pa_dtype)
        if mark is not None:
            request.node.add_marker(mark)

        op_name = all_arithmetic_operators
        ser = pd.Series(data)
        # pd.Series([ser.iloc[0]] * len(ser)) may not return ArrowExtensionArray
        # since ser.iloc[0] is a python scalar
        other = pd.Series(pd.array([ser.iloc[0]] * len(ser), dtype=data.dtype))

        if pa.types.is_floating(pa_dtype) or (
            pa.types.is_integer(pa_dtype) and all_arithmetic_operators != "__truediv__"
        ):
            monkeypatch.setattr(TestBaseArithmeticOps, "_combine", self._patch_combine)
        self.check_opname(ser, op_name, other, exc=self.series_array_exc)

    def test_add_series_with_extension_array(self, data, request):
        pa_dtype = data.dtype.pyarrow_dtype
        if not (
            pa.types.is_integer(pa_dtype)
            or pa.types.is_floating(pa_dtype)
            or (not pa_version_under8p0 and pa.types.is_duration(pa_dtype))
        ):
            request.node.add_marker(
                pytest.mark.xfail(
                    raises=NotImplementedError,
                    reason=f"add_checked not implemented for {pa_dtype}",
                )
            )
        elif pa_dtype.equals("int8"):
            request.node.add_marker(
                pytest.mark.xfail(
                    raises=pa.ArrowInvalid,
                    reason=f"raises on overflow for {pa_dtype}",
                )
            )
        super().test_add_series_with_extension_array(data)


class TestBaseComparisonOps(base.BaseComparisonOpsTests):
    def assert_series_equal(self, left, right, *args, **kwargs):
        # Series.combine for "expected" retains bool[pyarrow] dtype
        # While "result" return "boolean" dtype
        right = pd.Series(right._values.to_numpy(), dtype="boolean")
        super().assert_series_equal(left, right, *args, **kwargs)

    def test_compare_array(self, data, comparison_op, na_value, request):
        pa_dtype = data.dtype.pyarrow_dtype
        ser = pd.Series(data)
        # pd.Series([ser.iloc[0]] * len(ser)) may not return ArrowExtensionArray
        # since ser.iloc[0] is a python scalar
        other = pd.Series(pd.array([ser.iloc[0]] * len(ser), dtype=data.dtype))
        if comparison_op.__name__ in ["eq", "ne"]:
            # comparison should match point-wise comparisons
            result = comparison_op(ser, other)
            # Series.combine does not calculate the NA mask correctly
            # when comparing over an array
            assert result[8] is na_value
            assert result[97] is na_value
            expected = ser.combine(other, comparison_op)
            expected[8] = na_value
            expected[97] = na_value
            self.assert_series_equal(result, expected)

        else:
            exc = None
            try:
                result = comparison_op(ser, other)
            except Exception as err:
                exc = err

            if exc is None:
                # Didn't error, then should match point-wise behavior
                if pa.types.is_temporal(pa_dtype):
                    # point-wise comparison with pd.NA raises TypeError
                    assert result[8] is na_value
                    assert result[97] is na_value
                    result = result.drop([8, 97]).reset_index(drop=True)
                    ser = ser.drop([8, 97])
                    other = other.drop([8, 97])
                expected = ser.combine(other, comparison_op)
                self.assert_series_equal(result, expected)
            else:
                with pytest.raises(type(exc)):
                    ser.combine(other, comparison_op)

    def test_invalid_other_comp(self, data, comparison_op):
        # GH 48833
        with pytest.raises(
            NotImplementedError, match=".* not implemented for <class 'object'>"
        ):
            comparison_op(data, object())


def test_arrowdtype_construct_from_string_type_with_unsupported_parameters():
    with pytest.raises(NotImplementedError, match="Passing pyarrow type"):
        ArrowDtype.construct_from_string("timestamp[s, tz=UTC][pyarrow]")


@pytest.mark.parametrize(
    "interpolation", ["linear", "lower", "higher", "nearest", "midpoint"]
)
@pytest.mark.parametrize("quantile", [0.5, [0.5, 0.5]])
def test_quantile(data, interpolation, quantile, request):
    pa_dtype = data.dtype.pyarrow_dtype

    data = data.take([0, 0, 0])
    ser = pd.Series(data)

    if (
        pa.types.is_string(pa_dtype)
        or pa.types.is_binary(pa_dtype)
        or pa.types.is_boolean(pa_dtype)
    ):
        # For string, bytes, and bool, we don't *expect* to have quantile work
        # Note this matches the non-pyarrow behavior
        if pa_version_under7p0:
            msg = r"Function quantile has no kernel matching input types \(.*\)"
        else:
            msg = r"Function 'quantile' has no kernel matching input types \(.*\)"
        with pytest.raises(pa.ArrowNotImplementedError, match=msg):
            ser.quantile(q=quantile, interpolation=interpolation)
        return

    if pa.types.is_integer(pa_dtype) or pa.types.is_floating(pa_dtype):
        pass
    elif pa.types.is_temporal(data._data.type) and interpolation in ["lower", "higher"]:
        pass
    else:
        request.node.add_marker(
            pytest.mark.xfail(
                raises=pa.ArrowNotImplementedError,
                reason=f"quantile not supported by pyarrow for {pa_dtype}",
            )
        )
    data = data.take([0, 0, 0])
    ser = pd.Series(data)
    result = ser.quantile(q=quantile, interpolation=interpolation)
    if quantile == 0.5:
        assert result == data[0]
    else:
        # Just check the values
        expected = pd.Series(data.take([0, 0]), index=[0.5, 0.5])
        if pa.types.is_integer(pa_dtype) or pa.types.is_floating(pa_dtype):
            expected = expected.astype("float64[pyarrow]")
            result = result.astype("float64[pyarrow]")
        tm.assert_series_equal(result, expected)


@pytest.mark.xfail(
    pa_version_under6p0,
    raises=NotImplementedError,
    reason="mode only supported for pyarrow version >= 6.0",
)
@pytest.mark.parametrize("dropna", [True, False])
@pytest.mark.parametrize(
    "take_idx, exp_idx",
    [[[0, 0, 2, 2, 4, 4], [4, 0]], [[0, 0, 0, 2, 4, 4], [0]]],
    ids=["multi_mode", "single_mode"],
)
def test_mode(data_for_grouping, dropna, take_idx, exp_idx, request):
    pa_dtype = data_for_grouping.dtype.pyarrow_dtype
    if pa.types.is_string(pa_dtype) or pa.types.is_binary(pa_dtype):
        request.node.add_marker(
            pytest.mark.xfail(
                raises=pa.ArrowNotImplementedError,
                reason=f"mode not supported by pyarrow for {pa_dtype}",
            )
        )
    elif (
        pa.types.is_boolean(pa_dtype)
        and "multi_mode" in request.node.nodeid
        and pa_version_under9p0
    ):
        request.node.add_marker(
            pytest.mark.xfail(
                reason="https://issues.apache.org/jira/browse/ARROW-17096",
            )
        )
    data = data_for_grouping.take(take_idx)
    ser = pd.Series(data)
    result = ser.mode(dropna=dropna)
    expected = pd.Series(data_for_grouping.take(exp_idx))
    tm.assert_series_equal(result, expected)


def test_is_bool_dtype():
    # GH 22667
    data = ArrowExtensionArray(pa.array([True, False, True]))
    assert is_bool_dtype(data)
    assert pd.core.common.is_bool_indexer(data)
    s = pd.Series(range(len(data)))
    result = s[data]
    expected = s[np.asarray(data)]
    tm.assert_series_equal(result, expected)


def test_is_numeric_dtype(data):
    # GH 50563
    pa_type = data.dtype.pyarrow_dtype
    if (
        pa.types.is_floating(pa_type)
        or pa.types.is_integer(pa_type)
        or pa.types.is_decimal(pa_type)
    ):
        assert is_numeric_dtype(data)
    else:
        assert not is_numeric_dtype(data)


def test_is_integer_dtype(data):
    # GH 50667
    pa_type = data.dtype.pyarrow_dtype
    if pa.types.is_integer(pa_type):
        assert is_integer_dtype(data)
    else:
        assert not is_integer_dtype(data)


def test_is_any_integer_dtype(data):
    # GH 50667
    pa_type = data.dtype.pyarrow_dtype
    if pa.types.is_integer(pa_type):
        assert is_any_int_dtype(data)
    else:
        assert not is_any_int_dtype(data)


def test_is_signed_integer_dtype(data):
    pa_type = data.dtype.pyarrow_dtype
    if pa.types.is_signed_integer(pa_type):
        assert is_signed_integer_dtype(data)
    else:
        assert not is_signed_integer_dtype(data)


def test_is_unsigned_integer_dtype(data):
    pa_type = data.dtype.pyarrow_dtype
    if pa.types.is_unsigned_integer(pa_type):
        assert is_unsigned_integer_dtype(data)
    else:
        assert not is_unsigned_integer_dtype(data)


def test_is_float_dtype(data):
    pa_type = data.dtype.pyarrow_dtype
    if pa.types.is_floating(pa_type):
        assert is_float_dtype(data)
    else:
        assert not is_float_dtype(data)


def test_pickle_roundtrip(data):
    # GH 42600
    expected = pd.Series(data)
    expected_sliced = expected.head(2)
    full_pickled = pickle.dumps(expected)
    sliced_pickled = pickle.dumps(expected_sliced)

    assert len(full_pickled) > len(sliced_pickled)

    result = pickle.loads(full_pickled)
    tm.assert_series_equal(result, expected)

    result_sliced = pickle.loads(sliced_pickled)
    tm.assert_series_equal(result_sliced, expected_sliced)


def test_astype_from_non_pyarrow(data):
    # GH49795
    pd_array = data._data.to_pandas().array
    result = pd_array.astype(data.dtype)
    assert not isinstance(pd_array.dtype, ArrowDtype)
    assert isinstance(result.dtype, ArrowDtype)
    tm.assert_extension_array_equal(result, data)


def test_astype_float_from_non_pyarrow_str():
    # GH50430
    ser = pd.Series(["1.0"])
    result = ser.astype("float64[pyarrow]")
    expected = pd.Series([1.0], dtype="float64[pyarrow]")
    tm.assert_series_equal(result, expected)


def test_to_numpy_with_defaults(data):
    # GH49973
    result = data.to_numpy()

    pa_type = data._data.type
    if pa.types.is_duration(pa_type) or pa.types.is_timestamp(pa_type):
        expected = np.array(list(data))
    else:
        expected = np.array(data._data)

    if data._hasna:
        expected = expected.astype(object)
        expected[pd.isna(data)] = pd.NA

    tm.assert_numpy_array_equal(result, expected)


def test_setitem_null_slice(data):
    # GH50248
    orig = data.copy()

    result = orig.copy()
    result[:] = data[0]
    expected = ArrowExtensionArray(
        pa.array([data[0]] * len(data), type=data._data.type)
    )
    tm.assert_extension_array_equal(result, expected)

    result = orig.copy()
    result[:] = data[::-1]
    expected = data[::-1]
    tm.assert_extension_array_equal(result, expected)

    result = orig.copy()
    result[:] = data.tolist()
    expected = data
    tm.assert_extension_array_equal(result, expected)


def test_setitem_invalid_dtype(data):
    # GH50248
    pa_type = data._data.type
    if pa.types.is_string(pa_type) or pa.types.is_binary(pa_type):
        fill_value = 123
        err = TypeError
        msg = "Invalid value '123' for dtype"
    elif (
        pa.types.is_integer(pa_type)
        or pa.types.is_floating(pa_type)
        or pa.types.is_boolean(pa_type)
    ):
        fill_value = "foo"
        err = pa.ArrowInvalid
        msg = "Could not convert"
    else:
        fill_value = "foo"
        err = TypeError
        msg = "Invalid value 'foo' for dtype"
    with pytest.raises(err, match=msg):
        data[:] = fill_value


def test_round():
    dtype = "float64[pyarrow]"

    ser = pd.Series([0.0, 1.23, 2.56, pd.NA], dtype=dtype)
    result = ser.round(1)
    expected = pd.Series([0.0, 1.2, 2.6, pd.NA], dtype=dtype)
    tm.assert_series_equal(result, expected)

    ser = pd.Series([123.4, pd.NA, 56.78], dtype=dtype)
    result = ser.round(-1)
    expected = pd.Series([120.0, pd.NA, 60.0], dtype=dtype)
    tm.assert_series_equal(result, expected)


def test_searchsorted_with_na_raises(data_for_sorting, as_series):
    # GH50447
    b, c, a = data_for_sorting
    arr = data_for_sorting.take([2, 0, 1])  # to get [a, b, c]
    arr[-1] = pd.NA

    if as_series:
        arr = pd.Series(arr)

    msg = (
        "searchsorted requires array to be sorted, "
        "which is impossible with NAs present."
    )
    with pytest.raises(ValueError, match=msg):
        arr.searchsorted(b)
