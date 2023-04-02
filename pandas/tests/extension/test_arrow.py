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
from decimal import Decimal
from io import (
    BytesIO,
    StringIO,
)
import pickle
import re

import numpy as np
import pytest

from pandas.compat import (
    PY311,
    is_ci_environment,
    is_platform_windows,
    pa_version_under7p0,
    pa_version_under8p0,
    pa_version_under9p0,
    pa_version_under11p0,
)
from pandas.errors import PerformanceWarning

from pandas.core.dtypes.common import is_any_int_dtype
from pandas.core.dtypes.dtypes import CategoricalDtypeType

import pandas as pd
import pandas._testing as tm
from pandas.api.types import (
    is_bool_dtype,
    is_float_dtype,
    is_integer_dtype,
    is_numeric_dtype,
    is_signed_integer_dtype,
    is_string_dtype,
    is_unsigned_integer_dtype,
)
from pandas.tests.extension import base

pa = pytest.importorskip("pyarrow", minversion="7.0.0")

from pandas.core.arrays.arrow.array import ArrowExtensionArray

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
    elif pa.types.is_decimal(pa_dtype):
        data = (
            [Decimal("1"), Decimal("0.0")] * 4
            + [None]
            + [Decimal("-2.0"), Decimal("-1.0")] * 44
            + [None]
            + [Decimal("0.5"), Decimal("33.123")]
        )
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
    elif pa.types.is_decimal(pa_dtype):
        A = Decimal("-1.1")
        B = Decimal("0.0")
        C = Decimal("1.1")
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
        if pa.types.is_string(pa_dtype) or pa.types.is_decimal(pa_dtype):
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

    def test_from_sequence_pa_array(self, data):
        # https://github.com/pandas-dev/pandas/pull/47034#discussion_r955500784
        # data._data = pa.ChunkedArray
        result = type(data)._from_sequence(data._data)
        tm.assert_extension_array_equal(result, data)
        assert isinstance(result._data, pa.ChunkedArray)

        result = type(data)._from_sequence(data._data.combine_chunks())
        tm.assert_extension_array_equal(result, data)
        assert isinstance(result._data, pa.ChunkedArray)

    def test_from_sequence_pa_array_notimplemented(self, request):
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
        elif pa_version_under11p0 and (
            pa.types.is_duration(pa_dtype) or pa.types.is_decimal(pa_dtype)
        ):
            request.node.add_marker(
                pytest.mark.xfail(
                    raises=pa.ArrowNotImplementedError,
                    reason=f"pyarrow doesn't support parsing {pa_dtype}",
                )
            )
        elif pa.types.is_timestamp(pa_dtype) and pa_dtype.tz is not None:
            if is_platform_windows() and is_ci_environment():
                request.node.add_marker(
                    pytest.mark.xfail(
                        raises=pa.ArrowInvalid,
                        reason=(
                            "TODO: Set ARROW_TIMEZONE_DATABASE environment variable "
                            "on CI to path to the tzdata for pyarrow."
                        ),
                    )
                )
        pa_array = data._data.cast(pa.string())
        result = type(data)._from_sequence_of_strings(pa_array, dtype=data.dtype)
        tm.assert_extension_array_equal(result, data)

        pa_array = pa_array.combine_chunks()
        result = type(data)._from_sequence_of_strings(pa_array, dtype=data.dtype)
        tm.assert_extension_array_equal(result, data)


class TestGetitemTests(base.BaseGetitemTests):
    pass


class TestBaseAccumulateTests(base.BaseAccumulateTests):
    def check_accumulate(self, ser, op_name, skipna):
        result = getattr(ser, op_name)(skipna=skipna)

        if ser.dtype.kind == "m":
            # Just check that we match the integer behavior.
            ser = ser.astype("int64[pyarrow]")
            result = result.astype("int64[pyarrow]")

        result = result.astype("Float64")
        expected = getattr(ser.astype("Float64"), op_name)(skipna=skipna)
        self.assert_series_equal(result, expected, check_dtype=False)

    @pytest.mark.parametrize("skipna", [True, False])
    def test_accumulate_series_raises(self, data, all_numeric_accumulations, skipna):
        pa_type = data.dtype.pyarrow_dtype
        if (
            (
                pa.types.is_integer(pa_type)
                or pa.types.is_floating(pa_type)
                or pa.types.is_duration(pa_type)
            )
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
            if request.config.option.skip_slow:
                # equivalent to marking these cases with @pytest.mark.slow,
                #  these xfails take a long time to run because pytest
                #  renders the exception messages even when not showing them
                pytest.skip("pyarrow xfail slow")

            request.node.add_marker(
                pytest.mark.xfail(
                    reason=f"{all_numeric_accumulations} not implemented",
                    raises=NotImplementedError,
                )
            )
        elif all_numeric_accumulations == "cumsum" and (
            pa.types.is_boolean(pa_type) or pa.types.is_decimal(pa_type)
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
            all_numeric_reductions in {"var", "std", "median"}
            and pa_version_under7p0
            and pa.types.is_decimal(pa_dtype)
        ):
            request.node.add_marker(xfail_mark)
        elif all_numeric_reductions == "sem" and pa_version_under8p0:
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
        if pa.types.is_string(pa_dtype) or pa.types.is_binary(pa_dtype):
            # We *might* want to make this behave like the non-pyarrow cases,
            #  but have not yet decided.
            request.node.add_marker(xfail_mark)

        op_name = all_boolean_reductions
        ser = pd.Series(data)

        if pa.types.is_temporal(pa_dtype) and not pa.types.is_duration(pa_dtype):
            # xref GH#34479 we support this in our non-pyarrow datetime64 dtypes,
            #  but it isn't obvious we _should_.  For now, we keep the pyarrow
            #  behavior which does not support this.

            with pytest.raises(TypeError, match="does not support reduction"):
                getattr(ser, op_name)(skipna=skipna)

            return

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
        super().test_groupby_extension_transform(data_for_grouping)

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
        super().test_groupby_extension_agg(as_index, data_for_grouping)

    def test_in_numeric_groupby(self, data_for_grouping):
        if is_string_dtype(data_for_grouping.dtype):
            df = pd.DataFrame(
                {
                    "A": [1, 1, 2, 2, 3, 3, 1, 4],
                    "B": data_for_grouping,
                    "C": [1, 1, 1, 1, 1, 1, 1, 1],
                }
            )

            expected = pd.Index(["C"])
            with pytest.raises(TypeError, match="does not support"):
                df.groupby("A").sum().columns
            result = df.groupby("A").sum(numeric_only=True).columns
            tm.assert_index_equal(result, expected)
        else:
            super().test_in_numeric_groupby(data_for_grouping)


class TestBaseDtype(base.BaseDtypeTests):
    def test_check_dtype(self, data, request):
        pa_dtype = data.dtype.pyarrow_dtype
        if pa.types.is_decimal(pa_dtype) and pa_version_under8p0:
            request.node.add_marker(
                pytest.mark.xfail(
                    raises=ValueError,
                    reason="decimal string repr affects numpy comparison",
                )
            )
        super().test_check_dtype(data)

    def test_construct_from_string_own_name(self, dtype, request):
        pa_dtype = dtype.pyarrow_dtype
        if pa.types.is_decimal(pa_dtype):
            request.node.add_marker(
                pytest.mark.xfail(
                    raises=NotImplementedError,
                    reason=f"pyarrow.type_for_alias cannot infer {pa_dtype}",
                )
            )

        if pa.types.is_string(pa_dtype):
            # We still support StringDtype('pyarrow') over ArrowDtype(pa.string())
            msg = r"string\[pyarrow\] should be constructed by StringDtype"
            with pytest.raises(TypeError, match=msg):
                dtype.construct_from_string(dtype.name)

            return

        super().test_construct_from_string_own_name(dtype)

    def test_is_dtype_from_name(self, dtype, request):
        pa_dtype = dtype.pyarrow_dtype
        if pa.types.is_string(pa_dtype):
            # We still support StringDtype('pyarrow') over ArrowDtype(pa.string())
            assert not type(dtype).is_dtype(dtype.name)
        else:
            if pa.types.is_decimal(pa_dtype):
                request.node.add_marker(
                    pytest.mark.xfail(
                        raises=NotImplementedError,
                        reason=f"pyarrow.type_for_alias cannot infer {pa_dtype}",
                    )
                )
            super().test_is_dtype_from_name(dtype)

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
            or pa.types.is_binary(pa_dtype)
            or pa.types.is_decimal(pa_dtype)
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

    def test_is_not_string_type(self, dtype):
        pa_dtype = dtype.pyarrow_dtype
        if pa.types.is_string(pa_dtype):
            assert is_string_dtype(dtype)
        else:
            super().test_is_not_string_type(dtype)


class TestBaseIndex(base.BaseIndexTests):
    pass


class TestBaseInterface(base.BaseInterfaceTests):
    @pytest.mark.xfail(
        reason="GH 45419: pyarrow.ChunkedArray does not support views.", run=False
    )
    def test_view(self, data):
        super().test_view(data)


class TestBaseMissing(base.BaseMissingTests):
    def test_fillna_no_op_returns_copy(self, data):
        data = data[~data.isna()]

        valid = data[0]
        result = data.fillna(valid)
        assert result is not data
        self.assert_extension_array_equal(result, data)
        with tm.assert_produces_warning(PerformanceWarning):
            result = data.fillna(method="backfill")
        assert result is not data
        self.assert_extension_array_equal(result, data)

    def test_fillna_series_method(self, data_missing, fillna_method):
        with tm.maybe_produces_warning(
            PerformanceWarning, fillna_method is not None, check_stacklevel=False
        ):
            super().test_fillna_series_method(data_missing, fillna_method)


class TestBasePrinting(base.BasePrintingTests):
    pass


class TestBaseReshaping(base.BaseReshapingTests):
    @pytest.mark.xfail(
        reason="GH 45419: pyarrow.ChunkedArray does not support views", run=False
    )
    def test_transpose(self, data):
        super().test_transpose(data)


class TestBaseSetitem(base.BaseSetitemTests):
    @pytest.mark.xfail(
        reason="GH 45419: pyarrow.ChunkedArray does not support views", run=False
    )
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
        elif pa.types.is_decimal(pa_dtype):
            request.node.add_marker(
                pytest.mark.xfail(
                    raises=NotImplementedError,
                    reason=f"Parameterized types {pa_dtype} not supported.",
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

    def test_value_counts_returns_pyarrow_int64(self, data):
        # GH 51462
        data = data[:10]
        result = data.value_counts()
        assert result.dtype == ArrowDtype(pa.int64())

    def test_value_counts_with_normalize(self, data, request):
        data = data[:10].unique()
        values = np.array(data[~data.isna()])
        ser = pd.Series(data, dtype=data.dtype)

        result = ser.value_counts(normalize=True).sort_index()

        expected = pd.Series(
            [1 / len(values)] * len(values), index=result.index, name="proportion"
        )
        expected = expected.astype("double[pyarrow]")

        self.assert_series_equal(result, expected)

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
        elif pa.types.is_decimal(pa_dtype) and pa_version_under7p0:
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=f"No pyarrow kernel for {pa_dtype}",
                    raises=pa.ArrowNotImplementedError,
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
        pa_dtype = data_missing_for_sorting.dtype.pyarrow_dtype
        if pa.types.is_decimal(pa_dtype) and pa_version_under7p0 and skipna:
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=f"No pyarrow kernel for {pa_dtype}",
                    raises=pa.ArrowNotImplementedError,
                )
            )
        super().test_argreduce_series(
            data_missing_for_sorting, op_name, skipna, expected
        )

    def test_factorize(self, data_for_grouping, request):
        pa_dtype = data_for_grouping.dtype.pyarrow_dtype
        if pa.types.is_boolean(pa_dtype):
            request.node.add_marker(
                pytest.mark.xfail(
                    reason=f"{pa_dtype} only has 2 unique possible values",
                )
            )
        super().test_factorize(data_for_grouping)

    _combine_le_expected_dtype = "bool[pyarrow]"

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

    @classmethod
    def assert_equal(cls, left, right, **kwargs):
        if isinstance(left, pd.DataFrame):
            left_pa_type = left.iloc[:, 0].dtype.pyarrow_dtype
            right_pa_type = right.iloc[:, 0].dtype.pyarrow_dtype
        else:
            left_pa_type = left.dtype.pyarrow_dtype
            right_pa_type = right.dtype.pyarrow_dtype
        if pa.types.is_decimal(left_pa_type) or pa.types.is_decimal(right_pa_type):
            # decimal precision can resize in the result type depending on data
            # just compare the float values
            left = left.astype("float[pyarrow]")
            right = right.astype("float[pyarrow]")
        tm.assert_equal(left, right, **kwargs)

    def get_op_from_name(self, op_name):
        short_opname = op_name.strip("_")
        if short_opname == "rtruediv":
            # use the numpy version that won't raise on division by zero
            return lambda x, y: np.divide(y, x)
        elif short_opname == "rfloordiv":
            return lambda x, y: np.floor_divide(y, x)

        return tm.get_op_from_name(op_name)

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

        pa_expected = pa.array(expected_data._values)

        if pa.types.is_duration(pa_expected.type):
            # pyarrow sees sequence of datetime/timedelta objects and defaults
            #  to "us" but the non-pointwise op retains unit
            unit = original_dtype.pyarrow_dtype.unit
            if type(other) in [datetime, timedelta] and unit in ["s", "ms"]:
                # pydatetime/pytimedelta objects have microsecond reso, so we
                #  take the higher reso of the original and microsecond. Note
                #  this matches what we would do with DatetimeArray/TimedeltaArray
                unit = "us"
            pa_expected = pa_expected.cast(f"duration[{unit}]")
        else:
            pa_expected = pa_expected.cast(original_dtype.pyarrow_dtype)

        pd_expected = type(expected_data._values)(pa_expected)
        if was_frame:
            expected = pd.DataFrame(
                pd_expected, index=expected.index, columns=expected.columns
            )
        else:
            expected = pd.Series(pd_expected)
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
        elif opname in ["__add__", "__radd__"] and (
            pa.types.is_string(pa_dtype) or pa.types.is_binary(pa_dtype)
        ):
            exc = None
        elif not (
            pa.types.is_floating(pa_dtype)
            or pa.types.is_integer(pa_dtype)
            or pa.types.is_decimal(pa_dtype)
        ):
            exc = pa.ArrowNotImplementedError
        else:
            exc = None
        return exc

    def _get_arith_xfail_marker(self, opname, pa_dtype):
        mark = None

        arrow_temporal_supported = self._is_temporal_supported(opname, pa_dtype)

        if (
            opname == "__rpow__"
            and (
                pa.types.is_floating(pa_dtype)
                or pa.types.is_integer(pa_dtype)
                or pa.types.is_decimal(pa_dtype)
            )
            and not pa_version_under7p0
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
            opname == "__rfloordiv__"
            and (pa.types.is_integer(pa_dtype) or pa.types.is_decimal(pa_dtype))
            and not pa_version_under7p0
        ):
            mark = pytest.mark.xfail(
                raises=pa.ArrowInvalid,
                reason="divide by 0",
            )
        elif (
            opname == "__rtruediv__"
            and pa.types.is_decimal(pa_dtype)
            and not pa_version_under7p0
        ):
            mark = pytest.mark.xfail(
                raises=pa.ArrowInvalid,
                reason="divide by 0",
            )
        elif (
            opname == "__pow__"
            and pa.types.is_decimal(pa_dtype)
            and pa_version_under7p0
        ):
            mark = pytest.mark.xfail(
                raises=pa.ArrowInvalid,
                reason="Invalid decimal function: power_checked",
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

        if (
            (
                all_arithmetic_operators == "__floordiv__"
                and pa.types.is_integer(pa_dtype)
            )
            or pa.types.is_duration(pa_dtype)
            or pa.types.is_timestamp(pa_dtype)
        ):
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

        if (
            (
                all_arithmetic_operators == "__floordiv__"
                and pa.types.is_integer(pa_dtype)
            )
            or pa.types.is_duration(pa_dtype)
            or pa.types.is_timestamp(pa_dtype)
        ):
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
            and not pa_version_under7p0
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

        if (
            pa.types.is_floating(pa_dtype)
            or (
                pa.types.is_integer(pa_dtype)
                and all_arithmetic_operators not in ["__truediv__", "__rtruediv__"]
            )
            or pa.types.is_duration(pa_dtype)
            or pa.types.is_timestamp(pa_dtype)
        ):
            monkeypatch.setattr(TestBaseArithmeticOps, "_combine", self._patch_combine)
        self.check_opname(ser, op_name, other, exc=self.series_array_exc)

    def test_add_series_with_extension_array(self, data, request):
        pa_dtype = data.dtype.pyarrow_dtype

        if pa.types.is_temporal(pa_dtype) and not pa.types.is_duration(pa_dtype):
            # i.e. timestamp, date, time, but not timedelta; these *should*
            #  raise when trying to add
            ser = pd.Series(data)
            if pa_version_under7p0:
                msg = "Function add_checked has no kernel matching input types"
            else:
                msg = "Function 'add_checked' has no kernel matching input types"
            with pytest.raises(NotImplementedError, match=msg):
                # TODO: this is a pa.lib.ArrowNotImplementedError, might
                #  be better to reraise a TypeError; more consistent with
                #  non-pyarrow cases
                ser + data

            return

        if (pa_version_under8p0 and pa.types.is_duration(pa_dtype)) or (
            pa.types.is_boolean(pa_dtype)
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
    def test_compare_array(self, data, comparison_op, na_value, request):
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


class TestLogicalOps:
    """Various Series and DataFrame logical ops methods."""

    def test_kleene_or(self):
        a = pd.Series([True] * 3 + [False] * 3 + [None] * 3, dtype="boolean[pyarrow]")
        b = pd.Series([True, False, None] * 3, dtype="boolean[pyarrow]")
        result = a | b
        expected = pd.Series(
            [True, True, True, True, False, None, True, None, None],
            dtype="boolean[pyarrow]",
        )
        tm.assert_series_equal(result, expected)

        result = b | a
        tm.assert_series_equal(result, expected)

        # ensure we haven't mutated anything inplace
        tm.assert_series_equal(
            a,
            pd.Series([True] * 3 + [False] * 3 + [None] * 3, dtype="boolean[pyarrow]"),
        )
        tm.assert_series_equal(
            b, pd.Series([True, False, None] * 3, dtype="boolean[pyarrow]")
        )

    @pytest.mark.parametrize(
        "other, expected",
        [
            (None, [True, None, None]),
            (pd.NA, [True, None, None]),
            (True, [True, True, True]),
            (np.bool_(True), [True, True, True]),
            (False, [True, False, None]),
            (np.bool_(False), [True, False, None]),
        ],
    )
    def test_kleene_or_scalar(self, other, expected):
        a = pd.Series([True, False, None], dtype="boolean[pyarrow]")
        result = a | other
        expected = pd.Series(expected, dtype="boolean[pyarrow]")
        tm.assert_series_equal(result, expected)

        result = other | a
        tm.assert_series_equal(result, expected)

        # ensure we haven't mutated anything inplace
        tm.assert_series_equal(
            a, pd.Series([True, False, None], dtype="boolean[pyarrow]")
        )

    def test_kleene_and(self):
        a = pd.Series([True] * 3 + [False] * 3 + [None] * 3, dtype="boolean[pyarrow]")
        b = pd.Series([True, False, None] * 3, dtype="boolean[pyarrow]")
        result = a & b
        expected = pd.Series(
            [True, False, None, False, False, False, None, False, None],
            dtype="boolean[pyarrow]",
        )
        tm.assert_series_equal(result, expected)

        result = b & a
        tm.assert_series_equal(result, expected)

        # ensure we haven't mutated anything inplace
        tm.assert_series_equal(
            a,
            pd.Series([True] * 3 + [False] * 3 + [None] * 3, dtype="boolean[pyarrow]"),
        )
        tm.assert_series_equal(
            b, pd.Series([True, False, None] * 3, dtype="boolean[pyarrow]")
        )

    @pytest.mark.parametrize(
        "other, expected",
        [
            (None, [None, False, None]),
            (pd.NA, [None, False, None]),
            (True, [True, False, None]),
            (False, [False, False, False]),
            (np.bool_(True), [True, False, None]),
            (np.bool_(False), [False, False, False]),
        ],
    )
    def test_kleene_and_scalar(self, other, expected):
        a = pd.Series([True, False, None], dtype="boolean[pyarrow]")
        result = a & other
        expected = pd.Series(expected, dtype="boolean[pyarrow]")
        tm.assert_series_equal(result, expected)

        result = other & a
        tm.assert_series_equal(result, expected)

        # ensure we haven't mutated anything inplace
        tm.assert_series_equal(
            a, pd.Series([True, False, None], dtype="boolean[pyarrow]")
        )

    def test_kleene_xor(self):
        a = pd.Series([True] * 3 + [False] * 3 + [None] * 3, dtype="boolean[pyarrow]")
        b = pd.Series([True, False, None] * 3, dtype="boolean[pyarrow]")
        result = a ^ b
        expected = pd.Series(
            [False, True, None, True, False, None, None, None, None],
            dtype="boolean[pyarrow]",
        )
        tm.assert_series_equal(result, expected)

        result = b ^ a
        tm.assert_series_equal(result, expected)

        # ensure we haven't mutated anything inplace
        tm.assert_series_equal(
            a,
            pd.Series([True] * 3 + [False] * 3 + [None] * 3, dtype="boolean[pyarrow]"),
        )
        tm.assert_series_equal(
            b, pd.Series([True, False, None] * 3, dtype="boolean[pyarrow]")
        )

    @pytest.mark.parametrize(
        "other, expected",
        [
            (None, [None, None, None]),
            (pd.NA, [None, None, None]),
            (True, [False, True, None]),
            (np.bool_(True), [False, True, None]),
            (np.bool_(False), [True, False, None]),
        ],
    )
    def test_kleene_xor_scalar(self, other, expected):
        a = pd.Series([True, False, None], dtype="boolean[pyarrow]")
        result = a ^ other
        expected = pd.Series(expected, dtype="boolean[pyarrow]")
        tm.assert_series_equal(result, expected)

        result = other ^ a
        tm.assert_series_equal(result, expected)

        # ensure we haven't mutated anything inplace
        tm.assert_series_equal(
            a, pd.Series([True, False, None], dtype="boolean[pyarrow]")
        )


def test_arrowdtype_construct_from_string_type_with_unsupported_parameters():
    with pytest.raises(NotImplementedError, match="Passing pyarrow type"):
        ArrowDtype.construct_from_string("not_a_real_dype[s, tz=UTC][pyarrow]")

    # but as of GH#50689, timestamptz is supported
    dtype = ArrowDtype.construct_from_string("timestamp[s, tz=UTC][pyarrow]")
    expected = ArrowDtype(pa.timestamp("s", "UTC"))
    assert dtype == expected

    with pytest.raises(NotImplementedError, match="Passing pyarrow type"):
        ArrowDtype.construct_from_string("decimal(7, 2)[pyarrow]")


def test_arrowdtype_construct_from_string_type_only_one_pyarrow():
    # GH#51225
    invalid = "int64[pyarrow]foobar[pyarrow]"
    msg = (
        r"Passing pyarrow type specific parameters \(\[pyarrow\]\) in the "
        r"string is not supported\."
    )
    with pytest.raises(NotImplementedError, match=msg):
        pd.Series(range(3), dtype=invalid)


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

    if (
        pa.types.is_integer(pa_dtype)
        or pa.types.is_floating(pa_dtype)
        or (pa.types.is_decimal(pa_dtype) and not pa_version_under7p0)
    ):
        pass
    elif pa.types.is_temporal(data._data.type):
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

    if pa.types.is_timestamp(pa_dtype) and interpolation not in ["lower", "higher"]:
        # rounding error will make the check below fail
        #  (e.g. '2020-01-01 01:01:01.000001' vs '2020-01-01 01:01:01.000001024'),
        #  so we'll check for now that we match the numpy analogue
        if pa_dtype.tz:
            pd_dtype = f"M8[{pa_dtype.unit}, {pa_dtype.tz}]"
        else:
            pd_dtype = f"M8[{pa_dtype.unit}]"
        ser_np = ser.astype(pd_dtype)

        expected = ser_np.quantile(q=quantile, interpolation=interpolation)
        if quantile == 0.5:
            if pa_dtype.unit == "us":
                expected = expected.to_pydatetime(warn=False)
            assert result == expected
        else:
            if pa_dtype.unit == "us":
                expected = expected.dt.floor("us")
            tm.assert_series_equal(result, expected.astype(data.dtype))
        return

    if quantile == 0.5:
        assert result == data[0]
    else:
        # Just check the values
        expected = pd.Series(data.take([0, 0]), index=[0.5, 0.5])
        if (
            pa.types.is_integer(pa_dtype)
            or pa.types.is_floating(pa_dtype)
            or pa.types.is_decimal(pa_dtype)
        ):
            expected = expected.astype("float64[pyarrow]")
            result = result.astype("float64[pyarrow]")
        tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "take_idx, exp_idx",
    [[[0, 0, 2, 2, 4, 4], [0, 4]], [[0, 0, 0, 2, 4, 4], [0]]],
    ids=["multi_mode", "single_mode"],
)
def test_mode_dropna_true(data_for_grouping, take_idx, exp_idx):
    data = data_for_grouping.take(take_idx)
    ser = pd.Series(data)
    result = ser.mode(dropna=True)
    expected = pd.Series(data_for_grouping.take(exp_idx))
    tm.assert_series_equal(result, expected)


def test_mode_dropna_false_mode_na(data):
    # GH 50982
    more_nans = pd.Series([None, None, data[0]], dtype=data.dtype)
    result = more_nans.mode(dropna=False)
    expected = pd.Series([None], dtype=data.dtype)
    tm.assert_series_equal(result, expected)

    expected = pd.Series([None, data[0]], dtype=data.dtype)
    result = expected.mode(dropna=False)
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "arrow_dtype, expected_type",
    [
        [pa.binary(), bytes],
        [pa.binary(16), bytes],
        [pa.large_binary(), bytes],
        [pa.large_string(), str],
        [pa.list_(pa.int64()), list],
        [pa.large_list(pa.int64()), list],
        [pa.map_(pa.string(), pa.int64()), dict],
        [pa.dictionary(pa.int64(), pa.int64()), CategoricalDtypeType],
    ],
)
def test_arrow_dtype_type(arrow_dtype, expected_type):
    # GH 51845
    # TODO: Redundant with test_getitem_scalar once arrow_dtype exists in data fixture
    assert ArrowDtype(arrow_dtype).type == expected_type


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


def test_to_numpy_int_with_na():
    # GH51227: ensure to_numpy does not convert int to float
    data = [1, None]
    arr = pd.array(data, dtype="int64[pyarrow]")
    result = arr.to_numpy()
    expected = np.array([1, pd.NA], dtype=object)
    assert isinstance(result[0], int)
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


@pytest.mark.parametrize("pat", ["abc", "a[a-z]{2}"])
def test_str_count(pat):
    ser = pd.Series(["abc", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.count(pat)
    expected = pd.Series([1, None], dtype=ArrowDtype(pa.int32()))
    tm.assert_series_equal(result, expected)


def test_str_count_flags_unsupported():
    ser = pd.Series(["abc", None], dtype=ArrowDtype(pa.string()))
    with pytest.raises(NotImplementedError, match="count not"):
        ser.str.count("abc", flags=1)


@pytest.mark.parametrize(
    "side, str_func", [["left", "rjust"], ["right", "ljust"], ["both", "center"]]
)
def test_str_pad(side, str_func):
    ser = pd.Series(["a", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.pad(width=3, side=side, fillchar="x")
    expected = pd.Series(
        [getattr("a", str_func)(3, "x"), None], dtype=ArrowDtype(pa.string())
    )
    tm.assert_series_equal(result, expected)


def test_str_pad_invalid_side():
    ser = pd.Series(["a", None], dtype=ArrowDtype(pa.string()))
    with pytest.raises(ValueError, match="Invalid side: foo"):
        ser.str.pad(3, "foo", "x")


@pytest.mark.parametrize(
    "pat, case, na, regex, exp",
    [
        ["ab", False, None, False, [True, None]],
        ["Ab", True, None, False, [False, None]],
        ["ab", False, True, False, [True, True]],
        ["a[a-z]{1}", False, None, True, [True, None]],
        ["A[a-z]{1}", True, None, True, [False, None]],
    ],
)
def test_str_contains(pat, case, na, regex, exp):
    ser = pd.Series(["abc", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.contains(pat, case=case, na=na, regex=regex)
    expected = pd.Series(exp, dtype=ArrowDtype(pa.bool_()))
    tm.assert_series_equal(result, expected)


def test_str_contains_flags_unsupported():
    ser = pd.Series(["abc", None], dtype=ArrowDtype(pa.string()))
    with pytest.raises(NotImplementedError, match="contains not"):
        ser.str.contains("a", flags=1)


@pytest.mark.parametrize(
    "side, pat, na, exp",
    [
        ["startswith", "ab", None, [True, None]],
        ["startswith", "b", False, [False, False]],
        ["endswith", "b", True, [False, True]],
        ["endswith", "bc", None, [True, None]],
    ],
)
def test_str_start_ends_with(side, pat, na, exp):
    ser = pd.Series(["abc", None], dtype=ArrowDtype(pa.string()))
    result = getattr(ser.str, side)(pat, na=na)
    expected = pd.Series(exp, dtype=ArrowDtype(pa.bool_()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "arg_name, arg",
    [["pat", re.compile("b")], ["repl", str], ["case", False], ["flags", 1]],
)
def test_str_replace_unsupported(arg_name, arg):
    ser = pd.Series(["abc", None], dtype=ArrowDtype(pa.string()))
    kwargs = {"pat": "b", "repl": "x", "regex": True}
    kwargs[arg_name] = arg
    with pytest.raises(NotImplementedError, match="replace is not supported"):
        ser.str.replace(**kwargs)


@pytest.mark.parametrize(
    "pat, repl, n, regex, exp",
    [
        ["a", "x", -1, False, ["xbxc", None]],
        ["a", "x", 1, False, ["xbac", None]],
        ["[a-b]", "x", -1, True, ["xxxc", None]],
    ],
)
def test_str_replace(pat, repl, n, regex, exp):
    ser = pd.Series(["abac", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.replace(pat, repl, n=n, regex=regex)
    expected = pd.Series(exp, dtype=ArrowDtype(pa.string()))
    tm.assert_series_equal(result, expected)


def test_str_repeat_unsupported():
    ser = pd.Series(["abc", None], dtype=ArrowDtype(pa.string()))
    with pytest.raises(NotImplementedError, match="repeat is not"):
        ser.str.repeat([1, 2])


@pytest.mark.xfail(
    pa_version_under7p0,
    reason="Unsupported for pyarrow < 7",
    raises=NotImplementedError,
)
def test_str_repeat():
    ser = pd.Series(["abc", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.repeat(2)
    expected = pd.Series(["abcabc", None], dtype=ArrowDtype(pa.string()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "pat, case, na, exp",
    [
        ["ab", False, None, [True, None]],
        ["Ab", True, None, [False, None]],
        ["bc", True, None, [False, None]],
        ["ab", False, True, [True, True]],
        ["a[a-z]{1}", False, None, [True, None]],
        ["A[a-z]{1}", True, None, [False, None]],
    ],
)
def test_str_match(pat, case, na, exp):
    ser = pd.Series(["abc", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.match(pat, case=case, na=na)
    expected = pd.Series(exp, dtype=ArrowDtype(pa.bool_()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "pat, case, na, exp",
    [
        ["abc", False, None, [True, None]],
        ["Abc", True, None, [False, None]],
        ["bc", True, None, [False, None]],
        ["ab", False, True, [True, True]],
        ["a[a-z]{2}", False, None, [True, None]],
        ["A[a-z]{1}", True, None, [False, None]],
    ],
)
def test_str_fullmatch(pat, case, na, exp):
    ser = pd.Series(["abc", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.match(pat, case=case, na=na)
    expected = pd.Series(exp, dtype=ArrowDtype(pa.bool_()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "sub, start, end, exp, exp_typ",
    [["ab", 0, None, [0, None], pa.int32()], ["bc", 1, 3, [2, None], pa.int64()]],
)
def test_str_find(sub, start, end, exp, exp_typ):
    ser = pd.Series(["abc", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.find(sub, start=start, end=end)
    expected = pd.Series(exp, dtype=ArrowDtype(exp_typ))
    tm.assert_series_equal(result, expected)


def test_str_find_notimplemented():
    ser = pd.Series(["abc", None], dtype=ArrowDtype(pa.string()))
    with pytest.raises(NotImplementedError, match="find not implemented"):
        ser.str.find("ab", start=1)


@pytest.mark.parametrize(
    "i, exp",
    [
        [1, ["b", "e", None]],
        [-1, ["c", "e", None]],
        [2, ["c", None, None]],
        [-3, ["a", None, None]],
        [4, [None, None, None]],
    ],
)
def test_str_get(i, exp):
    ser = pd.Series(["abc", "de", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.get(i)
    expected = pd.Series(exp, dtype=ArrowDtype(pa.string()))
    tm.assert_series_equal(result, expected)


@pytest.mark.xfail(
    reason="TODO: StringMethods._validate should support Arrow list types",
    raises=AttributeError,
)
def test_str_join():
    ser = pd.Series(ArrowExtensionArray(pa.array([list("abc"), list("123"), None])))
    result = ser.str.join("=")
    expected = pd.Series(["a=b=c", "1=2=3", None], dtype=ArrowDtype(pa.string()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "start, stop, step, exp",
    [
        [None, 2, None, ["ab", None]],
        [None, 2, 1, ["ab", None]],
        [1, 3, 1, ["bc", None]],
    ],
)
def test_str_slice(start, stop, step, exp):
    ser = pd.Series(["abcd", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.slice(start, stop, step)
    expected = pd.Series(exp, dtype=ArrowDtype(pa.string()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "start, stop, repl, exp",
    [
        [1, 2, "x", ["axcd", None]],
        [None, 2, "x", ["xcd", None]],
        [None, 2, None, ["cd", None]],
    ],
)
def test_str_slice_replace(start, stop, repl, exp):
    ser = pd.Series(["abcd", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.slice_replace(start, stop, repl)
    expected = pd.Series(exp, dtype=ArrowDtype(pa.string()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "value, method, exp",
    [
        ["a1c", "isalnum", True],
        ["!|,", "isalnum", False],
        ["aaa", "isalpha", True],
        ["!!!", "isalpha", False],
        ["٠", "isdecimal", True],
        ["~!", "isdecimal", False],
        ["2", "isdigit", True],
        ["~", "isdigit", False],
        ["aaa", "islower", True],
        ["aaA", "islower", False],
        ["123", "isnumeric", True],
        ["11I", "isnumeric", False],
        [" ", "isspace", True],
        ["", "isspace", False],
        ["The That", "istitle", True],
        ["the That", "istitle", False],
        ["AAA", "isupper", True],
        ["AAc", "isupper", False],
    ],
)
def test_str_is_functions(value, method, exp):
    ser = pd.Series([value, None], dtype=ArrowDtype(pa.string()))
    result = getattr(ser.str, method)()
    expected = pd.Series([exp, None], dtype=ArrowDtype(pa.bool_()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "method, exp",
    [
        ["capitalize", "Abc def"],
        ["title", "Abc Def"],
        ["swapcase", "AbC Def"],
        ["lower", "abc def"],
        ["upper", "ABC DEF"],
    ],
)
def test_str_transform_functions(method, exp):
    ser = pd.Series(["aBc dEF", None], dtype=ArrowDtype(pa.string()))
    result = getattr(ser.str, method)()
    expected = pd.Series([exp, None], dtype=ArrowDtype(pa.string()))
    tm.assert_series_equal(result, expected)


def test_str_len():
    ser = pd.Series(["abcd", None], dtype=ArrowDtype(pa.string()))
    result = ser.str.len()
    expected = pd.Series([4, None], dtype=ArrowDtype(pa.int32()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "method, to_strip, val",
    [
        ["strip", None, " abc "],
        ["strip", "x", "xabcx"],
        ["lstrip", None, " abc"],
        ["lstrip", "x", "xabc"],
        ["rstrip", None, "abc "],
        ["rstrip", "x", "abcx"],
    ],
)
def test_str_strip(method, to_strip, val):
    ser = pd.Series([val, None], dtype=ArrowDtype(pa.string()))
    result = getattr(ser.str, method)(to_strip=to_strip)
    expected = pd.Series(["abc", None], dtype=ArrowDtype(pa.string()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("val", ["abc123", "abc"])
def test_str_removesuffix(val):
    ser = pd.Series([val, None], dtype=ArrowDtype(pa.string()))
    result = ser.str.removesuffix("123")
    expected = pd.Series(["abc", None], dtype=ArrowDtype(pa.string()))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "method, args",
    [
        ["partition", ("abc", False)],
        ["rpartition", ("abc", False)],
        ["removeprefix", ("abc",)],
        ["casefold", ()],
        ["encode", ("abc",)],
        ["extract", (r"[ab](\d)",)],
        ["findall", ("abc",)],
        ["get_dummies", ()],
        ["index", ("abc",)],
        ["rindex", ("abc",)],
        ["normalize", ("abc",)],
        ["rfind", ("abc",)],
        ["split", ()],
        ["rsplit", ()],
        ["translate", ("abc",)],
        ["wrap", ("abc",)],
    ],
)
def test_str_unsupported_methods(method, args):
    ser = pd.Series(["abc", None], dtype=ArrowDtype(pa.string()))
    with pytest.raises(
        NotImplementedError, match=f"str.{method} not supported with pd.ArrowDtype"
    ):
        getattr(ser.str, method)(*args)


@pytest.mark.parametrize("unit", ["ns", "us", "ms", "s"])
def test_duration_from_strings_with_nat(unit):
    # GH51175
    strings = ["1000", "NaT"]
    pa_type = pa.duration(unit)
    result = ArrowExtensionArray._from_sequence_of_strings(strings, dtype=pa_type)
    expected = ArrowExtensionArray(pa.array([1000, None], type=pa_type))
    tm.assert_extension_array_equal(result, expected)


def test_unsupported_dt(data):
    pa_dtype = data.dtype.pyarrow_dtype
    if not pa.types.is_temporal(pa_dtype):
        with pytest.raises(
            AttributeError, match="Can only use .dt accessor with datetimelike values"
        ):
            pd.Series(data).dt


@pytest.mark.parametrize(
    "prop, expected",
    [
        ["year", 2023],
        ["day", 2],
        ["day_of_week", 0],
        ["dayofweek", 0],
        ["weekday", 0],
        ["day_of_year", 2],
        ["dayofyear", 2],
        ["hour", 3],
        ["minute", 4],
        pytest.param(
            "is_leap_year",
            False,
            marks=pytest.mark.xfail(
                pa_version_under8p0,
                raises=NotImplementedError,
                reason="is_leap_year not implemented for pyarrow < 8.0",
            ),
        ),
        ["microsecond", 5],
        ["month", 1],
        ["nanosecond", 6],
        ["quarter", 1],
        ["second", 7],
        ["date", date(2023, 1, 2)],
        ["time", time(3, 4, 7, 5)],
    ],
)
def test_dt_properties(prop, expected):
    ser = pd.Series(
        [
            pd.Timestamp(
                year=2023,
                month=1,
                day=2,
                hour=3,
                minute=4,
                second=7,
                microsecond=5,
                nanosecond=6,
            ),
            None,
        ],
        dtype=ArrowDtype(pa.timestamp("ns")),
    )
    result = getattr(ser.dt, prop)
    exp_type = None
    if isinstance(expected, date):
        exp_type = pa.date32()
    elif isinstance(expected, time):
        exp_type = pa.time64("ns")
    expected = pd.Series(ArrowExtensionArray(pa.array([expected, None], type=exp_type)))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("unit", ["us", "ns"])
def test_dt_time_preserve_unit(unit):
    ser = pd.Series(
        [datetime(year=2023, month=1, day=2, hour=3), None],
        dtype=ArrowDtype(pa.timestamp(unit)),
    )
    result = ser.dt.time
    expected = pd.Series(
        ArrowExtensionArray(pa.array([time(3, 0), None], type=pa.time64(unit)))
    )
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("tz", [None, "UTC", "US/Pacific"])
def test_dt_tz(tz):
    ser = pd.Series(
        [datetime(year=2023, month=1, day=2, hour=3), None],
        dtype=ArrowDtype(pa.timestamp("ns", tz=tz)),
    )
    result = ser.dt.tz
    assert result == tz


def test_dt_isocalendar():
    ser = pd.Series(
        [datetime(year=2023, month=1, day=2, hour=3), None],
        dtype=ArrowDtype(pa.timestamp("ns")),
    )
    result = ser.dt.isocalendar()
    expected = pd.DataFrame(
        [[2023, 1, 1], [0, 0, 0]],
        columns=["year", "week", "day"],
        dtype="int64[pyarrow]",
    )
    tm.assert_frame_equal(result, expected)


def test_dt_strftime(request):
    if is_platform_windows() and is_ci_environment():
        request.node.add_marker(
            pytest.mark.xfail(
                raises=pa.ArrowInvalid,
                reason=(
                    "TODO: Set ARROW_TIMEZONE_DATABASE environment variable "
                    "on CI to path to the tzdata for pyarrow."
                ),
            )
        )
    ser = pd.Series(
        [datetime(year=2023, month=1, day=2, hour=3), None],
        dtype=ArrowDtype(pa.timestamp("ns")),
    )
    result = ser.dt.strftime("%Y-%m-%dT%H:%M:%S")
    expected = pd.Series(
        ["2023-01-02T03:00:00.000000000", None], dtype=ArrowDtype(pa.string())
    )
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("method", ["ceil", "floor", "round"])
def test_dt_roundlike_tz_options_not_supported(method):
    ser = pd.Series(
        [datetime(year=2023, month=1, day=2, hour=3), None],
        dtype=ArrowDtype(pa.timestamp("ns")),
    )
    with pytest.raises(NotImplementedError, match="ambiguous is not supported."):
        getattr(ser.dt, method)("1H", ambiguous="NaT")

    with pytest.raises(NotImplementedError, match="nonexistent is not supported."):
        getattr(ser.dt, method)("1H", nonexistent="NaT")


@pytest.mark.parametrize("method", ["ceil", "floor", "round"])
def test_dt_roundlike_unsupported_freq(method):
    ser = pd.Series(
        [datetime(year=2023, month=1, day=2, hour=3), None],
        dtype=ArrowDtype(pa.timestamp("ns")),
    )
    with pytest.raises(ValueError, match="freq='1B' is not supported"):
        getattr(ser.dt, method)("1B")

    with pytest.raises(ValueError, match="Must specify a valid frequency: None"):
        getattr(ser.dt, method)(None)


@pytest.mark.xfail(
    pa_version_under7p0, reason="Methods not supported for pyarrow < 7.0"
)
@pytest.mark.parametrize("freq", ["D", "H", "T", "S", "L", "U", "N"])
@pytest.mark.parametrize("method", ["ceil", "floor", "round"])
def test_dt_ceil_year_floor(freq, method):
    ser = pd.Series(
        [datetime(year=2023, month=1, day=1), None],
    )
    pa_dtype = ArrowDtype(pa.timestamp("ns"))
    expected = getattr(ser.dt, method)(f"1{freq}").astype(pa_dtype)
    result = getattr(ser.astype(pa_dtype).dt, method)(f"1{freq}")
    tm.assert_series_equal(result, expected)


def test_dt_to_pydatetime():
    # GH 51859
    data = [datetime(2022, 1, 1), datetime(2023, 1, 1)]
    ser = pd.Series(data, dtype=ArrowDtype(pa.timestamp("ns")))

    result = ser.dt.to_pydatetime()
    expected = np.array(data, dtype=object)
    tm.assert_numpy_array_equal(result, expected)
    assert all(type(res) is datetime for res in result)

    expected = ser.astype("datetime64[ns]").dt.to_pydatetime()
    tm.assert_numpy_array_equal(result, expected)


def test_dt_tz_localize_unsupported_tz_options():
    ser = pd.Series(
        [datetime(year=2023, month=1, day=2, hour=3), None],
        dtype=ArrowDtype(pa.timestamp("ns")),
    )
    with pytest.raises(NotImplementedError, match="ambiguous='NaT' is not supported"):
        ser.dt.tz_localize("UTC", ambiguous="NaT")

    with pytest.raises(NotImplementedError, match="nonexistent='NaT' is not supported"):
        ser.dt.tz_localize("UTC", nonexistent="NaT")


def test_dt_tz_localize_none():
    ser = pd.Series(
        [datetime(year=2023, month=1, day=2, hour=3), None],
        dtype=ArrowDtype(pa.timestamp("ns", tz="US/Pacific")),
    )
    result = ser.dt.tz_localize(None)
    expected = pd.Series(
        [datetime(year=2023, month=1, day=2, hour=3), None],
        dtype=ArrowDtype(pa.timestamp("ns")),
    )
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("unit", ["us", "ns"])
def test_dt_tz_localize(unit):
    ser = pd.Series(
        [datetime(year=2023, month=1, day=2, hour=3), None],
        dtype=ArrowDtype(pa.timestamp(unit)),
    )
    result = ser.dt.tz_localize("US/Pacific")
    expected = pd.Series(
        [datetime(year=2023, month=1, day=2, hour=3), None],
        dtype=ArrowDtype(pa.timestamp(unit, "US/Pacific")),
    )
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("skipna", [True, False])
def test_boolean_reduce_series_all_null(all_boolean_reductions, skipna):
    # GH51624
    ser = pd.Series([None], dtype="float64[pyarrow]")
    result = getattr(ser, all_boolean_reductions)(skipna=skipna)
    if skipna:
        expected = all_boolean_reductions == "all"
    else:
        expected = pd.NA
    assert result is expected


@pytest.mark.parametrize("dtype", ["string", "string[pyarrow]"])
def test_series_from_string_array(dtype):
    arr = pa.array("the quick brown fox".split())
    ser = pd.Series(arr, dtype=dtype)
    expected = pd.Series(ArrowExtensionArray(arr), dtype=dtype)
    tm.assert_series_equal(ser, expected)


def test_setitem_boolean_replace_with_mask_segfault():
    # GH#52059
    N = 145_000
    arr = ArrowExtensionArray(pa.chunked_array([np.ones((N,), dtype=np.bool_)]))
    expected = arr.copy()
    arr[np.zeros((N,), dtype=np.bool_)] = False
    assert arr._data == expected._data
