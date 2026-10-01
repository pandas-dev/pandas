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

from __future__ import annotations

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
import operator
import re

import numpy as np
import pytest

from pandas.compat import (
    PY312,
    is_platform_windows,
    pa_version_under20p0,
    pa_version_under21p0,
)
from pandas.compat.pyarrow import pa_version_under22p0
from pandas.errors import (
    Pandas4Warning,
)

from pandas.core.dtypes.common import pandas_dtype
from pandas.core.dtypes.dtypes import (
    ArrowDtype,
)

import pandas as pd
import pandas._testing as tm
from pandas.api.extensions import no_default
from pandas.api.types import (
    is_string_dtype,
)
from pandas.tests.extension import base

pa = pytest.importorskip("pyarrow")

from pandas.core.arrays.arrow.array import ArrowExtensionArray

# GH#62423; matched instead of the leading clause, whose wording differs per warn site
depr_msg = "In a future version these will be treated as scalar-like"


def _require_timezone_database(request):
    if is_platform_windows() and pa_version_under22p0:
        mark = pytest.mark.xfail(
            raises=pa.ArrowInvalid,
            reason=(
                "TODO: Set ARROW_TIMEZONE_DATABASE environment variable "
                "on CI to path to the tzdata for pyarrow."
            ),
        )
        request.applymarker(mark)


@pytest.fixture(params=tm.ALL_PYARROW_DTYPES, ids=str)
def dtype(request):
    return ArrowDtype(pyarrow_dtype=request.param)


@pytest.fixture
def data(dtype):
    pa_dtype = dtype.pyarrow_dtype
    if pa.types.is_boolean(pa_dtype):
        data = [True, False] * 2 + [None] + [True, False] + [None] + [True, False]
    elif pa.types.is_floating(pa_dtype):
        data = [1.0, 0.0] * 2 + [None] + [-2.0, -1.0] + [None] + [0.5, 99.5]
    elif pa.types.is_signed_integer(pa_dtype):
        data = [1, 0] * 2 + [None] + [-2, -1] + [None] + [1, 99]
    elif pa.types.is_unsigned_integer(pa_dtype):
        data = [1, 0] * 2 + [None] + [2, 1] + [None] + [1, 99]
    elif pa.types.is_decimal(pa_dtype):
        data = (
            [Decimal("1"), Decimal("0.0")] * 2
            + [None]
            + [Decimal("-2.0"), Decimal("-1.0")]
            + [None]
            + [Decimal("0.5"), Decimal("33.123")]
        )
    elif pa.types.is_date(pa_dtype):
        data = (
            [date(2022, 1, 1), date(1999, 12, 31)] * 2
            + [None]
            + [date(2022, 1, 1), date(2022, 1, 1)]
            + [None]
            + [date(1999, 12, 31), date(1999, 12, 31)]
        )
    elif pa.types.is_timestamp(pa_dtype):
        data = (
            [datetime(2020, 1, 1, 1, 1, 1, 1), datetime(1999, 1, 1, 1, 1, 1, 1)] * 2
            + [None]
            + [datetime(2020, 1, 1, 1), datetime(1999, 1, 1, 1)]
            + [None]
            + [datetime(2020, 1, 1), datetime(1999, 1, 1)]
        )
    elif pa.types.is_duration(pa_dtype):
        data = (
            [timedelta(1), timedelta(1, 1)] * 2
            + [None]
            + [timedelta(-1), timedelta(0)]
            + [None]
            + [timedelta(-10), timedelta(10)]
        )
    elif pa.types.is_time(pa_dtype):
        data = (
            [time(12, 0), time(0, 12)] * 2
            + [None]
            + [time(0, 0), time(1, 1)]
            + [None]
            + [time(0, 5), time(5, 0)]
        )
    elif pa.types.is_string(pa_dtype):
        data = ["a", "b"] * 2 + [None] + ["1", "2"] + [None] + ["!", ">"]
    elif pa.types.is_binary(pa_dtype):
        data = [b"a", b"b"] * 2 + [None] + [b"1", b"2"] + [None] + [b"!", b">"]
    else:
        raise NotImplementedError
    return pd.array(data, dtype=dtype)


@pytest.fixture
def data_missing(data):
    """Length-2 array with [NA, Valid]"""
    return type(data)._from_sequence([None, data[0]], dtype=data.dtype)


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
        [data_for_grouping[0], data_for_grouping[7], data_for_grouping[4]],
        dtype=data_for_grouping.dtype,
    )


@pytest.fixture
def data_missing_for_sorting(data_for_grouping):
    """
    Length-3 array with a known sort order.

    This should be three items [B, NA, A] with
    A < B and NA missing.
    """
    return type(data_for_grouping)._from_sequence(
        [data_for_grouping[0], data_for_grouping[2], data_for_grouping[4]],
        dtype=data_for_grouping.dtype,
    )


@pytest.fixture
def data_for_twos(data):
    """Length-100 array in which all the elements are two."""
    pa_dtype = data.dtype.pyarrow_dtype
    if (
        pa.types.is_integer(pa_dtype)
        or pa.types.is_floating(pa_dtype)
        or pa.types.is_decimal(pa_dtype)
        or pa.types.is_duration(pa_dtype)
    ):
        return pd.array([2] * 10, dtype=data.dtype)
    # tests will be xfailed where 2 is not a valid scalar for pa_dtype
    return data
    # TODO: skip otherwise?


class TestArrowArray(base.ExtensionTests):
    def _honors_copy_keyword(self, data) -> bool:
        return False

    def _construct_for_combine_add(self, left, right):
        dtype = left.dtype

        # in a couple cases, addition is not dtype-preserving
        if dtype == "bool[pyarrow]":
            dtype = pandas_dtype("int64[pyarrow]")
        elif dtype == "int8[pyarrow]" and isinstance(right, type(left)):
            dtype = pandas_dtype("int64[pyarrow]")

        if isinstance(right, type(left)):
            return left._from_sequence(
                [a + b for (a, b) in zip(list(left), list(right), strict=True)],
                dtype=dtype,
            )
        else:
            return left._from_sequence(
                [a + right for a in list(left)],
                dtype=dtype,
            )

    def test_compare_range_len(self, data, comparison_op):
        # GH#63429 a range compares elementwise like the equivalent list.
        #  Note we can't go through _compare_other here: its pointwise
        #  expectation uses Series.combine, which treats the range as a scalar
        #  and so matches an all-False result too.
        ser = pd.Series(data)
        rng = range(len(ser))

        try:
            expected = comparison_op(ser, list(rng))
        except Exception as err:
            with tm.assert_produces_warning(Pandas4Warning, match=depr_msg):
                with pytest.raises(type(err)):
                    comparison_op(ser, rng)
            return

        with tm.assert_produces_warning(Pandas4Warning, match=depr_msg):
            result = comparison_op(ser, rng)
        tm.assert_series_equal(result, expected)

    def test_compare_range_mismatched_len(self, data, comparison_op):
        # GH#63429 the length check must not be bypassed for a range
        ser = pd.Series(data)
        rng = range(len(ser) + 1)

        with tm.assert_produces_warning(Pandas4Warning, match=depr_msg):
            with pytest.raises(ValueError, match="Lengths must match to compare"):
                comparison_op(ser, rng)

    def test_astype_str(self, data, request, using_infer_string):
        pa_dtype = data.dtype.pyarrow_dtype
        if pa.types.is_binary(pa_dtype):
            request.applymarker(
                pytest.mark.xfail(
                    reason=f"For {pa_dtype} .astype(str) decodes.",
                )
            )
        elif not using_infer_string and (
            (pa.types.is_timestamp(pa_dtype) and pa_dtype.tz is None)
            or pa.types.is_duration(pa_dtype)
        ):
            request.applymarker(
                pytest.mark.xfail(
                    reason="pd.Timestamp/pd.Timedelta repr different from numpy repr",
                )
            )
        super().test_astype_str(data)

    def test_from_dtype(self, data, request):
        pa_dtype = data.dtype.pyarrow_dtype
        if pa.types.is_string(pa_dtype) or pa.types.is_decimal(pa_dtype):
            if pa.types.is_string(pa_dtype):
                reason = "ArrowDtype(pa.string()) != StringDtype('pyarrow')"
            else:
                reason = f"pyarrow.type_for_alias cannot infer {pa_dtype}"

            request.applymarker(
                pytest.mark.xfail(
                    reason=reason,
                )
            )
        super().test_from_dtype(data)

    def test_from_sequence_pa_array(self, data):
        # https://github.com/pandas-dev/pandas/pull/47034#discussion_r955500784
        # data._pa_array = pa.ChunkedArray
        result = type(data)._from_sequence(data._pa_array, dtype=data.dtype)
        tm.assert_extension_array_equal(result, data)
        assert isinstance(result._pa_array, pa.ChunkedArray)

        result = type(data)._from_sequence(
            data._pa_array.combine_chunks(), dtype=data.dtype
        )
        tm.assert_extension_array_equal(result, data)
        assert isinstance(result._pa_array, pa.ChunkedArray)

    def test_from_sequence_pa_array_notimplemented(self, request):
        dtype = ArrowDtype(pa.month_day_nano_interval())
        with pytest.raises(NotImplementedError, match="Converting strings to"):
            ArrowExtensionArray._from_sequence_of_strings(["12-1"], dtype=dtype)

    def test_from_sequence_of_strings_pa_array(self, data, request):
        pa_dtype = data.dtype.pyarrow_dtype
        if pa.types.is_timestamp(pa_dtype) and pa_dtype.tz is not None:
            _require_timezone_database(request)

        pa_array = data._pa_array.cast(pa.string())
        result = type(data)._from_sequence_of_strings(pa_array, dtype=data.dtype)
        tm.assert_extension_array_equal(result, data)

        pa_array = pa_array.combine_chunks()
        result = type(data)._from_sequence_of_strings(pa_array, dtype=data.dtype)
        tm.assert_extension_array_equal(result, data)

    def check_accumulate(self, ser, op_name, skipna):
        result = getattr(ser, op_name)(skipna=skipna)

        pa_type = ser.dtype.pyarrow_dtype
        if pa.types.is_temporal(pa_type):
            # Just check that we match the integer behavior.
            if pa_type.bit_width == 32:
                int_type = "int32[pyarrow]"
            else:
                int_type = "int64[pyarrow]"
            ser = ser.astype(int_type)
            result = result.astype(int_type)

        result = result.astype("Float64")
        expected = getattr(ser.astype("Float64"), op_name)(skipna=skipna)
        tm.assert_series_equal(result, expected, check_dtype=False)

    def _supports_accumulation(self, ser: pd.Series, op_name: str) -> bool:
        # error: Item "dtype[Any]" of "dtype[Any] | ExtensionDtype" has no
        # attribute "pyarrow_dtype"
        pa_type = ser.dtype.pyarrow_dtype  # type: ignore[union-attr]

        if pa.types.is_binary(pa_type) or pa.types.is_decimal(pa_type):
            if op_name in ["cumsum", "cumprod", "cummax", "cummin"]:
                return False
        elif pa.types.is_string(pa_type):
            if op_name == "cumprod":
                return False
        elif pa.types.is_boolean(pa_type):
            if op_name in ["cumprod", "cummax", "cummin"]:
                return False
        elif pa.types.is_temporal(pa_type):
            if op_name == "cumsum" and not pa.types.is_duration(pa_type):
                return False
            elif op_name == "cumprod":
                return False
        return True

    @pytest.mark.parametrize("skipna", [True, False])
    def test_accumulate_series(self, data, all_numeric_accumulations, skipna, request):
        pa_type = data.dtype.pyarrow_dtype
        op_name = all_numeric_accumulations

        if pa.types.is_string(pa_type) and op_name in ["cumsum", "cummin", "cummax"]:
            # https://github.com/pandas-dev/pandas/pull/60633
            # Doesn't fit test structure, tested in series/test_cumulative.py instead.
            return

        ser = pd.Series(data)

        if not self._supports_accumulation(ser, op_name):
            # The base class test will check that we raise
            return super().test_accumulate_series(
                data, all_numeric_accumulations, skipna
            )

        if all_numeric_accumulations == "cumsum" and (
            pa.types.is_boolean(pa_type) or pa.types.is_decimal(pa_type)
        ):
            request.applymarker(
                pytest.mark.xfail(
                    reason=f"{all_numeric_accumulations} not implemented for {pa_type}",
                    raises=TypeError,
                )
            )

        self.check_accumulate(ser, op_name, skipna)

    def _supports_reduction(self, ser: pd.Series, op_name: str) -> bool:
        if pa_version_under20p0 and op_name in ("skew", "kurt"):
            return False

        dtype = ser.dtype
        # error: Item "dtype[Any]" of "dtype[Any] | ExtensionDtype" has
        # no attribute "pyarrow_dtype"
        pa_dtype = dtype.pyarrow_dtype  # type: ignore[union-attr]
        if pa.types.is_temporal(pa_dtype) and op_name in [
            "sum",
            "var",
            "prod",
            "skew",
            "kurt",
        ]:
            if pa.types.is_duration(pa_dtype) and op_name in ["sum"]:
                # summing timedeltas is one case that *is* well-defined
                pass
            else:
                return False
        elif pa.types.is_binary(pa_dtype) and op_name in [
            "sum",
            "skew",
            "kurt",
            "any",
            "all",
        ]:
            return False
        elif (
            pa.types.is_string(pa_dtype) or pa.types.is_binary(pa_dtype)
        ) and op_name in [
            "mean",
            "median",
            "prod",
            "std",
            "sem",
            "var",
            "skew",
            "kurt",
            "any",
            "all",
        ]:
            return False

        if (
            pa.types.is_temporal(pa_dtype)
            and not pa.types.is_duration(pa_dtype)
            and op_name in ["any", "all"]
        ):
            # xref GH#34479 we support this in our non-pyarrow datetime64 dtypes,
            #  but it isn't obvious we _should_.  For now, we keep the pyarrow
            #  behavior which does not support this.
            return False

        if pa.types.is_boolean(pa_dtype) and op_name in [
            "median",
            "std",
            "var",
            "skew",
            "kurt",
            "sem",
        ]:
            return False

        return True

    def check_reduce(self, ser: pd.Series, op_name: str, skipna: bool):
        # error: Item "dtype[Any]" of "dtype[Any] | ExtensionDtype" has no
        # attribute "pyarrow_dtype"
        pa_dtype = ser.dtype.pyarrow_dtype  # type: ignore[union-attr]
        if pa.types.is_integer(pa_dtype) or pa.types.is_floating(pa_dtype):
            alt = ser.astype("Float64")
        else:
            # TODO: in the opposite case, aren't we testing... nothing? For
            # e.g. date/time dtypes trying to calculate 'expected' by converting
            # to object will raise for mean, std etc
            alt = ser

        # TODO: in the opposite case, aren't we testing... nothing?
        if op_name == "count":
            result = getattr(ser, op_name)()
            expected = getattr(alt, op_name)()
        else:
            result = getattr(ser, op_name)(skipna=skipna)
            expected = getattr(alt, op_name)(skipna=skipna)
        tm.assert_almost_equal(result, expected)

    def _get_expected_reduction_dtype(self, arr, op_name: str, skipna: bool):
        pa_type = arr._pa_array.type

        if op_name in ["max", "min"]:
            cmp_dtype = arr.dtype
        elif pa.types.is_temporal(pa_type):
            if op_name in ["std", "sem"]:
                if pa.types.is_duration(pa_type):
                    cmp_dtype = arr.dtype
                elif pa.types.is_date32(pa_type):
                    cmp_dtype = ArrowDtype(pa.duration("s"))
                elif pa.types.is_date64(pa_type):
                    cmp_dtype = ArrowDtype(pa.duration("ms"))
                elif pa.types.is_time(pa_type):
                    cmp_dtype = ArrowDtype(pa.duration(pa_type.unit))
                else:
                    cmp_dtype = ArrowDtype(pa.duration(pa_type.unit))
            else:
                cmp_dtype = arr.dtype
        elif arr.dtype.name == "decimal128(7, 3)[pyarrow]":
            if op_name == "sum" and not pa_version_under21p0:
                # https://github.com/apache/arrow/pull/44184
                cmp_dtype = ArrowDtype(pa.decimal128(38, 3))
            elif op_name not in ["median", "var", "std", "sem", "skew", "kurt"]:
                cmp_dtype = arr.dtype
            else:
                cmp_dtype = "float64[pyarrow]"
        elif op_name in ["median", "var", "std", "mean", "skew", "kurt", "sem"]:
            cmp_dtype = "float64[pyarrow]"
        elif op_name in ["sum", "prod"] and pa.types.is_boolean(pa_type):
            cmp_dtype = "uint64[pyarrow]"
        elif op_name == "sum" and pa.types.is_string(pa_type):
            cmp_dtype = arr.dtype
        else:
            cmp_dtype = {
                "i": "int64[pyarrow]",
                "u": "uint64[pyarrow]",
                "f": "float64[pyarrow]",
            }[arr.dtype.kind]
        return cmp_dtype

    @pytest.mark.parametrize("skipna", [True, False])
    def test_reduce_series_numeric(self, data, all_numeric_reductions, skipna):
        return super().test_reduce_series_numeric(data, all_numeric_reductions, skipna)

    @pytest.mark.parametrize("skipna", [True, False])
    def test_reduce_frame(self, data, all_numeric_reductions, skipna, request):
        op_name = all_numeric_reductions
        if op_name in ("skew", "kurt") and pa_version_under20p0:
            if data.dtype._is_numeric:
                mark = pytest.mark.xfail(reason=f"{op_name} not implemented")
                request.applymarker(mark)
        return super().test_reduce_frame(data, all_numeric_reductions, skipna)

    @pytest.mark.parametrize("typ", ["int64", "uint64", "float64"])
    def test_median_not_approximate(self, typ):
        # GH 52679
        result = pd.Series([1, 2], dtype=f"{typ}[pyarrow]").median()
        assert result == 1.5

    def test_construct_from_string_own_name(self, dtype, request):
        pa_dtype = dtype.pyarrow_dtype
        if pa.types.is_decimal(pa_dtype):
            request.applymarker(
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
                request.applymarker(
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

    def test_is_not_string_type(self, dtype):
        pa_dtype = dtype.pyarrow_dtype
        if pa.types.is_string(pa_dtype):
            assert is_string_dtype(dtype)
        else:
            super().test_is_not_string_type(dtype)

    @pytest.mark.xfail(
        reason="GH 45419: pyarrow.ChunkedArray does not support views.", run=False
    )
    def test_view(self, data):
        super().test_view(data)

    def test_fillna_no_op_returns_copy(self, data):
        data = data[~data.isna()]

        valid = data[0]
        result = data.fillna(valid)
        assert result is not data
        tm.assert_extension_array_equal(result, data)

    @pytest.mark.xfail(
        reason="GH 45419: pyarrow.ChunkedArray does not support views", run=False
    )
    def test_transpose(self, data):
        super().test_transpose(data)

    @pytest.mark.xfail(
        reason="GH 45419: pyarrow.ChunkedArray does not support views", run=False
    )
    def test_setitem_preserves_views(self, data):
        super().test_setitem_preserves_views(data)

    @pytest.mark.parametrize("dtype_backend", ["pyarrow", no_default])
    @pytest.mark.parametrize("engine", ["c", "python"])
    def test_EA_types(self, engine, data, dtype_backend, request, using_nan_is_na):
        pa_dtype = data.dtype.pyarrow_dtype
        if pa.types.is_decimal(pa_dtype):
            request.applymarker(
                pytest.mark.xfail(
                    raises=NotImplementedError,
                    reason=f"Parameterized types {pa_dtype} not supported.",
                )
            )
        elif pa.types.is_timestamp(pa_dtype) and pa_dtype.unit in ("us", "ns"):
            request.applymarker(
                pytest.mark.xfail(
                    raises=ValueError,
                    reason="https://github.com/pandas-dev/pandas/issues/49767",
                )
            )
        elif pa.types.is_binary(pa_dtype):
            request.applymarker(
                pytest.mark.xfail(reason="CSV parsers don't correctly handle binary")
            )
        df = pd.DataFrame({"with_dtype": pd.Series(data, dtype=str(data.dtype))})
        if not using_nan_is_na:
            csv_output = df.to_csv(index=False, na_rep="NA")
        else:
            csv_output = df.to_csv(index=False, na_rep=np.nan)
        if pa.types.is_binary(pa_dtype):
            csv_output = BytesIO(csv_output)
        else:
            csv_output = StringIO(csv_output)
        result = pd.read_csv(
            csv_output,
            dtype={"with_dtype": str(data.dtype)},
            engine=engine,
            dtype_backend=dtype_backend,
        )
        expected = df
        tm.assert_frame_equal(result, expected)

    def test_invert(self, data, request):
        pa_dtype = data.dtype.pyarrow_dtype
        if not (
            pa.types.is_boolean(pa_dtype)
            or pa.types.is_integer(pa_dtype)
            or pa.types.is_string(pa_dtype)
        ):
            request.applymarker(
                pytest.mark.xfail(
                    raises=pa.ArrowNotImplementedError,
                    reason=f"pyarrow.compute.invert does support {pa_dtype}",
                )
            )
        if PY312 and pa.types.is_boolean(pa_dtype):
            with tm.assert_produces_warning(
                DeprecationWarning, match="Bitwise inversion", check_stacklevel=False
            ):
                super().test_invert(data)
        else:
            super().test_invert(data)

    @pytest.mark.parametrize("periods", [1, -2])
    def test_diff(self, data, periods, request):
        pa_dtype = data.dtype.pyarrow_dtype
        if pa.types.is_unsigned_integer(pa_dtype) and periods == 1:
            request.applymarker(
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

    _combine_le_expected_dtype = "bool[pyarrow]"

    def get_op_from_name(self, op_name):
        short_opname = op_name.strip("_")
        if short_opname == "rtruediv":
            # use the numpy version that won't raise on division by zero

            def rtruediv(x, y):
                return np.divide(y, x)

            return rtruediv
        elif short_opname == "rfloordiv":
            return lambda x, y: np.floor_divide(y, x)

        return tm.get_op_from_name(op_name)

    # TODO: use EA._cast_pointwise_result, same with other test files that
    #  override this
    def _cast_pointwise_result(self, op_name: str, obj, other, pointwise_result):
        # BaseOpsUtil._combine can upcast expected dtype
        # (because it generates expected on python scalars)
        # while ArrowExtensionArray maintains original type
        expected = pointwise_result

        if op_name in ["eq", "ne", "lt", "le", "gt", "ge"]:
            return pointwise_result.astype("boolean[pyarrow]")

        original_dtype = tm.get_dtype(expected)

        was_frame = False
        if isinstance(expected, pd.DataFrame):
            was_frame = True
            expected_data = expected.iloc[:, 0]
        else:
            expected_data = expected

        # the pointwise method will have retained our original dtype, while
        #  the op(ser, other) version will have cast to 64bit
        if type(other) is int and op_name not in ["__floordiv__"]:
            if original_dtype.kind == "f":
                return expected.astype("float64[pyarrow]")
            else:
                return expected.astype("int64[pyarrow]")
        elif type(other) is float:
            return expected.astype("float64[pyarrow]")

        # error: Item "ExtensionDtype" of "dtype[Any] | ExtensionDtype" has
        #  no attribute "pyarrow_dtype"
        orig_pa_type = original_dtype.pyarrow_dtype  # type: ignore[union-attr]
        if not was_frame and isinstance(other, pd.Series):
            # i.e. test_arith_series_with_array
            if not (
                pa.types.is_floating(orig_pa_type)
                or (
                    pa.types.is_integer(orig_pa_type)
                    and op_name not in ["__truediv__", "__rtruediv__"]
                )
                or pa.types.is_duration(orig_pa_type)
                or pa.types.is_timestamp(orig_pa_type)
                or pa.types.is_date(orig_pa_type)
                or pa.types.is_decimal(orig_pa_type)
            ):
                # base class _combine always returns int64, while
                #  ArrowExtensionArray does not upcast
                return expected
        elif not (
            (op_name == "__floordiv__" and pa.types.is_integer(orig_pa_type))
            or pa.types.is_duration(orig_pa_type)
            or pa.types.is_timestamp(orig_pa_type)
            or pa.types.is_date(orig_pa_type)
            or pa.types.is_decimal(orig_pa_type)
        ):
            # base class _combine always returns int64, while
            #  ArrowExtensionArray does not upcast
            return expected

        pa_expected = pa.array(expected_data._values)

        if pa.types.is_decimal(pa_expected.type) and pa.types.is_decimal(orig_pa_type):
            # decimal precision can resize in the result type depending on data
            # just compare the float values
            alt = getattr(obj, op_name)(other)
            alt_dtype = tm.get_dtype(alt)
            assert isinstance(alt_dtype, ArrowDtype)
            if op_name == "__pow__" and isinstance(other, Decimal):
                # TODO: would it make more sense to retain Decimal here?
                alt_dtype = ArrowDtype(pa.float64())
            elif (
                op_name == "__pow__"
                and isinstance(other, pd.Series)
                and other.dtype == original_dtype
            ):
                # TODO: would it make more sense to retain Decimal here?
                alt_dtype = ArrowDtype(pa.float64())
            else:
                assert pa.types.is_decimal(alt_dtype.pyarrow_dtype)
            return expected.astype(alt_dtype)

        else:
            pa_expected = pa_expected.cast(orig_pa_type)

        pd_expected = type(expected_data._values)(pa_expected)
        if was_frame:
            expected = pd.DataFrame(
                pd_expected, index=expected.index, columns=expected.columns
            )
        else:
            expected = pd.Series(pd_expected)
        return expected

    def _is_temporal_supported(self, opname, pa_dtype):
        return (
            opname
            in (
                "__add__",
                "__radd__",
                "__truediv__",
                "__rtruediv__",
                "__floordiv__",
                "__rfloordiv__",
            )
            and pa.types.is_duration(pa_dtype)
        ) or (opname in ("__sub__", "__rsub__") and pa.types.is_temporal(pa_dtype))

    def _get_expected_exception(
        self, op_name: str, obj, other
    ) -> type[Exception] | tuple[type[Exception], ...] | None:
        if op_name in ("__divmod__", "__rdivmod__"):
            return (NotImplementedError, TypeError)

        exc: type[Exception] | tuple[type[Exception], ...] | None
        dtype = tm.get_dtype(obj)
        # error: Item "dtype[Any]" of "dtype[Any] | ExtensionDtype" has no
        # attribute "pyarrow_dtype"
        pa_dtype = dtype.pyarrow_dtype  # type: ignore[union-attr]

        arrow_temporal_supported = self._is_temporal_supported(op_name, pa_dtype)
        if op_name in {
            "__mod__",
            "__rmod__",
        }:
            exc = (NotImplementedError, TypeError)
        elif arrow_temporal_supported:
            exc = None
        elif op_name in ["__add__", "__radd__"] and (
            pa.types.is_string(pa_dtype) or pa.types.is_binary(pa_dtype)
        ):
            exc = None
        elif not (
            pa.types.is_floating(pa_dtype)
            or pa.types.is_integer(pa_dtype)
            or pa.types.is_decimal(pa_dtype)
        ):
            exc = TypeError
        else:
            exc = None
        return exc

    def _get_arith_xfail_marker(self, opname, pa_dtype):
        mark = None

        arrow_temporal_supported = self._is_temporal_supported(opname, pa_dtype)

        if opname == "__rpow__" and (
            pa.types.is_floating(pa_dtype)
            or pa.types.is_integer(pa_dtype)
            or pa.types.is_decimal(pa_dtype)
        ):
            mark = pytest.mark.xfail(
                reason=(
                    f"GH#29997: 1**pandas.NA == 1 while 1**pyarrow.NA == NULL "
                    f"for {pa_dtype}"
                )
            )
        elif arrow_temporal_supported and (
            pa.types.is_time(pa_dtype)
            or (
                opname
                in ("__truediv__", "__rtruediv__", "__floordiv__", "__rfloordiv__")
                and pa.types.is_duration(pa_dtype)
            )
        ):
            mark = pytest.mark.xfail(
                raises=TypeError,
                reason=(
                    f"{opname} not supported betweenpd.NA and {pa_dtype} Python scalar"
                ),
            )
        elif opname == "__rfloordiv__" and (
            pa.types.is_integer(pa_dtype) or pa.types.is_decimal(pa_dtype)
        ):
            mark = pytest.mark.xfail(
                raises=pa.ArrowInvalid,
                reason="divide by 0",
            )
        elif opname == "__rtruediv__" and pa.types.is_decimal(pa_dtype):
            mark = pytest.mark.xfail(
                raises=pa.ArrowInvalid,
                reason="divide by 0",
            )

        return mark

    def test_arith_series_with_scalar(self, data, all_arithmetic_operators, request):
        pa_dtype = data.dtype.pyarrow_dtype

        if all_arithmetic_operators == "__rmod__" and pa.types.is_binary(pa_dtype):
            pytest.skip("Skip testing Python string formatting")

        mark = self._get_arith_xfail_marker(all_arithmetic_operators, pa_dtype)
        if mark is not None:
            request.applymarker(mark)

        super().test_arith_series_with_scalar(data, all_arithmetic_operators)

    def test_arith_frame_with_scalar(self, data, all_arithmetic_operators, request):
        pa_dtype = data.dtype.pyarrow_dtype

        if all_arithmetic_operators == "__rmod__" and (
            pa.types.is_string(pa_dtype) or pa.types.is_binary(pa_dtype)
        ):
            pytest.skip("Skip testing Python string formatting")

        mark = self._get_arith_xfail_marker(all_arithmetic_operators, pa_dtype)
        if mark is not None:
            request.applymarker(mark)

        super().test_arith_frame_with_scalar(data, all_arithmetic_operators)

    def test_arith_series_with_array(self, data, all_arithmetic_operators, request):
        pa_dtype = data.dtype.pyarrow_dtype

        if all_arithmetic_operators in (
            "__sub__",
            "__rsub__",
        ) and pa.types.is_unsigned_integer(pa_dtype):
            request.applymarker(
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
            request.applymarker(mark)

        op_name = all_arithmetic_operators
        ser = pd.Series(data)
        # pd.Series([ser.iloc[0]] * len(ser)) may not return ArrowExtensionArray
        # since ser.iloc[0] is a python scalar
        other = pd.Series(pd.array([ser.iloc[0]] * len(ser), dtype=data.dtype))

        self.check_opname(ser, op_name, other)

    def test_add_series_with_extension_array(self, data, request):
        pa_dtype = data.dtype.pyarrow_dtype

        if pa_dtype.equals("int8"):
            request.applymarker(
                pytest.mark.xfail(
                    raises=pa.ArrowInvalid,
                    reason=f"raises on overflow for {pa_dtype}",
                )
            )
        super().test_add_series_with_extension_array(data)

    def test_invalid_other_comp(self, data, comparison_op):
        # GH 48833
        with pytest.raises(
            NotImplementedError, match=".* not implemented for <class 'object'>"
        ):
            comparison_op(data, object())

    @pytest.mark.parametrize("masked_dtype", ["boolean", "Int64", "Float64"])
    def test_comp_masked_numpy(self, masked_dtype, comparison_op):
        # GH 52625
        data = [1, 0, None]
        ser_masked = pd.Series(data, dtype=masked_dtype)
        ser_pa = pd.Series(data, dtype=f"{masked_dtype.lower()}[pyarrow]")
        result = comparison_op(ser_pa, ser_masked)
        if comparison_op in [operator.lt, operator.gt, operator.ne]:
            exp = [False, False, None]
        else:
            exp = [True, True, None]
        expected = pd.Series(exp, dtype=ArrowDtype(pa.bool_()))
        tm.assert_series_equal(result, expected)

    @pytest.mark.filterwarnings(
        "ignore:The default formatting of datetime/timedelta values:DeprecationWarning"
    )
    def test_values_for_json(self, data, request):
        # GH 65127
        # The date32 and date64 dtypes fail already in serialization due to as_unit not
        # implemented for them. Currently the json serialization relies on the default
        # 'epoch' format for datetimes, leading to the filtered Pandas4Warning.
        if data.dtype in [ArrowDtype(pa.date32()), ArrowDtype(pa.date64())]:
            try:
                super().test_values_for_json(data)
            # date32/date64
            except NotImplementedError as err:
                if "as_unit not implemented for date" in str(err):
                    request.applymarker(
                        pytest.mark.xfail(
                            raises=NotImplementedError,
                            reason="as_unit not implemented for date",
                        )
                    )
                raise
        else:
            super().test_values_for_json(data)

    @pytest.mark.filterwarnings(
        "ignore:The default formatting of datetime/timedelta values:DeprecationWarning"
    )
    def test_json_roundtrip(self, data, request):
        # GH 65127
        # All datetime and duration ArrowDtypes with non default resolution of ms fail
        # on roundtrip. The date32 and date64 dtypes fail already in serialization due
        # to as_unit not implemented for them. Currently the json serialization relies
        # on the default 'epoch' format for datetimes, leading to the filtered
        # Pandas4Warning.
        if ((data.dtype.kind in "Mm") and ("ms" not in str(data.dtype))) or (
            data.dtype in [ArrowDtype(pa.date32()), ArrowDtype(pa.date64())]
        ):
            try:
                super().test_json_roundtrip(data)
            except NotImplementedError as err:
                # date32/date64
                if "as_unit not implemented for date" in str(err):
                    request.applymarker(
                        pytest.mark.xfail(
                            raises=NotImplementedError,
                            reason="as_unit not implemented for date",
                        )
                    )
                # timestamp with s unit and US/Pacific or US/Eastern tz
                elif (
                    "Localizing Timestamps which are outside the range of Python"
                ) in str(err):
                    request.applymarker(
                        pytest.mark.xfail(
                            raises=NotImplementedError,
                            reason=(
                                "Localizing Timestamps which are outside the range of "
                                "Python's standard library datetime is not supported"
                            ),
                        )
                    )
                raise
            # all others
            except AssertionError as err:
                if "Series are different" in str(err):
                    request.applymarker(
                        pytest.mark.xfail(
                            raises=AssertionError, reason="Series are different"
                        )
                    )
                raise
        else:
            super().test_json_roundtrip(data)

    def test_plot_on_x_axis(self, plot_data):
        # GH 64535
        # Setup expected exception and message for expected failures
        err_cls = None
        err_msg = ""

        # Set expected values for certain dtypes
        pa_dtype = plot_data["Data"].dtype.pyarrow_dtype

        # Certain dtypes fail with NA values
        if plot_data["Data"].isna().any():
            if pa.types.is_timestamp(pa_dtype) and pa_dtype.tz is not None:
                err_cls = TypeError
                err_msg = re.escape("Failed to convert value(s) to axis units: array(")

        # Call test, errors should only be raised for unsupported dtypes set above
        if err_cls:
            with pytest.raises(err_cls, match=err_msg):
                super().test_plot_on_x_axis(plot_data)
        else:
            super().test_plot_on_x_axis(plot_data)

    def test_plot_on_y_axis(self, plot_data):
        # GH 64535
        # Setup expected exception and message for expected failures
        err_cls = None
        err_msg = ""
        wrn_cls = None
        wrn_msg = ""

        # Set expected values for certain dtypes
        pa_dtype = plot_data["Data"].dtype.pyarrow_dtype

        # str and binary are supported by matplotlib, but not pandas at the moment
        if pa.types.is_string(pa_dtype):
            err_cls = TypeError
            err_msg = "no numeric data to plot"
        elif pa.types.is_binary(pa_dtype):
            err_cls = TypeError
            err_msg = "no numeric data to plot"

        # Certain dtypes fail with NA values
        if plot_data["Data"].isna().any():
            if pa.types.is_time(pa_dtype):
                err_cls = TypeError
                err_msg = re.escape(
                    "float() argument must be a string or a real number, not 'NAType'"
                )
                # tm.assert_produces_warning raises if the warning is not produced even
                # if the test is skipped further down in the call stack, so skip here
                # already to avoid the test failing due to the warning not being raised
                pytest.importorskip("matplotlib")
                wrn_cls = UserWarning
                wrn_msg = "Warning: converting a masked element to nan."
            elif pa.types.is_timestamp(pa_dtype) and pa_dtype.tz is not None:
                err_cls = TypeError
                err_msg = re.escape(
                    "Failed to convert value(s) to axis units: masked_array(data="
                )

        # Call test, errors and warnings should only be raised for unsupported dtypes
        # set above
        if err_cls and wrn_cls:
            with tm.assert_produces_warning(
                wrn_cls, check_stacklevel=False, match=wrn_msg
            ):
                with pytest.raises(err_cls, match=err_msg):
                    super().test_plot_on_y_axis(plot_data)
        elif err_cls:
            with pytest.raises(err_cls, match=err_msg):
                super().test_plot_on_y_axis(plot_data)
        else:
            super().test_plot_on_y_axis(plot_data)
