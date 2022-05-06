"""
Tests for scalar Timedelta arithmetic ops
"""

from __future__ import annotations

from datetime import (
    datetime,
    timedelta,
)
import operator
import re

import numpy as np
import pytest

from pandas._libs.tslibs import OutOfBoundsTimedelta

import pandas as pd
from pandas import (
    NA,
    NaT,
    Timedelta,
    Timestamp,
    offsets,
)
import pandas._testing as tm
from pandas.core import ops


@pytest.fixture(name="tdlike_cls", params=(Timedelta, timedelta, np.timedelta64))
def fixture_tdlike_cls(request) -> type:
    return request.param


@pytest.fixture(
    name="tdlike_or_offset_cls",
    params=(Timedelta, timedelta, np.timedelta64, offsets.Nano),
)
def fixture_tdlike_or_offset_cls(request) -> type:
    return request.param


@pytest.fixture(name="ten_days")
def fixture_ten_days() -> Timedelta:
    return Timedelta(days=10)


@pytest.fixture(name="y2k", params=(Timestamp, np.datetime64, datetime.fromisoformat))
def fixture_y2k(request):
    return request.param("2000-01-01")


@pytest.fixture(name="one_day")
def fixture_one_day(tdlike_cls: type):
    if tdlike_cls is np.timedelta64:
        return np.timedelta64(1, "D")
    return tdlike_cls(days=1)


@pytest.fixture(
    name="na_value",
    params=(None, np.nan, np.float64("NaN"), NaT, NA),
    ids=("None", "np.nan", "np.float64('NaN')", "NaT", "NA"),
)
def fixture_na_value(request):
    return request.param


@pytest.fixture(name="add_op", params=(operator.add, ops.radd))
def fixture_add_op(request):
    return request.param


@pytest.fixture(name="sub_op", params=(operator.sub, ops.rsub))
def fixture_sub_op(request):
    return request.param


@pytest.fixture(
    name="add_or_sub",
    params=(operator.add, ops.radd, operator.sub, ops.rsub),
)
def fixture_add_or_sub(request):
    return request.param


@pytest.fixture(name="mul_op", params=(operator.mul, ops.rmul))
def fixture_mul_op(request):
    return request.param


@pytest.fixture(name="truediv_op", params=(operator.truediv, ops.rtruediv))
def fixture_truediv_op(request):
    return request.param


@pytest.fixture(
    name="floor_mod_divmod_op",
    params=(
        operator.floordiv,
        ops.rfloordiv,
        operator.mod,
        ops.rmod,
        divmod,
        ops.rdivmod,
    ),
)
def fixture_floor_mod_divmod_op(request):
    return request.param


@pytest.fixture(name="td_overflow_msg")
def fixture_td_overflow_msg() -> str:
    return re.escape(
        "outside allowed range [-9223372036854775807ns, 9223372036854775807ns]"
    )


@pytest.fixture(name="invalid_op_msg")
def fixture_invalid_op_msg() -> str:
    messages = (
        "cannot use operands with types",
        "Concatenation operation is not implemented for NumPy arrays",
        "cannot perform",
        "not supported between instances of 'Timedelta' and ",
        re.escape("unsupported operand type(s)"),
    )
    return "|".join(messages)


xfail_type_error = pytest.mark.xfail(
    reason="unsupported",
    raises=TypeError,
    strict=True,
)


def test_binary_ops_not_implemented_for_arbitrary_types(
    ten_days: Timedelta,
    invalid_op_msg: str,
    all_binary_operators,
):
    if all_binary_operators not in (operator.eq, operator.ne):
        with pytest.raises(TypeError, match=invalid_op_msg):
            all_binary_operators(ten_days, object())


class TestAdditionSubtractionScalar:
    """
    Tests for Timedelta.{__add__,__radd__,__sub__,__rsub__} where second operand is a
    scalar.
    """

    @pytest.mark.parametrize(
        "ten_seconds",
        [
            Timedelta(10, unit="s"),
            timedelta(seconds=10),
            np.timedelta64(10, "s"),
            np.timedelta64(10000000000, "ns"),
            offsets.Second(10),
        ],
    )
    def test_add_sub_ten_seconds(self, ten_seconds):
        # GH#6808
        base = Timestamp("20130101 09:01:12.123456")
        expected_add = Timestamp("20130101 09:01:22.123456")
        expected_sub = Timestamp("20130101 09:01:02.123456")

        result = base + ten_seconds
        assert result == expected_add

        result = base - ten_seconds
        assert result == expected_sub

    @pytest.mark.parametrize(
        "one_day_ten_secs",
        [
            Timedelta("1 day, 00:00:10"),
            Timedelta("1 days, 00:00:10"),
            timedelta(days=1, seconds=10),
            np.timedelta64(1, "D") + np.timedelta64(10, "s"),
            offsets.Day() + offsets.Second(10),
        ],
    )
    def test_add_sub_one_day_ten_seconds(self, one_day_ten_secs):
        # GH#6808
        base = Timestamp("20130102 09:01:12.123456")
        expected_add = Timestamp("20130103 09:01:22.123456")
        expected_sub = Timestamp("20130101 09:01:02.123456")

        result = base + one_day_ten_secs
        assert result == expected_add

        result = base - one_day_ten_secs
        assert result == expected_sub

    @pytest.mark.parametrize("value", (2, 2.0), ids=("int", "float"))
    def test_add_or_sub_numeric_raises(
        self,
        ten_days: Timedelta,
        add_or_sub,
        value,
        invalid_op_msg: str,
    ):
        with pytest.raises(TypeError, match=invalid_op_msg):
            add_or_sub(ten_days, value)

    def test_add_datetimelike(self, ten_days: Timedelta, add_op, y2k):
        # GH#19738
        result = add_op(y2k, ten_days)
        expected = Timestamp("2000-01-11")

        if type(y2k) != datetime and add_op != ops.radd:
            # datetime + Timedelta does _not_ call Timedelta.__radd__,
            # so we get a datetime back instead of a Timestamp
            assert isinstance(result, Timestamp)
        assert result == expected

    def test_sub_datetimelike(self, ten_days: Timedelta, y2k, invalid_op_msg: str):
        assert y2k - ten_days == Timestamp("1999-12-22")

        with pytest.raises(TypeError, match=invalid_op_msg):
            ten_days - y2k

    def test_add_timedeltalike(self, ten_days: Timedelta, add_op, one_day):
        result = add_op(ten_days, one_day)
        expected = Timedelta(days=11)
        assert isinstance(result, Timedelta)
        assert result == expected

    def test_sub_timedeltalike(self, ten_days: Timedelta, one_day, sub_op):
        result = sub_op(ten_days, one_day)
        expected = Timedelta(days=9) if sub_op is operator.sub else Timedelta(days=-9)
        assert isinstance(result, Timedelta)
        assert result == expected

    def test_add_offset(self, ten_days: Timedelta, add_op):
        result = add_op(ten_days, offsets.Hour(6))
        expected = Timedelta(days=10, hours=6)
        assert isinstance(result, Timedelta)
        assert result == expected

    def test_sub_offset(self, ten_days: Timedelta, sub_op):
        result = sub_op(ten_days, offsets.Hour(1))
        if sub_op is operator.sub:
            expected = Timedelta(hours=239)
        else:
            expected = Timedelta(hours=-239)

        assert isinstance(result, Timedelta)
        assert result == expected

    def test_with_timedeltadlike_raises_for_any_result_above_td_max(
        self,
        tdlike_or_offset_cls,
        td_overflow_msg: str,
    ):
        with pytest.raises(OutOfBoundsTimedelta, match=td_overflow_msg):
            Timedelta.max + tdlike_or_offset_cls(1)

        with pytest.raises(OutOfBoundsTimedelta, match=td_overflow_msg):
            Timedelta.max - (tdlike_or_offset_cls(-1))

    def test_no_error_for_result_1ns_below_td_min(self):
        assert Timedelta.min + Timedelta(-1, "ns") is NaT
        assert offsets.Nano(-1) + Timedelta.min is NaT
        assert Timedelta.min - np.timedelta64(1, "ns") is NaT

    def test_raises_for_any_result_2ns_below_td_min(
        self,
        tdlike_or_offset_cls: type,
        td_overflow_msg: str,
    ):
        with pytest.raises(OutOfBoundsTimedelta, match=td_overflow_msg):
            Timedelta.min + tdlike_or_offset_cls(-2)

        with pytest.raises(OutOfBoundsTimedelta, match=td_overflow_msg):
            Timedelta.min - tdlike_or_offset_cls(2)

    def test_add_or_sub_na(self, request, ten_days: Timedelta, add_or_sub, na_value):
        if na_value is NA:
            request.applymarker(xfail_type_error)
        result = add_or_sub(ten_days, na_value)
        assert result is NaT


class TestAdditionSubtractionBox:
    """
    Tests for Timedelta.{__add__,__radd__,__sub__,__rsub__} where second operand is a
    Array/Index/Series/DataFrame.
    """

    @pytest.mark.parametrize("value", (2, 2.0), ids=("int", "float"))
    def test_add_or_sub_numeric_raises(
        self,
        ten_days: Timedelta,
        add_or_sub,
        box_with_array,
        value,
        invalid_op_msg: str,
    ):
        other = tm.box_expected([value], box_with_array)
        with pytest.raises(TypeError, match=invalid_op_msg):
            add_or_sub(ten_days, other)

    def test_add_datetimelike(self):
        pass

    def test_sub_from_datetimelike(self, ten_days: Timedelta, box_with_array):
        # GH#21980
        other = tm.box_expected([np.datetime64("2000-01-11")], box_with_array)
        expected = tm.box_expected([np.datetime64("2000-01-01")], box_with_array)
        result = other - ten_days
        tm.assert_equal(result, expected)

    def test_sub_mixed_most_timedeltalike_object_dtype_array(self):
        # GH#21980
        now = Timestamp("2021-11-09 09:54:00")
        arr = np.array([now, Timedelta("1D"), np.timedelta64(2, "h")])
        exp = np.array(
            [
                now - Timedelta("1D"),
                Timedelta("0D"),
                np.timedelta64(2, "h") - Timedelta("1D"),
            ]
        )
        res = arr - Timedelta("1D")
        tm.assert_numpy_array_equal(res, exp)

    def test_rsub_mixed_most_timedeltalike_object_dtype_array(self, invalid_op_msg):
        # GH#21980
        now = Timestamp("2021-11-09 09:54:00")
        arr = np.array([now, Timedelta("1D"), np.timedelta64(2, "h")])
        with pytest.raises(TypeError, match=invalid_op_msg):
            Timedelta("1D") - arr

    def test_add_timedeltalike_object_dtype_array(self, add_op):
        # GH#21980
        arr = np.array([Timestamp("20130101 9:01"), Timestamp("20121230 9:02")])
        exp = np.array([Timestamp("20130102 9:01"), Timestamp("20121231 9:02")])
        res = add_op(arr, Timedelta("1D"))
        tm.assert_numpy_array_equal(res, exp)

    def test_add_mixed_timedeltalike_object_dtype_array(self, add_op):
        # GH#21980
        now = Timestamp("2021-11-09 09:54:00")
        arr = np.array([now, Timedelta("1D")])
        exp = np.array([now + Timedelta("1D"), Timedelta("2D")])
        res = add_op(arr, Timedelta("1D"))
        tm.assert_numpy_array_equal(res, exp)

    def test_add_td64_ndarray(self, ten_days: Timedelta, add_op):
        result = add_op(ten_days, np.array([np.timedelta64(1, "D")]))
        expected = np.array([Timedelta(days=11).to_timedelta64()])
        tm.assert_numpy_array_equal(result, expected)

    def test_sub_td64_ndarray(self, ten_days: Timedelta, sub_op):
        result = sub_op(ten_days, np.array([np.timedelta64(10, "D")]))
        expected = np.array([0], dtype="timedelta64[ns]")
        tm.assert_numpy_array_equal(result, expected)

    def test_add_sub_dt64_ndarray(self):
        td = Timedelta("1 day")
        other = pd.to_datetime(["2000-01-01"]).values

        expected = pd.to_datetime(["2000-01-02"]).values
        tm.assert_numpy_array_equal(td + other, expected)
        tm.assert_numpy_array_equal(other + td, expected)

        expected = pd.to_datetime(["1999-12-31"]).values
        tm.assert_numpy_array_equal(-td + other, expected)
        tm.assert_numpy_array_equal(other - td, expected)

    def test_na(self):
        pass


class TestMultiplicationScalar:
    """
    Tests for Timedelta.{__mul__,__rmul__} where second operand is a scalar.
    """

    @pytest.mark.parametrize(
        "factor,expected",
        ((2, 20), (1.5, 15), (-1, -10), (-1, -10)),
    )
    def test_numeric(self, ten_days: Timedelta, mul_op, factor, expected):
        # GH#19738
        result = mul_op(ten_days, factor)
        assert result == Timedelta(expected, "D")
        assert isinstance(result, Timedelta)

    @pytest.mark.parametrize("value", (Timestamp("2020-01-02"), Timedelta(1)))
    def test_datetimelike_or_timedeltalike_raises(
        self,
        ten_days: Timedelta,
        mul_op,
        value,
        invalid_op_msg: str,
    ):
        # timedelta * datetime is gibberish, as is multiplying by another timedelta
        with pytest.raises(TypeError, match=invalid_op_msg):
            mul_op(ten_days, value)

    def test_offset_raises(self):
        pass

    @pytest.mark.parametrize("value", (Timedelta.min, Timedelta.max))
    def test_raises_for_overflow(self, mul_op, td_overflow_msg: str, value: Timedelta):
        with pytest.raises(OutOfBoundsTimedelta, match=td_overflow_msg):
            mul_op(value, 2)

        with pytest.raises(OutOfBoundsTimedelta, match=td_overflow_msg):
            mul_op(value, 1.1)

    def test_na(self, request, ten_days: Timedelta, mul_op, na_value):
        if na_value is None or na_value is NaT or na_value is NA:
            request.applymarker(xfail_type_error)
        result = mul_op(ten_days, na_value)
        assert result is NaT


class TestMultiplicationBox:
    """
    Tests for Timedelta.{__mul__,__rmul__} where second operand is a
    Array/Index/Series/DataFrame.
    """

    @pytest.mark.parametrize("factor,expected", ((2, 20), (1.5, 15)))
    def test_numeric(self, ten_days, mul_op, factor, expected, box_with_array):
        other = tm.box_expected([factor], box_with_array)
        expected = tm.box_expected(
            [Timedelta(expected, "D").to_timedelta64()],
            box_with_array,
        )
        result = mul_op(ten_days, other)
        tm.assert_equal(result, expected)

    @pytest.mark.parametrize("value", (Timestamp.min, Timedelta.max))
    def test_datetimelike_or_timedeltalike_raises(
        self,
        ten_days: Timedelta,
        mul_op,
        value,
        box_with_array,
        invalid_op_msg: str,
    ):
        box = tm.box_expected([value], box_with_array)
        with pytest.raises(TypeError, match=invalid_op_msg):
            mul_op(ten_days, box)

    def test_offset_raises(self):
        pass

    def test_raises_for_overflow(self):
        pass

    def test_na(self):
        pass


class TestTrueDivisionScalar:
    """
    Tests for Timedelta.{__truediv__,__rtruediv__} where second operand is a scalar.
    """

    def test_truediv_numeric(self, ten_days: Timedelta, any_real_numpy_dtype):
        # GH#19738
        scalar = np.dtype(any_real_numpy_dtype).type(2.0)
        result = ten_days / scalar
        assert isinstance(result, Timedelta)
        assert result == Timedelta(days=5)

    def test_rtruediv_numeric_raises(
        self,
        ten_days: Timedelta,
        invalid_op_msg: str,
        any_real_numpy_dtype,
    ):
        scalar = np.dtype(any_real_numpy_dtype).type(2.0)
        with pytest.raises(TypeError, match=invalid_op_msg):
            scalar / ten_days

    def test_datetimelike_raises(
        self,
        ten_days: Timedelta,
        truediv_op,
        y2k,
        invalid_op_msg: str,
    ):
        with pytest.raises(TypeError, match=invalid_op_msg):
            truediv_op(ten_days, y2k)

    def test_timedeltalike(self, ten_days: Timedelta, one_day):
        # GH#19738
        assert ten_days / one_day == 10
        assert one_day / ten_days == 0.1

    def test_offset(self, ten_days: Timedelta):
        assert ten_days / offsets.Hour(12) == 20
        assert offsets.Hour(12) / ten_days == 0.05

    def test_na(self, request, ten_days: Timedelta, truediv_op, na_value):
        expected = NaT
        if na_value is NA or (
            truediv_op is ops.rtruediv and isinstance(na_value, float)
        ):
            request.applymarker(xfail_type_error)
        elif na_value is None or na_value is NaT:
            expected = np.nan
        result = truediv_op(ten_days, na_value)
        assert result is expected


class TestTrueDivisionBox:
    """
    Tests for Timedelta.{__floordiv__,__rfloordiv__,__truediv__,__rtruediv__} where
    second operand is a Array/Index/Series/DataFrame.
    """

    def test_truediv_numeric(self, ten_days: Timedelta, any_real_numpy_dtype):
        # GH#19738
        scalar = np.dtype(any_real_numpy_dtype).type(2.0)
        result = ten_days / scalar
        assert isinstance(result, Timedelta)
        assert result == Timedelta(days=5)

    def test_rtruediv_numeric_raises(
        self,
        ten_days: Timedelta,
        invalid_op_msg: str,
        any_real_numpy_dtype,
    ):
        scalar = np.dtype(any_real_numpy_dtype).type(2.0)
        with pytest.raises(TypeError, match=invalid_op_msg):
            scalar / ten_days

    def test_datetimelike_raises(
        self,
        ten_days: Timedelta,
        truediv_op,
        y2k,
        box_with_array,
        invalid_op_msg: str,
    ):
        other = tm.box_expected((y2k,), box_with_array)
        with pytest.raises(TypeError, match=invalid_op_msg):
            truediv_op(ten_days, other)

    def test_timedeltalike(
        self,
        ten_days: Timedelta,
        truediv_op,
        tdlike_cls,
        box_with_array,
    ):
        # TODO:
        elem = tdlike_cls(days=10) if tdlike_cls is timedelta else tdlike_cls(10, "D")
        other = tm.box_expected((elem,), box_with_array)

        if box_with_array is pd.array:
            expected = np.array((1.0,))
        else:
            expected = tm.box_expected((1.0,), box_with_array)

        result = truediv_op(ten_days, other)
        tm.assert_equal(result, expected)

    def test_offset(self, ten_days: Timedelta):
        ...

    def test_na(self, request, ten_days: Timedelta, truediv_op, na_value):
        ...


class TestFloorModuloDivisionScalar:
    """
    Timedelta.{__floordiv__,__rfloordiv__,__mod__,__rmod__,__divmod__,__rdivmod__} tests
    where second operand is a scalar.
    """

    def test_floordiv_numeric(self):
        pass

    def test_rfloordiv_numeric(
        self,
        ten_days: Timedelta,
        any_real_numpy_dtype,
        invalid_op_msg: str,
    ):
        # int32 deprecated GH#19761, enforced GH#29797
        scalar = np.dtype(any_real_numpy_dtype).type(1.0)
        assert ten_days.__rfloordiv__(scalar) is NotImplemented
        with pytest.raises(TypeError, match=invalid_op_msg):
            scalar // ten_days

    def test_mod_numeric(self):
        # GH#19365
        td = Timedelta(hours=37)

        # Numeric Others
        result = td % 2
        assert isinstance(result, Timedelta)
        assert result == Timedelta(0)

        result = td % 1e12
        assert isinstance(result, Timedelta)
        assert result == Timedelta(minutes=3, seconds=20)

        result = td % int(1e12)
        assert isinstance(result, Timedelta)
        assert result == Timedelta(minutes=3, seconds=20)

    def test_rmod_numeric(self):
        # GH#19365
        td = Timedelta(minutes=3)

        msg = "unsupported operand"
        with pytest.raises(TypeError, match=msg):
            Timestamp("2018-01-22") % td

        with pytest.raises(TypeError, match=msg):
            15 % td

        with pytest.raises(TypeError, match=msg):
            16.0 % td

        msg = "Invalid dtype int"
        with pytest.raises(TypeError, match=msg):
            np.array([22, 24]) % td

    def test_divmod_numeric(self):
        # GH#19365
        td = Timedelta(days=2, hours=6)

        result = divmod(td, 53 * 3600 * 1e9)
        assert result[0] == Timedelta(1, unit="ns")
        assert isinstance(result[1], Timedelta)
        assert result[1] == Timedelta(hours=1)

        assert result
        result = divmod(td, np.nan)
        assert result[0] is NaT
        assert result[1] is NaT

    def test_rdivmod_numeric(self):
        pass

    def test_datetimelike_raises(
        self,
        ten_days: Timedelta,
        floor_mod_divmod_op,
        y2k,
        invalid_op_msg: str,
    ):
        # GH#18846
        with pytest.raises(TypeError, match=invalid_op_msg):
            floor_mod_divmod_op(ten_days, y2k)

    def test_floordiv_timedeltalike(self):
        pass

    def test_rfloordiv_timedeltalike(self):
        # GH#18846
        td = Timedelta(hours=3, minutes=3)
        scalar = Timedelta(hours=3, minutes=4)

        # scalar others
        # x // Timedelta is defined only for timedelta-like x. int-like,
        # float-like, and date-like, in particular, should all either
        # a) raise TypeError directly or
        # b) return NotImplemented, following which the reversed
        #    operation will raise TypeError.
        assert td.__rfloordiv__(scalar) == 1
        assert (-td).__rfloordiv__(scalar.to_pytimedelta()) == -2
        assert (2 * td).__rfloordiv__(scalar.to_timedelta64()) == 0

    def test_mod_timedeltalike(self):
        # GH#19365
        td = Timedelta(hours=37)

        # Timedelta-like others
        result = td % Timedelta(hours=6)
        assert isinstance(result, Timedelta)
        assert result == Timedelta(hours=1)

        result = td % timedelta(minutes=60)
        assert isinstance(result, Timedelta)
        assert result == Timedelta(0)

        result = td % NaT
        assert result is NaT

    def test_mod_timedelta64(self):
        # GH#19365
        td = Timedelta(hours=37)

        result = td % np.timedelta64(2, "h")
        assert isinstance(result, Timedelta)
        assert result == Timedelta(hours=1)

    def test_rmod_timedeltalike(self):
        # GH#19365
        td = Timedelta(minutes=3)

        result = timedelta(minutes=4) % td
        assert isinstance(result, Timedelta)
        assert result == Timedelta(minutes=1)

    def test_rmod_timedelta64(self):
        # GH#19365
        td = Timedelta(minutes=3)
        result = np.timedelta64(5, "m") % td
        assert isinstance(result, Timedelta)
        assert result == Timedelta(minutes=2)

    def test_divmod_timedeltalike(self):
        # GH#19365
        td = Timedelta(days=2, hours=6)

        result = divmod(td, timedelta(days=1))
        assert result[0] == 2
        assert isinstance(result[1], Timedelta)
        assert result[1] == Timedelta(hours=6)

        result = divmod(td, 54)
        assert result[0] == Timedelta(hours=1)
        assert isinstance(result[1], Timedelta)
        assert result[1] == Timedelta(0)

        result = divmod(td, NaT)
        assert np.isnan(result[0])
        assert result[1] is NaT

    def test_rdivmod_pytimedelta(self):
        # GH#19365
        result = divmod(timedelta(days=2, hours=6), Timedelta(days=1))
        assert result[0] == 2
        assert isinstance(result[1], Timedelta)
        assert result[1] == Timedelta(hours=6)

    def test_floordiv_offsets(self):
        # GH#19738
        td = Timedelta(hours=3, minutes=4)
        assert td // offsets.Hour(1) == 3

        assert td // offsets.Minute(2) == 92

    def test_rfloordiv_offsets(self):
        # GH#19738
        assert offsets.Hour(1) // Timedelta(minutes=25) == 2

    def test_mod_offset(self):
        # GH#19365
        td = Timedelta(hours=37)

        result = td % offsets.Hour(5)
        assert isinstance(result, Timedelta)
        assert result == Timedelta(hours=2)

    def test_rmod_offset(self):
        pass

    def test_divmod_offset(self):
        # GH#19365
        td = Timedelta(days=2, hours=6)

        result = divmod(td, offsets.Hour(-4))
        assert result[0] == -14
        assert isinstance(result[1], Timedelta)
        assert result[1] == Timedelta(hours=-2)

    def test_rdivmod_offset(self):
        result = divmod(offsets.Hour(54), Timedelta(hours=-4))
        assert result[0] == -14
        assert isinstance(result[1], Timedelta)
        assert result[1] == Timedelta(hours=-2)

    def test_floordiv_na(self, request, ten_days: Timedelta, na_value):
        expected = NaT
        if na_value is NA:
            request.applymarker(xfail_type_error)
        elif na_value is None or na_value is NaT:
            expected = np.nan

        result = ten_days // na_value
        assert result is expected

    def test_rfloordiv_na(self, request, ten_days: Timedelta, na_value):
        expected = np.nan
        if na_value is NA or isinstance(na_value, float):
            request.applymarker(xfail_type_error)

        result = na_value // ten_days
        assert result is expected

    def test_mod_na(self, request, ten_days: Timedelta, na_value):
        expected = NaT
        if na_value is None or na_value is NA:
            request.applymarker(xfail_type_error)

        result = ten_days % na_value
        assert result is expected

    def test_rmod_na(self, request, ten_days: Timedelta, na_value):
        if na_value is not NaT:
            request.applymarker(xfail_type_error)

        result = na_value % ten_days
        assert result is NaT

    def test_divmod_na(self, request, ten_days: Timedelta, na_value):
        expected = (NaT, NaT)
        if na_value is None or na_value is NA:
            request.applymarker(xfail_type_error)
        elif na_value is NaT:
            expected = (np.nan, NaT)

        result = divmod(ten_days, na_value)
        assert result == expected

    def test_rdivmod_na(self, request, ten_days: Timedelta, na_value):
        expected = (np.nan, NaT)
        if na_value is not NaT:
            request.applymarker(xfail_type_error)

        result = ops.rdivmod(ten_days, na_value)
        assert result == expected


class TestFloorModuloDivisionBox:
    """
    Timedelta.{__floordiv__,__rfloordiv__,__mod__,__rmod__,__divmod__, __rdivmod__}
    tests where second operand is a Array/Index/Series/DataFrame.
    """

    def test_floordiv_numeric_series(self):
        # GH#18846
        td = Timedelta(hours=3, minutes=4)
        ser = pd.Series([1], dtype=np.int64)
        res = td // ser
        assert res.dtype.kind == "m"

    def test_rfloordiv_numeric_series(self):
        # GH#18846
        td = Timedelta(hours=3, minutes=3)
        ser = pd.Series([1], dtype=np.int64)
        res = td.__rfloordiv__(ser)
        assert res is NotImplemented

        msg = "Invalid dtype"
        with pytest.raises(TypeError, match=msg):
            # Deprecated GH#19761, enforced GH#29797
            ser // td

    def test_rfloordiv_intarray(self):
        # deprecated GH#19761, enforced GH#29797
        ints = np.array([1349654400, 1349740800, 1349827200, 1349913600]) * 10**9

        msg = "Invalid dtype"
        with pytest.raises(TypeError, match=msg):
            ints // Timedelta(1, unit="s")

    def test_floordiv_timedeltalike_array(self):
        # GH#18846
        td = Timedelta(hours=3, minutes=4)
        scalar = Timedelta(hours=3, minutes=3)

        # Array-like others
        assert td // np.array(scalar.to_timedelta64()) == 1

        res = (3 * td) // np.array([scalar.to_timedelta64()])
        expected = np.array([3], dtype=np.int64)
        tm.assert_numpy_array_equal(res, expected)

        res = (10 * td) // np.array([scalar.to_timedelta64(), np.timedelta64("NaT")])
        expected = np.array([10, np.nan])
        tm.assert_numpy_array_equal(res, expected)

    def test_rfloordiv_timedeltalike_array(self):
        # GH#18846
        td = Timedelta(hours=3, minutes=3)
        scalar = Timedelta(hours=3, minutes=4)

        # Array-like others
        assert td.__rfloordiv__(np.array(scalar.to_timedelta64())) == 1

        res = td.__rfloordiv__(np.array([(3 * scalar).to_timedelta64()]))
        expected = np.array([3], dtype=np.int64)
        tm.assert_numpy_array_equal(res, expected)

        arr = np.array([(10 * scalar).to_timedelta64(), np.timedelta64("NaT")])
        res = td.__rfloordiv__(arr)
        expected = np.array([10, np.nan])
        tm.assert_numpy_array_equal(res, expected)

    def test_na(self):
        ...


class TestComparison:
    def test_compare_tick(self, tick_classes):
        cls = tick_classes

        off = cls(4)
        td = off.delta
        assert isinstance(td, Timedelta)

        assert td == off
        assert not td != off
        assert td <= off
        assert td >= off
        assert not td < off
        assert not td > off

        assert not td == 2 * off
        assert td != 2 * off
        assert td <= 2 * off
        assert td < 2 * off
        assert not td >= 2 * off
        assert not td > 2 * off

    def test_comparison_object_array(self):
        # analogous to GH#15183
        td = Timedelta("2 days")
        other = Timedelta("3 hours")

        arr = np.array([other, td], dtype=object)
        res = arr == td
        expected = np.array([False, True], dtype=bool)
        assert (res == expected).all()

        # 2D case
        arr = np.array([[other, td], [td, other]], dtype=object)
        res = arr != td
        expected = np.array([[True, False], [False, True]], dtype=bool)
        assert res.shape == expected.shape
        assert (res == expected).all()

    def test_compare_timedelta_ndarray(self):
        # GH#11835
        periods = [Timedelta("0 days 01:00:00"), Timedelta("0 days 01:00:00")]
        arr = np.array(periods)
        result = arr[0] > arr
        expected = np.array([False, False])
        tm.assert_numpy_array_equal(result, expected)

    def test_compare_td64_ndarray(self):
        # GG#33441
        arr = np.arange(5).astype("timedelta64[ns]")
        td = Timedelta(arr[1])

        expected = np.array([False, True, False, False, False], dtype=bool)

        result = td == arr
        tm.assert_numpy_array_equal(result, expected)

        result = arr == td
        tm.assert_numpy_array_equal(result, expected)

        result = td != arr
        tm.assert_numpy_array_equal(result, ~expected)

        result = arr != td
        tm.assert_numpy_array_equal(result, ~expected)

    def test_compare_custom_object(self):
        """
        Make sure non supported operations on Timedelta returns NonImplemented
        and yields to other operand (GH#20829).
        """

        class CustomClass:
            def __init__(self, cmp_result=None) -> None:
                self.cmp_result = cmp_result

            def generic_result(self):
                if self.cmp_result is None:
                    return NotImplemented
                else:
                    return self.cmp_result

            def __eq__(self, other):
                return self.generic_result()

            def __gt__(self, other):
                return self.generic_result()

        t = Timedelta("1s")

        assert not (t == "string")
        assert not (t == 1)
        assert not (t == CustomClass())
        assert not (t == CustomClass(cmp_result=False))

        assert t < CustomClass(cmp_result=True)
        assert not (t < CustomClass(cmp_result=False))

        assert t == CustomClass(cmp_result=True)

    @pytest.mark.parametrize("val", ["string", 1])
    def test_compare_unknown_type(self, val):
        # GH#20829
        t = Timedelta("1s")
        msg = "not supported between instances of 'Timedelta' and '(int|str)'"
        with pytest.raises(TypeError, match=msg):
            t >= val
        with pytest.raises(TypeError, match=msg):
            t > val
        with pytest.raises(TypeError, match=msg):
            t <= val
        with pytest.raises(TypeError, match=msg):
            t < val

    def test_na(self):
        pass
