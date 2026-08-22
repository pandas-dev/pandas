"""
Tests for scalar Timedelta arithmetic ops
"""

from datetime import (
    datetime,
    timedelta,
)
import operator

import numpy as np
import pytest

from pandas.compat import PY313
from pandas.errors import (
    OutOfBoundsDatetime,
    OutOfBoundsTimedelta,
    Pandas4Warning,
)

import pandas as pd
from pandas import (
    NaT,
    Timedelta,
    Timestamp,
    offsets,
)
import pandas._testing as tm
from pandas.core import ops


class TestTimedeltaAdditionSubtraction:
    """
    Tests for Timedelta methods:

        __add__, __radd__,
        __sub__, __rsub__
    """

    def test_td_add_sub_pydatetime(self, unit):
        # GH#53643
        td = Timedelta(hours=23).as_unit(unit)
        dt = datetime(2016, 1, 1)

        expected = datetime(2016, 1, 1, 23)
        result = dt + td
        assert result == expected
        result = td + dt
        assert result == expected

        expected = datetime(2015, 12, 31, 1)
        result = dt - td
        assert result == expected

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
    def test_td_add_sub_ten_seconds(self, ten_seconds):
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
    def test_td_add_sub_one_day_ten_seconds(self, one_day_ten_secs):
        # GH#6808
        base = Timestamp("20130102 09:01:12.123456")
        expected_add = Timestamp("20130103 09:01:22.123456")
        expected_sub = Timestamp("20130101 09:01:02.123456")

        result = base + one_day_ten_secs
        assert result == expected_add

        result = base - one_day_ten_secs
        assert result == expected_sub

    @pytest.mark.parametrize("op", [operator.add, ops.radd])
    def test_td_add_datetimelike_scalar(self, op):
        # GH#19738
        td = Timedelta(10, unit="D")

        result = op(td, datetime(2016, 1, 1))
        if op is operator.add:
            # datetime + Timedelta does _not_ call Timedelta.__radd__,
            # so we get a datetime back instead of a Timestamp
            assert isinstance(result, Timestamp)
        assert result == Timestamp(2016, 1, 11)

        result = op(td, Timestamp("2018-01-12 18:09"))
        assert isinstance(result, Timestamp)
        assert result == Timestamp("2018-01-22 18:09")

        result = op(td, np.datetime64("2018-01-12"))
        assert isinstance(result, Timestamp)
        assert result == Timestamp("2018-01-22")

        result = op(td, NaT)
        assert result is NaT

    def test_td_add_timestamp_overflow(self):
        ts = Timestamp("1700-01-01").as_unit("ns")
        msg = "Cannot cast 259987 days 00:00:00 to unit='ns' without overflow."
        with pytest.raises(OutOfBoundsTimedelta, match=msg):
            ts + Timedelta(13 * 19999, unit="D")

        msg = "Cannot cast 259987 days 00:00:00 to unit='ns' without overflow"
        with pytest.raises(OutOfBoundsTimedelta, match=msg):
            ts + timedelta(days=13 * 19999)

    @pytest.mark.parametrize("op", [operator.add, ops.radd])
    def test_td_add_td(self, op):
        td = Timedelta(10, unit="D")

        result = op(td, Timedelta(days=10))
        assert isinstance(result, Timedelta)
        assert result == Timedelta(days=20)

    @pytest.mark.parametrize("op", [operator.add, ops.radd])
    def test_td_add_pytimedelta(self, op):
        td = Timedelta(10, unit="D")
        result = op(td, timedelta(days=9))
        assert isinstance(result, Timedelta)
        assert result == Timedelta(days=19)

    @pytest.mark.parametrize("op", [operator.add, ops.radd])
    def test_td_add_timedelta64(self, op):
        td = Timedelta(10, unit="D")
        result = op(td, np.timedelta64(-4, "D"))
        assert isinstance(result, Timedelta)
        assert result == Timedelta(days=6)

    @pytest.mark.parametrize("op", [operator.add, ops.radd])
    def test_td_add_offset(self, op):
        td = Timedelta(10, unit="D")

        result = op(td, offsets.Hour(6))
        assert isinstance(result, Timedelta)
        assert result == Timedelta(days=10, hours=6)

    def test_td_sub_td(self):
        td = Timedelta(10, unit="D")
        expected = Timedelta(0, unit="ns")
        result = td - td
        assert isinstance(result, Timedelta)
        assert result == expected

    def test_td_sub_pytimedelta(self):
        td = Timedelta(10, unit="D")
        expected = Timedelta(0, unit="ns")

        result = td - td.to_pytimedelta()
        assert isinstance(result, Timedelta)
        assert result == expected

        result = td.to_pytimedelta() - td
        assert isinstance(result, Timedelta)
        assert result == expected

    def test_td_sub_timedelta64(self):
        td = Timedelta(10, unit="D")
        expected = Timedelta(0, unit="ns")

        result = td - td.to_timedelta64()
        assert isinstance(result, Timedelta)
        assert result == expected

        result = td.to_timedelta64() - td
        assert isinstance(result, Timedelta)
        assert result == expected

    def test_td_sub_nat(self):
        # In this context pd.NaT is treated as timedelta-like
        td = Timedelta(10, unit="D")
        result = td - NaT
        assert result is NaT

    def test_td_sub_td64_nat(self):
        td = Timedelta(10, unit="D")
        td_nat = np.timedelta64("NaT", "ns")

        result = td - td_nat
        assert result is NaT

        result = td_nat - td
        assert result is NaT

    def test_td_sub_offset(self):
        td = Timedelta(10, unit="D")
        result = td - offsets.Hour(1)
        assert isinstance(result, Timedelta)
        assert result == Timedelta(239, unit="h")

    def test_td_add_sub_numeric_raises(self):
        td = Timedelta(10, unit="D")
        msg = "unsupported operand type"
        for other in [2, 2.0, np.int64(2), np.float64(2)]:
            with pytest.raises(TypeError, match=msg):
                td + other
            with pytest.raises(TypeError, match=msg):
                other + td
            with pytest.raises(TypeError, match=msg):
                td - other
            with pytest.raises(TypeError, match=msg):
                other - td

    def test_td_add_sub_int_ndarray(self):
        td = Timedelta("1 day")
        other = np.array([1])

        msg = r"unsupported operand type\(s\) for \+: 'Timedelta' and 'int'"
        with pytest.raises(TypeError, match=msg):
            td + np.array([1])

        msg = "|".join(
            [
                (
                    r"unsupported operand type\(s\) for \+: 'numpy.ndarray' "
                    "and 'Timedelta'"
                ),
                # This message goes on to say "Please do not rely on this error;
                #  it may not be given on all Python implementations"
                "Concatenation operation is not implemented for NumPy arrays",
            ]
        )
        with pytest.raises(TypeError, match=msg):
            other + td
        msg = r"unsupported operand type\(s\) for -: 'Timedelta' and 'int'"
        with pytest.raises(TypeError, match=msg):
            td - other
        msg = r"unsupported operand type\(s\) for -: 'numpy.ndarray' and 'Timedelta'"
        with pytest.raises(TypeError, match=msg):
            other - td

    def test_td_rsub_nat(self):
        td = Timedelta(10, unit="D")
        result = NaT - td
        assert result is NaT

        result = np.datetime64("NaT", "ns") - td
        assert result is NaT

    def test_td_rsub_offset(self):
        result = offsets.Hour(1) - Timedelta(10, unit="D")
        assert isinstance(result, Timedelta)
        assert result == Timedelta(-239, unit="h")

    def test_td_sub_timedeltalike_object_dtype_array(self):
        # GH#21980
        arr = np.array([Timestamp("20130101 9:01"), Timestamp("20121230 9:02")])
        exp = np.array([Timestamp("20121231 9:01"), Timestamp("20121229 9:02")])
        res = arr - Timedelta("1D")
        tm.assert_numpy_array_equal(res, exp)

    def test_td_sub_mixed_most_timedeltalike_object_dtype_array(self):
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

    def test_td_rsub_mixed_most_timedeltalike_object_dtype_array(self):
        # GH#21980
        now = Timestamp("2021-11-09 09:54:00")
        arr = np.array([now, Timedelta("1D"), np.timedelta64(2, "h")])
        msg = r"unsupported operand type\(s\) for \-: 'Timedelta' and 'Timestamp'"
        with pytest.raises(TypeError, match=msg):
            Timedelta("1D") - arr

    @pytest.mark.parametrize("op", [operator.add, ops.radd])
    def test_td_add_timedeltalike_object_dtype_array(self, op):
        # GH#21980
        arr = np.array([Timestamp("20130101 9:01"), Timestamp("20121230 9:02")])
        exp = np.array([Timestamp("20130102 9:01"), Timestamp("20121231 9:02")])
        res = op(arr, Timedelta("1D"))
        tm.assert_numpy_array_equal(res, exp)

    @pytest.mark.parametrize("op", [operator.add, ops.radd])
    def test_td_add_mixed_timedeltalike_object_dtype_array(self, op):
        # GH#21980
        now = Timestamp("2021-11-09 09:54:00")
        arr = np.array([now, Timedelta("1D")])
        exp = np.array([now + Timedelta("1D"), Timedelta("2D")])
        res = op(arr, Timedelta("1D"))
        tm.assert_numpy_array_equal(res, exp)

    def test_td_add_sub_td64_ndarray(self):
        td = Timedelta("1 day")

        other = np.array([td.to_timedelta64()])
        expected = np.array([Timedelta("2 Days").to_timedelta64()])

        result = td + other
        tm.assert_numpy_array_equal(result, expected)
        result = other + td
        tm.assert_numpy_array_equal(result, expected)

        result = td - other
        tm.assert_numpy_array_equal(result, expected * 0)
        result = other - td
        tm.assert_numpy_array_equal(result, expected * 0)

    def test_td_add_sub_td64_ndarray_lands_on_nat_sentinel(self, unit):
        # GH#66552 stepping one unit past the bottom of the range lands on
        #  iNaT, which is not NaT but is indistinguishable from it once
        #  stored, so it has to raise rather than come back as NaT
        td_min = Timedelta(np.timedelta64(-(2**63) + 1, unit))
        one = np.array([1], dtype=f"m8[{unit}]")

        attrname = {"s": "second", "ms": "millisecond", "us": "microsecond"}.get(
            unit, "nanosecond"
        )
        msg = f"Out of bounds {attrname} timedelta"
        with pytest.raises(OutOfBoundsTimedelta, match=msg):
            td_min - one
        with pytest.raises(OutOfBoundsTimedelta, match=msg):
            td_min + (-one)
        with pytest.raises(OutOfBoundsTimedelta, match=msg):
            (-one) + td_min

        # __rsub__ reaches the sentinel from the other end of the range
        td_max = Timedelta(np.timedelta64(2**63 - 1, unit))
        with pytest.raises(OutOfBoundsTimedelta, match=msg):
            (-one) - td_max

    def test_td_add_sub_td64_ndarray_overflow_in_cast(self):
        # GH#66552 the promotion to the finer of the two units used to wrap
        td = Timedelta(10**12, "s")
        other = np.array([1], dtype="m8[ns]")

        msg = "Cannot cast 11574074 days 01:46:40 to unit='ns' without overflow"
        with pytest.raises(OutOfBoundsTimedelta, match=msg):
            td + other
        with pytest.raises(OutOfBoundsTimedelta, match=msg):
            other - td

        # the coarser operand is the one cast when the Timedelta is finer
        td = Timedelta(1, "ns")
        other = np.array([2**62], dtype="m8[s]")

        msg = "Cannot convert 4611686018427387904 seconds to timedelta64"
        with pytest.raises(OutOfBoundsTimedelta, match=msg):
            td + other
        with pytest.raises(OutOfBoundsTimedelta, match=msg):
            other - td

    def test_td_add_sub_td64_ndarray_nat_operand(self):
        # GH#66552 a missing operand is still missing, not an overflow
        td = Timedelta(1, "ns")
        nat = np.timedelta64("NaT", "ns")
        other = np.array([1, nat], dtype="m8[ns]")

        expected = np.array([2, nat], dtype="m8[ns]")
        tm.assert_numpy_array_equal(td + other, expected)
        tm.assert_numpy_array_equal(other + td, expected)

        expected = np.array([0, nat], dtype="m8[ns]")
        tm.assert_numpy_array_equal(td - other, expected)
        tm.assert_numpy_array_equal(other - td, expected)

    def test_td_add_sub_td64_ndarray_byteswapped(self, unit):
        # GH#66552 the values are viewed as i8, which a non-native buffer
        #  would misread
        td = Timedelta(1, unit)
        other = np.array([1, 2], dtype=np.dtype(f"m8[{unit}]").newbyteorder(">"))

        expected = np.array([2, 3], dtype=f"m8[{unit}]")
        tm.assert_numpy_array_equal(td + other, expected)
        tm.assert_numpy_array_equal(other + td, expected)

    @pytest.mark.parametrize("freq", ["Y", "M"])
    def test_td_add_sub_td64_ndarray_year_month_unit(self, freq):
        # GH#66552 numpy itself refuses to combine these with a time unit; we
        #  leave them to it rather than inventing a common unit
        other = np.array([1], dtype=f"m8[{freq}]")

        msg = "Cannot get a common metadata divisor"
        with pytest.raises(TypeError, match=msg):
            Timedelta(1, "ns") + other
        with pytest.raises(TypeError, match=msg):
            other - Timedelta(1, "ns")

    def test_td_add_sub_dt64_ndarray(self):
        td = Timedelta("1 day")
        other = np.array(["2000-01-01"], dtype="M8[ns]")

        expected = np.array(["2000-01-02"], dtype="M8[ns]")
        tm.assert_numpy_array_equal(td + other, expected)
        tm.assert_numpy_array_equal(other + td, expected)

        expected = np.array(["1999-12-31"], dtype="M8[ns]")
        tm.assert_numpy_array_equal(-td + other, expected)
        tm.assert_numpy_array_equal(other - td, expected)

    def test_td_add_sub_dt64_ndarray_out_of_bounds(self, unit):
        # GH#66552 stepping past the top of the range used to wrap
        td = Timedelta(1, unit)
        attrname = {"s": "second", "ms": "millisecond", "us": "microsecond"}.get(
            unit, "nanosecond"
        )
        msg = f"Out of bounds {attrname} timestamp"

        dt_max = np.array([2**63 - 1], dtype=f"M8[{unit}]")
        with pytest.raises(OutOfBoundsDatetime, match=msg):
            td + dt_max
        with pytest.raises(OutOfBoundsDatetime, match=msg):
            dt_max + td

        # stepping one unit past the bottom lands on iNaT, which is not NaT but
        #  is indistinguishable from it once stored, so it has to raise too
        dt_min = np.array([-(2**63) + 1], dtype=f"M8[{unit}]")
        with pytest.raises(OutOfBoundsDatetime, match=msg):
            dt_min - td
        with pytest.raises(OutOfBoundsDatetime, match=msg):
            -td + dt_min

    def test_td_add_sub_dt64_ndarray_out_of_bounds_message_unit(self):
        # GH#66552 the message names the promoted resolution, not the
        #  Timedelta's own
        td = Timedelta(1, "s")

        msg = "Out of bounds nanosecond timestamp"
        with pytest.raises(OutOfBoundsDatetime, match=msg):
            td + np.array([2**63 - 1], dtype="M8[ns]")
        with pytest.raises(OutOfBoundsDatetime, match=msg):
            np.array([-(2**63) + 1], dtype="M8[ns]") - td

    def test_td_add_sub_dt64_ndarray_overflow_in_cast(self):
        # GH#66552 the promotion to the finer of the two units used to wrap
        td = Timedelta(1, "ns")
        other = np.array([2**40], dtype="M8[s]")

        msg = "Out of bounds nanosecond timestamp: 36812-02-20"
        with pytest.raises(OutOfBoundsDatetime, match=msg):
            td + other
        with pytest.raises(OutOfBoundsDatetime, match=msg):
            other - td

        # the Timedelta is the one cast when it is the coarser of the two
        td = Timedelta(10**12, "s")
        other = np.array([1], dtype="M8[ns]")

        msg = "Cannot cast 11574074 days 01:46:40 to unit='ns' without overflow"
        with pytest.raises(OutOfBoundsTimedelta, match=msg):
            td + other
        with pytest.raises(OutOfBoundsTimedelta, match=msg):
            other - td

    def test_td_add_sub_dt64_ndarray_nat_operand(self):
        # GH#66552 a missing operand is still missing, not an overflow
        td = Timedelta(1, "ns")
        nat = np.datetime64("NaT", "ns")
        other = np.array([1, nat], dtype="M8[ns]")

        expected = np.array([2, nat], dtype="M8[ns]")
        tm.assert_numpy_array_equal(td + other, expected)
        tm.assert_numpy_array_equal(other + td, expected)

        expected = np.array([0, nat], dtype="M8[ns]")
        tm.assert_numpy_array_equal(other - td, expected)

    def test_td_add_sub_dt64_ndarray_byteswapped(self, unit):
        # GH#66552 the values are viewed as i8, which a non-native buffer
        #  would misread
        td = Timedelta(1, unit)
        other = np.array([1, 2], dtype=np.dtype(f"M8[{unit}]").newbyteorder(">"))

        expected = np.array([2, 3], dtype=f"M8[{unit}]")
        tm.assert_numpy_array_equal(td + other, expected)
        tm.assert_numpy_array_equal(other + td, expected)

        # the cast to a common unit has to survive the non-native buffer too
        other = np.array(["2000-01-01"], dtype=np.dtype("M8[s]").newbyteorder(">"))
        expected = np.array(["2000-01-01"], dtype="M8[s]") + td.to_timedelta64()
        tm.assert_numpy_array_equal(td + other, expected)
        tm.assert_numpy_array_equal(other + td, expected)

    @pytest.mark.parametrize("freq", ["Y", "M", "W", "D", "h", "m", "s", "ms", "us"])
    def test_td_add_sub_dt64_ndarray_mixed_units(self, freq, unit):
        # GH#66552 matching numpy, the operands are cast to the finer of the
        #  two units; in-bounds results, NaT included, are unchanged
        td = Timedelta(1, unit)
        other = np.array(["2000-01-01", "NaT"], dtype=f"M8[{freq}]")
        m8 = td.to_timedelta64()

        tm.assert_numpy_array_equal(td + other, m8 + other)
        tm.assert_numpy_array_equal(other + td, m8 + other)
        tm.assert_numpy_array_equal(other - td, other - m8)

    @pytest.mark.parametrize("freq", ["Y", "M"])
    def test_td_add_sub_dt64_ndarray_year_month_unit_overflow(self, freq):
        # GH#66552 numpy casts these to the Timedelta's unit and wraps
        td = Timedelta(1, "ns")
        other = np.array(["4940-01-01"], dtype=f"M8[{freq}]")

        msg = "Out of bounds nanosecond timestamp: 4940-01-01"
        with pytest.raises(OutOfBoundsDatetime, match=msg):
            td + other
        with pytest.raises(OutOfBoundsDatetime, match=msg):
            other - td

    @pytest.mark.parametrize("freq", ["ps", "fs", "as"])
    def test_td_add_sub_dt64_ndarray_subnano_unit(self, freq):
        # GH#66552 we have no reso for these, so numpy keeps handling them
        td = Timedelta(1, "ns")
        other = np.array([1], dtype=f"M8[{freq}]")
        m8 = td.to_timedelta64()

        tm.assert_numpy_array_equal(td + other, m8 + other)
        tm.assert_numpy_array_equal(other + td, m8 + other)
        tm.assert_numpy_array_equal(other - td, other - m8)

    def test_td_add_sub_dt64_ndarray_generic_unit(self, unit):
        # GH#66552 numpy reads a generic datetime64 in the other operand's unit
        td = Timedelta(1, unit)
        other = np.zeros(2, dtype="M8")

        expected = np.array([1, 1], dtype=f"M8[{unit}]")
        tm.assert_numpy_array_equal(td + other, expected)
        tm.assert_numpy_array_equal(other + td, expected)

        expected = np.array([-1, -1], dtype=f"M8[{unit}]")
        tm.assert_numpy_array_equal(other - td, expected)

    def test_td_sub_dt64_ndarray_invalid(self):
        # GH#66552 numpy refuses this and so do we
        td = Timedelta(1, "ns")
        other = np.array(["2000-01-01"], dtype="M8[ns]")

        msg = "ufunc 'subtract' cannot use operands with types"
        with pytest.raises(TypeError, match=msg):
            td - other

    def test_td_add_sub_ndarray_0d(self):
        td = Timedelta("1 day")
        other = np.array(td.asm8)

        result = td + other
        assert isinstance(result, Timedelta)
        assert result == 2 * td

        result = other + td
        assert isinstance(result, Timedelta)
        assert result == 2 * td

        result = other - td
        assert isinstance(result, Timedelta)
        assert result == 0 * td

        result = td - other
        assert isinstance(result, Timedelta)
        assert result == 0 * td


class TestTimedeltaMultiplicationDivision:
    """
    Tests for Timedelta methods:

        __mul__, __rmul__,
        __div__, __rdiv__,
        __truediv__, __rtruediv__,
        __floordiv__, __rfloordiv__,
        __mod__, __rmod__,
        __divmod__, __rdivmod__
    """

    # ---------------------------------------------------------------
    # Timedelta.__mul__, __rmul__

    @pytest.mark.parametrize("td_nat", [NaT, np.timedelta64("NaT", "ns")])
    @pytest.mark.parametrize("op", [operator.mul, ops.rmul])
    def test_td_mul_nat(self, op, td_nat):
        # GH#19819
        td = Timedelta(10, unit="D")
        typs = "|".join(["numpy.timedelta64", "NaTType", "Timedelta"])
        msg = "|".join(
            [
                rf"unsupported operand type\(s\) for \*: '{typs}' and '{typs}'",
                r"ufunc '?multiply'? cannot use operands with types",
            ]
        )
        with pytest.raises(TypeError, match=msg):
            op(td, td_nat)

    @pytest.mark.parametrize("nan", [np.nan, np.float64("NaN"), float("nan")])
    @pytest.mark.parametrize("op", [operator.mul, ops.rmul])
    def test_td_mul_nan(self, op, nan):
        # np.float64('NaN') has a 'dtype' attr, avoid treating as array
        td = Timedelta(10, unit="D")
        result = op(td, nan)
        assert result is NaT

    @pytest.mark.parametrize("op", [operator.mul, ops.rmul])
    def test_td_mul_scalar(self, op):
        # GH#19738
        td = Timedelta(minutes=3)

        result = op(td, 2)
        assert result == Timedelta(minutes=6)

        result = op(td, 1.5)
        assert result == Timedelta(minutes=4, seconds=30)

        assert op(td, np.nan) is NaT

        assert op(-1, td)._value == -1 * td._value
        assert op(-1.0, td)._value == -1.0 * td._value

        msg = "unsupported operand type"
        with pytest.raises(TypeError, match=msg):
            # timedelta * datetime is gibberish
            op(td, Timestamp(2016, 1, 2))

        with pytest.raises(TypeError, match=msg):
            # invalid multiply with another timedelta
            op(td, td)

    @pytest.mark.parametrize("dtype", [np.int64, np.uint64, np.int32, np.int8])
    @pytest.mark.parametrize("op", [operator.mul, ops.rmul])
    def test_td_mul_numpy_integer_overflow(self, op, dtype):
        # GH#66551 numpy integer scalars used to keep the multiply in the
        #  operand's own dtype, so an out-of-bounds product wrapped silently
        #  instead of raising the way a Python int does.
        td = Timedelta(2**62, unit="ns")

        msg = "|".join(
            [
                "Python int too large to convert to C long",
                # windows, 32bit linux builds
                "int too big to convert",
            ]
        )
        with pytest.raises(OverflowError, match=msg):
            op(td, dtype(4))

    @pytest.mark.parametrize("dtype", [np.float32, np.float64])
    @pytest.mark.parametrize("op", [operator.mul, ops.rmul])
    def test_td_mul_numpy_float_precision(self, op, dtype):
        # GH#66551 a float32 multiplier used to drag the product down to
        #  float32 precision, losing 16 seconds here.
        td = Timedelta(10**18, unit="ns")

        result = op(td, dtype(1.0))
        assert result == td

    def test_td_mul_numeric_ndarray(self):
        td = Timedelta("1 day")
        other = np.array([2])
        expected = np.array([Timedelta("2 Days").to_timedelta64()])

        result = td * other
        tm.assert_numpy_array_equal(result, expected)

        result = other * td
        tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.parametrize("unit", ["s", "ms", "us", "ns"])
    def test_td_mul_int_ndarray_overflow(self, unit):
        # GH#66552 the product used to be formed in int64 and wrap silently,
        #  where the TimedeltaIndex equivalent raises
        td = Timedelta(4, unit)
        other = np.array([2**62])

        msg = "Overflow in int64 multiplication"
        with pytest.raises(OutOfBoundsTimedelta, match=msg):
            td * other
        with pytest.raises(OutOfBoundsTimedelta, match=msg):
            other * td

        # int64.min is representable but would be misread as NaT
        with pytest.raises(OutOfBoundsTimedelta, match=msg):
            Timedelta(1, unit) * np.array([-(2**63)])

        # one step in from the sentinel is fine
        result = Timedelta(1, unit) * np.array([-(2**63) + 1])
        expected = np.array([-(2**63) + 1], dtype=f"m8[{unit}]")
        tm.assert_numpy_array_equal(result, expected)

    def test_td_mul_uint64_ndarray_overflow(self):
        # GH#66552 a multiplier above int64.max wraps negative in the i8 cast
        td = Timedelta(4, "ns")
        other = np.array([2**63], dtype=np.uint64)

        msg = "Overflow in int64 multiplication"
        with pytest.raises(OutOfBoundsTimedelta, match=msg):
            td * other
        with pytest.raises(OutOfBoundsTimedelta, match=msg):
            other * td

    @pytest.mark.parametrize("factor", [1e30, np.inf, -np.inf])
    def test_td_mul_float_ndarray_overflow(self, factor):
        # GH#66552 the product used to saturate to int64.max on the cast, or
        #  come back as NaT for an infinite multiplier
        td = Timedelta(4, "ns")
        other = np.array([factor])

        msg = "Overflow in timedelta multiplication"
        with pytest.raises(OutOfBoundsTimedelta, match=msg):
            td * other
        with pytest.raises(OutOfBoundsTimedelta, match=msg):
            other * td

    def test_td_mul_float_ndarray_nan(self):
        # GH#66552 a NaN multiplier is still missing, not an overflow
        td = Timedelta(4, "ns")
        other = np.array([np.nan, 2.0])

        expected = np.array([np.timedelta64("NaT", "ns"), np.timedelta64(8, "ns")])
        tm.assert_numpy_array_equal(td * other, expected)
        tm.assert_numpy_array_equal(other * td, expected)

    def test_td_mul_float32_ndarray_precision(self):
        # GH#66552 a float32 multiplier used to drag the product down to
        #  float32 precision
        td = Timedelta(10**18, "ns")
        other = np.array([1.0], dtype=np.float32)

        expected = np.array([10**18], dtype="m8[ns]")
        tm.assert_numpy_array_equal(td * other, expected)
        tm.assert_numpy_array_equal(other * td, expected)

    def test_td_mul_numeric_ndarray_byteswapped(self):
        # GH#66552 a non-native buffer has to be converted before the multiply
        td = Timedelta(4, "ns")
        other = np.array([2, 3], dtype=np.dtype("i8").newbyteorder(">"))

        expected = np.array([8, 12], dtype="m8[ns]")
        tm.assert_numpy_array_equal(td * other, expected)
        tm.assert_numpy_array_equal(other * td, expected)

    def test_td_mul_numeric_ndarray_0d(self):
        td = Timedelta("1 day")
        other = np.array(2, dtype=np.int64)
        assert other.ndim == 0
        expected = Timedelta("2 days")

        res = td * other
        assert type(res) is Timedelta
        assert res == expected

        res = other * td
        assert type(res) is Timedelta
        assert res == expected

    def test_td_mul_td64_ndarray_invalid(self):
        td = Timedelta("1 day")
        other = np.array([Timedelta("2 Days").to_timedelta64()])

        msg = (
            "ufunc '?multiply'? cannot use operands with types "
            rf"dtype\('{tm.ENDIAN}m8\[us\]'\) and dtype\('{tm.ENDIAN}m8\[us\]'\)"
        )
        with pytest.raises(TypeError, match=msg):
            td * other
        with pytest.raises(TypeError, match=msg):
            other * td

    # ---------------------------------------------------------------
    # Timedelta.__div__, __truediv__

    def test_td_div_timedeltalike_scalar(self):
        # GH#19738
        td = Timedelta(10, unit="D")

        result = td / offsets.Hour(1)
        assert result == 240

        assert td / td == 1
        assert td / np.timedelta64(60, "h") == 4

        assert np.isnan(td / NaT)

    def test_td_div_td64_non_nano(self):
        # truediv
        td = Timedelta("1 days 2 hours 3 ns")
        result = td / np.timedelta64(1, "D")
        assert result == td._value / (86400 * 10**9)
        result = td / np.timedelta64(1, "s")
        assert result == td._value / 10**9
        result = td / np.timedelta64(1, "ns")
        assert result == td._value

        # floordiv
        td = Timedelta("1 days 2 hours 3 ns")
        result = td // np.timedelta64(1, "D")
        assert result == 1
        result = td // np.timedelta64(1, "s")
        assert result == 93600
        result = td // np.timedelta64(1, "ns")
        assert result == td._value

    def test_td_div_numeric_scalar(self):
        # GH#19738
        td = Timedelta(10, unit="D")

        result = td / 2
        assert isinstance(result, Timedelta)
        assert result == Timedelta(days=5)

        result = td / 5
        assert isinstance(result, Timedelta)
        assert result == Timedelta(days=2)

    @pytest.mark.parametrize(
        "value, divisor, expected",
        [
            (Timedelta.min._value, 1, Timedelta.min._value),
            (Timedelta.max._value, 1, Timedelta.max._value),
            (Timedelta.min._value, -1, Timedelta.max._value),
            (2**53 + 1, 1, 2**53 + 1),
            # truncation toward zero, matching numpy
            (36028797018963967, 2, 18014398509481983),
            (-36028797018963967, 2, -18014398509481983),
            (36028797018963967, -2, -18014398509481983),
        ],
    )
    def test_td_div_integer_exact(self, value, divisor, expected):
        # GH#66551 the quotient was computed in float64, whose 53-bit mantissa
        #  rounds values that int64 holds exactly, so dividing by 1 gave NaT at
        #  Timedelta.min and raised OverflowError at Timedelta.max.
        td = Timedelta(value, unit="ns")

        result = td / divisor
        assert result._value == expected

        # the vectorized path was already exact; the scalar disagreed with it
        assert (pd.array([td]) / divisor)[0]._value == expected

    @pytest.mark.parametrize(
        "nan",
        [
            np.nan,
            np.float64("NaN"),
            float("nan"),
        ],
    )
    def test_td_div_nan(self, nan):
        # np.float64('NaN') has a 'dtype' attr, avoid treating as array
        td = Timedelta(10, unit="D")
        result = td / nan
        assert result is NaT

        result = td // nan
        assert result is NaT

    def test_td_div_td64_ndarray(self):
        td = Timedelta("1 day")

        other = np.array([Timedelta("2 Days").to_timedelta64()])
        expected = np.array([0.5])

        result = td / other
        tm.assert_numpy_array_equal(result, expected)

        result = other / td
        tm.assert_numpy_array_equal(result, expected * 4)

    def test_td_div_ndarray_0d(self):
        td = Timedelta("1 day")

        other = np.array(1)
        res = td / other
        assert isinstance(res, Timedelta)
        assert res == td

    # ---------------------------------------------------------------
    # Timedelta.__rdiv__

    def test_td_rdiv_timedeltalike_scalar(self):
        # GH#19738
        td = Timedelta(10, unit="D")
        result = offsets.Hour(1) / td
        assert result == 1 / 240.0

        assert np.timedelta64(60, "h") / td == 0.25

    def test_td_rdiv_na_scalar(self):
        # GH#31869 None gets cast to NaT
        td = Timedelta(10, unit="D")

        result = NaT / td
        assert np.isnan(result)

        result = None / td
        assert np.isnan(result)

        result = np.timedelta64("NaT", "ns") / td
        assert np.isnan(result)

        msg = r"unsupported operand type\(s\) for /: 'numpy.datetime64' and 'Timedelta'"
        with pytest.raises(TypeError, match=msg):
            np.datetime64("NaT", "ns") / td

        msg = r"unsupported operand type\(s\) for /: 'float' and 'Timedelta'"
        with pytest.raises(TypeError, match=msg):
            np.nan / td

    def test_td_rdiv_ndarray(self):
        td = Timedelta(10, unit="D")

        arr = np.array([td], dtype=object)
        result = arr / td
        expected = np.array([1], dtype=np.float64)
        tm.assert_numpy_array_equal(result, expected)

        arr = np.array([None])
        result = arr / td
        expected = np.array([np.nan])
        tm.assert_numpy_array_equal(result, expected)

        arr = np.array([np.nan], dtype=object)
        msg = r"unsupported operand type\(s\) for /: 'float' and 'Timedelta'"
        with pytest.raises(TypeError, match=msg):
            arr / td

        arr = np.array([np.nan], dtype=np.float64)
        msg = "cannot use operands with types dtype"
        with pytest.raises(TypeError, match=msg):
            arr / td

    def test_td_rdiv_ndarray_0d(self):
        td = Timedelta(10, unit="D")

        arr = np.array(td.asm8)

        assert arr / td == 1

    # ---------------------------------------------------------------
    # Timedelta.__floordiv__

    def test_td_floordiv_timedeltalike_scalar(self):
        # GH#18846
        td = Timedelta(hours=3, minutes=4)
        scalar = Timedelta(hours=3, minutes=3)

        assert td // scalar == 1
        assert -td // scalar.to_pytimedelta() == -2
        assert (2 * td) // scalar.to_timedelta64() == 2

    def test_td_floordiv_null_scalar(self):
        # GH#18846
        td = Timedelta(hours=3, minutes=4)

        assert td // np.nan is NaT
        assert np.isnan(td // NaT)
        assert np.isnan(td // np.timedelta64("NaT", "ns"))

    def test_td_floordiv_offsets(self):
        # GH#19738
        td = Timedelta(hours=3, minutes=4)
        assert td // offsets.Hour(1) == 3
        assert td // offsets.Minute(2) == 92

    def test_td_floordiv_invalid_scalar(self):
        # GH#18846
        td = Timedelta(hours=3, minutes=4)

        # CPython 3.13 reworded the invalid-keyword error raised by the
        # np.datetime64 constructor (boundary confirmed on CI: 3.12 old,
        # 3.13 new; independent of NumPy version).
        if PY313:
            msg = "this function got an unexpected keyword argument 'dtype'"
        else:
            msg = "|".join(
                [
                    "'dtype' is an invalid keyword argument for this function",
                    r"ufunc '?floor_divide'? cannot use operands with types",
                ]
            )
        with pytest.raises(TypeError, match=msg):
            td // np.datetime64("2016-01-01", dtype="datetime64[us]")

    def test_td_floordiv_numeric_scalar(self):
        # GH#18846
        td = Timedelta(hours=3, minutes=4)

        expected = Timedelta(hours=1, minutes=32)
        assert td // 2 == expected
        assert td // 2.0 == expected
        assert td // np.float64(2.0) == expected
        assert td // np.int32(2.0) == expected
        assert td // np.uint8(2.0) == expected

    def test_td_floordiv_timedeltalike_array(self):
        # GH#18846
        td = Timedelta(hours=3, minutes=4)
        scalar = Timedelta(hours=3, minutes=3)

        # Array-like others
        assert td // np.array(scalar.to_timedelta64()) == 1

        res = (3 * td) // np.array([scalar.to_timedelta64()])
        expected = np.array([3], dtype=np.int64)
        tm.assert_numpy_array_equal(res, expected)

        res = (10 * td) // np.array(
            [scalar.to_timedelta64(), np.timedelta64("NaT", "ns")]
        )
        expected = np.array([10, np.nan])
        tm.assert_numpy_array_equal(res, expected)

    def test_td_floordiv_numeric_series(self):
        # GH#18846
        td = Timedelta(hours=3, minutes=4)
        ser = pd.Series([1], dtype=np.int64)
        res = td // ser
        assert res.dtype.kind == "m"

    # ---------------------------------------------------------------
    # Timedelta.__rfloordiv__

    def test_td_rfloordiv_timedeltalike_scalar(self):
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

    def test_td_rfloordiv_null_scalar(self):
        # GH#18846
        td = Timedelta(hours=3, minutes=3)

        assert np.isnan(td.__rfloordiv__(NaT))
        assert np.isnan(td.__rfloordiv__(np.timedelta64("NaT", "ns")))

    def test_td_rfloordiv_offsets(self):
        # GH#19738
        assert offsets.Hour(1) // Timedelta(minutes=25) == 2

    def test_td_rfloordiv_invalid_scalar(self):
        # GH#18846
        td = Timedelta(hours=3, minutes=3)

        dt64 = np.datetime64("2016-01-01", "us")

        assert td.__rfloordiv__(dt64) is NotImplemented

        msg = (
            r"unsupported operand type\(s\) for //: 'numpy.datetime64' and 'Timedelta'"
        )
        with pytest.raises(TypeError, match=msg):
            dt64 // td

    def test_td_rfloordiv_numeric_scalar(self):
        # GH#18846
        td = Timedelta(hours=3, minutes=3)

        assert td.__rfloordiv__(np.nan) is NotImplemented
        assert td.__rfloordiv__(3.5) is NotImplemented
        assert td.__rfloordiv__(2) is NotImplemented
        assert td.__rfloordiv__(np.float64(2.0)) is NotImplemented
        assert td.__rfloordiv__(np.uint8(9)) is NotImplemented
        assert td.__rfloordiv__(np.int32(2.0)) is NotImplemented

        msg = r"unsupported operand type\(s\) for //: '.*' and 'Timedelta"
        with pytest.raises(TypeError, match=msg):
            np.float64(2.0) // td
        with pytest.raises(TypeError, match=msg):
            np.uint8(9) // td
        with pytest.raises(TypeError, match=msg):
            # deprecated GH#19761, enforced GH#29797
            np.int32(2.0) // td

    def test_td_rfloordiv_timedeltalike_array(self):
        # GH#18846
        td = Timedelta(hours=3, minutes=3)
        scalar = Timedelta(hours=3, minutes=4)

        # Array-like others
        assert td.__rfloordiv__(np.array(scalar.to_timedelta64())) == 1

        res = td.__rfloordiv__(np.array([(3 * scalar).to_timedelta64()]))
        expected = np.array([3], dtype=np.int64)
        tm.assert_numpy_array_equal(res, expected)

        arr = np.array([(10 * scalar).to_timedelta64(), np.timedelta64("NaT", "ns")])
        res = td.__rfloordiv__(arr)
        expected = np.array([10, np.nan])
        tm.assert_numpy_array_equal(res, expected)

    def test_td_rfloordiv_intarray(self):
        # deprecated GH#19761, enforced GH#29797
        ints = np.array([1349654400, 1349740800, 1349827200, 1349913600]) * 10**9

        msg = "Invalid dtype"
        with pytest.raises(TypeError, match=msg):
            ints // Timedelta(1, unit="s")

    def test_td_rfloordiv_numeric_series(self):
        # GH#18846
        td = Timedelta(hours=3, minutes=3)
        ser = pd.Series([1], dtype=np.int64)
        res = td.__rfloordiv__(ser)
        assert res is NotImplemented

        msg = "Invalid dtype"
        with pytest.raises(TypeError, match=msg):
            # Deprecated GH#19761, enforced GH#29797
            ser // td

    # ----------------------------------------------------------------
    # Timedelta.__mod__, __rmod__

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

    def test_mod_timedelta64_nat(self):
        # GH#19365
        td = Timedelta(hours=37)

        result = td % np.timedelta64("NaT", "ns")
        assert result is NaT

    def test_mod_timedelta64(self):
        # GH#19365
        td = Timedelta(hours=37)

        result = td % np.timedelta64(2, "h")
        assert isinstance(result, Timedelta)
        assert result == Timedelta(hours=1)

    def test_mod_offset(self):
        # GH#19365
        td = Timedelta(hours=37)

        result = td % offsets.Hour(5)
        assert isinstance(result, Timedelta)
        assert result == Timedelta(hours=2)

    def test_mod_numeric(self):
        # GH#19365
        td = Timedelta(hours=37)

        # Numeric Others
        result = td % 2
        assert isinstance(result, Timedelta)
        assert result == Timedelta(0)

        result = td % 1e9
        assert isinstance(result, Timedelta)
        assert result == Timedelta(minutes=3, seconds=20)

        result = td % int(1e9)
        assert isinstance(result, Timedelta)
        assert result == Timedelta(minutes=3, seconds=20)

    def test_mod_invalid(self):
        # GH#19365
        td = Timedelta(hours=37)
        msg = "unsupported operand type"
        with pytest.raises(TypeError, match=msg):
            td % Timestamp("2018-01-22")

        with pytest.raises(TypeError, match=msg):
            td % []

    def test_rmod_pytimedelta(self):
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

    def test_rmod_invalid(self):
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

    # ----------------------------------------------------------------
    # Timedelta.__divmod__, __rdivmod__

    def test_divmod_numeric(self):
        # GH#19365
        td = Timedelta(days=2, hours=6)

        result = divmod(td, 53 * 3600 * 1e6)
        assert result[0] == Timedelta(1, unit="us").as_unit("us")
        assert isinstance(result[1], Timedelta)
        assert result[1] == Timedelta(hours=1)

        assert result
        result = divmod(td, np.nan)
        assert result[0] is NaT
        assert result[1] is NaT

    def test_divmod(self):
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

    def test_divmod_offset(self):
        # GH#19365
        td = Timedelta(days=2, hours=6)

        result = divmod(td, offsets.Hour(-4))
        assert result[0] == -14
        assert isinstance(result[1], Timedelta)
        assert result[1] == Timedelta(hours=-2)

    def test_divmod_invalid(self):
        # GH#19365
        td = Timedelta(days=2, hours=6)

        msg = r"unsupported operand type\(s\) for //: 'Timedelta' and 'Timestamp'"
        with pytest.raises(TypeError, match=msg):
            divmod(td, Timestamp("2018-01-22"))

    def test_rdivmod_pytimedelta(self):
        # GH#19365
        result = divmod(timedelta(days=2, hours=6), Timedelta(days=1))
        assert result[0] == 2
        assert isinstance(result[1], Timedelta)
        assert result[1] == Timedelta(hours=6)

    def test_rdivmod_offset(self):
        result = divmod(offsets.Hour(54), Timedelta(hours=-4))
        assert result[0] == -14
        assert isinstance(result[1], Timedelta)
        assert result[1] == Timedelta(hours=-2)

    def test_rdivmod_invalid(self):
        # GH#19365
        td = Timedelta(minutes=3)
        msg = "unsupported operand type"

        with pytest.raises(TypeError, match=msg):
            divmod(Timestamp("2018-01-22"), td)

        with pytest.raises(TypeError, match=msg):
            divmod(15, td)

        with pytest.raises(TypeError, match=msg):
            divmod(16.0, td)

        msg = "Invalid dtype int"
        with pytest.raises(TypeError, match=msg):
            divmod(np.array([22, 24]), td)

    # ----------------------------------------------------------------

    @pytest.mark.parametrize(
        "op", [operator.mul, ops.rmul, operator.truediv, ops.rdiv, ops.rsub]
    )
    @pytest.mark.parametrize(
        "arr",
        [
            [Timestamp("20130101 9:01"), Timestamp("20121230 9:02")],
            [Timestamp("2021-11-09 09:54:00"), Timedelta("1D")],
        ],
    )
    def test_td_op_timedelta_timedeltalike_array(self, op, arr):
        arr = np.array(arr)
        msg = "|".join(["unsupported operand type", "cannot use operands with types"])
        with pytest.raises(TypeError, match=msg):
            op(arr, Timedelta("1D"))

    def test_mul_bool_invalid(self):
        # GH#62316
        msg = "Cannot multiply Timedelta by bool. Explicitly cast to integer"
        with pytest.raises(TypeError, match=msg):
            Timedelta("1 day") * True


class TestTimedeltaComparison:
    def test_compare_pytimedelta_bounds(self):
        # GH#49021 don't overflow on comparison with very large pytimedeltas

        for unit in ["ns", "us"]:
            tdmax = Timedelta.max.as_unit(unit).max
            tdmin = Timedelta.min.as_unit(unit).min

            assert tdmax < timedelta.max
            assert tdmax <= timedelta.max
            assert not tdmax > timedelta.max
            assert not tdmax >= timedelta.max
            assert tdmax != timedelta.max
            assert not tdmax == timedelta.max

            assert tdmin > timedelta.min
            assert tdmin >= timedelta.min
            assert not tdmin < timedelta.min
            assert not tdmin <= timedelta.min
            assert tdmin != timedelta.min
            assert not tdmin == timedelta.min

        # But the "ms" and "s"-reso bounds extend pass pytimedelta
        for unit in ["ms", "s"]:
            tdmax = Timedelta.max.as_unit(unit).max
            tdmin = Timedelta.min.as_unit(unit).min

            assert tdmax > timedelta.max
            assert tdmax >= timedelta.max
            assert not tdmax < timedelta.max
            assert not tdmax <= timedelta.max
            assert tdmax != timedelta.max
            assert not tdmax == timedelta.max

            assert tdmin < timedelta.min
            assert tdmin <= timedelta.min
            assert not tdmin > timedelta.min
            assert not tdmin >= timedelta.min
            assert tdmin != timedelta.min
            assert not tdmin == timedelta.min

    def test_compare_pytimedelta_bounds2(self):
        # a pytimedelta outside the microsecond bounds
        pytd = timedelta(days=999999999, seconds=86399)
        # NB: np.timedelta64(td, "s"") incorrectly overflows
        td64 = np.timedelta64(pytd.days, "D") + np.timedelta64(pytd.seconds, "s")
        td = Timedelta(td64)
        assert td.days == pytd.days
        assert td.seconds == pytd.seconds

        assert td == pytd
        assert not td != pytd
        assert not td < pytd
        assert not td > pytd
        assert td <= pytd
        assert td >= pytd

        td2 = td - Timedelta(seconds=1).as_unit("s")
        assert td2 != pytd
        assert not td2 == pytd
        assert td2 < pytd
        assert td2 <= pytd
        assert not td2 > pytd
        assert not td2 >= pytd

    def test_compare_tick(self, tick_classes):
        cls = tick_classes

        off = cls(4)
        td = off._as_pd_timedelta
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

        assert t != "string"
        assert t != 1
        assert t != CustomClass()
        assert t != CustomClass(cmp_result=False)

        assert t < CustomClass(cmp_result=True)
        assert not t < CustomClass(cmp_result=False)

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


def test_ops_notimplemented():
    class Other:
        pass

    other = Other()

    td = Timedelta("1 day")
    assert td.__add__(other) is NotImplemented
    assert td.__sub__(other) is NotImplemented
    assert td.__truediv__(other) is NotImplemented
    assert td.__mul__(other) is NotImplemented
    assert td.__floordiv__(other) is NotImplemented


def test_ops_error_str():
    # GH#13624
    td = Timedelta("1 day")

    for left, right in [(td, "a"), ("a", td)]:
        msg = "|".join(
            [
                "unsupported operand type",
                r'can only concatenate str \(not "Timedelta"\) to str',
                "must be str, not Timedelta",
            ]
        )
        with pytest.raises(TypeError, match=msg):
            left + right

        msg = "not supported between instances of"
        with pytest.raises(TypeError, match=msg):
            left > right

        assert not left == right
        assert left != right


@pytest.mark.parametrize("box", [True, False])
def test_ops_str_deprecated(box):
    # GH#59653
    td = Timedelta("1 day")
    item = "1"
    if box:
        item = np.array([item], dtype=object)

    msg = "Scalar operations between Timedelta and string are deprecated"
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        td + item
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        item + td
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        td - item
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        item - td
    with tm.assert_produces_warning(Pandas4Warning, match=msg):
        item / td
    if not box:
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            td / item
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            item // td
        with tm.assert_produces_warning(Pandas4Warning, match=msg):
            td // item
    else:
        # true division dispatches to NumPy; older NumPy raised via the ufunc
        div_msg = "|".join(
            [
                r"unsupported operand type\(s\) for /: 'datetime.timedelta' and 'str'",
                "ufunc 'divide' cannot use operands",  # older NumPy
            ]
        )
        with pytest.raises(TypeError, match=div_msg):
            td / item
        # floor division (either operand order) raises on the object dtype
        floordiv_msg = "|".join(
            [
                "Invalid dtype object for __floordiv__",
                r"unsupported operand type\(s\) for /: 'int' and 'str'",  # older NumPy
            ]
        )
        with pytest.raises(TypeError, match=floordiv_msg):
            item // td
        with pytest.raises(TypeError, match=floordiv_msg):
            td // item


@pytest.mark.parametrize("unit", ["s", "ms", "us", "ns"])
@pytest.mark.parametrize("factor", [-(2**63), float(-(2**63))])
def test_td_mul_lands_on_nat_sentinel(unit, factor):
    # GH#66551 the result is iNaT, which is not NaT but is indistinguishable
    #  from it once stored, so it has to raise rather than be constructed.
    #  The guard used to be a bare `assert`, which `python -O` strips.
    td = Timedelta(1, unit).as_unit(unit)

    attrname = {"s": "second", "ms": "millisecond", "us": "microsecond"}.get(
        unit, "nanosecond"
    )
    msg = f"Out of bounds {attrname} timedelta: {-(2**63)}"
    with pytest.raises(OutOfBoundsTimedelta, match=msg):
        td * factor
    with pytest.raises(OutOfBoundsTimedelta, match=msg):
        factor * td


@pytest.mark.parametrize("unit", ["s", "ms", "us", "ns"])
@pytest.mark.parametrize("dtype", ["f8", "f4"])
def test_td_div_float_ndarray_overflow(unit, dtype):
    # GH#66552 a quotient outside the int64 range used to saturate on the cast,
    #  where the TimedeltaIndex equivalent raises
    td = Timedelta(4, unit)
    other = np.array([1e-30], dtype=dtype)

    msg = "Overflow in timedelta division"
    with pytest.raises(OutOfBoundsTimedelta, match=msg):
        td / other
    with pytest.raises(OutOfBoundsTimedelta, match=msg):
        td // other


@pytest.mark.parametrize("other", [np.array([0]), np.array([0.0]), np.array([np.nan])])
def test_td_div_ndarray_zero_or_nan_still_nat(other):
    # GH#66552 a zero or nan divisor is not an overflow: numpy calls it NaT and
    #  TimedeltaIndex keeps that, so the overflow check has to exempt it
    td = Timedelta(4, "ns")
    nat = np.array([np.timedelta64("NaT", "ns")])

    # a float zero divisor also trips numpy's own zero-division warning, which
    #  the TimedeltaIndex path suppresses and this one does not
    with np.errstate(divide="ignore", invalid="ignore"):
        tm.assert_numpy_array_equal(td / other, nat)
        tm.assert_numpy_array_equal(td // other, nat)


def test_td_div_ndarray_in_bounds():
    # GH#66552 the overflow check leaves representable quotients alone,
    #  including the ndim > 1 and empty cases
    td = Timedelta(12, "ns")

    expected = np.array([[6], [4]], dtype="m8[ns]")
    tm.assert_numpy_array_equal(td / np.array([[2.0], [3.0]]), expected)
    tm.assert_numpy_array_equal(td // np.array([[2.0], [3.0]]), expected)

    empty = np.array([], dtype="m8[ns]")
    tm.assert_numpy_array_equal(td / np.array([], dtype="f8"), empty)
    tm.assert_numpy_array_equal(td // np.array([], dtype="f8"), empty)


@pytest.mark.parametrize("unit", ["s", "ms", "us", "ns"])
def test_td_add_sub_lands_on_nat_sentinel(unit):
    # GH#66552 stepping one unit past Timedelta.min lands on iNaT, which is not
    #  NaT but is indistinguishable from it once stored, so it has to raise
    #  rather than come back as NaT.
    td_min = Timedelta(np.timedelta64(-(2**63) + 1, unit))

    attrname = {"s": "second", "ms": "millisecond", "us": "microsecond"}.get(
        unit, "nanosecond"
    )
    msg = f"Out of bounds {attrname} timedelta: {-(2**63)}"
    with pytest.raises(OutOfBoundsTimedelta, match=msg):
        td_min - Timedelta(1, unit)
    with pytest.raises(OutOfBoundsTimedelta, match=msg):
        td_min + Timedelta(-1, unit)
    with pytest.raises(OutOfBoundsTimedelta, match=msg):
        Timedelta(-1, unit) + td_min
    with pytest.raises(OutOfBoundsTimedelta, match=msg):
        td_min + np.timedelta64(-1, unit)


@pytest.mark.parametrize("op", [operator.mul, operator.truediv, operator.floordiv])
def test_td_integral_float_op_is_exact(op):
    # GH#66551 an integral float operand has an exact int equivalent, so it
    #  takes the int64 path rather than float64.  Timedelta.min is iNaT + 1,
    #  which float64 used to round down onto iNaT, and Timedelta.max used to
    #  round up out of the int64 range.
    assert op(Timedelta.min, 1.0) == Timedelta.min
    assert op(Timedelta.max, 1.0) == Timedelta.max

    # above 2**53 the float64 mantissa runs out, so the operand had to be
    #  applied exactly for these to round-trip
    td = Timedelta(2**53 + 1, "ns")
    assert op(td, 1.0)._value == 2**53 + 1
    assert op(td, 1.0) == op(td, 1)


def assert_masked_array_equal(result, expected):
    # what a MaskedArray holds underneath a masked slot is unspecified, so only
    #  the mask itself and the visible values are meaningful
    assert isinstance(result, np.ma.MaskedArray)
    mask = np.ma.getmaskarray(expected)
    tm.assert_numpy_array_equal(np.ma.getmaskarray(result), mask)
    tm.assert_numpy_array_equal(result.data[~mask], expected.data[~mask])


@pytest.mark.parametrize("kind", ["m", "M"])
def test_td_add_sub_ndarray_subclass(kind):
    # GH#66552 the overflow guard views the operand as i8 and rebuilds a plain
    #  ndarray, which drops an ndarray subclass' own semantics; a MaskedArray
    #  used to come back unmasked, exposing the payload under the mask
    td = Timedelta(5, "ns")
    other = np.ma.MaskedArray(
        np.array([10, 20], dtype=f"{kind}8[ns]"), mask=[False, True]
    )
    m8 = td.to_timedelta64()

    assert_masked_array_equal(td + other, m8 + other)
    assert_masked_array_equal(other + td, other + m8)
    assert_masked_array_equal(other - td, other - m8)
    if kind == "m":
        assert_masked_array_equal(td - other, m8 - other)


@pytest.mark.parametrize("dtype", ["i8", "f8"])
def test_td_mul_ndarray_subclass(dtype):
    # GH#66552 the float overflow guard reduces with np.max(..., initial=0.0),
    #  which dispatches to a subclass' own max() and raised there
    td = Timedelta(4, "ns")
    other = np.ma.MaskedArray(np.array([2, 3], dtype=dtype), mask=[False, True])
    m8 = td.to_timedelta64()

    assert_masked_array_equal(td * other, other * m8)
    assert_masked_array_equal(other * td, other * m8)


@pytest.mark.filterwarnings("ignore:__array_wrap__:DeprecationWarning")
def test_td_div_ndarray_subclass():
    # GH#66552 the overflow guard reduces with np.max(..., initial=, where=),
    #  which dispatches to a subclass' own max() and raised there.
    # NB: numpy's own MaskedArray.__array_wrap__ is broken on this path, with
    #  an m8 divisor or not: it drops the mask and emits the DeprecationWarning
    #  filtered above. So the expected values are not meaningful on their own;
    #  they are taken from numpy so that what is pinned is only that Timedelta
    #  does not diverge from timedelta64.
    td = Timedelta(4, "ns")
    other = np.ma.MaskedArray(np.array([2.0, 4.0]), mask=[False, True])
    m8 = td.to_timedelta64()

    assert_masked_array_equal(td / other, m8 / other)
    assert_masked_array_equal(td // other, m8 // other)


@pytest.mark.parametrize("kind", ["m", "M"])
def test_td_add_sub_ndarray_unit_multiplier(kind):
    # GH#66552 a dtype such as m8[10s] carries a multiplier that our
    #  resolutions cannot express; it used to be silently read as m8[s]
    td = Timedelta(1, "s")
    other = np.array([1, 2], dtype=f"{kind}8[10s]")

    expected = np.array([11, 21], dtype=f"{kind}8[s]")
    tm.assert_numpy_array_equal(td + other, expected)
    tm.assert_numpy_array_equal(other + td, expected)

    expected = np.array([9, 19], dtype=f"{kind}8[s]")
    tm.assert_numpy_array_equal(other - td, expected)

    if kind == "m":
        expected = np.array([-9, -19], dtype="m8[s]")
        tm.assert_numpy_array_equal(td - other, expected)
