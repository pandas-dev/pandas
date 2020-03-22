import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm
from pandas.api.types import is_float, is_float_dtype, is_scalar
from pandas.core.arrays import IntegerArray, integer_array
from pandas.core.arrays.integer import (
    Int8Dtype,
    Int16Dtype,
    Int32Dtype,
    Int64Dtype,
    UInt8Dtype,
    UInt16Dtype,
    UInt32Dtype,
    UInt64Dtype,
)
from pandas.tests.extension.base import BaseOpsUtil


@pytest.fixture(
    params=[
        Int8Dtype,
        Int16Dtype,
        Int32Dtype,
        Int64Dtype,
        UInt8Dtype,
        UInt16Dtype,
        UInt32Dtype,
        UInt64Dtype,
    ]
)
def dtype(request):
    return request.param()


@pytest.fixture
def data(dtype):
    return integer_array(
        list(range(8)) + [np.nan] + list(range(10, 98)) + [np.nan] + [99, 100],
        dtype=dtype,
    )


class TestArithmeticOps(BaseOpsUtil):
    def _check_divmod_op(self, s, op, other, exc=None):
        super()._check_divmod_op(s, op, other, None)

    def _check_op(self, s, op_name, other, exc=None):
        op = self.get_op_from_name(op_name)
        result = op(s, other)

        # compute expected
        mask = s.isna()

        # if s is a DataFrame, squeeze to a Series
        # for comparison
        if isinstance(s, pd.DataFrame):
            result = result.squeeze()
            s = s.squeeze()
            mask = mask.squeeze()

        # other array is an Integer
        if isinstance(other, IntegerArray):
            omask = getattr(other, "mask", None)
            mask = getattr(other, "data", other)
            if omask is not None:
                mask |= omask

        # 1 ** na is na, so need to unmask those
        if op_name == "__pow__":
            mask = np.where(~s.isna() & (s == 1), False, mask)

        elif op_name == "__rpow__":
            other_is_one = other == 1
            if isinstance(other_is_one, pd.Series):
                other_is_one = other_is_one.fillna(False)
            mask = np.where(other_is_one, False, mask)

        # float result type or float op
        if (
            is_float_dtype(other)
            or is_float(other)
            or op_name in ["__rtruediv__", "__truediv__", "__rdiv__", "__div__"]
        ):
            rs = s.astype("float")
            expected = op(rs, other)
            self._check_op_float(result, expected, mask, s, op_name, other)

        # integer result type
        else:
            rs = pd.Series(s.values._data, name=s.name)
            expected = op(rs, other)
            self._check_op_integer(result, expected, mask, s, op_name, other)

    def _check_op_float(self, result, expected, mask, s, op_name, other):
        # check comparisons that are resulting in float dtypes

        expected[mask] = np.nan
        if "floordiv" in op_name:
            # Series op sets 1//0 to np.inf, which IntegerArray does not do (yet)
            mask2 = np.isinf(expected) & np.isnan(result)
            expected[mask2] = np.nan
        tm.assert_series_equal(result, expected)

    def _check_op_integer(self, result, expected, mask, s, op_name, other):
        # check comparisons that are resulting in integer dtypes

        # to compare properly, we convert the expected
        # to float, mask to nans and convert infs
        # if we have uints then we process as uints
        # then convert to float
        # and we ultimately want to create a IntArray
        # for comparisons

        fill_value = 0

        # mod/rmod turn floating 0 into NaN while
        # integer works as expected (no nan)
        if op_name in ["__mod__", "__rmod__"]:
            if is_scalar(other):
                if other == 0:
                    expected[s.values == 0] = 0
                else:
                    expected = expected.fillna(0)
            else:
                expected[
                    (s.values == 0).fillna(False)
                    & ((expected == 0).fillna(False) | expected.isna())
                ] = 0
        try:
            expected[
                ((expected == np.inf) | (expected == -np.inf)).fillna(False)
            ] = fill_value
            original = expected
            expected = expected.astype(s.dtype)

        except ValueError:

            expected = expected.astype(float)
            expected[
                ((expected == np.inf) | (expected == -np.inf)).fillna(False)
            ] = fill_value
            original = expected
            expected = expected.astype(s.dtype)

        expected[mask] = pd.NA

        # assert that the expected astype is ok
        # (skip for unsigned as they have wrap around)
        if not s.dtype.is_unsigned_integer:
            original = pd.Series(original)

            # we need to fill with 0's to emulate what an astype('int') does
            # (truncation) for certain ops
            if op_name in ["__rtruediv__", "__rdiv__"]:
                mask |= original.isna()
                original = original.fillna(0).astype("int")

            original = original.astype("float")
            original[mask] = np.nan
            tm.assert_series_equal(original, expected.astype("float"))

        # assert our expected result
        tm.assert_series_equal(result, expected)

    def test_arith_integer_array(self, data, all_arithmetic_operators):
        # we operate with a rhs of an integer array

        op = all_arithmetic_operators

        s = pd.Series(data)
        rhs = pd.Series([1] * len(data), dtype=data.dtype)
        rhs.iloc[-1] = np.nan

        self._check_op(s, op, rhs)

    def test_arith_series_with_scalar(self, data, all_arithmetic_operators):
        # scalar
        op = all_arithmetic_operators
        s = pd.Series(data)
        self._check_op(s, op, 1, exc=TypeError)

    def test_arith_frame_with_scalar(self, data, all_arithmetic_operators):
        # frame & scalar
        op = all_arithmetic_operators
        df = pd.DataFrame({"A": data})
        self._check_op(df, op, 1, exc=TypeError)

    def test_arith_series_with_array(self, data, all_arithmetic_operators):
        # ndarray & other series
        op = all_arithmetic_operators
        s = pd.Series(data)
        other = np.ones(len(s), dtype=s.dtype.type)
        self._check_op(s, op, other, exc=TypeError)

    def test_arith_coerce_scalar(self, data, all_arithmetic_operators):

        op = all_arithmetic_operators
        s = pd.Series(data)

        other = 0.01
        self._check_op(s, op, other)

    @pytest.mark.parametrize("other", [1.0, np.array(1.0)])
    def test_arithmetic_conversion(self, all_arithmetic_operators, other):
        # if we have a float operand we should have a float result
        # if that is equal to an integer
        op = self.get_op_from_name(all_arithmetic_operators)

        s = pd.Series([1, 2, 3], dtype="Int64")
        result = op(s, other)
        assert result.dtype is np.dtype("float")

    def test_arith_len_mismatch(self, all_arithmetic_operators):
        # operating with a list-like with non-matching length raises
        op = self.get_op_from_name(all_arithmetic_operators)
        other = np.array([1.0])

        s = pd.Series([1, 2, 3], dtype="Int64")
        with pytest.raises(ValueError, match="Lengths must match"):
            op(s, other)

    @pytest.mark.parametrize("other", [0, 0.5])
    def test_arith_zero_dim_ndarray(self, other):
        arr = integer_array([1, None, 2])
        result = arr + np.array(other)
        expected = arr + other
        tm.assert_equal(result, expected)

    def test_error(self, data, all_arithmetic_operators):
        # invalid ops

        op = all_arithmetic_operators
        s = pd.Series(data)
        ops = getattr(s, op)
        opa = getattr(data, op)

        # invalid scalars
        msg = (
            r"(:?can only perform ops with numeric values)"
            r"|(:?IntegerArray cannot perform the operation mod)"
        )
        with pytest.raises(TypeError, match=msg):
            ops("foo")
        with pytest.raises(TypeError, match=msg):
            ops(pd.Timestamp("20180101"))

        # invalid array-likes
        with pytest.raises(TypeError, match=msg):
            ops(pd.Series("foo", index=s.index))

        if op != "__rpow__":
            # TODO(extension)
            # rpow with a datetimelike coerces the integer array incorrectly
            msg = (
                "can only perform ops with numeric values|"
                "cannot perform .* with this index type: DatetimeArray|"
                "Addition/subtraction of integers and integer-arrays "
                "with DatetimeArray is no longer supported. *"
            )
            with pytest.raises(TypeError, match=msg):
                ops(pd.Series(pd.date_range("20180101", periods=len(s))))

        # 2d
        result = opa(pd.DataFrame({"A": s}))
        assert result is NotImplemented

        msg = r"can only perform ops with 1-d structures"
        with pytest.raises(NotImplementedError, match=msg):
            opa(np.arange(len(s)).reshape(-1, len(s)))

    @pytest.mark.parametrize("zero, negative", [(0, False), (0.0, False), (-0.0, True)])
    def test_divide_by_zero(self, zero, negative):
        # https://github.com/pandas-dev/pandas/issues/27398
        a = pd.array([0, 1, -1, None], dtype="Int64")
        result = a / zero
        expected = np.array([np.nan, np.inf, -np.inf, np.nan])
        if negative:
            expected *= -1
        tm.assert_numpy_array_equal(result, expected)

    def test_pow_scalar(self):
        a = pd.array([-1, 0, 1, None, 2], dtype="Int64")
        result = a ** 0
        expected = pd.array([1, 1, 1, 1, 1], dtype="Int64")
        tm.assert_extension_array_equal(result, expected)

        result = a ** 1
        expected = pd.array([-1, 0, 1, None, 2], dtype="Int64")
        tm.assert_extension_array_equal(result, expected)

        result = a ** pd.NA
        expected = pd.array([None, None, 1, None, None], dtype="Int64")
        tm.assert_extension_array_equal(result, expected)

        result = a ** np.nan
        expected = np.array([np.nan, np.nan, 1, np.nan, np.nan], dtype="float64")
        tm.assert_numpy_array_equal(result, expected)

        # reversed
        a = a[1:]  # Can't raise integers to negative powers.

        result = 0 ** a
        expected = pd.array([1, 0, None, 0], dtype="Int64")
        tm.assert_extension_array_equal(result, expected)

        result = 1 ** a
        expected = pd.array([1, 1, 1, 1], dtype="Int64")
        tm.assert_extension_array_equal(result, expected)

        result = pd.NA ** a
        expected = pd.array([1, None, None, None], dtype="Int64")
        tm.assert_extension_array_equal(result, expected)

        result = np.nan ** a
        expected = np.array([1, np.nan, np.nan, np.nan], dtype="float64")
        tm.assert_numpy_array_equal(result, expected)

    def test_pow_array(self):
        a = integer_array([0, 0, 0, 1, 1, 1, None, None, None])
        b = integer_array([0, 1, None, 0, 1, None, 0, 1, None])
        result = a ** b
        expected = integer_array([1, 0, None, 1, 1, 1, 1, None, None])
        tm.assert_extension_array_equal(result, expected)

    def test_rpow_one_to_na(self):
        # https://github.com/pandas-dev/pandas/issues/22022
        # https://github.com/pandas-dev/pandas/issues/29997
        arr = integer_array([np.nan, np.nan])
        result = np.array([1.0, 2.0]) ** arr
        expected = np.array([1.0, np.nan])
        tm.assert_numpy_array_equal(result, expected)


def test_cross_type_arithmetic():

    df = pd.DataFrame(
        {
            "A": pd.Series([1, 2, np.nan], dtype="Int64"),
            "B": pd.Series([1, np.nan, 3], dtype="UInt8"),
            "C": [1, 2, 3],
        }
    )

    result = df.A + df.C
    expected = pd.Series([2, 4, np.nan], dtype="Int64")
    tm.assert_series_equal(result, expected)

    result = (df.A + df.C) * 3 == 12
    expected = pd.Series([False, True, None], dtype="boolean")
    tm.assert_series_equal(result, expected)

    result = df.A + df.B
    expected = pd.Series([2, np.nan, np.nan], dtype="Int64")
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("op", ["mean"])
def test_reduce_to_float(op):
    # some reduce ops always return float, even if the result
    # is a rounded number
    df = pd.DataFrame(
        {
            "A": ["a", "b", "b"],
            "B": [1, None, 3],
            "C": integer_array([1, None, 3], dtype="Int64"),
        }
    )

    # op
    result = getattr(df.C, op)()
    assert isinstance(result, float)

    # groupby
    result = getattr(df.groupby("A"), op)()

    expected = pd.DataFrame(
        {"B": np.array([1.0, 3.0]), "C": integer_array([1, 3], dtype="Int64")},
        index=pd.Index(["a", "b"], name="A"),
    )
    tm.assert_frame_equal(result, expected)
