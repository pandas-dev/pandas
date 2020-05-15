import numpy as np
import pytest

import pandas as pd
import pandas._testing as tm
from pandas.core.arrays import FloatingArray, IntegerArray
from pandas.tests.extension.base import BaseOpsUtil


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
        if isinstance(other, (IntegerArray, FloatingArray)):
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

        rs = s.astype(s.dtype.numpy_dtype)
        expected = op(rs, other).astype(s.dtype)
        expected[mask] = np.nan
        if "floordiv" in op_name:
            # Series op sets 1//0 to np.inf, which IntegerArray does not do (yet)
            mask2 = np.isinf(expected) & np.isnan(result)
            expected[mask2] = np.nan
        tm.assert_series_equal(result, expected)

    def test_arith_floating_array(self, data, all_arithmetic_operators):
        # we operate with a rhs of an floating array

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

    def test_arith_len_mismatch(self, all_arithmetic_operators):
        # operating with a list-like with non-matching length raises
        op = self.get_op_from_name(all_arithmetic_operators)
        other = np.array([1.0])

        s = pd.Series([1, 2, 3], dtype="Float64")
        with pytest.raises(ValueError, match="Lengths must match"):
            op(s, other)

    @pytest.mark.parametrize("other", [0, 0.5])
    def test_arith_zero_dim_ndarray(self, other):
        arr = pd.array([1, None, 2], dtype="Float64")
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
            r"|(:?FloatingArray cannot perform the operation mod)"
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
        # TODO pending NA/NaN discussion
        # https://github.com/pandas-dev/pandas/issues/32265/
        a = pd.array([0, 1, -1, None], dtype="Float64")
        result = a / zero
        expected = FloatingArray(
            np.array([np.nan, np.inf, -np.inf, np.nan]),
            np.array([False, False, False, True]),
        )
        if negative:
            expected *= -1
        tm.assert_extension_array_equal(result, expected)

    def test_pow_scalar(self):
        a = pd.array([-1, 0, 1, None, 2], dtype="Float64")
        result = a ** 0
        expected = pd.array([1, 1, 1, 1, 1], dtype="Float64")
        tm.assert_extension_array_equal(result, expected)

        result = a ** 1
        expected = pd.array([-1, 0, 1, None, 2], dtype="Float64")
        tm.assert_extension_array_equal(result, expected)

        result = a ** pd.NA
        expected = pd.array([None, None, 1, None, None], dtype="Float64")
        tm.assert_extension_array_equal(result, expected)

        result = a ** np.nan
        # TODO np.nan should be converted to pd.NA / missing before operation?
        expected = FloatingArray(
            np.array([np.nan, np.nan, 1, np.nan, np.nan], dtype="float64"), mask=a._mask
        )
        tm.assert_extension_array_equal(result, expected)

        # reversed
        a = a[1:]  # Can't raise integers to negative powers.

        result = 0 ** a
        expected = pd.array([1, 0, None, 0], dtype="Float64")
        tm.assert_extension_array_equal(result, expected)

        result = 1 ** a
        expected = pd.array([1, 1, 1, 1], dtype="Float64")
        tm.assert_extension_array_equal(result, expected)

        result = pd.NA ** a
        expected = pd.array([1, None, None, None], dtype="Float64")
        tm.assert_extension_array_equal(result, expected)

        result = np.nan ** a
        expected = FloatingArray(
            np.array([1, np.nan, np.nan, np.nan], dtype="float64"), mask=a._mask
        )
        tm.assert_extension_array_equal(result, expected)

    def test_pow_array(self):
        a = pd.array([0, 0, 0, 1, 1, 1, None, None, None], dtype="Float64")
        b = pd.array([0, 1, None, 0, 1, None, 0, 1, None], dtype="Float64")
        result = a ** b
        expected = pd.array([1, 0, None, 1, 1, 1, 1, None, None], dtype="Float64")
        tm.assert_extension_array_equal(result, expected)

    def test_rpow_one_to_na(self):
        # https://github.com/pandas-dev/pandas/issues/22022
        # https://github.com/pandas-dev/pandas/issues/29997
        arr = pd.array([np.nan, np.nan], dtype="Float64")
        result = np.array([1.0, 2.0]) ** arr
        expected = pd.array([1.0, np.nan], dtype="Float64")
        tm.assert_extension_array_equal(result, expected)


def test_cross_type_arithmetic():

    df = pd.DataFrame(
        {
            "A": pd.array([1, 2, np.nan], dtype="Float64"),
            "B": pd.array([1, np.nan, 3], dtype="Float32"),
            "C": np.array([1, 2, 3], dtype="float64"),
        }
    )

    result = df.A + df.C
    expected = pd.Series([2, 4, np.nan], dtype="Float64")
    tm.assert_series_equal(result, expected)

    result = (df.A + df.C) * 3 == 12
    expected = pd.Series([False, True, None], dtype="boolean")
    tm.assert_series_equal(result, expected)

    result = df.A + df.B
    expected = pd.Series([2, np.nan, np.nan], dtype="Float64")
    tm.assert_series_equal(result, expected)
