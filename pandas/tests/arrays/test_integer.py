import numpy as np
import pytest

import pandas.util._test_decorators as td

from pandas.core.dtypes.generic import ABCIndexClass

import pandas as pd
import pandas._testing as tm
from pandas.api.types import is_float, is_float_dtype, is_integer, is_scalar
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


def make_data():
    return list(range(8)) + [np.nan] + list(range(10, 98)) + [np.nan] + [99, 100]


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
    return integer_array(make_data(), dtype=dtype)


@pytest.fixture
def data_missing(dtype):
    return integer_array([np.nan, 1], dtype=dtype)


@pytest.fixture(params=["data", "data_missing"])
def all_data(request, data, data_missing):
    """Parametrized fixture giving 'data' and 'data_missing'"""
    if request.param == "data":
        return data
    elif request.param == "data_missing":
        return data_missing


def test_dtypes(dtype):
    # smoke tests on auto dtype construction

    if dtype.is_signed_integer:
        assert np.dtype(dtype.type).kind == "i"
    else:
        assert np.dtype(dtype.type).kind == "u"
    assert dtype.name is not None


@pytest.mark.parametrize(
    "dtype, expected",
    [
        (Int8Dtype(), "Int8Dtype()"),
        (Int16Dtype(), "Int16Dtype()"),
        (Int32Dtype(), "Int32Dtype()"),
        (Int64Dtype(), "Int64Dtype()"),
        (UInt8Dtype(), "UInt8Dtype()"),
        (UInt16Dtype(), "UInt16Dtype()"),
        (UInt32Dtype(), "UInt32Dtype()"),
        (UInt64Dtype(), "UInt64Dtype()"),
    ],
)
def test_repr_dtype(dtype, expected):
    assert repr(dtype) == expected


def test_repr_array():
    result = repr(integer_array([1, None, 3]))
    expected = "<IntegerArray>\n[1, <NA>, 3]\nLength: 3, dtype: Int64"
    assert result == expected


def test_repr_array_long():
    data = integer_array([1, 2, None] * 1000)
    expected = (
        "<IntegerArray>\n"
        "[   1,    2, <NA>,    1,    2, <NA>,    1,    2, <NA>,    1,\n"
        " ...\n"
        " <NA>,    1,    2, <NA>,    1,    2, <NA>,    1,    2, <NA>]\n"
        "Length: 3000, dtype: Int64"
    )
    result = repr(data)
    assert result == expected


class TestConstructors:
    def test_uses_pandas_na(self):
        a = pd.array([1, None], dtype=pd.Int64Dtype())
        assert a[1] is pd.NA

    def test_from_dtype_from_float(self, data):
        # construct from our dtype & string dtype
        dtype = data.dtype

        # from float
        expected = pd.Series(data)
        result = pd.Series(
            data.to_numpy(na_value=np.nan, dtype="float"), dtype=str(dtype)
        )
        tm.assert_series_equal(result, expected)

        # from int / list
        expected = pd.Series(data)
        result = pd.Series(np.array(data).tolist(), dtype=str(dtype))
        tm.assert_series_equal(result, expected)

        # from int / array
        expected = pd.Series(data).dropna().reset_index(drop=True)
        dropped = np.array(data.dropna()).astype(np.dtype((dtype.type)))
        result = pd.Series(dropped, dtype=str(dtype))
        tm.assert_series_equal(result, expected)


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
        with pytest.raises(TypeError):
            ops("foo")
        with pytest.raises(TypeError):
            ops(pd.Timestamp("20180101"))

        # invalid array-likes
        with pytest.raises(TypeError):
            ops(pd.Series("foo", index=s.index))

        if op != "__rpow__":
            # TODO(extension)
            # rpow with a datetimelike coerces the integer array incorrectly
            with pytest.raises(TypeError):
                ops(pd.Series(pd.date_range("20180101", periods=len(s))))

        # 2d
        result = opa(pd.DataFrame({"A": s}))
        assert result is NotImplemented

        with pytest.raises(NotImplementedError):
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
        a = pd.array([0, 1, None, 2], dtype="Int64")
        result = a ** 0
        expected = pd.array([1, 1, 1, 1], dtype="Int64")
        tm.assert_extension_array_equal(result, expected)

        result = a ** 1
        expected = pd.array([0, 1, None, 2], dtype="Int64")
        tm.assert_extension_array_equal(result, expected)

        result = a ** pd.NA
        expected = pd.array([None, 1, None, None], dtype="Int64")
        tm.assert_extension_array_equal(result, expected)

        result = a ** np.nan
        expected = np.array([np.nan, 1, np.nan, np.nan], dtype="float64")
        tm.assert_numpy_array_equal(result, expected)

        # reversed
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


class TestComparisonOps(BaseOpsUtil):
    def _compare_other(self, data, op_name, other):
        op = self.get_op_from_name(op_name)

        # array
        result = pd.Series(op(data, other))
        expected = pd.Series(op(data._data, other), dtype="boolean")

        # fill the nan locations
        expected[data._mask] = pd.NA

        tm.assert_series_equal(result, expected)

        # series
        s = pd.Series(data)
        result = op(s, other)

        expected = op(pd.Series(data._data), other)

        # fill the nan locations
        expected[data._mask] = pd.NA
        expected = expected.astype("boolean")

        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize("other", [True, False, pd.NA, -1, 0, 1])
    def test_scalar(self, other, all_compare_operators):
        op = self.get_op_from_name(all_compare_operators)
        a = pd.array([1, 0, None], dtype="Int64")

        result = op(a, other)

        if other is pd.NA:
            expected = pd.array([None, None, None], dtype="boolean")
        else:
            values = op(a._data, other)
            expected = pd.arrays.BooleanArray(values, a._mask, copy=True)
        tm.assert_extension_array_equal(result, expected)

        # ensure we haven't mutated anything inplace
        result[0] = pd.NA
        tm.assert_extension_array_equal(a, pd.array([1, 0, None], dtype="Int64"))

    def test_array(self, all_compare_operators):
        op = self.get_op_from_name(all_compare_operators)
        a = pd.array([0, 1, 2, None, None, None], dtype="Int64")
        b = pd.array([0, 1, None, 0, 1, None], dtype="Int64")

        result = op(a, b)
        values = op(a._data, b._data)
        mask = a._mask | b._mask

        expected = pd.arrays.BooleanArray(values, mask)
        tm.assert_extension_array_equal(result, expected)

        # ensure we haven't mutated anything inplace
        result[0] = pd.NA
        tm.assert_extension_array_equal(
            a, pd.array([0, 1, 2, None, None, None], dtype="Int64")
        )
        tm.assert_extension_array_equal(
            b, pd.array([0, 1, None, 0, 1, None], dtype="Int64")
        )

    def test_compare_with_booleanarray(self, all_compare_operators):
        op = self.get_op_from_name(all_compare_operators)
        a = pd.array([True, False, None] * 3, dtype="boolean")
        b = pd.array([0] * 3 + [1] * 3 + [None] * 3, dtype="Int64")
        other = pd.array([False] * 3 + [True] * 3 + [None] * 3, dtype="boolean")
        expected = op(a, other)
        result = op(a, b)
        tm.assert_extension_array_equal(result, expected)

    def test_no_shared_mask(self, data):
        result = data + 1
        assert np.shares_memory(result._mask, data._mask) is False

    def test_compare_to_string(self, any_nullable_int_dtype):
        # GH 28930
        s = pd.Series([1, None], dtype=any_nullable_int_dtype)
        result = s == "a"
        expected = pd.Series([False, pd.NA], dtype="boolean")

        self.assert_series_equal(result, expected)

    def test_compare_to_int(self, any_nullable_int_dtype, all_compare_operators):
        # GH 28930
        s1 = pd.Series([1, None, 3], dtype=any_nullable_int_dtype)
        s2 = pd.Series([1, None, 3], dtype="float")

        method = getattr(s1, all_compare_operators)
        result = method(2)

        method = getattr(s2, all_compare_operators)
        expected = method(2).astype("boolean")
        expected[s2.isna()] = pd.NA

        self.assert_series_equal(result, expected)


class TestCasting:
    @pytest.mark.parametrize("dropna", [True, False])
    def test_construct_index(self, all_data, dropna):
        # ensure that we do not coerce to Float64Index, rather
        # keep as Index

        all_data = all_data[:10]
        if dropna:
            other = np.array(all_data[~all_data.isna()])
        else:
            other = all_data

        result = pd.Index(integer_array(other, dtype=all_data.dtype))
        expected = pd.Index(other, dtype=object)

        tm.assert_index_equal(result, expected)

    @pytest.mark.parametrize("dropna", [True, False])
    def test_astype_index(self, all_data, dropna):
        # as an int/uint index to Index

        all_data = all_data[:10]
        if dropna:
            other = all_data[~all_data.isna()]
        else:
            other = all_data

        dtype = all_data.dtype
        idx = pd.Index(np.array(other))
        assert isinstance(idx, ABCIndexClass)

        result = idx.astype(dtype)
        expected = idx.astype(object).astype(dtype)
        tm.assert_index_equal(result, expected)

    def test_astype(self, all_data):
        all_data = all_data[:10]

        ints = all_data[~all_data.isna()]
        mixed = all_data
        dtype = Int8Dtype()

        # coerce to same type - ints
        s = pd.Series(ints)
        result = s.astype(all_data.dtype)
        expected = pd.Series(ints)
        tm.assert_series_equal(result, expected)

        # coerce to same other - ints
        s = pd.Series(ints)
        result = s.astype(dtype)
        expected = pd.Series(ints, dtype=dtype)
        tm.assert_series_equal(result, expected)

        # coerce to same numpy_dtype - ints
        s = pd.Series(ints)
        result = s.astype(all_data.dtype.numpy_dtype)
        expected = pd.Series(ints._data.astype(all_data.dtype.numpy_dtype))
        tm.assert_series_equal(result, expected)

        # coerce to same type - mixed
        s = pd.Series(mixed)
        result = s.astype(all_data.dtype)
        expected = pd.Series(mixed)
        tm.assert_series_equal(result, expected)

        # coerce to same other - mixed
        s = pd.Series(mixed)
        result = s.astype(dtype)
        expected = pd.Series(mixed, dtype=dtype)
        tm.assert_series_equal(result, expected)

        # coerce to same numpy_dtype - mixed
        s = pd.Series(mixed)
        with pytest.raises(ValueError):
            s.astype(all_data.dtype.numpy_dtype)

        # coerce to object
        s = pd.Series(mixed)
        result = s.astype("object")
        expected = pd.Series(np.asarray(mixed))
        tm.assert_series_equal(result, expected)

    def test_astype_to_larger_numpy(self):
        a = pd.array([1, 2], dtype="Int32")
        result = a.astype("int64")
        expected = np.array([1, 2], dtype="int64")
        tm.assert_numpy_array_equal(result, expected)

        a = pd.array([1, 2], dtype="UInt32")
        result = a.astype("uint64")
        expected = np.array([1, 2], dtype="uint64")
        tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.parametrize("dtype", [Int8Dtype(), "Int8", UInt32Dtype(), "UInt32"])
    def test_astype_specific_casting(self, dtype):
        s = pd.Series([1, 2, 3], dtype="Int64")
        result = s.astype(dtype)
        expected = pd.Series([1, 2, 3], dtype=dtype)
        tm.assert_series_equal(result, expected)

        s = pd.Series([1, 2, 3, None], dtype="Int64")
        result = s.astype(dtype)
        expected = pd.Series([1, 2, 3, None], dtype=dtype)
        tm.assert_series_equal(result, expected)

    def test_construct_cast_invalid(self, dtype):

        msg = "cannot safely"
        arr = [1.2, 2.3, 3.7]
        with pytest.raises(TypeError, match=msg):
            integer_array(arr, dtype=dtype)

        with pytest.raises(TypeError, match=msg):
            pd.Series(arr).astype(dtype)

        arr = [1.2, 2.3, 3.7, np.nan]
        with pytest.raises(TypeError, match=msg):
            integer_array(arr, dtype=dtype)

        with pytest.raises(TypeError, match=msg):
            pd.Series(arr).astype(dtype)

    @pytest.mark.parametrize("in_series", [True, False])
    def test_to_numpy_na_nan(self, in_series):
        a = pd.array([0, 1, None], dtype="Int64")
        if in_series:
            a = pd.Series(a)

        result = a.to_numpy(dtype="float64", na_value=np.nan)
        expected = np.array([0.0, 1.0, np.nan], dtype="float64")
        tm.assert_numpy_array_equal(result, expected)

        result = a.to_numpy(dtype="int64", na_value=-1)
        expected = np.array([0, 1, -1], dtype="int64")
        tm.assert_numpy_array_equal(result, expected)

        result = a.to_numpy(dtype="bool", na_value=False)
        expected = np.array([False, True, False], dtype="bool")
        tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.parametrize("in_series", [True, False])
    @pytest.mark.parametrize("dtype", ["int32", "int64", "bool"])
    def test_to_numpy_dtype(self, dtype, in_series):
        a = pd.array([0, 1], dtype="Int64")
        if in_series:
            a = pd.Series(a)

        result = a.to_numpy(dtype=dtype)
        expected = np.array([0, 1], dtype=dtype)
        tm.assert_numpy_array_equal(result, expected)

    @pytest.mark.parametrize("dtype", ["float64", "int64", "bool"])
    def test_to_numpy_na_raises(self, dtype):
        a = pd.array([0, 1, None], dtype="Int64")
        with pytest.raises(ValueError, match=dtype):
            a.to_numpy(dtype=dtype)

    def test_astype_str(self):
        a = pd.array([1, 2, None], dtype="Int64")
        expected = np.array(["1", "2", "<NA>"], dtype=object)

        tm.assert_numpy_array_equal(a.astype(str), expected)
        tm.assert_numpy_array_equal(a.astype("str"), expected)


def test_frame_repr(data_missing):

    df = pd.DataFrame({"A": data_missing})
    result = repr(df)
    expected = "      A\n0  <NA>\n1     1"
    assert result == expected


def test_conversions(data_missing):

    # astype to object series
    df = pd.DataFrame({"A": data_missing})
    result = df["A"].astype("object")
    expected = pd.Series(np.array([np.nan, 1], dtype=object), name="A")
    tm.assert_series_equal(result, expected)

    # convert to object ndarray
    # we assert that we are exactly equal
    # including type conversions of scalars
    result = df["A"].astype("object").values
    expected = np.array([pd.NA, 1], dtype=object)
    tm.assert_numpy_array_equal(result, expected)

    for r, e in zip(result, expected):
        if pd.isnull(r):
            assert pd.isnull(e)
        elif is_integer(r):
            assert r == e
            assert is_integer(e)
        else:
            assert r == e
            assert type(r) == type(e)


def test_integer_array_constructor():
    values = np.array([1, 2, 3, 4], dtype="int64")
    mask = np.array([False, False, False, True], dtype="bool")

    result = IntegerArray(values, mask)
    expected = integer_array([1, 2, 3, np.nan], dtype="int64")
    tm.assert_extension_array_equal(result, expected)

    with pytest.raises(TypeError):
        IntegerArray(values.tolist(), mask)

    with pytest.raises(TypeError):
        IntegerArray(values, mask.tolist())

    with pytest.raises(TypeError):
        IntegerArray(values.astype(float), mask)

    with pytest.raises(TypeError):
        IntegerArray(values)


@pytest.mark.parametrize(
    "a, b",
    [
        ([1, None], [1, np.nan]),
        ([None], [np.nan]),
        ([None, np.nan], [np.nan, np.nan]),
        ([np.nan, np.nan], [np.nan, np.nan]),
    ],
)
def test_integer_array_constructor_none_is_nan(a, b):
    result = integer_array(a)
    expected = integer_array(b)
    tm.assert_extension_array_equal(result, expected)


def test_integer_array_constructor_copy():
    values = np.array([1, 2, 3, 4], dtype="int64")
    mask = np.array([False, False, False, True], dtype="bool")

    result = IntegerArray(values, mask)
    assert result._data is values
    assert result._mask is mask

    result = IntegerArray(values, mask, copy=True)
    assert result._data is not values
    assert result._mask is not mask


@pytest.mark.parametrize(
    "values",
    [
        ["foo", "bar"],
        ["1", "2"],
        "foo",
        1,
        1.0,
        pd.date_range("20130101", periods=2),
        np.array(["foo"]),
        [[1, 2], [3, 4]],
        [np.nan, {"a": 1}],
    ],
)
def test_to_integer_array_error(values):
    # error in converting existing arrays to IntegerArrays
    with pytest.raises(TypeError):
        integer_array(values)


def test_to_integer_array_inferred_dtype():
    # if values has dtype -> respect it
    result = integer_array(np.array([1, 2], dtype="int8"))
    assert result.dtype == Int8Dtype()
    result = integer_array(np.array([1, 2], dtype="int32"))
    assert result.dtype == Int32Dtype()

    # if values have no dtype -> always int64
    result = integer_array([1, 2])
    assert result.dtype == Int64Dtype()


def test_to_integer_array_dtype_keyword():
    result = integer_array([1, 2], dtype="int8")
    assert result.dtype == Int8Dtype()

    # if values has dtype -> override it
    result = integer_array(np.array([1, 2], dtype="int8"), dtype="int32")
    assert result.dtype == Int32Dtype()


def test_to_integer_array_float():
    result = integer_array([1.0, 2.0])
    expected = integer_array([1, 2])
    tm.assert_extension_array_equal(result, expected)

    with pytest.raises(TypeError, match="cannot safely cast non-equivalent"):
        integer_array([1.5, 2.0])

    # for float dtypes, the itemsize is not preserved
    result = integer_array(np.array([1.0, 2.0], dtype="float32"))
    assert result.dtype == Int64Dtype()


@pytest.mark.parametrize(
    "bool_values, int_values, target_dtype, expected_dtype",
    [
        ([False, True], [0, 1], Int64Dtype(), Int64Dtype()),
        ([False, True], [0, 1], "Int64", Int64Dtype()),
        ([False, True, np.nan], [0, 1, np.nan], Int64Dtype(), Int64Dtype()),
    ],
)
def test_to_integer_array_bool(bool_values, int_values, target_dtype, expected_dtype):
    result = integer_array(bool_values, dtype=target_dtype)
    assert result.dtype == expected_dtype
    expected = integer_array(int_values, dtype=target_dtype)
    tm.assert_extension_array_equal(result, expected)


@pytest.mark.parametrize(
    "values, to_dtype, result_dtype",
    [
        (np.array([1], dtype="int64"), None, Int64Dtype),
        (np.array([1, np.nan]), None, Int64Dtype),
        (np.array([1, np.nan]), "int8", Int8Dtype),
    ],
)
def test_to_integer_array(values, to_dtype, result_dtype):
    # convert existing arrays to IntegerArrays
    result = integer_array(values, dtype=to_dtype)
    assert result.dtype == result_dtype()
    expected = integer_array(values, dtype=result_dtype())
    tm.assert_extension_array_equal(result, expected)


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


@pytest.mark.parametrize("op", ["sum", "min", "max", "prod"])
def test_preserve_dtypes(op):
    # TODO(#22346): preserve Int64 dtype
    # for ops that enable (mean would actually work here
    # but generally it is a float return value)
    df = pd.DataFrame(
        {
            "A": ["a", "b", "b"],
            "B": [1, None, 3],
            "C": integer_array([1, None, 3], dtype="Int64"),
        }
    )

    # op
    result = getattr(df.C, op)()
    assert isinstance(result, int)

    # groupby
    result = getattr(df.groupby("A"), op)()

    expected = pd.DataFrame(
        {"B": np.array([1.0, 3.0]), "C": integer_array([1, 3], dtype="Int64")},
        index=pd.Index(["a", "b"], name="A"),
    )
    tm.assert_frame_equal(result, expected)


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


def test_astype_nansafe():
    # see gh-22343
    arr = integer_array([np.nan, 1, 2], dtype="Int8")
    msg = "cannot convert to 'uint32'-dtype NumPy array with missing values."

    with pytest.raises(ValueError, match=msg):
        arr.astype("uint32")


@pytest.mark.parametrize("ufunc", [np.abs, np.sign])
def test_ufuncs_single_int(ufunc):
    a = integer_array([1, 2, -3, np.nan])
    result = ufunc(a)
    expected = integer_array(ufunc(a.astype(float)))
    tm.assert_extension_array_equal(result, expected)

    s = pd.Series(a)
    result = ufunc(s)
    expected = pd.Series(integer_array(ufunc(a.astype(float))))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("ufunc", [np.log, np.exp, np.sin, np.cos, np.sqrt])
def test_ufuncs_single_float(ufunc):
    a = integer_array([1, 2, -3, np.nan])
    with np.errstate(invalid="ignore"):
        result = ufunc(a)
        expected = ufunc(a.astype(float))
    tm.assert_numpy_array_equal(result, expected)

    s = pd.Series(a)
    with np.errstate(invalid="ignore"):
        result = ufunc(s)
        expected = ufunc(s.astype(float))
    tm.assert_series_equal(result, expected)


@pytest.mark.parametrize("ufunc", [np.add, np.subtract])
def test_ufuncs_binary_int(ufunc):
    # two IntegerArrays
    a = integer_array([1, 2, -3, np.nan])
    result = ufunc(a, a)
    expected = integer_array(ufunc(a.astype(float), a.astype(float)))
    tm.assert_extension_array_equal(result, expected)

    # IntegerArray with numpy array
    arr = np.array([1, 2, 3, 4])
    result = ufunc(a, arr)
    expected = integer_array(ufunc(a.astype(float), arr))
    tm.assert_extension_array_equal(result, expected)

    result = ufunc(arr, a)
    expected = integer_array(ufunc(arr, a.astype(float)))
    tm.assert_extension_array_equal(result, expected)

    # IntegerArray with scalar
    result = ufunc(a, 1)
    expected = integer_array(ufunc(a.astype(float), 1))
    tm.assert_extension_array_equal(result, expected)

    result = ufunc(1, a)
    expected = integer_array(ufunc(1, a.astype(float)))
    tm.assert_extension_array_equal(result, expected)


@pytest.mark.parametrize("values", [[0, 1], [0, None]])
def test_ufunc_reduce_raises(values):
    a = integer_array(values)
    with pytest.raises(NotImplementedError):
        np.add.reduce(a)


@td.skip_if_no("pyarrow", min_version="0.15.0")
def test_arrow_array(data):
    # protocol added in 0.15.0
    import pyarrow as pa

    arr = pa.array(data)
    expected = np.array(data, dtype=object)
    expected[data.isna()] = None
    expected = pa.array(expected, type=data.dtype.name.lower(), from_pandas=True)
    assert arr.equals(expected)


@td.skip_if_no("pyarrow", min_version="0.15.1.dev")
def test_arrow_roundtrip(data):
    # roundtrip possible from arrow 1.0.0
    import pyarrow as pa

    df = pd.DataFrame({"a": data})
    table = pa.table(df)
    assert table.field("a").type == str(data.dtype.numpy_dtype)
    result = table.to_pandas()
    tm.assert_frame_equal(result, df)


@pytest.mark.parametrize(
    "pandasmethname, kwargs",
    [
        ("var", {"ddof": 0}),
        ("var", {"ddof": 1}),
        ("kurtosis", {}),
        ("skew", {}),
        ("sem", {}),
    ],
)
def test_stat_method(pandasmethname, kwargs):
    s = pd.Series(data=[1, 2, 3, 4, 5, 6, np.nan, np.nan], dtype="Int64")
    pandasmeth = getattr(s, pandasmethname)
    result = pandasmeth(**kwargs)
    s2 = pd.Series(data=[1, 2, 3, 4, 5, 6], dtype="Int64")
    pandasmeth = getattr(s2, pandasmethname)
    expected = pandasmeth(**kwargs)
    assert expected == result


# TODO(jreback) - these need testing / are broken

# shift

# set_index (destroys type)
