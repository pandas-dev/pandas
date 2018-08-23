# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import pandas.util.testing as tm
import pytest

from pandas.api.types import is_integer, is_float, is_float_dtype, is_scalar
from pandas.core.dtypes.generic import ABCIndexClass

from pandas.core.arrays import (
    integer_array, IntegerArray)
from pandas.core.arrays.integer import (
    Int8Dtype, Int16Dtype, Int32Dtype, Int64Dtype,
    UInt8Dtype, UInt16Dtype, UInt32Dtype, UInt64Dtype)

from ..extension.base import BaseOpsUtil


def make_data():
    return (list(range(8)) +
            [np.nan] +
            list(range(10, 98)) +
            [np.nan] +
            [99, 100])


@pytest.fixture(params=[Int8Dtype, Int16Dtype, Int32Dtype, Int64Dtype,
                        UInt8Dtype, UInt16Dtype, UInt32Dtype, UInt64Dtype])
def dtype(request):
    return request.param()


@pytest.fixture
def data(dtype):
    return integer_array(make_data(), dtype=dtype)


@pytest.fixture
def data_missing(dtype):
    return integer_array([np.nan, 1], dtype=dtype)


@pytest.fixture(params=['data', 'data_missing'])
def all_data(request, data, data_missing):
    """Parametrized fixture giving 'data' and 'data_missing'"""
    if request.param == 'data':
        return data
    elif request.param == 'data_missing':
        return data_missing


def test_dtypes(dtype):
    # smoke tests on auto dtype construction

    if dtype.is_signed_integer:
        assert np.dtype(dtype.type).kind == 'i'
    else:
        assert np.dtype(dtype.type).kind == 'u'
    assert dtype.name is not None


class TestInterface(object):

    def test_repr_array(self, data):
        result = repr(data)

        # not long
        assert '...' not in result

        assert 'dtype=' in result
        assert 'IntegerArray' in result

    def test_repr_array_long(self, data):
        # some arrays may be able to assert a ... in the repr
        with pd.option_context('display.max_seq_items', 1):
            result = repr(data)

            assert '...' in result
            assert 'length' in result


class TestConstructors(object):

    def test_from_dtype_from_float(self, data):
        # construct from our dtype & string dtype
        dtype = data.dtype

        # from float
        expected = pd.Series(data)
        result = pd.Series(np.array(data).astype('float'), dtype=str(dtype))
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
        super(TestArithmeticOps, self)._check_divmod_op(s, op, other, None)

    def _check_op(self, s, op_name, other, exc=None):
        op = self.get_op_from_name(op_name)
        result = op(s, other)

        # compute expected
        mask = s.isna()

        # other array is an Integer
        if isinstance(other, IntegerArray):
            omask = getattr(other, 'mask', None)
            mask = getattr(other, 'data', other)
            if omask is not None:
                mask |= omask

        # float result type or float op
        if ((is_float_dtype(other) or is_float(other) or
             op_name in ['__rtruediv__', '__truediv__',
                         '__rdiv__', '__div__'])):
            rs = s.astype('float')
            expected = op(rs, other)
            self._check_op_float(result, expected, mask, s, op_name, other)

        # integer result type
        else:
            rs = pd.Series(s.values._data)
            expected = op(rs, other)
            self._check_op_integer(result, expected, mask, s, op_name, other)

    def _check_op_float(self, result, expected, mask, s, op_name, other):
        # check comparisions that are resulting in float dtypes

        expected[mask] = np.nan
        tm.assert_series_equal(result, expected)

    def _check_op_integer(self, result, expected, mask, s, op_name, other):
        # check comparisions that are resulting in integer dtypes

        # to compare properly, we convert the expected
        # to float, mask to nans and convert infs
        # if we have uints then we process as uints
        # then conert to float
        # and we ultimately want to create a IntArray
        # for comparisons

        fill_value = 0

        # mod/rmod turn floating 0 into NaN while
        # integer works as expected (no nan)
        if op_name in ['__mod__', '__rmod__']:
            if is_scalar(other):
                if other == 0:
                    expected[s.values == 0] = 0
                else:
                    expected = expected.fillna(0)
            else:
                expected[(s.values == 0) &
                         ((expected == 0) | expected.isna())] = 0

        try:
            expected[(expected == np.inf) | (expected == -np.inf)] = fill_value
            original = expected
            expected = expected.astype(s.dtype)

        except ValueError:

            expected = expected.astype(float)
            expected[(expected == np.inf) | (expected == -np.inf)] = fill_value
            original = expected
            expected = expected.astype(s.dtype)

        expected[mask] = np.nan

        # assert that the expected astype is ok
        # (skip for unsigned as they have wrap around)
        if not s.dtype.is_unsigned_integer:
            original = pd.Series(original)

            # we need to fill with 0's to emulate what an astype('int') does
            # (truncation) for certain ops
            if op_name in ['__rtruediv__', '__rdiv__']:
                mask |= original.isna()
                original = original.fillna(0).astype('int')

            original = original.astype('float')
            original[mask] = np.nan
            tm.assert_series_equal(original, expected.astype('float'))

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

    @pytest.mark.xfail(run=False, reason="_reduce needs implementation")
    def test_arith_frame_with_scalar(self, data, all_arithmetic_operators):
        # frame & scalar
        op = all_arithmetic_operators

        df = pd.DataFrame({'A': data})
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

    @pytest.mark.parametrize("other", [1., 1.0, np.array(1.), np.array([1.])])
    def test_arithmetic_conversion(self, all_arithmetic_operators, other):
        # if we have a float operand we should have a float result
        # if if that is equal to an integer
        op = self.get_op_from_name(all_arithmetic_operators)

        s = pd.Series([1, 2, 3], dtype='Int64')
        result = op(s, other)
        assert result.dtype is np.dtype('float')

    def test_error(self, data, all_arithmetic_operators):
        # invalid ops

        op = all_arithmetic_operators
        s = pd.Series(data)
        ops = getattr(s, op)
        opa = getattr(data, op)

        # invalid scalars
        with pytest.raises(TypeError):
            ops('foo')
        with pytest.raises(TypeError):
            ops(pd.Timestamp('20180101'))

        # invalid array-likes
        with pytest.raises(TypeError):
            ops(pd.Series('foo', index=s.index))

        if op != '__rpow__':
            # TODO(extension)
            # rpow with a datetimelike coerces the integer array incorrectly
            with pytest.raises(TypeError):
                ops(pd.Series(pd.date_range('20180101', periods=len(s))))

        # 2d
        with pytest.raises(NotImplementedError):
            opa(pd.DataFrame({'A': s}))
        with pytest.raises(NotImplementedError):
            opa(np.arange(len(s)).reshape(-1, len(s)))


class TestComparisonOps(BaseOpsUtil):

    def _compare_other(self, s, data, op_name, other):
        op = self.get_op_from_name(op_name)

        # array
        result = op(s, other)
        expected = pd.Series(op(data._data, other))

        # fill the nan locations
        expected[data._mask] = True if op_name == '__ne__' else False

        tm.assert_series_equal(result, expected)

        # series
        s = pd.Series(data)
        result = op(s, other)

        expected = pd.Series(data._data)
        expected = op(expected, other)

        # fill the nan locations
        expected[data._mask] = True if op_name == '__ne__' else False

        tm.assert_series_equal(result, expected)

    def test_compare_scalar(self, data, all_compare_operators):
        op_name = all_compare_operators
        s = pd.Series(data)
        self._compare_other(s, data, op_name, 0)

    def test_compare_array(self, data, all_compare_operators):
        op_name = all_compare_operators
        s = pd.Series(data)
        other = pd.Series([0] * len(data))
        self._compare_other(s, data, op_name, other)


class TestCasting(object):
    pass

    @pytest.mark.parametrize('dropna', [True, False])
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

    @pytest.mark.parametrize('dropna', [True, False])
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
        expected = pd.Series(ints._data.astype(
            all_data.dtype.numpy_dtype))
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
        result = s.astype('object')
        expected = pd.Series(np.asarray(mixed))
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize('dtype', [Int8Dtype(), 'Int8',
                                       UInt32Dtype(), 'UInt32'])
    def test_astype_specific_casting(self, dtype):
        s = pd.Series([1, 2, 3], dtype='Int64')
        result = s.astype(dtype)
        expected = pd.Series([1, 2, 3], dtype=dtype)
        tm.assert_series_equal(result, expected)

        s = pd.Series([1, 2, 3, None], dtype='Int64')
        result = s.astype(dtype)
        expected = pd.Series([1, 2, 3, None], dtype=dtype)
        tm.assert_series_equal(result, expected)

    def test_construct_cast_invalid(self, dtype):

        msg = "cannot safely"
        arr = [1.2, 2.3, 3.7]
        with tm.assert_raises_regex(TypeError, msg):
            integer_array(arr, dtype=dtype)

        with tm.assert_raises_regex(TypeError, msg):
            pd.Series(arr).astype(dtype)

        arr = [1.2, 2.3, 3.7, np.nan]
        with tm.assert_raises_regex(TypeError, msg):
            integer_array(arr, dtype=dtype)

        with tm.assert_raises_regex(TypeError, msg):
            pd.Series(arr).astype(dtype)


def test_frame_repr(data_missing):

    df = pd.DataFrame({'A': data_missing})
    result = repr(df)
    expected = '     A\n0  NaN\n1    1'
    assert result == expected


def test_conversions(data_missing):

    # astype to object series
    df = pd.DataFrame({'A': data_missing})
    result = df['A'].astype('object')
    expected = pd.Series(np.array([np.nan, 1], dtype=object), name='A')
    tm.assert_series_equal(result, expected)

    # convert to object ndarray
    # we assert that we are exactly equal
    # including type conversions of scalars
    result = df['A'].astype('object').values
    expected = np.array([np.nan, 1], dtype=object)
    tm.assert_numpy_array_equal(result, expected)

    for r, e in zip(result, expected):
        if pd.isnull(r):
            assert pd.isnull(e)
        elif is_integer(r):
            # PY2 can be int or long
            assert r == e
            assert is_integer(e)
        else:
            assert r == e
            assert type(r) == type(e)


def test_integer_array_constructor():
    values = np.array([1, 2, 3, 4], dtype='int64')
    mask = np.array([False, False, False, True], dtype='bool')

    result = IntegerArray(values, mask)
    expected = integer_array([1, 2, 3, np.nan], dtype='int64')
    tm.assert_extension_array_equal(result, expected)

    with pytest.raises(TypeError):
        IntegerArray(values.tolist(), mask)

    with pytest.raises(TypeError):
        IntegerArray(values, mask.tolist())

    with pytest.raises(TypeError):
        IntegerArray(values.astype(float), mask)

    with pytest.raises(TypeError):
        IntegerArray(values)


def test_integer_array_constructor_copy():
    values = np.array([1, 2, 3, 4], dtype='int64')
    mask = np.array([False, False, False, True], dtype='bool')

    result = IntegerArray(values, mask)
    assert result._data is values
    assert result._mask is mask

    result = IntegerArray(values, mask, copy=True)
    assert result._data is not values
    assert result._mask is not mask


@pytest.mark.parametrize(
    'values',
    [
        ['foo', 'bar'],
        ['1', '2'],
        'foo',
        1,
        1.0,
        pd.date_range('20130101', periods=2),
        np.array(['foo'])])
def test_to_integer_array_error(values):
    # error in converting existing arrays to IntegerArrays
    with pytest.raises(TypeError):
        integer_array(values)


def test_to_integer_array_inferred_dtype():
    # if values has dtype -> respect it
    result = integer_array(np.array([1, 2], dtype='int8'))
    assert result.dtype == Int8Dtype()
    result = integer_array(np.array([1, 2], dtype='int32'))
    assert result.dtype == Int32Dtype()

    # if values have no dtype -> always int64
    result = integer_array([1, 2])
    assert result.dtype == Int64Dtype()


def test_to_integer_array_dtype_keyword():
    result = integer_array([1, 2], dtype='int8')
    assert result.dtype == Int8Dtype()

    # if values has dtype -> override it
    result = integer_array(np.array([1, 2], dtype='int8'), dtype='int32')
    assert result.dtype == Int32Dtype()


def test_to_integer_array_float():
    result = integer_array([1., 2.])
    expected = integer_array([1, 2])
    tm.assert_extension_array_equal(result, expected)

    with pytest.raises(TypeError, match="cannot safely cast non-equivalent"):
        integer_array([1.5, 2.])

    # for float dtypes, the itemsize is not preserved
    result = integer_array(np.array([1., 2.], dtype='float32'))
    assert result.dtype == Int64Dtype()


@pytest.mark.parametrize(
    'values, to_dtype, result_dtype',
    [
        (np.array([1], dtype='int64'), None, Int64Dtype),
        (np.array([1, np.nan]), None, Int64Dtype),
        (np.array([1, np.nan]), 'int8', Int8Dtype)])
def test_to_integer_array(values, to_dtype, result_dtype):
    # convert existing arrays to IntegerArrays
    result = integer_array(values, dtype=to_dtype)
    assert result.dtype == result_dtype()
    expected = integer_array(values, dtype=result_dtype())
    tm.assert_extension_array_equal(result, expected)


def test_cross_type_arithmetic():

    df = pd.DataFrame({'A': pd.Series([1, 2, np.nan], dtype='Int64'),
                       'B': pd.Series([1, np.nan, 3], dtype='UInt8'),
                       'C': [1, 2, 3]})

    result = df.A + df.C
    expected = pd.Series([2, 4, np.nan], dtype='Int64')
    tm.assert_series_equal(result, expected)

    result = (df.A + df.C) * 3 == 12
    expected = pd.Series([False, True, False])
    tm.assert_series_equal(result, expected)

    result = df.A + df.B
    expected = pd.Series([2, np.nan, np.nan], dtype='Int64')
    tm.assert_series_equal(result, expected)


def test_groupby_mean_included():
    df = pd.DataFrame({
        "A": ['a', 'b', 'b'],
        "B": [1, None, 3],
        "C": integer_array([1, None, 3], dtype='Int64'),
    })

    result = df.groupby("A").sum()
    # TODO(#22346): preserve Int64 dtype
    expected = pd.DataFrame({
        "B": np.array([1.0, 3.0]),
        "C": np.array([1, 3], dtype="int64")
    }, index=pd.Index(['a', 'b'], name='A'))
    tm.assert_frame_equal(result, expected)


def test_astype_nansafe():
    # https://github.com/pandas-dev/pandas/pull/22343
    arr = integer_array([np.nan, 1, 2], dtype="Int8")

    with tm.assert_raises_regex(
            ValueError, 'cannot convert float NaN to integer'):
        arr.astype('uint32')


# TODO(jreback) - these need testing / are broken

# shift

# set_index (destroys type)
