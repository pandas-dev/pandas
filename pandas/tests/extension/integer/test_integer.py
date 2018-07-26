import numpy as np
import pandas as pd
import pandas.util.testing as tm
import pytest

from pandas.tests.extension import base
from pandas.api.types import (
    is_integer, is_scalar, is_float, is_float_dtype)
from pandas.core.dtypes.generic import ABCIndexClass

from pandas.core.arrays import (
    integer_array, IntegerArray)
from pandas.core.arrays.integer import (
    Int8Dtype, Int16Dtype, Int32Dtype, Int64Dtype,
    UInt8Dtype, UInt16Dtype, UInt32Dtype, UInt64Dtype)


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


@pytest.fixture
def data_repeated(data):
    def gen(count):
        for _ in range(count):
            yield data
    yield gen


@pytest.fixture
def data_for_sorting(dtype):
    return integer_array([1, 2, 0], dtype=dtype)


@pytest.fixture
def data_missing_for_sorting(dtype):
    return integer_array([1, np.nan, 0], dtype=dtype)


@pytest.fixture
def na_cmp():
    # we are np.nan
    return lambda x, y: np.isnan(x) and np.isnan(y)


@pytest.fixture
def na_value():
    return np.nan


@pytest.fixture
def data_for_grouping(dtype):
    b = 1
    a = 0
    c = 2
    na = np.nan
    return integer_array([b, b, na, na, a, a, b, c], dtype=dtype)


def test_dtypes(dtype):
    # smoke tests on auto dtype construction

    if dtype.is_signed_integer:
        assert np.dtype(dtype.type).kind == 'i'
    else:
        assert np.dtype(dtype.type).kind == 'u'
    assert dtype.name is not None


class BaseInteger(object):

    def assert_index_equal(self, left, right, *args, **kwargs):

        left_na = left.isna()
        right_na = right.isna()

        tm.assert_numpy_array_equal(left_na, right_na)
        return tm.assert_index_equal(left[~left_na],
                                     right[~right_na],
                                     *args, **kwargs)

    def assert_series_equal(self, left, right, *args, **kwargs):

        left_na = left.isna()
        right_na = right.isna()

        tm.assert_series_equal(left_na, right_na)
        return tm.assert_series_equal(left[~left_na],
                                      right[~right_na],
                                      *args, **kwargs)

    def assert_frame_equal(self, left, right, *args, **kwargs):
        # TODO(EA): select_dtypes
        tm.assert_index_equal(
            left.columns, right.columns,
            exact=kwargs.get('check_column_type', 'equiv'),
            check_names=kwargs.get('check_names', True),
            check_exact=kwargs.get('check_exact', False),
            check_categorical=kwargs.get('check_categorical', True),
            obj='{obj}.columns'.format(obj=kwargs.get('obj', 'DataFrame')))

        integers = (left.dtypes == 'integer').index

        for col in integers:
            self.assert_series_equal(left[col], right[col],
                                     *args, **kwargs)

        left = left.drop(columns=integers)
        right = right.drop(columns=integers)
        tm.assert_frame_equal(left, right, *args, **kwargs)


class TestDtype(BaseInteger, base.BaseDtypeTests):

    @pytest.mark.skip(reason="using multiple dtypes")
    def test_is_dtype_unboxes_dtype(self):
        # we have multiple dtypes, so skip
        pass

    def test_array_type_with_arg(self, data, dtype):
        assert dtype.construct_array_type() is IntegerArray


class TestArithmeticOps(BaseInteger, base.BaseArithmeticOpsTests):

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
        self.assert_series_equal(result, expected)

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
            self.assert_series_equal(original, expected.astype('float'))

        # assert our expected result
        self.assert_series_equal(result, expected)

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


class TestComparisonOps(BaseInteger, base.BaseComparisonOpsTests):

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


class TestInterface(BaseInteger, base.BaseInterfaceTests):

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


class TestConstructors(BaseInteger, base.BaseConstructorsTests):

    def test_from_dtype_from_float(self, data):
        # construct from our dtype & string dtype
        dtype = data.dtype

        # from float
        expected = pd.Series(data)
        result = pd.Series(np.array(data).astype('float'), dtype=str(dtype))
        self.assert_series_equal(result, expected)

        # from int / list
        expected = pd.Series(data)
        result = pd.Series(np.array(data).tolist(), dtype=str(dtype))
        self.assert_series_equal(result, expected)

        # from int / array
        expected = pd.Series(data).dropna().reset_index(drop=True)
        dropped = np.array(data.dropna()).astype(np.dtype((dtype.type)))
        result = pd.Series(dropped, dtype=str(dtype))
        self.assert_series_equal(result, expected)


class TestReshaping(BaseInteger, base.BaseReshapingTests):

    def test_concat_mixed_dtypes(self, data):
        # https://github.com/pandas-dev/pandas/issues/20762
        df1 = pd.DataFrame({'A': data[:3]})
        df2 = pd.DataFrame({"A": [1, 2, 3]})
        df3 = pd.DataFrame({"A": ['a', 'b', 'c']}).astype('category')
        df4 = pd.DataFrame({"A": pd.SparseArray([1, 2, 3])})
        dfs = [df1, df2, df3, df4]

        # dataframes
        result = pd.concat(dfs)
        expected = pd.concat([x.astype(object) for x in dfs])
        self.assert_frame_equal(result, expected)

        # series
        result = pd.concat([x['A'] for x in dfs])
        expected = pd.concat([x['A'].astype(object) for x in dfs])
        self.assert_series_equal(result, expected)

        result = pd.concat([df1, df2])
        expected = pd.concat([df1.astype('object'), df2.astype('object')])
        self.assert_frame_equal(result, expected)

        # concat of an Integer and Int coerces to object dtype
        # TODO(jreback) once integrated this would
        # be a result of Integer
        result = pd.concat([df1['A'], df2['A']])
        expected = pd.concat([df1['A'].astype('object'),
                              df2['A'].astype('object')])
        self.assert_series_equal(result, expected)


class TestGetitem(BaseInteger, base.BaseGetitemTests):
    pass


class TestMissing(BaseInteger, base.BaseMissingTests):
    pass


class TestMethods(BaseInteger, base.BaseMethodsTests):

    @pytest.mark.parametrize('dropna', [True, False])
    def test_value_counts(self, all_data, dropna):
        all_data = all_data[:10]
        if dropna:
            other = np.array(all_data[~all_data.isna()])
        else:
            other = all_data

        result = pd.Series(all_data).value_counts(dropna=dropna).sort_index()
        expected = pd.Series(other).value_counts(
            dropna=dropna).sort_index()
        expected.index = expected.index.astype(all_data.dtype)

        self.assert_series_equal(result, expected)

    def test_combine_add(self, data_repeated):
        # GH 20825
        orig_data1, orig_data2 = data_repeated(2)
        s1 = pd.Series(orig_data1)
        s2 = pd.Series(orig_data2)

        # fundamentally this is not a great operation
        # as overflow / underflow can easily happen here
        # e.g. int8 + int8
        def scalar_add(a, b):

            # TODO; should really be a type specific NA
            if pd.isna(a) or pd.isna(b):
                return np.nan
            if is_integer(a):
                a = int(a)
            elif is_integer(b):
                b = int(b)
            return a + b

        result = s1.combine(s2, scalar_add)
        expected = pd.Series(
            orig_data1._from_sequence([scalar_add(a, b) for (a, b) in
                                       zip(orig_data1,
                                           orig_data2)]))
        self.assert_series_equal(result, expected)

        val = s1.iloc[0]
        result = s1.combine(val, lambda x1, x2: x1 + x2)
        expected = pd.Series(
            orig_data1._from_sequence([a + val for a in list(orig_data1)]))
        self.assert_series_equal(result, expected)


class TestCasting(BaseInteger, base.BaseCastingTests):

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

        self.assert_index_equal(result, expected)

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
        self.assert_index_equal(result, expected)

    def test_astype(self, all_data):
        all_data = all_data[:10]

        ints = all_data[~all_data.isna()]
        mixed = all_data
        dtype = Int8Dtype()

        # coerce to same type - ints
        s = pd.Series(ints)
        result = s.astype(all_data.dtype)
        expected = pd.Series(ints)
        self.assert_series_equal(result, expected)

        # coerce to same other - ints
        s = pd.Series(ints)
        result = s.astype(dtype)
        expected = pd.Series(ints, dtype=dtype)
        self.assert_series_equal(result, expected)

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
        self.assert_series_equal(result, expected)

        # coerce to same other - mixed
        s = pd.Series(mixed)
        result = s.astype(dtype)
        expected = pd.Series(mixed, dtype=dtype)
        self.assert_series_equal(result, expected)

        # coerce to same numpy_dtype - mixed
        s = pd.Series(mixed)
        with pytest.raises(ValueError):
            s.astype(all_data.dtype.numpy_dtype)

        # coerce to object
        s = pd.Series(mixed)
        result = s.astype('object')
        expected = pd.Series(np.asarray(mixed))
        tm.assert_series_equal(result, expected)

    @pytest.mark.parametrize('dtype', [Int8Dtype(), 'Int8'])
    def test_astype_specific_casting(self, dtype):
        s = pd.Series([1, 2, 3], dtype='Int64')
        result = s.astype(dtype)
        expected = pd.Series([1, 2, 3], dtype='Int8')
        self.assert_series_equal(result, expected)

        s = pd.Series([1, 2, 3, None], dtype='Int64')
        result = s.astype(dtype)
        expected = pd.Series([1, 2, 3, None], dtype='Int8')
        self.assert_series_equal(result, expected)

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


class TestGroupby(BaseInteger, base.BaseGroupbyTests):

    @pytest.mark.xfail(reason="groupby not working")
    def test_groupby_extension_no_sort(self, data_for_grouping):
        super(TestGroupby, self).test_groupby_extension_no_sort(
            data_for_grouping)

    @pytest.mark.xfail(reason="groupby not working")
    @pytest.mark.parametrize('as_index', [True, False])
    def test_groupby_extension_agg(self, as_index, data_for_grouping):
        super(TestGroupby, self).test_groupby_extension_agg(
            as_index, data_for_grouping)


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


@pytest.mark.parametrize(
    'values',
    [
        ['foo', 'bar'],
        'foo',
        1,
        1.0,
        pd.date_range('20130101', periods=2),
        np.array(['foo'])])
def test_to_integer_array_error(values):
    # error in converting existing arrays to IntegerArrays
    with pytest.raises(TypeError):
        integer_array(values)


@pytest.mark.parametrize(
    'values, to_dtype, result_dtype',
    [
        (np.array([1], dtype='int64'), None, Int64Dtype),
        (np.array([1, np.nan]), None, Int64Dtype),
        (np.array([1, np.nan]), 'int8', Int8Dtype)])
def test_to_integer_array(values, to_dtype, result_dtype):
    # convert existing arrays to IntegerArrays
    result = integer_array(values, dtype=to_dtype)
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


# TODO(jreback) - these need testing / are broken

# shift

# set_index (destroys type)
