import numpy as np
import pandas as pd
import pandas.util.testing as tm
import pytest

from pandas.tests.extension import base
from pandas.core.dtypes.common import is_extension_array_dtype

from pandas.core.arrays import IntegerArray, integer_array
from pandas.core.arrays.integer import (
    Int8Dtype, Int16Dtype, Int32Dtype, Int64Dtype,
    UInt8Dtype, UInt16Dtype, UInt32Dtype, UInt64Dtype)


def make_data():
    return (list(range(1, 9)) +
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

    def check_opname(self, s, op_name, other, exc=None):
        # overwriting to indicate ops don't raise an error
        super(TestArithmeticOps, self).check_opname(s, op_name,
                                                    other, exc=None)

    def _check_op(self, s, op, other, op_name, exc=NotImplementedError):
        if exc is None:
            if s.dtype.is_unsigned_integer and (op_name == '__rsub__'):
                # TODO see https://github.com/pandas-dev/pandas/issues/22023
                pytest.skip("unsigned subtraction gives negative values")

            if (hasattr(other, 'dtype')
                    and not is_extension_array_dtype(other.dtype)
                    and pd.api.types.is_integer_dtype(other.dtype)):
                # other is np.int64 and would therefore always result in
                # upcasting, so keeping other as same numpy_dtype
                other = other.astype(s.dtype.numpy_dtype)

            result = op(s, other)
            expected = s.combine(other, op)

            if op_name == '__rdiv__':
                # combine is not giving the correct result for this case
                pytest.skip("skipping reverse div in python 2")
            elif op_name in ('__rtruediv__', '__truediv__', '__div__'):
                expected = expected.astype(float)
                if op_name == '__rtruediv__':
                    # TODO reverse operators result in object dtype
                    result = result.astype(float)
            elif op_name.startswith('__r'):
                # TODO reverse operators result in object dtype
                # see https://github.com/pandas-dev/pandas/issues/22024
                expected = expected.astype(s.dtype)
                result = result.astype(s.dtype)
            else:
                # combine method result in 'biggest' (int64) dtype
                expected = expected.astype(s.dtype)
                pass
            if (op_name == '__rpow__') and isinstance(other, pd.Series):
                # TODO pow on Int arrays gives different result with NA
                # see https://github.com/pandas-dev/pandas/issues/22022
                result = result.fillna(1)

            self.assert_series_equal(result, expected)
        else:
            with pytest.raises(exc):
                op(s, other)

    def _check_divmod_op(self, s, op, other, exc=None):
        super(TestArithmeticOps, self)._check_divmod_op(s, op, other, None)

    @pytest.mark.skip(reason="intNA does not error on ops")
    def test_error(self, data, all_arithmetic_operators):
        # other specific errors tested in the integer array specific tests
        pass


class TestComparisonOps(BaseInteger, base.BaseComparisonOpsTests):

    def check_opname(self, s, op_name, other, exc=None):
        super(TestComparisonOps, self).check_opname(s, op_name,
                                                    other, exc=None)

    def _compare_other(self, s, data, op_name, other):
        self.check_opname(s, op_name, other)


class TestInterface(BaseInteger, base.BaseInterfaceTests):
    pass


class TestConstructors(BaseInteger, base.BaseConstructorsTests):
    pass


class TestReshaping(BaseInteger, base.BaseReshapingTests):
    pass

    # for test_concat_mixed_dtypes test
    # concat of an Integer and Int coerces to object dtype
    # TODO(jreback) once integrated this would


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


class TestCasting(BaseInteger, base.BaseCastingTests):
    pass


class TestGroupby(BaseInteger, base.BaseGroupbyTests):

    @pytest.mark.xfail(reason="groupby not working", strict=True)
    def test_groupby_extension_no_sort(self, data_for_grouping):
        super(TestGroupby, self).test_groupby_extension_no_sort(
            data_for_grouping)

    @pytest.mark.parametrize('as_index', [
        pytest.param(True,
                     marks=pytest.mark.xfail(reason="groupby not working",
                                             strict=True)),
        False
    ])
    def test_groupby_extension_agg(self, as_index, data_for_grouping):
        super(TestGroupby, self).test_groupby_extension_agg(
            as_index, data_for_grouping)
