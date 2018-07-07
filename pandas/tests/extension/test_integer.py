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
import numpy as np
import pandas as pd
import pytest

from pandas.tests.extension import base
from pandas.core.dtypes.common import is_extension_array_dtype

from pandas.core.arrays import integer_array
from pandas.core.arrays.integer import (
    Int8Dtype, Int16Dtype, Int32Dtype, Int64Dtype,
    UInt8Dtype, UInt16Dtype, UInt32Dtype, UInt64Dtype)


def make_data():
    return (list(range(1, 9)) + [np.nan] + list(range(10, 98))
            + [np.nan] + [99, 100])


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


class TestDtype(base.BaseDtypeTests):

    @pytest.mark.skip(reason="using multiple dtypes")
    def test_is_dtype_unboxes_dtype(self):
        # we have multiple dtypes, so skip
        pass


class TestArithmeticOps(base.BaseArithmeticOpsTests):

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

    @pytest.mark.xfail(reason="EA is listified. GH-22922", strict=True)
    def test_add_series_with_extension_array(self, data):
        super(TestArithmeticOps, self).test_add_series_with_extension_array(
            data
        )


class TestComparisonOps(base.BaseComparisonOpsTests):

    def check_opname(self, s, op_name, other, exc=None):
        super(TestComparisonOps, self).check_opname(s, op_name,
                                                    other, exc=None)

    def _compare_other(self, s, data, op_name, other):
        self.check_opname(s, op_name, other)


class TestInterface(base.BaseInterfaceTests):
    pass


class TestConstructors(base.BaseConstructorsTests):
    pass


class TestReshaping(base.BaseReshapingTests):
    pass

    # for test_concat_mixed_dtypes test
    # concat of an Integer and Int coerces to object dtype
    # TODO(jreback) once integrated this would


class TestGetitem(base.BaseGetitemTests):
    pass


class TestMissing(base.BaseMissingTests):
    pass


class TestMethods(base.BaseMethodsTests):

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


class TestCasting(base.BaseCastingTests):
    pass


class TestGroupby(base.BaseGroupbyTests):

    @pytest.mark.parametrize('as_index', [True, False])
    def test_groupby_extension_agg(self, as_index, data_for_grouping):
        df = pd.DataFrame({"A": [1, 1, 2, 2, 3, 3, 1, 4],
                           "B": data_for_grouping})
        result = df.groupby("B", as_index=as_index).A.mean()
        _, index = pd.factorize(data_for_grouping, sort=True)

        # TODO(ExtensionIndex): remove coercion to object
        # we don't have an easy way to represent an EA as an Index object
        index = pd.Index(index, name="B", dtype=object)
        expected = pd.Series([3, 1, 4], index=index, name="A")
        if as_index:
            self.assert_series_equal(result, expected)
        else:
            expected = expected.reset_index()
            self.assert_frame_equal(result, expected)

    def test_groupby_extension_no_sort(self, data_for_grouping):
        df = pd.DataFrame({"A": [1, 1, 2, 2, 3, 3, 1, 4],
                           "B": data_for_grouping})
        result = df.groupby("B", sort=False).A.mean()
        _, index = pd.factorize(data_for_grouping, sort=False)

        # TODO(ExtensionIndex): remove coercion to object
        # we don't have an easy way to represent an EA as an Index object
        index = pd.Index(index, name="B", dtype=object)
        expected = pd.Series([1, 3, 4], index=index, name="A")
        self.assert_series_equal(result, expected)


class TestNumericReduce(base.BaseNumericReduceTests):
    pass


class TestBooleanReduce(base.BaseBooleanReduceTests):
    pass
