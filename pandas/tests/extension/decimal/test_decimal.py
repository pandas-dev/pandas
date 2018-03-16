import decimal

import numpy as np
import pandas as pd
import pandas.util.testing as tm
import pytest

from pandas.tests.extension import base

from .array import DecimalDtype, DecimalArray, make_data


@pytest.fixture
def dtype():
    return DecimalDtype()


@pytest.fixture
def data():
    return DecimalArray(make_data())


@pytest.fixture
def data_missing():
    return DecimalArray([decimal.Decimal('NaN'), decimal.Decimal(1)])


@pytest.fixture
def na_cmp():
    return lambda x, y: x.is_nan() and y.is_nan()


@pytest.fixture
def na_value():
    return decimal.Decimal("NaN")


class BaseDecimal(object):

    def assert_series_equal(self, left, right, *args, **kwargs):

        left_na = left.isna()
        right_na = right.isna()

        tm.assert_series_equal(left_na, right_na)
        return tm.assert_series_equal(left[~left_na],
                                      right[~right_na],
                                      *args, **kwargs)

    def assert_frame_equal(self, left, right, *args, **kwargs):
        self.assert_series_equal(left.dtypes, right.dtypes)
        for col in left.columns:
            self.assert_series_equal(left[col], right[col],
                                     *args, **kwargs)


class TestDtype(BaseDecimal, base.BaseDtypeTests):
    pass


class TestInterface(BaseDecimal, base.BaseInterfaceTests):
    pass


class TestConstructors(BaseDecimal, base.BaseConstructorsTests):
    pass


class TestReshaping(BaseDecimal, base.BaseReshapingTests):

    def test_align(self, data, na_value):
        # Have to override since assert_series_equal doesn't
        # compare Decimal(NaN) properly.
        a = data[:3]
        b = data[2:5]
        r1, r2 = pd.Series(a).align(pd.Series(b, index=[1, 2, 3]))

        # NaN handling
        e1 = pd.Series(type(data)(list(a) + [na_value]))
        e2 = pd.Series(type(data)([na_value] + list(b)))
        tm.assert_series_equal(r1.iloc[:3], e1.iloc[:3])
        assert r1[3].is_nan()
        assert e1[3].is_nan()

        tm.assert_series_equal(r2.iloc[1:], e2.iloc[1:])
        assert r2[0].is_nan()
        assert e2[0].is_nan()

    def test_align_frame(self, data, na_value):
        # Override for Decimal(NaN) comparison
        a = data[:3]
        b = data[2:5]
        r1, r2 = pd.DataFrame({'A': a}).align(
            pd.DataFrame({'A': b}, index=[1, 2, 3])
        )

        # Assumes that the ctor can take a list of scalars of the type
        e1 = pd.DataFrame({'A': type(data)(list(a) + [na_value])})
        e2 = pd.DataFrame({'A': type(data)([na_value] + list(b))})

        tm.assert_frame_equal(r1.iloc[:3], e1.iloc[:3])
        assert r1.loc[3, 'A'].is_nan()
        assert e1.loc[3, 'A'].is_nan()

        tm.assert_frame_equal(r2.iloc[1:], e2.iloc[1:])
        assert r2.loc[0, 'A'].is_nan()
        assert e2.loc[0, 'A'].is_nan()


class TestGetitem(BaseDecimal, base.BaseGetitemTests):
    pass


class TestMissing(BaseDecimal, base.BaseMissingTests):
    pass


class TestMethods(BaseDecimal, base.BaseMethodsTests):
    @pytest.mark.parametrize('dropna', [True, False])
    @pytest.mark.xfail(reason="value_counts not implemented yet.")
    def test_value_counts(self, all_data, dropna):
        all_data = all_data[:10]
        if dropna:
            other = np.array(all_data[~all_data.isna()])
        else:
            other = all_data

        result = pd.Series(all_data).value_counts(dropna=dropna).sort_index()
        expected = pd.Series(other).value_counts(dropna=dropna).sort_index()

        tm.assert_series_equal(result, expected)


class TestCasting(BaseDecimal, base.BaseCastingTests):
    pass


def test_series_constructor_coerce_data_to_extension_dtype_raises():
    xpr = ("Cannot cast data to extension dtype 'decimal'. Pass the "
           "extension array directly.")
    with tm.assert_raises_regex(ValueError, xpr):
        pd.Series([0, 1, 2], dtype=DecimalDtype())


def test_series_constructor_with_same_dtype_ok():
    arr = DecimalArray([decimal.Decimal('10.0')])
    result = pd.Series(arr, dtype=DecimalDtype())
    expected = pd.Series(arr)
    tm.assert_series_equal(result, expected)


def test_series_constructor_coerce_extension_array_to_dtype_raises():
    arr = DecimalArray([decimal.Decimal('10.0')])
    xpr = "Cannot specify a dtype 'int64' .* \('decimal'\)."

    with tm.assert_raises_regex(ValueError, xpr):
        pd.Series(arr, dtype='int64')


def test_dataframe_constructor_with_same_dtype_ok():
    arr = DecimalArray([decimal.Decimal('10.0')])

    result = pd.DataFrame({"A": arr}, dtype=DecimalDtype())
    expected = pd.DataFrame({"A": arr})
    tm.assert_frame_equal(result, expected)


def test_dataframe_constructor_with_different_dtype_raises():
    arr = DecimalArray([decimal.Decimal('10.0')])

    xpr = "Cannot coerce extension array to dtype 'int64'. "
    with tm.assert_raises_regex(ValueError, xpr):
        pd.DataFrame({"A": arr}, dtype='int64')
