import decimal

import numpy as np
import pandas as pd
import pandas.util.testing as tm
import pytest

from pandas.tests.extension import base

from .array import DecimalDtype, DecimalNoNaArray, make_data
from .test_decimal import (
    BaseDecimal, TestDtype, TestInterface, TestConstructors,
    TestReshaping, TestGetitem, TestMissing, TestCasting,
    TestGroupby)


@pytest.fixture
def dtype():
    return DecimalDtype()


@pytest.fixture
def data():
    return DecimalNoNaArray(make_data())


@pytest.fixture
def data_missing():
    pytest.skip("No missing data tests for _can_hold_na=False")


@pytest.fixture
def data_for_sorting():
    return DecimalNoNaArray([decimal.Decimal('1'),
                             decimal.Decimal('2'),
                             decimal.Decimal('0')])


@pytest.fixture
def data_missing_for_sorting():
    pytest.skip("No missing data tests for _can_hold_na=False")


@pytest.fixture
def na_cmp():
    pytest.skip("No missing data tests for _can_hold_na=False")


@pytest.fixture
def na_value():
    pytest.skip("No missing data tests for _can_hold_na=False")


@pytest.fixture
def data_for_grouping():
    b = decimal.Decimal('1.0')
    a = decimal.Decimal('0.0')
    c = decimal.Decimal('2.0')
    d = decimal.Decimal('-1.0')
    return DecimalNoNaArray([b, b, d, d, a, a, b, c])


class TestNoNaDtype(TestDtype):
    pass


class TestNoNaInterface(TestInterface):
    pass


class TestNoNaConstructors(TestConstructors):
    pass


class TestNoNaReshaping(TestReshaping):
    pass


class TestNoNaGetitem(TestGetitem):
    pass


class TestNoNaMissing(TestMissing):
    pass


class TestNoNaMethods(BaseDecimal, base.BaseMethodsTests):
    def test_factorize(self, data_for_grouping, na_sentinel=None):
        labels, uniques = pd.factorize(data_for_grouping)
        expected_labels = np.array([0, 0, 1,
                                   1, 2, 2, 0, 3],
                                   dtype=np.intp)
        expected_uniques = data_for_grouping.take([0, 2, 4, 7])

        tm.assert_numpy_array_equal(labels, expected_labels)
        self.assert_extension_array_equal(uniques, expected_uniques)


class TestNoNaCasting(TestCasting):
    pass


class TestNoNaGroupby(TestGroupby):
    pass


def test_series_constructor_with_same_dtype_ok():
    arr = DecimalNoNaArray([decimal.Decimal('10.0')])
    result = pd.Series(arr, dtype=DecimalDtype())
    expected = pd.Series(arr)
    tm.assert_series_equal(result, expected)


def test_series_constructor_coerce_extension_array_to_dtype_raises():
    arr = DecimalNoNaArray([decimal.Decimal('10.0')])
    xpr = r"Cannot specify a dtype 'int64' .* \('decimal'\)."

    with tm.assert_raises_regex(ValueError, xpr):
        pd.Series(arr, dtype='int64')


def test_dataframe_constructor_with_same_dtype_ok():
    arr = DecimalNoNaArray([decimal.Decimal('10.0')])

    result = pd.DataFrame({"A": arr}, dtype=DecimalDtype())
    expected = pd.DataFrame({"A": arr})
    tm.assert_frame_equal(result, expected)


def test_dataframe_constructor_with_different_dtype_raises():
    arr = DecimalNoNaArray([decimal.Decimal('10.0')])

    xpr = "Cannot coerce extension array to dtype 'int64'. "
    with tm.assert_raises_regex(ValueError, xpr):
        pd.DataFrame({"A": arr}, dtype='int64')
