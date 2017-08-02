# -*- coding: utf-8 -*-

import pytest

import numpy as np
from numpy.testing import assert_allclose

import pandas as pd
from pandas.core.dtypes.units import DimensionedFloatDtype
from pandas.core.dtypes.common import is_dtype_equal

from .test_dtypes import Base


def test_construction_string():
    """
    Assert that ainstances of UnitDType have
    """
    a = DimensionedFloatDtype("meter")
    assert isinstance(a, DimensionedFloatDtype)


def test_equality():
    assert DimensionedFloatDtype("meter") == DimensionedFloatDtype("meter")
    a = DimensionedFloatDtype("meter")
    assert a == a


def test_series_construction():
    a = pd.Series([1, 2, 3, 4, 5], dtype=DimensionedFloatDtype("meter"))
    assert a.dtype == DimensionedFloatDtype("meter")

def test_series_addition_same_unit():
    a = pd.Series([1, 2, 3, 4, 5], dtype=DimensionedFloatDtype("meter"))
    b = a + a
    assert b[0] == 2
    assert b.dtype == DimensionedFloatDtype("meter")
    assert_allclose(b, pd.Series([2, 4, 6, 8, 10.],
              dtype=DimensionedFloatDtype("meter")))

def test_series_addition_compatible_units():
    a = pd.Series([1, 2, 3, 4, 5], dtype=DimensionedFloatDtype("meter"))
    c = pd.Series([5, 10, 50, 100, 500],
                  dtype=DimensionedFloatDtype("centimeter"))
    assert_allclose(a + c,
                    pd.Series([1.05, 2.1, 3.5, 5, 10.],
                              dtype=DimensionedFloatDtype("meter")))
                              
def test_series_addition_compatible_units():
    a = pd.Series([1, 2, 3, 4, 5], dtype=DimensionedFloatDtype("meter"))
    c = pd.Series([5, 10, 50, 100, 500],
                  dtype=DimensionedFloatDtype("centimeter"))
    product = a * c
    assert_allclose(product,
                    pd.Series([5, 20, 150, 400, 2500.],
                              dtype=DimensionedFloatDtype("meter*centimeter")))
    assert product.dtype == DimensionedFloatDtype("meter*centimeter")

class TestUnitDtype(Base):
    def create(self):
        return DimensionedFloatDtype("meter")

    def test_hash_vs_equality(self):
        # make sure that we satisfy is semantics
        dtype = self.dtype
        dtype2 = DimensionedFloatDtype('meter')
        assert dtype == dtype2
        assert dtype2 == dtype
        assert dtype is dtype2
        assert dtype2 is dtype
        assert hash(dtype) == hash(dtype2)

    def test_construction(self):
        pytest.raises(Exception,  # pint.UndefinedUnitError
                      lambda: DimensionedFloatDtype('thisIsNotAUnit'))

    def test_construction_from_string(self):
        result = DimensionedFloatDtype.construct_from_string(
            'dimensionedFloat[meter]')
        assert is_dtype_equal(self.dtype, result)
        pytest.raises(TypeError,
                    lambda: DimensionedFloatDtype.construct_from_string('foo'))

    def test_is_dtype(self):
        assert not DimensionedFloatDtype.is_dtype(None)
        assert DimensionedFloatDtype.is_dtype(self.dtype)
        assert DimensionedFloatDtype.is_dtype('dimensionedFloat[meter]')
        assert not DimensionedFloatDtype.is_dtype('foo')
        assert DimensionedFloatDtype.is_dtype(DimensionedFloatDtype('hours'))
        assert not DimensionedFloatDtype.is_dtype(np.float64)

    def test_equality(self):
        assert is_dtype_equal(self.dtype, 'dimensionedFloat[meter]')
        assert not is_dtype_equal(self.dtype, 'foo')
