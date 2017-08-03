# -*- coding: utf-8 -*-

import pytest

import numpy as np
from numpy.testing import assert_allclose

import pandas as pd
from pandas.core.dtypes.costum_dtypes import NumpyDtypeWithMetadata, AlwaysSame
from pandas.core.dtypes.common import is_dtype_equal

from .test_dtypes import Base

def test_alwys_same_series_construction():
    a = pd.Series([1,2,3,4,5,6], dtype = AlwaysSame(42))
    # The '==' operator for the whole series will be tested elsewhere
    for i in range(0,6):
        assert a[i]==42

def test_always_same_with_series_eq():
    a = pd.Series([1,2,3,4,5,6], dtype = AlwaysSame(42))
    assert np.all(a == a)
    assert np.all(a == pd.Series([42,42,42,42,42,42], dtype = int))
    assert np.all(a != pd.Series([1,2,4,4,4,4], dtype = AlwaysSame(41)))
    assert np.all(a == pd.Series([1,2,4,4,4,4], dtype = AlwaysSame(42)))


def test_always_same_with_series_add():
    a = pd.Series([1,2,3,4,5,6], dtype = AlwaysSame(42))
    assert (a == a + a).all()


def test_always_same_with_series_add_mixing_types():
    a = pd.Series([1,2,3,4,5,6], dtype = AlwaysSame(42))
    b = pd.Series([1,2,3,4,5,6], dtype = AlwaysSame(43))
    c = pd.Series([1,2,3,4,5,6], dtype = int)
    pytest.raises(TypeError,
                  lambda: a + b)
    pytest.raises(TypeError,
                  lambda: b + a)
    assert (a + c == a).all()
    assert (c + a == a).all()


class TestAlwaysSame(Base):
    def create(self):
        return AlwaysSame(42)

    def test_hash_vs_equality(self):
        # make sure that we satisfy is semantics
        dtype = self.dtype
        dtype2 = AlwaysSame(42)
        assert dtype == dtype2
        assert dtype2 == dtype
        assert dtype is dtype2
        assert dtype2 is dtype
        assert hash(dtype) == hash(dtype2)


class TestNumpyDtypeWithMetadata:

    def test_is_not_dtype(self):
        assert not NumpyDtypeWithMetadata.is_dtype(None)
        assert not NumpyDtypeWithMetadata.is_dtype(np.float64)
        assert NumpyDtypeWithMetadata.is_dtype(AlwaysSame('Beer'))
