import string

import pytest
import pandas as pd
import numpy as np

from pandas.core.sparse.dtype import SparseDtype
from pandas import SparseArray
from pandas.tests.extension import base


def make_data():
    data = np.random.uniform(size=100)
    data[1::3] = np.nan
    return data


@pytest.fixture
def dtype():
    return SparseDtype()


@pytest.fixture
def data():
    """Length-100 PeriodArray for semantics test."""
    res = SparseArray(make_data())
    return res


@pytest.fixture
def data_missing():
    """Length 2 array with [NA, Valid]"""
    return SparseArray([np.nan, 1.0])


@pytest.fixture
def data_repeated():
    """Return different versions of data for count times"""
    def gen(count):
        for _ in range(count):
            yield SparseArray(make_data())
    yield gen


@pytest.fixture
def data_for_sorting():
    return SparseArray([1, 2, 3])


@pytest.fixture
def data_missing_for_sorting():
    return SparseArray([1, np.nan, 2])


@pytest.fixture
def na_value():
    return np.nan


@pytest.fixture
def data_for_grouping():
    return SparseArray([1, 1, np.nan, np.nan, 2, 2, 1, 3])


class TestDtype(base.BaseDtypeTests):

    def test_array_type_with_arg(self, data, dtype):
        assert dtype.construct_array_type() is SparseArray


class TestInterface(base.BaseInterfaceTests):
    def test_no_values_attribute(self, data):
        pytest.skip("Welp")


class TestConstructors(base.BaseConstructorsTests):
    pass


class TestReshaping(base.BaseReshapingTests):
    pass


class TestGetitem(base.BaseGetitemTests):
    pass


class TestSetitem(base.BaseSetitemTests):
    pass


class TestMissing(base.BaseMissingTests):
    pass


class TestMethods(base.BaseMethodsTests):
    pass


class TestCasting(base.BaseCastingTests):
    pass


class TestArithmeticOps(base.BaseArithmeticOpsTests):
    pass


class TestComparisonOps(base.BaseComparisonOpsTests):
    pass
