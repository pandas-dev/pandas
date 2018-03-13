import operator
import sys

import pytest


from pandas.tests.extension import base

from .array import JSONArray, JSONDtype, make_data

pytestmark = pytest.mark.skipif(sys.version_info[0] == 2,
                                reason="Py2 doesn't have a UserDict")


@pytest.fixture
def dtype():
    return JSONDtype()


@pytest.fixture
def data():
    """Length-100 PeriodArray for semantics test."""
    return JSONArray(make_data())


@pytest.fixture
def data_missing():
    """Length 2 array with [NA, Valid]"""
    return JSONArray([{}, {'a': 10}])


@pytest.fixture
def na_value():
    return {}


@pytest.fixture
def na_cmp():
    return operator.eq


@pytest.fixture
def data_for_grouping():
    return JSONArray([
        {'b': 1}, {'b': 1},
        {}, {},
        {'a': 0, 'c': 2}, {'a': 0, 'c': 2},
        {'b': 1},
        {'c': 2},
    ])


class TestDtype(base.BaseDtypeTests):
    pass


class TestInterface(base.BaseInterfaceTests):
    pass


class TestConstructors(base.BaseConstructorsTests):
    pass


class TestReshaping(base.BaseReshapingTests):
    pass


class TestGetitem(base.BaseGetitemTests):
    pass


class TestMissing(base.BaseMissingTests):
    pass


class TestMethods(base.BaseMethodsTests):
    unhashable = pytest.mark.skip(reason="Unhashable")

    @unhashable
    def test_factorize(self):
        pass


class TestCasting(base.BaseCastingTests):
    pass
