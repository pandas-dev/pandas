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
def data_for_sorting():
    return JSONArray([{'b': 1}, {'c': 4}, {'a': 2, 'c': 3}])


@pytest.fixture
def data_missing_for_sorting():
    return JSONArray([{'b': 1}, {}, {'c': 4}])


@pytest.fixture
def na_value():
    return {}


@pytest.fixture
def na_cmp():
    return operator.eq


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
    @pytest.mark.xfail(reason="Setting a dict as a scalar")
    def test_fillna_series(self):
        """We treat dictionaries as a mapping in fillna, not a scalar."""

    @pytest.mark.xfail(reason="Setting a dict as a scalar")
    def test_fillna_frame(self):
        """We treat dictionaries as a mapping in fillna, not a scalar."""


class TestMethods(base.BaseMethodsTests):
    @pytest.mark.skip(reason="Unhashable")
    def test_value_counts(self, all_data, dropna):
        pass

    @pytest.mark.skip(reason="Dictionaries are not orderable.")
    def test_argsort(self):
        pass

    @pytest.mark.skip(reason="Dictionaries are not orderable.")
    def test_argsort_missing(self):
        pass

    @pytest.mark.skip(reason="Dictionaries are not orderable.")
    def test_sort_values(self):
        pass

    @pytest.mark.skip(reason="Dictionaries are not orderable.")
    def test_sort_values_missing(self):
        pass

    @pytest.mark.skip(reason="Dictionaries are not orderable.")
    def test_sort_values_frame(self):
        pass


class TestCasting(base.BaseCastingTests):
    pass
