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
    return JSONArray([{'b': 1}, {}, {'a': 4}])


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
    unhashable = pytest.mark.skip(reason="Unhashable")
    unstable = pytest.mark.skipif(sys.version_info <= (3, 5),
                                  reason="Dictionary order unstable")

    @unhashable
    def test_value_counts(self, all_data, dropna):
        pass

    @unhashable
    def test_sort_values_frame(self):
        # TODO (EA.factorize): see if _values_for_factorize allows this.
        pass

    @unstable
    def test_argsort(self, data_for_sorting):
        super(TestMethods, self).test_argsort(data_for_sorting)

    @unstable
    def test_argsort_missing(self, data_missing_for_sorting):
        super(TestMethods, self).test_argsort_missing(
            data_missing_for_sorting)

    @unstable
    @pytest.mark.parametrize('ascending', [True, False])
    def test_sort_values(self, data_for_sorting, ascending):
        super(TestMethods, self).test_sort_values(
            data_for_sorting, ascending)

    @pytest.mark.parametrize('ascending', [True, False])
    def test_sort_values_missing(self, data_missing_for_sorting, ascending):
        super(TestMethods, self).test_sort_values_missing(
            data_missing_for_sorting, ascending)


class TestCasting(base.BaseCastingTests):
    pass
