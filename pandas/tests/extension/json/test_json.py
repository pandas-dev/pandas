import operator

import pytest


from pandas.compat import PY2, PY36
from pandas.tests.extension import base

from .array import JSONArray, JSONDtype, make_data

pytestmark = pytest.mark.skipif(PY2, reason="Py2 doesn't have a UserDict")


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
    @pytest.mark.xfail(reason="Setting a dict as a scalar")
    def test_fillna_series(self):
        """We treat dictionaries as a mapping in fillna, not a scalar."""

    @pytest.mark.xfail(reason="Setting a dict as a scalar")
    def test_fillna_frame(self):
        """We treat dictionaries as a mapping in fillna, not a scalar."""


unhashable = pytest.mark.skip(reason="Unhashable")
unstable = pytest.mark.skipif(not PY36,  # 3.6 or higher
                              reason="Dictionary order unstable")


class TestMethods(base.BaseMethodsTests):
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

    @unstable
    @pytest.mark.parametrize('ascending', [True, False])
    def test_sort_values_missing(self, data_missing_for_sorting, ascending):
        super(TestMethods, self).test_sort_values_missing(
            data_missing_for_sorting, ascending)


class TestCasting(base.BaseCastingTests):
    pass


class TestGroupby(base.BaseGroupbyTests):

    @unhashable
    def test_groupby_extension_transform(self):
        """
        This currently fails in Series.name.setter, since the
        name must be hashable, but the value is a dictionary.
        I think this is what we want, i.e. `.name` should be the original
        values, and not the values for factorization.
        """

    @unhashable
    def test_groupby_extension_apply(self):
        """
        This fails in Index._do_unique_check with

        >   hash(val)
        E   TypeError: unhashable type: 'UserDict' with

        I suspect that once we support Index[ExtensionArray],
        we'll be able to dispatch unique.
        """

    @unstable
    @pytest.mark.parametrize('as_index', [True, False])
    def test_groupby_extension_agg(self, as_index, data_for_grouping):
        super(TestGroupby, self).test_groupby_extension_agg(
            as_index, data_for_grouping
        )
