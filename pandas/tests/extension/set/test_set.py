import numpy as np
import pandas.util.testing as tm
import pytest

from pandas.tests.extension import base

from pandas.core.arrays.set import SetDtype, SetArray


def make_string_sets():
    s = tm.makeStringSeries()
    return s.index.map(set).values


def make_int_sets():
    s = tm.makeFloatSeries().astype(str).str.replace(r'\D', '')
    return s.map(lambda x: set(map(int, x))).values


def make_data():
    return (list(make_string_sets()) +
            [np.nan] +
            list(make_int_sets()) +
            [np.nan, None, set()])


@pytest.fixture
def dtype():
    return SetDtype()


@pytest.fixture
def data():
    return SetArray(make_int_sets())


@pytest.fixture
def data_missing():
    return SetArray([np.nan, {1}])


@pytest.fixture
def data_repeated(data):
    def gen(count):
        for _ in range(count):
            yield data
    yield gen


# @pytest.fixture
# def data_for_sorting(dtype):
#     return SetArray(...)


# @pytest.fixture
# def data_missing_for_sorting(dtype):
#     return SetArray(...)


@pytest.fixture
def na_cmp():
    # we are np.nan
    return lambda x, y: np.isnan(x) and np.isnan(y)


@pytest.fixture
def na_value():
    return np.nan

# @pytest.fixture
# def data_for_grouping(dtype):
#     return SetArray(...)


class TestDtype(base.BaseDtypeTests):

    def test_array_type_with_arg(self, data, dtype):
        assert dtype.construct_array_type() is SetArray


class TestInterface(base.BaseInterfaceTests):

    def test_len(self, data):
        assert len(data) == 30


class TestConstructors(base.BaseConstructorsTests):
    pass


class TestReshaping(base.BaseReshapingTests):
    pass


class TestGetitem(base.BaseGetitemTests):

    @pytest.mark.skip(reason="Need to think about it.")
    def test_take_non_na_fill_value(self, data_missing):
        pass


class TestSetitem(base.BaseGetitemTests):
    pass


class TestMissing(base.BaseMissingTests):

    def test_fillna_frame(self, data_missing):
        pytest.skip('df.fillna does not dispatch to EA')

    def test_fillna_limit_pad(self):
        pytest.skip('TODO')

    def test_fillna_limit_backfill(self):
        pytest.skip('TODO')

    def test_fillna_series_method(self):
        pytest.skip('TODO')

    def test_fillna_series(self):
        pytest.skip('series.fillna does not dispatch to EA')


# # most methods (value_counts, unique, factorize) will not be for SetArray
# # rest still buggy
class TestMethods(base.BaseMethodsTests):
    pass


class TestCasting(base.BaseCastingTests):
    pass


class TestArithmeticOps(base.BaseArithmeticOpsTests):

    def check_opname(self, s, op_name, other, exc='ignored'):
        op = self.get_op_from_name(op_name)

        self._check_op(s, op, other,
                       None if op_name == '__sub__' else NotImplementedError)

    def test_divmod(self, data):
        pytest.skip('Not relevant')

    def test_error(self, data, all_arithmetic_operators):
        pytest.skip('TODO')


class TestComparisonOps(base.BaseComparisonOpsTests):

    def _compare_other(self, s, data, op_name, other):
        op = self.get_op_from_name(op_name)
        result = op(s, other)
        expected = s.combine(other, op)
        self.assert_series_equal(result, expected)

# # GroupBy won't be implemented for SetArray
# class TestGroupby(base.BaseGroupbyTests):
#     pass
