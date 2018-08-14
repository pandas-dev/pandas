import numpy as np
import pandas as pd
import pandas.util.testing as tm
import pytest

from pandas.tests.extension import base
from pandas.api.types import (
    is_integer, is_scalar, is_float, is_float_dtype)
from pandas.core.dtypes.generic import ABCIndexClass

from pandas.core.arrays.set import (SetDtype,
                                    to_set_array, SetArray)

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
            [np.nan] +
            [set()] + [None])


@pytest.fixture
def dtype():
    return SetDtype()


@pytest.fixture
def data():
    return SetArray(make_int_sets())


@pytest.fixture
def data_missing():
    return SetArray(make_data())


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

# class BaseInteger(object):
# 
#     def assert_index_equal(self, left, right, *args, **kwargs):
# 
#         left_na = left.isna()
#         right_na = right.isna()
# 
#         tm.assert_numpy_array_equal(left_na, right_na)
#         return tm.assert_index_equal(left[~left_na],
#                                      right[~right_na],
#                                      *args, **kwargs)
# 
#     def assert_series_equal(self, left, right, *args, **kwargs):
# 
#         left_na = left.isna()
#         right_na = right.isna()
# 
#         tm.assert_series_equal(left_na, right_na)
#         return tm.assert_series_equal(left[~left_na],
#                                       right[~right_na],
#                                       *args, **kwargs)
# 
#     def assert_frame_equal(self, left, right, *args, **kwargs):
#         # TODO(EA): select_dtypes
#         tm.assert_index_equal(
#             left.columns, right.columns,
#             exact=kwargs.get('check_column_type', 'equiv'),
#             check_names=kwargs.get('check_names', True),
#             check_exact=kwargs.get('check_exact', False),
#             check_categorical=kwargs.get('check_categorical', True),
#             obj='{obj}.columns'.format(obj=kwargs.get('obj', 'DataFrame')))
# 
#         integers = (left.dtypes == 'integer').index
# 
#         for col in integers:
#             self.assert_series_equal(left[col], right[col],
#                                      *args, **kwargs)
# 
#         left = left.drop(columns=integers)
#         right = right.drop(columns=integers)
#         tm.assert_frame_equal(left, right, *args, **kwargs)


class TestDtype(base.BaseDtypeTests):

    def test_array_type_with_arg(self, data, dtype):
        assert dtype.construct_array_type() is SetArray


class TestInterface(base.BaseInterfaceTests):

    def test_no_values_attribute(self, data):
        pytest.skip("Welp")


class TestConstructors(base.BaseConstructorsTests):
    pass


class TestReshaping(base.BaseReshapingTests):
    pass


class TestGetitem(base.BaseGetitemTests):

    @pytest.mark.skip(reason="Need to think about it.")
    def test_take_non_na_fill_value(self, data_missing):
        pass

    def test_get(self, data):
        s = pd.Series(data, index=[2 * i for i in range(len(data))])
        assert np.isnan(s.get(4)) and np.isnan(s.iloc[2])
        assert s.get(2) == s.iloc[1]


class TestGetitem(base.BaseGetitemTests):
    pass


class TestMissing(base.BaseMissingTests):

    def test_fillna_limit_pad(self):
        pass

    def test_fillna_limit_backfill(self):
        pass

    def test_fillna_series_method(self):
        pass

    def test_fillna_series(self):
        # this one looks doable.
        pass


class TestMethods(base.BaseMethodsTests):
    pass


class TestCasting(base.BaseCastingTests):
    pass


class TestArithmeticOps(base.BaseArithmeticOpsTests):
    pass


class TestComparisonOps(base.BaseComparisonOpsTests):
    pass


class TestInterface(base.BaseInterfaceTests):

    def test_repr_array(self, data):
        result = repr(data)

        # not long
        assert '...' not in result

        assert 'dtype=' in result
        assert 'SetArray' in result

    def test_repr_array_long(self, data):
        # some arrays may be able to assert a ... in the repr
        with pd.option_context('display.max_seq_items', 1):
            result = repr(data)

            assert '...' in result
            assert 'length' in result


class TestGroupby(base.BaseGroupbyTests):

    @pytest.mark.xfail(reason="groupby not working", strict=True)
    def test_groupby_extension_no_sort(self, data_for_grouping):
        super(TestGroupby, self).test_groupby_extension_no_sort(
            data_for_grouping)

    @pytest.mark.parametrize('as_index', [
        pytest.param(True,
                     marks=pytest.mark.xfail(reason="groupby not working",
                                             strict=True)),
        False
    ])
    def test_groupby_extension_agg(self, as_index, data_for_grouping):
        super(TestGroupby, self).test_groupby_extension_agg(
            as_index, data_for_grouping)

