import operator

import numpy as np
import pandas as pd
import pandas.util.testing as tm
import pytest

from pandas.tests.extension import base

from .array import RelObjectDtype, RelObjectArray, RelObj, make_data


@pytest.fixture
def dtype():
    return RelObjectDtype()


@pytest.fixture
def data():
    return RelObjectArray(make_data())


@pytest.fixture
def data_missing():
    # Since _can_hold_na is False, we don't have missing values
    return RelObjectArray([RelObj(10), RelObj(5)])


@pytest.fixture
def data_for_sorting():
    return RelObjectArray([RelObj(4), RelObj(5), RelObj(3)])


@pytest.fixture
def data_missing_for_sorting():
    # Since _can_hold_na is False, we don't have missing values
    # Tests assume middle value is smallest
    return RelObjectArray([RelObj(1), RelObj(-1), RelObj(0)])


@pytest.fixture
def na_cmp():
    return lambda x, y: x is np.nan and y is np.nan


@pytest.fixture
def na_value():
    return np.nan


@pytest.fixture
def data_for_grouping():
    b = RelObj(1)
    a = RelObj(0)
    c = RelObj(2)
    d = RelObj(-1)
    return RelObjectArray([b, b, d, d, a, a, b, c])


class BaseRelObject(object):

    def assert_series_equal(self, left, right, *args, **kwargs):

        result = tm.assert_series_equal(left, right,
                                        *args, **kwargs)
        if result:
            diff = 0
            for l, r in zip(left.values, right.values):
                if l is not r:
                    diff += 1
            if diff > 0:
                diff = diff * 100.0 / len(left)
                obj = 'RelObjArray'
                msg = '{obj} values are different ({pct} %)'.format(
                    obj='RelObjArray', pct=np.round(diff, 5))
                tm.raise_assert_detail(obj, msg, left, right)
        return result

    def assert_frame_equal(self, left, right, *args, **kwargs):
        relobjs = (left.dtypes == 'relobj').index

        for col in relobjs:
            self.assert_series_equal(left[col], right[col],
                                     *args, **kwargs)

        left = left.drop(columns=relobjs)
        right = right.drop(columns=relobjs)
        tm.assert_frame_equal(left, right, *args, **kwargs)


class TestDtype(BaseRelObject, base.BaseDtypeTests):
    pass


class TestInterface(BaseRelObject, base.BaseInterfaceTests):
    pass


class TestConstructors(BaseRelObject, base.BaseConstructorsTests):
    pass


class TestReshaping(BaseRelObject, base.BaseReshapingTests):
    pass


class TestGetitem(BaseRelObject, base.BaseGetitemTests):
    pass


class TestMissing(BaseRelObject, base.BaseMissingTests):
    pass


class TestMethods(BaseRelObject, base.BaseMethodsTests):
    @pytest.mark.parametrize('dropna', [True, False])
    @pytest.mark.xfail(reason="value_counts not implemented yet.")
    def test_value_counts(self, all_data, dropna):
        all_data = all_data[:10]
        other = all_data

        result = pd.Series(all_data).value_counts(dropna=dropna).sort_index()
        expected = pd.Series(other).value_counts(dropna=dropna).sort_index()

        tm.assert_series_equal(result, expected)

    @pytest.mark.xfail(reason="sorting not appropriate")
    def test_argsort(self, data_for_sorting):
        pass

    @pytest.mark.xfail(reason="sorting not appropriate")
    def test_argsort_missing(self, data_missing_for_sorting):
        pass

    @pytest.mark.parametrize('ascending', [True, False])
    @pytest.mark.xfail(reason="sorting not appropriate")
    def test_sort_values(self, data_for_sorting, ascending):
        pass

    @pytest.mark.parametrize('ascending', [True, False])
    @pytest.mark.xfail(reason="sorting not appropriate")
    def test_sort_values_missing(self, data_missing_for_sorting, ascending):
        pass

    @pytest.mark.parametrize('ascending', [True, False])
    @pytest.mark.xfail(reason="sorting not appropriate")
    def test_sort_values_frame(self, data_for_sorting, ascending):
        pass

    def test_factorize(self, data_for_grouping, na_sentinel=None):
        labels, uniques = pd.factorize(data_for_grouping,
                                       na_sentinel=na_sentinel)
        expected_labels = np.array([0, 0, 1,
                                   1, 2, 2, 0, 3],
                                   dtype=np.intp)
        expected_uniques = data_for_grouping.take([0, 2, 4, 7])

        tm.assert_numpy_array_equal(labels, expected_labels)
        self.assert_extension_array_equal(uniques, expected_uniques)


class TestCasting(BaseRelObject, base.BaseCastingTests):
    pass


class TestGroupby(BaseRelObject, base.BaseGroupbyTests):

    @pytest.mark.xfail(reason="transform fails when __eq__ returns obj")
    def test_groupby_extension_transform(self, data_for_grouping):
        pass

    @pytest.mark.xfail(reason="apply fails when __eq__ returns obj")
    def test_groupby_extension_apply(self, data_for_grouping, op):
        pass


def test_series_constructor_coerce_data_to_extension_dtype_raises():
    xpr = ("Cannot cast data to extension dtype 'relobj'. Pass the "
           "extension array directly.")
    with tm.assert_raises_regex(ValueError, xpr):
        pd.Series([0, 1, 2], dtype=RelObjectDtype())


def test_series_constructor_with_same_dtype_ok():
    arr = RelObjectArray([10])
    result = pd.Series(arr, dtype=RelObjectDtype())
    expected = pd.Series(arr)
    tm.assert_series_equal(result, expected)


def test_series_constructor_coerce_extension_array_to_dtype_raises():
    arr = RelObjectArray([10])
    xpr = r"Cannot specify a dtype 'float.* \('relobj'\)."

    with tm.assert_raises_regex(ValueError, xpr):
        pd.Series(arr, dtype='float')


def test_dataframe_constructor_with_same_dtype_ok():
    arr = RelObjectArray([10])

    result = pd.DataFrame({"A": arr}, dtype=RelObjectDtype())
    expected = pd.DataFrame({"A": arr})
    tm.assert_frame_equal(result, expected)


def test_dataframe_constructor_with_different_dtype_raises():
    arr = RelObjectArray([10])

    xpr = "Cannot coerce extension array to dtype 'float"
    with tm.assert_raises_regex(ValueError, xpr):
        pd.DataFrame({"A": arr}, dtype='float')


@pytest.mark.parametrize(
    'op, supported',
    [
        ('lt', False),
        ('le', True),
        ('gt', False),
        ('ge', True),
        ('eq', True),
        ('ne', False)])
def test_comparisons(op, supported):
    arr1 = RelObjectArray(make_data())
    arr2 = RelObjectArray(make_data())
    ser1 = pd.Series(arr1)
    ser2 = pd.Series(arr2)
    func = getattr(operator, op)

    nsuppmsg = op + " not supported"
    nocomparemsg = "Cannot compare RelObj to object of different type"

    if supported:
        result = func(ser1, ser2)
        expected = pd.Series([func(a, b) for (a, b) in zip(arr1, arr2)])
        tm.assert_series_equal(result, expected)
    else:
        with tm.assert_raises_regex(Exception, nsuppmsg):
            result = func(ser1, ser2)

    oneval = 10
    if supported:
        etype = TypeError
        msg = nocomparemsg
    else:
        etype = Exception
        msg = nsuppmsg
    with tm.assert_raises_regex(etype, msg):
        result = func(ser1, oneval)

    alist = [i for i in arr2]
    if supported:
        result = func(ser1, alist)
        expected = pd.Series([func(a, b) for (a, b) in zip(arr1, alist)])
        tm.assert_series_equal(result, expected)
    else:
        with tm.assert_raises_regex(Exception, nsuppmsg):
            result = func(ser1, alist)

    if op not in ['eq', 'ne']:
        l2 = list(arr2)
        l2[5] = 'abc'
        with tm.assert_raises_regex(etype, msg):
            func(ser1, "abc")

        with tm.assert_raises_regex(etype, msg):
            func(ser1, l2)
