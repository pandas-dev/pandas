import pytest
import numpy as np

from pandas import Interval, IntervalIndex, date_range, timedelta_range
from pandas.core.arrays import IntervalArray
from pandas.core.dtypes.dtypes import IntervalDtype
from pandas.tests.extension import base
import pandas.util.testing as tm


def make_data():
    N = 100
    left = np.random.uniform(size=N).cumsum()
    right = left + np.random.uniform(size=N)
    return [Interval(l, r) for l, r in zip(left, right)]


@pytest.fixture
def dtype():
    return IntervalDtype()


@pytest.fixture
def data():
    """Length-100 PeriodArray for semantics test."""
    return IntervalArray(make_data())


@pytest.fixture
def data_missing():
    """Length 2 array with [NA, Valid]"""
    return IntervalArray.from_tuples([None, (0, 1)])


@pytest.fixture
def data_for_sorting():
    return IntervalArray.from_tuples([(1, 2), (2, 3), (0, 1)])


@pytest.fixture
def data_missing_for_sorting():
    return IntervalArray.from_tuples([(1, 2), None, (0, 1)])


@pytest.fixture
def na_value():
    return np.nan


@pytest.fixture
def data_for_grouping():
    a = (0, 1)
    b = (1, 2)
    c = (2, 3)
    return IntervalArray.from_tuples([b, b, None, None, a, a, b, c])


class BaseInterval(object):
    pass


class TestDtype(BaseInterval, base.BaseDtypeTests):
    pass


class TestCasting(BaseInterval, base.BaseCastingTests):
    pass


class TestConstructors(BaseInterval, base.BaseConstructorsTests):
    pass


class TestGetitem(BaseInterval, base.BaseGetitemTests):
    pass


class TestGrouping(BaseInterval, base.BaseGroupbyTests):
    pass


class TestInterface(BaseInterval, base.BaseInterfaceTests):
    pass


class TestMethods(BaseInterval, base.BaseMethodsTests):
    pass


class TestMissing(BaseInterval, base.BaseMissingTests):
    # Index.fillna only accepts scalar `value`, so we have to skip all
    # non-scalar fill tests.
    unsupported_fill = pytest.mark.skip("Unsupported fillna option.")

    @unsupported_fill
    def test_fillna_limit_pad(self):
        pass

    @unsupported_fill
    def test_fillna_series_method(self):
        pass

    @unsupported_fill
    def test_fillna_limit_backfill(self):
        pass

    @unsupported_fill
    def test_fillna_series(self):
        pass

    def test_non_scalar_raises(self, data_missing):
        msg = "Got a 'list' instead."
        with tm.assert_raises_regex(TypeError, msg):
            data_missing.fillna([1, 1])


class TestReshaping(BaseInterval, base.BaseReshapingTests):
    pass


class TestSetitem(BaseInterval, base.BaseSetitemTests):

    @pytest.mark.parametrize('left, right', [
        (np.arange(3.0), np.arange(1.0, 4.0)),
        ([0, 2, 4], [1, 3, 5]),
        (timedelta_range('0 days', periods=3),
         timedelta_range('1 day', periods=3)),
        (date_range('20170101', periods=3), date_range('20170102', periods=3)),
        pytest.param(date_range('20170101', periods=3, tz='US/Eastern'),
                     date_range('20170102', periods=3, tz='US/Eastern'),
                     marks=pytest.mark.xfail(reason='fixed after rebase?'))])
    def test_set_na(self, left, right):
        result = IntervalArray.from_arrays(left, right)
        result[0] = np.nan

        expected = IntervalArray.from_arrays(
            [np.nan] + list(left[1:]), [np.nan] + list(right[1:]))

        self.assert_extension_array_equal(result, expected)


def test_repr_matches():
    idx = IntervalIndex.from_breaks([1, 2, 3])
    a = repr(idx)
    b = repr(idx.values)
    assert a.replace("Index", "Array") == b
