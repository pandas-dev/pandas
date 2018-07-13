import pytest
import numpy as np

from pandas import Index, Interval, IntervalIndex, date_range, timedelta_range
from pandas.core.arrays import IntervalArray
from pandas.core.dtypes.dtypes import IntervalDtype
from pandas.tests.extension import base
import pandas.util.testing as tm


def make_data():
    N = 100
    left = np.random.uniform(size=N).cumsum()
    right = left + np.random.uniform(size=N)
    return [Interval(l, r) for l, r in zip(left, right)]


@pytest.fixture(params=[
    (Index([0, 2, 4]), Index([1, 3, 5])),
    (Index([0., 1., 2.]), Index([1., 2., 3.])),
    (timedelta_range('0 days', periods=3),
     timedelta_range('1 day', periods=3)),
    (date_range('20170101', periods=3), date_range('20170102', periods=3)),
    (date_range('20170101', periods=3, tz='US/Eastern'),
     date_range('20170102', periods=3, tz='US/Eastern'))],
    ids=lambda x: str(x[0].dtype))
def left_right_dtypes(request):
    """
    Fixture for building an IntervalArray from various dtypes
    """
    return request.param


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
def data_repeated():
    """Return different versions of data for count times"""
    def gen(count):
        for _ in range(count):
            yield IntervalArray(make_data())
    yield gen


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

    def test_array_type_with_arg(self, data, dtype):
        assert dtype.construct_array_type() is IntervalArray


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
    @pytest.mark.parametrize('repeats', [0, 1, 5])
    def test_repeat(self, left_right_dtypes, repeats):
        left, right = left_right_dtypes
        result = IntervalArray.from_arrays(left, right).repeat(repeats)
        expected = IntervalArray.from_arrays(
            left.repeat(repeats), right.repeat(repeats))
        tm.assert_extension_array_equal(result, expected)

    @pytest.mark.parametrize('bad_repeats, msg', [
        (-1, 'negative dimensions are not allowed'),
        ('foo', r'invalid literal for (int|long)\(\) with base 10')])
    def test_repeat_errors(self, bad_repeats, msg):
        array = IntervalArray.from_breaks(range(4))
        with tm.assert_raises_regex(ValueError, msg):
            array.repeat(bad_repeats)

    @pytest.mark.parametrize('new_closed', [
        'left', 'right', 'both', 'neither'])
    def test_set_closed(self, closed, new_closed):
        # GH 21670
        array = IntervalArray.from_breaks(range(10), closed=closed)
        result = array.set_closed(new_closed)
        expected = IntervalArray.from_breaks(range(10), closed=new_closed)
        tm.assert_extension_array_equal(result, expected)

    @pytest.mark.skip(reason='addition is not defined for intervals')
    def test_combine_add(self, data_repeated):
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

    def test_set_na(self, left_right_dtypes):
        left, right = left_right_dtypes
        result = IntervalArray.from_arrays(left, right)
        result[0] = np.nan

        expected_left = Index([left._na_value] + list(left[1:]))
        expected_right = Index([right._na_value] + list(right[1:]))
        expected = IntervalArray.from_arrays(expected_left, expected_right)

        self.assert_extension_array_equal(result, expected)


def test_repr_matches():
    idx = IntervalIndex.from_breaks([1, 2, 3])
    a = repr(idx)
    b = repr(idx.values)
    assert a.replace("Index", "Array") == b
