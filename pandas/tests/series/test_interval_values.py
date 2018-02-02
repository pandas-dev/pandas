import pytest
from pandas import (
    Categorical,
    IntervalIndex,
    interval_range,
    Series)
from pandas.core.dtypes.common import is_categorical_dtype
from pandas.core.indexes.accessors import IntervalAccessor
import pandas.util.testing as tm


class TestSeriesIntervalAccessor(object):

    @pytest.mark.parametrize('prop', IntervalIndex._interval_ops)
    @pytest.mark.parametrize('data', [
        IntervalIndex.from_breaks([0, 1, 3, 6, 10]),
        Categorical(interval_range(0, 3).repeat(2))])
    def test_iv_properties(self, prop, data):
        s = Series(data)
        if is_categorical_dtype(data):
            ii = IntervalIndex(data.get_values())
        else:
            ii = data

        # check values
        result = getattr(s.iv, prop)
        expected = Series(getattr(ii, prop))
        tm.assert_series_equal(result, expected)

        # no modifications
        msg = ('modifications to a property of an IntervalIndex object are '
               'not supported. Change values on the original.')
        with tm.assert_raises_regex(ValueError, msg):
            setattr(s.iv, prop, 1)

    def test_iv_accessor_api(self):
        assert Series.iv is IntervalAccessor

        s = Series(interval_range(0, 5))
        assert isinstance(s.iv, IntervalAccessor)

        invalid = Series(list('abcde'))
        assert not hasattr(invalid, 'iv')

        with tm.assert_raises_regex(AttributeError, "only use .iv accessor"):
            invalid.iv

    def test_no_new_attributes(self):
        s = Series(interval_range(0, 5))
        msg = 'You cannot add any new attribute'
        with tm.assert_raises_regex(AttributeError, msg):
            s.iv.new_attribute = 'foo'
