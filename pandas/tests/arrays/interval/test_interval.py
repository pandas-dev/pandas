# -*- coding: utf-8 -*-
import numpy as np

from pandas import Index, date_range, option_context, timedelta_range
from pandas.core.arrays import IntervalArray
import pandas.util.testing as tm
import pytest


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


class TestMethods(object):

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


class TestSetitem(object):

    def test_set_na(self, left_right_dtypes):
        left, right = left_right_dtypes
        result = IntervalArray.from_arrays(left, right)
        result[0] = np.nan

        expected_left = Index([left._na_value] + list(left[1:]))
        expected_right = Index([right._na_value] + list(right[1:]))
        expected = IntervalArray.from_arrays(expected_left, expected_right)

        tm.assert_extension_array_equal(result, expected)


def test_repr_small():
    arr = IntervalArray.from_breaks([1, 2, 3])
    result = repr(arr)
    expected = (
        '<IntervalArray>\n'
        '[(1, 2], (2, 3]]\n'
        'Length: 2, dtype: interval[int64]'
    )
    assert result == expected


def test_repr_large():
    arr = IntervalArray.from_breaks([1, 2, 3, 4, 5, 6])
    with option_context('display.max_seq_items', 2):
        result = repr(arr)
    expected = (
        '<IntervalArray>\n'
        '[(1, 2],\n'
        ' ...\n'
        ' (5, 6]] \n'
        'Length: 5, dtype: interval[int64]'
    )
    assert result == expected
