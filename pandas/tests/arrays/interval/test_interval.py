import numpy as np
import pytest

import pandas as pd
from pandas import Index, Interval, IntervalIndex, date_range, timedelta_range
from pandas.core.arrays import IntervalArray
import pandas.util.testing as tm


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


class TestMethods:

    @pytest.mark.parametrize('new_closed', [
        'left', 'right', 'both', 'neither'])
    def test_set_closed(self, closed, new_closed):
        # GH 21670
        array = IntervalArray.from_breaks(range(10), closed=closed)
        result = array.set_closed(new_closed)
        expected = IntervalArray.from_breaks(range(10), closed=new_closed)
        tm.assert_extension_array_equal(result, expected)

    @pytest.mark.parametrize('other', [
        Interval(0, 1, closed='right'),
        IntervalArray.from_breaks([1, 2, 3, 4], closed='right'),
    ])
    def test_where_raises(self, other):
        ser = pd.Series(IntervalArray.from_breaks([1, 2, 3, 4],
                                                  closed='left'))
        match = "'value.closed' is 'right', expected 'left'."
        with pytest.raises(ValueError, match=match):
            ser.where([True, False, True], other=other)


class TestSetitem:

    def test_set_na(self, left_right_dtypes):
        left, right = left_right_dtypes
        result = IntervalArray.from_arrays(left, right)
        result[0] = np.nan

        expected_left = Index([left._na_value] + list(left[1:]))
        expected_right = Index([right._na_value] + list(right[1:]))
        expected = IntervalArray.from_arrays(expected_left, expected_right)

        tm.assert_extension_array_equal(result, expected)


def test_repr_matches():
    idx = IntervalIndex.from_breaks([1, 2, 3])
    a = repr(idx)
    b = repr(idx.values)
    assert a.replace("Index", "Array") == b


def test_point_interval():
    match = "both/neither sides must be closed when left == right"
    for closed in ('left', 'right'):
        with pytest.raises(ValueError, match=match):
            pd.Interval(0, 0, closed)

    pd.Interval(0, 0, 'neither')  # no exception
    pd.Interval(0, 0, 'both')  # no exception

    assert not pd.Interval(0, 1, "left").overlaps(pd.Interval(0, 0, "neither"))
    assert pd.Interval(0, 1, "left").overlaps(pd.Interval(0, 0, "both"))
    assert not pd.Interval(0, 1, "right").overlaps(
        pd.Interval(1, 1, "neither"))
    assert pd.Interval(0, 1, "right").overlaps(pd.Interval(1, 1, "both"))
    assert not pd.Interval(0, 1, "both").overlaps(pd.Interval(0, 0, "neither"))
    assert not pd.Interval(0, 1, "both").overlaps(pd.Interval(0, 0, "neither"))
    assert not pd.Interval(0, 1, "neither").overlaps(pd.Interval(1, 1, "both"))
    assert not pd.Interval(0, 1, "neither").overlaps(pd.Interval(1, 1, "both"))
