""" test Series-boxed versions of the scalar Timestamp """

import operator
import numpy as np

import pandas.util.testing as tm

from pandas.util.testing import assert_series_equal
from pandas import Timestamp, date_range, Series


class TestSeriesTimestamps(object):
    def test_series_box_timestamp(self):
        rng = date_range('20090415', '20090519', freq='B')
        s = Series(rng)

        assert isinstance(s[5], Timestamp)

        rng = date_range('20090415', '20090519', freq='B')
        s = Series(rng, index=rng)
        assert isinstance(s[5], Timestamp)

        assert isinstance(s.iat[5], Timestamp)

    def test_timestamp_equality(self):
        # GH 11034
        s = Series([Timestamp('2000-01-29 01:59:00'), 'NaT'])
        result = s != s
        assert_series_equal(result, Series([False, True]))
        result = s != s[0]
        assert_series_equal(result, Series([False, True]))
        result = s != s[1]
        assert_series_equal(result, Series([True, True]))

        result = s == s
        assert_series_equal(result, Series([True, False]))
        result = s == s[0]
        assert_series_equal(result, Series([True, False]))
        result = s == s[1]
        assert_series_equal(result, Series([False, False]))

    def test_timestamp_and_series(self):
        timestamp_series = Series(date_range('2014-03-17', periods=2, freq='D',
                                             tz='US/Eastern'))
        first_timestamp = timestamp_series[0]

        delta_series = Series([np.timedelta64(0, 'D'), np.timedelta64(1, 'D')])
        assert_series_equal(timestamp_series - first_timestamp, delta_series)
        assert_series_equal(first_timestamp - timestamp_series, -delta_series)

    def test_timestamp_compare_series(self):
        # make sure we can compare Timestamps on the right AND left hand side
        # GH4982
        s = Series(date_range('20010101', periods=10), name='dates')
        s_nat = s.copy(deep=True)

        s[0] = Timestamp('nat')
        s[3] = Timestamp('nat')

        ops = {'lt': 'gt', 'le': 'ge', 'eq': 'eq', 'ne': 'ne'}

        for left, right in ops.items():
            left_f = getattr(operator, left)
            right_f = getattr(operator, right)

            # no nats
            expected = left_f(s, Timestamp('20010109'))
            result = right_f(Timestamp('20010109'), s)
            tm.assert_series_equal(result, expected)

            # nats
            expected = left_f(s, Timestamp('nat'))
            result = right_f(Timestamp('nat'), s)
            tm.assert_series_equal(result, expected)

            # compare to timestamp with series containing nats
            expected = left_f(s_nat, Timestamp('20010109'))
            result = right_f(Timestamp('20010109'), s_nat)
            tm.assert_series_equal(result, expected)

            # compare to nat with series containing nats
            expected = left_f(s_nat, Timestamp('nat'))
            result = right_f(Timestamp('nat'), s_nat)
            tm.assert_series_equal(result, expected)
