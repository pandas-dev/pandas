import numpy as np

import pandas.lib as lib
import pandas.util.testing as tm
from pandas import Float64Index, date_range, Timestamp
from pandas import (Index, DatetimeIndex, datetime, offsets, to_datetime,
                    Series, DataFrame)


class TestDateTimeIndexToJulianDate(tm.TestCase):

    def test_1700(self):
        r1 = Float64Index([2345897.5, 2345898.5, 2345899.5, 2345900.5,
                           2345901.5])
        r2 = date_range(start=Timestamp('1710-10-01'), periods=5,
                        freq='D').to_julian_date()
        self.assertIsInstance(r2, Float64Index)
        tm.assert_index_equal(r1, r2)

    def test_2000(self):
        r1 = Float64Index([2451601.5, 2451602.5, 2451603.5, 2451604.5,
                           2451605.5])
        r2 = date_range(start=Timestamp('2000-02-27'), periods=5,
                        freq='D').to_julian_date()
        self.assertIsInstance(r2, Float64Index)
        tm.assert_index_equal(r1, r2)

    def test_hour(self):
        r1 = Float64Index(
            [2451601.5, 2451601.5416666666666666, 2451601.5833333333333333,
             2451601.625, 2451601.6666666666666666])
        r2 = date_range(start=Timestamp('2000-02-27'), periods=5,
                        freq='H').to_julian_date()
        self.assertIsInstance(r2, Float64Index)
        tm.assert_index_equal(r1, r2)

    def test_minute(self):
        r1 = Float64Index(
            [2451601.5, 2451601.5006944444444444, 2451601.5013888888888888,
             2451601.5020833333333333, 2451601.5027777777777777])
        r2 = date_range(start=Timestamp('2000-02-27'), periods=5,
                        freq='T').to_julian_date()
        self.assertIsInstance(r2, Float64Index)
        tm.assert_index_equal(r1, r2)

    def test_second(self):
        r1 = Float64Index(
            [2451601.5, 2451601.500011574074074, 2451601.5000231481481481,
             2451601.5000347222222222, 2451601.5000462962962962])
        r2 = date_range(start=Timestamp('2000-02-27'), periods=5,
                        freq='S').to_julian_date()
        self.assertIsInstance(r2, Float64Index)
        tm.assert_index_equal(r1, r2)


class TestTimeSeries(tm.TestCase):
    _multiprocess_can_split_ = True

    def test_pass_datetimeindex_to_index(self):
        # Bugs in #1396
        rng = date_range('1/1/2000', '3/1/2000')
        idx = Index(rng, dtype=object)

        expected = Index(rng.to_pydatetime(), dtype=object)

        self.assert_numpy_array_equal(idx.values, expected.values)

    def test_range_edges(self):
        # GH 13672
        idx = DatetimeIndex(start=Timestamp('1970-01-01 00:00:00.000000001'),
                            end=Timestamp('1970-01-01 00:00:00.000000004'),
                            freq='N')
        exp = DatetimeIndex(['1970-01-01 00:00:00.000000001',
                             '1970-01-01 00:00:00.000000002',
                             '1970-01-01 00:00:00.000000003',
                             '1970-01-01 00:00:00.000000004'])
        tm.assert_index_equal(idx, exp)

        idx = DatetimeIndex(start=Timestamp('1970-01-01 00:00:00.000000004'),
                            end=Timestamp('1970-01-01 00:00:00.000000001'),
                            freq='N')
        exp = DatetimeIndex([])
        tm.assert_index_equal(idx, exp)

        idx = DatetimeIndex(start=Timestamp('1970-01-01 00:00:00.000000001'),
                            end=Timestamp('1970-01-01 00:00:00.000000001'),
                            freq='N')
        exp = DatetimeIndex(['1970-01-01 00:00:00.000000001'])
        tm.assert_index_equal(idx, exp)

        idx = DatetimeIndex(start=Timestamp('1970-01-01 00:00:00.000001'),
                            end=Timestamp('1970-01-01 00:00:00.000004'),
                            freq='U')
        exp = DatetimeIndex(['1970-01-01 00:00:00.000001',
                             '1970-01-01 00:00:00.000002',
                             '1970-01-01 00:00:00.000003',
                             '1970-01-01 00:00:00.000004'])
        tm.assert_index_equal(idx, exp)

        idx = DatetimeIndex(start=Timestamp('1970-01-01 00:00:00.001'),
                            end=Timestamp('1970-01-01 00:00:00.004'),
                            freq='L')
        exp = DatetimeIndex(['1970-01-01 00:00:00.001',
                             '1970-01-01 00:00:00.002',
                             '1970-01-01 00:00:00.003',
                             '1970-01-01 00:00:00.004'])
        tm.assert_index_equal(idx, exp)

        idx = DatetimeIndex(start=Timestamp('1970-01-01 00:00:01'),
                            end=Timestamp('1970-01-01 00:00:04'), freq='S')
        exp = DatetimeIndex(['1970-01-01 00:00:01', '1970-01-01 00:00:02',
                             '1970-01-01 00:00:03', '1970-01-01 00:00:04'])
        tm.assert_index_equal(idx, exp)

        idx = DatetimeIndex(start=Timestamp('1970-01-01 00:01'),
                            end=Timestamp('1970-01-01 00:04'), freq='T')
        exp = DatetimeIndex(['1970-01-01 00:01', '1970-01-01 00:02',
                             '1970-01-01 00:03', '1970-01-01 00:04'])
        tm.assert_index_equal(idx, exp)

        idx = DatetimeIndex(start=Timestamp('1970-01-01 01:00'),
                            end=Timestamp('1970-01-01 04:00'), freq='H')
        exp = DatetimeIndex(['1970-01-01 01:00', '1970-01-01 02:00',
                             '1970-01-01 03:00', '1970-01-01 04:00'])
        tm.assert_index_equal(idx, exp)

        idx = DatetimeIndex(start=Timestamp('1970-01-01'),
                            end=Timestamp('1970-01-04'), freq='D')
        exp = DatetimeIndex(['1970-01-01', '1970-01-02',
                             '1970-01-03', '1970-01-04'])
        tm.assert_index_equal(idx, exp)

    def test_date_range_businesshour(self):
        idx = DatetimeIndex(['2014-07-04 09:00', '2014-07-04 10:00',
                             '2014-07-04 11:00',
                             '2014-07-04 12:00', '2014-07-04 13:00',
                             '2014-07-04 14:00',
                             '2014-07-04 15:00', '2014-07-04 16:00'],
                            freq='BH')
        rng = date_range('2014-07-04 09:00', '2014-07-04 16:00', freq='BH')
        tm.assert_index_equal(idx, rng)

        idx = DatetimeIndex(
            ['2014-07-04 16:00', '2014-07-07 09:00'], freq='BH')
        rng = date_range('2014-07-04 16:00', '2014-07-07 09:00', freq='BH')
        tm.assert_index_equal(idx, rng)

        idx = DatetimeIndex(['2014-07-04 09:00', '2014-07-04 10:00',
                             '2014-07-04 11:00',
                             '2014-07-04 12:00', '2014-07-04 13:00',
                             '2014-07-04 14:00',
                             '2014-07-04 15:00', '2014-07-04 16:00',
                             '2014-07-07 09:00', '2014-07-07 10:00',
                             '2014-07-07 11:00',
                             '2014-07-07 12:00', '2014-07-07 13:00',
                             '2014-07-07 14:00',
                             '2014-07-07 15:00', '2014-07-07 16:00',
                             '2014-07-08 09:00', '2014-07-08 10:00',
                             '2014-07-08 11:00',
                             '2014-07-08 12:00', '2014-07-08 13:00',
                             '2014-07-08 14:00',
                             '2014-07-08 15:00', '2014-07-08 16:00'],
                            freq='BH')
        rng = date_range('2014-07-04 09:00', '2014-07-08 16:00', freq='BH')
        tm.assert_index_equal(idx, rng)

    def test_datetimeindex_integers_shift(self):
        rng = date_range('1/1/2000', periods=20)

        result = rng + 5
        expected = rng.shift(5)
        tm.assert_index_equal(result, expected)

        result = rng - 5
        expected = rng.shift(-5)
        tm.assert_index_equal(result, expected)

    def test_datetimeindex_repr_short(self):
        dr = date_range(start='1/1/2012', periods=1)
        repr(dr)

        dr = date_range(start='1/1/2012', periods=2)
        repr(dr)

        dr = date_range(start='1/1/2012', periods=3)
        repr(dr)


class TestDatetime64(tm.TestCase):

    def test_datetimeindex_accessors(self):
        dti = DatetimeIndex(freq='D', start=datetime(1998, 1, 1), periods=365)

        self.assertEqual(dti.year[0], 1998)
        self.assertEqual(dti.month[0], 1)
        self.assertEqual(dti.day[0], 1)
        self.assertEqual(dti.hour[0], 0)
        self.assertEqual(dti.minute[0], 0)
        self.assertEqual(dti.second[0], 0)
        self.assertEqual(dti.microsecond[0], 0)
        self.assertEqual(dti.dayofweek[0], 3)

        self.assertEqual(dti.dayofyear[0], 1)
        self.assertEqual(dti.dayofyear[120], 121)

        self.assertEqual(dti.weekofyear[0], 1)
        self.assertEqual(dti.weekofyear[120], 18)

        self.assertEqual(dti.quarter[0], 1)
        self.assertEqual(dti.quarter[120], 2)

        self.assertEqual(dti.days_in_month[0], 31)
        self.assertEqual(dti.days_in_month[90], 30)

        self.assertEqual(dti.is_month_start[0], True)
        self.assertEqual(dti.is_month_start[1], False)
        self.assertEqual(dti.is_month_start[31], True)
        self.assertEqual(dti.is_quarter_start[0], True)
        self.assertEqual(dti.is_quarter_start[90], True)
        self.assertEqual(dti.is_year_start[0], True)
        self.assertEqual(dti.is_year_start[364], False)
        self.assertEqual(dti.is_month_end[0], False)
        self.assertEqual(dti.is_month_end[30], True)
        self.assertEqual(dti.is_month_end[31], False)
        self.assertEqual(dti.is_month_end[364], True)
        self.assertEqual(dti.is_quarter_end[0], False)
        self.assertEqual(dti.is_quarter_end[30], False)
        self.assertEqual(dti.is_quarter_end[89], True)
        self.assertEqual(dti.is_quarter_end[364], True)
        self.assertEqual(dti.is_year_end[0], False)
        self.assertEqual(dti.is_year_end[364], True)

        # GH 11128
        self.assertEqual(dti.weekday_name[4], u'Monday')
        self.assertEqual(dti.weekday_name[5], u'Tuesday')
        self.assertEqual(dti.weekday_name[6], u'Wednesday')
        self.assertEqual(dti.weekday_name[7], u'Thursday')
        self.assertEqual(dti.weekday_name[8], u'Friday')
        self.assertEqual(dti.weekday_name[9], u'Saturday')
        self.assertEqual(dti.weekday_name[10], u'Sunday')

        self.assertEqual(Timestamp('2016-04-04').weekday_name, u'Monday')
        self.assertEqual(Timestamp('2016-04-05').weekday_name, u'Tuesday')
        self.assertEqual(Timestamp('2016-04-06').weekday_name, u'Wednesday')
        self.assertEqual(Timestamp('2016-04-07').weekday_name, u'Thursday')
        self.assertEqual(Timestamp('2016-04-08').weekday_name, u'Friday')
        self.assertEqual(Timestamp('2016-04-09').weekday_name, u'Saturday')
        self.assertEqual(Timestamp('2016-04-10').weekday_name, u'Sunday')

        self.assertEqual(len(dti.year), 365)
        self.assertEqual(len(dti.month), 365)
        self.assertEqual(len(dti.day), 365)
        self.assertEqual(len(dti.hour), 365)
        self.assertEqual(len(dti.minute), 365)
        self.assertEqual(len(dti.second), 365)
        self.assertEqual(len(dti.microsecond), 365)
        self.assertEqual(len(dti.dayofweek), 365)
        self.assertEqual(len(dti.dayofyear), 365)
        self.assertEqual(len(dti.weekofyear), 365)
        self.assertEqual(len(dti.quarter), 365)
        self.assertEqual(len(dti.is_month_start), 365)
        self.assertEqual(len(dti.is_month_end), 365)
        self.assertEqual(len(dti.is_quarter_start), 365)
        self.assertEqual(len(dti.is_quarter_end), 365)
        self.assertEqual(len(dti.is_year_start), 365)
        self.assertEqual(len(dti.is_year_end), 365)
        self.assertEqual(len(dti.weekday_name), 365)

        dti = DatetimeIndex(freq='BQ-FEB', start=datetime(1998, 1, 1),
                            periods=4)

        self.assertEqual(sum(dti.is_quarter_start), 0)
        self.assertEqual(sum(dti.is_quarter_end), 4)
        self.assertEqual(sum(dti.is_year_start), 0)
        self.assertEqual(sum(dti.is_year_end), 1)

        # Ensure is_start/end accessors throw ValueError for CustomBusinessDay,
        # CBD requires np >= 1.7
        bday_egypt = offsets.CustomBusinessDay(weekmask='Sun Mon Tue Wed Thu')
        dti = date_range(datetime(2013, 4, 30), periods=5, freq=bday_egypt)
        self.assertRaises(ValueError, lambda: dti.is_month_start)

        dti = DatetimeIndex(['2000-01-01', '2000-01-02', '2000-01-03'])

        self.assertEqual(dti.is_month_start[0], 1)

        tests = [
            (Timestamp('2013-06-01', freq='M').is_month_start, 1),
            (Timestamp('2013-06-01', freq='BM').is_month_start, 0),
            (Timestamp('2013-06-03', freq='M').is_month_start, 0),
            (Timestamp('2013-06-03', freq='BM').is_month_start, 1),
            (Timestamp('2013-02-28', freq='Q-FEB').is_month_end, 1),
            (Timestamp('2013-02-28', freq='Q-FEB').is_quarter_end, 1),
            (Timestamp('2013-02-28', freq='Q-FEB').is_year_end, 1),
            (Timestamp('2013-03-01', freq='Q-FEB').is_month_start, 1),
            (Timestamp('2013-03-01', freq='Q-FEB').is_quarter_start, 1),
            (Timestamp('2013-03-01', freq='Q-FEB').is_year_start, 1),
            (Timestamp('2013-03-31', freq='QS-FEB').is_month_end, 1),
            (Timestamp('2013-03-31', freq='QS-FEB').is_quarter_end, 0),
            (Timestamp('2013-03-31', freq='QS-FEB').is_year_end, 0),
            (Timestamp('2013-02-01', freq='QS-FEB').is_month_start, 1),
            (Timestamp('2013-02-01', freq='QS-FEB').is_quarter_start, 1),
            (Timestamp('2013-02-01', freq='QS-FEB').is_year_start, 1),
            (Timestamp('2013-06-30', freq='BQ').is_month_end, 0),
            (Timestamp('2013-06-30', freq='BQ').is_quarter_end, 0),
            (Timestamp('2013-06-30', freq='BQ').is_year_end, 0),
            (Timestamp('2013-06-28', freq='BQ').is_month_end, 1),
            (Timestamp('2013-06-28', freq='BQ').is_quarter_end, 1),
            (Timestamp('2013-06-28', freq='BQ').is_year_end, 0),
            (Timestamp('2013-06-30', freq='BQS-APR').is_month_end, 0),
            (Timestamp('2013-06-30', freq='BQS-APR').is_quarter_end, 0),
            (Timestamp('2013-06-30', freq='BQS-APR').is_year_end, 0),
            (Timestamp('2013-06-28', freq='BQS-APR').is_month_end, 1),
            (Timestamp('2013-06-28', freq='BQS-APR').is_quarter_end, 1),
            (Timestamp('2013-03-29', freq='BQS-APR').is_year_end, 1),
            (Timestamp('2013-11-01', freq='AS-NOV').is_year_start, 1),
            (Timestamp('2013-10-31', freq='AS-NOV').is_year_end, 1),
            (Timestamp('2012-02-01').days_in_month, 29),
            (Timestamp('2013-02-01').days_in_month, 28)]

        for ts, value in tests:
            self.assertEqual(ts, value)

    def test_datetimeindex_diff(self):
        dti1 = DatetimeIndex(freq='Q-JAN', start=datetime(1997, 12, 31),
                             periods=100)
        dti2 = DatetimeIndex(freq='Q-JAN', start=datetime(1997, 12, 31),
                             periods=98)
        self.assertEqual(len(dti1.difference(dti2)), 2)

    def test_nanosecond_field(self):
        dti = DatetimeIndex(np.arange(10))

        self.assert_numpy_array_equal(dti.nanosecond,
                                      np.arange(10, dtype=np.int32))

    def test_datetimeindex_constructor(self):
        arr = ['1/1/2005', '1/2/2005', 'Jn 3, 2005', '2005-01-04']
        self.assertRaises(Exception, DatetimeIndex, arr)

        arr = ['1/1/2005', '1/2/2005', '1/3/2005', '2005-01-04']
        idx1 = DatetimeIndex(arr)

        arr = [datetime(2005, 1, 1), '1/2/2005', '1/3/2005', '2005-01-04']
        idx2 = DatetimeIndex(arr)

        arr = [lib.Timestamp(datetime(2005, 1, 1)), '1/2/2005', '1/3/2005',
               '2005-01-04']
        idx3 = DatetimeIndex(arr)

        arr = np.array(['1/1/2005', '1/2/2005', '1/3/2005',
                        '2005-01-04'], dtype='O')
        idx4 = DatetimeIndex(arr)

        arr = to_datetime(['1/1/2005', '1/2/2005', '1/3/2005', '2005-01-04'])
        idx5 = DatetimeIndex(arr)

        arr = to_datetime(['1/1/2005', '1/2/2005', 'Jan 3, 2005', '2005-01-04'
                           ])
        idx6 = DatetimeIndex(arr)

        idx7 = DatetimeIndex(['12/05/2007', '25/01/2008'], dayfirst=True)
        idx8 = DatetimeIndex(['2007/05/12', '2008/01/25'], dayfirst=False,
                             yearfirst=True)
        tm.assert_index_equal(idx7, idx8)

        for other in [idx2, idx3, idx4, idx5, idx6]:
            self.assertTrue((idx1.values == other.values).all())

        sdate = datetime(1999, 12, 25)
        edate = datetime(2000, 1, 1)
        idx = DatetimeIndex(start=sdate, freq='1B', periods=20)
        self.assertEqual(len(idx), 20)
        self.assertEqual(idx[0], sdate + 0 * offsets.BDay())
        self.assertEqual(idx.freq, 'B')

        idx = DatetimeIndex(end=edate, freq=('D', 5), periods=20)
        self.assertEqual(len(idx), 20)
        self.assertEqual(idx[-1], edate)
        self.assertEqual(idx.freq, '5D')

        idx1 = DatetimeIndex(start=sdate, end=edate, freq='W-SUN')
        idx2 = DatetimeIndex(start=sdate, end=edate,
                             freq=offsets.Week(weekday=6))
        self.assertEqual(len(idx1), len(idx2))
        self.assertEqual(idx1.offset, idx2.offset)

        idx1 = DatetimeIndex(start=sdate, end=edate, freq='QS')
        idx2 = DatetimeIndex(start=sdate, end=edate,
                             freq=offsets.QuarterBegin(startingMonth=1))
        self.assertEqual(len(idx1), len(idx2))
        self.assertEqual(idx1.offset, idx2.offset)

        idx1 = DatetimeIndex(start=sdate, end=edate, freq='BQ')
        idx2 = DatetimeIndex(start=sdate, end=edate,
                             freq=offsets.BQuarterEnd(startingMonth=12))
        self.assertEqual(len(idx1), len(idx2))
        self.assertEqual(idx1.offset, idx2.offset)

    def test_dayfirst(self):
        # GH 5917
        arr = ['10/02/2014', '11/02/2014', '12/02/2014']
        expected = DatetimeIndex([datetime(2014, 2, 10), datetime(2014, 2, 11),
                                  datetime(2014, 2, 12)])
        idx1 = DatetimeIndex(arr, dayfirst=True)
        idx2 = DatetimeIndex(np.array(arr), dayfirst=True)
        idx3 = to_datetime(arr, dayfirst=True)
        idx4 = to_datetime(np.array(arr), dayfirst=True)
        idx5 = DatetimeIndex(Index(arr), dayfirst=True)
        idx6 = DatetimeIndex(Series(arr), dayfirst=True)
        tm.assert_index_equal(expected, idx1)
        tm.assert_index_equal(expected, idx2)
        tm.assert_index_equal(expected, idx3)
        tm.assert_index_equal(expected, idx4)
        tm.assert_index_equal(expected, idx5)
        tm.assert_index_equal(expected, idx6)

    def test_dti_set_index_reindex(self):
        # GH 6631
        df = DataFrame(np.random.random(6))
        idx1 = date_range('2011/01/01', periods=6, freq='M', tz='US/Eastern')
        idx2 = date_range('2013', periods=6, freq='A', tz='Asia/Tokyo')

        df = df.set_index(idx1)
        tm.assert_index_equal(df.index, idx1)
        df = df.reindex(idx2)
        tm.assert_index_equal(df.index, idx2)

        # 11314
        # with tz
        index = date_range(datetime(2015, 10, 1),
                           datetime(2015, 10, 1, 23),
                           freq='H', tz='US/Eastern')
        df = DataFrame(np.random.randn(24, 1), columns=['a'], index=index)
        new_index = date_range(datetime(2015, 10, 2),
                               datetime(2015, 10, 2, 23),
                               freq='H', tz='US/Eastern')

        # TODO: unused?
        result = df.set_index(new_index)  # noqa

        self.assertEqual(new_index.freq, index.freq)

    def test_datetimeindex_union_join_empty(self):
        dti = DatetimeIndex(start='1/1/2001', end='2/1/2001', freq='D')
        empty = Index([])

        result = dti.union(empty)
        tm.assertIsInstance(result, DatetimeIndex)
        self.assertIs(result, result)

        result = dti.join(empty)
        tm.assertIsInstance(result, DatetimeIndex)


class TestTimeSeriesDuplicates(tm.TestCase):
    _multiprocess_can_split_ = True

    def test_recreate_from_data(self):
        freqs = ['M', 'Q', 'A', 'D', 'B', 'BH', 'T', 'S', 'L', 'U', 'H', 'N',
                 'C']

        for f in freqs:
            org = DatetimeIndex(start='2001/02/01 09:00', freq=f, periods=1)
            idx = DatetimeIndex(org, freq=f)
            tm.assert_index_equal(idx, org)

            org = DatetimeIndex(start='2001/02/01 09:00', freq=f,
                                tz='US/Pacific', periods=1)
            idx = DatetimeIndex(org, freq=f, tz='US/Pacific')
            tm.assert_index_equal(idx, org)
