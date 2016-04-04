# coding=utf-8
# pylint: disable-msg=E1101,W0612

from datetime import timedelta

from numpy import nan
import numpy as np
import pandas as pd

from pandas import Series, isnull
from pandas.tseries.index import Timestamp

from pandas.compat import range
from pandas.util.testing import assert_series_equal
import pandas.util.testing as tm

from .common import TestData


class TestSeriesMissingData(TestData, tm.TestCase):

    _multiprocess_can_split_ = True

    def test_timedelta_fillna(self):
        # GH 3371
        s = Series([Timestamp('20130101'), Timestamp('20130101'), Timestamp(
            '20130102'), Timestamp('20130103 9:01:01')])
        td = s.diff()

        # reg fillna
        result = td.fillna(0)
        expected = Series([timedelta(0), timedelta(0), timedelta(1), timedelta(
            days=1, seconds=9 * 3600 + 60 + 1)])
        assert_series_equal(result, expected)

        # interprested as seconds
        result = td.fillna(1)
        expected = Series([timedelta(seconds=1), timedelta(0), timedelta(1),
                           timedelta(days=1, seconds=9 * 3600 + 60 + 1)])
        assert_series_equal(result, expected)

        result = td.fillna(timedelta(days=1, seconds=1))
        expected = Series([timedelta(days=1, seconds=1), timedelta(
            0), timedelta(1), timedelta(days=1, seconds=9 * 3600 + 60 + 1)])
        assert_series_equal(result, expected)

        result = td.fillna(np.timedelta64(int(1e9)))
        expected = Series([timedelta(seconds=1), timedelta(0), timedelta(1),
                           timedelta(days=1, seconds=9 * 3600 + 60 + 1)])
        assert_series_equal(result, expected)

        from pandas import tslib
        result = td.fillna(tslib.NaT)
        expected = Series([tslib.NaT, timedelta(0), timedelta(1),
                           timedelta(days=1, seconds=9 * 3600 + 60 + 1)],
                          dtype='m8[ns]')
        assert_series_equal(result, expected)

        # ffill
        td[2] = np.nan
        result = td.ffill()
        expected = td.fillna(0)
        expected[0] = np.nan
        assert_series_equal(result, expected)

        # bfill
        td[2] = np.nan
        result = td.bfill()
        expected = td.fillna(0)
        expected[2] = timedelta(days=1, seconds=9 * 3600 + 60 + 1)
        assert_series_equal(result, expected)

    def test_datetime64_fillna(self):

        s = Series([Timestamp('20130101'), Timestamp('20130101'), Timestamp(
            '20130102'), Timestamp('20130103 9:01:01')])
        s[2] = np.nan

        # reg fillna
        result = s.fillna(Timestamp('20130104'))
        expected = Series([Timestamp('20130101'), Timestamp(
            '20130101'), Timestamp('20130104'), Timestamp('20130103 9:01:01')])
        assert_series_equal(result, expected)

        from pandas import tslib
        result = s.fillna(tslib.NaT)
        expected = s
        assert_series_equal(result, expected)

        # ffill
        result = s.ffill()
        expected = Series([Timestamp('20130101'), Timestamp(
            '20130101'), Timestamp('20130101'), Timestamp('20130103 9:01:01')])
        assert_series_equal(result, expected)

        # bfill
        result = s.bfill()
        expected = Series([Timestamp('20130101'), Timestamp('20130101'),
                           Timestamp('20130103 9:01:01'), Timestamp(
                               '20130103 9:01:01')])
        assert_series_equal(result, expected)

        # GH 6587
        # make sure that we are treating as integer when filling
        # this also tests inference of a datetime-like with NaT's
        s = Series([pd.NaT, pd.NaT, '2013-08-05 15:30:00.000001'])
        expected = Series(
            ['2013-08-05 15:30:00.000001', '2013-08-05 15:30:00.000001',
             '2013-08-05 15:30:00.000001'], dtype='M8[ns]')
        result = s.fillna(method='backfill')
        assert_series_equal(result, expected)

    def test_datetime64_tz_fillna(self):
        for tz in ['US/Eastern', 'Asia/Tokyo']:
            # DatetimeBlock
            s = Series([Timestamp('2011-01-01 10:00'), pd.NaT, Timestamp(
                '2011-01-03 10:00'), pd.NaT])
            result = s.fillna(pd.Timestamp('2011-01-02 10:00'))
            expected = Series([Timestamp('2011-01-01 10:00'), Timestamp(
                '2011-01-02 10:00'), Timestamp('2011-01-03 10:00'), Timestamp(
                    '2011-01-02 10:00')])
            self.assert_series_equal(expected, result)

            result = s.fillna(pd.Timestamp('2011-01-02 10:00', tz=tz))
            expected = Series([Timestamp('2011-01-01 10:00'), Timestamp(
                '2011-01-02 10:00', tz=tz), Timestamp('2011-01-03 10:00'),
                Timestamp('2011-01-02 10:00', tz=tz)])
            self.assert_series_equal(expected, result)

            result = s.fillna('AAA')
            expected = Series([Timestamp('2011-01-01 10:00'), 'AAA',
                               Timestamp('2011-01-03 10:00'), 'AAA'],
                              dtype=object)
            self.assert_series_equal(expected, result)

            result = s.fillna({1: pd.Timestamp('2011-01-02 10:00', tz=tz),
                               3: pd.Timestamp('2011-01-04 10:00')})
            expected = Series([Timestamp('2011-01-01 10:00'), Timestamp(
                '2011-01-02 10:00', tz=tz), Timestamp('2011-01-03 10:00'),
                Timestamp('2011-01-04 10:00')])
            self.assert_series_equal(expected, result)

            result = s.fillna({1: pd.Timestamp('2011-01-02 10:00'),
                               3: pd.Timestamp('2011-01-04 10:00')})
            expected = Series([Timestamp('2011-01-01 10:00'), Timestamp(
                '2011-01-02 10:00'), Timestamp('2011-01-03 10:00'), Timestamp(
                    '2011-01-04 10:00')])
            self.assert_series_equal(expected, result)

            # DatetimeBlockTZ
            idx = pd.DatetimeIndex(['2011-01-01 10:00', pd.NaT,
                                    '2011-01-03 10:00', pd.NaT], tz=tz)
            s = pd.Series(idx)
            result = s.fillna(pd.Timestamp('2011-01-02 10:00'))
            expected = Series([Timestamp('2011-01-01 10:00', tz=tz), Timestamp(
                '2011-01-02 10:00'), Timestamp('2011-01-03 10:00', tz=tz),
                Timestamp('2011-01-02 10:00')])
            self.assert_series_equal(expected, result)

            result = s.fillna(pd.Timestamp('2011-01-02 10:00', tz=tz))
            idx = pd.DatetimeIndex(['2011-01-01 10:00', '2011-01-02 10:00',
                                    '2011-01-03 10:00', '2011-01-02 10:00'],
                                   tz=tz)
            expected = Series(idx)
            self.assert_series_equal(expected, result)

            result = s.fillna(pd.Timestamp(
                '2011-01-02 10:00', tz=tz).to_pydatetime())
            idx = pd.DatetimeIndex(['2011-01-01 10:00', '2011-01-02 10:00',
                                    '2011-01-03 10:00', '2011-01-02 10:00'],
                                   tz=tz)
            expected = Series(idx)
            self.assert_series_equal(expected, result)

            result = s.fillna('AAA')
            expected = Series([Timestamp('2011-01-01 10:00', tz=tz), 'AAA',
                               Timestamp('2011-01-03 10:00', tz=tz), 'AAA'],
                              dtype=object)
            self.assert_series_equal(expected, result)

            result = s.fillna({1: pd.Timestamp('2011-01-02 10:00', tz=tz),
                               3: pd.Timestamp('2011-01-04 10:00')})
            expected = Series([Timestamp('2011-01-01 10:00', tz=tz), Timestamp(
                '2011-01-02 10:00', tz=tz), Timestamp(
                    '2011-01-03 10:00', tz=tz), Timestamp('2011-01-04 10:00')])
            self.assert_series_equal(expected, result)

            result = s.fillna({1: pd.Timestamp('2011-01-02 10:00', tz=tz),
                               3: pd.Timestamp('2011-01-04 10:00', tz=tz)})
            expected = Series([Timestamp('2011-01-01 10:00', tz=tz), Timestamp(
                '2011-01-02 10:00', tz=tz), Timestamp(
                    '2011-01-03 10:00', tz=tz), Timestamp('2011-01-04 10:00',
                                                          tz=tz)])
            self.assert_series_equal(expected, result)

            # filling with a naive/other zone, coerce to object
            result = s.fillna(Timestamp('20130101'))
            expected = Series([Timestamp('2011-01-01 10:00', tz=tz), Timestamp(
                '2013-01-01'), Timestamp('2011-01-03 10:00', tz=tz), Timestamp(
                    '2013-01-01')])
            self.assert_series_equal(expected, result)

            result = s.fillna(Timestamp('20130101', tz='US/Pacific'))
            expected = Series([Timestamp('2011-01-01 10:00', tz=tz),
                               Timestamp('2013-01-01', tz='US/Pacific'),
                               Timestamp('2011-01-03 10:00', tz=tz),
                               Timestamp('2013-01-01', tz='US/Pacific')])
            self.assert_series_equal(expected, result)

    def test_fillna_int(self):
        s = Series(np.random.randint(-100, 100, 50))
        s.fillna(method='ffill', inplace=True)
        assert_series_equal(s.fillna(method='ffill', inplace=False), s)

    def test_fillna_raise(self):
        s = Series(np.random.randint(-100, 100, 50))
        self.assertRaises(TypeError, s.fillna, [1, 2])
        self.assertRaises(TypeError, s.fillna, (1, 2))

    def test_isnull_for_inf(self):
        s = Series(['a', np.inf, np.nan, 1.0])
        with pd.option_context('mode.use_inf_as_null', True):
            r = s.isnull()
            dr = s.dropna()
        e = Series([False, True, True, False])
        de = Series(['a', 1.0], index=[0, 3])
        tm.assert_series_equal(r, e)
        tm.assert_series_equal(dr, de)

    def test_fillna(self):
        ts = Series([0., 1., 2., 3., 4.], index=tm.makeDateIndex(5))

        self.assert_numpy_array_equal(ts, ts.fillna(method='ffill'))

        ts[2] = np.NaN

        self.assert_numpy_array_equal(ts.fillna(method='ffill'),
                                      [0., 1., 1., 3., 4.])
        self.assert_numpy_array_equal(ts.fillna(method='backfill'),
                                      [0., 1., 3., 3., 4.])

        self.assert_numpy_array_equal(ts.fillna(value=5), [0., 1., 5., 3., 4.])

        self.assertRaises(ValueError, ts.fillna)
        self.assertRaises(ValueError, self.ts.fillna, value=0, method='ffill')

        # GH 5703
        s1 = Series([np.nan])
        s2 = Series([1])
        result = s1.fillna(s2)
        expected = Series([1.])
        assert_series_equal(result, expected)
        result = s1.fillna({})
        assert_series_equal(result, s1)
        result = s1.fillna(Series(()))
        assert_series_equal(result, s1)
        result = s2.fillna(s1)
        assert_series_equal(result, s2)
        result = s1.fillna({0: 1})
        assert_series_equal(result, expected)
        result = s1.fillna({1: 1})
        assert_series_equal(result, Series([np.nan]))
        result = s1.fillna({0: 1, 1: 1})
        assert_series_equal(result, expected)
        result = s1.fillna(Series({0: 1, 1: 1}))
        assert_series_equal(result, expected)
        result = s1.fillna(Series({0: 1, 1: 1}, index=[4, 5]))
        assert_series_equal(result, s1)

        s1 = Series([0, 1, 2], list('abc'))
        s2 = Series([0, np.nan, 2], list('bac'))
        result = s2.fillna(s1)
        expected = Series([0, 0, 2.], list('bac'))
        assert_series_equal(result, expected)

        # limit
        s = Series(np.nan, index=[0, 1, 2])
        result = s.fillna(999, limit=1)
        expected = Series([999, np.nan, np.nan], index=[0, 1, 2])
        assert_series_equal(result, expected)

        result = s.fillna(999, limit=2)
        expected = Series([999, 999, np.nan], index=[0, 1, 2])
        assert_series_equal(result, expected)

        # GH 9043
        # make sure a string representation of int/float values can be filled
        # correctly without raising errors or being converted
        vals = ['0', '1.5', '-0.3']
        for val in vals:
            s = Series([0, 1, np.nan, np.nan, 4], dtype='float64')
            result = s.fillna(val)
            expected = Series([0, 1, val, val, 4], dtype='object')
            assert_series_equal(result, expected)

    def test_fillna_bug(self):
        x = Series([nan, 1., nan, 3., nan], ['z', 'a', 'b', 'c', 'd'])
        filled = x.fillna(method='ffill')
        expected = Series([nan, 1., 1., 3., 3.], x.index)
        assert_series_equal(filled, expected)

        filled = x.fillna(method='bfill')
        expected = Series([1., 1., 3., 3., nan], x.index)
        assert_series_equal(filled, expected)

    def test_fillna_inplace(self):
        x = Series([nan, 1., nan, 3., nan], ['z', 'a', 'b', 'c', 'd'])
        y = x.copy()

        y.fillna(value=0, inplace=True)

        expected = x.fillna(value=0)
        assert_series_equal(y, expected)

    def test_fillna_invalid_method(self):
        try:
            self.ts.fillna(method='ffil')
        except ValueError as inst:
            self.assertIn('ffil', str(inst))

    def test_ffill(self):
        ts = Series([0., 1., 2., 3., 4.], index=tm.makeDateIndex(5))
        ts[2] = np.NaN
        assert_series_equal(ts.ffill(), ts.fillna(method='ffill'))

    def test_bfill(self):
        ts = Series([0., 1., 2., 3., 4.], index=tm.makeDateIndex(5))
        ts[2] = np.NaN
        assert_series_equal(ts.bfill(), ts.fillna(method='bfill'))

    def test_timedelta64_nan(self):

        from pandas import tslib
        td = Series([timedelta(days=i) for i in range(10)])

        # nan ops on timedeltas
        td1 = td.copy()
        td1[0] = np.nan
        self.assertTrue(isnull(td1[0]))
        self.assertEqual(td1[0].value, tslib.iNaT)
        td1[0] = td[0]
        self.assertFalse(isnull(td1[0]))

        td1[1] = tslib.iNaT
        self.assertTrue(isnull(td1[1]))
        self.assertEqual(td1[1].value, tslib.iNaT)
        td1[1] = td[1]
        self.assertFalse(isnull(td1[1]))

        td1[2] = tslib.NaT
        self.assertTrue(isnull(td1[2]))
        self.assertEqual(td1[2].value, tslib.iNaT)
        td1[2] = td[2]
        self.assertFalse(isnull(td1[2]))

        # boolean setting
        # this doesn't work, not sure numpy even supports it
        # result = td[(td>np.timedelta64(timedelta(days=3))) &
        # td<np.timedelta64(timedelta(days=7)))] = np.nan
        # self.assertEqual(isnull(result).sum(), 7)

        # NumPy limitiation =(

        # def test_logical_range_select(self):
        #     np.random.seed(12345)
        #     selector = -0.5 <= self.ts <= 0.5
        #     expected = (self.ts >= -0.5) & (self.ts <= 0.5)
        #     assert_series_equal(selector, expected)

    def test_dropna_empty(self):
        s = Series([])
        self.assertEqual(len(s.dropna()), 0)
        s.dropna(inplace=True)
        self.assertEqual(len(s), 0)

        # invalid axis
        self.assertRaises(ValueError, s.dropna, axis=1)

    def test_datetime64_tz_dropna(self):
        # DatetimeBlock
        s = Series([Timestamp('2011-01-01 10:00'), pd.NaT, Timestamp(
            '2011-01-03 10:00'), pd.NaT])
        result = s.dropna()
        expected = Series([Timestamp('2011-01-01 10:00'),
                           Timestamp('2011-01-03 10:00')], index=[0, 2])
        self.assert_series_equal(result, expected)

        # DatetimeBlockTZ
        idx = pd.DatetimeIndex(['2011-01-01 10:00', pd.NaT,
                                '2011-01-03 10:00', pd.NaT],
                               tz='Asia/Tokyo')
        s = pd.Series(idx)
        self.assertEqual(s.dtype, 'datetime64[ns, Asia/Tokyo]')
        result = s.dropna()
        expected = Series([Timestamp('2011-01-01 10:00', tz='Asia/Tokyo'),
                           Timestamp('2011-01-03 10:00', tz='Asia/Tokyo')],
                          index=[0, 2])
        self.assertEqual(result.dtype, 'datetime64[ns, Asia/Tokyo]')
        self.assert_series_equal(result, expected)

    def test_dropna_no_nan(self):
        for s in [Series([1, 2, 3], name='x'), Series(
                [False, True, False], name='x')]:

            result = s.dropna()
            self.assert_series_equal(result, s)
            self.assertFalse(result is s)

            s2 = s.copy()
            s2.dropna(inplace=True)
            self.assert_series_equal(s2, s)

    def test_valid(self):
        ts = self.ts.copy()
        ts[::2] = np.NaN

        result = ts.valid()
        self.assertEqual(len(result), ts.count())

        tm.assert_dict_equal(result, ts, compare_keys=False)

    def test_isnull(self):
        ser = Series([0, 5.4, 3, nan, -0.001])
        np.array_equal(ser.isnull(),
                       Series([False, False, False, True, False]).values)
        ser = Series(["hi", "", nan])
        np.array_equal(ser.isnull(), Series([False, False, True]).values)

    def test_notnull(self):
        ser = Series([0, 5.4, 3, nan, -0.001])
        np.array_equal(ser.notnull(),
                       Series([True, True, True, False, True]).values)
        ser = Series(["hi", "", nan])
        np.array_equal(ser.notnull(), Series([True, True, False]).values)

    def test_pad_nan(self):
        x = Series([np.nan, 1., np.nan, 3., np.nan], ['z', 'a', 'b', 'c', 'd'],
                   dtype=float)

        x.fillna(method='pad', inplace=True)

        expected = Series([np.nan, 1.0, 1.0, 3.0, 3.0],
                          ['z', 'a', 'b', 'c', 'd'], dtype=float)
        assert_series_equal(x[1:], expected[1:])
        self.assertTrue(np.isnan(x[0]), np.isnan(expected[0]))

    def test_dropna_preserve_name(self):
        self.ts[:5] = np.nan
        result = self.ts.dropna()
        self.assertEqual(result.name, self.ts.name)
        name = self.ts.name
        ts = self.ts.copy()
        ts.dropna(inplace=True)
        self.assertEqual(ts.name, name)

    def test_fill_value_when_combine_const(self):
        # GH12723
        s = Series([0, 1, np.nan, 3, 4, 5])

        exp = s.fillna(0).add(2)
        res = s.add(2, fill_value=0)
        assert_series_equal(res, exp)


if __name__ == '__main__':
    import nose
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   # '--with-coverage', '--cover-package=pandas.core']
                   exit=False)
