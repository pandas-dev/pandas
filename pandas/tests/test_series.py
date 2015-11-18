# coding=utf-8
# pylint: disable-msg=E1101,W0612

import re
import sys
from datetime import datetime, timedelta
import operator
import string
from inspect import getargspec
from itertools import product, starmap
from distutils.version import LooseVersion
import warnings
import random

import nose

from numpy import nan, inf
import numpy as np
import numpy.ma as ma
import pandas as pd

from pandas import (Index, Series, DataFrame, isnull, notnull, bdate_range,
                    date_range, period_range, timedelta_range, _np_version_under1p8)
from pandas.core.index import MultiIndex
from pandas.core.indexing import IndexingError
from pandas.tseries.period import PeriodIndex
from pandas.tseries.index import Timestamp, DatetimeIndex
from pandas.tseries.tdi import Timedelta, TimedeltaIndex
import pandas.core.common as com
import pandas.core.config as cf
import pandas.lib as lib

import pandas.core.datetools as datetools
import pandas.core.nanops as nanops

from pandas.compat import StringIO, lrange, range, zip, u, OrderedDict, long
from pandas import compat
from pandas.util.testing import (assert_series_equal,
                                 assert_almost_equal,
                                 assert_frame_equal,
                                 assert_index_equal,
                                 ensure_clean)
import pandas.util.testing as tm


#------------------------------------------------------------------------------
# Series test cases

JOIN_TYPES = ['inner', 'outer', 'left', 'right']


class CheckNameIntegration(object):

    _multiprocess_can_split_ = True

    def test_scalarop_preserve_name(self):
        result = self.ts * 2
        self.assertEqual(result.name, self.ts.name)

    def test_copy_name(self):
        result = self.ts.copy()
        self.assertEqual(result.name, self.ts.name)

    def test_copy_index_name_checking(self):
        # don't want to be able to modify the index stored elsewhere after
        # making a copy

        self.ts.index.name = None
        self.assertIsNone(self.ts.index.name)
        self.assertIs(self.ts, self.ts)

        cp = self.ts.copy()
        cp.index.name = 'foo'
        com.pprint_thing(self.ts.index.name)
        self.assertIsNone(self.ts.index.name)

    def test_append_preserve_name(self):
        result = self.ts[:5].append(self.ts[5:])
        self.assertEqual(result.name, self.ts.name)

    def test_dt_namespace_accessor(self):

        # GH 7207
        # test .dt namespace accessor

        ok_for_base = ['year','month','day','hour','minute','second','weekofyear','week','dayofweek','weekday','dayofyear','quarter','freq','days_in_month','daysinmonth']
        ok_for_period = ok_for_base + ['qyear']
        ok_for_period_methods = ['strftime']
        ok_for_dt = ok_for_base + ['date','time','microsecond','nanosecond', 'is_month_start', 'is_month_end', 'is_quarter_start',
                                   'is_quarter_end', 'is_year_start', 'is_year_end', 'tz']
        ok_for_dt_methods = ['to_period','to_pydatetime','tz_localize','tz_convert', 'normalize', 'strftime']
        ok_for_td = ['days','seconds','microseconds','nanoseconds']
        ok_for_td_methods = ['components','to_pytimedelta','total_seconds']

        def get_expected(s, name):
            result = getattr(Index(s._values),prop)
            if isinstance(result, np.ndarray):
                if com.is_integer_dtype(result):
                    result = result.astype('int64')
            elif not com.is_list_like(result):
                return result
            return Series(result,index=s.index)

        def compare(s, name):
            a = getattr(s.dt,prop)
            b = get_expected(s,prop)
            if not (com.is_list_like(a) and com.is_list_like(b)):
                self.assertEqual(a,b)
            else:
                tm.assert_series_equal(a,b)

        # datetimeindex
        for s in [Series(date_range('20130101',periods=5)),
                  Series(date_range('20130101',periods=5,freq='s')),
                  Series(date_range('20130101 00:00:00',periods=5,freq='ms'))]:
            for prop in ok_for_dt:
                # we test freq below
                if prop != 'freq':
                    compare(s, prop)

            for prop in ok_for_dt_methods:
                getattr(s.dt, prop)

            result = s.dt.to_pydatetime()
            self.assertIsInstance(result,np.ndarray)
            self.assertTrue(result.dtype == object)

            result = s.dt.tz_localize('US/Eastern')
            expected = Series(DatetimeIndex(s.values).tz_localize('US/Eastern'),index=s.index)
            tm.assert_series_equal(result, expected)

            tz_result = result.dt.tz
            self.assertEqual(str(tz_result), 'US/Eastern')
            freq_result = s.dt.freq
            self.assertEqual(freq_result, DatetimeIndex(s.values, freq='infer').freq)

            # let's localize, then convert
            result = s.dt.tz_localize('UTC').dt.tz_convert('US/Eastern')
            expected = Series(DatetimeIndex(s.values).tz_localize('UTC').tz_convert('US/Eastern'),index=s.index)
            tm.assert_series_equal(result, expected)

        # datetimeindex with tz
        s = Series(date_range('20130101',periods=5,tz='US/Eastern'))
        for prop in ok_for_dt:

            # we test freq below
            if prop != 'freq':
                compare(s, prop)

        for prop in ok_for_dt_methods:
            getattr(s.dt,prop)

        result = s.dt.to_pydatetime()
        self.assertIsInstance(result,np.ndarray)
        self.assertTrue(result.dtype == object)

        result = s.dt.tz_convert('CET')
        expected = Series(s._values.tz_convert('CET'),index=s.index)
        tm.assert_series_equal(result, expected)

        tz_result = result.dt.tz
        self.assertEqual(str(tz_result), 'CET')
        freq_result = s.dt.freq
        self.assertEqual(freq_result, DatetimeIndex(s.values, freq='infer').freq)

        # timedeltaindex
        for s in [Series(timedelta_range('1 day',periods=5),index=list('abcde')),
                  Series(timedelta_range('1 day 01:23:45',periods=5,freq='s')),
                  Series(timedelta_range('2 days 01:23:45.012345',periods=5,freq='ms'))]:
            for prop in ok_for_td:
                # we test freq below
                if prop != 'freq':
                    compare(s, prop)

            for prop in ok_for_td_methods:
                getattr(s.dt, prop)

            result = s.dt.components
            self.assertIsInstance(result,DataFrame)
            tm.assert_index_equal(result.index,s.index)

            result = s.dt.to_pytimedelta()
            self.assertIsInstance(result,np.ndarray)
            self.assertTrue(result.dtype == object)

            result = s.dt.total_seconds()
            self.assertIsInstance(result,pd.Series)
            self.assertTrue(result.dtype == 'float64')

            freq_result = s.dt.freq
            self.assertEqual(freq_result, TimedeltaIndex(s.values, freq='infer').freq)

        # both
        index = date_range('20130101',periods=3,freq='D')
        s = Series(date_range('20140204',periods=3,freq='s'),index=index)
        tm.assert_series_equal(s.dt.year,Series(np.array([2014,2014,2014],dtype='int64'),index=index))
        tm.assert_series_equal(s.dt.month,Series(np.array([2,2,2],dtype='int64'),index=index))
        tm.assert_series_equal(s.dt.second,Series(np.array([0,1,2],dtype='int64'),index=index))
        tm.assert_series_equal(s.dt.normalize(), pd.Series([s[0]] * 3, index=index))

        # periodindex
        for s in [Series(period_range('20130101',periods=5,freq='D'))]:
            for prop in ok_for_period:
                # we test freq below
                if prop != 'freq':
                    compare(s, prop)

            for prop in ok_for_period_methods:
                getattr(s.dt, prop)

            freq_result = s.dt.freq
            self.assertEqual(freq_result, PeriodIndex(s.values).freq)

        # test limited display api
        def get_dir(s):
            results = [ r for r in s.dt.__dir__() if not r.startswith('_') ]
            return list(sorted(set(results)))

        s = Series(date_range('20130101',periods=5,freq='D'))
        results = get_dir(s)
        tm.assert_almost_equal(results,list(sorted(set(ok_for_dt + ok_for_dt_methods))))

        s = Series(period_range('20130101',periods=5,freq='D').asobject)
        results = get_dir(s)
        tm.assert_almost_equal(results, list(sorted(set(ok_for_period + ok_for_period_methods))))

        # 11295
        # ambiguous time error on the conversions
        s = Series(pd.date_range('2015-01-01', '2016-01-01', freq='T'))
        s = s.dt.tz_localize('UTC').dt.tz_convert('America/Chicago')
        results = get_dir(s)
        tm.assert_almost_equal(results, list(sorted(set(ok_for_dt + ok_for_dt_methods))))
        expected = Series(pd.date_range('2015-01-01',
                                        '2016-01-01',
                                        freq='T',
                                        tz='UTC').tz_convert('America/Chicago'))
        tm.assert_series_equal(s, expected)

        # no setting allowed
        s = Series(date_range('20130101',periods=5,freq='D'))
        with tm.assertRaisesRegexp(ValueError, "modifications"):
            s.dt.hour = 5

        # trying to set a copy
        with pd.option_context('chained_assignment','raise'):
            def f():
                s.dt.hour[0] = 5
            self.assertRaises(com.SettingWithCopyError, f)

    def test_dt_accessor_no_new_attributes(self):
        # https://github.com/pydata/pandas/issues/10673
        s = Series(date_range('20130101',periods=5,freq='D'))
        with tm.assertRaisesRegexp(AttributeError, "You cannot add any new attribute"):
            s.dt.xlabel = "a"

    def test_strftime(self):
        # GH 10086
        s = Series(date_range('20130101', periods=5))
        result = s.dt.strftime('%Y/%m/%d')
        expected = Series(['2013/01/01', '2013/01/02', '2013/01/03', '2013/01/04', '2013/01/05'])
        tm.assert_series_equal(result, expected)

        s = Series(date_range('2015-02-03 11:22:33.4567', periods=5))
        result = s.dt.strftime('%Y/%m/%d %H-%M-%S')
        expected = Series(['2015/02/03 11-22-33', '2015/02/04 11-22-33', '2015/02/05 11-22-33',
                           '2015/02/06 11-22-33', '2015/02/07 11-22-33'])
        tm.assert_series_equal(result, expected)

        s = Series(period_range('20130101', periods=5))
        result = s.dt.strftime('%Y/%m/%d')
        expected = Series(['2013/01/01', '2013/01/02', '2013/01/03', '2013/01/04', '2013/01/05'])
        tm.assert_series_equal(result, expected)

        s = Series(period_range('2015-02-03 11:22:33.4567', periods=5, freq='s'))
        result = s.dt.strftime('%Y/%m/%d %H-%M-%S')
        expected = Series(['2015/02/03 11-22-33', '2015/02/03 11-22-34', '2015/02/03 11-22-35',
                           '2015/02/03 11-22-36', '2015/02/03 11-22-37'])
        tm.assert_series_equal(result, expected)

        s = Series(date_range('20130101', periods=5))
        s.iloc[0] = pd.NaT
        result = s.dt.strftime('%Y/%m/%d')
        expected = Series(['NaT', '2013/01/02', '2013/01/03', '2013/01/04', '2013/01/05'])
        tm.assert_series_equal(result, expected)

        datetime_index = date_range('20150301', periods=5)
        result = datetime_index.strftime("%Y/%m/%d")
        expected = np.array(['2015/03/01', '2015/03/02', '2015/03/03', '2015/03/04', '2015/03/05'], dtype=object)
        self.assert_numpy_array_equal(result, expected)

        period_index = period_range('20150301', periods=5)
        result = period_index.strftime("%Y/%m/%d")
        expected = np.array(['2015/03/01', '2015/03/02', '2015/03/03', '2015/03/04', '2015/03/05'], dtype=object)
        self.assert_numpy_array_equal(result, expected)

        s = Series([datetime(2013, 1, 1, 2, 32, 59), datetime(2013, 1, 2, 14, 32, 1)])
        result = s.dt.strftime('%Y-%m-%d %H:%M:%S')
        expected = Series(["2013-01-01 02:32:59", "2013-01-02 14:32:01"])
        tm.assert_series_equal(result, expected)

        s = Series(period_range('20130101', periods=4, freq='H'))
        result = s.dt.strftime('%Y/%m/%d %H:%M:%S')
        expected = Series(["2013/01/01 00:00:00", "2013/01/01 01:00:00",
                           "2013/01/01 02:00:00", "2013/01/01 03:00:00"])

        s = Series(period_range('20130101', periods=4, freq='L'))
        result = s.dt.strftime('%Y/%m/%d %H:%M:%S.%l')
        expected = Series(["2013/01/01 00:00:00.000", "2013/01/01 00:00:00.001",
                           "2013/01/01 00:00:00.002", "2013/01/01 00:00:00.003"])
        tm.assert_series_equal(result, expected)

    def test_valid_dt_with_missing_values(self):

        from datetime import date, time

        # GH 8689
        s = Series(date_range('20130101',periods=5,freq='D'))
        s.iloc[2] = pd.NaT

        for attr in ['microsecond','nanosecond','second','minute','hour','day']:
            expected = getattr(s.dt,attr).copy()
            expected.iloc[2] = np.nan
            result = getattr(s.dt,attr)
            tm.assert_series_equal(result, expected)

        result = s.dt.date
        expected = Series([date(2013,1,1),date(2013,1,2),np.nan,date(2013,1,4),date(2013,1,5)],dtype='object')
        tm.assert_series_equal(result, expected)

        result = s.dt.time
        expected = Series([time(0),time(0),np.nan,time(0),time(0)],dtype='object')
        tm.assert_series_equal(result, expected)

    def test_dt_accessor_api(self):
        # GH 9322
        from pandas.tseries.common import (CombinedDatetimelikeProperties,
                                           DatetimeProperties)
        self.assertIs(Series.dt, CombinedDatetimelikeProperties)

        s = Series(date_range('2000-01-01', periods=3))
        self.assertIsInstance(s.dt, DatetimeProperties)

        for s in [Series(np.arange(5)),
                  Series(list('abcde')),
                  Series(np.random.randn(5))]:
            with tm.assertRaisesRegexp(AttributeError,
                                       "only use .dt accessor"):
                s.dt
            self.assertFalse(hasattr(s, 'dt'))

    def test_tab_completion(self):
        # GH 9910
        s = Series(list('abcd'))
        # Series of str values should have .str but not .dt/.cat in __dir__
        self.assertTrue('str' in dir(s))
        self.assertTrue('dt' not in dir(s))
        self.assertTrue('cat' not in dir(s))

        # similiarly for .dt
        s = Series(date_range('1/1/2015', periods=5))
        self.assertTrue('dt' in dir(s))
        self.assertTrue('str' not in dir(s))
        self.assertTrue('cat' not in dir(s))

        # similiarly for .cat, but with the twist that str and dt should be there
        # if the categories are of that type
        # first cat and str
        s = Series(list('abbcd'), dtype="category")
        self.assertTrue('cat' in dir(s))
        self.assertTrue('str' in dir(s)) # as it is a string categorical
        self.assertTrue('dt' not in dir(s))

        # similar to cat and str
        s = Series(date_range('1/1/2015', periods=5)).astype("category")
        self.assertTrue('cat' in dir(s))
        self.assertTrue('str' not in dir(s))
        self.assertTrue('dt' in dir(s)) # as it is a datetime categorical

    def test_binop_maybe_preserve_name(self):
        # names match, preserve
        result = self.ts * self.ts
        self.assertEqual(result.name, self.ts.name)
        result = self.ts.mul(self.ts)
        self.assertEqual(result.name, self.ts.name)

        result = self.ts * self.ts[:-2]
        self.assertEqual(result.name, self.ts.name)

        # names don't match, don't preserve
        cp = self.ts.copy()
        cp.name = 'something else'
        result = self.ts + cp
        self.assertIsNone(result.name)
        result = self.ts.add(cp)
        self.assertIsNone(result.name)

        ops = ['add', 'sub', 'mul', 'div', 'truediv', 'floordiv', 'mod', 'pow']
        ops = ops + ['r' + op for op in ops]
        for op in ops:
            # names match, preserve
            s = self.ts.copy()
            result = getattr(s, op)(s)
            self.assertEqual(result.name, self.ts.name)

            # names don't match, don't preserve
            cp = self.ts.copy()
            cp.name = 'changed'
            result = getattr(s, op)(cp)
            self.assertIsNone(result.name)

    def test_combine_first_name(self):
        result = self.ts.combine_first(self.ts[:5])
        self.assertEqual(result.name, self.ts.name)

    def test_combine_first_dt64(self):
        from pandas.tseries.tools import to_datetime
        s0 = to_datetime(Series(["2010", np.NaN]))
        s1 = to_datetime(Series([np.NaN, "2011"]))
        rs = s0.combine_first(s1)
        xp = to_datetime(Series(['2010', '2011']))
        assert_series_equal(rs, xp)

        s0 = to_datetime(Series(["2010", np.NaN]))
        s1 = Series([np.NaN, "2011"])
        rs = s0.combine_first(s1)
        xp = Series([datetime(2010, 1, 1), '2011'])
        assert_series_equal(rs, xp)

    def test_get(self):

        # GH 6383
        s = Series(np.array([43, 48, 60, 48, 50, 51, 50, 45, 57, 48, 56,
                             45, 51, 39, 55, 43, 54, 52, 51, 54]))

        result = s.get(25, 0)
        expected = 0
        self.assertEqual(result,expected)

        s = Series(np.array([43, 48, 60, 48, 50, 51, 50, 45, 57, 48, 56,
                             45, 51, 39, 55, 43, 54, 52, 51, 54]),
                   index=pd.Float64Index([25.0, 36.0, 49.0, 64.0, 81.0, 100.0,
                                          121.0, 144.0, 169.0, 196.0, 1225.0,
                                          1296.0, 1369.0, 1444.0, 1521.0, 1600.0,
                                          1681.0, 1764.0, 1849.0, 1936.0],
                                         dtype='object'))

        result = s.get(25, 0)
        expected = 43
        self.assertEqual(result,expected)

        # GH 7407
        # with a boolean accessor
        df = pd.DataFrame({'i':[0]*3, 'b':[False]*3})
        vc = df.i.value_counts()
        result = vc.get(99,default='Missing')
        self.assertEqual(result,'Missing')

        vc = df.b.value_counts()
        result = vc.get(False,default='Missing')
        self.assertEqual(result,3)

        result = vc.get(True,default='Missing')
        self.assertEqual(result,'Missing')

    def test_delitem(self):

        # GH 5542
        # should delete the item inplace
        s = Series(lrange(5))
        del s[0]

        expected = Series(lrange(1,5),index=lrange(1,5))
        assert_series_equal(s, expected)

        del s[1]
        expected = Series(lrange(2,5),index=lrange(2,5))
        assert_series_equal(s, expected)

        # empty
        s = Series()
        def f():
            del s[0]
        self.assertRaises(KeyError, f)

        # only 1 left, del, add, del
        s = Series(1)
        del s[0]
        assert_series_equal(s, Series(dtype='int64', index=Index([], dtype='int64')))
        s[0] = 1
        assert_series_equal(s, Series(1))
        del s[0]
        assert_series_equal(s, Series(dtype='int64', index=Index([], dtype='int64')))

        # Index(dtype=object)
        s = Series(1, index=['a'])
        del s['a']
        assert_series_equal(s, Series(dtype='int64', index=Index([], dtype='object')))
        s['a'] = 1
        assert_series_equal(s, Series(1, index=['a']))
        del s['a']
        assert_series_equal(s, Series(dtype='int64', index=Index([], dtype='object')))

    def test_getitem_preserve_name(self):
        result = self.ts[self.ts > 0]
        self.assertEqual(result.name, self.ts.name)

        result = self.ts[[0, 2, 4]]
        self.assertEqual(result.name, self.ts.name)

        result = self.ts[5:10]
        self.assertEqual(result.name, self.ts.name)

    def test_getitem_setitem_ellipsis(self):
        s = Series(np.random.randn(10))

        np.fix(s)

        result = s[...]
        assert_series_equal(result, s)

        s[...] = 5
        self.assertTrue((result == 5).all())

    def test_getitem_negative_out_of_bounds(self):
        s = Series(tm.rands_array(5, 10), index=tm.rands_array(10, 10))

        self.assertRaises(IndexError, s.__getitem__, -11)
        self.assertRaises(IndexError, s.__setitem__, -11, 'foo')

    def test_multilevel_name_print(self):
        index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'],
                                   ['one', 'two', 'three']],
                           labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                                   [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                           names=['first', 'second'])
        s = Series(lrange(0, len(index)), index=index, name='sth')
        expected = ["first  second",
                    "foo    one       0",
                    "       two       1",
                    "       three     2",
                    "bar    one       3",
                    "       two       4",
                    "baz    two       5",
                    "       three     6",
                    "qux    one       7",
                    "       two       8",
                    "       three     9",
                    "Name: sth, dtype: int64"]
        expected = "\n".join(expected)
        self.assertEqual(repr(s), expected)

    def test_multilevel_preserve_name(self):
        index = MultiIndex(levels=[['foo', 'bar', 'baz', 'qux'],
                                   ['one', 'two', 'three']],
                           labels=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3],
                                   [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
                           names=['first', 'second'])
        s = Series(np.random.randn(len(index)), index=index, name='sth')

        result = s['foo']
        result2 = s.ix['foo']
        self.assertEqual(result.name, s.name)
        self.assertEqual(result2.name, s.name)

    def test_name_printing(self):
        # test small series
        s = Series([0, 1, 2])
        s.name = "test"
        self.assertIn("Name: test", repr(s))
        s.name = None
        self.assertNotIn("Name:", repr(s))
        # test big series (diff code path)
        s = Series(lrange(0, 1000))
        s.name = "test"
        self.assertIn("Name: test", repr(s))
        s.name = None
        self.assertNotIn("Name:", repr(s))

        s = Series(index=date_range('20010101', '20020101'), name='test')
        self.assertIn("Name: test", repr(s))

    def test_pickle_preserve_name(self):
        unpickled = self._pickle_roundtrip_name(self.ts)
        self.assertEqual(unpickled.name, self.ts.name)

    def _pickle_roundtrip_name(self, obj):

        with ensure_clean() as path:
            obj.to_pickle(path)
            unpickled = pd.read_pickle(path)
            return unpickled

    def test_argsort_preserve_name(self):
        result = self.ts.argsort()
        self.assertEqual(result.name, self.ts.name)

    def test_sort_index_name(self):
        result = self.ts.sort_index(ascending=False)
        self.assertEqual(result.name, self.ts.name)

    def test_to_sparse_pass_name(self):
        result = self.ts.to_sparse()
        self.assertEqual(result.name, self.ts.name)


class TestNanops(tm.TestCase):

    _multiprocess_can_split_ = True

    def test_comparisons(self):
        left = np.random.randn(10)
        right = np.random.randn(10)
        left[:3] = np.nan

        result = nanops.nangt(left, right)
        expected = (left > right).astype('O')
        expected[:3] = np.nan

        assert_almost_equal(result, expected)

        s = Series(['a', 'b', 'c'])
        s2 = Series([False, True, False])

        # it works!
        s == s2
        s2 == s

    def test_sum_zero(self):
        arr = np.array([])
        self.assertEqual(nanops.nansum(arr), 0)

        arr = np.empty((10, 0))
        self.assertTrue((nanops.nansum(arr, axis=1) == 0).all())

        # GH #844
        s = Series([], index=[])
        self.assertEqual(s.sum(), 0)

        df = DataFrame(np.empty((10, 0)))
        self.assertTrue((df.sum(1) == 0).all())

    def test_nansum_buglet(self):
        s = Series([1.0, np.nan], index=[0, 1])
        result = np.nansum(s)
        assert_almost_equal(result, 1)

    def test_overflow(self):
        # GH 6915
        # overflowing on the smaller int dtypes
        for dtype in ['int32','int64']:
            v = np.arange(5000000,dtype=dtype)
            s = Series(v)

            # no bottleneck
            result = s.sum(skipna=False)
            self.assertEqual(int(result),v.sum(dtype='int64'))
            result = s.min(skipna=False)
            self.assertEqual(int(result),0)
            result = s.max(skipna=False)
            self.assertEqual(int(result),v[-1])

            # use bottleneck if available
            result = s.sum()
            self.assertEqual(int(result),v.sum(dtype='int64'))
            result = s.min()
            self.assertEqual(int(result),0)
            result = s.max()
            self.assertEqual(int(result),v[-1])

        for dtype in ['float32', 'float64']:
            v = np.arange(5000000, dtype=dtype)
            s = Series(v)

            # no bottleneck
            result = s.sum(skipna=False)
            self.assertEqual(result, v.sum(dtype=dtype))
            result = s.min(skipna=False)
            self.assertTrue(np.allclose(float(result), 0.0))
            result = s.max(skipna=False)
            self.assertTrue(np.allclose(float(result), v[-1]))

            # use bottleneck if available
            result = s.sum()
            self.assertEqual(result, v.sum(dtype=dtype))
            result = s.min()
            self.assertTrue(np.allclose(float(result), 0.0))
            result = s.max()
            self.assertTrue(np.allclose(float(result), v[-1]))

class SafeForSparse(object):
    pass

_ts = tm.makeTimeSeries()

class TestSeries(tm.TestCase, CheckNameIntegration):

    _multiprocess_can_split_ = True

    def setUp(self):
        import warnings

        self.ts = _ts.copy()
        self.ts.name = 'ts'

        self.series = tm.makeStringSeries()
        self.series.name = 'series'

        self.objSeries = tm.makeObjectSeries()
        self.objSeries.name = 'objects'

        self.empty = Series([], index=[])

    def test_scalar_conversion(self):

        # Pass in scalar is disabled
        scalar = Series(0.5)
        self.assertNotIsInstance(scalar, float)

        # coercion
        self.assertEqual(float(Series([1.])), 1.0)
        self.assertEqual(int(Series([1.])), 1)
        self.assertEqual(long(Series([1.])), 1)


    def test_astype(self):
        s = Series(np.random.randn(5),name='foo')

        for dtype in ['float32','float64','int64','int32']:
            astyped = s.astype(dtype)
            self.assertEqual(astyped.dtype, dtype)
            self.assertEqual(astyped.name, s.name)

    def test_TimeSeries_deprecation(self):

        # deprecation TimeSeries, #10890
        with tm.assert_produces_warning(FutureWarning):
            pd.TimeSeries(1,index=date_range('20130101',periods=3))

    def test_constructor(self):
        # Recognize TimeSeries
        with tm.assert_produces_warning(FutureWarning):
            self.assertTrue(self.ts.is_time_series)
        self.assertTrue(self.ts.index.is_all_dates)

        # Pass in Series
        derived = Series(self.ts)
        with tm.assert_produces_warning(FutureWarning):
            self.assertTrue(derived.is_time_series)
        self.assertTrue(derived.index.is_all_dates)

        self.assertTrue(tm.equalContents(derived.index, self.ts.index))
        # Ensure new index is not created
        self.assertEqual(id(self.ts.index), id(derived.index))

        # Mixed type Series
        mixed = Series(['hello', np.NaN], index=[0, 1])
        self.assertEqual(mixed.dtype, np.object_)
        self.assertIs(mixed[1], np.NaN)

        with tm.assert_produces_warning(FutureWarning):
            self.assertFalse(self.empty.is_time_series)
        self.assertFalse(self.empty.index.is_all_dates)
        with tm.assert_produces_warning(FutureWarning):
            self.assertFalse(Series({}).is_time_series)
        self.assertFalse(Series({}).index.is_all_dates)
        self.assertRaises(Exception, Series, np.random.randn(3, 3),
                          index=np.arange(3))

        mixed.name = 'Series'
        rs = Series(mixed).name
        xp = 'Series'
        self.assertEqual(rs, xp)

        # raise on MultiIndex GH4187
        m = MultiIndex.from_arrays([[1, 2], [3, 4]])
        self.assertRaises(NotImplementedError, Series, m)

    def test_constructor_empty(self):
        empty = Series()
        empty2 = Series([])
        assert_series_equal(empty, empty2, check_index_type=False)

        empty = Series(index=lrange(10))
        empty2 = Series(np.nan, index=lrange(10))
        assert_series_equal(empty, empty2)

    def test_constructor_series(self):
        index1 = ['d', 'b', 'a', 'c']
        index2 = sorted(index1)
        s1 = Series([4, 7, -5, 3], index=index1)
        s2 = Series(s1, index=index2)

        assert_series_equal(s2, s1.sort_index())

    def test_constructor_iterator(self):

        expected = Series(list(range(10)),dtype='int64')
        result = Series(range(10),dtype='int64')
        assert_series_equal(result, expected)

    def test_constructor_generator(self):
        gen = (i for i in range(10))

        result = Series(gen)
        exp = Series(lrange(10))
        assert_series_equal(result, exp)

        gen = (i for i in range(10))
        result = Series(gen, index=lrange(10, 20))
        exp.index = lrange(10, 20)
        assert_series_equal(result, exp)

    def test_constructor_map(self):
        # GH8909
        m = map(lambda x: x, range(10))

        result = Series(m)
        exp = Series(lrange(10))
        assert_series_equal(result, exp)

        m = map(lambda x: x, range(10))
        result = Series(m, index=lrange(10, 20))
        exp.index = lrange(10, 20)
        assert_series_equal(result, exp)

    def test_constructor_categorical(self):
        cat = pd.Categorical([0, 1, 2, 0, 1, 2], ['a', 'b', 'c'], fastpath=True)
        res = Series(cat)
        self.assertTrue(res.values.equals(cat))

    def test_constructor_maskedarray(self):
        data = ma.masked_all((3,), dtype=float)
        result = Series(data)
        expected = Series([nan, nan, nan])
        assert_series_equal(result, expected)

        data[0] = 0.0
        data[2] = 2.0
        index = ['a', 'b', 'c']
        result = Series(data, index=index)
        expected = Series([0.0, nan, 2.0], index=index)
        assert_series_equal(result, expected)

        data[1] = 1.0
        result = Series(data, index=index)
        expected = Series([0.0, 1.0, 2.0], index=index)
        assert_series_equal(result, expected)

        data = ma.masked_all((3,), dtype=int)
        result = Series(data)
        expected = Series([nan, nan, nan], dtype=float)
        assert_series_equal(result, expected)

        data[0] = 0
        data[2] = 2
        index = ['a', 'b', 'c']
        result = Series(data, index=index)
        expected = Series([0, nan, 2], index=index, dtype=float)
        assert_series_equal(result, expected)

        data[1] = 1
        result = Series(data, index=index)
        expected = Series([0, 1, 2], index=index, dtype=int)
        assert_series_equal(result, expected)

        data = ma.masked_all((3,), dtype=bool)
        result = Series(data)
        expected = Series([nan, nan, nan], dtype=object)
        assert_series_equal(result, expected)

        data[0] = True
        data[2] = False
        index = ['a', 'b', 'c']
        result = Series(data, index=index)
        expected = Series([True, nan, False], index=index, dtype=object)
        assert_series_equal(result, expected)

        data[1] = True
        result = Series(data, index=index)
        expected = Series([True, True, False], index=index, dtype=bool)
        assert_series_equal(result, expected)

        from pandas import tslib
        data = ma.masked_all((3,), dtype='M8[ns]')
        result = Series(data)
        expected = Series([tslib.iNaT, tslib.iNaT, tslib.iNaT], dtype='M8[ns]')
        assert_series_equal(result, expected)

        data[0] = datetime(2001, 1, 1)
        data[2] = datetime(2001, 1, 3)
        index = ['a', 'b', 'c']
        result = Series(data, index=index)
        expected = Series([datetime(2001, 1, 1), tslib.iNaT,
                           datetime(2001, 1, 3)], index=index, dtype='M8[ns]')
        assert_series_equal(result, expected)

        data[1] = datetime(2001, 1, 2)
        result = Series(data, index=index)
        expected = Series([datetime(2001, 1, 1), datetime(2001, 1, 2),
                           datetime(2001, 1, 3)], index=index, dtype='M8[ns]')
        assert_series_equal(result, expected)

    def test_constructor_default_index(self):
        s = Series([0, 1, 2])
        assert_almost_equal(s.index, np.arange(3))

    def test_constructor_corner(self):
        df = tm.makeTimeDataFrame()
        objs = [df, df]
        s = Series(objs, index=[0, 1])
        tm.assertIsInstance(s, Series)

    def test_constructor_sanitize(self):
        s = Series(np.array([1., 1., 8.]), dtype='i8')
        self.assertEqual(s.dtype, np.dtype('i8'))

        s = Series(np.array([1., 1., np.nan]), copy=True, dtype='i8')
        self.assertEqual(s.dtype, np.dtype('f8'))

    def test_constructor_pass_none(self):
        s = Series(None, index=lrange(5))
        self.assertEqual(s.dtype, np.float64)

        s = Series(None, index=lrange(5), dtype=object)
        self.assertEqual(s.dtype, np.object_)

        # GH 7431
        # inference on the index
        s = Series(index=np.array([None]))
        expected = Series(index=Index([None]))
        assert_series_equal(s,expected)

    def test_constructor_cast(self):
        self.assertRaises(ValueError, Series, ['a', 'b', 'c'], dtype=float)

    def test_constructor_dtype_nocast(self):
        # 1572
        s = Series([1, 2, 3])

        s2 = Series(s, dtype=np.int64)

        s2[1] = 5
        self.assertEqual(s[1], 5)

    def test_constructor_datelike_coercion(self):

        # GH 9477
        # incorrectly infering on dateimelike looking when object dtype is specified
        s = Series([Timestamp('20130101'),'NOV'],dtype=object)
        self.assertEqual(s.iloc[0],Timestamp('20130101'))
        self.assertEqual(s.iloc[1],'NOV')
        self.assertTrue(s.dtype == object)

        # the dtype was being reset on the slicing and re-inferred to datetime even
        # thought the blocks are mixed
        belly = '216 3T19'.split()
        wing1 = '2T15 4H19'.split()
        wing2 = '416 4T20'.split()
        mat = pd.to_datetime('2016-01-22 2019-09-07'.split())
        df = pd.DataFrame({'wing1':wing1, 'wing2':wing2, 'mat':mat}, index=belly)

        result = df.loc['3T19']
        self.assertTrue(result.dtype == object)
        result = df.loc['216']
        self.assertTrue(result.dtype == object)

    def test_constructor_dtype_datetime64(self):
        import pandas.tslib as tslib

        s = Series(tslib.iNaT, dtype='M8[ns]', index=lrange(5))
        self.assertTrue(isnull(s).all())

        # in theory this should be all nulls, but since
        # we are not specifying a dtype is ambiguous
        s = Series(tslib.iNaT, index=lrange(5))
        self.assertFalse(isnull(s).all())

        s = Series(nan, dtype='M8[ns]', index=lrange(5))
        self.assertTrue(isnull(s).all())

        s = Series([datetime(2001, 1, 2, 0, 0), tslib.iNaT], dtype='M8[ns]')
        self.assertTrue(isnull(s[1]))
        self.assertEqual(s.dtype, 'M8[ns]')

        s = Series([datetime(2001, 1, 2, 0, 0), nan], dtype='M8[ns]')
        self.assertTrue(isnull(s[1]))
        self.assertEqual(s.dtype, 'M8[ns]')

        # GH3416
        dates = [
            np.datetime64(datetime(2013, 1, 1)),
            np.datetime64(datetime(2013, 1, 2)),
            np.datetime64(datetime(2013, 1, 3)),
            ]

        s = Series(dates)
        self.assertEqual(s.dtype, 'M8[ns]')

        s.ix[0] = np.nan
        self.assertEqual(s.dtype, 'M8[ns]')

        # invalid astypes
        for t in ['s', 'D', 'us', 'ms']:
            self.assertRaises(TypeError, s.astype, 'M8[%s]' % t)

        # GH3414 related
        self.assertRaises(TypeError, lambda x: Series(
            Series(dates).astype('int') / 1000000, dtype='M8[ms]'))
        self.assertRaises(
            TypeError, lambda x: Series(dates, dtype='datetime64'))

        # invalid dates can be help as object
        result = Series([datetime(2,1,1)])
        self.assertEqual(result[0], datetime(2,1,1,0,0))

        result = Series([datetime(3000,1,1)])
        self.assertEqual(result[0], datetime(3000,1,1,0,0))

        # don't mix types
        result = Series([ Timestamp('20130101'), 1],index=['a','b'])
        self.assertEqual(result['a'], Timestamp('20130101'))
        self.assertEqual(result['b'], 1)

        # GH6529
        # coerce datetime64 non-ns properly
        dates = date_range('01-Jan-2015', '01-Dec-2015', freq='M')
        values2 = dates.view(np.ndarray).astype('datetime64[ns]')
        expected = Series(values2, dates)

        for dtype in ['s', 'D', 'ms', 'us', 'ns']:
            values1 = dates.view(np.ndarray).astype('M8[{0}]'.format(dtype))
            result = Series(values1, dates)
            assert_series_equal(result,expected)

        # leave datetime.date alone
        dates2 = np.array([d.date() for d in dates.to_pydatetime()],
                          dtype=object)
        series1 = Series(dates2, dates)
        self.assert_numpy_array_equal(series1.values,dates2)
        self.assertEqual(series1.dtype,object)

        # these will correctly infer a datetime
        s = Series([None, pd.NaT, '2013-08-05 15:30:00.000001'])
        self.assertEqual(s.dtype,'datetime64[ns]')
        s = Series([np.nan, pd.NaT, '2013-08-05 15:30:00.000001'])
        self.assertEqual(s.dtype,'datetime64[ns]')
        s = Series([pd.NaT, None, '2013-08-05 15:30:00.000001'])
        self.assertEqual(s.dtype,'datetime64[ns]')
        s = Series([pd.NaT, np.nan, '2013-08-05 15:30:00.000001'])
        self.assertEqual(s.dtype,'datetime64[ns]')

        # tz-aware (UTC and other tz's)
        # GH 8411
        dr = date_range('20130101',periods=3)
        self.assertTrue(Series(dr).iloc[0].tz is None)
        dr = date_range('20130101',periods=3,tz='UTC')
        self.assertTrue(str(Series(dr).iloc[0].tz) == 'UTC')
        dr = date_range('20130101',periods=3,tz='US/Eastern')
        self.assertTrue(str(Series(dr).iloc[0].tz) == 'US/Eastern')

        # non-convertible
        s = Series([1479596223000, -1479590, pd.NaT])
        self.assertTrue(s.dtype == 'object')
        self.assertTrue(s[2] is pd.NaT)
        self.assertTrue('NaT' in str(s))

        # if we passed a NaT it remains
        s = Series([datetime(2010, 1, 1), datetime(2, 1, 1), pd.NaT])
        self.assertTrue(s.dtype == 'object')
        self.assertTrue(s[2] is pd.NaT)
        self.assertTrue('NaT' in str(s))

        # if we passed a nan it remains
        s = Series([datetime(2010, 1, 1), datetime(2, 1, 1), np.nan])
        self.assertTrue(s.dtype == 'object')
        self.assertTrue(s[2] is np.nan)
        self.assertTrue('NaN' in str(s))

    def test_constructor_with_datetime_tz(self):

        # 8260
        # support datetime64 with tz

        dr = date_range('20130101',periods=3,tz='US/Eastern')
        s = Series(dr)
        self.assertTrue(s.dtype.name == 'datetime64[ns, US/Eastern]')
        self.assertTrue(s.dtype == 'datetime64[ns, US/Eastern]')
        self.assertTrue(com.is_datetime64tz_dtype(s.dtype))
        self.assertTrue('datetime64[ns, US/Eastern]' in str(s))

        # export
        result = s.values
        self.assertIsInstance(result, np.ndarray)
        self.assertTrue(result.dtype == 'datetime64[ns]')
        self.assertTrue(dr.equals(pd.DatetimeIndex(result).tz_localize('UTC').tz_convert(tz=s.dt.tz)))

        # indexing
        result = s.iloc[0]
        self.assertEqual(result,Timestamp('2013-01-01 00:00:00-0500', tz='US/Eastern', offset='D'))
        result = s[0]
        self.assertEqual(result,Timestamp('2013-01-01 00:00:00-0500', tz='US/Eastern', offset='D'))

        result = s[Series([True,True,False],index=s.index)]
        assert_series_equal(result,s[0:2])

        result = s.iloc[0:1]
        assert_series_equal(result,Series(dr[0:1]))

        # concat
        result = pd.concat([s.iloc[0:1],s.iloc[1:]])
        assert_series_equal(result,s)

        # astype
        result = s.astype(object)
        expected = Series(DatetimeIndex(s._values).asobject)
        assert_series_equal(result, expected)

        result = Series(s.values).dt.tz_localize('UTC').dt.tz_convert(s.dt.tz)
        assert_series_equal(result, s)

        # astype - datetime64[ns, tz]
        result = Series(s.values).astype('datetime64[ns, US/Eastern]')
        assert_series_equal(result, s)

        result = Series(s.values).astype(s.dtype)
        assert_series_equal(result, s)

        result = s.astype('datetime64[ns, CET]')
        expected = Series(date_range('20130101 06:00:00',periods=3,tz='CET'))
        assert_series_equal(result, expected)

        # short str
        self.assertTrue('datetime64[ns, US/Eastern]' in str(s))

        # formatting with NaT
        result = s.shift()
        self.assertTrue('datetime64[ns, US/Eastern]' in str(result))
        self.assertTrue('NaT' in str(result))

        # long str
        t = Series(date_range('20130101',periods=1000,tz='US/Eastern'))
        self.assertTrue('datetime64[ns, US/Eastern]' in str(t))

        result = pd.DatetimeIndex(s,freq='infer')
        tm.assert_index_equal(result, dr)

        # inference
        s = Series([pd.Timestamp('2013-01-01 13:00:00-0800', tz='US/Pacific'),pd.Timestamp('2013-01-02 14:00:00-0800', tz='US/Pacific')])
        self.assertTrue(s.dtype == 'datetime64[ns, US/Pacific]')
        self.assertTrue(lib.infer_dtype(s) == 'datetime64')

        s = Series([pd.Timestamp('2013-01-01 13:00:00-0800', tz='US/Pacific'),pd.Timestamp('2013-01-02 14:00:00-0800', tz='US/Eastern')])
        self.assertTrue(s.dtype == 'object')
        self.assertTrue(lib.infer_dtype(s) == 'datetime')

    def test_constructor_periodindex(self):
        # GH7932
        # converting a PeriodIndex when put in a Series

        pi = period_range('20130101',periods=5,freq='D')
        s = Series(pi)
        expected = Series(pi.asobject)
        assert_series_equal(s, expected)

    def test_constructor_dict(self):
        d = {'a': 0., 'b': 1., 'c': 2.}
        result = Series(d, index=['b', 'c', 'd', 'a'])
        expected = Series([1, 2, nan, 0], index=['b', 'c', 'd', 'a'])
        assert_series_equal(result, expected)

        pidx = tm.makePeriodIndex(100)
        d = {pidx[0]: 0, pidx[1]: 1}
        result = Series(d, index=pidx)
        expected = Series(np.nan, pidx)
        expected.ix[0] = 0
        expected.ix[1] = 1
        assert_series_equal(result, expected)

    def test_constructor_dict_multiindex(self):
        check = lambda result, expected: tm.assert_series_equal(
            result, expected, check_dtype=True, check_index_type=True,
            check_series_type=True)
        d = {('a', 'a'): 0., ('b', 'a'): 1., ('b', 'c'): 2.}
        _d = sorted(d.items())
        ser = Series(d)
        expected = Series([x[1] for x in _d],
                          index=MultiIndex.from_tuples([x[0] for x in _d]))
        check(ser, expected)

        d['z'] = 111.
        _d.insert(0, ('z', d['z']))
        ser = Series(d)
        expected = Series(
            [x[1] for x in _d],
            index=Index([x[0] for x in _d], tupleize_cols=False))
        ser = ser.reindex(index=expected.index)
        check(ser, expected)

    def test_constructor_subclass_dict(self):
        data = tm.TestSubDict((x, 10.0 * x) for x in range(10))
        series = Series(data)
        refseries = Series(dict(compat.iteritems(data)))
        assert_series_equal(refseries, series)

    def test_constructor_dict_datetime64_index(self):
        # GH 9456

        dates_as_str = ['1984-02-19', '1988-11-06', '1989-12-03', '1990-03-15']
        values = [42544017.198965244, 1234565, 40512335.181958228, -1]

        def create_data(constructor):
            return dict(zip((constructor(x) for x in dates_as_str), values))

        data_datetime64 = create_data(np.datetime64)
        data_datetime = create_data(lambda x: datetime.strptime(x, '%Y-%m-%d'))
        data_Timestamp = create_data(Timestamp)

        expected = Series(values, (Timestamp(x) for x in dates_as_str))

        result_datetime64 = Series(data_datetime64)
        result_datetime = Series(data_datetime)
        result_Timestamp = Series(data_Timestamp)

        assert_series_equal(result_datetime64, expected)
        assert_series_equal(result_datetime, expected)
        assert_series_equal(result_Timestamp, expected)

    def test_orderedDict_ctor(self):
        # GH3283
        import pandas
        import random
        data = OrderedDict([('col%s' % i, random.random()) for i in range(12)])
        s = pandas.Series(data)
        self.assertTrue(all(s.values == list(data.values())))

    def test_orderedDict_subclass_ctor(self):
        # GH3283
        import pandas
        import random

        class A(OrderedDict):
            pass
        data = A([('col%s' % i, random.random()) for i in range(12)])
        s = pandas.Series(data)
        self.assertTrue(all(s.values == list(data.values())))

    def test_constructor_list_of_tuples(self):
        data = [(1, 1), (2, 2), (2, 3)]
        s = Series(data)
        self.assertEqual(list(s), data)

    def test_constructor_tuple_of_tuples(self):
        data = ((1, 1), (2, 2), (2, 3))
        s = Series(data)
        self.assertEqual(tuple(s), data)

    def test_constructor_set(self):
        values = set([1, 2, 3, 4, 5])
        self.assertRaises(TypeError, Series, values)
        values = frozenset(values)
        self.assertRaises(TypeError, Series, values)

    def test_fromDict(self):
        data = {'a': 0, 'b': 1, 'c': 2, 'd': 3}

        series = Series(data)
        self.assertTrue(tm.is_sorted(series.index))

        data = {'a': 0, 'b': '1', 'c': '2', 'd': datetime.now()}
        series = Series(data)
        self.assertEqual(series.dtype, np.object_)

        data = {'a': 0, 'b': '1', 'c': '2', 'd': '3'}
        series = Series(data)
        self.assertEqual(series.dtype, np.object_)

        data = {'a': '0', 'b': '1'}
        series = Series(data, dtype=float)
        self.assertEqual(series.dtype, np.float64)

    def test_setindex(self):
        # wrong type
        series = self.series.copy()
        self.assertRaises(TypeError, setattr, series, 'index', None)

        # wrong length
        series = self.series.copy()
        self.assertRaises(Exception, setattr, series, 'index',
                          np.arange(len(series) - 1))

        # works
        series = self.series.copy()
        series.index = np.arange(len(series))
        tm.assertIsInstance(series.index, Index)

    def test_array_finalize(self):
        pass

    def test_pop(self):
        # GH 6600
        df = DataFrame({
            'A': 0,
            'B': np.arange(5,dtype='int64'),
            'C': 0,
            })
        k = df.iloc[4]

        result = k.pop('B')
        self.assertEqual(result, 4)

        expected = Series([0, 0], index=['A', 'C'], name=4)
        assert_series_equal(k, expected)

    def test_not_hashable(self):
        s_empty = Series()
        s = Series([1])
        self.assertRaises(TypeError, hash, s_empty)
        self.assertRaises(TypeError, hash, s)

    def test_fromValue(self):

        nans = Series(np.NaN, index=self.ts.index)
        self.assertEqual(nans.dtype, np.float_)
        self.assertEqual(len(nans), len(self.ts))

        strings = Series('foo', index=self.ts.index)
        self.assertEqual(strings.dtype, np.object_)
        self.assertEqual(len(strings), len(self.ts))

        d = datetime.now()
        dates = Series(d, index=self.ts.index)
        self.assertEqual(dates.dtype, 'M8[ns]')
        self.assertEqual(len(dates), len(self.ts))

    def test_contains(self):
        tm.assert_contains_all(self.ts.index, self.ts)

    def test_pickle(self):
        unp_series = self._pickle_roundtrip(self.series)
        unp_ts = self._pickle_roundtrip(self.ts)
        assert_series_equal(unp_series, self.series)
        assert_series_equal(unp_ts, self.ts)

    def _pickle_roundtrip(self, obj):

        with ensure_clean() as path:
            obj.to_pickle(path)
            unpickled = pd.read_pickle(path)
            return unpickled

    def test_getitem_get(self):
        idx1 = self.series.index[5]
        idx2 = self.objSeries.index[5]

        self.assertEqual(self.series[idx1], self.series.get(idx1))
        self.assertEqual(self.objSeries[idx2], self.objSeries.get(idx2))

        self.assertEqual(self.series[idx1], self.series[5])
        self.assertEqual(self.objSeries[idx2], self.objSeries[5])

        self.assertEqual(
            self.series.get(-1), self.series.get(self.series.index[-1]))
        self.assertEqual(self.series[5], self.series.get(self.series.index[5]))

        # missing
        d = self.ts.index[0] - datetools.bday
        self.assertRaises(KeyError, self.ts.__getitem__, d)

        # None
        # GH 5652
        for s in [Series(), Series(index=list('abc'))]:
            result = s.get(None)
            self.assertIsNone(result)

    def test_iget(self):

        s = Series(np.random.randn(10), index=lrange(0, 20, 2))

        # 10711, deprecated
        with tm.assert_produces_warning(FutureWarning):
            s.iget(1)

        # 10711, deprecated
        with tm.assert_produces_warning(FutureWarning):
            s.irow(1)

        # 10711, deprecated
        with tm.assert_produces_warning(FutureWarning):
            s.iget_value(1)

        for i in range(len(s)):
            result = s.iloc[i]
            exp = s[s.index[i]]
            assert_almost_equal(result, exp)

        # pass a slice
        result = s.iloc[slice(1, 3)]
        expected = s.ix[2:4]
        assert_series_equal(result, expected)

        # test slice is a view
        result[:] = 0
        self.assertTrue((s[1:3] == 0).all())

        # list of integers
        result = s.iloc[[0, 2, 3, 4, 5]]
        expected = s.reindex(s.index[[0, 2, 3, 4, 5]])
        assert_series_equal(result, expected)

    def test_iget_nonunique(self):
        s = Series([0, 1, 2], index=[0, 1, 0])
        self.assertEqual(s.iloc[2], 2)

    def test_getitem_regression(self):
        s = Series(lrange(5), index=lrange(5))
        result = s[lrange(5)]
        assert_series_equal(result, s)

    def test_getitem_setitem_slice_bug(self):
        s = Series(lrange(10), lrange(10))
        result = s[-12:]
        assert_series_equal(result, s)

        result = s[-7:]
        assert_series_equal(result, s[3:])

        result = s[:-12]
        assert_series_equal(result, s[:0])

        s = Series(lrange(10), lrange(10))
        s[-12:] = 0
        self.assertTrue((s == 0).all())

        s[:-12] = 5
        self.assertTrue((s == 0).all())

    def test_getitem_int64(self):
        idx = np.int64(5)
        self.assertEqual(self.ts[idx], self.ts[5])

    def test_getitem_fancy(self):
        slice1 = self.series[[1, 2, 3]]
        slice2 = self.objSeries[[1, 2, 3]]
        self.assertEqual(self.series.index[2], slice1.index[1])
        self.assertEqual(self.objSeries.index[2], slice2.index[1])
        self.assertEqual(self.series[2], slice1[1])
        self.assertEqual(self.objSeries[2], slice2[1])

    def test_getitem_boolean(self):
        s = self.series
        mask = s > s.median()

        # passing list is OK
        result = s[list(mask)]
        expected = s[mask]
        assert_series_equal(result, expected)
        self.assert_numpy_array_equal(result.index, s.index[mask])

    def test_getitem_boolean_empty(self):
        s = Series([], dtype=np.int64)
        s.index.name = 'index_name'
        s = s[s.isnull()]
        self.assertEqual(s.index.name, 'index_name')
        self.assertEqual(s.dtype, np.int64)

        # GH5877
        # indexing with empty series
        s = Series(['A', 'B'])
        expected = Series(np.nan,index=['C'],dtype=object)
        result = s[Series(['C'], dtype=object)]
        assert_series_equal(result, expected)

        s = Series(['A', 'B'])
        expected = Series(dtype=object, index=Index([], dtype='int64'))
        result = s[Series([], dtype=object)]
        assert_series_equal(result, expected)

        # invalid because of the boolean indexer
        # that's empty or not-aligned
        def f():
            s[Series([], dtype=bool)]
        self.assertRaises(IndexingError, f)

        def f():
            s[Series([True], dtype=bool)]
        self.assertRaises(IndexingError, f)

    def test_getitem_generator(self):
        gen = (x > 0 for x in self.series)
        result = self.series[gen]
        result2 = self.series[iter(self.series > 0)]
        expected = self.series[self.series > 0]
        assert_series_equal(result, expected)
        assert_series_equal(result2, expected)

    def test_getitem_boolean_object(self):
        # using column from DataFrame

        s = self.series
        mask = s > s.median()
        omask = mask.astype(object)

        # getitem
        result = s[omask]
        expected = s[mask]
        assert_series_equal(result, expected)

        # setitem
        s2  = s.copy()
        cop = s.copy()
        cop[omask] = 5
        s2[mask] = 5
        assert_series_equal(cop, s2)

        # nans raise exception
        omask[5:10] = np.nan
        self.assertRaises(Exception, s.__getitem__, omask)
        self.assertRaises(Exception, s.__setitem__, omask, 5)

    def test_getitem_setitem_boolean_corner(self):
        ts = self.ts
        mask_shifted = ts.shift(1, freq=datetools.bday) > ts.median()

        # these used to raise...??

        self.assertRaises(Exception, ts.__getitem__, mask_shifted)
        self.assertRaises(Exception, ts.__setitem__, mask_shifted, 1)
        #ts[mask_shifted]
        #ts[mask_shifted] = 1

        self.assertRaises(Exception, ts.ix.__getitem__, mask_shifted)
        self.assertRaises(Exception, ts.ix.__setitem__, mask_shifted, 1)
        #ts.ix[mask_shifted]
        #ts.ix[mask_shifted] = 2

    def test_getitem_setitem_slice_integers(self):
        s = Series(np.random.randn(8), index=[2, 4, 6, 8, 10, 12, 14, 16])

        result = s[:4]
        expected = s.reindex([2, 4, 6, 8])
        assert_series_equal(result, expected)

        s[:4] = 0
        self.assertTrue((s[:4] == 0).all())
        self.assertTrue(not (s[4:] == 0).any())

    def test_getitem_out_of_bounds(self):
        # don't segfault, GH #495
        self.assertRaises(IndexError, self.ts.__getitem__, len(self.ts))

        # GH #917
        s = Series([])
        self.assertRaises(IndexError, s.__getitem__, -1)

    def test_getitem_setitem_integers(self):
        # caused bug without test
        s = Series([1, 2, 3], ['a', 'b', 'c'])

        self.assertEqual(s.ix[0], s['a'])
        s.ix[0] = 5
        self.assertAlmostEqual(s['a'], 5)

    def test_getitem_box_float64(self):
        value = self.ts[5]
        tm.assertIsInstance(value, np.float64)

    def test_getitem_ambiguous_keyerror(self):
        s = Series(lrange(10), index=lrange(0, 20, 2))
        self.assertRaises(KeyError, s.__getitem__, 1)
        self.assertRaises(KeyError, s.ix.__getitem__, 1)

    def test_getitem_unordered_dup(self):
        obj = Series(lrange(5), index=['c', 'a', 'a', 'b', 'b'])
        self.assertTrue(np.isscalar(obj['c']))
        self.assertEqual(obj['c'], 0)

    def test_getitem_dups_with_missing(self):

        # breaks reindex, so need to use .ix internally
        # GH 4246
        s = Series([1, 2, 3, 4], ['foo', 'bar', 'foo', 'bah'])
        expected = s.ix[['foo', 'bar', 'bah', 'bam']]
        result = s[['foo', 'bar', 'bah', 'bam']]
        assert_series_equal(result, expected)

    def test_getitem_dups(self):
        s = Series(range(5),index=['A','A','B','C','C'],dtype=np.int64)
        expected = Series([3,4],index=['C','C'],dtype=np.int64)
        result = s['C']
        assert_series_equal(result, expected)

    def test_getitem_dataframe(self):
        rng = list(range(10))
        s   = pd.Series(10, index=rng)
        df  = pd.DataFrame(rng, index=rng)
        self.assertRaises(TypeError, s.__getitem__, df>5)

    def test_setitem_ambiguous_keyerror(self):
        s = Series(lrange(10), index=lrange(0, 20, 2))

        # equivalent of an append
        s2 = s.copy()
        s2[1] = 5
        expected = s.append(Series([5],index=[1]))
        assert_series_equal(s2,expected)

        s2 = s.copy()
        s2.ix[1] = 5
        expected = s.append(Series([5],index=[1]))
        assert_series_equal(s2,expected)

    def test_setitem_float_labels(self):
        # note labels are floats
        s = Series(['a', 'b', 'c'], index=[0, 0.5, 1])
        tmp = s.copy()

        s.ix[1] = 'zoo'
        tmp.iloc[2] = 'zoo'

        assert_series_equal(s, tmp)

    def test_slice(self):
        numSlice = self.series[10:20]
        numSliceEnd = self.series[-10:]
        objSlice = self.objSeries[10:20]

        self.assertNotIn(self.series.index[9], numSlice.index)
        self.assertNotIn(self.objSeries.index[9], objSlice.index)

        self.assertEqual(len(numSlice), len(numSlice.index))
        self.assertEqual(self.series[numSlice.index[0]],
                         numSlice[numSlice.index[0]])

        self.assertEqual(numSlice.index[1], self.series.index[11])

        self.assertTrue(tm.equalContents(numSliceEnd,
                                         np.array(self.series)[-10:]))

        # test return view
        sl = self.series[10:20]
        sl[:] = 0
        self.assertTrue((self.series[10:20] == 0).all())

    def test_slice_can_reorder_not_uniquely_indexed(self):
        s = Series(1, index=['a', 'a', 'b', 'b', 'c'])
        result = s[::-1]  # it works!

    def test_slice_float_get_set(self):

        self.assertRaises(TypeError, lambda : self.ts[4.0:10.0])
        def f():
            self.ts[4.0:10.0] = 0
        self.assertRaises(TypeError, f)

        self.assertRaises(TypeError, self.ts.__getitem__, slice(4.5, 10.0))
        self.assertRaises(TypeError, self.ts.__setitem__, slice(4.5, 10.0), 0)

    def test_slice_floats2(self):
        s = Series(np.random.rand(10), index=np.arange(10, 20, dtype=float))

        self.assertEqual(len(s.ix[12.0:]), 8)
        self.assertEqual(len(s.ix[12.5:]), 7)

        i = np.arange(10, 20, dtype=float)
        i[2] = 12.2
        s.index = i
        self.assertEqual(len(s.ix[12.0:]), 8)
        self.assertEqual(len(s.ix[12.5:]), 7)

    def test_slice_float64(self):

        values = np.arange(10., 50., 2)
        index = Index(values)

        start, end = values[[5, 15]]

        s = Series(np.random.randn(20), index=index)

        result = s[start:end]
        expected = s.iloc[5:16]
        assert_series_equal(result, expected)

        result = s.loc[start:end]
        assert_series_equal(result, expected)

        df = DataFrame(np.random.randn(20, 3), index=index)

        result = df[start:end]
        expected = df.iloc[5:16]
        tm.assert_frame_equal(result, expected)

        result = df.loc[start:end]
        tm.assert_frame_equal(result, expected)

    def test_setitem(self):
        self.ts[self.ts.index[5]] = np.NaN
        self.ts[[1, 2, 17]] = np.NaN
        self.ts[6] = np.NaN
        self.assertTrue(np.isnan(self.ts[6]))
        self.assertTrue(np.isnan(self.ts[2]))
        self.ts[np.isnan(self.ts)] = 5
        self.assertFalse(np.isnan(self.ts[2]))

        # caught this bug when writing tests
        series = Series(tm.makeIntIndex(20).astype(float),
                        index=tm.makeIntIndex(20))

        series[::2] = 0
        self.assertTrue((series[::2] == 0).all())

        # set item that's not contained
        s = self.series.copy()
        s['foobar'] = 1

        app = Series([1], index=['foobar'], name='series')
        expected = self.series.append(app)
        assert_series_equal(s, expected)

        # Test for issue #10193
        key = pd.Timestamp('2012-01-01')
        series = pd.Series()
        series[key] = 47
        expected = pd.Series(47, [key])
        assert_series_equal(series, expected)

        series = pd.Series([], pd.DatetimeIndex([], freq='D'))
        series[key] = 47
        expected = pd.Series(47, pd.DatetimeIndex([key], freq='D'))
        assert_series_equal(series, expected)

    def test_setitem_dtypes(self):

        # change dtypes
        # GH 4463
        expected = Series([np.nan,2,3])

        s = Series([1,2,3])
        s.iloc[0] = np.nan
        assert_series_equal(s,expected)

        s = Series([1,2,3])
        s.loc[0] = np.nan
        assert_series_equal(s,expected)

        s = Series([1,2,3])
        s[0] = np.nan
        assert_series_equal(s,expected)

        s = Series([False])
        s.loc[0] = np.nan
        assert_series_equal(s,Series([np.nan]))

        s = Series([False,True])
        s.loc[0] = np.nan
        assert_series_equal(s,Series([np.nan,1.0]))

    def test_set_value(self):
        idx = self.ts.index[10]
        res = self.ts.set_value(idx, 0)
        self.assertIs(res, self.ts)
        self.assertEqual(self.ts[idx], 0)

        # equiv
        s = self.series.copy()
        res = s.set_value('foobar', 0)
        self.assertIs(res, s)
        self.assertEqual(res.index[-1], 'foobar')
        self.assertEqual(res['foobar'], 0)

        s = self.series.copy()
        s.loc['foobar'] = 0
        self.assertEqual(s.index[-1], 'foobar')
        self.assertEqual(s['foobar'], 0)

    def test_setslice(self):
        sl = self.ts[5:20]
        self.assertEqual(len(sl), len(sl.index))
        self.assertTrue(sl.index.is_unique)

    def test_basic_getitem_setitem_corner(self):
        # invalid tuples, e.g. self.ts[:, None] vs. self.ts[:, 2]
        with tm.assertRaisesRegexp(ValueError, 'tuple-index'):
            self.ts[:, 2]
        with tm.assertRaisesRegexp(ValueError, 'tuple-index'):
            self.ts[:, 2] = 2

        # weird lists. [slice(0, 5)] will work but not two slices
        result = self.ts[[slice(None, 5)]]
        expected = self.ts[:5]
        assert_series_equal(result, expected)

        # OK
        self.assertRaises(Exception, self.ts.__getitem__,
                          [5, slice(None, None)])
        self.assertRaises(Exception, self.ts.__setitem__,
                          [5, slice(None, None)], 2)

    def test_reshape_non_2d(self):
        # GH 4554
        x = Series(np.random.random(201), name='x')
        self.assertTrue(x.reshape(x.shape,) is x)

        # GH 2719
        a = Series([1, 2, 3, 4])
        result = a.reshape(2, 2)
        expected = a.values.reshape(2, 2)
        tm.assert_numpy_array_equal(result, expected)
        self.assertTrue(type(result) is type(expected))

    def test_reshape_2d_return_array(self):
        x = Series(np.random.random(201), name='x')
        result = x.reshape((-1, 1))
        self.assertNotIsInstance(result, Series)

        result2 = np.reshape(x, (-1, 1))
        self.assertNotIsInstance(result2, Series)

        result = x[:, None]
        expected = x.reshape((-1, 1))
        assert_almost_equal(result, expected)

    def test_basic_getitem_with_labels(self):
        indices = self.ts.index[[5, 10, 15]]

        result = self.ts[indices]
        expected = self.ts.reindex(indices)
        assert_series_equal(result, expected)

        result = self.ts[indices[0]:indices[2]]
        expected = self.ts.ix[indices[0]:indices[2]]
        assert_series_equal(result, expected)

        # integer indexes, be careful
        s = Series(np.random.randn(10), index=lrange(0, 20, 2))
        inds = [0, 2, 5, 7, 8]
        arr_inds = np.array([0, 2, 5, 7, 8])
        result = s[inds]
        expected = s.reindex(inds)
        assert_series_equal(result, expected)

        result = s[arr_inds]
        expected = s.reindex(arr_inds)
        assert_series_equal(result, expected)

    def test_basic_setitem_with_labels(self):
        indices = self.ts.index[[5, 10, 15]]

        cp = self.ts.copy()
        exp = self.ts.copy()
        cp[indices] = 0
        exp.ix[indices] = 0
        assert_series_equal(cp, exp)

        cp = self.ts.copy()
        exp = self.ts.copy()
        cp[indices[0]:indices[2]] = 0
        exp.ix[indices[0]:indices[2]] = 0
        assert_series_equal(cp, exp)

        # integer indexes, be careful
        s = Series(np.random.randn(10), index=lrange(0, 20, 2))
        inds = [0, 4, 6]
        arr_inds = np.array([0, 4, 6])

        cp = s.copy()
        exp = s.copy()
        s[inds] = 0
        s.ix[inds] = 0
        assert_series_equal(cp, exp)

        cp = s.copy()
        exp = s.copy()
        s[arr_inds] = 0
        s.ix[arr_inds] = 0
        assert_series_equal(cp, exp)

        inds_notfound = [0, 4, 5, 6]
        arr_inds_notfound = np.array([0, 4, 5, 6])
        self.assertRaises(Exception, s.__setitem__, inds_notfound, 0)
        self.assertRaises(Exception, s.__setitem__, arr_inds_notfound, 0)

    def test_ix_getitem(self):
        inds = self.series.index[[3, 4, 7]]
        assert_series_equal(self.series.ix[inds], self.series.reindex(inds))
        assert_series_equal(self.series.ix[5::2], self.series[5::2])

        # slice with indices
        d1, d2 = self.ts.index[[5, 15]]
        result = self.ts.ix[d1:d2]
        expected = self.ts.truncate(d1, d2)
        assert_series_equal(result, expected)

        # boolean
        mask = self.series > self.series.median()
        assert_series_equal(self.series.ix[mask], self.series[mask])

        # ask for index value
        self.assertEqual(self.ts.ix[d1], self.ts[d1])
        self.assertEqual(self.ts.ix[d2], self.ts[d2])

    def test_ix_getitem_not_monotonic(self):
        d1, d2 = self.ts.index[[5, 15]]

        ts2 = self.ts[::2][[1, 2, 0]]

        self.assertRaises(KeyError, ts2.ix.__getitem__, slice(d1, d2))
        self.assertRaises(KeyError, ts2.ix.__setitem__, slice(d1, d2), 0)

    def test_ix_getitem_setitem_integer_slice_keyerrors(self):
        s = Series(np.random.randn(10), index=lrange(0, 20, 2))

        # this is OK
        cp = s.copy()
        cp.ix[4:10] = 0
        self.assertTrue((cp.ix[4:10] == 0).all())

        # so is this
        cp = s.copy()
        cp.ix[3:11] = 0
        self.assertTrue((cp.ix[3:11] == 0).values.all())

        result = s.ix[4:10]
        result2 = s.ix[3:11]
        expected = s.reindex([4, 6, 8, 10])

        assert_series_equal(result, expected)
        assert_series_equal(result2, expected)

        # non-monotonic, raise KeyError
        s2 = s.iloc[lrange(5) + lrange(5, 10)[::-1]]
        self.assertRaises(KeyError, s2.ix.__getitem__, slice(3, 11))
        self.assertRaises(KeyError, s2.ix.__setitem__, slice(3, 11), 0)

    def test_ix_getitem_iterator(self):
        idx = iter(self.series.index[:10])
        result = self.series.ix[idx]
        assert_series_equal(result, self.series[:10])

    def test_where(self):
        s = Series(np.random.randn(5))
        cond = s > 0

        rs = s.where(cond).dropna()
        rs2 = s[cond]
        assert_series_equal(rs, rs2)

        rs = s.where(cond, -s)
        assert_series_equal(rs, s.abs())

        rs = s.where(cond)
        assert(s.shape == rs.shape)
        assert(rs is not s)

        # test alignment
        cond = Series([True,False,False,True,False],index=s.index)
        s2 = -(s.abs())

        expected = s2[cond].reindex(s2.index[:3]).reindex(s2.index)
        rs = s2.where(cond[:3])
        assert_series_equal(rs, expected)

        expected = s2.abs()
        expected.ix[0] = s2[0]
        rs = s2.where(cond[:3], -s2)
        assert_series_equal(rs, expected)

        self.assertRaises(ValueError, s.where, 1)
        self.assertRaises(ValueError, s.where, cond[:3].values, -s)

        # GH 2745
        s = Series([1, 2])
        s[[True, False]] = [0, 1]
        expected = Series([0, 2])
        assert_series_equal(s, expected)

        # failures
        self.assertRaises(
            ValueError, s.__setitem__, tuple([[[True, False]]]), [0, 2, 3])
        self.assertRaises(
            ValueError, s.__setitem__, tuple([[[True, False]]]), [])

        # unsafe dtype changes
        for dtype in [np.int8, np.int16, np.int32, np.int64, np.float16, np.float32, np.float64]:
            s = Series(np.arange(10), dtype=dtype)
            mask = s < 5
            s[mask] = lrange(2, 7)
            expected = Series(lrange(2, 7) + lrange(5, 10), dtype=dtype)
            assert_series_equal(s, expected)
            self.assertEqual(s.dtype, expected.dtype)

        # these are allowed operations, but are upcasted
        for dtype in [np.int64, np.float64]:
            s = Series(np.arange(10), dtype=dtype)
            mask = s < 5
            values = [2.5, 3.5, 4.5, 5.5, 6.5]
            s[mask] = values
            expected = Series(values + lrange(5, 10), dtype='float64')
            assert_series_equal(s, expected)
            self.assertEqual(s.dtype, expected.dtype)

        # GH 9731
        s = Series(np.arange(10), dtype='int64')
        mask = s > 5
        values = [2.5, 3.5, 4.5, 5.5]
        s[mask] = values
        expected = Series(lrange(6) + values, dtype='float64')
        assert_series_equal(s, expected)

        # can't do these as we are forced to change the itemsize of the input
        # to something we cannot
        for dtype in [np.int8, np.int16, np.int32, np.float16, np.float32]:
            s = Series(np.arange(10), dtype=dtype)
            mask = s < 5
            values = [2.5, 3.5, 4.5, 5.5, 6.5]
            self.assertRaises(Exception, s.__setitem__, tuple(mask), values)

        # GH3235
        s = Series(np.arange(10), dtype='int64')
        mask = s < 5
        s[mask] = lrange(2, 7)
        expected = Series(lrange(2, 7) + lrange(5, 10), dtype='int64')
        assert_series_equal(s, expected)
        self.assertEqual(s.dtype, expected.dtype)

        s = Series(np.arange(10), dtype='int64')
        mask = s > 5
        s[mask] = [0] * 4
        expected = Series([0, 1, 2, 3, 4, 5] + [0] * 4, dtype='int64')
        assert_series_equal(s, expected)

        s = Series(np.arange(10))
        mask = s > 5
        def f():
            s[mask] = [5,4,3,2,1]
        self.assertRaises(ValueError, f)
        def f():
            s[mask] = [0] * 5
        self.assertRaises(ValueError, f)

        # dtype changes
        s = Series([1,2,3,4])
        result = s.where(s>2,np.nan)
        expected = Series([np.nan,np.nan,3,4])
        assert_series_equal(result, expected)

        # GH 4667
        # setting with None changes dtype
        s = Series(range(10)).astype(float)
        s[8] = None
        result = s[8]
        self.assertTrue(isnull(result))

        s = Series(range(10)).astype(float)
        s[s > 8] = None
        result = s[isnull(s)]
        expected = Series(np.nan,index=[9])
        assert_series_equal(result, expected)

    def test_where_setitem_invalid(self):

        # GH 2702
        # make sure correct exceptions are raised on invalid list assignment

        # slice
        s = Series(list('abc'))
        def f():
            s[0:3] = list(range(27))
        self.assertRaises(ValueError, f)

        s[0:3] = list(range(3))
        expected = Series([0,1,2])
        assert_series_equal(s.astype(np.int64), expected, )

        # slice with step
        s = Series(list('abcdef'))
        def f():
            s[0:4:2] = list(range(27))
        self.assertRaises(ValueError, f)

        s = Series(list('abcdef'))
        s[0:4:2] = list(range(2))
        expected = Series([0,'b',1,'d','e','f'])
        assert_series_equal(s, expected)

        # neg slices
        s = Series(list('abcdef'))
        def f():
            s[:-1] = list(range(27))
        self.assertRaises(ValueError, f)

        s[-3:-1] = list(range(2))
        expected = Series(['a','b','c',0,1,'f'])
        assert_series_equal(s, expected)

        # list
        s = Series(list('abc'))
        def f():
            s[[0,1,2]] = list(range(27))
        self.assertRaises(ValueError, f)

        s = Series(list('abc'))
        def f():
            s[[0,1,2]] = list(range(2))
        self.assertRaises(ValueError, f)

        # scalar
        s = Series(list('abc'))
        s[0] = list(range(10))
        expected = Series([list(range(10)),'b','c'])
        assert_series_equal(s, expected)

    def test_where_broadcast(self):
        # Test a variety of differently sized series
        for size in range(2, 6):
            # Test a variety of boolean indices
            for selection in [np.resize([True, False, False, False, False], size),  # First element should be set
                              # Set alternating elements]
                              np.resize([True, False], size),
                              np.resize([False], size)]:  # No element should be set
                # Test a variety of different numbers as content
                for item in [2.0, np.nan, np.finfo(np.float).max, np.finfo(np.float).min]:
                    # Test numpy arrays, lists and tuples as the input to be
                    # broadcast
                    for arr in [np.array([item]), [item], (item,)]:
                        data = np.arange(size, dtype=float)
                        s = Series(data)
                        s[selection] = arr
                        # Construct the expected series by taking the source
                        # data or item based on the selection
                        expected = Series([item if use_item else data[i]
                                          for i, use_item in enumerate(selection)])
                        assert_series_equal(s, expected)

                        s = Series(data)
                        result = s.where(~selection, arr)
                        assert_series_equal(result, expected)

    def test_where_inplace(self):
        s = Series(np.random.randn(5))
        cond = s > 0

        rs = s.copy()

        rs.where(cond, inplace=True)
        assert_series_equal(rs.dropna(), s[cond])
        assert_series_equal(rs, s.where(cond))

        rs = s.copy()
        rs.where(cond, -s, inplace=True)
        assert_series_equal(rs, s.where(cond, -s))

    def test_where_dups(self):
        # GH 4550
        # where crashes with dups in index
        s1 = Series(list(range(3)))
        s2 = Series(list(range(3)))
        comb = pd.concat([s1,s2])
        result = comb.where(comb < 2)
        expected = Series([0,1,np.nan,0,1,np.nan],index=[0,1,2,0,1,2])
        assert_series_equal(result, expected)

        # GH 4548
        # inplace updating not working with dups
        comb[comb<1] = 5
        expected = Series([5,1,2,5,1,2],index=[0,1,2,0,1,2])
        assert_series_equal(comb, expected)

        comb[comb<2] += 10
        expected = Series([5,11,2,5,11,2],index=[0,1,2,0,1,2])
        assert_series_equal(comb, expected)

    def test_where_datetime(self):
        s = Series(date_range('20130102', periods=2))
        expected = Series([10, 10], dtype='datetime64[ns]')
        mask = np.array([False, False])

        rs = s.where(mask, [10, 10])
        assert_series_equal(rs, expected)

        rs = s.where(mask, 10)
        assert_series_equal(rs, expected)

        rs = s.where(mask, 10.0)
        assert_series_equal(rs, expected)

        rs = s.where(mask, [10.0, 10.0])
        assert_series_equal(rs, expected)

        rs = s.where(mask, [10.0, np.nan])
        expected = Series([10, None], dtype='datetime64[ns]')
        assert_series_equal(rs, expected)

    def test_where_timedelta(self):
        s = Series([1, 2], dtype='timedelta64[ns]')
        expected = Series([10, 10], dtype='timedelta64[ns]')
        mask = np.array([False, False])

        rs = s.where(mask, [10, 10])
        assert_series_equal(rs, expected)

        rs = s.where(mask, 10)
        assert_series_equal(rs, expected)

        rs = s.where(mask, 10.0)
        assert_series_equal(rs, expected)

        rs = s.where(mask, [10.0, 10.0])
        assert_series_equal(rs, expected)

        rs = s.where(mask, [10.0, np.nan])
        expected = Series([10, None], dtype='timedelta64[ns]')
        assert_series_equal(rs, expected)

    def test_mask(self):
        # compare with tested results in test_where
        s = Series(np.random.randn(5))
        cond = s > 0

        rs = s.where(~cond, np.nan)
        assert_series_equal(rs, s.mask(cond))

        rs = s.where(~cond)
        rs2 = s.mask(cond)
        assert_series_equal(rs, rs2)

        rs = s.where(~cond, -s)
        rs2 = s.mask(cond, -s)
        assert_series_equal(rs, rs2)

        cond = Series([True, False, False, True, False], index=s.index)
        s2 = -(s.abs())
        rs = s2.where(~cond[:3])
        rs2 = s2.mask(cond[:3])
        assert_series_equal(rs, rs2)

        rs = s2.where(~cond[:3], -s2)
        rs2 = s2.mask(cond[:3], -s2)
        assert_series_equal(rs, rs2)

        self.assertRaises(ValueError, s.mask, 1)
        self.assertRaises(ValueError, s.mask, cond[:3].values, -s)

        # dtype changes
        s = Series([1,2,3,4])
        result = s.mask(s>2, np.nan)
        expected = Series([1, 2, np.nan, np.nan])
        assert_series_equal(result, expected)

    def test_mask_broadcast(self):
        # GH 8801
        # copied from test_where_broadcast
        for size in range(2, 6):
            for selection in [np.resize([True, False, False, False, False], size),  # First element should be set
                              # Set alternating elements]
                              np.resize([True, False], size),
                              np.resize([False], size)]:  # No element should be set
                for item in [2.0, np.nan, np.finfo(np.float).max, np.finfo(np.float).min]:
                    for arr in [np.array([item]), [item], (item,)]:
                        data = np.arange(size, dtype=float)
                        s = Series(data)
                        result = s.mask(selection, arr)
                        expected = Series([item if use_item else data[i]
                                          for i, use_item in enumerate(selection)])
                        assert_series_equal(result, expected)

    def test_mask_inplace(self):
        s = Series(np.random.randn(5))
        cond = s > 0

        rs = s.copy()
        rs.mask(cond, inplace=True)
        assert_series_equal(rs.dropna(), s[~cond])
        assert_series_equal(rs, s.mask(cond))

        rs = s.copy()
        rs.mask(cond, -s, inplace=True)
        assert_series_equal(rs, s.mask(cond, -s))

    def test_drop(self):

        # unique
        s = Series([1,2],index=['one','two'])
        expected = Series([1],index=['one'])
        result = s.drop(['two'])
        assert_series_equal(result,expected)
        result = s.drop('two', axis='rows')
        assert_series_equal(result,expected)

        # non-unique
        # GH 5248
        s = Series([1,1,2],index=['one','two','one'])
        expected = Series([1,2],index=['one','one'])
        result = s.drop(['two'], axis=0)
        assert_series_equal(result,expected)
        result = s.drop('two')
        assert_series_equal(result,expected)

        expected = Series([1],index=['two'])
        result = s.drop(['one'])
        assert_series_equal(result,expected)
        result = s.drop('one')
        assert_series_equal(result,expected)

        # single string/tuple-like
        s = Series(range(3),index=list('abc'))
        self.assertRaises(ValueError, s.drop, 'bc')
        self.assertRaises(ValueError, s.drop, ('a',))

        # errors='ignore'
        s = Series(range(3),index=list('abc'))
        result = s.drop('bc', errors='ignore')
        assert_series_equal(result, s)
        result = s.drop(['a', 'd'], errors='ignore')
        expected = s.ix[1:]
        assert_series_equal(result, expected)

        # bad axis
        self.assertRaises(ValueError, s.drop, 'one', axis='columns')

        # GH 8522
        s = Series([2,3], index=[True, False])
        self.assertTrue(s.index.is_object())
        result = s.drop(True)
        expected = Series([3],index=[False])
        assert_series_equal(result,expected)

    def test_ix_setitem(self):
        inds = self.series.index[[3, 4, 7]]

        result = self.series.copy()
        result.ix[inds] = 5

        expected = self.series.copy()
        expected[[3, 4, 7]] = 5
        assert_series_equal(result, expected)

        result.ix[5:10] = 10
        expected[5:10] = 10
        assert_series_equal(result, expected)

        # set slice with indices
        d1, d2 = self.series.index[[5, 15]]
        result.ix[d1:d2] = 6
        expected[5:16] = 6  # because it's inclusive
        assert_series_equal(result, expected)

        # set index value
        self.series.ix[d1] = 4
        self.series.ix[d2] = 6
        self.assertEqual(self.series[d1], 4)
        self.assertEqual(self.series[d2], 6)

    def test_where_numeric_with_string(self):
        # GH 9280
        s = pd.Series([1, 2, 3])
        w = s.where(s>1, 'X')

        self.assertFalse(com.is_integer(w[0]))
        self.assertTrue(com.is_integer(w[1]))
        self.assertTrue(com.is_integer(w[2]))
        self.assertTrue(isinstance(w[0], str))
        self.assertTrue(w.dtype == 'object')

        w = s.where(s>1, ['X', 'Y', 'Z'])
        self.assertFalse(com.is_integer(w[0]))
        self.assertTrue(com.is_integer(w[1]))
        self.assertTrue(com.is_integer(w[2]))
        self.assertTrue(isinstance(w[0], str))
        self.assertTrue(w.dtype == 'object')

        w = s.where(s>1, np.array(['X', 'Y', 'Z']))
        self.assertFalse(com.is_integer(w[0]))
        self.assertTrue(com.is_integer(w[1]))
        self.assertTrue(com.is_integer(w[2]))
        self.assertTrue(isinstance(w[0], str))
        self.assertTrue(w.dtype == 'object')

    def test_setitem_boolean(self):
        mask = self.series > self.series.median()

        # similiar indexed series
        result = self.series.copy()
        result[mask] = self.series * 2
        expected = self.series * 2
        assert_series_equal(result[mask], expected[mask])

        # needs alignment
        result = self.series.copy()
        result[mask] = (self.series * 2)[0:5]
        expected = (self.series * 2)[0:5].reindex_like(self.series)
        expected[-mask] = self.series[mask]
        assert_series_equal(result[mask], expected[mask])

    def test_ix_setitem_boolean(self):
        mask = self.series > self.series.median()

        result = self.series.copy()
        result.ix[mask] = 0
        expected = self.series
        expected[mask] = 0
        assert_series_equal(result, expected)

    def test_ix_setitem_corner(self):
        inds = list(self.series.index[[5, 8, 12]])
        self.series.ix[inds] = 5
        self.assertRaises(Exception, self.series.ix.__setitem__,
                          inds + ['foo'], 5)

    def test_get_set_boolean_different_order(self):
        ordered = self.series.sort_values()

        # setting
        copy = self.series.copy()
        copy[ordered > 0] = 0

        expected = self.series.copy()
        expected[expected > 0] = 0

        assert_series_equal(copy, expected)

        # getting
        sel = self.series[ordered > 0]
        exp = self.series[self.series > 0]
        assert_series_equal(sel, exp)

    def test_repr(self):
        str(self.ts)
        str(self.series)
        str(self.series.astype(int))
        str(self.objSeries)

        str(Series(tm.randn(1000), index=np.arange(1000)))
        str(Series(tm.randn(1000), index=np.arange(1000, 0, step=-1)))

        # empty
        str(self.empty)

        # with NaNs
        self.series[5:7] = np.NaN
        str(self.series)

        # with Nones
        ots = self.ts.astype('O')
        ots[::2] = None
        repr(ots)

        # various names
        for name in ['', 1, 1.2, 'foo', u('\u03B1\u03B2\u03B3'),
                     'loooooooooooooooooooooooooooooooooooooooooooooooooooong',
                     ('foo', 'bar', 'baz'),
                     (1, 2),
                     ('foo', 1, 2.3),
                     (u('\u03B1'), u('\u03B2'), u('\u03B3')),
                     (u('\u03B1'), 'bar')]:
            self.series.name = name
            repr(self.series)

        biggie = Series(tm.randn(1000), index=np.arange(1000),
                        name=('foo', 'bar', 'baz'))
        repr(biggie)

        # 0 as name
        ser = Series(np.random.randn(100), name=0)
        rep_str = repr(ser)
        self.assertIn("Name: 0", rep_str)

        # tidy repr
        ser = Series(np.random.randn(1001), name=0)
        rep_str = repr(ser)
        self.assertIn("Name: 0", rep_str)

        ser = Series(["a\n\r\tb"], name=["a\n\r\td"], index=["a\n\r\tf"])
        self.assertFalse("\t" in repr(ser))
        self.assertFalse("\r" in repr(ser))
        self.assertFalse("a\n" in repr(ser))

        # with empty series (#4651)
        s = Series([], dtype=np.int64, name='foo')
        self.assertEqual(repr(s), 'Series([], Name: foo, dtype: int64)')

        s = Series([], dtype=np.int64, name=None)
        self.assertEqual(repr(s), 'Series([], dtype: int64)')

    def test_tidy_repr(self):
        a = Series([u("\u05d0")] * 1000)
        a.name = 'title1'
        repr(a)         # should not raise exception

    def test_repr_bool_fails(self):
        s = Series([DataFrame(np.random.randn(2, 2)) for i in range(5)])

        import sys

        buf = StringIO()
        tmp = sys.stderr
        sys.stderr = buf
        try:
        # it works (with no Cython exception barf)!
            repr(s)
        finally:
            sys.stderr = tmp
        self.assertEqual(buf.getvalue(), '')

    def test_repr_name_iterable_indexable(self):
        s = Series([1, 2, 3], name=np.int64(3))

        # it works!
        repr(s)

        s.name = (u("\u05d0"),) * 2
        repr(s)

    def test_repr_should_return_str(self):
        # http://docs.python.org/py3k/reference/datamodel.html#object.__repr__
        # http://docs.python.org/reference/datamodel.html#object.__repr__
        # ...The return value must be a string object.

        # (str on py2.x, str (unicode) on py3)

        data = [8, 5, 3, 5]
        index1 = [u("\u03c3"), u("\u03c4"), u("\u03c5"), u("\u03c6")]
        df = Series(data, index=index1)
        self.assertTrue(type(df.__repr__() == str))  # both py2 / 3

    def test_repr_max_rows(self):
        # GH 6863
        with pd.option_context('max_rows', None):
            str(Series(range(1001))) # should not raise exception

    def test_unicode_string_with_unicode(self):
        df = Series([u("\u05d0")], name=u("\u05d1"))
        if compat.PY3:
            str(df)
        else:
            compat.text_type(df)

    def test_bytestring_with_unicode(self):
        df = Series([u("\u05d0")], name=u("\u05d1"))
        if compat.PY3:
            bytes(df)
        else:
            str(df)

    def test_timeseries_repr_object_dtype(self):
        index = Index([datetime(2000, 1, 1) + timedelta(i)
                       for i in range(1000)], dtype=object)
        ts = Series(np.random.randn(len(index)), index)
        repr(ts)

        ts = tm.makeTimeSeries(1000)
        self.assertTrue(repr(ts).splitlines()[-1].startswith('Freq:'))

        ts2 = ts.ix[np.random.randint(0, len(ts) - 1, 400)]
        repr(ts2).splitlines()[-1]

    def test_timeseries_periodindex(self):
        # GH2891
        from pandas import period_range
        prng = period_range('1/1/2011', '1/1/2012', freq='M')
        ts = Series(np.random.randn(len(prng)), prng)
        new_ts = self.round_trip_pickle(ts)
        self.assertEqual(new_ts.index.freq, 'M')

    def test_iter(self):
        for i, val in enumerate(self.series):
            self.assertEqual(val, self.series[i])

        for i, val in enumerate(self.ts):
            self.assertEqual(val, self.ts[i])

    def test_keys(self):
        # HACK: By doing this in two stages, we avoid 2to3 wrapping the call
        # to .keys() in a list()
        getkeys = self.ts.keys
        self.assertIs(getkeys(), self.ts.index)

    def test_values(self):
        self.assert_numpy_array_equal(self.ts, self.ts.values)

    def test_iteritems(self):
        for idx, val in compat.iteritems(self.series):
            self.assertEqual(val, self.series[idx])

        for idx, val in compat.iteritems(self.ts):
            self.assertEqual(val, self.ts[idx])

        # assert is lazy (genrators don't define reverse, lists do)
        self.assertFalse(hasattr(self.series.iteritems(), 'reverse'))

    def test_sum(self):
        self._check_stat_op('sum', np.sum, check_allna=True)

    def test_sum_inf(self):
        import pandas.core.nanops as nanops

        s = Series(np.random.randn(10))
        s2 = s.copy()

        s[5:8] = np.inf
        s2[5:8] = np.nan

        self.assertTrue(np.isinf(s.sum()))

        arr = np.random.randn(100, 100).astype('f4')
        arr[:, 2] = np.inf

        with cf.option_context("mode.use_inf_as_null", True):
            assert_almost_equal(s.sum(), s2.sum())

        res = nanops.nansum(arr, axis=1)
        self.assertTrue(np.isinf(res).all())

    def test_mean(self):
        self._check_stat_op('mean', np.mean)

    def test_median(self):
        self._check_stat_op('median', np.median)

        # test with integers, test failure
        int_ts = Series(np.ones(10, dtype=int), index=lrange(10))
        self.assertAlmostEqual(np.median(int_ts), int_ts.median())

    def test_mode(self):
        s = Series([12, 12, 11, 10, 19, 11])
        exp = Series([11, 12])
        assert_series_equal(s.mode(), exp)

        assert_series_equal(Series([1, 2, 3]).mode(), Series([], dtype='int64'))

        lst = [5] * 20 + [1] * 10 + [6] * 25
        np.random.shuffle(lst)
        s = Series(lst)
        assert_series_equal(s.mode(), Series([6]))

        s = Series([5] * 10)
        assert_series_equal(s.mode(), Series([5]))

        s = Series(lst)
        s[0] = np.nan
        assert_series_equal(s.mode(), Series([6.]))

        s = Series(list('adfasbasfwewefwefweeeeasdfasnbam'))
        assert_series_equal(s.mode(), Series(['e']))

        s = Series(['2011-01-03', '2013-01-02', '1900-05-03'], dtype='M8[ns]')
        assert_series_equal(s.mode(), Series([], dtype="M8[ns]"))
        s = Series(['2011-01-03', '2013-01-02', '1900-05-03', '2011-01-03',
                    '2013-01-02'], dtype='M8[ns]')
        assert_series_equal(s.mode(), Series(['2011-01-03', '2013-01-02'],
                                             dtype='M8[ns]'))

    def test_prod(self):
        self._check_stat_op('prod', np.prod)

    def test_min(self):
        self._check_stat_op('min', np.min, check_objects=True)

    def test_max(self):
        self._check_stat_op('max', np.max, check_objects=True)

    def test_var_std(self):
        alt = lambda x: np.std(x, ddof=1)
        self._check_stat_op('std', alt)

        alt = lambda x: np.var(x, ddof=1)
        self._check_stat_op('var', alt)

        result = self.ts.std(ddof=4)
        expected = np.std(self.ts.values, ddof=4)
        assert_almost_equal(result, expected)

        result = self.ts.var(ddof=4)
        expected = np.var(self.ts.values, ddof=4)
        assert_almost_equal(result, expected)

        # 1 - element series with ddof=1
        s = self.ts.iloc[[0]]
        result = s.var(ddof=1)
        self.assertTrue(isnull(result))

        result = s.std(ddof=1)
        self.assertTrue(isnull(result))

    def test_sem(self):
        alt = lambda x: np.std(x, ddof=1)/np.sqrt(len(x))
        self._check_stat_op('sem', alt)

        result = self.ts.sem(ddof=4)
        expected = np.std(self.ts.values, ddof=4)/np.sqrt(len(self.ts.values))
        assert_almost_equal(result, expected)

        # 1 - element series with ddof=1
        s = self.ts.iloc[[0]]
        result = s.sem(ddof=1)
        self.assertTrue(isnull(result))

    def test_skew(self):
        tm._skip_if_no_scipy()

        from scipy.stats import skew
        alt = lambda x: skew(x, bias=False)
        self._check_stat_op('skew', alt)

        # test corner cases, skew() returns NaN unless there's at least 3 values
        min_N = 3
        for i in range(1, min_N + 1):
            s = Series(np.ones(i))
            df = DataFrame(np.ones((i, i)))
            if i < min_N:
                self.assertTrue(np.isnan(s.skew()))
                self.assertTrue(np.isnan(df.skew()).all())
            else:
                self.assertEqual(0, s.skew())
                self.assertTrue((df.skew() == 0).all())

    def test_kurt(self):
        tm._skip_if_no_scipy()

        from scipy.stats import kurtosis
        alt = lambda x: kurtosis(x, bias=False)
        self._check_stat_op('kurt', alt)

        index = MultiIndex(levels=[['bar'], ['one', 'two', 'three'], [0, 1]],
                           labels=[[0, 0, 0, 0, 0, 0],
                                   [0, 1, 2, 0, 1, 2],
                                   [0, 1, 0, 1, 0, 1]])
        s = Series(np.random.randn(6), index=index)
        self.assertAlmostEqual(s.kurt(), s.kurt(level=0)['bar'])

        # test corner cases, kurt() returns NaN unless there's at least 4 values
        min_N = 4
        for i in range(1, min_N + 1):
            s = Series(np.ones(i))
            df = DataFrame(np.ones((i, i)))
            if i < min_N:
                self.assertTrue(np.isnan(s.kurt()))
                self.assertTrue(np.isnan(df.kurt()).all())
            else:
                self.assertEqual(0, s.kurt())
                self.assertTrue((df.kurt() == 0).all())

    def test_argsort(self):
        self._check_accum_op('argsort')
        argsorted = self.ts.argsort()
        self.assertTrue(issubclass(argsorted.dtype.type, np.integer))

        # GH 2967 (introduced bug in 0.11-dev I think)
        s = Series([Timestamp('201301%02d' % (i + 1)) for i in range(5)])
        self.assertEqual(s.dtype, 'datetime64[ns]')
        shifted = s.shift(-1)
        self.assertEqual(shifted.dtype, 'datetime64[ns]')
        self.assertTrue(isnull(shifted[4]))

        result = s.argsort()
        expected = Series(lrange(5), dtype='int64')
        assert_series_equal(result, expected)

        result = shifted.argsort()
        expected = Series(lrange(4) + [-1], dtype='int64')
        assert_series_equal(result, expected)

    def test_argsort_stable(self):
        s = Series(np.random.randint(0, 100, size=10000))
        mindexer = s.argsort(kind='mergesort')
        qindexer = s.argsort()

        mexpected = np.argsort(s.values, kind='mergesort')
        qexpected = np.argsort(s.values, kind='quicksort')

        self.assert_numpy_array_equal(mindexer, mexpected)
        self.assert_numpy_array_equal(qindexer, qexpected)
        self.assertFalse(np.array_equal(qindexer, mindexer))

    def test_reorder_levels(self):
        index = MultiIndex(levels=[['bar'], ['one', 'two', 'three'], [0, 1]],
                           labels=[[0, 0, 0, 0, 0, 0],
                                   [0, 1, 2, 0, 1, 2],
                                   [0, 1, 0, 1, 0, 1]],
                           names=['L0', 'L1', 'L2'])
        s = Series(np.arange(6), index=index)

        # no change, position
        result = s.reorder_levels([0, 1, 2])
        assert_series_equal(s, result)

        # no change, labels
        result = s.reorder_levels(['L0', 'L1', 'L2'])
        assert_series_equal(s, result)

        # rotate, position
        result = s.reorder_levels([1, 2, 0])
        e_idx = MultiIndex(levels=[['one', 'two', 'three'], [0, 1], ['bar']],
                           labels=[[0, 1, 2, 0, 1, 2],
                                   [0, 1, 0, 1, 0, 1],
                                   [0, 0, 0, 0, 0, 0]],
                           names=['L1', 'L2', 'L0'])
        expected = Series(np.arange(6), index=e_idx)
        assert_series_equal(result, expected)

        result = s.reorder_levels([0, 0, 0])
        e_idx = MultiIndex(levels=[['bar'], ['bar'], ['bar']],
                           labels=[[0, 0, 0, 0, 0, 0],
                                   [0, 0, 0, 0, 0, 0],
                                   [0, 0, 0, 0, 0, 0]],
                           names=['L0', 'L0', 'L0'])
        expected = Series(range(6), index=e_idx)
        assert_series_equal(result, expected)

        result = s.reorder_levels(['L0', 'L0', 'L0'])
        assert_series_equal(result, expected)

    def test_cumsum(self):
        self._check_accum_op('cumsum')

    def test_cumprod(self):
        self._check_accum_op('cumprod')

    def test_cummin(self):
        self.assert_numpy_array_equal(self.ts.cummin(),
                                      np.minimum.accumulate(np.array(self.ts)))
        ts = self.ts.copy()
        ts[::2] = np.NaN
        result = ts.cummin()[1::2]
        expected = np.minimum.accumulate(ts.valid())

        self.assert_numpy_array_equal(result, expected)

    def test_cummax(self):
        self.assert_numpy_array_equal(self.ts.cummax(),
                                      np.maximum.accumulate(np.array(self.ts)))
        ts = self.ts.copy()
        ts[::2] = np.NaN
        result = ts.cummax()[1::2]
        expected = np.maximum.accumulate(ts.valid())

        self.assert_numpy_array_equal(result, expected)

    def test_cummin_datetime64(self):
        s = pd.Series(pd.to_datetime(
            ['NaT', '2000-1-2', 'NaT', '2000-1-1', 'NaT', '2000-1-3']))

        expected = pd.Series(pd.to_datetime(
            ['NaT', '2000-1-2', 'NaT', '2000-1-1', 'NaT', '2000-1-1']))
        result = s.cummin(skipna=True)
        self.assert_series_equal(expected, result)

        expected = pd.Series(pd.to_datetime(
            ['NaT', '2000-1-2', '2000-1-2', '2000-1-1', '2000-1-1', '2000-1-1']))
        result = s.cummin(skipna=False)
        self.assert_series_equal(expected, result)

    def test_cummax_datetime64(self):
        s = pd.Series(pd.to_datetime(
            ['NaT', '2000-1-2', 'NaT', '2000-1-1', 'NaT', '2000-1-3']))

        expected = pd.Series(pd.to_datetime(
            ['NaT', '2000-1-2', 'NaT', '2000-1-2', 'NaT', '2000-1-3']))
        result = s.cummax(skipna=True)
        self.assert_series_equal(expected, result)

        expected = pd.Series(pd.to_datetime(
            ['NaT', '2000-1-2', '2000-1-2', '2000-1-2', '2000-1-2', '2000-1-3']))
        result = s.cummax(skipna=False)
        self.assert_series_equal(expected, result)

    def test_cummin_timedelta64(self):
        s = pd.Series(pd.to_timedelta(
            ['NaT', '2 min', 'NaT', '1 min', 'NaT', '3 min', ]))

        expected = pd.Series(pd.to_timedelta(
            ['NaT', '2 min', 'NaT', '1 min', 'NaT', '1 min', ]))
        result = s.cummin(skipna=True)
        self.assert_series_equal(expected, result)

        expected = pd.Series(pd.to_timedelta(
            ['NaT', '2 min', '2 min', '1 min', '1 min', '1 min', ]))
        result = s.cummin(skipna=False)
        self.assert_series_equal(expected, result)

    def test_cummax_timedelta64(self):
        s = pd.Series(pd.to_timedelta(
            ['NaT', '2 min', 'NaT', '1 min', 'NaT', '3 min', ]))

        expected = pd.Series(pd.to_timedelta(
            ['NaT', '2 min', 'NaT', '2 min', 'NaT', '3 min', ]))
        result = s.cummax(skipna=True)
        self.assert_series_equal(expected, result)

        expected = pd.Series(pd.to_timedelta(
            ['NaT', '2 min', '2 min', '2 min', '2 min', '3 min', ]))
        result = s.cummax(skipna=False)
        self.assert_series_equal(expected, result)

    def test_npdiff(self):
        raise nose.SkipTest("skipping due to Series no longer being an "
                            "ndarray")

        # no longer works as the return type of np.diff is now nd.array
        s = Series(np.arange(5))

        r = np.diff(s)
        assert_series_equal(Series([nan, 0, 0, 0, nan]), r)

    def _check_stat_op(self, name, alternate, check_objects=False, check_allna=False):
        import pandas.core.nanops as nanops

        def testit():
            f = getattr(Series, name)

            # add some NaNs
            self.series[5:15] = np.NaN

            # idxmax, idxmin, min, and max are valid for dates
            if name not in ['max','min']:
                ds = Series(date_range('1/1/2001', periods=10))
                self.assertRaises(TypeError, f, ds)

            # skipna or no
            self.assertTrue(notnull(f(self.series)))
            self.assertTrue(isnull(f(self.series, skipna=False)))

            # check the result is correct
            nona = self.series.dropna()
            assert_almost_equal(f(nona), alternate(nona.values))
            assert_almost_equal(f(self.series), alternate(nona.values))

            allna = self.series * nan

            if check_allna:
                # xref 9422
                # bottleneck >= 1.0 give 0.0 for an allna Series sum
                try:
                    self.assertTrue(nanops._USE_BOTTLENECK)
                    import bottleneck as bn
                    self.assertTrue(bn.__version__ >= LooseVersion('1.0'))
                    self.assertEqual(f(allna),0.0)
                except:
                    self.assertTrue(np.isnan(f(allna)))

            # dtype=object with None, it works!
            s = Series([1, 2, 3, None, 5])
            f(s)

            # 2888
            l = [0]
            l.extend(lrange(2 ** 40, 2 ** 40+1000))
            s = Series(l, dtype='int64')
            assert_almost_equal(float(f(s)), float(alternate(s.values)))

            # check date range
            if check_objects:
                s = Series(bdate_range('1/1/2000', periods=10))
                res = f(s)
                exp = alternate(s)
                self.assertEqual(res, exp)

            # check on string data
            if name not in ['sum','min','max']:
                self.assertRaises(TypeError, f, Series(list('abc')))

            # Invalid axis.
            self.assertRaises(ValueError, f, self.series, axis=1)

            # Unimplemented numeric_only parameter.
            if 'numeric_only' in getargspec(f).args:
                self.assertRaisesRegexp(NotImplementedError, name, f,
                                        self.series, numeric_only=True)

        testit()

        try:
            import bottleneck as bn
            nanops._USE_BOTTLENECK = False
            testit()
            nanops._USE_BOTTLENECK = True
        except ImportError:
            pass

    def _check_accum_op(self, name):
        func = getattr(np, name)
        self.assert_numpy_array_equal(func(self.ts), func(np.array(self.ts)))

        # with missing values
        ts = self.ts.copy()
        ts[::2] = np.NaN

        result = func(ts)[1::2]
        expected = func(np.array(ts.valid()))

        self.assert_numpy_array_equal(result, expected)

    def test_round(self):
        # numpy.round doesn't preserve metadata, probably a numpy bug,
        # re: GH #314
        result = np.round(self.ts, 2)
        expected = Series(np.round(self.ts.values, 2), index=self.ts.index,
                          name='ts')
        assert_series_equal(result, expected)
        self.assertEqual(result.name, self.ts.name)

    def test_prod_numpy16_bug(self):
        s = Series([1., 1., 1.], index=lrange(3))
        result = s.prod()
        self.assertNotIsInstance(result, Series)

    def test_quantile(self):
        from numpy import percentile

        q = self.ts.quantile(0.1)
        self.assertEqual(q, percentile(self.ts.valid(), 10))

        q = self.ts.quantile(0.9)
        self.assertEqual(q, percentile(self.ts.valid(), 90))

        # object dtype
        q = Series(self.ts,dtype=object).quantile(0.9)
        self.assertEqual(q, percentile(self.ts.valid(), 90))

        # datetime64[ns] dtype
        dts = self.ts.index.to_series()
        q = dts.quantile(.2)
        self.assertEqual(q, Timestamp('2000-01-10 19:12:00'))

        # timedelta64[ns] dtype
        tds = dts.diff()
        q = tds.quantile(.25)
        self.assertEqual(q, pd.to_timedelta('24:00:00'))

        # GH7661
        result = Series([np.timedelta64('NaT')]).sum()
        self.assertTrue(result is pd.NaT)

        msg = 'percentiles should all be in the interval \\[0, 1\\]'
        for invalid in [-1, 2, [0.5, -1], [0.5, 2]]:
            with tm.assertRaisesRegexp(ValueError, msg):
                self.ts.quantile(invalid)

    def test_quantile_multi(self):
        from numpy import percentile

        qs = [.1, .9]
        result = self.ts.quantile(qs)
        expected = pd.Series([percentile(self.ts.valid(), 10),
                              percentile(self.ts.valid(), 90)],
                             index=qs, name=self.ts.name)
        assert_series_equal(result, expected)

        dts = self.ts.index.to_series()
        dts.name = 'xxx'
        result = dts.quantile((.2, .2))
        expected = Series([Timestamp('2000-01-10 19:12:00'),
                           Timestamp('2000-01-10 19:12:00')],
                          index=[.2, .2], name='xxx')
        assert_series_equal(result, expected)

        result = self.ts.quantile([])
        expected = pd.Series([], name=self.ts.name, index=Index([], dtype=float))
        assert_series_equal(result, expected)

    def test_append(self):
        appendedSeries = self.series.append(self.objSeries)
        for idx, value in compat.iteritems(appendedSeries):
            if idx in self.series.index:
                self.assertEqual(value, self.series[idx])
            elif idx in self.objSeries.index:
                self.assertEqual(value, self.objSeries[idx])
            else:
                self.fail("orphaned index!")

        self.assertRaises(ValueError, self.ts.append, self.ts,
                          verify_integrity=True)

    def test_append_many(self):
        pieces = [self.ts[:5], self.ts[5:10], self.ts[10:]]

        result = pieces[0].append(pieces[1:])
        assert_series_equal(result, self.ts)

    def test_all_any(self):
        ts = tm.makeTimeSeries()
        bool_series = ts > 0
        self.assertFalse(bool_series.all())
        self.assertTrue(bool_series.any())

        # Alternative types, with implicit 'object' dtype.
        s = Series(['abc', True])
        self.assertEqual('abc', s.any())  # 'abc' || True => 'abc'

    def test_all_any_params(self):
        # Check skipna, with implicit 'object' dtype.
        s1 = Series([np.nan, True])
        s2 = Series([np.nan, False])
        self.assertTrue(s1.all(skipna=False))  # nan && True => True
        self.assertTrue(s1.all(skipna=True))
        self.assertTrue(np.isnan(s2.any(skipna=False)))  # nan || False => nan
        self.assertFalse(s2.any(skipna=True))

        # Check level.
        s = pd.Series([False, False, True, True, False, True],
                      index=[0, 0, 1, 1, 2, 2])
        assert_series_equal(s.all(level=0), Series([False, True, False]))
        assert_series_equal(s.any(level=0), Series([False, True, True]))

        # bool_only is not implemented with level option.
        self.assertRaises(NotImplementedError, s.any, bool_only=True, level=0)
        self.assertRaises(NotImplementedError, s.all, bool_only=True, level=0)

        # bool_only is not implemented alone.
        self.assertRaises(NotImplementedError, s.any, bool_only=True)
        self.assertRaises(NotImplementedError, s.all, bool_only=True)

    def test_op_method(self):
        def check(series, other, check_reverse=False):
            simple_ops = ['add', 'sub', 'mul', 'floordiv', 'truediv', 'pow']
            if not compat.PY3:
                simple_ops.append('div')

            for opname in simple_ops:
                op = getattr(Series, opname)

                if op == 'div':
                    alt = operator.truediv
                else:
                    alt = getattr(operator, opname)

                result = op(series, other)
                expected = alt(series, other)
                tm.assert_almost_equal(result, expected)
                if check_reverse:
                    rop = getattr(Series, "r" + opname)
                    result = rop(series, other)
                    expected = alt(other, series)
                    tm.assert_almost_equal(result, expected)

        check(self.ts, self.ts * 2)
        check(self.ts, self.ts[::2])
        check(self.ts, 5, check_reverse=True)
        check(tm.makeFloatSeries(), tm.makeFloatSeries(), check_reverse=True)

    def test_neg(self):
        assert_series_equal(-self.series, -1 * self.series)

    def test_invert(self):
        assert_series_equal(-(self.series < 0), ~(self.series < 0))

    def test_modulo(self):

        # GH3590, modulo as ints
        p = DataFrame({'first': [3, 4, 5, 8], 'second': [0, 0, 0, 3]})
        result = p['first'] % p['second']
        expected = Series(p['first'].values %
                          p['second'].values, dtype='float64')
        expected.iloc[0:3] = np.nan
        assert_series_equal(result, expected)

        result = p['first'] % 0
        expected = Series(np.nan, index=p.index, name='first')
        assert_series_equal(result, expected)

        p = p.astype('float64')
        result = p['first'] % p['second']
        expected = Series(p['first'].values % p['second'].values)
        assert_series_equal(result, expected)

        p = p.astype('float64')
        result = p['first'] % p['second']
        result2 = p['second'] % p['first']
        self.assertFalse(np.array_equal(result, result2))

        # GH 9144
        s = Series([0, 1])

        result = s % 0
        expected = Series([nan, nan])
        assert_series_equal(result, expected)

        result = 0 % s
        expected = Series([nan, 0.0])
        assert_series_equal(result, expected)

    def test_div(self):

        # no longer do integer div for any ops, but deal with the 0's
        p = DataFrame({'first': [3, 4, 5, 8], 'second': [0, 0, 0, 3]})
        result = p['first'] / p['second']
        expected = Series(p['first'].values.astype(float) / p['second'].values,
                          dtype='float64')
        expected.iloc[0:3] = np.inf
        assert_series_equal(result, expected)

        result = p['first'] / 0
        expected = Series(np.inf, index=p.index, name='first')
        assert_series_equal(result, expected)

        p = p.astype('float64')
        result = p['first'] / p['second']
        expected = Series(p['first'].values / p['second'].values)
        assert_series_equal(result, expected)

        p = DataFrame({'first': [3, 4, 5, 8], 'second': [1, 1, 1, 1]})
        result = p['first'] / p['second']
        assert_series_equal(result, p['first'].astype('float64'), check_names=False)
        self.assertTrue(result.name is None)
        self.assertFalse(np.array_equal(result, p['second'] / p['first']))

        # inf signing
        s = Series([np.nan,1.,-1.])
        result = s / 0
        expected = Series([np.nan,np.inf,-np.inf])
        assert_series_equal(result, expected)

        # float/integer issue
        # GH 7785
        p = DataFrame({'first': (1,0), 'second': (-0.01,-0.02)})
        expected = Series([-0.01,-np.inf])

        result = p['second'].div(p['first'])
        assert_series_equal(result, expected, check_names=False)

        result = p['second'] / p['first']
        assert_series_equal(result, expected)

        # GH 9144
        s = Series([-1, 0, 1])

        result = 0 / s
        expected = Series([0.0, nan, 0.0])
        assert_series_equal(result, expected)

        result = s / 0
        expected = Series([-inf, nan, inf])
        assert_series_equal(result, expected)

        result = s // 0
        expected = Series([-inf, nan, inf])
        assert_series_equal(result, expected)

    def test_operators(self):

        def _check_op(series, other, op, pos_only=False):
            left = np.abs(series) if pos_only else series
            right = np.abs(other) if pos_only else other

            cython_or_numpy = op(left, right)
            python = left.combine(right, op)
            tm.assert_almost_equal(cython_or_numpy, python)

        def check(series, other):
            simple_ops = ['add', 'sub', 'mul', 'truediv', 'floordiv', 'mod']

            for opname in simple_ops:
                _check_op(series, other, getattr(operator, opname))

            _check_op(series, other, operator.pow, pos_only=True)

            _check_op(series, other, lambda x, y: operator.add(y, x))
            _check_op(series, other, lambda x, y: operator.sub(y, x))
            _check_op(series, other, lambda x, y: operator.truediv(y, x))
            _check_op(series, other, lambda x, y: operator.floordiv(y, x))
            _check_op(series, other, lambda x, y: operator.mul(y, x))
            _check_op(series, other, lambda x, y: operator.pow(y, x),
                      pos_only=True)
            _check_op(series, other, lambda x, y: operator.mod(y, x))

        check(self.ts, self.ts * 2)
        check(self.ts, self.ts * 0)
        check(self.ts, self.ts[::2])
        check(self.ts, 5)

        def check_comparators(series, other):
            _check_op(series, other, operator.gt)
            _check_op(series, other, operator.ge)
            _check_op(series, other, operator.eq)
            _check_op(series, other, operator.lt)
            _check_op(series, other, operator.le)

        check_comparators(self.ts, 5)
        check_comparators(self.ts, self.ts + 1)

    def test_operators_empty_int_corner(self):
        s1 = Series([], [], dtype=np.int32)
        s2 = Series({'x': 0.})
        tm.assert_series_equal(s1 * s2, Series([np.nan], index=['x']))

    def test_constructor_dtype_timedelta64(self):

        # basic
        td = Series([timedelta(days=i) for i in range(3)])
        self.assertEqual(td.dtype, 'timedelta64[ns]')

        td = Series([timedelta(days=1)])
        self.assertEqual(td.dtype, 'timedelta64[ns]')

        td = Series([timedelta(days=1),timedelta(days=2),np.timedelta64(1,'s')])
        self.assertEqual(td.dtype, 'timedelta64[ns]')

        # mixed with NaT
        from pandas import tslib
        td = Series([timedelta(days=1),tslib.NaT ], dtype='m8[ns]' )
        self.assertEqual(td.dtype, 'timedelta64[ns]')

        td = Series([timedelta(days=1),np.nan ], dtype='m8[ns]' )
        self.assertEqual(td.dtype, 'timedelta64[ns]')

        td = Series([np.timedelta64(300000000), pd.NaT],dtype='m8[ns]')
        self.assertEqual(td.dtype, 'timedelta64[ns]')

        # improved inference
        # GH5689
        td = Series([np.timedelta64(300000000), pd.NaT])
        self.assertEqual(td.dtype, 'timedelta64[ns]')

        td = Series([np.timedelta64(300000000), tslib.iNaT])
        self.assertEqual(td.dtype, 'timedelta64[ns]')

        td = Series([np.timedelta64(300000000), np.nan])
        self.assertEqual(td.dtype, 'timedelta64[ns]')

        td = Series([pd.NaT, np.timedelta64(300000000)])
        self.assertEqual(td.dtype, 'timedelta64[ns]')

        td = Series([np.timedelta64(1,'s')])
        self.assertEqual(td.dtype, 'timedelta64[ns]')

        # these are frequency conversion astypes
        #for t in ['s', 'D', 'us', 'ms']:
        #    self.assertRaises(TypeError, td.astype, 'm8[%s]' % t)

        # valid astype
        td.astype('int64')

        # invalid casting
        self.assertRaises(TypeError, td.astype, 'int32')

        # this is an invalid casting
        def f():
            Series([timedelta(days=1), 'foo'],dtype='m8[ns]')
        self.assertRaises(Exception, f)

        # leave as object here
        td = Series([timedelta(days=i) for i in range(3)] + ['foo'])
        self.assertEqual(td.dtype, 'object')

        # these will correctly infer a timedelta
        s = Series([None, pd.NaT, '1 Day'])
        self.assertEqual(s.dtype,'timedelta64[ns]')
        s = Series([np.nan, pd.NaT, '1 Day'])
        self.assertEqual(s.dtype,'timedelta64[ns]')
        s = Series([pd.NaT, None, '1 Day'])
        self.assertEqual(s.dtype,'timedelta64[ns]')
        s = Series([pd.NaT, np.nan, '1 Day'])
        self.assertEqual(s.dtype,'timedelta64[ns]')

    def test_operators_timedelta64(self):

        # invalid ops
        self.assertRaises(Exception, self.objSeries.__add__, 1)
        self.assertRaises(
            Exception, self.objSeries.__add__, np.array(1, dtype=np.int64))
        self.assertRaises(Exception, self.objSeries.__sub__, 1)
        self.assertRaises(
            Exception, self.objSeries.__sub__, np.array(1, dtype=np.int64))

        # seriese ops
        v1 = date_range('2012-1-1', periods=3, freq='D')
        v2 = date_range('2012-1-2', periods=3, freq='D')
        rs = Series(v2) - Series(v1)
        xp = Series(1e9 * 3600 * 24, rs.index).astype(
            'int64').astype('timedelta64[ns]')
        assert_series_equal(rs, xp)
        self.assertEqual(rs.dtype, 'timedelta64[ns]')

        df = DataFrame(dict(A=v1))
        td = Series([timedelta(days=i) for i in range(3)])
        self.assertEqual(td.dtype, 'timedelta64[ns]')

        # series on the rhs
        result = df['A'] - df['A'].shift()
        self.assertEqual(result.dtype, 'timedelta64[ns]')

        result = df['A'] + td
        self.assertEqual(result.dtype, 'M8[ns]')

        # scalar Timestamp on rhs
        maxa = df['A'].max()
        tm.assertIsInstance(maxa, Timestamp)

        resultb = df['A'] - df['A'].max()
        self.assertEqual(resultb.dtype, 'timedelta64[ns]')

        # timestamp on lhs
        result = resultb + df['A']
        values = [Timestamp('20111230'), Timestamp('20120101'), Timestamp('20120103')]
        expected = Series(values, name='A')
        assert_series_equal(result, expected)

        # datetimes on rhs
        result = df['A'] - datetime(2001, 1, 1)
        expected = Series([timedelta(days=4017 + i) for i in range(3)], name='A')
        assert_series_equal(result, expected)
        self.assertEqual(result.dtype, 'm8[ns]')

        d = datetime(2001, 1, 1, 3, 4)
        resulta = df['A'] - d
        self.assertEqual(resulta.dtype, 'm8[ns]')

        # roundtrip
        resultb = resulta + d
        assert_series_equal(df['A'], resultb)

        # timedeltas on rhs
        td = timedelta(days=1)
        resulta = df['A'] + td
        resultb = resulta - td
        assert_series_equal(resultb, df['A'])
        self.assertEqual(resultb.dtype, 'M8[ns]')

        # roundtrip
        td = timedelta(minutes=5, seconds=3)
        resulta = df['A'] + td
        resultb = resulta - td
        assert_series_equal(df['A'], resultb)
        self.assertEqual(resultb.dtype, 'M8[ns]')

        # inplace
        value = rs[2] + np.timedelta64(timedelta(minutes=5,seconds=1))
        rs[2] += np.timedelta64(timedelta(minutes=5,seconds=1))
        self.assertEqual(rs[2], value)

    def test_timedeltas_with_DateOffset(self):

        # GH 4532
        # operate with pd.offsets
        s = Series([Timestamp('20130101 9:01'), Timestamp('20130101 9:02')])

        result = s + pd.offsets.Second(5)
        result2 = pd.offsets.Second(5) + s
        expected = Series(
            [Timestamp('20130101 9:01:05'), Timestamp('20130101 9:02:05')])
        assert_series_equal(result, expected)
        assert_series_equal(result2, expected)

        result = s + pd.offsets.Milli(5)
        result2 = pd.offsets.Milli(5) + s
        expected = Series(
            [Timestamp('20130101 9:01:00.005'), Timestamp('20130101 9:02:00.005')])
        assert_series_equal(result, expected)
        assert_series_equal(result2, expected)

        result = s + pd.offsets.Minute(5) + pd.offsets.Milli(5)
        expected = Series(
            [Timestamp('20130101 9:06:00.005'), Timestamp('20130101 9:07:00.005')])
        assert_series_equal(result, expected)

        # operate with np.timedelta64 correctly
        result = s + np.timedelta64(1, 's')
        result2 = np.timedelta64(1, 's') + s
        expected = Series(
            [Timestamp('20130101 9:01:01'), Timestamp('20130101 9:02:01')])
        assert_series_equal(result, expected)
        assert_series_equal(result2, expected)

        result = s + np.timedelta64(5, 'ms')
        result2 = np.timedelta64(5, 'ms') + s
        expected = Series(
            [Timestamp('20130101 9:01:00.005'), Timestamp('20130101 9:02:00.005')])
        assert_series_equal(result, expected)
        assert_series_equal(result2, expected)

        # valid DateOffsets
        for do in [ 'Hour', 'Minute', 'Second', 'Day', 'Micro',
                    'Milli', 'Nano' ]:
            op = getattr(pd.offsets,do)
            s + op(5)
            op(5) + s


    def test_timedelta64_operations_with_DateOffset(self):
        # GH 10699
        td = Series([timedelta(minutes=5, seconds=3)] * 3)
        result = td + pd.offsets.Minute(1)
        expected = Series([timedelta(minutes=6, seconds=3)] * 3)
        assert_series_equal(result, expected)

        result = td - pd.offsets.Minute(1)
        expected = Series([timedelta(minutes=4, seconds=3)] * 3)
        assert_series_equal(result, expected)

        result = td + Series([pd.offsets.Minute(1), pd.offsets.Second(3),
                              pd.offsets.Hour(2)])
        expected = Series([timedelta(minutes=6, seconds=3),
                           timedelta(minutes=5, seconds=6),
                           timedelta(hours=2, minutes=5, seconds=3)])
        assert_series_equal(result, expected)

        result = td + pd.offsets.Minute(1) + pd.offsets.Second(12)
        expected = Series([timedelta(minutes=6, seconds=15)] * 3)
        assert_series_equal(result, expected)

        # valid DateOffsets
        for do in [ 'Hour', 'Minute', 'Second', 'Day', 'Micro',
                    'Milli', 'Nano' ]:
            op = getattr(pd.offsets,do)
            td + op(5)
            op(5) + td
            td - op(5)
            op(5) - td

    def test_timedelta64_operations_with_timedeltas(self):

        # td operate with td
        td1 = Series([timedelta(minutes=5, seconds=3)] * 3)
        td2 = timedelta(minutes=5, seconds=4)
        result = td1 - td2
        expected = Series([timedelta(seconds=0)] * 3) -Series(
            [timedelta(seconds=1)] * 3)
        self.assertEqual(result.dtype, 'm8[ns]')
        assert_series_equal(result, expected)

        result2 = td2 - td1
        expected = (Series([timedelta(seconds=1)] * 3) -
                    Series([timedelta(seconds=0)] * 3))
        assert_series_equal(result2, expected)

        # roundtrip
        assert_series_equal(result + td2,td1)

        # Now again, using pd.to_timedelta, which should build
        # a Series or a scalar, depending on input.
        td1 = Series(pd.to_timedelta(['00:05:03'] * 3))
        td2 = pd.to_timedelta('00:05:04')
        result = td1 - td2
        expected = Series([timedelta(seconds=0)] * 3) -Series(
            [timedelta(seconds=1)] * 3)
        self.assertEqual(result.dtype, 'm8[ns]')
        assert_series_equal(result, expected)

        result2 = td2 - td1
        expected = (Series([timedelta(seconds=1)] * 3) -
                    Series([timedelta(seconds=0)] * 3))
        assert_series_equal(result2, expected)

        # roundtrip
        assert_series_equal(result + td2,td1)

    def test_timedelta64_operations_with_integers(self):

        # GH 4521
        # divide/multiply by integers
        startdate = Series(date_range('2013-01-01', '2013-01-03'))
        enddate = Series(date_range('2013-03-01', '2013-03-03'))

        s1 = enddate - startdate
        s1[2] = np.nan
        s2 = Series([2, 3, 4])
        expected = Series(s1.values.astype(np.int64) / s2, dtype='m8[ns]')
        expected[2] = np.nan
        result = s1 / s2
        assert_series_equal(result,expected)

        s2 = Series([20, 30, 40])
        expected = Series(s1.values.astype(np.int64) / s2, dtype='m8[ns]')
        expected[2] = np.nan
        result = s1 / s2
        assert_series_equal(result,expected)

        result = s1 / 2
        expected = Series(s1.values.astype(np.int64) / 2, dtype='m8[ns]')
        expected[2] = np.nan
        assert_series_equal(result,expected)

        s2 = Series([20, 30, 40])
        expected = Series(s1.values.astype(np.int64) * s2, dtype='m8[ns]')
        expected[2] = np.nan
        result = s1 * s2
        assert_series_equal(result,expected)

        for dtype in ['int32','int16','uint32','uint64','uint32','uint16','uint8']:
            s2 = Series([20, 30, 40],dtype=dtype)
            expected = Series(s1.values.astype(np.int64) * s2.astype(np.int64), dtype='m8[ns]')
            expected[2] = np.nan
            result = s1 * s2
            assert_series_equal(result,expected)

        result = s1 * 2
        expected = Series(s1.values.astype(np.int64) * 2, dtype='m8[ns]')
        expected[2] = np.nan
        assert_series_equal(result,expected)

        result = s1 * -1
        expected = Series(s1.values.astype(np.int64) * -1, dtype='m8[ns]')
        expected[2] = np.nan
        assert_series_equal(result,expected)

        # invalid ops
        for op in ['__true_div__','__div__','__mul__']:
            sop = getattr(s1,op,None)
            if sop is not None:
                self.assertRaises(TypeError, sop, s2.astype(float))
                self.assertRaises(TypeError, sop, 2.)

        for op in ['__add__','__sub__']:
            sop = getattr(s1,op,None)
            if sop is not None:
                self.assertRaises(TypeError, sop, 1)
                self.assertRaises(TypeError, sop, s2.values)

    def test_timedelta64_conversions(self):
        startdate = Series(date_range('2013-01-01', '2013-01-03'))
        enddate = Series(date_range('2013-03-01', '2013-03-03'))

        s1 = enddate - startdate
        s1[2] = np.nan

        for m in [1, 3, 10]:
            for unit in ['D','h','m','s','ms','us','ns']:

                # op
                expected = s1.apply(lambda x: x / np.timedelta64(m,unit))
                result = s1 / np.timedelta64(m,unit)
                assert_series_equal(result, expected)

                if m == 1 and unit != 'ns':

                    # astype
                    result = s1.astype("timedelta64[{0}]".format(unit))
                    assert_series_equal(result, expected)

                # reverse op
                expected = s1.apply(lambda x: np.timedelta64(m,unit) / x)
                result = np.timedelta64(m,unit) / s1

        # astype
        s = Series(date_range('20130101',periods=3))
        result = s.astype(object)
        self.assertIsInstance(result.iloc[0],datetime)
        self.assertTrue(result.dtype == np.object_)

        result = s1.astype(object)
        self.assertIsInstance(result.iloc[0],timedelta)
        self.assertTrue(result.dtype == np.object_)

    def test_timedelta64_equal_timedelta_supported_ops(self):
        ser = Series([Timestamp('20130301'), Timestamp('20130228 23:00:00'),
                      Timestamp('20130228 22:00:00'),
                      Timestamp('20130228 21:00:00')])

        intervals = 'D', 'h', 'm', 's', 'us'
        npy16_mappings = {'D': 24 * 60 * 60 * 1000000, 'h': 60 * 60 * 1000000,
                          'm': 60 * 1000000, 's': 1000000, 'us': 1}

        def timedelta64(*args):
            return sum(starmap(np.timedelta64, zip(args, intervals)))

        for op, d, h, m, s, us in product([operator.add, operator.sub],
                                          *([range(2)] * 5)):
            nptd = timedelta64(d, h, m, s, us)
            pytd = timedelta(days=d, hours=h, minutes=m, seconds=s,
                             microseconds=us)
            lhs = op(ser, nptd)
            rhs = op(ser, pytd)

            try:
                assert_series_equal(lhs, rhs)
            except:
                raise AssertionError(
                    "invalid comparsion [op->{0},d->{1},h->{2},m->{3},s->{4},us->{5}]\n{6}\n{7}\n".format(op, d, h, m, s, us, lhs, rhs))

    def test_timedelta_assignment(self):
        # GH 8209
        s = Series([])
        s.loc['B'] = timedelta(1)
        tm.assert_series_equal(s,Series(Timedelta('1 days'),index=['B']))

        s = s.reindex(s.index.insert(0, 'A'))
        tm.assert_series_equal(s,Series([np.nan,Timedelta('1 days')],index=['A','B']))

        result = s.fillna(timedelta(1))
        expected = Series(Timedelta('1 days'),index=['A','B'])
        tm.assert_series_equal(result, expected)

        s.loc['A'] = timedelta(1)
        tm.assert_series_equal(s, expected)

    def test_operators_datetimelike(self):

        def run_ops(ops, get_ser, test_ser):

            # check that we are getting a TypeError
            # with 'operate' (from core/ops.py) for the ops that are not defined
            for op_str in ops:
                op = getattr(get_ser, op_str, None)
                with tm.assertRaisesRegexp(TypeError, 'operate'):
                    op(test_ser)

        ### timedelta64 ###
        td1 = Series([timedelta(minutes=5,seconds=3)]*3)
        td1.iloc[2] = np.nan
        td2 = timedelta(minutes=5,seconds=4)
        ops = ['__mul__','__floordiv__','__pow__',
               '__rmul__','__rfloordiv__','__rpow__']
        run_ops(ops, td1, td2)
        td1 + td2
        td2 + td1
        td1 - td2
        td2 - td1
        td1 / td2
        td2 / td1

        ### datetime64 ###
        dt1 = Series([Timestamp('20111230'), Timestamp('20120101'),
                      Timestamp('20120103')])
        dt1.iloc[2] = np.nan
        dt2 = Series([Timestamp('20111231'), Timestamp('20120102'),
                      Timestamp('20120104')])
        ops = ['__add__', '__mul__', '__floordiv__', '__truediv__', '__div__',
               '__pow__', '__radd__', '__rmul__', '__rfloordiv__',
               '__rtruediv__', '__rdiv__', '__rpow__']
        run_ops(ops, dt1, dt2)
        dt1 - dt2
        dt2 - dt1

        ### datetime64 with timetimedelta ###
        ops = ['__mul__', '__floordiv__', '__truediv__', '__div__', '__pow__',
               '__rmul__', '__rfloordiv__', '__rtruediv__', '__rdiv__',
               '__rpow__']
        run_ops(ops, dt1, td1)
        dt1 + td1
        td1 + dt1
        dt1 - td1
        # TODO: Decide if this ought to work.
        # td1 - dt1

        ### timetimedelta with datetime64 ###
        ops = ['__sub__', '__mul__', '__floordiv__', '__truediv__', '__div__',
               '__pow__', '__rsub__', '__rmul__', '__rfloordiv__',
               '__rtruediv__', '__rdiv__', '__rpow__']
        run_ops(ops, td1, dt1)
        td1 + dt1
        dt1 + td1

        # 8260, 10763
        # datetime64 with tz
        ops = ['__mul__', '__floordiv__', '__truediv__', '__div__', '__pow__',
               '__rmul__', '__rfloordiv__', '__rtruediv__', '__rdiv__',
               '__rpow__']
        dt1 = Series(date_range('2000-01-01 09:00:00',periods=5,tz='US/Eastern'),name='foo')
        dt2 = dt1.copy()
        dt2.iloc[2] = np.nan
        td1 = Series(timedelta_range('1 days 1 min',periods=5, freq='H'))
        td2 = td1.copy()
        td2.iloc[1] = np.nan
        run_ops(ops, dt1, td1)

        result = dt1 + td1[0]
        expected = (dt1.dt.tz_localize(None) + td1[0]).dt.tz_localize('US/Eastern')
        assert_series_equal(result, expected)

        result = dt2 + td2[0]
        expected = (dt2.dt.tz_localize(None) + td2[0]).dt.tz_localize('US/Eastern')
        assert_series_equal(result, expected)

        # odd numpy behavior with scalar timedeltas
        if not _np_version_under1p8:
            result = td1[0] + dt1
            expected = (dt1.dt.tz_localize(None) + td1[0]).dt.tz_localize('US/Eastern')
            assert_series_equal(result, expected)

            result = td2[0] + dt2
            expected = (dt2.dt.tz_localize(None) + td2[0]).dt.tz_localize('US/Eastern')
            assert_series_equal(result, expected)

        result = dt1 - td1[0]
        expected = (dt1.dt.tz_localize(None) - td1[0]).dt.tz_localize('US/Eastern')
        assert_series_equal(result, expected)
        self.assertRaises(TypeError, lambda: td1[0] - dt1)

        result = dt2 - td2[0]
        expected = (dt2.dt.tz_localize(None) - td2[0]).dt.tz_localize('US/Eastern')
        assert_series_equal(result, expected)
        self.assertRaises(TypeError, lambda: td2[0] - dt2)

        result = dt1 + td1
        expected = (dt1.dt.tz_localize(None) + td1).dt.tz_localize('US/Eastern')
        assert_series_equal(result, expected)

        result = dt2 + td2
        expected = (dt2.dt.tz_localize(None) + td2).dt.tz_localize('US/Eastern')
        assert_series_equal(result, expected)

        result = dt1 - td1
        expected = (dt1.dt.tz_localize(None) - td1).dt.tz_localize('US/Eastern')
        assert_series_equal(result, expected)

        result = dt2 - td2
        expected = (dt2.dt.tz_localize(None) - td2).dt.tz_localize('US/Eastern')
        assert_series_equal(result, expected)

        self.assertRaises(TypeError, lambda: td1 - dt1)
        self.assertRaises(TypeError, lambda: td2 - dt2)

    def test_ops_datetimelike_align(self):
        # GH 7500
        # datetimelike ops need to align
        dt = Series(date_range('2012-1-1', periods=3, freq='D'))
        dt.iloc[2] = np.nan
        dt2 = dt[::-1]

        expected = Series([timedelta(0), timedelta(0), pd.NaT])
        # name is reset
        result = dt2 - dt
        assert_series_equal(result, expected)

        expected = Series(expected, name=0)
        result = (dt2.to_frame() - dt.to_frame())[0]
        assert_series_equal(result, expected)

    def test_timedelta64_functions(self):

        from datetime import timedelta
        from pandas import date_range

        # index min/max
        td = Series(date_range('2012-1-1', periods=3, freq='D')) - \
                    Timestamp('20120101')

        result = td.idxmin()
        self.assertEqual(result, 0)

        result = td.idxmax()
        self.assertEqual(result, 2)

        # GH 2982
        # with NaT
        td[0] = np.nan

        result = td.idxmin()
        self.assertEqual(result, 1)

        result = td.idxmax()
        self.assertEqual(result, 2)

        # abs
        s1 = Series(date_range('20120101', periods=3))
        s2 = Series(date_range('20120102', periods=3))
        expected = Series(s2 - s1)

        # this fails as numpy returns timedelta64[us]
        #result = np.abs(s1-s2)
        # assert_frame_equal(result,expected)

        result = (s1 - s2).abs()
        assert_series_equal(result, expected)

        # max/min
        result = td.max()
        expected = Timedelta('2 days')
        self.assertEqual(result, expected)

        result = td.min()
        expected = Timedelta('1 days')
        self.assertEqual(result, expected)

    def test_ops_consistency_on_empty(self):

        # GH 7869
        # consistency on empty

        # float
        result = Series(dtype=float).sum()
        self.assertEqual(result,0)

        result = Series(dtype=float).mean()
        self.assertTrue(isnull(result))

        result = Series(dtype=float).median()
        self.assertTrue(isnull(result))

        # timedelta64[ns]
        result = Series(dtype='m8[ns]').sum()
        self.assertEqual(result, Timedelta(0))

        result = Series(dtype='m8[ns]').mean()
        self.assertTrue(result is pd.NaT)

        result = Series(dtype='m8[ns]').median()
        self.assertTrue(result is pd.NaT)

    def test_timedelta_fillna(self):
        #GH 3371
        s = Series([Timestamp('20130101'), Timestamp('20130101'),
                    Timestamp('20130102'), Timestamp('20130103 9:01:01')])
        td = s.diff()

        # reg fillna
        result = td.fillna(0)
        expected = Series([timedelta(0), timedelta(0), timedelta(1),
                           timedelta(days=1, seconds=9*3600+60+1)])
        assert_series_equal(result, expected)

        # interprested as seconds
        result = td.fillna(1)
        expected = Series([timedelta(seconds=1), timedelta(0),
                           timedelta(1), timedelta(days=1, seconds=9*3600+60+1)])
        assert_series_equal(result, expected)

        result = td.fillna(timedelta(days=1, seconds=1))
        expected = Series([timedelta(days=1, seconds=1), timedelta(0),
                           timedelta(1), timedelta(days=1, seconds=9*3600+60+1)])
        assert_series_equal(result, expected)

        result = td.fillna(np.timedelta64(int(1e9)))
        expected = Series([timedelta(seconds=1), timedelta(0), timedelta(1),
                           timedelta(days=1, seconds=9*3600+60+1)])
        assert_series_equal(result, expected)

        from pandas import tslib
        result = td.fillna(tslib.NaT)
        expected = Series([tslib.NaT, timedelta(0), timedelta(1),
                           timedelta(days=1, seconds=9*3600+60+1)], dtype='m8[ns]')
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
        expected[2] = timedelta(days=1, seconds=9*3600+60+1)
        assert_series_equal(result, expected)

    def test_datetime64_fillna(self):

        s = Series([Timestamp('20130101'), Timestamp('20130101'),
                    Timestamp('20130102'), Timestamp('20130103 9:01:01')])
        s[2] = np.nan

        # reg fillna
        result = s.fillna(Timestamp('20130104'))
        expected = Series([Timestamp('20130101'), Timestamp('20130101'),
                           Timestamp('20130104'), Timestamp('20130103 9:01:01')])
        assert_series_equal(result, expected)

        from pandas import tslib
        result = s.fillna(tslib.NaT)
        expected = s
        assert_series_equal(result, expected)

        # ffill
        result = s.ffill()
        expected = Series([Timestamp('20130101'), Timestamp('20130101'),
                           Timestamp('20130101'), Timestamp('20130103 9:01:01')])
        assert_series_equal(result, expected)

        # bfill
        result = s.bfill()
        expected = Series([Timestamp('20130101'), Timestamp('20130101'),
                           Timestamp('20130103 9:01:01'),
                           Timestamp('20130103 9:01:01')])
        assert_series_equal(result, expected)

        # GH 6587
        # make sure that we are treating as integer when filling
        # this also tests inference of a datetime-like with NaT's
        s = Series([pd.NaT, pd.NaT, '2013-08-05 15:30:00.000001'])
        expected = Series(['2013-08-05 15:30:00.000001', '2013-08-05 15:30:00.000001', '2013-08-05 15:30:00.000001'], dtype='M8[ns]')
        result = s.fillna(method='backfill')
        assert_series_equal(result, expected)

    def test_datetime64_tz_fillna(self):
        for tz in ['US/Eastern', 'Asia/Tokyo']:
            # DatetimeBlock
            s = Series([Timestamp('2011-01-01 10:00'), pd.NaT,
                        Timestamp('2011-01-03 10:00'), pd.NaT])
            result = s.fillna(pd.Timestamp('2011-01-02 10:00'))
            expected = Series([Timestamp('2011-01-01 10:00'), Timestamp('2011-01-02 10:00'),
                               Timestamp('2011-01-03 10:00'), Timestamp('2011-01-02 10:00')])
            self.assert_series_equal(expected, result)

            result = s.fillna(pd.Timestamp('2011-01-02 10:00', tz=tz))
            expected = Series([Timestamp('2011-01-01 10:00'),
                               Timestamp('2011-01-02 10:00', tz=tz),
                               Timestamp('2011-01-03 10:00'),
                               Timestamp('2011-01-02 10:00', tz=tz)])
            self.assert_series_equal(expected, result)

            result = s.fillna('AAA')
            expected = Series([Timestamp('2011-01-01 10:00'), 'AAA',
                               Timestamp('2011-01-03 10:00'), 'AAA'], dtype=object)
            self.assert_series_equal(expected, result)

            result = s.fillna({1: pd.Timestamp('2011-01-02 10:00', tz=tz),
                               3: pd.Timestamp('2011-01-04 10:00')})
            expected = Series([Timestamp('2011-01-01 10:00'),
                               Timestamp('2011-01-02 10:00', tz=tz),
                               Timestamp('2011-01-03 10:00'),
                               Timestamp('2011-01-04 10:00')])
            self.assert_series_equal(expected, result)

            result = s.fillna({1: pd.Timestamp('2011-01-02 10:00'),
                               3: pd.Timestamp('2011-01-04 10:00')})
            expected = Series([Timestamp('2011-01-01 10:00'), Timestamp('2011-01-02 10:00'),
                               Timestamp('2011-01-03 10:00'), Timestamp('2011-01-04 10:00')])
            self.assert_series_equal(expected, result)

            # DatetimeBlockTZ
            idx = pd.DatetimeIndex(['2011-01-01 10:00', pd.NaT,
                                    '2011-01-03 10:00', pd.NaT], tz=tz)
            s = pd.Series(idx)
            result = s.fillna(pd.Timestamp('2011-01-02 10:00'))
            expected = Series([Timestamp('2011-01-01 10:00', tz=tz),
                               Timestamp('2011-01-02 10:00'),
                               Timestamp('2011-01-03 10:00', tz=tz),
                               Timestamp('2011-01-02 10:00')])
            self.assert_series_equal(expected, result)

            result = s.fillna(pd.Timestamp('2011-01-02 10:00', tz=tz))
            idx = pd.DatetimeIndex(['2011-01-01 10:00', '2011-01-02 10:00',
                                    '2011-01-03 10:00', '2011-01-02 10:00'],
                                   tz=tz)
            expected = Series(idx)
            self.assert_series_equal(expected, result)

            result = s.fillna(pd.Timestamp('2011-01-02 10:00', tz=tz).to_pydatetime())
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
            expected = Series([Timestamp('2011-01-01 10:00', tz=tz),
                               Timestamp('2011-01-02 10:00', tz=tz),
                               Timestamp('2011-01-03 10:00', tz=tz),
                               Timestamp('2011-01-04 10:00')])
            self.assert_series_equal(expected, result)

            result = s.fillna({1: pd.Timestamp('2011-01-02 10:00', tz=tz),
                               3: pd.Timestamp('2011-01-04 10:00', tz=tz)})
            expected = Series([Timestamp('2011-01-01 10:00', tz=tz),
                               Timestamp('2011-01-02 10:00', tz=tz),
                               Timestamp('2011-01-03 10:00', tz=tz),
                               Timestamp('2011-01-04 10:00', tz=tz)])
            self.assert_series_equal(expected, result)

            # filling with a naive/other zone, coerce to object
            result = s.fillna(Timestamp('20130101'))
            expected = Series([Timestamp('2011-01-01 10:00', tz=tz),
                               Timestamp('2013-01-01'),
                               Timestamp('2011-01-03 10:00', tz=tz),
                               Timestamp('2013-01-01')])
            self.assert_series_equal(expected, result)

            result = s.fillna(Timestamp('20130101',tz='US/Pacific'))
            expected = Series([Timestamp('2011-01-01 10:00', tz=tz),
                               Timestamp('2013-01-01',tz='US/Pacific'),
                               Timestamp('2011-01-03 10:00', tz=tz),
                               Timestamp('2013-01-01',tz='US/Pacific')])
            self.assert_series_equal(expected, result)

    def test_fillna_int(self):
        s = Series(np.random.randint(-100, 100, 50))
        s.fillna(method='ffill', inplace=True)
        assert_series_equal(s.fillna(method='ffill', inplace=False), s)

    def test_fillna_raise(self):
        s = Series(np.random.randint(-100, 100, 50))
        self.assertRaises(TypeError, s.fillna, [1, 2])
        self.assertRaises(TypeError, s.fillna, (1, 2))

    def test_raise_on_info(self):
        s = Series(np.random.randn(10))
        with tm.assertRaises(AttributeError):
            s.info()

    def test_isnull_for_inf(self):
        s = Series(['a', np.inf, np.nan, 1.0])
        with pd.option_context('mode.use_inf_as_null', True):
            r = s.isnull()
            dr = s.dropna()
        e = Series([False, True, True, False])
        de = Series(['a', 1.0], index=[0, 3])
        tm.assert_series_equal(r, e)
        tm.assert_series_equal(dr, de)


# TimeSeries-specific

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
        assert_series_equal(result,expected)
        result = s1.fillna({})
        assert_series_equal(result,s1)
        result = s1.fillna(Series(()))
        assert_series_equal(result,s1)
        result = s2.fillna(s1)
        assert_series_equal(result,s2)
        result = s1.fillna({ 0 : 1})
        assert_series_equal(result,expected)
        result = s1.fillna({ 1 : 1})
        assert_series_equal(result,Series([np.nan]))
        result = s1.fillna({ 0 : 1, 1 : 1})
        assert_series_equal(result,expected)
        result = s1.fillna(Series({ 0 : 1, 1 : 1}))
        assert_series_equal(result,expected)
        result = s1.fillna(Series({ 0 : 1, 1 : 1},index=[4,5]))
        assert_series_equal(result,s1)

        s1 = Series([0, 1, 2], list('abc'))
        s2 = Series([0, np.nan, 2], list('bac'))
        result = s2.fillna(s1)
        expected = Series([0,0,2.], list('bac'))
        assert_series_equal(result,expected)

        # limit
        s = Series(np.nan,index=[0,1,2])
        result = s.fillna(999,limit=1)
        expected = Series([999,np.nan,np.nan],index=[0,1,2])
        assert_series_equal(result,expected)

        result = s.fillna(999,limit=2)
        expected = Series([999,999,np.nan],index=[0,1,2])
        assert_series_equal(result,expected)

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

    def test_sub_of_datetime_from_TimeSeries(self):
        from pandas.tseries.timedeltas import to_timedelta
        from datetime import datetime
        a = Timestamp(datetime(1993, 0o1, 0o7, 13, 30, 00))
        b = datetime(1993, 6, 22, 13, 30)
        a = Series([a])
        result = to_timedelta(np.abs(a - b))
        self.assertEqual(result.dtype, 'timedelta64[ns]')

    def test_datetime64_with_index(self):

        # arithmetic integer ops with an index
        s = Series(np.random.randn(5))
        expected = s - s.index.to_series()
        result = s - s.index
        assert_series_equal(result, expected)

        # GH 4629
        # arithmetic datetime64 ops with an index
        s = Series(date_range('20130101', periods=5),
                   index=date_range('20130101', periods=5))
        expected = s - s.index.to_series()
        result = s - s.index
        assert_series_equal(result, expected)

        result = s - s.index.to_period()
        assert_series_equal(result, expected)

        df = DataFrame(np.random.randn(5,2),
                      index=date_range('20130101', periods=5))
        df['date'] = Timestamp('20130102')
        df['expected'] = df['date'] - df.index.to_series()
        df['result'] = df['date'] - df.index
        assert_series_equal(df['result'], df['expected'], check_names=False)

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
        #result = td[(td>np.timedelta64(timedelta(days=3))) & (td<np.timedelta64(timedelta(days=7)))] = np.nan
        #self.assertEqual(isnull(result).sum(), 7)

    # NumPy limitiation =(

    # def test_logical_range_select(self):
    #     np.random.seed(12345)
    #     selector = -0.5 <= self.ts <= 0.5
    #     expected = (self.ts >= -0.5) & (self.ts <= 0.5)
    #     assert_series_equal(selector, expected)

    def test_operators_na_handling(self):
        from decimal import Decimal
        from datetime import date
        s = Series([Decimal('1.3'), Decimal('2.3')],
                   index=[date(2012, 1, 1), date(2012, 1, 2)])

        result = s + s.shift(1)
        result2 = s.shift(1) + s
        self.assertTrue(isnull(result[0]))
        self.assertTrue(isnull(result2[0]))

        s = Series(['foo', 'bar', 'baz', np.nan])
        result = 'prefix_' + s
        expected = Series(['prefix_foo', 'prefix_bar', 'prefix_baz', np.nan])
        assert_series_equal(result, expected)

        result = s + '_suffix'
        expected = Series(['foo_suffix', 'bar_suffix', 'baz_suffix', np.nan])
        assert_series_equal(result, expected)

    def test_object_comparisons(self):
        s = Series(['a', 'b', np.nan, 'c', 'a'])

        result = s == 'a'
        expected = Series([True, False, False, False, True])
        assert_series_equal(result, expected)

        result = s < 'a'
        expected = Series([False, False, False, False, False])
        assert_series_equal(result, expected)

        result = s != 'a'
        expected = -(s == 'a')
        assert_series_equal(result, expected)

    def test_comparison_tuples(self):
        # GH11339
        # comparisons vs tuple
        s = Series([(1,1),(1,2)])

        result = s == (1,2)
        expected = Series([False,True])
        assert_series_equal(result, expected)

        result = s != (1,2)
        expected = Series([True, False])
        assert_series_equal(result, expected)

        result = s == (0,0)
        expected = Series([False, False])
        assert_series_equal(result, expected)

        result = s != (0,0)
        expected = Series([True, True])
        assert_series_equal(result, expected)

        s = Series([(1,1),(1,1)])

        result = s == (1,1)
        expected = Series([True, True])
        assert_series_equal(result, expected)

        result = s != (1,1)
        expected = Series([False, False])
        assert_series_equal(result, expected)

        s = Series([frozenset([1]),frozenset([1,2])])

        result = s == frozenset([1])
        expected = Series([True, False])
        assert_series_equal(result, expected)

    def test_comparison_operators_with_nas(self):
        s = Series(bdate_range('1/1/2000', periods=10), dtype=object)
        s[::2] = np.nan

        # test that comparisons work
        ops = ['lt', 'le', 'gt', 'ge', 'eq', 'ne']
        for op in ops:
            val = s[5]

            f = getattr(operator, op)
            result = f(s, val)

            expected = f(s.dropna(), val).reindex(s.index)

            if op == 'ne':
                expected = expected.fillna(True).astype(bool)
            else:
                expected = expected.fillna(False).astype(bool)

            assert_series_equal(result, expected)

            # fffffffuuuuuuuuuuuu
            # result = f(val, s)
            # expected = f(val, s.dropna()).reindex(s.index)
            # assert_series_equal(result, expected)

        # boolean &, |, ^ should work with object arrays and propagate NAs

        ops = ['and_', 'or_', 'xor']
        mask = s.isnull()
        for bool_op in ops:
            f = getattr(operator, bool_op)

            filled = s.fillna(s[0])

            result = f(s < s[9], s > s[3])

            expected = f(filled < filled[9], filled > filled[3])
            expected[mask] = False
            assert_series_equal(result, expected)

    def test_comparison_object_numeric_nas(self):
        s = Series(np.random.randn(10), dtype=object)
        shifted = s.shift(2)

        ops = ['lt', 'le', 'gt', 'ge', 'eq', 'ne']
        for op in ops:
            f = getattr(operator, op)

            result = f(s, shifted)
            expected = f(s.astype(float), shifted.astype(float))
            assert_series_equal(result, expected)

    def test_comparison_invalid(self):

        # GH4968
        # invalid date/int comparisons
        s = Series(range(5))
        s2 = Series(date_range('20010101', periods=5))

        for (x, y) in [(s,s2),(s2,s)]:
            self.assertRaises(TypeError, lambda : x == y)
            self.assertRaises(TypeError, lambda : x != y)
            self.assertRaises(TypeError, lambda : x >= y)
            self.assertRaises(TypeError, lambda : x > y)
            self.assertRaises(TypeError, lambda : x < y)
            self.assertRaises(TypeError, lambda : x <= y)

    def test_more_na_comparisons(self):
        left = Series(['a', np.nan, 'c'])
        right = Series(['a', np.nan, 'd'])

        result = left == right
        expected = Series([True, False, False])
        assert_series_equal(result, expected)

        result = left != right
        expected = Series([False, True, True])
        assert_series_equal(result, expected)

        result = left == np.nan
        expected = Series([False, False, False])
        assert_series_equal(result, expected)

        result = left != np.nan
        expected = Series([True, True, True])
        assert_series_equal(result, expected)

    def test_comparison_different_length(self):
        a = Series(['a', 'b', 'c'])
        b = Series(['b', 'a'])
        self.assertRaises(ValueError, a.__lt__, b)

        a = Series([1, 2])
        b = Series([2, 3, 4])
        self.assertRaises(ValueError, a.__eq__, b)

    def test_comparison_label_based(self):

        # GH 4947
        # comparisons should be label based

        a = Series([True, False, True], list('bca'))
        b = Series([False, True, False], list('abc'))

        expected = Series([True, False, False], list('bca'))
        result = a & b
        assert_series_equal(result,expected)

        expected = Series([True, False, True], list('bca'))
        result = a | b
        assert_series_equal(result,expected)

        expected = Series([False, False, True], list('bca'))
        result = a ^ b
        assert_series_equal(result,expected)

        # rhs is bigger
        a = Series([True, False, True], list('bca'))
        b = Series([False, True, False, True], list('abcd'))

        expected = Series([True, False, False], list('bca'))
        result = a & b
        assert_series_equal(result,expected)

        expected = Series([True, False, True], list('bca'))
        result = a | b
        assert_series_equal(result,expected)

        # filling

        # vs empty
        result = a & Series([])
        expected = Series([False, False, False], list('bca'))
        assert_series_equal(result,expected)

        result = a | Series([])
        expected = Series([True, False, True], list('bca'))
        assert_series_equal(result,expected)

        # vs non-matching
        result = a & Series([1],['z'])
        expected = Series([False, False, False], list('bca'))
        assert_series_equal(result,expected)

        result = a | Series([1],['z'])
        expected = Series([True, False, True], list('bca'))
        assert_series_equal(result,expected)

        # identity
        # we would like s[s|e] == s to hold for any e, whether empty or not
        for e in [Series([]),Series([1],['z']),Series(['z']),Series(np.nan,b.index),Series(np.nan,a.index)]:
            result = a[a | e]
            assert_series_equal(result,a[a])

        # vs scalars
        index = list('bca')
        t = Series([True,False,True])

        for v in [True,1,2]:
            result = Series([True,False,True],index=index) | v
            expected = Series([True,True,True],index=index)
            assert_series_equal(result,expected)

        for v in [np.nan,'foo']:
            self.assertRaises(TypeError, lambda : t | v)

        for v in [False,0]:
            result = Series([True,False,True],index=index) | v
            expected = Series([True,False,True],index=index)
            assert_series_equal(result,expected)

        for v in [True,1]:
            result = Series([True,False,True],index=index) & v
            expected = Series([True,False,True],index=index)
            assert_series_equal(result,expected)

        for v in [False,0]:
            result = Series([True,False,True],index=index) & v
            expected = Series([False,False,False],index=index)
            assert_series_equal(result,expected)
        for v in [np.nan]:
            self.assertRaises(TypeError, lambda : t & v)

    def test_operators_bitwise(self):
        # GH 9016: support bitwise op for integer types
        index = list('bca')

        s_tft = Series([True, False, True], index=index)
        s_fff = Series([False, False, False], index=index)
        s_tff = Series([True, False, False], index=index)
        s_empty = Series([])
        s_0101 = Series([0,1,0,1])
        s_0123 = Series(range(4),dtype='int64')
        s_3333 = Series([3] * 4)
        s_4444 = Series([4] * 4)

        res = s_tft & s_empty
        expected = s_fff
        assert_series_equal(res, expected)

        res = s_tft | s_empty
        expected = s_tft
        assert_series_equal(res, expected)

        res = s_0123 & s_3333
        expected = Series(range(4),dtype='int64')
        assert_series_equal(res, expected)

        res = s_0123 | s_4444
        expected = Series(range(4, 8),dtype='int64')
        assert_series_equal(res, expected)

        s_a0b1c0 = Series([1], list('b'))

        res = s_tft & s_a0b1c0
        expected = s_tff
        assert_series_equal(res, expected)

        res = s_tft | s_a0b1c0
        expected = s_tft
        assert_series_equal(res, expected)

        n0 = 0
        res = s_tft & n0
        expected = s_fff
        assert_series_equal(res, expected)

        res = s_0123 & n0
        expected = Series([0] * 4)
        assert_series_equal(res, expected)

        n1 = 1
        res = s_tft & n1
        expected = s_tft
        assert_series_equal(res, expected)

        res = s_0123 & n1
        expected = Series([0, 1, 0, 1])
        assert_series_equal(res, expected)

        s_1111 = Series([1]*4, dtype='int8')
        res = s_0123 & s_1111
        expected = Series([0, 1, 0, 1], dtype='int64')
        assert_series_equal(res, expected)

        res = s_0123.astype(np.int16) | s_1111.astype(np.int32)
        expected = Series([1, 1, 3, 3], dtype='int32')
        assert_series_equal(res, expected)

        self.assertRaises(TypeError, lambda: s_1111 & 'a')
        self.assertRaises(TypeError, lambda: s_1111 & ['a','b','c','d'])
        self.assertRaises(TypeError, lambda: s_0123 & np.NaN)
        self.assertRaises(TypeError, lambda: s_0123 & 3.14)
        self.assertRaises(TypeError, lambda: s_0123 & [0.1, 4, 3.14, 2])

        # s_0123 will be all false now because of reindexing like s_tft
        assert_series_equal(s_tft & s_0123, Series([False] * 3, list('bca')))
        # s_tft will be all false now because of reindexing like s_0123
        assert_series_equal(s_0123 & s_tft, Series([False] * 4))
        assert_series_equal(s_0123 & False, Series([False] * 4))
        assert_series_equal(s_0123 ^ False, Series([False, True, True, True]))
        assert_series_equal(s_0123 & [False], Series([False] * 4))
        assert_series_equal(s_0123 & (False), Series([False] * 4))
        assert_series_equal(s_0123 & Series([False, np.NaN, False, False]), Series([False] * 4))

        s_ftft = Series([False, True, False, True])
        assert_series_equal(s_0123 & Series([0.1, 4, -3.14, 2]), s_ftft)

        s_abNd = Series(['a','b',np.NaN,'d'])
        res = s_0123 & s_abNd
        expected = s_ftft
        assert_series_equal(res, expected)

    def test_between(self):
        s = Series(bdate_range('1/1/2000', periods=20).asobject)
        s[::2] = np.nan

        result = s[s.between(s[3], s[17])]
        expected = s[3:18].dropna()
        assert_series_equal(result, expected)

        result = s[s.between(s[3], s[17], inclusive=False)]
        expected = s[5:16].dropna()
        assert_series_equal(result, expected)

    def test_setitem_na(self):
        # these induce dtype changes
        expected = Series([np.nan, 3, np.nan, 5, np.nan, 7, np.nan, 9, np.nan])
        s = Series([2, 3, 4, 5, 6, 7, 8, 9, 10])
        s[::2] = np.nan
        assert_series_equal(s, expected)

        # get's coerced to float, right?
        expected = Series([np.nan, 1, np.nan, 0])
        s = Series([True, True, False, False])
        s[::2] = np.nan
        assert_series_equal(s, expected)

        expected = Series([np.nan, np.nan, np.nan, np.nan, np.nan, 5, 6, 7, 8, 9])
        s = Series(np.arange(10))
        s[:5] = np.nan
        assert_series_equal(s, expected)

    def test_scalar_na_cmp_corners(self):
        s = Series([2, 3, 4, 5, 6, 7, 8, 9, 10])

        def tester(a, b):
            return a & b

        self.assertRaises(TypeError, tester, s, datetime(2005, 1, 1))

        s = Series([2, 3, 4, 5, 6, 7, 8, 9, datetime(2005, 1, 1)])
        s[::2] = np.nan

        expected = Series(True,index=s.index)
        expected[::2] = False
        assert_series_equal(tester(s, list(s)), expected)

        d = DataFrame({'A': s})
        # TODO: Fix this exception - needs to be fixed! (see GH5035)
        # (previously this was a TypeError because series returned
        # NotImplemented
        self.assertRaises(ValueError, tester, s, d)

    def test_idxmin(self):
        # test idxmin
        # _check_stat_op approach can not be used here because of isnull check.

        # add some NaNs
        self.series[5:15] = np.NaN

        # skipna or no
        self.assertEqual(self.series[self.series.idxmin()], self.series.min())
        self.assertTrue(isnull(self.series.idxmin(skipna=False)))

        # no NaNs
        nona = self.series.dropna()
        self.assertEqual(nona[nona.idxmin()], nona.min())
        self.assertEqual(nona.index.values.tolist().index(nona.idxmin()),
                         nona.values.argmin())

        # all NaNs
        allna = self.series * nan
        self.assertTrue(isnull(allna.idxmin()))

        # datetime64[ns]
        from pandas import date_range
        s = Series(date_range('20130102', periods=6))
        result = s.idxmin()
        self.assertEqual(result, 0)

        s[0] = np.nan
        result = s.idxmin()
        self.assertEqual(result, 1)

    def test_idxmax(self):
        # test idxmax
        # _check_stat_op approach can not be used here because of isnull check.

        # add some NaNs
        self.series[5:15] = np.NaN

        # skipna or no
        self.assertEqual(self.series[self.series.idxmax()], self.series.max())
        self.assertTrue(isnull(self.series.idxmax(skipna=False)))

        # no NaNs
        nona = self.series.dropna()
        self.assertEqual(nona[nona.idxmax()], nona.max())
        self.assertEqual(nona.index.values.tolist().index(nona.idxmax()),
                         nona.values.argmax())

        # all NaNs
        allna = self.series * nan
        self.assertTrue(isnull(allna.idxmax()))

        from pandas import date_range
        s = Series(date_range('20130102', periods=6))
        result = s.idxmax()
        self.assertEqual(result, 5)

        s[5] = np.nan
        result = s.idxmax()
        self.assertEqual(result, 4)

        # Float64Index
        # GH 5914
        s = pd.Series([1,2,3],[1.1,2.1,3.1])
        result = s.idxmax()
        self.assertEqual(result, 3.1)
        result = s.idxmin()
        self.assertEqual(result, 1.1)

        s = pd.Series(s.index, s.index)
        result = s.idxmax()
        self.assertEqual(result, 3.1)
        result = s.idxmin()
        self.assertEqual(result, 1.1)

    def test_ndarray_compat(self):

        # test numpy compat with Series as sub-class of NDFrame
        tsdf = DataFrame(np.random.randn(1000, 3), columns=['A', 'B', 'C'],
                         index=date_range('1/1/2000', periods=1000))

        def f(x):
            return x[x.argmax()]
        result = tsdf.apply(f)
        expected = tsdf.max()
        assert_series_equal(result,expected)

        # .item()
        s = Series([1])
        result = s.item()
        self.assertEqual(result, 1)
        self.assertEqual(s.item(), s.iloc[0])

        # using an ndarray like function
        s = Series(np.random.randn(10))
        result = np.ones_like(s)
        expected = Series(1,index=range(10),dtype='float64')
        #assert_series_equal(result,expected)

        # ravel
        s = Series(np.random.randn(10))
        tm.assert_almost_equal(s.ravel(order='F'),s.values.ravel(order='F'))

        # compress
        # GH 6658
        s = Series([0, 1., -1], index=list('abc'))
        result = np.compress(s > 0, s)
        assert_series_equal(result, Series([1.], index=['b']))

        result = np.compress(s < -1, s)
        # result empty Index(dtype=object) as the same as original
        exp = Series([], dtype='float64', index=Index([], dtype='object'))
        assert_series_equal(result, exp)

        s = Series([0, 1., -1], index=[.1, .2, .3])
        result = np.compress(s > 0, s)
        assert_series_equal(result, Series([1.], index=[.2]))

        result = np.compress(s < -1, s)
        # result empty Float64Index as the same as original
        exp = Series([], dtype='float64', index=Index([], dtype='float64'))
        assert_series_equal(result, exp)

    def test_complexx(self):

        # GH4819
        # complex access for ndarray compat
        a = np.arange(5)
        b = Series(a + 4j*a)
        tm.assert_almost_equal(a,b.real)
        tm.assert_almost_equal(4*a,b.imag)

        b.real = np.arange(5)+5
        tm.assert_almost_equal(a+5,b.real)
        tm.assert_almost_equal(4*a,b.imag)

    def test_underlying_data_conversion(self):

        # GH 4080
        df = DataFrame(dict((c, [1,2,3]) for c in ['a', 'b', 'c']))
        df.set_index(['a', 'b', 'c'], inplace=True)
        s = Series([1], index=[(2,2,2)])
        df['val'] = 0
        df
        df['val'].update(s)

        expected = DataFrame(dict(a = [1,2,3], b = [1,2,3], c = [1,2,3], val = [0,1,0]))
        expected.set_index(['a', 'b', 'c'], inplace=True)
        tm.assert_frame_equal(df,expected)

        # GH 3970
        # these are chained assignments as well
        pd.set_option('chained_assignment',None)
        df = DataFrame({ "aa":range(5), "bb":[2.2]*5})
        df["cc"] = 0.0
        ck = [True]*len(df)
        df["bb"].iloc[0] = .13
        df_tmp = df.iloc[ck]
        df["bb"].iloc[0] = .15
        self.assertEqual(df['bb'].iloc[0], 0.15)
        pd.set_option('chained_assignment','raise')

        # GH 3217
        df = DataFrame(dict(a = [1,3], b = [np.nan, 2]))
        df['c'] = np.nan
        df['c'].update(pd.Series(['foo'],index=[0]))

        expected = DataFrame(dict(a = [1,3], b = [np.nan, 2], c = ['foo',np.nan]))
        tm.assert_frame_equal(df,expected)

    def test_operators_corner(self):
        series = self.ts

        empty = Series([], index=Index([]))

        result = series + empty
        self.assertTrue(np.isnan(result).all())

        result = empty + Series([], index=Index([]))
        self.assertEqual(len(result), 0)

        # TODO: this returned NotImplemented earlier, what to do?
        # deltas = Series([timedelta(1)] * 5, index=np.arange(5))
        # sub_deltas = deltas[::2]
        # deltas5 = deltas * 5
        # deltas = deltas + sub_deltas

        # float + int
        int_ts = self.ts.astype(int)[:-5]
        added = self.ts + int_ts
        expected = self.ts.values[:-5] + int_ts.values
        self.assert_numpy_array_equal(added[:-5], expected)

    def test_operators_reverse_object(self):
        # GH 56
        arr = Series(np.random.randn(10), index=np.arange(10),
                     dtype=object)

        def _check_op(arr, op):
            result = op(1., arr)
            expected = op(1., arr.astype(float))
            assert_series_equal(result.astype(float), expected)

        _check_op(arr, operator.add)
        _check_op(arr, operator.sub)
        _check_op(arr, operator.mul)
        _check_op(arr, operator.truediv)
        _check_op(arr, operator.floordiv)

    def test_series_frame_radd_bug(self):
        import operator

        # GH 353
        vals = Series(tm.rands_array(5, 10))
        result = 'foo_' + vals
        expected = vals.map(lambda x: 'foo_' + x)
        assert_series_equal(result, expected)

        frame = DataFrame({'vals': vals})
        result = 'foo_' + frame
        expected = DataFrame({'vals': vals.map(lambda x: 'foo_' + x)})
        tm.assert_frame_equal(result, expected)

        # really raise this time
        self.assertRaises(TypeError, operator.add, datetime.now(), self.ts)

    def test_operators_frame(self):
        # rpow does not work with DataFrame
        df = DataFrame({'A': self.ts})

        tm.assert_almost_equal(self.ts + self.ts, self.ts + df['A'])
        tm.assert_almost_equal(self.ts ** self.ts, self.ts ** df['A'])
        tm.assert_almost_equal(self.ts < self.ts, self.ts < df['A'])
        tm.assert_almost_equal(self.ts / self.ts, self.ts / df['A'])

    def test_operators_combine(self):
        def _check_fill(meth, op, a, b, fill_value=0):
            exp_index = a.index.union(b.index)
            a = a.reindex(exp_index)
            b = b.reindex(exp_index)

            amask = isnull(a)
            bmask = isnull(b)

            exp_values = []
            for i in range(len(exp_index)):
                if amask[i]:
                    if bmask[i]:
                        exp_values.append(nan)
                        continue
                    exp_values.append(op(fill_value, b[i]))
                elif bmask[i]:
                    if amask[i]:
                        exp_values.append(nan)
                        continue
                    exp_values.append(op(a[i], fill_value))
                else:
                    exp_values.append(op(a[i], b[i]))

            result = meth(a, b, fill_value=fill_value)
            expected = Series(exp_values, exp_index)
            assert_series_equal(result, expected)

        a = Series([nan, 1., 2., 3., nan], index=np.arange(5))
        b = Series([nan, 1, nan, 3, nan, 4.], index=np.arange(6))

        pairings = []
        for op in ['add', 'sub', 'mul', 'pow', 'truediv', 'floordiv']:
            fv = 0
            lop = getattr(Series, op)
            lequiv = getattr(operator, op)
            rop = getattr(Series, 'r' + op)
            # bind op at definition time...
            requiv = lambda x, y, op=op: getattr(operator, op)(y, x)
            pairings.append((lop, lequiv, fv))
            pairings.append((rop, requiv, fv))

        if compat.PY3:
            pairings.append((Series.div, operator.truediv, 1))
            pairings.append((Series.rdiv, lambda x, y: operator.truediv(y, x), 1))
        else:
            pairings.append((Series.div, operator.div, 1))
            pairings.append((Series.rdiv, lambda x, y: operator.div(y, x), 1))

        for op, equiv_op, fv in pairings:
            result = op(a, b)
            exp = equiv_op(a, b)
            assert_series_equal(result, exp)
            _check_fill(op, equiv_op, a, b, fill_value=fv)
            # should accept axis=0 or axis='rows'
            op(a, b, axis=0)

    def test_combine_first(self):
        values = tm.makeIntIndex(20).values.astype(float)
        series = Series(values, index=tm.makeIntIndex(20))

        series_copy = series * 2
        series_copy[::2] = np.NaN

        # nothing used from the input
        combined = series.combine_first(series_copy)

        self.assert_numpy_array_equal(combined, series)

        # Holes filled from input
        combined = series_copy.combine_first(series)
        self.assertTrue(np.isfinite(combined).all())

        self.assert_numpy_array_equal(combined[::2], series[::2])
        self.assert_numpy_array_equal(combined[1::2], series_copy[1::2])

        # mixed types
        index = tm.makeStringIndex(20)
        floats = Series(tm.randn(20), index=index)
        strings = Series(tm.makeStringIndex(10), index=index[::2])

        combined = strings.combine_first(floats)

        tm.assert_dict_equal(strings, combined, compare_keys=False)
        tm.assert_dict_equal(floats[1::2], combined, compare_keys=False)

        # corner case
        s = Series([1., 2, 3], index=[0, 1, 2])
        result = s.combine_first(Series([], index=[]))
        assert_series_equal(s, result)

    def test_update(self):
        s = Series([1.5, nan, 3., 4., nan])
        s2 = Series([nan, 3.5, nan, 5.])
        s.update(s2)

        expected = Series([1.5, 3.5, 3., 5., np.nan])
        assert_series_equal(s, expected)

        # GH 3217
        df = DataFrame([{"a": 1}, {"a": 3, "b": 2}])
        df['c'] = np.nan

        # this will fail as long as series is a sub-class of ndarray
        # df['c'].update(Series(['foo'],index=[0])) #####

    def test_corr(self):
        tm._skip_if_no_scipy()

        import scipy.stats as stats

        # full overlap
        self.assertAlmostEqual(self.ts.corr(self.ts), 1)

        # partial overlap
        self.assertAlmostEqual(self.ts[:15].corr(self.ts[5:]), 1)

        self.assertTrue(isnull(self.ts[:15].corr(self.ts[5:], min_periods=12)))

        ts1 = self.ts[:15].reindex(self.ts.index)
        ts2 = self.ts[5:].reindex(self.ts.index)
        self.assertTrue(isnull(ts1.corr(ts2, min_periods=12)))

        # No overlap
        self.assertTrue(np.isnan(self.ts[::2].corr(self.ts[1::2])))

        # all NA
        cp = self.ts[:10].copy()
        cp[:] = np.nan
        self.assertTrue(isnull(cp.corr(cp)))

        A = tm.makeTimeSeries()
        B = tm.makeTimeSeries()
        result = A.corr(B)
        expected, _ = stats.pearsonr(A, B)
        self.assertAlmostEqual(result, expected)

    def test_corr_rank(self):
        tm._skip_if_no_scipy()

        import scipy
        import scipy.stats as stats

        # kendall and spearman
        A = tm.makeTimeSeries()
        B = tm.makeTimeSeries()
        A[-5:] = A[:5]
        result = A.corr(B, method='kendall')
        expected = stats.kendalltau(A, B)[0]
        self.assertAlmostEqual(result, expected)

        result = A.corr(B, method='spearman')
        expected = stats.spearmanr(A, B)[0]
        self.assertAlmostEqual(result, expected)

        # these methods got rewritten in 0.8
        if scipy.__version__ < LooseVersion('0.9'):
            raise nose.SkipTest("skipping corr rank because of scipy version "
                                "{0}".format(scipy.__version__))

        # results from R
        A = Series([-0.89926396, 0.94209606, -1.03289164, -0.95445587,
                    0.76910310, -0.06430576, -2.09704447, 0.40660407,
                    -0.89926396, 0.94209606])
        B = Series([-1.01270225, -0.62210117, -1.56895827, 0.59592943,
                    -0.01680292, 1.17258718, -1.06009347, -0.10222060,
                    -0.89076239, 0.89372375])
        kexp = 0.4319297
        sexp = 0.5853767
        self.assertAlmostEqual(A.corr(B, method='kendall'), kexp)
        self.assertAlmostEqual(A.corr(B, method='spearman'), sexp)

    def test_cov(self):
        # full overlap
        self.assertAlmostEqual(self.ts.cov(self.ts), self.ts.std() ** 2)

        # partial overlap
        self.assertAlmostEqual(
            self.ts[:15].cov(self.ts[5:]), self.ts[5:15].std() ** 2)

        # No overlap
        self.assertTrue(np.isnan(self.ts[::2].cov(self.ts[1::2])))

        # all NA
        cp = self.ts[:10].copy()
        cp[:] = np.nan
        self.assertTrue(isnull(cp.cov(cp)))

        # min_periods
        self.assertTrue(isnull(self.ts[:15].cov(self.ts[5:], min_periods=12)))

        ts1 = self.ts[:15].reindex(self.ts.index)
        ts2 = self.ts[5:].reindex(self.ts.index)
        self.assertTrue(isnull(ts1.cov(ts2, min_periods=12)))

    def test_copy(self):
        ts = self.ts.copy()

        ts[::2] = np.NaN

        # Did not modify original Series
        self.assertFalse(np.isnan(self.ts[0]))

    def test_count(self):
        self.assertEqual(self.ts.count(), len(self.ts))

        self.ts[::2] = np.NaN

        self.assertEqual(self.ts.count(), np.isfinite(self.ts).sum())

        mi = MultiIndex.from_arrays([list('aabbcc'), [1, 2, 2, nan, 1, 2]])
        ts = Series(np.arange(len(mi)), index=mi)

        left = ts.count(level=1)
        right = Series([2, 3, 1], index=[1, 2, nan])
        assert_series_equal(left, right)

        ts.iloc[[0, 3, 5]] = nan
        assert_series_equal(ts.count(level=1), right - 1)

    def test_dtype(self):

        self.assertEqual(self.ts.dtype, np.dtype('float64'))
        self.assertEqual(self.ts.dtypes, np.dtype('float64'))
        self.assertEqual(self.ts.ftype, 'float64:dense')
        self.assertEqual(self.ts.ftypes, 'float64:dense')
        assert_series_equal(self.ts.get_dtype_counts(),Series(1,['float64']))
        assert_series_equal(self.ts.get_ftype_counts(),Series(1,['float64:dense']))

    def test_dot(self):
        a = Series(np.random.randn(4), index=['p', 'q', 'r', 's'])
        b = DataFrame(np.random.randn(3, 4), index=['1', '2', '3'],
                      columns=['p', 'q', 'r', 's']).T

        result = a.dot(b)
        expected = Series(np.dot(a.values, b.values),
                             index=['1', '2', '3'])
        assert_series_equal(result, expected)

        # Check index alignment
        b2 = b.reindex(index=reversed(b.index))
        result = a.dot(b)
        assert_series_equal(result, expected)

        # Check ndarray argument
        result = a.dot(b.values)
        self.assertTrue(np.all(result == expected.values))
        assert_almost_equal(a.dot(b['2'].values), expected['2'])

        # Check series argument
        assert_almost_equal(a.dot(b['1']), expected['1'])
        assert_almost_equal(a.dot(b2['1']), expected['1'])

        self.assertRaises(Exception, a.dot, a.values[:3])
        self.assertRaises(ValueError, a.dot, b.T)

    def test_value_counts_nunique(self):

        # basics.rst doc example
        series = Series(np.random.randn(500))
        series[20:500] = np.nan
        series[10:20]  = 5000
        result = series.nunique()
        self.assertEqual(result, 11)

    def test_unique(self):

        # 714 also, dtype=float
        s = Series([1.2345] * 100)
        s[::2] = np.nan
        result = s.unique()
        self.assertEqual(len(result), 2)

        s = Series([1.2345] * 100, dtype='f4')
        s[::2] = np.nan
        result = s.unique()
        self.assertEqual(len(result), 2)

        # NAs in object arrays #714
        s = Series(['foo'] * 100, dtype='O')
        s[::2] = np.nan
        result = s.unique()
        self.assertEqual(len(result), 2)

        # decision about None
        s = Series([1, 2, 3, None, None, None], dtype=object)
        result = s.unique()
        expected = np.array([1, 2, 3, None], dtype=object)
        self.assert_numpy_array_equal(result, expected)

    def test_dropna_empty(self):
        s = Series([])
        self.assertEqual(len(s.dropna()), 0)
        s.dropna(inplace=True)
        self.assertEqual(len(s), 0)

        # invalid axis
        self.assertRaises(ValueError, s.dropna, axis=1)

    def test_datetime64_tz_dropna(self):
        # DatetimeBlock
        s = Series([Timestamp('2011-01-01 10:00'), pd.NaT,
                    Timestamp('2011-01-03 10:00'), pd.NaT])
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
        for s in [Series([1, 2, 3], name='x'),
                  Series([False, True, False], name='x')]:

            result = s.dropna()
            self.assert_series_equal(result, s)
            self.assertFalse(result is s)

            s2 = s.copy()
            s2.dropna(inplace=True)
            self.assert_series_equal(s2, s)

    def test_axis_alias(self):
        s = Series([1, 2, np.nan])
        assert_series_equal(s.dropna(axis='rows'), s.dropna(axis='index'))
        self.assertEqual(s.dropna().sum('rows'), 3)
        self.assertEqual(s._get_axis_number('rows'), 0)
        self.assertEqual(s._get_axis_name('rows'), 'index')

    def test_drop_duplicates(self):
        # check both int and object
        for s in [Series([1, 2, 3, 3]), Series(['1', '2', '3', '3'])]:
            expected = Series([False, False, False, True])
            assert_series_equal(s.duplicated(), expected)
            assert_series_equal(s.drop_duplicates(), s[~expected])
            sc = s.copy()
            sc.drop_duplicates(inplace=True)
            assert_series_equal(sc, s[~expected])

            expected = Series([False, False, True, False])
            assert_series_equal(s.duplicated(keep='last'), expected)
            assert_series_equal(s.drop_duplicates(keep='last'), s[~expected])
            sc = s.copy()
            sc.drop_duplicates(keep='last', inplace=True)
            assert_series_equal(sc, s[~expected])

            # deprecate take_last
            with tm.assert_produces_warning(FutureWarning):
                assert_series_equal(s.duplicated(take_last=True), expected)
            with tm.assert_produces_warning(FutureWarning):
                assert_series_equal(s.drop_duplicates(take_last=True), s[~expected])
            sc = s.copy()
            with tm.assert_produces_warning(FutureWarning):
                sc.drop_duplicates(take_last=True, inplace=True)
            assert_series_equal(sc, s[~expected])

            expected = Series([False, False, True, True])
            assert_series_equal(s.duplicated(keep=False), expected)
            assert_series_equal(s.drop_duplicates(keep=False), s[~expected])
            sc = s.copy()
            sc.drop_duplicates(keep=False, inplace=True)
            assert_series_equal(sc, s[~expected])

        for s in [Series([1, 2, 3, 5, 3, 2, 4]),
                  Series(['1', '2', '3', '5', '3', '2', '4'])]:
            expected = Series([False, False, False, False, True, True, False])
            assert_series_equal(s.duplicated(), expected)
            assert_series_equal(s.drop_duplicates(), s[~expected])
            sc = s.copy()
            sc.drop_duplicates(inplace=True)
            assert_series_equal(sc, s[~expected])

            expected = Series([False, True, True, False, False, False, False])
            assert_series_equal(s.duplicated(keep='last'), expected)
            assert_series_equal(s.drop_duplicates(keep='last'), s[~expected])
            sc = s.copy()
            sc.drop_duplicates(keep='last', inplace=True)
            assert_series_equal(sc, s[~expected])

            # deprecate take_last
            with tm.assert_produces_warning(FutureWarning):
                assert_series_equal(s.duplicated(take_last=True), expected)
            with tm.assert_produces_warning(FutureWarning):
                assert_series_equal(s.drop_duplicates(take_last=True), s[~expected])
            sc = s.copy()
            with tm.assert_produces_warning(FutureWarning):
                sc.drop_duplicates(take_last=True, inplace=True)
            assert_series_equal(sc, s[~expected])

            expected = Series([False, True, True, False, True, True, False])
            assert_series_equal(s.duplicated(keep=False), expected)
            assert_series_equal(s.drop_duplicates(keep=False), s[~expected])
            sc = s.copy()
            sc.drop_duplicates(keep=False, inplace=True)
            assert_series_equal(sc, s[~expected])

    def test_sort_values(self):

        ts = self.ts.copy()

        # 9816 deprecated
        with tm.assert_produces_warning(FutureWarning):
            ts.sort()

        self.assert_numpy_array_equal(ts, self.ts.sort_values())
        self.assert_numpy_array_equal(ts.index, self.ts.sort_values().index)

        ts.sort_values(ascending=False, inplace=True)
        self.assert_numpy_array_equal(ts, self.ts.sort_values(ascending=False))
        self.assert_numpy_array_equal(ts.index,
                                      self.ts.sort_values(ascending=False).index)

        # GH 5856/5853
        # Series.sort_values operating on a view
        df = DataFrame(np.random.randn(10,4))
        s = df.iloc[:,0]
        def f():
            s.sort_values(inplace=True)
        self.assertRaises(ValueError, f)

        # test order/sort inplace
        # GH6859
        ts1 = self.ts.copy()
        ts1.sort_values(ascending=False, inplace=True)
        ts2 = self.ts.copy()
        ts2.sort_values(ascending=False, inplace=True)
        assert_series_equal(ts1,ts2)

        ts1 = self.ts.copy()
        ts1 = ts1.sort_values(ascending=False, inplace=False)
        ts2 = self.ts.copy()
        ts2 = ts.sort_values(ascending=False)
        assert_series_equal(ts1,ts2)

    def test_sort_index(self):
        rindex = list(self.ts.index)
        random.shuffle(rindex)

        random_order = self.ts.reindex(rindex)
        sorted_series = random_order.sort_index()
        assert_series_equal(sorted_series, self.ts)

        # descending
        sorted_series = random_order.sort_index(ascending=False)
        assert_series_equal(sorted_series,
                            self.ts.reindex(self.ts.index[::-1]))

    def test_sort_index_inplace(self):

        # For #11402
        rindex = list(self.ts.index)
        random.shuffle(rindex)

        # descending
        random_order = self.ts.reindex(rindex)
        result = random_order.sort_index(ascending=False, inplace=True)
        self.assertIs(result, None,
                      msg='sort_index() inplace should return None')
        assert_series_equal(random_order,
                            self.ts.reindex(self.ts.index[::-1]))

        # ascending
        random_order = self.ts.reindex(rindex)
        result = random_order.sort_index(ascending=True, inplace=True)
        self.assertIs(result, None,
                      msg='sort_index() inplace should return None')
        assert_series_equal(random_order, self.ts)

    def test_sort_API(self):

        # API for 9816

        # sortlevel
        mi = MultiIndex.from_tuples([[1, 1, 3], [1, 1, 1]], names=list('ABC'))
        s = Series([1, 2], mi)
        backwards = s.iloc[[1, 0]]

        res = s.sort_index(level='A')
        assert_series_equal(backwards, res)

        # sort_index
        rindex = list(self.ts.index)
        random.shuffle(rindex)

        random_order = self.ts.reindex(rindex)
        sorted_series = random_order.sort_index(level=0)
        assert_series_equal(sorted_series, self.ts)

        # compat on axis
        sorted_series = random_order.sort_index(axis=0)
        assert_series_equal(sorted_series, self.ts)

        self.assertRaises(ValueError, lambda : random_order.sort_values(axis=1))

        sorted_series = random_order.sort_index(level=0, axis=0)
        assert_series_equal(sorted_series, self.ts)

        self.assertRaises(ValueError, lambda : random_order.sort_index(level=0, axis=1))

    def test_order(self):

        # 9816 deprecated
        with tm.assert_produces_warning(FutureWarning):
            self.ts.order()

        ts = self.ts.copy()
        ts[:5] = np.NaN
        vals = ts.values

        result = ts.sort_values()
        self.assertTrue(np.isnan(result[-5:]).all())
        self.assert_numpy_array_equal(result[:-5], np.sort(vals[5:]))

        result = ts.sort_values(na_position='first')
        self.assertTrue(np.isnan(result[:5]).all())
        self.assert_numpy_array_equal(result[5:], np.sort(vals[5:]))

        # something object-type
        ser = Series(['A', 'B'], [1, 2])
        # no failure
        ser.sort_values()

        # ascending=False
        ordered = ts.sort_values(ascending=False)
        expected = np.sort(ts.valid().values)[::-1]
        assert_almost_equal(expected, ordered.valid().values)
        ordered = ts.sort_values(ascending=False, na_position='first')
        assert_almost_equal(expected, ordered.valid().values)

    def test_nsmallest_nlargest(self):
        # float, int, datetime64 (use i8), timedelts64 (same),
        # object that are numbers, object that are strings

        base = [3, 2, 1, 2, 5]

        s_list = [
            Series(base, dtype='int8'),
            Series(base, dtype='int16'),
            Series(base, dtype='int32'),
            Series(base, dtype='int64'),
            Series(base, dtype='float32'),
            Series(base, dtype='float64'),
            Series(base, dtype='uint8'),
            Series(base, dtype='uint16'),
            Series(base, dtype='uint32'),
            Series(base, dtype='uint64'),
            Series(base).astype('timedelta64[ns]'),
            Series(pd.to_datetime(['2003', '2002', '2001', '2002', '2005'])),
        ]

        raising = [
            Series([3., 2, 1, 2, '5'], dtype='object'),
            Series([3., 2, 1, 2, 5], dtype='object'),
            # not supported on some archs
            # Series([3., 2, 1, 2, 5], dtype='complex256'),
            Series([3., 2, 1, 2, 5], dtype='complex128'),
        ]

        for r in raising:
            dt = r.dtype
            msg = "Cannot use method 'n(larg|small)est' with dtype %s" % dt
            args = 2, len(r), 0, -1
            methods = r.nlargest, r.nsmallest
            for method, arg in product(methods, args):
                with tm.assertRaisesRegexp(TypeError, msg):
                    method(arg)

        for s in s_list:

            assert_series_equal(s.nsmallest(2), s.iloc[[2, 1]])

            assert_series_equal(s.nsmallest(2, keep='last'), s.iloc[[2, 3]])
            with tm.assert_produces_warning(FutureWarning):
                assert_series_equal(s.nsmallest(2, take_last=True), s.iloc[[2, 3]])

            assert_series_equal(s.nlargest(3), s.iloc[[4, 0, 1]])

            assert_series_equal(s.nlargest(3, keep='last'), s.iloc[[4, 0, 3]])
            with tm.assert_produces_warning(FutureWarning):
                assert_series_equal(s.nlargest(3, take_last=True), s.iloc[[4, 0, 3]])

            empty = s.iloc[0:0]
            assert_series_equal(s.nsmallest(0), empty)
            assert_series_equal(s.nsmallest(-1), empty)
            assert_series_equal(s.nlargest(0), empty)
            assert_series_equal(s.nlargest(-1), empty)

            assert_series_equal(s.nsmallest(len(s)), s.sort_values())
            assert_series_equal(s.nsmallest(len(s) + 1), s.sort_values())
            assert_series_equal(s.nlargest(len(s)), s.iloc[[4, 0, 1, 3, 2]])
            assert_series_equal(s.nlargest(len(s) + 1),
                                s.iloc[[4, 0, 1, 3, 2]])

        s = Series([3., np.nan, 1, 2, 5])
        assert_series_equal(s.nlargest(), s.iloc[[4, 0, 3, 2]])
        assert_series_equal(s.nsmallest(), s.iloc[[2, 3, 0, 4]])

        msg = 'keep must be either "first", "last"'
        with tm.assertRaisesRegexp(ValueError, msg):
            s.nsmallest(keep='invalid')
        with tm.assertRaisesRegexp(ValueError, msg):
            s.nlargest(keep='invalid')

    def test_rank(self):
        tm._skip_if_no_scipy()
        from scipy.stats import rankdata

        self.ts[::2] = np.nan
        self.ts[:10][::3] = 4.

        ranks = self.ts.rank()
        oranks = self.ts.astype('O').rank()

        assert_series_equal(ranks, oranks)

        mask = np.isnan(self.ts)
        filled = self.ts.fillna(np.inf)

        # rankdata returns a ndarray
        exp = Series(rankdata(filled),index=filled.index)
        exp[mask] = np.nan

        assert_almost_equal(ranks, exp)

        iseries = Series(np.arange(5).repeat(2))

        iranks = iseries.rank()
        exp = iseries.astype(float).rank()
        assert_series_equal(iranks, exp)
        iseries = Series(np.arange(5)) + 1.0
        exp = iseries / 5.0
        iranks = iseries.rank(pct=True)

        assert_series_equal(iranks, exp)

        iseries = Series(np.repeat(1, 100))
        exp = Series(np.repeat(0.505, 100))
        iranks = iseries.rank(pct=True)
        assert_series_equal(iranks, exp)

        iseries[1] = np.nan
        exp = Series(np.repeat(50.0 / 99.0, 100))
        exp[1] =  np.nan
        iranks = iseries.rank(pct=True)
        assert_series_equal(iranks, exp)

        iseries = Series(np.arange(5)) + 1.0
        iseries[4] = np.nan
        exp = iseries / 4.0
        iranks = iseries.rank(pct=True)
        assert_series_equal(iranks, exp)

        iseries = Series(np.repeat(np.nan, 100))
        exp = iseries.copy()
        iranks = iseries.rank(pct=True)
        assert_series_equal(iranks, exp)

        iseries = Series(np.arange(5)) + 1
        iseries[4] = np.nan
        exp = iseries / 4.0
        iranks = iseries.rank(pct=True)
        assert_series_equal(iranks, exp)

        rng = date_range('1/1/1990', periods=5)
        iseries = Series(np.arange(5), rng) + 1
        iseries.ix[4] = np.nan
        exp = iseries / 4.0
        iranks = iseries.rank(pct=True)
        assert_series_equal(iranks, exp)

        iseries = Series([1e-50, 1e-100, 1e-20, 1e-2, 1e-20+1e-30, 1e-1])
        exp = Series([2, 1, 3, 5, 4, 6.0])
        iranks = iseries.rank()
        assert_series_equal(iranks, exp)

        values = np.array([-50, -1, -1e-20, -1e-25, -1e-50, 0, 1e-40, 1e-20, 1e-10, 2, 40], dtype='float64')
        random_order = np.random.permutation(len(values))
        iseries = Series(values[random_order])
        exp = Series(random_order + 1.0, dtype='float64')
        iranks = iseries.rank()
        assert_series_equal(iranks, exp)

    def test_rank_inf(self):
        raise nose.SkipTest('DataFrame.rank does not currently rank np.inf and -np.inf properly')

        values = np.array([-np.inf, -50, -1, -1e-20, -1e-25, -1e-50, 0, 1e-40, 1e-20, 1e-10, 2, 40, np.inf], dtype='float64')
        random_order = np.random.permutation(len(values))
        iseries = Series(values[random_order])
        exp = Series(random_order + 1.0, dtype='float64')
        iranks = iseries.rank()
        assert_series_equal(iranks, exp)


    def test_from_csv(self):

        with ensure_clean() as path:
            self.ts.to_csv(path)
            ts = Series.from_csv(path)
            assert_series_equal(self.ts, ts, check_names=False)
            self.assertTrue(ts.name is None)
            self.assertTrue(ts.index.name is None)

            # GH10483
            self.ts.to_csv(path, header=True)
            ts_h = Series.from_csv(path, header=0)
            self.assertTrue(ts_h.name == 'ts')

            self.series.to_csv(path)
            series = Series.from_csv(path)
            self.assertIsNone(series.name)
            self.assertIsNone(series.index.name)
            assert_series_equal(self.series, series, check_names=False)
            self.assertTrue(series.name is None)
            self.assertTrue(series.index.name is None)

            self.series.to_csv(path, header=True)
            series_h = Series.from_csv(path, header=0)
            self.assertTrue(series_h.name == 'series')

            outfile = open(path, 'w')
            outfile.write('1998-01-01|1.0\n1999-01-01|2.0')
            outfile.close()
            series = Series.from_csv(path, sep='|')
            checkseries = Series(
                 {datetime(1998, 1, 1): 1.0, datetime(1999, 1, 1): 2.0})
            assert_series_equal(checkseries, series)

            series = Series.from_csv(path, sep='|', parse_dates=False)
            checkseries = Series({'1998-01-01': 1.0, '1999-01-01': 2.0})
            assert_series_equal(checkseries, series)

    def test_to_csv(self):
        import io

        with ensure_clean() as path:
            self.ts.to_csv(path)

            lines = io.open(path, newline=None).readlines()
            assert(lines[1] != '\n')

            self.ts.to_csv(path, index=False)
            arr = np.loadtxt(path)
            assert_almost_equal(arr, self.ts.values)

    def test_to_csv_unicode_index(self):
        buf = StringIO()
        s = Series([u("\u05d0"), "d2"], index=[u("\u05d0"), u("\u05d1")])

        s.to_csv(buf, encoding='UTF-8')
        buf.seek(0)

        s2 = Series.from_csv(buf, index_col=0, encoding='UTF-8')

        assert_series_equal(s, s2)

    def test_tolist(self):
        rs = self.ts.tolist()
        xp = self.ts.values.tolist()
        assert_almost_equal(rs, xp)

        # datetime64
        s = Series(self.ts.index)
        rs = s.tolist()
        self.assertEqual(self.ts.index[0], rs[0])

    def test_to_frame(self):
        self.ts.name = None
        rs = self.ts.to_frame()
        xp = pd.DataFrame(self.ts.values, index=self.ts.index)
        assert_frame_equal(rs, xp)

        self.ts.name = 'testname'
        rs = self.ts.to_frame()
        xp = pd.DataFrame(dict(testname=self.ts.values), index=self.ts.index)
        assert_frame_equal(rs, xp)

        rs = self.ts.to_frame(name='testdifferent')
        xp = pd.DataFrame(dict(testdifferent=self.ts.values), index=self.ts.index)
        assert_frame_equal(rs, xp)

    def test_to_dict(self):
        self.assert_numpy_array_equal(Series(self.ts.to_dict()), self.ts)

    def test_to_csv_float_format(self):

        with ensure_clean() as filename:
            ser = Series([0.123456, 0.234567, 0.567567])
            ser.to_csv(filename, float_format='%.2f')

            rs = Series.from_csv(filename)
            xp = Series([0.12, 0.23, 0.57])
            assert_series_equal(rs, xp)

    def test_to_csv_list_entries(self):
        s = Series(['jack and jill', 'jesse and frank'])

        split = s.str.split(r'\s+and\s+')

        buf = StringIO()
        split.to_csv(buf)

    def test_to_csv_path_is_none(self):
        # GH 8215
        # Series.to_csv() was returning None, inconsistent with
        # DataFrame.to_csv() which returned string
        s = Series([1, 2, 3])
        csv_str = s.to_csv(path=None)
        self.assertIsInstance(csv_str, str)

    def test_str_attribute(self):
        # GH9068
        methods = ['strip', 'rstrip', 'lstrip']
        s = Series([' jack', 'jill ', ' jesse ', 'frank'])
        for method in methods:
            expected = Series([getattr(str, method)(x) for x in s.values])
            assert_series_equal(getattr(Series.str, method)(s.str), expected)

        # str accessor only valid with string values
        s = Series(range(5))
        with self.assertRaisesRegexp(AttributeError, 'only use .str accessor'):
            s.str.repeat(2)

    def test_clip(self):
        val = self.ts.median()

        self.assertEqual(self.ts.clip_lower(val).min(), val)
        self.assertEqual(self.ts.clip_upper(val).max(), val)

        self.assertEqual(self.ts.clip(lower=val).min(), val)
        self.assertEqual(self.ts.clip(upper=val).max(), val)

        result = self.ts.clip(-0.5, 0.5)
        expected = np.clip(self.ts, -0.5, 0.5)
        assert_series_equal(result, expected)
        tm.assertIsInstance(expected, Series)

    def test_clip_types_and_nulls(self):

        sers = [Series([np.nan, 1.0, 2.0, 3.0]),
                Series([None, 'a', 'b', 'c']),
                Series(pd.to_datetime([np.nan, 1, 2, 3], unit='D'))]

        for s in sers:
            thresh = s[2]
            l = s.clip_lower(thresh)
            u = s.clip_upper(thresh)
            self.assertEqual(l[notnull(l)].min(), thresh)
            self.assertEqual(u[notnull(u)].max(), thresh)
            self.assertEqual(list(isnull(s)), list(isnull(l)))
            self.assertEqual(list(isnull(s)), list(isnull(u)))

    def test_clip_against_series(self):
        # GH #6966

        s = Series([1.0, 1.0, 4.0])
        threshold = Series([1.0, 2.0, 3.0])

        assert_series_equal(s.clip_lower(threshold), Series([1.0, 2.0, 4.0]))
        assert_series_equal(s.clip_upper(threshold), Series([1.0, 1.0, 3.0]))

        lower = Series([1.0, 2.0, 3.0])
        upper = Series([1.5, 2.5, 3.5])
        assert_series_equal(s.clip(lower, upper), Series([1.0, 2.0, 3.5]))
        assert_series_equal(s.clip(1.5, upper), Series([1.5, 1.5, 3.5]))

    def test_valid(self):
        ts = self.ts.copy()
        ts[::2] = np.NaN

        result = ts.valid()
        self.assertEqual(len(result), ts.count())

        tm.assert_dict_equal(result, ts, compare_keys=False)

    def test_isnull(self):
        ser = Series([0, 5.4, 3, nan, -0.001])
        np.array_equal(
            ser.isnull(), Series([False, False, False, True, False]).values)
        ser = Series(["hi", "", nan])
        np.array_equal(ser.isnull(), Series([False, False, True]).values)

    def test_notnull(self):
        ser = Series([0, 5.4, 3, nan, -0.001])
        np.array_equal(
            ser.notnull(), Series([True, True, True, False, True]).values)
        ser = Series(["hi", "", nan])
        np.array_equal(ser.notnull(), Series([True, True, False]).values)

    def test_shift(self):
        shifted = self.ts.shift(1)
        unshifted = shifted.shift(-1)

        tm.assert_dict_equal(unshifted.valid(), self.ts, compare_keys=False)

        offset = datetools.bday
        shifted = self.ts.shift(1, freq=offset)
        unshifted = shifted.shift(-1, freq=offset)

        assert_series_equal(unshifted, self.ts)

        unshifted = self.ts.shift(0, freq=offset)
        assert_series_equal(unshifted, self.ts)

        shifted = self.ts.shift(1, freq='B')
        unshifted = shifted.shift(-1, freq='B')

        assert_series_equal(unshifted, self.ts)

        # corner case
        unshifted = self.ts.shift(0)
        assert_series_equal(unshifted, self.ts)

        # Shifting with PeriodIndex
        ps = tm.makePeriodSeries()
        shifted = ps.shift(1)
        unshifted = shifted.shift(-1)
        tm.assert_dict_equal(unshifted.valid(), ps, compare_keys=False)

        shifted2 = ps.shift(1, 'B')
        shifted3 = ps.shift(1, datetools.bday)
        assert_series_equal(shifted2, shifted3)
        assert_series_equal(ps, shifted2.shift(-1, 'B'))

        self.assertRaises(ValueError, ps.shift, freq='D')

        # legacy support
        shifted4 = ps.shift(1, freq='B')
        assert_series_equal(shifted2, shifted4)

        shifted5 = ps.shift(1, freq=datetools.bday)
        assert_series_equal(shifted5, shifted4)

        # 32-bit taking
        # GH 8129
        index=date_range('2000-01-01',periods=5)
        for dtype in ['int32','int64']:
            s1 = Series(np.arange(5,dtype=dtype),index=index)
            p = s1.iloc[1]
            result = s1.shift(periods=p)
            expected = Series([np.nan,0,1,2,3],index=index)
            assert_series_equal(result,expected)

        # xref 8260
        # with tz
        s = Series(date_range('2000-01-01 09:00:00',periods=5,tz='US/Eastern'),name='foo')
        result = s-s.shift()
        assert_series_equal(result,Series(TimedeltaIndex(['NaT'] + ['1 days']*4),name='foo'))

        # incompat tz
        s2 = Series(date_range('2000-01-01 09:00:00',periods=5,tz='CET'),name='foo')
        self.assertRaises(ValueError, lambda : s-s2)

    def test_tshift(self):
        # PeriodIndex
        ps = tm.makePeriodSeries()
        shifted = ps.tshift(1)
        unshifted = shifted.tshift(-1)

        assert_series_equal(unshifted, ps)

        shifted2 = ps.tshift(freq='B')
        assert_series_equal(shifted, shifted2)

        shifted3 = ps.tshift(freq=datetools.bday)
        assert_series_equal(shifted, shifted3)

        self.assertRaises(ValueError, ps.tshift, freq='M')

        # DatetimeIndex
        shifted = self.ts.tshift(1)
        unshifted = shifted.tshift(-1)

        assert_series_equal(self.ts, unshifted)

        shifted2 = self.ts.tshift(freq=self.ts.index.freq)
        assert_series_equal(shifted, shifted2)

        inferred_ts = Series(self.ts.values, Index(np.asarray(self.ts.index)),
                             name='ts')
        shifted = inferred_ts.tshift(1)
        unshifted = shifted.tshift(-1)
        assert_series_equal(shifted, self.ts.tshift(1))
        assert_series_equal(unshifted, inferred_ts)

        no_freq = self.ts[[0, 5, 7]]
        self.assertRaises(ValueError, no_freq.tshift)

    def test_shift_int(self):
        ts = self.ts.astype(int)
        shifted = ts.shift(1)
        expected = ts.astype(float).shift(1)
        assert_series_equal(shifted, expected)

    def test_shift_categorical(self):
        # GH 9416
        s = pd.Series(['a', 'b', 'c', 'd'], dtype='category')

        assert_series_equal(s.iloc[:-1], s.shift(1).shift(-1).valid())

        sp1 = s.shift(1)
        assert_index_equal(s.index, sp1.index)
        self.assertTrue(np.all(sp1.values.codes[:1] == -1))
        self.assertTrue(np.all(s.values.codes[:-1] == sp1.values.codes[1:]))

        sn2 = s.shift(-2)
        assert_index_equal(s.index, sn2.index)
        self.assertTrue(np.all(sn2.values.codes[-2:] == -1))
        self.assertTrue(np.all(s.values.codes[2:] == sn2.values.codes[:-2]))

        assert_index_equal(s.values.categories, sp1.values.categories)
        assert_index_equal(s.values.categories, sn2.values.categories)

    def test_truncate(self):
        offset = datetools.bday

        ts = self.ts[::3]

        start, end = self.ts.index[3], self.ts.index[6]
        start_missing, end_missing = self.ts.index[2], self.ts.index[7]

        # neither specified
        truncated = ts.truncate()
        assert_series_equal(truncated, ts)

        # both specified
        expected = ts[1:3]

        truncated = ts.truncate(start, end)
        assert_series_equal(truncated, expected)

        truncated = ts.truncate(start_missing, end_missing)
        assert_series_equal(truncated, expected)

        # start specified
        expected = ts[1:]

        truncated = ts.truncate(before=start)
        assert_series_equal(truncated, expected)

        truncated = ts.truncate(before=start_missing)
        assert_series_equal(truncated, expected)

        # end specified
        expected = ts[:3]

        truncated = ts.truncate(after=end)
        assert_series_equal(truncated, expected)

        truncated = ts.truncate(after=end_missing)
        assert_series_equal(truncated, expected)

        # corner case, empty series returned
        truncated = ts.truncate(after=self.ts.index[0] - offset)
        assert(len(truncated) == 0)

        truncated = ts.truncate(before=self.ts.index[-1] + offset)
        assert(len(truncated) == 0)

        self.assertRaises(ValueError, ts.truncate,
                          before=self.ts.index[-1] + offset,
                          after=self.ts.index[0] - offset)

    def test_ptp(self):
        N = 1000
        arr = np.random.randn(N)
        ser = Series(arr)
        self.assertEqual(np.ptp(ser), np.ptp(arr))

        # GH11163
        s = Series([3, 5, np.nan, -3, 10])
        self.assertEqual(s.ptp(), 13)
        self.assertTrue(pd.isnull(s.ptp(skipna=False)))

        mi = pd.MultiIndex.from_product([['a','b'], [1,2,3]])
        s = pd.Series([1, np.nan, 7, 3, 5, np.nan], index=mi)

        expected = pd.Series([6, 2], index=['a', 'b'], dtype=np.float64)
        self.assert_series_equal(s.ptp(level=0), expected)

        expected = pd.Series([np.nan, np.nan], index=['a', 'b'])
        self.assert_series_equal(s.ptp(level=0, skipna=False), expected)

        with self.assertRaises(ValueError):
            s.ptp(axis=1)

        s = pd.Series(['a', 'b', 'c', 'd', 'e'])
        with self.assertRaises(TypeError):
            s.ptp()

        with self.assertRaises(NotImplementedError):
            s.ptp(numeric_only=True)


    def test_asof(self):
        # array or list or dates
        N = 50
        rng = date_range('1/1/1990', periods=N, freq='53s')
        ts = Series(np.random.randn(N), index=rng)
        ts[15:30] = np.nan
        dates = date_range('1/1/1990', periods=N * 3, freq='25s')

        result = ts.asof(dates)
        self.assertTrue(notnull(result).all())
        lb = ts.index[14]
        ub = ts.index[30]

        result = ts.asof(list(dates))
        self.assertTrue(notnull(result).all())
        lb = ts.index[14]
        ub = ts.index[30]

        mask = (result.index >= lb) & (result.index < ub)
        rs = result[mask]
        self.assertTrue((rs == ts[lb]).all())

        val = result[result.index[result.index >= ub][0]]
        self.assertEqual(ts[ub], val)

        self.ts[5:10] = np.NaN
        self.ts[15:20] = np.NaN

        val1 = self.ts.asof(self.ts.index[7])
        val2 = self.ts.asof(self.ts.index[19])

        self.assertEqual(val1, self.ts[4])
        self.assertEqual(val2, self.ts[14])

        # accepts strings
        val1 = self.ts.asof(str(self.ts.index[7]))
        self.assertEqual(val1, self.ts[4])

        # in there
        self.assertEqual(self.ts.asof(self.ts.index[3]), self.ts[3])

        # no as of value
        d = self.ts.index[0] - datetools.bday
        self.assertTrue(np.isnan(self.ts.asof(d)))

    def test_getitem_setitem_datetimeindex(self):
        from pandas import date_range
        N = 50
        # testing with timezone, GH #2785
        rng = date_range('1/1/1990', periods=N, freq='H', tz='US/Eastern')
        ts = Series(np.random.randn(N), index=rng)

        result = ts["1990-01-01 04:00:00"]
        expected = ts[4]
        self.assertEqual(result, expected)

        result = ts.copy()
        result["1990-01-01 04:00:00"] = 0
        result["1990-01-01 04:00:00"] = ts[4]
        assert_series_equal(result, ts)

        result = ts["1990-01-01 04:00:00":"1990-01-01 07:00:00"]
        expected = ts[4:8]
        assert_series_equal(result, expected)

        result = ts.copy()
        result["1990-01-01 04:00:00":"1990-01-01 07:00:00"] = 0
        result["1990-01-01 04:00:00":"1990-01-01 07:00:00"] = ts[4:8]
        assert_series_equal(result, ts)

        lb = "1990-01-01 04:00:00"
        rb = "1990-01-01 07:00:00"
        result = ts[(ts.index >= lb) & (ts.index <= rb)]
        expected = ts[4:8]
        assert_series_equal(result, expected)

        # repeat all the above with naive datetimes
        result = ts[datetime(1990, 1, 1, 4)]
        expected = ts[4]
        self.assertEqual(result, expected)

        result = ts.copy()
        result[datetime(1990, 1, 1, 4)] = 0
        result[datetime(1990, 1, 1, 4)] = ts[4]
        assert_series_equal(result, ts)

        result = ts[datetime(1990, 1, 1, 4):datetime(1990, 1, 1, 7)]
        expected = ts[4:8]
        assert_series_equal(result, expected)

        result = ts.copy()
        result[datetime(1990, 1, 1, 4):datetime(1990, 1, 1, 7)] = 0
        result[datetime(1990, 1, 1, 4):datetime(1990, 1, 1, 7)] = ts[4:8]
        assert_series_equal(result, ts)

        lb = datetime(1990, 1, 1, 4)
        rb = datetime(1990, 1, 1, 7)
        result = ts[(ts.index >= lb) & (ts.index <= rb)]
        expected = ts[4:8]
        assert_series_equal(result, expected)

        result = ts[ts.index[4]]
        expected = ts[4]
        self.assertEqual(result, expected)

        result = ts[ts.index[4:8]]
        expected = ts[4:8]
        assert_series_equal(result, expected)

        result = ts.copy()
        result[ts.index[4:8]] = 0
        result[4:8] = ts[4:8]
        assert_series_equal(result, ts)

        # also test partial date slicing
        result = ts["1990-01-02"]
        expected = ts[24:48]
        assert_series_equal(result, expected)

        result = ts.copy()
        result["1990-01-02"] = 0
        result["1990-01-02"] = ts[24:48]
        assert_series_equal(result, ts)

    def test_getitem_setitem_datetime_tz_pytz(self):
        tm._skip_if_no_pytz()
        from pytz import timezone as tz

        from pandas import date_range
        N = 50
        # testing with timezone, GH #2785
        rng = date_range('1/1/1990', periods=N, freq='H', tz='US/Eastern')
        ts = Series(np.random.randn(N), index=rng)

        # also test Timestamp tz handling, GH #2789
        result = ts.copy()
        result["1990-01-01 09:00:00+00:00"] = 0
        result["1990-01-01 09:00:00+00:00"] = ts[4]
        assert_series_equal(result, ts)

        result = ts.copy()
        result["1990-01-01 03:00:00-06:00"] = 0
        result["1990-01-01 03:00:00-06:00"] = ts[4]
        assert_series_equal(result, ts)

        # repeat with datetimes
        result = ts.copy()
        result[datetime(1990, 1, 1, 9, tzinfo=tz('UTC'))] = 0
        result[datetime(1990, 1, 1, 9, tzinfo=tz('UTC'))] = ts[4]
        assert_series_equal(result, ts)

        result = ts.copy()

        # comparison dates with datetime MUST be localized!
        date = tz('US/Central').localize(datetime(1990, 1, 1, 3))
        result[date] = 0
        result[date] = ts[4]
        assert_series_equal(result, ts)


    def test_getitem_setitem_datetime_tz_dateutil(self):
        tm._skip_if_no_dateutil()
        from dateutil.tz import tzutc
        from pandas.tslib import _dateutil_gettz as gettz

        tz = lambda x: tzutc() if x == 'UTC' else gettz(x)  # handle special case for utc in dateutil

        from pandas import date_range
        N = 50
        # testing with timezone, GH #2785
        rng = date_range('1/1/1990', periods=N, freq='H', tz='US/Eastern')
        ts = Series(np.random.randn(N), index=rng)

        # also test Timestamp tz handling, GH #2789
        result = ts.copy()
        result["1990-01-01 09:00:00+00:00"] = 0
        result["1990-01-01 09:00:00+00:00"] = ts[4]
        assert_series_equal(result, ts)

        result = ts.copy()
        result["1990-01-01 03:00:00-06:00"] = 0
        result["1990-01-01 03:00:00-06:00"] = ts[4]
        assert_series_equal(result, ts)

        # repeat with datetimes
        result = ts.copy()
        result[datetime(1990, 1, 1, 9, tzinfo=tz('UTC'))] = 0
        result[datetime(1990, 1, 1, 9, tzinfo=tz('UTC'))] = ts[4]
        assert_series_equal(result, ts)

        result = ts.copy()
        result[datetime(1990, 1, 1, 3, tzinfo=tz('US/Central'))] = 0
        result[datetime(1990, 1, 1, 3, tzinfo=tz('US/Central'))] = ts[4]
        assert_series_equal(result, ts)

    def test_getitem_setitem_periodindex(self):
        from pandas import period_range
        N = 50
        rng = period_range('1/1/1990', periods=N, freq='H')
        ts = Series(np.random.randn(N), index=rng)

        result = ts["1990-01-01 04"]
        expected = ts[4]
        self.assertEqual(result, expected)

        result = ts.copy()
        result["1990-01-01 04"] = 0
        result["1990-01-01 04"] = ts[4]
        assert_series_equal(result, ts)

        result = ts["1990-01-01 04":"1990-01-01 07"]
        expected = ts[4:8]
        assert_series_equal(result, expected)

        result = ts.copy()
        result["1990-01-01 04":"1990-01-01 07"] = 0
        result["1990-01-01 04":"1990-01-01 07"] = ts[4:8]
        assert_series_equal(result, ts)

        lb = "1990-01-01 04"
        rb = "1990-01-01 07"
        result = ts[(ts.index >= lb) & (ts.index <= rb)]
        expected = ts[4:8]
        assert_series_equal(result, expected)

        # GH 2782
        result = ts[ts.index[4]]
        expected = ts[4]
        self.assertEqual(result, expected)

        result = ts[ts.index[4:8]]
        expected = ts[4:8]
        assert_series_equal(result, expected)

        result = ts.copy()
        result[ts.index[4:8]] = 0
        result[4:8] = ts[4:8]
        assert_series_equal(result, ts)

    def test_asof_periodindex(self):
        from pandas import period_range, PeriodIndex
        # array or list or dates
        N = 50
        rng = period_range('1/1/1990', periods=N, freq='H')
        ts = Series(np.random.randn(N), index=rng)
        ts[15:30] = np.nan
        dates = date_range('1/1/1990', periods=N * 3, freq='37min')

        result = ts.asof(dates)
        self.assertTrue(notnull(result).all())
        lb = ts.index[14]
        ub = ts.index[30]

        result = ts.asof(list(dates))
        self.assertTrue(notnull(result).all())
        lb = ts.index[14]
        ub = ts.index[30]

        pix = PeriodIndex(result.index.values, freq='H')
        mask = (pix >= lb) & (pix < ub)
        rs = result[mask]
        self.assertTrue((rs == ts[lb]).all())

        ts[5:10] = np.NaN
        ts[15:20] = np.NaN

        val1 = ts.asof(ts.index[7])
        val2 = ts.asof(ts.index[19])

        self.assertEqual(val1, ts[4])
        self.assertEqual(val2, ts[14])

        # accepts strings
        val1 = ts.asof(str(ts.index[7]))
        self.assertEqual(val1, ts[4])

        # in there
        self.assertEqual(ts.asof(ts.index[3]), ts[3])

        # no as of value
        d = ts.index[0].to_timestamp() - datetools.bday
        self.assertTrue(np.isnan(ts.asof(d)))

    def test_asof_more(self):
        from pandas import date_range
        s = Series([nan, nan, 1, 2, nan, nan, 3, 4, 5],
                   index=date_range('1/1/2000', periods=9))

        dates = s.index[[4, 5, 6, 2, 1]]

        result = s.asof(dates)
        expected = Series([2, 2, 3, 1, np.nan], index=dates)

        assert_series_equal(result, expected)

        s = Series([1.5, 2.5, 1, 2, nan, nan, 3, 4, 5],
                   index=date_range('1/1/2000', periods=9))
        result = s.asof(s.index[0])
        self.assertEqual(result, s[0])

    def test_cast_on_putmask(self):

        # GH 2746

        # need to upcast
        s = Series([1, 2], index=[1, 2], dtype='int64')
        s[[True, False]] = Series([0], index=[1], dtype='int64')
        expected = Series([0, 2], index=[1, 2], dtype='int64')

        assert_series_equal(s, expected)

    def test_type_promote_putmask(self):

        # GH8387: test that changing types does not break alignment
        ts = Series(np.random.randn(100), index=np.arange(100,0,-1)).round(5)
        left, mask = ts.copy(), ts > 0
        right = ts[mask].copy().map(str)
        left[mask] = right
        assert_series_equal(left, ts.map(lambda t: str(t) if t > 0 else t))

        s = Series([0, 1, 2, 0 ])
        mask = s > 0
        s2 = s[ mask ].map( str )
        s[mask] = s2
        assert_series_equal(s, Series([0, '1', '2', 0]))

        s = Series([0, 'foo', 'bar', 0 ])
        mask = Series([False, True, True, False])
        s2 = s[ mask ]
        s[mask] = s2
        assert_series_equal(s, Series([0, 'foo','bar', 0]))

    def test_astype_cast_nan_int(self):
        df = Series([1.0, 2.0, 3.0, np.nan])
        self.assertRaises(ValueError, df.astype, np.int64)

    def test_astype_cast_object_int(self):
        arr = Series(["car", "house", "tree", "1"])

        self.assertRaises(ValueError, arr.astype, int)
        self.assertRaises(ValueError, arr.astype, np.int64)
        self.assertRaises(ValueError, arr.astype, np.int8)

        arr = Series(['1', '2', '3', '4'], dtype=object)
        result = arr.astype(int)
        self.assert_numpy_array_equal(result, np.arange(1, 5))

    def test_astype_datetimes(self):
        import pandas.tslib as tslib

        s = Series(tslib.iNaT, dtype='M8[ns]', index=lrange(5))
        s = s.astype('O')
        self.assertEqual(s.dtype, np.object_)

        s = Series([datetime(2001, 1, 2, 0, 0)])
        s = s.astype('O')
        self.assertEqual(s.dtype, np.object_)

        s = Series([datetime(2001, 1, 2, 0, 0) for i in range(3)])
        s[1] = np.nan
        self.assertEqual(s.dtype, 'M8[ns]')
        s = s.astype('O')
        self.assertEqual(s.dtype, np.object_)

    def test_astype_str(self):
        # GH4405
        digits = string.digits
        s1 = Series([digits * 10, tm.rands(63), tm.rands(64),
                    tm.rands(1000)])
        s2 = Series([digits * 10, tm.rands(63), tm.rands(64), nan, 1.0])
        types = (compat.text_type, np.str_)
        for typ in types:
            for s in (s1, s2):
                res = s.astype(typ)
                expec = s.map(compat.text_type)
                assert_series_equal(res, expec)

        # GH9757
        # Test str and unicode on python 2.x and just str on python 3.x
        for tt in set([str, compat.text_type]):
            ts = Series([Timestamp('2010-01-04 00:00:00')])
            s = ts.astype(tt)
            expected = Series([tt('2010-01-04')])
            assert_series_equal(s, expected)

            ts = Series([Timestamp('2010-01-04 00:00:00', tz='US/Eastern')])
            s = ts.astype(tt)
            expected = Series([tt('2010-01-04 00:00:00-05:00')])
            assert_series_equal(s, expected)

            td = Series([Timedelta(1, unit='d')])
            s = td.astype(tt)
            expected = Series([tt('1 days 00:00:00.000000000')])
            assert_series_equal(s, expected)

    def test_astype_unicode(self):

        # GH7758
        # a bit of magic is required to set default encoding encoding to utf-8
        digits = string.digits
        test_series = [
            Series([digits * 10, tm.rands(63), tm.rands(64), tm.rands(1000)]),
            Series([u('データーサイエンス、お前はもう死んでいる')]),

        ]

        former_encoding = None
        if not compat.PY3:
            # in python we can force the default encoding
            # for this test
            former_encoding = sys.getdefaultencoding()
            reload(sys)
            sys.setdefaultencoding("utf-8")
        if sys.getdefaultencoding() == "utf-8":
            test_series.append(Series([u('野菜食べないとやばい').encode("utf-8")]))
        for s in test_series:
            res = s.astype("unicode")
            expec = s.map(compat.text_type)
            assert_series_equal(res, expec)
        # restore the former encoding
        if former_encoding is not None and former_encoding != "utf-8":
            reload(sys)
            sys.setdefaultencoding(former_encoding)


    def test_map(self):
        index, data = tm.getMixedTypeDict()

        source = Series(data['B'], index=data['C'])
        target = Series(data['C'][:4], index=data['D'][:4])

        merged = target.map(source)

        for k, v in compat.iteritems(merged):
            self.assertEqual(v, source[target[k]])

        # input could be a dict
        merged = target.map(source.to_dict())

        for k, v in compat.iteritems(merged):
            self.assertEqual(v, source[target[k]])

        # function
        result = self.ts.map(lambda x: x * 2)
        self.assert_numpy_array_equal(result, self.ts * 2)

        # GH 10324
        a = Series([1, 2, 3, 4])
        b = Series(["even", "odd", "even", "odd"], dtype="category")
        c = Series(["even", "odd", "even", "odd"])

        exp = Series(["odd", "even", "odd", np.nan], dtype="category")
        self.assert_series_equal(a.map(b), exp)
        exp = Series(["odd", "even", "odd", np.nan])
        self.assert_series_equal(a.map(c), exp)

        a = Series(['a', 'b', 'c', 'd'])
        b = Series([1, 2, 3, 4], index=pd.CategoricalIndex(['b', 'c', 'd', 'e']))
        c = Series([1, 2, 3, 4], index=Index(['b', 'c', 'd', 'e']))

        exp = Series([np.nan, 1, 2, 3])
        self.assert_series_equal(a.map(b), exp)
        exp = Series([np.nan, 1, 2, 3])
        self.assert_series_equal(a.map(c), exp)

        a = Series(['a', 'b', 'c', 'd'])
        b = Series(['B', 'C', 'D', 'E'], dtype='category',
                   index=pd.CategoricalIndex(['b', 'c', 'd', 'e']))
        c = Series(['B', 'C', 'D', 'E'], index=Index(['b', 'c', 'd', 'e']))

        exp = Series([np.nan, 'B', 'C', 'D'], dtype='category')
        self.assert_series_equal(a.map(b), exp)
        exp = Series([np.nan, 'B', 'C', 'D'])
        self.assert_series_equal(a.map(c), exp)

    def test_map_compat(self):
        # related GH 8024
        s = Series([True,True,False],index=[1,2,3])
        result = s.map({ True : 'foo', False : 'bar' })
        expected = Series(['foo','foo','bar'],index=[1,2,3])
        assert_series_equal(result,expected)

    def test_map_int(self):
        left = Series({'a': 1., 'b': 2., 'c': 3., 'd': 4})
        right = Series({1: 11, 2: 22, 3: 33})

        self.assertEqual(left.dtype, np.float_)
        self.assertTrue(issubclass(right.dtype.type, np.integer))

        merged = left.map(right)
        self.assertEqual(merged.dtype, np.float_)
        self.assertTrue(isnull(merged['d']))
        self.assertTrue(not isnull(merged['c']))

    def test_map_type_inference(self):
        s = Series(lrange(3))
        s2 = s.map(lambda x: np.where(x == 0, 0, 1))
        self.assertTrue(issubclass(s2.dtype.type, np.integer))

    def test_divide_decimal(self):
        ''' resolves issue #9787 '''
        from decimal import Decimal

        expected = Series([Decimal(5)])

        s =  Series([Decimal(10)])
        s =  s/Decimal(2)

        tm.assert_series_equal(expected, s)

        s =  Series([Decimal(10)])
        s =  s//Decimal(2)

        tm.assert_series_equal(expected, s)

    def test_map_decimal(self):
        from decimal import Decimal

        result = self.series.map(lambda x: Decimal(str(x)))
        self.assertEqual(result.dtype, np.object_)
        tm.assertIsInstance(result[0], Decimal)

    def test_map_na_exclusion(self):
        s = Series([1.5, np.nan, 3, np.nan, 5])

        result = s.map(lambda x: x * 2, na_action='ignore')
        exp = s * 2
        assert_series_equal(result, exp)

    def test_map_dict_with_tuple_keys(self):
        '''
        Due to new MultiIndex-ing behaviour in v0.14.0,
        dicts with tuple keys passed to map were being
        converted to a multi-index, preventing tuple values
        from being mapped properly.
        '''
        df = pd.DataFrame({'a': [(1,), (2,), (3, 4), (5, 6)]})
        label_mappings = {
            (1,): 'A',
            (2,): 'B',
            (3, 4): 'A',
            (5, 6): 'B'
        }
        df['labels'] = df['a'].map(label_mappings)
        df['expected_labels'] = pd.Series(['A', 'B', 'A', 'B'], index=df.index)
        # All labels should be filled now
        tm.assert_series_equal(df['labels'], df['expected_labels'], check_names=False)

    def test_apply(self):
        assert_series_equal(self.ts.apply(np.sqrt), np.sqrt(self.ts))

        # elementwise-apply
        import math
        assert_series_equal(self.ts.apply(math.exp), np.exp(self.ts))

        # how to handle Series result, #2316
        result = self.ts.apply(lambda x: Series([x, x ** 2],
                                                index=['x', 'x^2']))
        expected = DataFrame({'x': self.ts, 'x^2': self.ts ** 2})
        tm.assert_frame_equal(result, expected)

        # empty series
        s = Series(dtype=object, name='foo', index=pd.Index([], name='bar'))
        rs = s.apply(lambda x: x)
        tm.assert_series_equal(s, rs)
        # check all metadata (GH 9322)
        self.assertIsNot(s, rs)
        self.assertIs(s.index, rs.index)
        self.assertEqual(s.dtype, rs.dtype)
        self.assertEqual(s.name, rs.name)

        # index but no data
        s = Series(index=[1, 2, 3])
        rs = s.apply(lambda x: x)
        tm.assert_series_equal(s, rs)

    def test_apply_same_length_inference_bug(self):
        s = Series([1, 2])
        f = lambda x: (x, x + 1)

        result = s.apply(f)
        expected = s.map(f)
        assert_series_equal(result, expected)

        s = Series([1, 2, 3])
        result = s.apply(f)
        expected = s.map(f)
        assert_series_equal(result, expected)

    def test_apply_dont_convert_dtype(self):
        s = Series(np.random.randn(10))

        f = lambda x: x if x > 0 else np.nan
        result = s.apply(f, convert_dtype=False)
        self.assertEqual(result.dtype, object)

    def test_convert_objects(self):

        s = Series([1., 2, 3], index=['a', 'b', 'c'])
        with tm.assert_produces_warning(FutureWarning):
            result = s.convert_objects(convert_dates=False, convert_numeric=True)
        assert_series_equal(result, s)

        # force numeric conversion
        r = s.copy().astype('O')
        r['a'] = '1'
        with tm.assert_produces_warning(FutureWarning):
            result = r.convert_objects(convert_dates=False, convert_numeric=True)
        assert_series_equal(result, s)

        r = s.copy().astype('O')
        r['a'] = '1.'
        with tm.assert_produces_warning(FutureWarning):
            result = r.convert_objects(convert_dates=False, convert_numeric=True)
        assert_series_equal(result, s)

        r = s.copy().astype('O')
        r['a'] = 'garbled'
        expected = s.copy()
        expected['a'] = np.nan
        with tm.assert_produces_warning(FutureWarning):
            result = r.convert_objects(convert_dates=False, convert_numeric=True)
        assert_series_equal(result, expected)

        # GH 4119, not converting a mixed type (e.g.floats and object)
        s = Series([1, 'na', 3, 4])
        with tm.assert_produces_warning(FutureWarning):
            result = s.convert_objects(convert_numeric=True)
        expected = Series([1, np.nan, 3, 4])
        assert_series_equal(result, expected)

        s = Series([1, '', 3, 4])
        with tm.assert_produces_warning(FutureWarning):
            result = s.convert_objects(convert_numeric=True)
        expected = Series([1, np.nan, 3, 4])
        assert_series_equal(result, expected)

        # dates
        s = Series(
            [datetime(2001, 1, 1, 0, 0), datetime(2001, 1, 2, 0, 0), datetime(2001, 1, 3, 0, 0)])
        s2 = Series([datetime(2001, 1, 1, 0, 0), datetime(2001, 1, 2, 0, 0), datetime(
            2001, 1, 3, 0, 0), 'foo', 1.0, 1, Timestamp('20010104'), '20010105'], dtype='O')
        with tm.assert_produces_warning(FutureWarning):
            result = s.convert_objects(convert_dates=True, convert_numeric=False)
        expected = Series(
            [Timestamp('20010101'), Timestamp('20010102'), Timestamp('20010103')], dtype='M8[ns]')
        assert_series_equal(result, expected)

        with tm.assert_produces_warning(FutureWarning):
            result = s.convert_objects(convert_dates='coerce',
                                       convert_numeric=False)
        with tm.assert_produces_warning(FutureWarning):
            result = s.convert_objects(convert_dates='coerce',
                                       convert_numeric=True)
        assert_series_equal(result, expected)

        expected = Series(
            [Timestamp(
                '20010101'), Timestamp('20010102'), Timestamp('20010103'),
                          lib.NaT, lib.NaT, lib.NaT, Timestamp('20010104'), Timestamp('20010105')], dtype='M8[ns]')
        with tm.assert_produces_warning(FutureWarning):
            result = s2.convert_objects(convert_dates='coerce',
                                        convert_numeric=False)
        assert_series_equal(result, expected)
        with tm.assert_produces_warning(FutureWarning):
            result = s2.convert_objects(convert_dates='coerce',
                                        convert_numeric=True)
        assert_series_equal(result, expected)

        # preserver all-nans (if convert_dates='coerce')
        s = Series(['foo', 'bar', 1, 1.0], dtype='O')
        with tm.assert_produces_warning(FutureWarning):
            result = s.convert_objects(convert_dates='coerce',
                                       convert_numeric=False)
        assert_series_equal(result, s)

        # preserver if non-object
        s = Series([1], dtype='float32')
        with tm.assert_produces_warning(FutureWarning):
            result = s.convert_objects(convert_dates='coerce',
                                       convert_numeric=False)
        assert_series_equal(result, s)

        #r = s.copy()
        #r[0] = np.nan
        #result = r.convert_objects(convert_dates=True,convert_numeric=False)
        #self.assertEqual(result.dtype, 'M8[ns]')

        # dateutil parses some single letters into today's value as a date
        for x in 'abcdefghijklmnopqrstuvwxyz':
            s = Series([x])
            with tm.assert_produces_warning(FutureWarning):
                result = s.convert_objects(convert_dates='coerce')
            assert_series_equal(result, s)
            s = Series([x.upper()])
            with tm.assert_produces_warning(FutureWarning):
                result = s.convert_objects(convert_dates='coerce')
            assert_series_equal(result, s)

    def test_convert_objects_preserve_bool(self):
        s = Series([1, True, 3, 5], dtype=object)
        with tm.assert_produces_warning(FutureWarning):
            r = s.convert_objects(convert_numeric=True)
        e = Series([1, 1, 3, 5], dtype='i8')
        tm.assert_series_equal(r, e)

    def test_convert_objects_preserve_all_bool(self):
        s = Series([False, True, False, False], dtype=object)
        with tm.assert_produces_warning(FutureWarning):
            r = s.convert_objects(convert_numeric=True)
        e = Series([False, True, False, False], dtype=bool)
        tm.assert_series_equal(r, e)

    # GH 10265
    def test_convert(self):
        # Tests: All to nans, coerce, true
        # Test coercion returns correct type
        s = Series(['a', 'b', 'c'])
        results = s._convert(datetime=True, coerce=True)
        expected = Series([lib.NaT] * 3)
        assert_series_equal(results, expected)

        results = s._convert(numeric=True, coerce=True)
        expected = Series([np.nan] * 3)
        assert_series_equal(results, expected)

        expected = Series([lib.NaT] * 3, dtype=np.dtype('m8[ns]'))
        results = s._convert(timedelta=True, coerce=True)
        assert_series_equal(results, expected)

        dt = datetime(2001, 1, 1, 0, 0)
        td = dt - datetime(2000, 1, 1, 0, 0)

        # Test coercion with mixed types
        s = Series(['a', '3.1415', dt, td])
        results = s._convert(datetime=True, coerce=True)
        expected = Series([lib.NaT, lib.NaT, dt, lib.NaT])
        assert_series_equal(results, expected)

        results = s._convert(numeric=True, coerce=True)
        expected = Series([nan, 3.1415, nan, nan])
        assert_series_equal(results, expected)

        results = s._convert(timedelta=True, coerce=True)
        expected = Series([lib.NaT, lib.NaT, lib.NaT, td],
                          dtype=np.dtype('m8[ns]'))
        assert_series_equal(results, expected)

        # Test standard conversion returns original
        results = s._convert(datetime=True)
        assert_series_equal(results, s)
        results = s._convert(numeric=True)
        expected = Series([nan, 3.1415, nan, nan])
        assert_series_equal(results, expected)
        results = s._convert(timedelta=True)
        assert_series_equal(results, s)

        # test pass-through and non-conversion when other types selected
        s = Series(['1.0','2.0','3.0'])
        results = s._convert(datetime=True, numeric=True, timedelta=True)
        expected = Series([1.0,2.0,3.0])
        assert_series_equal(results, expected)
        results = s._convert(True,False,True)
        assert_series_equal(results, s)

        s = Series([datetime(2001, 1, 1, 0, 0),datetime(2001, 1, 1, 0, 0)],
                   dtype='O')
        results = s._convert(datetime=True, numeric=True, timedelta=True)
        expected = Series([datetime(2001, 1, 1, 0, 0),datetime(2001, 1, 1, 0, 0)])
        assert_series_equal(results, expected)
        results = s._convert(datetime=False,numeric=True,timedelta=True)
        assert_series_equal(results, s)

        td = datetime(2001, 1, 1, 0, 0) - datetime(2000, 1, 1, 0, 0)
        s = Series([td, td], dtype='O')
        results = s._convert(datetime=True, numeric=True, timedelta=True)
        expected = Series([td, td])
        assert_series_equal(results, expected)
        results = s._convert(True,True,False)
        assert_series_equal(results, s)


        s = Series([1., 2, 3], index=['a', 'b', 'c'])
        result = s._convert(numeric=True)
        assert_series_equal(result, s)

        # force numeric conversion
        r = s.copy().astype('O')
        r['a'] = '1'
        result = r._convert(numeric=True)
        assert_series_equal(result, s)

        r = s.copy().astype('O')
        r['a'] = '1.'
        result = r._convert(numeric=True)
        assert_series_equal(result, s)

        r = s.copy().astype('O')
        r['a'] = 'garbled'
        result = r._convert(numeric=True)
        expected = s.copy()
        expected['a'] = nan
        assert_series_equal(result, expected)

        # GH 4119, not converting a mixed type (e.g.floats and object)
        s = Series([1, 'na', 3, 4])
        result = s._convert(datetime=True, numeric=True)
        expected = Series([1, nan, 3, 4])
        assert_series_equal(result, expected)

        s = Series([1, '', 3, 4])
        result = s._convert(datetime=True, numeric=True)
        assert_series_equal(result, expected)

        # dates
        s = Series(
            [datetime(2001, 1, 1, 0, 0), datetime(2001, 1, 2, 0, 0), datetime(2001, 1, 3, 0, 0)])
        s2 = Series([datetime(2001, 1, 1, 0, 0), datetime(2001, 1, 2, 0, 0), datetime(
            2001, 1, 3, 0, 0), 'foo', 1.0, 1, Timestamp('20010104'), '20010105'], dtype='O')

        result = s._convert(datetime=True)
        expected = Series(
            [Timestamp('20010101'), Timestamp('20010102'), Timestamp('20010103')], dtype='M8[ns]')
        assert_series_equal(result, expected)

        result = s._convert(datetime=True, coerce=True)
        assert_series_equal(result, expected)

        expected = Series(
            [Timestamp(
                '20010101'), Timestamp('20010102'), Timestamp('20010103'),
                          lib.NaT, lib.NaT, lib.NaT, Timestamp('20010104'), Timestamp('20010105')], dtype='M8[ns]')
        result = s2._convert(datetime=True,
                                    numeric=False,
                                    timedelta=False,
                                    coerce=True)
        assert_series_equal(result, expected)
        result = s2._convert(datetime=True, coerce=True)
        assert_series_equal(result, expected)

        s = Series(['foo', 'bar', 1, 1.0], dtype='O')
        result = s._convert(datetime=True, coerce=True)
        expected = Series([lib.NaT]*4)
        assert_series_equal(result, expected)

        # preserver if non-object
        s = Series([1], dtype='float32')
        result = s._convert(datetime=True, coerce=True)
        assert_series_equal(result, s)

        #r = s.copy()
        #r[0] = np.nan
        #result = r._convert(convert_dates=True,convert_numeric=False)
        #self.assertEqual(result.dtype, 'M8[ns]')

        # dateutil parses some single letters into today's value as a date
        expected = Series([lib.NaT])
        for x in 'abcdefghijklmnopqrstuvwxyz':
            s = Series([x])
            result = s._convert(datetime=True, coerce=True)
            assert_series_equal(result, expected)
            s = Series([x.upper()])
            result = s._convert(datetime=True, coerce=True)
            assert_series_equal(result, expected)

    def test_convert_no_arg_error(self):
        s = Series(['1.0','2'])
        self.assertRaises(ValueError, s._convert)

    def test_convert_preserve_bool(self):
        s = Series([1, True, 3, 5], dtype=object)
        r = s._convert(datetime=True, numeric=True)
        e = Series([1, 1, 3, 5], dtype='i8')
        tm.assert_series_equal(r, e)

    def test_convert_preserve_all_bool(self):
        s = Series([False, True, False, False], dtype=object)
        r = s._convert(datetime=True, numeric=True)
        e = Series([False, True, False, False], dtype=bool)
        tm.assert_series_equal(r, e)

    def test_apply_args(self):
        s = Series(['foo,bar'])

        result = s.apply(str.split, args=(',',))
        self.assertEqual(result[0], ['foo', 'bar'])
        tm.assertIsInstance(result[0], list)

    def test_align(self):
        def _check_align(a, b, how='left', fill=None):
            aa, ab = a.align(b, join=how, fill_value=fill)

            join_index = a.index.join(b.index, how=how)
            if fill is not None:
                diff_a = aa.index.difference(join_index)
                diff_b = ab.index.difference(join_index)
                if len(diff_a) > 0:
                    self.assertTrue((aa.reindex(diff_a) == fill).all())
                if len(diff_b) > 0:
                    self.assertTrue((ab.reindex(diff_b) == fill).all())

            ea = a.reindex(join_index)
            eb = b.reindex(join_index)

            if fill is not None:
                ea = ea.fillna(fill)
                eb = eb.fillna(fill)

            assert_series_equal(aa, ea)
            assert_series_equal(ab, eb)
            self.assertEqual(aa.name, 'ts')
            self.assertEqual(ea.name, 'ts')
            self.assertEqual(ab.name, 'ts')
            self.assertEqual(eb.name, 'ts')

        for kind in JOIN_TYPES:
            _check_align(self.ts[2:], self.ts[:-5], how=kind)
            _check_align(self.ts[2:], self.ts[:-5], how=kind, fill=-1)

            # empty left
            _check_align(self.ts[:0], self.ts[:-5], how=kind)
            _check_align(self.ts[:0], self.ts[:-5], how=kind, fill=-1)

            # empty right
            _check_align(self.ts[:-5], self.ts[:0], how=kind)
            _check_align(self.ts[:-5], self.ts[:0], how=kind, fill=-1)

            # both empty
            _check_align(self.ts[:0], self.ts[:0], how=kind)
            _check_align(self.ts[:0], self.ts[:0], how=kind, fill=-1)

    def test_align_fill_method(self):
        def _check_align(a, b, how='left', method='pad', limit=None):
            aa, ab = a.align(b, join=how, method=method, limit=limit)

            join_index = a.index.join(b.index, how=how)
            ea = a.reindex(join_index)
            eb = b.reindex(join_index)

            ea = ea.fillna(method=method, limit=limit)
            eb = eb.fillna(method=method, limit=limit)

            assert_series_equal(aa, ea)
            assert_series_equal(ab, eb)

        for kind in JOIN_TYPES:
            for meth in ['pad', 'bfill']:
                _check_align(self.ts[2:], self.ts[:-5], how=kind, method=meth)
                _check_align(self.ts[2:], self.ts[:-5], how=kind,
                             method=meth, limit=1)

                # empty left
                _check_align(self.ts[:0], self.ts[:-5], how=kind, method=meth)
                _check_align(self.ts[:0], self.ts[:-5], how=kind, method=meth,
                             limit=1)

                # empty right
                _check_align(self.ts[:-5], self.ts[:0], how=kind, method=meth)
                _check_align(self.ts[:-5], self.ts[:0], how=kind, method=meth,
                             limit=1)

                # both empty
                _check_align(self.ts[:0], self.ts[:0], how=kind, method=meth)
                _check_align(self.ts[:0], self.ts[:0], how=kind, method=meth,
                             limit=1)

    def test_align_nocopy(self):
        b = self.ts[:5].copy()

        # do copy
        a = self.ts.copy()
        ra, _ = a.align(b, join='left')
        ra[:5] = 5
        self.assertFalse((a[:5] == 5).any())

        # do not copy
        a = self.ts.copy()
        ra, _ = a.align(b, join='left', copy=False)
        ra[:5] = 5
        self.assertTrue((a[:5] == 5).all())

        # do copy
        a = self.ts.copy()
        b = self.ts[:5].copy()
        _, rb = a.align(b, join='right')
        rb[:3] = 5
        self.assertFalse((b[:3] == 5).any())

        # do not copy
        a = self.ts.copy()
        b = self.ts[:5].copy()
        _, rb = a.align(b, join='right', copy=False)
        rb[:2] = 5
        self.assertTrue((b[:2] == 5).all())

    def test_align_sameindex(self):
        a, b = self.ts.align(self.ts, copy=False)
        self.assertIs(a.index, self.ts.index)
        self.assertIs(b.index, self.ts.index)

        # a, b = self.ts.align(self.ts, copy=True)
        # self.assertIsNot(a.index, self.ts.index)
        # self.assertIsNot(b.index, self.ts.index)

    def test_align_multiindex(self):
        # GH 10665

        midx = pd.MultiIndex.from_product([range(2), range(3), range(2)],
                                 names=('a', 'b', 'c'))
        idx = pd.Index(range(2), name='b')
        s1 = pd.Series(np.arange(12,dtype='int64'), index=midx)
        s2 = pd.Series(np.arange(2,dtype='int64'), index=idx)

        # these must be the same results (but flipped)
        res1l, res1r = s1.align(s2, join='left')
        res2l, res2r = s2.align(s1, join='right')

        expl = s1
        tm.assert_series_equal(expl, res1l)
        tm.assert_series_equal(expl, res2r)
        expr = pd.Series([0, 0, 1, 1, np.nan, np.nan] * 2, index=midx)
        tm.assert_series_equal(expr, res1r)
        tm.assert_series_equal(expr, res2l)

        res1l, res1r = s1.align(s2, join='right')
        res2l, res2r = s2.align(s1, join='left')

        exp_idx = pd.MultiIndex.from_product([range(2), range(2), range(2)],
                                             names=('a', 'b', 'c'))
        expl = pd.Series([0, 1, 2, 3, 6, 7, 8, 9], index=exp_idx)
        tm.assert_series_equal(expl, res1l)
        tm.assert_series_equal(expl, res2r)
        expr = pd.Series([0, 0, 1, 1] * 2, index=exp_idx)
        tm.assert_series_equal(expr, res1r)
        tm.assert_series_equal(expr, res2l)

    def test_reindex(self):

        identity = self.series.reindex(self.series.index)

        # __array_interface__ is not defined for older numpies
        # and on some pythons
        try:
            self.assertTrue(np.may_share_memory(self.series.index, identity.index))
        except (AttributeError):
            pass

        self.assertTrue(identity.index.is_(self.series.index))
        self.assertTrue(identity.index.identical(self.series.index))

        subIndex = self.series.index[10:20]
        subSeries = self.series.reindex(subIndex)

        for idx, val in compat.iteritems(subSeries):
            self.assertEqual(val, self.series[idx])

        subIndex2 = self.ts.index[10:20]
        subTS = self.ts.reindex(subIndex2)

        for idx, val in compat.iteritems(subTS):
            self.assertEqual(val, self.ts[idx])
        stuffSeries = self.ts.reindex(subIndex)

        self.assertTrue(np.isnan(stuffSeries).all())

        # This is extremely important for the Cython code to not screw up
        nonContigIndex = self.ts.index[::2]
        subNonContig = self.ts.reindex(nonContigIndex)
        for idx, val in compat.iteritems(subNonContig):
            self.assertEqual(val, self.ts[idx])

        # return a copy the same index here
        result = self.ts.reindex()
        self.assertFalse((result is self.ts))

    def test_reindex_nan(self):
        ts = Series([2, 3, 5, 7], index=[1, 4, nan, 8])

        i, j = [nan, 1, nan, 8, 4, nan], [2, 0, 2, 3, 1, 2]
        assert_series_equal(ts.reindex(i), ts.iloc[j])

        ts.index = ts.index.astype('object')
        # reindex coerces index.dtype to float, loc/iloc doesn't
        assert_series_equal(ts.reindex(i), ts.iloc[j], check_index_type=False)

    def test_reindex_corner(self):
        # (don't forget to fix this) I think it's fixed
        reindexed_dep = self.empty.reindex(self.ts.index, method='pad')

        # corner case: pad empty series
        reindexed = self.empty.reindex(self.ts.index, method='pad')

        # pass non-Index
        reindexed = self.ts.reindex(list(self.ts.index))
        assert_series_equal(self.ts, reindexed)

        # bad fill method
        ts = self.ts[::2]
        self.assertRaises(Exception, ts.reindex, self.ts.index, method='foo')

    def test_reindex_pad(self):

        s = Series(np.arange(10),dtype='int64')
        s2 = s[::2]

        reindexed = s2.reindex(s.index, method='pad')
        reindexed2 = s2.reindex(s.index, method='ffill')
        assert_series_equal(reindexed, reindexed2)

        expected = Series([0, 0, 2, 2, 4, 4, 6, 6, 8, 8], index=np.arange(10))
        assert_series_equal(reindexed, expected)

        # GH4604
        s = Series([1,2,3,4,5], index=['a', 'b', 'c', 'd', 'e'])
        new_index = ['a','g','c','f']
        expected = Series([1,1,3,3],index=new_index)

        # this changes dtype because the ffill happens after
        result = s.reindex(new_index).ffill()
        assert_series_equal(result, expected.astype('float64'))

        result = s.reindex(new_index).ffill(downcast='infer')
        assert_series_equal(result, expected)

        expected = Series([1, 5, 3, 5], index=new_index)
        result = s.reindex(new_index, method='ffill')
        assert_series_equal(result, expected)

        # inferrence of new dtype
        s = Series([True,False,False,True],index=list('abcd'))
        new_index='agc'
        result = s.reindex(list(new_index)).ffill()
        expected = Series([True,True,False],index=list(new_index))
        assert_series_equal(result, expected)

        # GH4618 shifted series downcasting
        s = Series(False,index=lrange(0,5))
        result = s.shift(1).fillna(method='bfill')
        expected = Series(False,index=lrange(0,5))
        assert_series_equal(result, expected)

    def test_reindex_nearest(self):
        s = Series(np.arange(10, dtype='int64'))
        target = [0.1, 0.9, 1.5, 2.0]
        actual = s.reindex(target, method='nearest')
        expected = Series(np.around(target).astype('int64'), target)
        assert_series_equal(expected, actual)

        actual = s.reindex_like(actual, method='nearest')
        assert_series_equal(expected, actual)

        actual = s.reindex_like(actual, method='nearest', tolerance=1)
        assert_series_equal(expected, actual)

        actual = s.reindex(target, method='nearest', tolerance=0.2)
        expected = Series([0, 1, np.nan, 2], target)
        assert_series_equal(expected, actual)

    def test_reindex_backfill(self):
        pass

    def test_reindex_int(self):
        ts = self.ts[::2]
        int_ts = Series(np.zeros(len(ts), dtype=int), index=ts.index)

        # this should work fine
        reindexed_int = int_ts.reindex(self.ts.index)

        # if NaNs introduced
        self.assertEqual(reindexed_int.dtype, np.float_)

        # NO NaNs introduced
        reindexed_int = int_ts.reindex(int_ts.index[::2])
        self.assertEqual(reindexed_int.dtype, np.int_)

    def test_reindex_bool(self):

        # A series other than float, int, string, or object
        ts = self.ts[::2]
        bool_ts = Series(np.zeros(len(ts), dtype=bool), index=ts.index)

        # this should work fine
        reindexed_bool = bool_ts.reindex(self.ts.index)

        # if NaNs introduced
        self.assertEqual(reindexed_bool.dtype, np.object_)

        # NO NaNs introduced
        reindexed_bool = bool_ts.reindex(bool_ts.index[::2])
        self.assertEqual(reindexed_bool.dtype, np.bool_)

    def test_reindex_bool_pad(self):
        # fail
        ts = self.ts[5:]
        bool_ts = Series(np.zeros(len(ts), dtype=bool), index=ts.index)
        filled_bool = bool_ts.reindex(self.ts.index, method='pad')
        self.assertTrue(isnull(filled_bool[:5]).all())

    def test_reindex_like(self):
        other = self.ts[::2]
        assert_series_equal(self.ts.reindex(other.index),
                            self.ts.reindex_like(other))

        # GH 7179
        day1 = datetime(2013,3,5)
        day2 = datetime(2013,5,5)
        day3 = datetime(2014,3,5)

        series1 = Series([5, None, None],[day1, day2, day3])
        series2 = Series([None, None], [day1, day3])

        result = series1.reindex_like(series2, method='pad')
        expected = Series([5, np.nan], index=[day1, day3])
        assert_series_equal(result, expected)

    def test_reindex_fill_value(self):
        #------------------------------------------------------------
        # floats
        floats = Series([1., 2., 3.])
        result = floats.reindex([1, 2, 3])
        expected = Series([2., 3., np.nan], index=[1, 2, 3])
        assert_series_equal(result, expected)

        result = floats.reindex([1, 2, 3], fill_value=0)
        expected = Series([2., 3., 0], index=[1, 2, 3])
        assert_series_equal(result, expected)

        #------------------------------------------------------------
        # ints
        ints = Series([1, 2, 3])

        result = ints.reindex([1, 2, 3])
        expected = Series([2., 3., np.nan], index=[1, 2, 3])
        assert_series_equal(result, expected)

        # don't upcast
        result = ints.reindex([1, 2, 3], fill_value=0)
        expected = Series([2, 3, 0], index=[1, 2, 3])
        self.assertTrue(issubclass(result.dtype.type, np.integer))
        assert_series_equal(result, expected)

        #------------------------------------------------------------
        # objects
        objects = Series([1, 2, 3], dtype=object)

        result = objects.reindex([1, 2, 3])
        expected = Series([2, 3, np.nan], index=[1, 2, 3], dtype=object)
        assert_series_equal(result, expected)

        result = objects.reindex([1, 2, 3], fill_value='foo')
        expected = Series([2, 3, 'foo'], index=[1, 2, 3], dtype=object)
        assert_series_equal(result, expected)

        #------------------------------------------------------------
        # bools
        bools = Series([True, False, True])

        result = bools.reindex([1, 2, 3])
        expected = Series([False, True, np.nan], index=[1, 2, 3], dtype=object)
        assert_series_equal(result, expected)

        result = bools.reindex([1, 2, 3], fill_value=False)
        expected = Series([False, True, False], index=[1, 2, 3])
        assert_series_equal(result, expected)

    def test_rename(self):
        renamer = lambda x: x.strftime('%Y%m%d')
        renamed = self.ts.rename(renamer)
        self.assertEqual(renamed.index[0], renamer(self.ts.index[0]))

        # dict
        rename_dict = dict(zip(self.ts.index, renamed.index))
        renamed2 = self.ts.rename(rename_dict)
        assert_series_equal(renamed, renamed2)

        # partial dict
        s = Series(np.arange(4), index=['a', 'b', 'c', 'd'], dtype='int64')
        renamed = s.rename({'b': 'foo', 'd': 'bar'})
        self.assert_numpy_array_equal(renamed.index, ['a', 'foo', 'c', 'bar'])

        # index with name
        renamer = Series(
            np.arange(4), index=Index(['a', 'b', 'c', 'd'], name='name'), dtype='int64')
        renamed = renamer.rename({})
        self.assertEqual(renamed.index.name, renamer.index.name)

    def test_rename_inplace(self):
        renamer = lambda x: x.strftime('%Y%m%d')
        expected = renamer(self.ts.index[0])

        self.ts.rename(renamer, inplace=True)
        self.assertEqual(self.ts.index[0], expected)

    def test_preserveRefs(self):
        seq = self.ts[[5, 10, 15]]
        seq[1] = np.NaN
        self.assertFalse(np.isnan(self.ts[10]))

    def test_ne(self):
        ts = Series([3, 4, 5, 6, 7], [3, 4, 5, 6, 7], dtype=float)
        expected = [True, True, False, True, True]
        self.assertTrue(tm.equalContents(ts.index != 5, expected))
        self.assertTrue(tm.equalContents(~(ts.index == 5), expected))

    def test_pad_nan(self):
        x = Series([np.nan, 1., np.nan, 3., np.nan],
                   ['z', 'a', 'b', 'c', 'd'], dtype=float)

        x.fillna(method='pad', inplace=True)

        expected = Series([np.nan, 1.0, 1.0, 3.0, 3.0],
                          ['z', 'a', 'b', 'c', 'd'], dtype=float)
        assert_series_equal(x[1:], expected[1:])
        self.assertTrue(np.isnan(x[0]), np.isnan(expected[0]))

    def test_unstack(self):
        from numpy import nan
        from pandas.util.testing import assert_frame_equal

        index = MultiIndex(levels=[['bar', 'foo'], ['one', 'three', 'two']],
                           labels=[[1, 1, 0, 0], [0, 1, 0, 2]])

        s = Series(np.arange(4.), index=index)
        unstacked = s.unstack()

        expected = DataFrame([[2., nan, 3.], [0., 1., nan]],
                             index=['bar', 'foo'],
                             columns=['one', 'three', 'two'])

        assert_frame_equal(unstacked, expected)

        unstacked = s.unstack(level=0)
        assert_frame_equal(unstacked, expected.T)

        index = MultiIndex(levels=[['bar'], ['one', 'two', 'three'], [0, 1]],
                           labels=[[0, 0, 0, 0, 0, 0],
                                   [0, 1, 2, 0, 1, 2],
                                   [0, 1, 0, 1, 0, 1]])
        s = Series(np.random.randn(6), index=index)
        exp_index = MultiIndex(levels=[['one', 'two', 'three'], [0, 1]],
                               labels=[[0, 1, 2, 0, 1, 2],
                                       [0, 1, 0, 1, 0, 1]])
        expected = DataFrame({'bar': s.values}, index=exp_index).sortlevel(0)
        unstacked = s.unstack(0)
        assert_frame_equal(unstacked, expected)

        # GH5873
        idx = pd.MultiIndex.from_arrays([[101, 102], [3.5, np.nan]])
        ts = pd.Series([1,2], index=idx)
        left = ts.unstack()
        right = DataFrame([[nan, 1], [2, nan]], index=[101, 102],
                          columns=[nan, 3.5])
        print(left)
        print(right)
        assert_frame_equal(left, right)

        idx = pd.MultiIndex.from_arrays([['cat', 'cat', 'cat', 'dog', 'dog'],
                ['a', 'a', 'b', 'a', 'b'], [1, 2, 1, 1, np.nan]])
        ts = pd.Series([1.0, 1.1, 1.2, 1.3, 1.4], index=idx)
        right = DataFrame([[1.0, 1.3], [1.1, nan], [nan, 1.4], [1.2, nan]],
                          columns=['cat', 'dog'])
        tpls = [('a', 1), ('a', 2), ('b', nan), ('b', 1)]
        right.index = pd.MultiIndex.from_tuples(tpls)
        assert_frame_equal(ts.unstack(level=0), right)

    def test_sortlevel(self):
        mi = MultiIndex.from_tuples([[1, 1, 3], [1, 1, 1]], names=list('ABC'))
        s = Series([1, 2], mi)
        backwards = s.iloc[[1, 0]]

        res = s.sortlevel('A')
        assert_series_equal(backwards, res)

        res = s.sortlevel(['A', 'B'])
        assert_series_equal(backwards, res)

        res = s.sortlevel('A', sort_remaining=False)
        assert_series_equal(s, res)

        res = s.sortlevel(['A', 'B'], sort_remaining=False)
        assert_series_equal(s, res)

    def test_head_tail(self):
        assert_series_equal(self.series.head(), self.series[:5])
        assert_series_equal(self.series.tail(), self.series[-5:])

    def test_isin(self):
        s = Series(['A', 'B', 'C', 'a', 'B', 'B', 'A', 'C'])

        result = s.isin(['A', 'C'])
        expected = Series([True, False, True, False, False, False, True, True])
        assert_series_equal(result, expected)

    def test_isin_with_string_scalar(self):
        # GH4763
        s = Series(['A', 'B', 'C', 'a', 'B', 'B', 'A', 'C'])
        with tm.assertRaises(TypeError):
            s.isin('a')

        with tm.assertRaises(TypeError):
            s = Series(['aaa', 'b', 'c'])
            s.isin('aaa')

    def test_isin_with_i8(self):
        # GH 5021

        expected = Series([True,True,False,False,False])
        expected2 = Series([False,True,False,False,False])

        # datetime64[ns]
        s = Series(date_range('jan-01-2013','jan-05-2013'))

        result = s.isin(s[0:2])
        assert_series_equal(result, expected)

        result = s.isin(s[0:2].values)
        assert_series_equal(result, expected)

        # fails on dtype conversion in the first place
        result = s.isin(s[0:2].values.astype('datetime64[D]'))
        assert_series_equal(result, expected)

        result = s.isin([s[1]])
        assert_series_equal(result, expected2)

        result = s.isin([np.datetime64(s[1])])
        assert_series_equal(result, expected2)

        # timedelta64[ns]
        s = Series(pd.to_timedelta(lrange(5),unit='d'))
        result = s.isin(s[0:2])
        assert_series_equal(result, expected)

#------------------------------------------------------------------------------
# TimeSeries-specific
    def test_cummethods_bool(self):
        # GH 6270
        # looks like a buggy np.maximum.accumulate for numpy 1.6.1, py 3.2
        def cummin(x):
            return np.minimum.accumulate(x)

        def cummax(x):
            return np.maximum.accumulate(x)

        a = pd.Series([False, False, False, True, True, False, False])
        b = ~a
        c = pd.Series([False] * len(b))
        d = ~c
        methods = {'cumsum': np.cumsum, 'cumprod': np.cumprod,
                   'cummin': cummin, 'cummax': cummax}
        args = product((a, b, c, d), methods)
        for s, method in args:
            expected = Series(methods[method](s.values))
            result = getattr(s, method)()
            assert_series_equal(result, expected)

        e = pd.Series([False, True, nan, False])
        cse = pd.Series([0, 1, nan, 1], dtype=object)
        cpe = pd.Series([False, 0, nan, 0])
        cmin = pd.Series([False, False, nan, False])
        cmax = pd.Series([False, True, nan, True])
        expecteds = {'cumsum': cse, 'cumprod': cpe, 'cummin': cmin,
                     'cummax': cmax}

        for method in methods:
            res = getattr(e, method)()
            assert_series_equal(res, expecteds[method])

    def test_replace(self):
        N = 100
        ser = Series(np.random.randn(N))
        ser[0:4] = np.nan
        ser[6:10] = 0

        # replace list with a single value
        ser.replace([np.nan], -1, inplace=True)

        exp = ser.fillna(-1)
        assert_series_equal(ser, exp)

        rs = ser.replace(0., np.nan)
        ser[ser == 0.] = np.nan
        assert_series_equal(rs, ser)

        ser = Series(np.fabs(np.random.randn(N)), tm.makeDateIndex(N),
                     dtype=object)
        ser[:5] = np.nan
        ser[6:10] = 'foo'
        ser[20:30] = 'bar'

        # replace list with a single value
        rs = ser.replace([np.nan, 'foo', 'bar'], -1)

        self.assertTrue((rs[:5] == -1).all())
        self.assertTrue((rs[6:10] == -1).all())
        self.assertTrue((rs[20:30] == -1).all())
        self.assertTrue((isnull(ser[:5])).all())

        # replace with different values
        rs = ser.replace({np.nan: -1, 'foo': -2, 'bar': -3})

        self.assertTrue((rs[:5] == -1).all())
        self.assertTrue((rs[6:10] == -2).all())
        self.assertTrue((rs[20:30] == -3).all())
        self.assertTrue((isnull(ser[:5])).all())

        # replace with different values with 2 lists
        rs2 = ser.replace([np.nan, 'foo', 'bar'], [-1, -2, -3])
        assert_series_equal(rs, rs2)

        # replace inplace
        ser.replace([np.nan, 'foo', 'bar'], -1, inplace=True)

        self.assertTrue((ser[:5] == -1).all())
        self.assertTrue((ser[6:10] == -1).all())
        self.assertTrue((ser[20:30] == -1).all())

        ser = Series([np.nan, 0, np.inf])
        assert_series_equal(ser.replace(np.nan, 0), ser.fillna(0))

        ser = Series([np.nan, 0, 'foo', 'bar', np.inf, None, lib.NaT])
        assert_series_equal(ser.replace(np.nan, 0), ser.fillna(0))
        filled = ser.copy()
        filled[4] = 0
        assert_series_equal(ser.replace(np.inf, 0), filled)

        ser = Series(self.ts.index)
        assert_series_equal(ser.replace(np.nan, 0), ser.fillna(0))

        # malformed
        self.assertRaises(ValueError, ser.replace, [1, 2, 3], [np.nan, 0])

        # make sure that we aren't just masking a TypeError because bools don't
        # implement indexing
        with tm.assertRaisesRegexp(TypeError, 'Cannot compare types .+'):
            ser.replace([1, 2], [np.nan, 0])

        ser = Series([0, 1, 2, 3, 4])
        result = ser.replace([0, 1, 2, 3, 4], [4, 3, 2, 1, 0])
        assert_series_equal(result, Series([4, 3, 2, 1, 0]))

        # API change from 0.12?
        # GH 5319
        ser = Series([0, np.nan, 2, 3, 4])
        expected = ser.ffill()
        result = ser.replace([np.nan])
        assert_series_equal(result, expected)

        ser = Series([0, np.nan, 2, 3, 4])
        expected = ser.ffill()
        result = ser.replace(np.nan)
        assert_series_equal(result, expected)
        #GH 5797
        ser = Series(date_range('20130101', periods=5))
        expected = ser.copy()
        expected.loc[2] = Timestamp('20120101')
        result = ser.replace({Timestamp('20130103'):
                              Timestamp('20120101')})
        assert_series_equal(result, expected)
        result = ser.replace(Timestamp('20130103'), Timestamp('20120101'))
        assert_series_equal(result, expected)

    def test_replace_with_single_list(self):
        ser = Series([0, 1, 2, 3, 4])
        result = ser.replace([1,2,3])
        assert_series_equal(result, Series([0,0,0,0,4]))

        s = ser.copy()
        s.replace([1,2,3],inplace=True)
        assert_series_equal(s, Series([0,0,0,0,4]))

        # make sure things don't get corrupted when fillna call fails
        s = ser.copy()
        with tm.assertRaises(ValueError):
            s.replace([1,2,3],inplace=True,method='crash_cymbal')
        assert_series_equal(s, ser)

    def test_replace_mixed_types(self):
        s = Series(np.arange(5),dtype='int64')

        def check_replace(to_rep, val, expected):
            sc = s.copy()
            r = s.replace(to_rep, val)
            sc.replace(to_rep, val, inplace=True)
            assert_series_equal(expected, r)
            assert_series_equal(expected, sc)

        # should NOT upcast to float
        e = Series([0,1,2,3,4])
        tr, v = [3], [3.0]
        check_replace(tr, v, e)

        # MUST upcast to float
        e = Series([0,1,2,3.5,4])
        tr, v = [3], [3.5]
        check_replace(tr, v, e)

        # casts to object
        e = Series([0,1,2,3.5,'a'])
        tr, v = [3,4], [3.5,'a']
        check_replace(tr, v, e)

        # again casts to object
        e = Series([0,1,2,3.5,Timestamp('20130101')])
        tr, v = [3,4],[3.5,Timestamp('20130101')]
        check_replace(tr, v, e)

        # casts to float
        e = Series([0,1,2,3.5,1])
        tr, v = [3,4],[3.5,True]
        check_replace(tr, v, e)

        # test an object with dates + floats + integers + strings
        dr = date_range('1/1/2001', '1/10/2001',
                        freq='D').to_series().reset_index(drop=True)
        result = dr.astype(object).replace([dr[0],dr[1],dr[2]], [1.0,2,'a'])
        expected = Series([1.0,2,'a'] + dr[3:].tolist(),dtype=object)
        assert_series_equal(result, expected)

    def test_replace_bool_with_string_no_op(self):
        s = Series([True, False, True])
        result = s.replace('fun', 'in-the-sun')
        tm.assert_series_equal(s, result)

    def test_replace_bool_with_string(self):
        # nonexistent elements
        s = Series([True, False, True])
        result = s.replace(True, '2u')
        expected = Series(['2u', False, '2u'])
        tm.assert_series_equal(expected, result)

    def test_replace_bool_with_bool(self):
        s = Series([True, False, True])
        result = s.replace(True, False)
        expected = Series([False] * len(s))
        tm.assert_series_equal(expected, result)

    def test_replace_with_dict_with_bool_keys(self):
        s = Series([True, False, True])
        with tm.assertRaisesRegexp(TypeError, 'Cannot compare types .+'):
            s.replace({'asdf': 'asdb', True: 'yes'})

    def test_asfreq(self):
        ts = Series([0., 1., 2.], index=[datetime(2009, 10, 30),
                                         datetime(2009, 11, 30),
                                         datetime(2009, 12, 31)])

        daily_ts = ts.asfreq('B')
        monthly_ts = daily_ts.asfreq('BM')
        self.assert_numpy_array_equal(monthly_ts, ts)

        daily_ts = ts.asfreq('B', method='pad')
        monthly_ts = daily_ts.asfreq('BM')
        self.assert_numpy_array_equal(monthly_ts, ts)

        daily_ts = ts.asfreq(datetools.bday)
        monthly_ts = daily_ts.asfreq(datetools.bmonthEnd)
        self.assert_numpy_array_equal(monthly_ts, ts)

        result = ts[:0].asfreq('M')
        self.assertEqual(len(result), 0)
        self.assertIsNot(result, ts)

    def test_diff(self):
        # Just run the function
        self.ts.diff()

        # int dtype
        a = 10000000000000000
        b = a + 1
        s = Series([a, b])

        rs = s.diff()
        self.assertEqual(rs[1], 1)

        # neg n
        rs = self.ts.diff(-1)
        xp = self.ts - self.ts.shift(-1)
        assert_series_equal(rs, xp)

        # 0
        rs = self.ts.diff(0)
        xp = self.ts - self.ts
        assert_series_equal(rs, xp)

        # datetime diff (GH3100)
        s = Series(date_range('20130102', periods=5))
        rs = s - s.shift(1)
        xp = s.diff()
        assert_series_equal(rs, xp)

        # timedelta diff
        nrs = rs - rs.shift(1)
        nxp = xp.diff()
        assert_series_equal(nrs, nxp)

        # with tz
        s = Series(date_range('2000-01-01 09:00:00',periods=5,tz='US/Eastern'), name='foo')
        result = s.diff()
        assert_series_equal(result,Series(TimedeltaIndex(['NaT'] + ['1 days']*4),name='foo'))

    def test_pct_change(self):
        rs = self.ts.pct_change(fill_method=None)
        assert_series_equal(rs, self.ts / self.ts.shift(1) - 1)

        rs = self.ts.pct_change(2)
        filled = self.ts.fillna(method='pad')
        assert_series_equal(rs, filled / filled.shift(2) - 1)

        rs = self.ts.pct_change(fill_method='bfill', limit=1)
        filled = self.ts.fillna(method='bfill', limit=1)
        assert_series_equal(rs, filled / filled.shift(1) - 1)

        rs = self.ts.pct_change(freq='5D')
        filled = self.ts.fillna(method='pad')
        assert_series_equal(rs, filled / filled.shift(freq='5D') - 1)

    def test_pct_change_shift_over_nas(self):
        s = Series([1., 1.5, np.nan, 2.5, 3.])

        chg = s.pct_change()
        expected = Series([np.nan, 0.5, np.nan, 2.5 / 1.5 - 1, .2])
        assert_series_equal(chg, expected)

    def test_autocorr(self):
        # Just run the function
        corr1 = self.ts.autocorr()

        # Now run it with the lag parameter
        corr2 = self.ts.autocorr(lag=1)

        # corr() with lag needs Series of at least length 2
        if len(self.ts) <= 2:
            self.assertTrue(np.isnan(corr1))
            self.assertTrue(np.isnan(corr2))
        else:
            self.assertEqual(corr1, corr2)

        # Choose a random lag between 1 and length of Series - 2
        # and compare the result with the Series corr() function
        n = 1 + np.random.randint(max(1, len(self.ts) - 2))
        corr1 = self.ts.corr(self.ts.shift(n))
        corr2 = self.ts.autocorr(lag=n)

        # corr() with lag needs Series of at least length 2
        if len(self.ts) <= 2:
            self.assertTrue(np.isnan(corr1))
            self.assertTrue(np.isnan(corr2))
        else:
            self.assertEqual(corr1, corr2)

    def test_first_last_valid(self):
        ts = self.ts.copy()
        ts[:5] = np.NaN

        index = ts.first_valid_index()
        self.assertEqual(index, ts.index[5])

        ts[-5:] = np.NaN
        index = ts.last_valid_index()
        self.assertEqual(index, ts.index[-6])

        ts[:] = np.nan
        self.assertIsNone(ts.last_valid_index())
        self.assertIsNone(ts.first_valid_index())

        ser = Series([], index=[])
        self.assertIsNone(ser.last_valid_index())
        self.assertIsNone(ser.first_valid_index())

    def test_mpl_compat_hack(self):
        result = self.ts[:, np.newaxis]
        expected = self.ts.values[:, np.newaxis]
        assert_almost_equal(result, expected)

#------------------------------------------------------------------------------
# GroupBy

    def test_select(self):
        n = len(self.ts)
        result = self.ts.select(lambda x: x >= self.ts.index[n // 2])
        expected = self.ts.reindex(self.ts.index[n // 2:])
        assert_series_equal(result, expected)

        result = self.ts.select(lambda x: x.weekday() == 2)
        expected = self.ts[self.ts.index.weekday == 2]
        assert_series_equal(result, expected)

#------------------------------------------------------------------------------
# Misc not safe for sparse

    def test_dropna_preserve_name(self):
        self.ts[:5] = np.nan
        result = self.ts.dropna()
        self.assertEqual(result.name, self.ts.name)
        name = self.ts.name
        ts = self.ts.copy()
        ts.dropna(inplace=True)
        self.assertEqual(ts.name, name)

    def test_numpy_unique(self):
        # it works!
        result = np.unique(self.ts)

    def test_concat_empty_series_dtypes_roundtrips(self):

        # round-tripping with self & like self
        dtypes = map(np.dtype,['float64','int8','uint8','bool','m8[ns]','M8[ns]'])

        for dtype in dtypes:
            self.assertEqual(pd.concat([Series(dtype=dtype)]).dtype, dtype)
            self.assertEqual(pd.concat([Series(dtype=dtype),
                                        Series(dtype=dtype)]).dtype, dtype)

        def int_result_type(dtype, dtype2):
            typs = set([dtype.kind,dtype2.kind])
            if not len(typs-set(['i','u','b'])) and (dtype.kind == 'i' or dtype2.kind == 'i'):
                return 'i'
            elif not len(typs-set(['u','b'])) and (dtype.kind == 'u' or dtype2.kind == 'u'):
                 return 'u'
            return None

        def float_result_type(dtype, dtype2):
            typs = set([dtype.kind,dtype2.kind])
            if not len(typs-set(['f','i','u'])) and (dtype.kind == 'f' or dtype2.kind == 'f'):
                return 'f'
            return None

        def get_result_type(dtype, dtype2):
            result = float_result_type(dtype, dtype2)
            if result is not None:
                return result
            result = int_result_type(dtype, dtype2)
            if result is not None:
                return result
            return 'O'

        for dtype in dtypes:
            for dtype2 in dtypes:
                if dtype == dtype2:
                    continue

                expected = get_result_type(dtype, dtype2)
                result = pd.concat([Series(dtype=dtype),
                                    Series(dtype=dtype2)]).dtype
                self.assertEqual(result.kind, expected)

    def test_concat_empty_series_dtypes(self):

        # bools
        self.assertEqual(pd.concat([Series(dtype=np.bool_),
                                    Series(dtype=np.int32)]).dtype, np.int32)
        self.assertEqual(pd.concat([Series(dtype=np.bool_),
                                    Series(dtype=np.float32)]).dtype, np.object_)

        # datetimelike
        self.assertEqual(pd.concat([Series(dtype='m8[ns]'),
                                    Series(dtype=np.bool)]).dtype, np.object_)
        self.assertEqual(pd.concat([Series(dtype='m8[ns]'),
                                    Series(dtype=np.int64)]).dtype, np.object_)
        self.assertEqual(pd.concat([Series(dtype='M8[ns]'),
                                    Series(dtype=np.bool)]).dtype, np.object_)
        self.assertEqual(pd.concat([Series(dtype='M8[ns]'),
                                    Series(dtype=np.int64)]).dtype, np.object_)
        self.assertEqual(pd.concat([Series(dtype='M8[ns]'),
                                    Series(dtype=np.bool_),
                                    Series(dtype=np.int64)]).dtype, np.object_)

        # categorical
        self.assertEqual(pd.concat([Series(dtype='category'),
                                    Series(dtype='category')]).dtype, 'category')
        self.assertEqual(pd.concat([Series(dtype='category'),
                                    Series(dtype='float64')]).dtype, np.object_)
        self.assertEqual(pd.concat([Series(dtype='category'),
                                    Series(dtype='object')]).dtype, 'category')

        # sparse
        result = pd.concat([Series(dtype='float64').to_sparse(),
                            Series(dtype='float64').to_sparse()])
        self.assertEqual(result.dtype,np.float64)
        self.assertEqual(result.ftype,'float64:sparse')

        result = pd.concat([Series(dtype='float64').to_sparse(),
                            Series(dtype='float64')])
        self.assertEqual(result.dtype,np.float64)
        self.assertEqual(result.ftype,'float64:sparse')

        result = pd.concat([Series(dtype='float64').to_sparse(),
                            Series(dtype='object')])
        self.assertEqual(result.dtype,np.object_)
        self.assertEqual(result.ftype,'object:dense')

    def test_searchsorted_numeric_dtypes_scalar(self):
        s = Series([1, 2, 90, 1000, 3e9])
        r = s.searchsorted(30)
        e = 2
        tm.assert_equal(r, e)

        r = s.searchsorted([30])
        e = np.array([2])
        tm.assert_numpy_array_equal(r, e)

    def test_searchsorted_numeric_dtypes_vector(self):
        s = Series([1, 2, 90, 1000, 3e9])
        r = s.searchsorted([91, 2e6])
        e = np.array([3, 4])
        tm.assert_numpy_array_equal(r, e)

    def test_search_sorted_datetime64_scalar(self):
        s = Series(pd.date_range('20120101', periods=10, freq='2D'))
        v = pd.Timestamp('20120102')
        r = s.searchsorted(v)
        e = 1
        tm.assert_equal(r, e)

    def test_search_sorted_datetime64_list(self):
        s = Series(pd.date_range('20120101', periods=10, freq='2D'))
        v = [pd.Timestamp('20120102'), pd.Timestamp('20120104')]
        r = s.searchsorted(v)
        e = np.array([1, 2])
        tm.assert_numpy_array_equal(r, e)

    def test_searchsorted_sorter(self):
        # GH8490
        s = Series([3, 1, 2])
        r = s.searchsorted([0, 3], sorter=np.argsort(s))
        e = np.array([0, 2])
        tm.assert_numpy_array_equal(r, e)

    def test_to_frame_expanddim(self):
        # GH 9762

        class SubclassedSeries(Series):
            @property
            def _constructor_expanddim(self):
                return SubclassedFrame

        class SubclassedFrame(DataFrame):
            pass

        s = SubclassedSeries([1, 2, 3], name='X')
        result = s.to_frame()
        self.assertTrue(isinstance(result, SubclassedFrame))
        expected = SubclassedFrame({'X': [1, 2, 3]})
        assert_frame_equal(result, expected)


class TestSeriesNonUnique(tm.TestCase):

    _multiprocess_can_split_ = True

    def setUp(self):
        pass

    def test_basic_indexing(self):
        s = Series(np.random.randn(5), index=['a', 'b', 'a', 'a', 'b'])

        self.assertRaises(IndexError, s.__getitem__, 5)
        self.assertRaises(IndexError, s.__setitem__, 5, 0)

        self.assertRaises(KeyError, s.__getitem__, 'c')

        s = s.sort_index()

        self.assertRaises(IndexError, s.__getitem__, 5)
        self.assertRaises(IndexError, s.__setitem__, 5, 0)


    def test_int_indexing(self):
        s = Series(np.random.randn(6), index=[0, 0, 1, 1, 2, 2])

        self.assertRaises(KeyError, s.__getitem__, 5)

        self.assertRaises(KeyError, s.__getitem__, 'c')

        # not monotonic
        s = Series(np.random.randn(6), index=[2, 2, 0, 0, 1, 1])

        self.assertRaises(KeyError, s.__getitem__, 5)

        self.assertRaises(KeyError, s.__getitem__, 'c')

    def test_datetime_indexing(self):
        from pandas import date_range

        index = date_range('1/1/2000', '1/7/2000')
        index = index.repeat(3)

        s = Series(len(index), index=index)
        stamp = Timestamp('1/8/2000')

        self.assertRaises(KeyError, s.__getitem__, stamp)
        s[stamp] = 0
        self.assertEqual(s[stamp], 0)

        # not monotonic
        s = Series(len(index), index=index)
        s = s[::-1]

        self.assertRaises(KeyError, s.__getitem__, stamp)
        s[stamp] = 0
        self.assertEqual(s[stamp], 0)

    def test_reset_index(self):
        df = tm.makeDataFrame()[:5]
        ser = df.stack()
        ser.index.names = ['hash', 'category']

        ser.name = 'value'
        df = ser.reset_index()
        self.assertIn('value', df)

        df = ser.reset_index(name='value2')
        self.assertIn('value2', df)

        # check inplace
        s = ser.reset_index(drop=True)
        s2 = ser
        s2.reset_index(drop=True, inplace=True)
        assert_series_equal(s, s2)

        # level
        index = MultiIndex(levels=[['bar'], ['one', 'two', 'three'], [0, 1]],
                           labels=[[0, 0, 0, 0, 0, 0],
                                   [0, 1, 2, 0, 1, 2],
                                   [0, 1, 0, 1, 0, 1]])
        s = Series(np.random.randn(6), index=index)
        rs = s.reset_index(level=1)
        self.assertEqual(len(rs.columns), 2)

        rs = s.reset_index(level=[0, 2], drop=True)
        self.assertTrue(rs.index.equals(Index(index.get_level_values(1))))
        tm.assertIsInstance(rs, Series)

    def test_set_index_makes_timeseries(self):
        idx = tm.makeDateIndex(10)

        s = Series(lrange(10))
        s.index = idx

        with tm.assert_produces_warning(FutureWarning):
            self.assertTrue(s.is_time_series == True)
        self.assertTrue(s.index.is_all_dates == True)

    def test_timeseries_coercion(self):
        idx = tm.makeDateIndex(10000)
        ser = Series(np.random.randn(len(idx)), idx.astype(object))
        with tm.assert_produces_warning(FutureWarning):
            self.assertTrue(ser.is_time_series)
        self.assertTrue(ser.index.is_all_dates)
        self.assertIsInstance(ser.index, DatetimeIndex)

    def test_replace(self):
        N = 100
        ser = Series(np.fabs(np.random.randn(N)), tm.makeDateIndex(N),
                     dtype=object)
        ser[:5] = np.nan
        ser[6:10] = 'foo'
        ser[20:30] = 'bar'

        # replace list with a single value
        rs = ser.replace([np.nan, 'foo', 'bar'], -1)

        self.assertTrue((rs[:5] == -1).all())
        self.assertTrue((rs[6:10] == -1).all())
        self.assertTrue((rs[20:30] == -1).all())
        self.assertTrue((isnull(ser[:5])).all())

        # replace with different values
        rs = ser.replace({np.nan: -1, 'foo': -2, 'bar': -3})

        self.assertTrue((rs[:5] == -1).all())
        self.assertTrue((rs[6:10] == -2).all())
        self.assertTrue((rs[20:30] == -3).all())
        self.assertTrue((isnull(ser[:5])).all())

        # replace with different values with 2 lists
        rs2 = ser.replace([np.nan, 'foo', 'bar'], [-1, -2, -3])
        assert_series_equal(rs, rs2)

        # replace inplace
        ser.replace([np.nan, 'foo', 'bar'], -1, inplace=True)
        self.assertTrue((ser[:5] == -1).all())
        self.assertTrue((ser[6:10] == -1).all())
        self.assertTrue((ser[20:30] == -1).all())

    def test_repeat(self):
        s = Series(np.random.randn(3), index=['a', 'b', 'c'])

        reps = s.repeat(5)
        exp = Series(s.values.repeat(5), index=s.index.values.repeat(5))
        assert_series_equal(reps, exp)

        to_rep = [2, 3, 4]
        reps = s.repeat(to_rep)
        exp = Series(s.values.repeat(to_rep),
                     index=s.index.values.repeat(to_rep))
        assert_series_equal(reps, exp)

    def test_unique_data_ownership(self):
        # it works! #1807
        Series(Series(["a", "c", "b"]).unique()).sort_values()

    def test_datetime_timedelta_quantiles(self):
        # covers #9694
        self.assertTrue(pd.isnull(Series([],dtype='M8[ns]').quantile(.5)))
        self.assertTrue(pd.isnull(Series([],dtype='m8[ns]').quantile(.5)))

    def test_empty_timeseries_redections_return_nat(self):
        # covers #11245
        for dtype in ('m8[ns]', 'm8[ns]', 'M8[ns]', 'M8[ns, UTC]'):
            self.assertIs(Series([], dtype=dtype).min(), pd.NaT)
            self.assertIs(Series([], dtype=dtype).max(), pd.NaT)

if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
