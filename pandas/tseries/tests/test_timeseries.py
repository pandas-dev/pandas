# pylint: disable-msg=E1101,W0612
from datetime import datetime, time, timedelta, date
import sys
import os
import operator

from distutils.version import LooseVersion

import nose

import numpy as np
randn = np.random.randn

from pandas import (Index, Series, TimeSeries, DataFrame,
                    isnull, date_range, Timestamp, Period, DatetimeIndex,
                    Int64Index, to_datetime, bdate_range)

from pandas.core.daterange import DateRange
import pandas.core.datetools as datetools
import pandas.tseries.offsets as offsets
import pandas.tseries.tools as tools
import pandas.tseries.frequencies as fmod
import pandas as pd

from pandas.util.testing import assert_series_equal, assert_almost_equal
import pandas.util.testing as tm

from pandas.tslib import NaT, iNaT
import pandas.lib as lib
import pandas.tslib as tslib

import pandas.index as _index

from pandas.compat import range, long, StringIO, lrange, lmap, zip, product
import pandas.core.datetools as dt
from numpy.random import rand
from numpy.testing import assert_array_equal
from pandas.util.testing import assert_frame_equal
import pandas.compat as compat
import pandas.core.common as com
from pandas import concat
from pandas import _np_version_under1p7

from numpy.testing.decorators import slow


def _skip_if_no_pytz():
    try:
        import pytz
    except ImportError:
        raise nose.SkipTest("pytz not installed")

def _skip_if_has_locale():
    import locale
    lang, _ = locale.getlocale()
    if lang is not None:
        raise nose.SkipTest("Specific locale is set {0}".format(lang))

class TestTimeSeriesDuplicates(tm.TestCase):
    _multiprocess_can_split_ = True

    def setUp(self):
        dates = [datetime(2000, 1, 2), datetime(2000, 1, 2),
                 datetime(2000, 1, 2), datetime(2000, 1, 3),
                 datetime(2000, 1, 3), datetime(2000, 1, 3),
                 datetime(2000, 1, 4), datetime(2000, 1, 4),
                 datetime(2000, 1, 4), datetime(2000, 1, 5)]

        self.dups = Series(np.random.randn(len(dates)), index=dates)

    def test_constructor(self):
        tm.assert_isinstance(self.dups, TimeSeries)
        tm.assert_isinstance(self.dups.index, DatetimeIndex)

    def test_is_unique_monotonic(self):
        self.assert_(not self.dups.index.is_unique)

    def test_index_unique(self):
        uniques = self.dups.index.unique()
        self.assert_(uniques.dtype == 'M8[ns]')  # sanity

        # #2563
        self.assertTrue(isinstance(uniques, DatetimeIndex))

        dups_local = self.dups.index.tz_localize('US/Eastern')
        dups_local.name = 'foo'
        result = dups_local.unique()
        self.assertTrue(result.tz is not None)
        self.assertEquals(result.name, 'foo')

    def test_index_dupes_contains(self):
        d = datetime(2011, 12, 5, 20, 30)
        ix = DatetimeIndex([d, d])
        self.assertTrue(d in ix)

    def test_duplicate_dates_indexing(self):
        ts = self.dups

        uniques = ts.index.unique()
        for date in uniques:
            result = ts[date]

            mask = ts.index == date
            total = (ts.index == date).sum()
            expected = ts[mask]
            if total > 1:
                assert_series_equal(result, expected)
            else:
                assert_almost_equal(result, expected[0])

            cp = ts.copy()
            cp[date] = 0
            expected = Series(np.where(mask, 0, ts), index=ts.index)
            assert_series_equal(cp, expected)

        self.assertRaises(KeyError, ts.__getitem__, datetime(2000, 1, 6))

        # new index
        ts[datetime(2000,1,6)] = 0
        self.assert_(ts[datetime(2000,1,6)] == 0)

    def test_range_slice(self):
        idx = DatetimeIndex(['1/1/2000', '1/2/2000', '1/2/2000', '1/3/2000',
                             '1/4/2000'])

        ts = Series(np.random.randn(len(idx)), index=idx)

        result = ts['1/2/2000':]
        expected = ts[1:]
        assert_series_equal(result, expected)

        result = ts['1/2/2000':'1/3/2000']
        expected = ts[1:4]
        assert_series_equal(result, expected)

    def test_groupby_average_dup_values(self):
        result = self.dups.groupby(level=0).mean()
        expected = self.dups.groupby(self.dups.index).mean()
        assert_series_equal(result, expected)

    def test_indexing_over_size_cutoff(self):
        import datetime
        # #1821

        old_cutoff = _index._SIZE_CUTOFF
        try:
            _index._SIZE_CUTOFF = 1000

            # create large list of non periodic datetime
            dates = []
            sec = datetime.timedelta(seconds=1)
            half_sec = datetime.timedelta(microseconds=500000)
            d = datetime.datetime(2011, 12, 5, 20, 30)
            n = 1100
            for i in range(n):
                dates.append(d)
                dates.append(d + sec)
                dates.append(d + sec + half_sec)
                dates.append(d + sec + sec + half_sec)
                d += 3 * sec

            # duplicate some values in the list
            duplicate_positions = np.random.randint(0, len(dates) - 1, 20)
            for p in duplicate_positions:
                dates[p + 1] = dates[p]

            df = DataFrame(np.random.randn(len(dates), 4),
                           index=dates,
                           columns=list('ABCD'))

            pos = n * 3
            timestamp = df.index[pos]
            self.assert_(timestamp in df.index)

            # it works!
            df.ix[timestamp]
            self.assert_(len(df.ix[[timestamp]]) > 0)
        finally:
            _index._SIZE_CUTOFF = old_cutoff

    def test_indexing_unordered(self):

        # GH 2437
        rng = date_range(start='2011-01-01', end='2011-01-15')
        ts  = Series(randn(len(rng)), index=rng)
        ts2 = concat([ts[0:4],ts[-4:],ts[4:-4]])

        for t in ts.index:
            s = str(t)
            expected = ts[t]
            result = ts2[t]
            self.assertTrue(expected == result)

        # GH 3448 (ranges)
        def compare(slobj):
            result = ts2[slobj].copy()
            result = result.sort_index()
            expected = ts[slobj]
            assert_series_equal(result,expected)

        compare(slice('2011-01-01','2011-01-15'))
        compare(slice('2010-12-30','2011-01-15'))
        compare(slice('2011-01-01','2011-01-16'))

        # partial ranges
        compare(slice('2011-01-01','2011-01-6'))
        compare(slice('2011-01-06','2011-01-8'))
        compare(slice('2011-01-06','2011-01-12'))

        # single values
        result = ts2['2011'].sort_index()
        expected = ts['2011']
        assert_series_equal(result,expected)

        # diff freq
        rng = date_range(datetime(2005, 1, 1), periods=20, freq='M')
        ts = Series(np.arange(len(rng)), index=rng)
        ts = ts.take(np.random.permutation(20))

        result = ts['2005']
        for t in result.index:
            self.assertTrue(t.year == 2005)

    def test_indexing(self):

        idx = date_range("2001-1-1", periods=20, freq='M')
        ts = Series(np.random.rand(len(idx)),index=idx)

        # getting

        # GH 3070, make sure semantics work on Series/Frame
        expected = ts['2001']

        df = DataFrame(dict(A = ts))
        result = df['2001']['A']
        assert_series_equal(expected,result)

        # setting
        ts['2001'] = 1
        expected = ts['2001']

        df.loc['2001','A'] = 1

        result = df['2001']['A']
        assert_series_equal(expected,result)

        # GH3546 (not including times on the last day)
        idx = date_range(start='2013-05-31 00:00', end='2013-05-31 23:00', freq='H')
        ts  = Series(lrange(len(idx)), index=idx)
        expected = ts['2013-05']
        assert_series_equal(expected,ts)

        idx = date_range(start='2013-05-31 00:00', end='2013-05-31 23:59', freq='S')
        ts  = Series(lrange(len(idx)), index=idx)
        expected = ts['2013-05']
        assert_series_equal(expected,ts)

        idx = [ Timestamp('2013-05-31 00:00'), Timestamp(datetime(2013,5,31,23,59,59,999999))]
        ts  = Series(lrange(len(idx)), index=idx)
        expected = ts['2013']
        assert_series_equal(expected,ts)

        # GH 3925, indexing with a seconds resolution string / datetime object
        df = DataFrame(randn(5,5),columns=['open','high','low','close','volume'],index=date_range('2012-01-02 18:01:00',periods=5,tz='US/Central',freq='s'))
        expected = df.loc[[df.index[2]]]
        result = df['2012-01-02 18:01:02']
        assert_frame_equal(result,expected)

        # this is a single date, so will raise
        self.assertRaises(KeyError, df.__getitem__, df.index[2],)


def assert_range_equal(left, right):
    assert(left.equals(right))
    assert(left.freq == right.freq)
    assert(left.tz == right.tz)


class TestTimeSeries(tm.TestCase):
    _multiprocess_can_split_ = True

    def test_is_(self):
        dti = DatetimeIndex(start='1/1/2005', end='12/1/2005', freq='M')
        self.assertTrue(dti.is_(dti))
        self.assertTrue(dti.is_(dti.view()))
        self.assertFalse(dti.is_(dti.copy()))

    def test_dti_slicing(self):
        dti = DatetimeIndex(start='1/1/2005', end='12/1/2005', freq='M')
        dti2 = dti[[1, 3, 5]]

        v1 = dti2[0]
        v2 = dti2[1]
        v3 = dti2[2]

        self.assertEquals(v1, Timestamp('2/28/2005'))
        self.assertEquals(v2, Timestamp('4/30/2005'))
        self.assertEquals(v3, Timestamp('6/30/2005'))

        # don't carry freq through irregular slicing
        self.assert_(dti2.freq is None)

    def test_pass_datetimeindex_to_index(self):
        # Bugs in #1396

        rng = date_range('1/1/2000', '3/1/2000')
        idx = Index(rng, dtype=object)

        expected = Index(rng.to_pydatetime(), dtype=object)

        self.assert_(np.array_equal(idx.values, expected.values))

    def test_contiguous_boolean_preserve_freq(self):
        rng = date_range('1/1/2000', '3/1/2000', freq='B')

        mask = np.zeros(len(rng), dtype=bool)
        mask[10:20] = True

        masked = rng[mask]
        expected = rng[10:20]
        self.assert_(expected.freq is not None)
        assert_range_equal(masked, expected)

        mask[22] = True
        masked = rng[mask]
        self.assert_(masked.freq is None)

    def test_getitem_median_slice_bug(self):
        index = date_range('20090415', '20090519', freq='2B')
        s = Series(np.random.randn(13), index=index)

        indexer = [slice(6, 7, None)]
        result = s[indexer]
        expected = s[indexer[0]]
        assert_series_equal(result, expected)

    def test_series_box_timestamp(self):
        rng = date_range('20090415', '20090519', freq='B')
        s = Series(rng)

        tm.assert_isinstance(s[5], Timestamp)

        rng = date_range('20090415', '20090519', freq='B')
        s = Series(rng, index=rng)
        tm.assert_isinstance(s[5], Timestamp)

        tm.assert_isinstance(s.iget_value(5), Timestamp)

    def test_date_range_ambiguous_arguments(self):
        # #2538
        start = datetime(2011, 1, 1, 5, 3, 40)
        end = datetime(2011, 1, 1, 8, 9, 40)

        self.assertRaises(ValueError, date_range, start, end,
                          freq='s', periods=10)

    def test_timestamp_to_datetime(self):
        _skip_if_no_pytz()
        rng = date_range('20090415', '20090519',
                         tz='US/Eastern')

        stamp = rng[0]
        dtval = stamp.to_pydatetime()
        self.assertEquals(stamp, dtval)
        self.assertEquals(stamp.tzinfo, dtval.tzinfo)

    def test_index_convert_to_datetime_array(self):
        _skip_if_no_pytz()

        def _check_rng(rng):
            converted = rng.to_pydatetime()
            tm.assert_isinstance(converted, np.ndarray)
            for x, stamp in zip(converted, rng):
                tm.assert_isinstance(x, datetime)
                self.assertEquals(x, stamp.to_pydatetime())
                self.assertEquals(x.tzinfo, stamp.tzinfo)

        rng = date_range('20090415', '20090519')
        rng_eastern = date_range('20090415', '20090519', tz='US/Eastern')
        rng_utc = date_range('20090415', '20090519', tz='utc')

        _check_rng(rng)
        _check_rng(rng_eastern)
        _check_rng(rng_utc)

    def test_ctor_str_intraday(self):
        rng = DatetimeIndex(['1-1-2000 00:00:01'])
        self.assert_(rng[0].second == 1)

    def test_series_ctor_plus_datetimeindex(self):
        rng = date_range('20090415', '20090519', freq='B')
        data = dict((k, 1) for k in rng)

        result = Series(data, index=rng)
        self.assert_(result.index is rng)

    def test_series_pad_backfill_limit(self):
        index = np.arange(10)
        s = Series(np.random.randn(10), index=index)

        result = s[:2].reindex(index, method='pad', limit=5)

        expected = s[:2].reindex(index).fillna(method='pad')
        expected[-3:] = np.nan
        assert_series_equal(result, expected)

        result = s[-2:].reindex(index, method='backfill', limit=5)

        expected = s[-2:].reindex(index).fillna(method='backfill')
        expected[:3] = np.nan
        assert_series_equal(result, expected)

    def test_series_fillna_limit(self):
        index = np.arange(10)
        s = Series(np.random.randn(10), index=index)

        result = s[:2].reindex(index)
        result = result.fillna(method='pad', limit=5)

        expected = s[:2].reindex(index).fillna(method='pad')
        expected[-3:] = np.nan
        assert_series_equal(result, expected)

        result = s[-2:].reindex(index)
        result = result.fillna(method='bfill', limit=5)

        expected = s[-2:].reindex(index).fillna(method='backfill')
        expected[:3] = np.nan
        assert_series_equal(result, expected)

    def test_frame_pad_backfill_limit(self):
        index = np.arange(10)
        df = DataFrame(np.random.randn(10, 4), index=index)

        result = df[:2].reindex(index, method='pad', limit=5)

        expected = df[:2].reindex(index).fillna(method='pad')
        expected.values[-3:] = np.nan
        tm.assert_frame_equal(result, expected)

        result = df[-2:].reindex(index, method='backfill', limit=5)

        expected = df[-2:].reindex(index).fillna(method='backfill')
        expected.values[:3] = np.nan
        tm.assert_frame_equal(result, expected)

    def test_frame_fillna_limit(self):
        index = np.arange(10)
        df = DataFrame(np.random.randn(10, 4), index=index)

        result = df[:2].reindex(index)
        result = result.fillna(method='pad', limit=5)

        expected = df[:2].reindex(index).fillna(method='pad')
        expected.values[-3:] = np.nan
        tm.assert_frame_equal(result, expected)

        result = df[-2:].reindex(index)
        result = result.fillna(method='backfill', limit=5)

        expected = df[-2:].reindex(index).fillna(method='backfill')
        expected.values[:3] = np.nan
        tm.assert_frame_equal(result, expected)

    def test_frame_setitem_timestamp(self):
        # 2155
        columns = DatetimeIndex(start='1/1/2012', end='2/1/2012',
                                freq=datetools.bday)
        index = lrange(10)
        data = DataFrame(columns=columns, index=index)
        t = datetime(2012, 11, 1)
        ts = Timestamp(t)
        data[ts] = np.nan  # works

    def test_sparse_series_fillna_limit(self):
        index = np.arange(10)
        s = Series(np.random.randn(10), index=index)

        ss = s[:2].reindex(index).to_sparse()
        result = ss.fillna(method='pad', limit=5)
        expected = ss.fillna(method='pad', limit=5)
        expected = expected.to_dense()
        expected[-3:] = np.nan
        expected = expected.to_sparse()
        assert_series_equal(result, expected)

        ss = s[-2:].reindex(index).to_sparse()
        result = ss.fillna(method='backfill', limit=5)
        expected = ss.fillna(method='backfill')
        expected = expected.to_dense()
        expected[:3] = np.nan
        expected = expected.to_sparse()
        assert_series_equal(result, expected)

    def test_sparse_series_pad_backfill_limit(self):
        index = np.arange(10)
        s = Series(np.random.randn(10), index=index)
        s = s.to_sparse()

        result = s[:2].reindex(index, method='pad', limit=5)
        expected = s[:2].reindex(index).fillna(method='pad')
        expected = expected.to_dense()
        expected[-3:] = np.nan
        expected = expected.to_sparse()
        assert_series_equal(result, expected)

        result = s[-2:].reindex(index, method='backfill', limit=5)
        expected = s[-2:].reindex(index).fillna(method='backfill')
        expected = expected.to_dense()
        expected[:3] = np.nan
        expected = expected.to_sparse()
        assert_series_equal(result, expected)

    def test_sparse_frame_pad_backfill_limit(self):
        index = np.arange(10)
        df = DataFrame(np.random.randn(10, 4), index=index)
        sdf = df.to_sparse()

        result = sdf[:2].reindex(index, method='pad', limit=5)

        expected = sdf[:2].reindex(index).fillna(method='pad')
        expected = expected.to_dense()
        expected.values[-3:] = np.nan
        expected = expected.to_sparse()
        tm.assert_frame_equal(result, expected)

        result = sdf[-2:].reindex(index, method='backfill', limit=5)

        expected = sdf[-2:].reindex(index).fillna(method='backfill')
        expected = expected.to_dense()
        expected.values[:3] = np.nan
        expected = expected.to_sparse()
        tm.assert_frame_equal(result, expected)

    def test_sparse_frame_fillna_limit(self):
        index = np.arange(10)
        df = DataFrame(np.random.randn(10, 4), index=index)
        sdf = df.to_sparse()

        result = sdf[:2].reindex(index)
        result = result.fillna(method='pad', limit=5)

        expected = sdf[:2].reindex(index).fillna(method='pad')
        expected = expected.to_dense()
        expected.values[-3:] = np.nan
        expected = expected.to_sparse()
        tm.assert_frame_equal(result, expected)

        result = sdf[-2:].reindex(index)
        result = result.fillna(method='backfill', limit=5)

        expected = sdf[-2:].reindex(index).fillna(method='backfill')
        expected = expected.to_dense()
        expected.values[:3] = np.nan
        expected = expected.to_sparse()
        tm.assert_frame_equal(result, expected)

    def test_pad_require_monotonicity(self):
        rng = date_range('1/1/2000', '3/1/2000', freq='B')

        rng2 = rng[::2][::-1]

        self.assertRaises(ValueError, rng2.get_indexer, rng,
                          method='pad')

    def test_frame_ctor_datetime64_column(self):
        rng = date_range('1/1/2000 00:00:00', '1/1/2000 1:59:50',
                         freq='10s')
        dates = np.asarray(rng)

        df = DataFrame({'A': np.random.randn(len(rng)), 'B': dates})
        self.assert_(np.issubdtype(df['B'].dtype, np.dtype('M8[ns]')))

    def test_frame_add_datetime64_column(self):
        rng = date_range('1/1/2000 00:00:00', '1/1/2000 1:59:50',
                         freq='10s')
        df = DataFrame(index=np.arange(len(rng)))

        df['A'] = rng
        self.assert_(np.issubdtype(df['A'].dtype, np.dtype('M8[ns]')))

    def test_frame_datetime64_pre1900_repr(self):
        df = DataFrame({'year': date_range('1/1/1700', periods=50,
                                           freq='A-DEC')})
        # it works!
        repr(df)

    def test_frame_add_datetime64_col_other_units(self):
        n = 100

        units = ['h', 'm', 's', 'ms', 'D', 'M', 'Y']

        ns_dtype = np.dtype('M8[ns]')

        for unit in units:
            dtype = np.dtype('M8[%s]' % unit)
            vals = np.arange(n, dtype=np.int64).view(dtype)

            df = DataFrame({'ints': np.arange(n)}, index=np.arange(n))
            df[unit] = vals

            ex_vals = to_datetime(vals.astype('O'))

            self.assert_(df[unit].dtype == ns_dtype)
            self.assert_((df[unit].values == ex_vals).all())

        # Test insertion into existing datetime64 column
        df = DataFrame({'ints': np.arange(n)}, index=np.arange(n))
        df['dates'] = np.arange(n, dtype=np.int64).view(ns_dtype)

        for unit in units:
            dtype = np.dtype('M8[%s]' % unit)
            vals = np.arange(n, dtype=np.int64).view(dtype)

            tmp = df.copy()

            tmp['dates'] = vals
            ex_vals = to_datetime(vals.astype('O'))

            self.assert_((tmp['dates'].values == ex_vals).all())

    def test_to_datetime_unit(self):

        epoch = 1370745748
        s = Series([ epoch + t for t in range(20) ])
        result = to_datetime(s,unit='s')
        expected = Series([ Timestamp('2013-06-09 02:42:28') + timedelta(seconds=t) for t in range(20) ])
        assert_series_equal(result,expected)

        s = Series([ epoch + t for t in range(20) ]).astype(float)
        result = to_datetime(s,unit='s')
        expected = Series([ Timestamp('2013-06-09 02:42:28') + timedelta(seconds=t) for t in range(20) ])
        assert_series_equal(result,expected)

        s = Series([ epoch + t for t in range(20) ] + [iNaT])
        result = to_datetime(s,unit='s')
        expected = Series([ Timestamp('2013-06-09 02:42:28') + timedelta(seconds=t) for t in range(20) ] + [NaT])
        assert_series_equal(result,expected)

        s = Series([ epoch + t for t in range(20) ] + [iNaT]).astype(float)
        result = to_datetime(s,unit='s')
        expected = Series([ Timestamp('2013-06-09 02:42:28') + timedelta(seconds=t) for t in range(20) ] + [NaT])
        assert_series_equal(result,expected)

        s = concat([Series([ epoch + t for t in range(20) ]).astype(float),Series([np.nan])],ignore_index=True)
        result = to_datetime(s,unit='s')
        expected = Series([ Timestamp('2013-06-09 02:42:28') + timedelta(seconds=t) for t in range(20) ] + [NaT])
        assert_series_equal(result,expected)

    def test_series_ctor_datetime64(self):
        rng = date_range('1/1/2000 00:00:00', '1/1/2000 1:59:50',
                         freq='10s')
        dates = np.asarray(rng)

        series = Series(dates)
        self.assert_(np.issubdtype(series.dtype, np.dtype('M8[ns]')))

    def test_index_cast_datetime64_other_units(self):
        arr = np.arange(0, 100, 10, dtype=np.int64).view('M8[D]')

        idx = Index(arr)

        self.assert_((idx.values == tslib.cast_to_nanoseconds(arr)).all())

    def test_index_astype_datetime64(self):
        idx = Index([datetime(2012, 1, 1)], dtype=object)

        if not _np_version_under1p7:
            raise nose.SkipTest("test only valid in numpy < 1.7")

        casted = idx.astype(np.dtype('M8[D]'))
        expected = DatetimeIndex(idx.values)
        tm.assert_isinstance(casted, DatetimeIndex)
        self.assert_(casted.equals(expected))

    def test_reindex_series_add_nat(self):
        rng = date_range('1/1/2000 00:00:00', periods=10, freq='10s')
        series = Series(rng)

        result = series.reindex(lrange(15))
        self.assert_(np.issubdtype(result.dtype, np.dtype('M8[ns]')))

        mask = result.isnull()
        self.assert_(mask[-5:].all())
        self.assert_(not mask[:-5].any())

    def test_reindex_frame_add_nat(self):
        rng = date_range('1/1/2000 00:00:00', periods=10, freq='10s')
        df = DataFrame({'A': np.random.randn(len(rng)), 'B': rng})

        result = df.reindex(lrange(15))
        self.assert_(np.issubdtype(result['B'].dtype, np.dtype('M8[ns]')))

        mask = com.isnull(result)['B']
        self.assert_(mask[-5:].all())
        self.assert_(not mask[:-5].any())

    def test_series_repr_nat(self):
        series = Series([0, 1000, 2000, iNaT], dtype='M8[ns]')

        result = repr(series)
        expected = ('0          1970-01-01 00:00:00\n'
                    '1   1970-01-01 00:00:00.000001\n'
                    '2   1970-01-01 00:00:00.000002\n'
                    '3                          NaT\n'
                    'dtype: datetime64[ns]')
        self.assertEquals(result, expected)

    def test_fillna_nat(self):
        series = Series([0, 1, 2, iNaT], dtype='M8[ns]')

        filled = series.fillna(method='pad')
        filled2 = series.fillna(value=series.values[2])

        expected = series.copy()
        expected.values[3] = expected.values[2]

        assert_series_equal(filled, expected)
        assert_series_equal(filled2, expected)

        df = DataFrame({'A': series})
        filled = df.fillna(method='pad')
        filled2 = df.fillna(value=series.values[2])
        expected = DataFrame({'A': expected})
        assert_frame_equal(filled, expected)
        assert_frame_equal(filled2, expected)

        series = Series([iNaT, 0, 1, 2], dtype='M8[ns]')

        filled = series.fillna(method='bfill')
        filled2 = series.fillna(value=series[1])

        expected = series.copy()
        expected[0] = expected[1]

        assert_series_equal(filled, expected)
        assert_series_equal(filled2, expected)

        df = DataFrame({'A': series})
        filled = df.fillna(method='bfill')
        filled2 = df.fillna(value=series[1])
        expected = DataFrame({'A': expected})
        assert_frame_equal(filled, expected)
        assert_frame_equal(filled2, expected)

    def test_string_na_nat_conversion(self):
        # GH #999, #858

        from pandas.compat import parse_date

        strings = np.array(['1/1/2000', '1/2/2000', np.nan,
                            '1/4/2000, 12:34:56'], dtype=object)

        expected = np.empty(4, dtype='M8[ns]')
        for i, val in enumerate(strings):
            if com.isnull(val):
                expected[i] = iNaT
            else:
                expected[i] = parse_date(val)

        result = tslib.array_to_datetime(strings)
        assert_almost_equal(result, expected)

        result2 = to_datetime(strings)
        tm.assert_isinstance(result2, DatetimeIndex)
        assert_almost_equal(result, result2)

        malformed = np.array(['1/100/2000', np.nan], dtype=object)
        result = to_datetime(malformed)
        assert_almost_equal(result, malformed)

        self.assertRaises(ValueError, to_datetime, malformed,
                          errors='raise')

        idx = ['a', 'b', 'c', 'd', 'e']
        series = Series(['1/1/2000', np.nan, '1/3/2000', np.nan,
                         '1/5/2000'], index=idx, name='foo')
        dseries = Series([to_datetime('1/1/2000'), np.nan,
                          to_datetime('1/3/2000'), np.nan,
                          to_datetime('1/5/2000')], index=idx, name='foo')

        result = to_datetime(series)
        dresult = to_datetime(dseries)

        expected = Series(np.empty(5, dtype='M8[ns]'), index=idx)
        for i in range(5):
            x = series[i]
            if isnull(x):
                expected[i] = iNaT
            else:
                expected[i] = to_datetime(x)

        assert_series_equal(result, expected)
        self.assertEquals(result.name, 'foo')

        assert_series_equal(dresult, expected)
        self.assertEquals(dresult.name, 'foo')

    def test_to_datetime_iso8601(self):
        result = to_datetime(["2012-01-01 00:00:00"])
        exp = Timestamp("2012-01-01 00:00:00")
        self.assert_(result[0] == exp)

        result = to_datetime(['20121001'])  # bad iso 8601
        exp = Timestamp('2012-10-01')
        self.assert_(result[0] == exp)

    def test_to_datetime_default(self):
        rs = to_datetime('2001')
        xp = datetime(2001, 1, 1)
        self.assert_(rs, xp)

        #### dayfirst is essentially broken
        #### to_datetime('01-13-2012', dayfirst=True)
        #### self.assertRaises(ValueError, to_datetime('01-13-2012', dayfirst=True))

    def test_to_datetime_on_datetime64_series(self):
        # #2699
        s = Series(date_range('1/1/2000', periods=10))

        result = to_datetime(s)
        self.assertEquals(result[0], s[0])

    def test_to_datetime_with_apply(self):
        # this is only locale tested with US/None locales
        _skip_if_has_locale()

        # GH 5195
        # with a format and coerce a single item to_datetime fails
        td = Series(['May 04', 'Jun 02', 'Dec 11'], index=[1,2,3])
        expected = pd.to_datetime(td, format='%b %y')
        result = td.apply(pd.to_datetime, format='%b %y')
        assert_series_equal(result, expected)

        td = pd.Series(['May 04', 'Jun 02', ''], index=[1,2,3])
        self.assertRaises(ValueError, lambda : pd.to_datetime(td,format='%b %y'))
        self.assertRaises(ValueError, lambda : td.apply(pd.to_datetime, format='%b %y'))
        expected = pd.to_datetime(td, format='%b %y', coerce=True)

        result = td.apply(lambda x: pd.to_datetime(x, format='%b %y', coerce=True))
        assert_series_equal(result, expected)

    def test_nat_vector_field_access(self):
        idx = DatetimeIndex(['1/1/2000', None, None, '1/4/2000'])

        fields = ['year', 'quarter', 'month', 'day', 'hour',
                  'minute', 'second', 'microsecond', 'nanosecond',
                  'week', 'dayofyear']
        for field in fields:
            result = getattr(idx, field)
            expected = [getattr(x, field) if x is not NaT else -1
                        for x in idx]
            self.assert_(np.array_equal(result, expected))

    def test_nat_scalar_field_access(self):
        fields = ['year', 'quarter', 'month', 'day', 'hour',
                  'minute', 'second', 'microsecond', 'nanosecond',
                  'week', 'dayofyear']
        for field in fields:
            result = getattr(NaT, field)
            self.assertEquals(result, -1)

        self.assertEquals(NaT.weekday(), -1)

    def test_to_datetime_types(self):

        # empty string
        result = to_datetime('')
        self.assert_(result is NaT)

        result = to_datetime(['', ''])
        self.assert_(isnull(result).all())

        # ints
        result = Timestamp(0)
        expected = to_datetime(0)
        self.assert_(result == expected)

        # GH 3888 (strings)
        expected = to_datetime(['2012'])[0]
        result = to_datetime('2012')
        self.assert_(result == expected)

        ### array = ['2012','20120101','20120101 12:01:01']
        array = ['20120101','20120101 12:01:01']
        expected = list(to_datetime(array))
        result = lmap(Timestamp,array)
        tm.assert_almost_equal(result,expected)

        ### currently fails ###
        ### result = Timestamp('2012')
        ### expected = to_datetime('2012')
        ### self.assert_(result == expected)

    def test_to_datetime_unprocessable_input(self):
        # GH 4928
        self.assert_(
            np.array_equal(
                to_datetime([1, '1']),
                np.array([1, '1'], dtype='O')
            )
        )
        self.assertRaises(TypeError, to_datetime, [1, '1'], errors='raise')

    def test_to_datetime_other_datetime64_units(self):
        # 5/25/2012
        scalar = np.int64(1337904000000000).view('M8[us]')
        as_obj = scalar.astype('O')

        index = DatetimeIndex([scalar])
        self.assertEquals(index[0], scalar.astype('O'))

        value = Timestamp(scalar)
        self.assertEquals(value, as_obj)

    def test_to_datetime_list_of_integers(self):
        rng = date_range('1/1/2000', periods=20)
        rng = DatetimeIndex(rng.values)

        ints = list(rng.asi8)

        result = DatetimeIndex(ints)

        self.assert_(rng.equals(result))

    def test_to_datetime_dt64s(self):
        in_bound_dts = [
            np.datetime64('2000-01-01'),
            np.datetime64('2000-01-02'),
        ]

        for dt in in_bound_dts:
            self.assertEqual(
                pd.to_datetime(dt),
                Timestamp(dt)
            )

        oob_dts = [
            np.datetime64('1000-01-01'),
            np.datetime64('5000-01-02'),
        ]

        for dt in oob_dts:
            self.assertRaises(ValueError, pd.to_datetime, dt, errors='raise')
            self.assertRaises(ValueError, tslib.Timestamp, dt)
            self.assert_(pd.to_datetime(dt, coerce=True) is NaT)

    def test_to_datetime_array_of_dt64s(self):
        dts = [
            np.datetime64('2000-01-01'),
            np.datetime64('2000-01-02'),
        ]

        # Assuming all datetimes are in bounds, to_datetime() returns
        # an array that is equal to Timestamp() parsing
        self.assert_(
            np.array_equal(
                pd.to_datetime(dts, box=False),
                np.array([Timestamp(x).asm8 for x in dts])
            )
        )

        # A list of datetimes where the last one is out of bounds
        dts_with_oob = dts + [np.datetime64('9999-01-01')]

        self.assertRaises(
            ValueError,
            pd.to_datetime,
            dts_with_oob,
            coerce=False,
            errors='raise'
        )

        self.assert_(
            np.array_equal(
                pd.to_datetime(dts_with_oob, box=False, coerce=True),
                np.array(
                    [
                        Timestamp(dts_with_oob[0]).asm8,
                        Timestamp(dts_with_oob[1]).asm8,
                        iNaT,
                    ],
                    dtype='M8'
                )
            )
        )

        # With coerce=False and errors='ignore', out of bounds datetime64s
        # are converted to their .item(), which depending on the version of
        # numpy is either a python datetime.datetime or datetime.date
        self.assert_(
            np.array_equal(
                pd.to_datetime(dts_with_oob, box=False, coerce=False),
                np.array(
                    [dt.item() for dt in dts_with_oob],
                    dtype='O'
                )
            )
        )

    def test_index_to_datetime(self):
        idx = Index(['1/1/2000', '1/2/2000', '1/3/2000'])

        result = idx.to_datetime()
        expected = DatetimeIndex(datetools.to_datetime(idx.values))
        self.assert_(result.equals(expected))

        today = datetime.today()
        idx = Index([today], dtype=object)
        result = idx.to_datetime()
        expected = DatetimeIndex([today])
        self.assert_(result.equals(expected))

    def test_to_datetime_freq(self):
        xp = bdate_range('2000-1-1', periods=10, tz='UTC')
        rs = xp.to_datetime()
        self.assert_(xp.freq == rs.freq)
        self.assert_(xp.tzinfo == rs.tzinfo)

    def test_range_misspecified(self):
        # GH #1095

        self.assertRaises(ValueError, date_range, '1/1/2000')
        self.assertRaises(ValueError, date_range, end='1/1/2000')
        self.assertRaises(ValueError, date_range, periods=10)

        self.assertRaises(ValueError, date_range, '1/1/2000', freq='H')
        self.assertRaises(ValueError, date_range, end='1/1/2000', freq='H')
        self.assertRaises(ValueError, date_range, periods=10, freq='H')

    def test_reasonable_keyerror(self):
        # GH #1062
        index = DatetimeIndex(['1/3/2000'])
        try:
            index.get_loc('1/1/2000')
        except KeyError as e:
            self.assert_('2000' in str(e))

    def test_reindex_with_datetimes(self):
        rng = date_range('1/1/2000', periods=20)
        ts = Series(np.random.randn(20), index=rng)

        result = ts.reindex(list(ts.index[5:10]))
        expected = ts[5:10]
        tm.assert_series_equal(result, expected)

        result = ts[list(ts.index[5:10])]
        tm.assert_series_equal(result, expected)

    def test_promote_datetime_date(self):
        rng = date_range('1/1/2000', periods=20)
        ts = Series(np.random.randn(20), index=rng)

        ts_slice = ts[5:]
        ts2 = ts_slice.copy()
        ts2.index = [x.date() for x in ts2.index]

        result = ts + ts2
        result2 = ts2 + ts
        expected = ts + ts[5:]
        assert_series_equal(result, expected)
        assert_series_equal(result2, expected)

        # test asfreq
        result = ts2.asfreq('4H', method='ffill')
        expected = ts[5:].asfreq('4H', method='ffill')
        assert_series_equal(result, expected)

        result = rng.get_indexer(ts2.index)
        expected = rng.get_indexer(ts_slice.index)
        self.assert_(np.array_equal(result, expected))

    def test_asfreq_normalize(self):
        rng = date_range('1/1/2000 09:30', periods=20)
        norm = date_range('1/1/2000', periods=20)
        vals = np.random.randn(20)
        ts = Series(vals, index=rng)

        result = ts.asfreq('D', normalize=True)
        norm = date_range('1/1/2000', periods=20)
        expected = Series(vals, index=norm)

        assert_series_equal(result, expected)

        vals = np.random.randn(20, 3)
        ts = DataFrame(vals, index=rng)

        result = ts.asfreq('D', normalize=True)
        expected = DataFrame(vals, index=norm)

        assert_frame_equal(result, expected)

    def test_date_range_gen_error(self):
        rng = date_range('1/1/2000 00:00', '1/1/2000 00:18', freq='5min')
        self.assertEquals(len(rng), 4)

    def test_first_subset(self):
        ts = _simple_ts('1/1/2000', '1/1/2010', freq='12h')
        result = ts.first('10d')
        self.assert_(len(result) == 20)

        ts = _simple_ts('1/1/2000', '1/1/2010')
        result = ts.first('10d')
        self.assert_(len(result) == 10)

        result = ts.first('3M')
        expected = ts[:'3/31/2000']
        assert_series_equal(result, expected)

        result = ts.first('21D')
        expected = ts[:21]
        assert_series_equal(result, expected)

        result = ts[:0].first('3M')
        assert_series_equal(result, ts[:0])

    def test_last_subset(self):
        ts = _simple_ts('1/1/2000', '1/1/2010', freq='12h')
        result = ts.last('10d')
        self.assert_(len(result) == 20)

        ts = _simple_ts('1/1/2000', '1/1/2010')
        result = ts.last('10d')
        self.assert_(len(result) == 10)

        result = ts.last('21D')
        expected = ts['12/12/2009':]
        assert_series_equal(result, expected)

        result = ts.last('21D')
        expected = ts[-21:]
        assert_series_equal(result, expected)

        result = ts[:0].last('3M')
        assert_series_equal(result, ts[:0])

    def test_add_offset(self):
        rng = date_range('1/1/2000', '2/1/2000')

        result = rng + offsets.Hour(2)
        expected = date_range('1/1/2000 02:00', '2/1/2000 02:00')
        self.assert_(result.equals(expected))

    def test_format_pre_1900_dates(self):
        rng = date_range('1/1/1850', '1/1/1950', freq='A-DEC')
        rng.format()
        ts = Series(1, index=rng)
        repr(ts)

    def test_repeat(self):
        rng = date_range('1/1/2000', '1/1/2001')

        result = rng.repeat(5)
        self.assert_(result.freq is None)
        self.assert_(len(result) == 5 * len(rng))

    def test_at_time(self):
        rng = date_range('1/1/2000', '1/5/2000', freq='5min')
        ts = Series(np.random.randn(len(rng)), index=rng)
        rs = ts.at_time(rng[1])
        self.assert_((rs.index.hour == rng[1].hour).all())
        self.assert_((rs.index.minute == rng[1].minute).all())
        self.assert_((rs.index.second == rng[1].second).all())

        result = ts.at_time('9:30')
        expected = ts.at_time(time(9, 30))
        assert_series_equal(result, expected)

        df = DataFrame(np.random.randn(len(rng), 3), index=rng)

        result = ts[time(9, 30)]
        result_df = df.ix[time(9, 30)]
        expected = ts[(rng.hour == 9) & (rng.minute == 30)]
        exp_df = df[(rng.hour == 9) & (rng.minute == 30)]

        # expected.index = date_range('1/1/2000', '1/4/2000')

        assert_series_equal(result, expected)
        tm.assert_frame_equal(result_df, exp_df)

        chunk = df.ix['1/4/2000':]
        result = chunk.ix[time(9, 30)]
        expected = result_df[-1:]
        tm.assert_frame_equal(result, expected)

        # midnight, everything
        rng = date_range('1/1/2000', '1/31/2000')
        ts = Series(np.random.randn(len(rng)), index=rng)

        result = ts.at_time(time(0, 0))
        assert_series_equal(result, ts)

        # time doesn't exist
        rng = date_range('1/1/2012', freq='23Min', periods=384)
        ts = Series(np.random.randn(len(rng)), rng)
        rs = ts.at_time('16:00')
        self.assert_(len(rs) == 0)

    def test_at_time_frame(self):
        rng = date_range('1/1/2000', '1/5/2000', freq='5min')
        ts = DataFrame(np.random.randn(len(rng), 2), index=rng)
        rs = ts.at_time(rng[1])
        self.assert_((rs.index.hour == rng[1].hour).all())
        self.assert_((rs.index.minute == rng[1].minute).all())
        self.assert_((rs.index.second == rng[1].second).all())

        result = ts.at_time('9:30')
        expected = ts.at_time(time(9, 30))
        assert_frame_equal(result, expected)

        result = ts.ix[time(9, 30)]
        expected = ts.ix[(rng.hour == 9) & (rng.minute == 30)]

        assert_frame_equal(result, expected)

        # midnight, everything
        rng = date_range('1/1/2000', '1/31/2000')
        ts = DataFrame(np.random.randn(len(rng), 3), index=rng)

        result = ts.at_time(time(0, 0))
        assert_frame_equal(result, ts)

        # time doesn't exist
        rng = date_range('1/1/2012', freq='23Min', periods=384)
        ts = DataFrame(np.random.randn(len(rng), 2), rng)
        rs = ts.at_time('16:00')
        self.assert_(len(rs) == 0)

    def test_between_time(self):
        rng = date_range('1/1/2000', '1/5/2000', freq='5min')
        ts = Series(np.random.randn(len(rng)), index=rng)
        stime = time(0, 0)
        etime = time(1, 0)

        close_open = product([True, False], [True, False])
        for inc_start, inc_end in close_open:
            filtered = ts.between_time(stime, etime, inc_start, inc_end)
            exp_len = 13 * 4 + 1
            if not inc_start:
                exp_len -= 5
            if not inc_end:
                exp_len -= 4

            self.assert_(len(filtered) == exp_len)
            for rs in filtered.index:
                t = rs.time()
                if inc_start:
                    self.assert_(t >= stime)
                else:
                    self.assert_(t > stime)

                if inc_end:
                    self.assert_(t <= etime)
                else:
                    self.assert_(t < etime)

        result = ts.between_time('00:00', '01:00')
        expected = ts.between_time(stime, etime)
        assert_series_equal(result, expected)

        # across midnight
        rng = date_range('1/1/2000', '1/5/2000', freq='5min')
        ts = Series(np.random.randn(len(rng)), index=rng)
        stime = time(22, 0)
        etime = time(9, 0)

        close_open = product([True, False], [True, False])
        for inc_start, inc_end in close_open:
            filtered = ts.between_time(stime, etime, inc_start, inc_end)
            exp_len = (12 * 11 + 1) * 4 + 1
            if not inc_start:
                exp_len -= 4
            if not inc_end:
                exp_len -= 4

            self.assert_(len(filtered) == exp_len)
            for rs in filtered.index:
                t = rs.time()
                if inc_start:
                    self.assert_((t >= stime) or (t <= etime))
                else:
                    self.assert_((t > stime) or (t <= etime))

                if inc_end:
                    self.assert_((t <= etime) or (t >= stime))
                else:
                    self.assert_((t < etime) or (t >= stime))

    def test_between_time_frame(self):
        rng = date_range('1/1/2000', '1/5/2000', freq='5min')
        ts = DataFrame(np.random.randn(len(rng), 2), index=rng)
        stime = time(0, 0)
        etime = time(1, 0)

        close_open = product([True, False], [True, False])
        for inc_start, inc_end in close_open:
            filtered = ts.between_time(stime, etime, inc_start, inc_end)
            exp_len = 13 * 4 + 1
            if not inc_start:
                exp_len -= 5
            if not inc_end:
                exp_len -= 4

            self.assert_(len(filtered) == exp_len)
            for rs in filtered.index:
                t = rs.time()
                if inc_start:
                    self.assert_(t >= stime)
                else:
                    self.assert_(t > stime)

                if inc_end:
                    self.assert_(t <= etime)
                else:
                    self.assert_(t < etime)

        result = ts.between_time('00:00', '01:00')
        expected = ts.between_time(stime, etime)
        assert_frame_equal(result, expected)

        # across midnight
        rng = date_range('1/1/2000', '1/5/2000', freq='5min')
        ts = DataFrame(np.random.randn(len(rng), 2), index=rng)
        stime = time(22, 0)
        etime = time(9, 0)

        close_open = product([True, False], [True, False])
        for inc_start, inc_end in close_open:
            filtered = ts.between_time(stime, etime, inc_start, inc_end)
            exp_len = (12 * 11 + 1) * 4 + 1
            if not inc_start:
                exp_len -= 4
            if not inc_end:
                exp_len -= 4

            self.assert_(len(filtered) == exp_len)
            for rs in filtered.index:
                t = rs.time()
                if inc_start:
                    self.assert_((t >= stime) or (t <= etime))
                else:
                    self.assert_((t > stime) or (t <= etime))

                if inc_end:
                    self.assert_((t <= etime) or (t >= stime))
                else:
                    self.assert_((t < etime) or (t >= stime))

    def test_dti_constructor_preserve_dti_freq(self):
        rng = date_range('1/1/2000', '1/2/2000', freq='5min')

        rng2 = DatetimeIndex(rng)
        self.assert_(rng.freq == rng2.freq)

    def test_normalize(self):
        rng = date_range('1/1/2000 9:30', periods=10, freq='D')

        result = rng.normalize()
        expected = date_range('1/1/2000', periods=10, freq='D')
        self.assert_(result.equals(expected))

        rng_ns = pd.DatetimeIndex(np.array([1380585623454345752, 1380585612343234312]).astype("datetime64[ns]"))
        rng_ns_normalized = rng_ns.normalize()
        expected = pd.DatetimeIndex(np.array([1380585600000000000, 1380585600000000000]).astype("datetime64[ns]"))
        self.assert_(rng_ns_normalized.equals(expected))

        self.assert_(result.is_normalized)
        self.assert_(not rng.is_normalized)

    def test_to_period(self):
        from pandas.tseries.period import period_range

        ts = _simple_ts('1/1/2000', '1/1/2001')

        pts = ts.to_period()
        exp = ts.copy()
        exp.index = period_range('1/1/2000', '1/1/2001')
        assert_series_equal(pts, exp)

        pts = ts.to_period('M')
        self.assert_(pts.index.equals(exp.index.asfreq('M')))

    def create_dt64_based_index(self):
        data = [Timestamp('2007-01-01 10:11:12.123456Z'),
                Timestamp('2007-01-01 10:11:13.789123Z')]
        index = DatetimeIndex(data)
        return index

    def test_to_period_millisecond(self):
        index = self.create_dt64_based_index()

        period = index.to_period(freq='L')
        self.assertEqual(2, len(period))
        self.assertEqual(period[0], Period('2007-01-01 10:11:12.123Z', 'L'))
        self.assertEqual(period[1], Period('2007-01-01 10:11:13.789Z', 'L'))

    def test_to_period_microsecond(self):
        index = self.create_dt64_based_index()

        period = index.to_period(freq='U')
        self.assertEqual(2, len(period))
        self.assertEqual(period[0], Period('2007-01-01 10:11:12.123456Z', 'U'))
        self.assertEqual(period[1], Period('2007-01-01 10:11:13.789123Z', 'U'))

    def test_to_period_tz(self):
        _skip_if_no_pytz()
        from dateutil.tz import tzlocal
        from pytz import utc as UTC

        xp = date_range('1/1/2000', '4/1/2000').to_period()

        ts = date_range('1/1/2000', '4/1/2000', tz='US/Eastern')

        result = ts.to_period()[0]
        expected = ts[0].to_period()

        self.assert_(result == expected)
        self.assert_(ts.to_period().equals(xp))

        ts = date_range('1/1/2000', '4/1/2000', tz=UTC)

        result = ts.to_period()[0]
        expected = ts[0].to_period()

        self.assert_(result == expected)
        self.assert_(ts.to_period().equals(xp))

        ts = date_range('1/1/2000', '4/1/2000', tz=tzlocal())

        result = ts.to_period()[0]
        expected = ts[0].to_period()

        self.assert_(result == expected)
        self.assert_(ts.to_period().equals(xp))

    def test_frame_to_period(self):
        K = 5
        from pandas.tseries.period import period_range

        dr = date_range('1/1/2000', '1/1/2001')
        pr = period_range('1/1/2000', '1/1/2001')
        df = DataFrame(randn(len(dr), K), index=dr)
        df['mix'] = 'a'

        pts = df.to_period()
        exp = df.copy()
        exp.index = pr
        assert_frame_equal(pts, exp)

        pts = df.to_period('M')
        self.assert_(pts.index.equals(exp.index.asfreq('M')))

        df = df.T
        pts = df.to_period(axis=1)
        exp = df.copy()
        exp.columns = pr
        assert_frame_equal(pts, exp)

        pts = df.to_period('M', axis=1)
        self.assert_(pts.columns.equals(exp.columns.asfreq('M')))

        self.assertRaises(ValueError, df.to_period, axis=2)

    def test_timestamp_fields(self):
        # extra fields from DatetimeIndex like quarter and week
        idx = tm.makeDateIndex(100)

        fields = ['dayofweek', 'dayofyear', 'week', 'weekofyear', 'quarter']
        for f in fields:
            expected = getattr(idx, f)[-1]
            result = getattr(Timestamp(idx[-1]), f)
            self.assertEqual(result, expected)

        self.assertEqual(idx.freq, Timestamp(idx[-1], idx.freq).freq)
        self.assertEqual(idx.freqstr, Timestamp(idx[-1], idx.freq).freqstr)

    def test_woy_boundary(self):
        # make sure weeks at year boundaries are correct
        d = datetime(2013,12,31)
        result = Timestamp(d).week
        expected = 1 # ISO standard
        self.assertEqual(result, expected)

        d = datetime(2008,12,28)
        result = Timestamp(d).week
        expected = 52 # ISO standard
        self.assertEqual(result, expected)

        d = datetime(2009,12,31)
        result = Timestamp(d).week
        expected = 53 # ISO standard
        self.assertEqual(result, expected)

        d = datetime(2010,1,1)
        result = Timestamp(d).week
        expected = 53 # ISO standard
        self.assertEqual(result, expected)

        d = datetime(2010,1,3)
        result = Timestamp(d).week
        expected = 53 # ISO standard
        self.assertEqual(result, expected)

        result = np.array([Timestamp(datetime(*args)).week for args in
                           [(2000,1,1),(2000,1,2),(2005,1,1),(2005,1,2)]])
        self.assertTrue((result == [52, 52, 53, 53]).all())

    def test_timestamp_date_out_of_range(self):
        self.assertRaises(ValueError, Timestamp, '1676-01-01')
        self.assertRaises(ValueError, Timestamp, '2263-01-01')

        # 1475
        self.assertRaises(ValueError, DatetimeIndex, ['1400-01-01'])
        self.assertRaises(ValueError, DatetimeIndex, [datetime(1400, 1, 1)])

    def test_timestamp_repr(self):
        # pre-1900
        stamp = Timestamp('1850-01-01', tz='US/Eastern')
        repr(stamp)

        iso8601 = '1850-01-01 01:23:45.012345'
        stamp = Timestamp(iso8601, tz='US/Eastern')
        result = repr(stamp)
        self.assert_(iso8601 in result)

    def test_timestamp_from_ordinal(self):

        # GH 3042
        dt = datetime(2011, 4, 16, 0, 0)
        ts = Timestamp.fromordinal(dt.toordinal())
        self.assert_(ts.to_pydatetime() == dt)

        # with a tzinfo
        stamp = Timestamp('2011-4-16', tz='US/Eastern')
        dt_tz = stamp.to_pydatetime()
        ts = Timestamp.fromordinal(dt_tz.toordinal(),tz='US/Eastern')
        self.assert_(ts.to_pydatetime() == dt_tz)

    def test_datetimeindex_integers_shift(self):
        rng = date_range('1/1/2000', periods=20)

        result = rng + 5
        expected = rng.shift(5)
        self.assert_(result.equals(expected))

        result = rng - 5
        expected = rng.shift(-5)
        self.assert_(result.equals(expected))

    def test_astype_object(self):
        # NumPy 1.6.1 weak ns support
        rng = date_range('1/1/2000', periods=20)

        casted = rng.astype('O')
        exp_values = list(rng)

        self.assert_(np.array_equal(casted, exp_values))

    def test_catch_infinite_loop(self):
        offset = datetools.DateOffset(minute=5)
        # blow up, don't loop forever
        self.assertRaises(Exception, date_range, datetime(2011, 11, 11),
                          datetime(2011, 11, 12), freq=offset)

    def test_append_concat(self):
        rng = date_range('5/8/2012 1:45', periods=10, freq='5T')
        ts = Series(np.random.randn(len(rng)), rng)
        df = DataFrame(np.random.randn(len(rng), 4), index=rng)

        result = ts.append(ts)
        result_df = df.append(df)
        ex_index = DatetimeIndex(np.tile(rng.values, 2))
        self.assert_(result.index.equals(ex_index))
        self.assert_(result_df.index.equals(ex_index))

        appended = rng.append(rng)
        self.assert_(appended.equals(ex_index))

        appended = rng.append([rng, rng])
        ex_index = DatetimeIndex(np.tile(rng.values, 3))
        self.assert_(appended.equals(ex_index))

        # different index names
        rng1 = rng.copy()
        rng2 = rng.copy()
        rng1.name = 'foo'
        rng2.name = 'bar'
        self.assert_(rng1.append(rng1).name == 'foo')
        self.assert_(rng1.append(rng2).name is None)

    def test_append_concat_tz(self):
        #GH 2938
        _skip_if_no_pytz()

        rng = date_range('5/8/2012 1:45', periods=10, freq='5T',
                         tz='US/Eastern')
        rng2 = date_range('5/8/2012 2:35', periods=10, freq='5T',
                         tz='US/Eastern')
        rng3 = date_range('5/8/2012 1:45', periods=20, freq='5T',
                         tz='US/Eastern')
        ts = Series(np.random.randn(len(rng)), rng)
        df = DataFrame(np.random.randn(len(rng), 4), index=rng)
        ts2 = Series(np.random.randn(len(rng2)), rng2)
        df2 = DataFrame(np.random.randn(len(rng2), 4), index=rng2)

        result = ts.append(ts2)
        result_df = df.append(df2)
        self.assert_(result.index.equals(rng3))
        self.assert_(result_df.index.equals(rng3))

        appended = rng.append(rng2)
        self.assert_(appended.equals(rng3))

    def test_set_dataframe_column_ns_dtype(self):
        x = DataFrame([datetime.now(), datetime.now()])
        self.assert_(x[0].dtype == np.dtype('M8[ns]'))

    def test_groupby_count_dateparseerror(self):
        dr = date_range(start='1/1/2012', freq='5min', periods=10)

        # BAD Example, datetimes first
        s = Series(np.arange(10), index=[dr, lrange(10)])
        grouped = s.groupby(lambda x: x[1] % 2 == 0)
        result = grouped.count()

        s = Series(np.arange(10), index=[lrange(10), dr])
        grouped = s.groupby(lambda x: x[0] % 2 == 0)
        expected = grouped.count()

        assert_series_equal(result, expected)

    def test_datetimeindex_repr_short(self):
        dr = date_range(start='1/1/2012', periods=1)
        repr(dr)

        dr = date_range(start='1/1/2012', periods=2)
        repr(dr)

        dr = date_range(start='1/1/2012', periods=3)
        repr(dr)

    def test_constructor_int64_nocopy(self):
        # #1624
        arr = np.arange(1000, dtype=np.int64)
        index = DatetimeIndex(arr)

        arr[50:100] = -1
        self.assert_((index.asi8[50:100] == -1).all())

        arr = np.arange(1000, dtype=np.int64)
        index = DatetimeIndex(arr, copy=True)

        arr[50:100] = -1
        self.assert_((index.asi8[50:100] != -1).all())

    def test_series_interpolate_method_values(self):
        # #1646
        ts = _simple_ts('1/1/2000', '1/20/2000')
        ts[::2] = np.nan

        result = ts.interpolate(method='values')
        exp = ts.interpolate()
        assert_series_equal(result, exp)

    def test_frame_datetime64_handling_groupby(self):
        # it works!
        df = DataFrame([(3, np.datetime64('2012-07-03')),
                        (3, np.datetime64('2012-07-04'))],
                       columns=['a', 'date'])
        result = df.groupby('a').first()
        self.assertEqual(result['date'][3], Timestamp('2012-07-03'))

    def test_series_interpolate_intraday(self):
        # #1698
        index = pd.date_range('1/1/2012', periods=4, freq='12D')
        ts = pd.Series([0, 12, 24, 36], index)
        new_index = index.append(index + pd.DateOffset(days=1)).order()

        exp = ts.reindex(new_index).interpolate(method='time')

        index = pd.date_range('1/1/2012', periods=4, freq='12H')
        ts = pd.Series([0, 12, 24, 36], index)
        new_index = index.append(index + pd.DateOffset(hours=1)).order()
        result = ts.reindex(new_index).interpolate(method='time')

        self.assert_(np.array_equal(result.values, exp.values))

    def test_frame_dict_constructor_datetime64_1680(self):
        dr = date_range('1/1/2012', periods=10)
        s = Series(dr, index=dr)

        # it works!
        DataFrame({'a': 'foo', 'b': s}, index=dr)
        DataFrame({'a': 'foo', 'b': s.values}, index=dr)

    def test_frame_datetime64_mixed_index_ctor_1681(self):
        dr = date_range('2011/1/1', '2012/1/1', freq='W-FRI')
        ts = Series(dr)

        # it works!
        d = DataFrame({'A': 'foo', 'B': ts}, index=dr)
        self.assert_(d['B'].isnull().all())

    def test_frame_timeseries_to_records(self):
        index = date_range('1/1/2000', periods=10)
        df = DataFrame(np.random.randn(10, 3), index=index,
                       columns=['a', 'b', 'c'])

        result = df.to_records()
        result['index'].dtype == 'M8[ns]'

        result = df.to_records(index=False)

    def test_frame_datetime64_duplicated(self):
        dates = date_range('2010-07-01', end='2010-08-05')

        tst = DataFrame({'symbol': 'AAA', 'date': dates})
        result = tst.duplicated(['date', 'symbol'])
        self.assert_((-result).all())

        tst = DataFrame({'date': dates})
        result = tst.duplicated()
        self.assert_((-result).all())

    def test_timestamp_compare_with_early_datetime(self):
        # e.g. datetime.min
        stamp = Timestamp('2012-01-01')

        self.assertFalse(stamp == datetime.min)
        self.assertFalse(stamp == datetime(1600, 1, 1))
        self.assertFalse(stamp == datetime(2700, 1, 1))
        self.assert_(stamp != datetime.min)
        self.assert_(stamp != datetime(1600, 1, 1))
        self.assert_(stamp != datetime(2700, 1, 1))
        self.assert_(stamp > datetime(1600, 1, 1))
        self.assert_(stamp >= datetime(1600, 1, 1))
        self.assert_(stamp < datetime(2700, 1, 1))
        self.assert_(stamp <= datetime(2700, 1, 1))

    def test_to_html_timestamp(self):
        rng = date_range('2000-01-01', periods=10)
        df = DataFrame(np.random.randn(10, 4), index=rng)

        result = df.to_html()
        self.assert_('2000-01-01' in result)

    def test_to_csv_numpy_16_bug(self):
        frame = DataFrame({'a': date_range('1/1/2000', periods=10)})

        buf = StringIO()
        frame.to_csv(buf)

        result = buf.getvalue()
        self.assert_('2000-01-01' in result)

    def test_series_map_box_timestamps(self):
        # #2689, #2627
        s = Series(date_range('1/1/2000', periods=10))

        def f(x):
            return (x.hour, x.day, x.month)

        # it works!
        s.map(f)
        s.apply(f)
        DataFrame(s).applymap(f)

    def test_concat_datetime_datetime64_frame(self):
        # #2624
        rows = []
        rows.append([datetime(2010, 1, 1), 1])
        rows.append([datetime(2010, 1, 2), 'hi'])

        df2_obj = DataFrame.from_records(rows, columns=['date', 'test'])

        ind = date_range(start="2000/1/1", freq="D", periods=10)
        df1 = DataFrame({'date': ind, 'test':lrange(10)})

        # it works!
        pd.concat([df1, df2_obj])

    def test_period_resample(self):
        # GH3609
        s = Series(range(100),index=date_range('20130101', freq='s', periods=100), dtype='float')
        s[10:30] = np.nan
        expected = Series([34.5, 79.5], index=[Period('2013-01-01 00:00', 'T'), Period('2013-01-01 00:01', 'T')])
        result = s.to_period().resample('T', kind='period')
        assert_series_equal(result, expected)
        result2 = s.resample('T', kind='period')
        assert_series_equal(result2, expected)

    def test_period_resample_with_local_timezone(self):
        # GH5430
        _skip_if_no_pytz()
        import pytz

        local_timezone = pytz.timezone('America/Los_Angeles')

        start = datetime(year=2013, month=11, day=1, hour=0, minute=0, tzinfo=pytz.utc)
        # 1 day later
        end = datetime(year=2013, month=11, day=2, hour=0, minute=0, tzinfo=pytz.utc)

        index = pd.date_range(start, end, freq='H')

        series = pd.Series(1, index=index)
        series = series.tz_convert(local_timezone)
        result = series.resample('D', kind='period')
        # Create the expected series
        expected_index = (pd.period_range(start=start, end=end, freq='D') - 1)  # Index is moved back a day with the timezone conversion from UTC to Pacific
        expected = pd.Series(1, index=expected_index)
        assert_series_equal(result, expected)


def _simple_ts(start, end, freq='D'):
    rng = date_range(start, end, freq=freq)
    return Series(np.random.randn(len(rng)), index=rng)


class TestDatetimeIndex(tm.TestCase):
    _multiprocess_can_split_ = True

    def test_hash_error(self):
        index = date_range('20010101', periods=10)
        with tm.assertRaisesRegexp(TypeError,
                                   "unhashable type: %r" %
                                   type(index).__name__):
            hash(index)

    def test_stringified_slice_with_tz(self):
        #GH2658
        import datetime
        start=datetime.datetime.now()
        idx=DatetimeIndex(start=start,freq="1d",periods=10)
        df=DataFrame(lrange(10),index=idx)
        df["2013-01-14 23:44:34.437768-05:00":] # no exception here

    def test_append_join_nondatetimeindex(self):
        rng = date_range('1/1/2000', periods=10)
        idx = Index(['a', 'b', 'c', 'd'])

        result = rng.append(idx)
        tm.assert_isinstance(result[0], Timestamp)

        # it works
        rng.join(idx, how='outer')

    def test_astype(self):
        rng = date_range('1/1/2000', periods=10)

        result = rng.astype('i8')
        self.assert_(np.array_equal(result, rng.asi8))

    def test_to_period_nofreq(self):
        idx = DatetimeIndex(['2000-01-01', '2000-01-02', '2000-01-04'])
        self.assertRaises(ValueError, idx.to_period)

        idx = DatetimeIndex(['2000-01-01', '2000-01-02', '2000-01-03'],
                            freq='infer')
        idx.to_period()

    def test_000constructor_resolution(self):
        # 2252
        t1 = Timestamp((1352934390 * 1000000000) + 1000000 + 1000 + 1)
        idx = DatetimeIndex([t1])

        self.assert_(idx.nanosecond[0] == t1.nanosecond)

    def test_constructor_coverage(self):
        rng = date_range('1/1/2000', periods=10.5)
        exp = date_range('1/1/2000', periods=10)
        self.assert_(rng.equals(exp))

        self.assertRaises(ValueError, DatetimeIndex, start='1/1/2000',
                          periods='foo', freq='D')

        self.assertRaises(ValueError, DatetimeIndex, start='1/1/2000',
                          end='1/10/2000')

        self.assertRaises(ValueError, DatetimeIndex, '1/1/2000')

        # generator expression
        gen = (datetime(2000, 1, 1) + timedelta(i) for i in range(10))
        result = DatetimeIndex(gen)
        expected = DatetimeIndex([datetime(2000, 1, 1) + timedelta(i)
                                  for i in range(10)])
        self.assert_(result.equals(expected))

        # NumPy string array
        strings = np.array(['2000-01-01', '2000-01-02', '2000-01-03'])
        result = DatetimeIndex(strings)
        expected = DatetimeIndex(strings.astype('O'))
        self.assert_(result.equals(expected))

        from_ints = DatetimeIndex(expected.asi8)
        self.assert_(from_ints.equals(expected))

        # non-conforming
        self.assertRaises(ValueError, DatetimeIndex,
                          ['2000-01-01', '2000-01-02', '2000-01-04'],
                          freq='D')

        self.assertRaises(ValueError, DatetimeIndex,
                          start='2011-01-01', freq='b')
        self.assertRaises(ValueError, DatetimeIndex,
                          end='2011-01-01', freq='B')
        self.assertRaises(ValueError, DatetimeIndex, periods=10, freq='D')

    def test_constructor_name(self):
        idx = DatetimeIndex(start='2000-01-01', periods=1, freq='A',
                            name='TEST')
        self.assertEquals(idx.name, 'TEST')

    def test_comparisons_coverage(self):
        rng = date_range('1/1/2000', periods=10)

        # raise TypeError for now
        self.assertRaises(TypeError, rng.__lt__, rng[3].value)

        result = rng == list(rng)
        exp = rng == rng
        self.assert_(np.array_equal(result, exp))

    def test_map(self):
        rng = date_range('1/1/2000', periods=10)

        f = lambda x: x.strftime('%Y%m%d')
        result = rng.map(f)
        exp = [f(x) for x in rng]
        self.assert_(np.array_equal(result, exp))

    def test_add_union(self):
        rng = date_range('1/1/2000', periods=5)
        rng2 = date_range('1/6/2000', periods=5)

        result = rng + rng2
        expected = rng.union(rng2)
        self.assert_(result.equals(expected))

    def test_misc_coverage(self):
        rng = date_range('1/1/2000', periods=5)
        result = rng.groupby(rng.day)
        tm.assert_isinstance(list(result.values())[0][0], Timestamp)

        idx = DatetimeIndex(['2000-01-03', '2000-01-01', '2000-01-02'])
        self.assert_(idx.equals(list(idx)))

        non_datetime = Index(list('abc'))
        self.assert_(not idx.equals(list(non_datetime)))

    def test_union_coverage(self):
        idx = DatetimeIndex(['2000-01-03', '2000-01-01', '2000-01-02'])
        ordered = DatetimeIndex(idx.order(), freq='infer')
        result = ordered.union(idx)
        self.assert_(result.equals(ordered))

        result = ordered[:0].union(ordered)
        self.assert_(result.equals(ordered))
        self.assert_(result.freq == ordered.freq)

    def test_union_bug_1730(self):
        rng_a = date_range('1/1/2012', periods=4, freq='3H')
        rng_b = date_range('1/1/2012', periods=4, freq='4H')

        result = rng_a.union(rng_b)
        exp = DatetimeIndex(sorted(set(list(rng_a)) | set(list(rng_b))))
        self.assert_(result.equals(exp))

    def test_union_bug_1745(self):
        left = DatetimeIndex(['2012-05-11 15:19:49.695000'])
        right = DatetimeIndex(['2012-05-29 13:04:21.322000',
                               '2012-05-11 15:27:24.873000',
                               '2012-05-11 15:31:05.350000'])

        result = left.union(right)
        exp = DatetimeIndex(sorted(set(list(left)) | set(list(right))))
        self.assert_(result.equals(exp))

    def test_union_bug_4564(self):
        from pandas import DateOffset
        left = date_range("2013-01-01", "2013-02-01")
        right = left + DateOffset(minutes=15)

        result = left.union(right)
        exp = DatetimeIndex(sorted(set(list(left)) | set(list(right))))
        self.assert_(result.equals(exp))

    def test_intersection_bug_1708(self):
        from pandas import DateOffset
        index_1 = date_range('1/1/2012', periods=4, freq='12H')
        index_2 = index_1 + DateOffset(hours=1)

        result = index_1 & index_2
        self.assertEqual(len(result), 0)

    # def test_add_timedelta64(self):
    #     rng = date_range('1/1/2000', periods=5)
    #     delta = rng.values[3] - rng.values[1]

    #     result = rng + delta
    #     expected = rng + timedelta(2)
    #     self.assert_(result.equals(expected))

    def test_get_duplicates(self):
        idx = DatetimeIndex(['2000-01-01', '2000-01-02', '2000-01-02',
                             '2000-01-03', '2000-01-03', '2000-01-04'])

        result = idx.get_duplicates()
        ex = DatetimeIndex(['2000-01-02', '2000-01-03'])
        self.assert_(result.equals(ex))

    def test_argmin_argmax(self):
        idx = DatetimeIndex(['2000-01-04', '2000-01-01', '2000-01-02'])
        self.assertEqual(idx.argmin(), 1)
        self.assertEqual(idx.argmax(), 0)

    def test_order(self):
        idx = DatetimeIndex(['2000-01-04', '2000-01-01', '2000-01-02'])

        ordered = idx.order()
        self.assert_(ordered.is_monotonic)

        ordered = idx.order(ascending=False)
        self.assert_(ordered[::-1].is_monotonic)

        ordered, dexer = idx.order(return_indexer=True)
        self.assert_(ordered.is_monotonic)
        self.assert_(np.array_equal(dexer, [1, 2, 0]))

        ordered, dexer = idx.order(return_indexer=True, ascending=False)
        self.assert_(ordered[::-1].is_monotonic)
        self.assert_(np.array_equal(dexer, [0, 2, 1]))

    def test_insert(self):
        idx = DatetimeIndex(['2000-01-04', '2000-01-01', '2000-01-02'])

        result = idx.insert(2, datetime(2000, 1, 5))
        exp = DatetimeIndex(['2000-01-04', '2000-01-01', '2000-01-05',
                             '2000-01-02'])
        self.assert_(result.equals(exp))

        # insertion of non-datetime should coerce to object index
        result = idx.insert(1, 'inserted')
        expected = Index([datetime(2000, 1, 4), 'inserted', datetime(2000, 1, 1),
                          datetime(2000, 1, 2)])
        self.assert_(not isinstance(result, DatetimeIndex))
        tm.assert_index_equal(result, expected)

        idx = date_range('1/1/2000', periods=3, freq='M')
        result = idx.insert(3, datetime(2000, 4, 30))
        self.assert_(result.freqstr == 'M')

    def test_map_bug_1677(self):
        index = DatetimeIndex(['2012-04-25 09:30:00.393000'])
        f = index.asof

        result = index.map(f)
        expected = np.array([f(index[0])])
        self.assert_(np.array_equal(result, expected))

    def test_groupby_function_tuple_1677(self):
        df = DataFrame(np.random.rand(100),
                       index=date_range("1/1/2000", periods=100))
        monthly_group = df.groupby(lambda x: (x.year, x.month))

        result = monthly_group.mean()
        tm.assert_isinstance(result.index[0], tuple)

    def test_append_numpy_bug_1681(self):
        # another datetime64 bug
        dr = date_range('2011/1/1', '2012/1/1', freq='W-FRI')
        a = DataFrame()
        c = DataFrame({'A': 'foo', 'B': dr}, index=dr)

        result = a.append(c)
        self.assert_((result['B'] == dr).all())

    def test_isin(self):
        index = tm.makeDateIndex(4)
        result = index.isin(index)
        self.assert_(result.all())

        result = index.isin(list(index))
        self.assert_(result.all())

        assert_almost_equal(index.isin([index[2], 5]),
                            [False, False, True, False])

    def test_union(self):
        i1 = Int64Index(np.arange(0, 20, 2))
        i2 = Int64Index(np.arange(10, 30, 2))
        result = i1.union(i2)
        expected = Int64Index(np.arange(0, 30, 2))
        self.assert_(np.array_equal(result, expected))

    def test_union_with_DatetimeIndex(self):
        i1 = Int64Index(np.arange(0, 20, 2))
        i2 = DatetimeIndex(start='2012-01-03 00:00:00', periods=10, freq='D')
        i1.union(i2)  # Works
        i2.union(i1)  # Fails with "AttributeError: can't set attribute"

    def test_time(self):
        rng = pd.date_range('1/1/2000', freq='12min', periods=10)
        result = pd.Index(rng).time
        expected = [t.time() for t in rng]
        self.assert_((result == expected).all())

    def test_date(self):
        rng = pd.date_range('1/1/2000', freq='12H', periods=10)
        result = pd.Index(rng).date
        expected = [t.date() for t in rng]
        self.assert_((result == expected).all())

    def test_does_not_convert_mixed_integer(self):
        df = tm.makeCustomDataframe(10, 10, data_gen_f=lambda *args, **kwargs:
                                    randn(), r_idx_type='i', c_idx_type='dt')
        cols = df.columns.join(df.index, how='outer')
        joined = cols.join(df.columns)
        self.assertEqual(cols.dtype, np.dtype('O'))
        self.assertEqual(cols.dtype, joined.dtype)
        assert_array_equal(cols.values, joined.values)

    def test_slice_keeps_name(self):
        # GH4226
        st = pd.Timestamp('2013-07-01 00:00:00', tz='America/Los_Angeles')
        et = pd.Timestamp('2013-07-02 00:00:00', tz='America/Los_Angeles')
        dr = pd.date_range(st, et, freq='H', name='timebucket')
        self.assertEqual(dr[1:].name, dr.name)

    def test_join_self(self):
        index = date_range('1/1/2000', periods=10)
        kinds = 'outer', 'inner', 'left', 'right'
        for kind in kinds:
            joined = index.join(index, how=kind)
            self.assert_(index is joined)

    def assert_index_parameters(self, index):
        assert index.freq == '40960N'
        assert index.inferred_freq == '40960N'

    def test_ns_index(self):

        if _np_version_under1p7:
            raise nose.SkipTest

        nsamples = 400
        ns = int(1e9 / 24414)
        dtstart = np.datetime64('2012-09-20T00:00:00')

        dt = dtstart + np.arange(nsamples) * np.timedelta64(ns, 'ns')
        freq = ns * pd.datetools.Nano()
        index = pd.DatetimeIndex(dt, freq=freq, name='time')
        self.assert_index_parameters(index)

        new_index = pd.DatetimeIndex(start=index[0], end=index[-1], freq=index.freq)
        self.assert_index_parameters(new_index)

    def test_join_with_period_index(self):
        df = tm.makeCustomDataframe(10, 10, data_gen_f=lambda *args:
                                    np.random.randint(2), c_idx_type='p',
                                    r_idx_type='dt')
        s = df.iloc[:5, 0]
        joins = 'left', 'right', 'inner', 'outer'

        for join in joins:
            with tm.assertRaisesRegexp(ValueError, 'can only call with other '
                                       'PeriodIndex-ed objects'):
                df.columns.join(s.index, how=join)


class TestDatetime64(tm.TestCase):
    """
    Also test supoprt for datetime64[ns] in Series / DataFrame
    """

    def setUp(self):
        dti = DatetimeIndex(start=datetime(2005, 1, 1),
                            end=datetime(2005, 1, 10), freq='Min')
        self.series = Series(rand(len(dti)), dti)

    def test_datetimeindex_accessors(self):
        dti = DatetimeIndex(
            freq='Q-JAN', start=datetime(1997, 12, 31), periods=100)

        self.assertEquals(dti.year[0], 1998)
        self.assertEquals(dti.month[0], 1)
        self.assertEquals(dti.day[0], 31)
        self.assertEquals(dti.hour[0], 0)
        self.assertEquals(dti.minute[0], 0)
        self.assertEquals(dti.second[0], 0)
        self.assertEquals(dti.microsecond[0], 0)
        self.assertEquals(dti.dayofweek[0], 5)

        self.assertEquals(dti.dayofyear[0], 31)
        self.assertEquals(dti.dayofyear[1], 120)

        self.assertEquals(dti.weekofyear[0], 5)
        self.assertEquals(dti.weekofyear[1], 18)

        self.assertEquals(dti.quarter[0], 1)
        self.assertEquals(dti.quarter[1], 2)

        self.assertEquals(len(dti.year), 100)
        self.assertEquals(len(dti.month), 100)
        self.assertEquals(len(dti.day), 100)
        self.assertEquals(len(dti.hour), 100)
        self.assertEquals(len(dti.minute), 100)
        self.assertEquals(len(dti.second), 100)
        self.assertEquals(len(dti.microsecond), 100)
        self.assertEquals(len(dti.dayofweek), 100)
        self.assertEquals(len(dti.dayofyear), 100)
        self.assertEquals(len(dti.weekofyear), 100)
        self.assertEquals(len(dti.quarter), 100)

    def test_nanosecond_field(self):
        dti = DatetimeIndex(np.arange(10))

        self.assert_(np.array_equal(dti.nanosecond, np.arange(10)))

    def test_datetimeindex_diff(self):
        dti1 = DatetimeIndex(freq='Q-JAN', start=datetime(1997, 12, 31),
                             periods=100)
        dti2 = DatetimeIndex(freq='Q-JAN', start=datetime(1997, 12, 31),
                             periods=98)
        self.assert_(len(dti1.diff(dti2)) == 2)

    def test_fancy_getitem(self):
        dti = DatetimeIndex(freq='WOM-1FRI', start=datetime(2005, 1, 1),
                            end=datetime(2010, 1, 1))

        s = Series(np.arange(len(dti)), index=dti)

        self.assertEquals(s[48], 48)
        self.assertEquals(s['1/2/2009'], 48)
        self.assertEquals(s['2009-1-2'], 48)
        self.assertEquals(s[datetime(2009, 1, 2)], 48)
        self.assertEquals(s[lib.Timestamp(datetime(2009, 1, 2))], 48)
        self.assertRaises(KeyError, s.__getitem__, '2009-1-3')

        assert_series_equal(s['3/6/2009':'2009-06-05'],
                            s[datetime(2009, 3, 6):datetime(2009, 6, 5)])

    def test_fancy_setitem(self):
        dti = DatetimeIndex(freq='WOM-1FRI', start=datetime(2005, 1, 1),
                            end=datetime(2010, 1, 1))

        s = Series(np.arange(len(dti)), index=dti)
        s[48] = -1
        self.assertEquals(s[48], -1)
        s['1/2/2009'] = -2
        self.assertEquals(s[48], -2)
        s['1/2/2009':'2009-06-05'] = -3
        self.assert_((s[48:54] == -3).all())

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

        arr = to_datetime(
            ['1/1/2005', '1/2/2005', 'Jan 3, 2005', '2005-01-04'])
        idx6 = DatetimeIndex(arr)

        idx7 = DatetimeIndex(['12/05/2007', '25/01/2008'], dayfirst=True)
        idx8 = DatetimeIndex(['2007/05/12', '2008/01/25'], dayfirst=False,
                             yearfirst=True)
        self.assert_(idx7.equals(idx8))

        for other in [idx2, idx3, idx4, idx5, idx6]:
            self.assert_((idx1.values == other.values).all())

        sdate = datetime(1999, 12, 25)
        edate = datetime(2000, 1, 1)
        idx = DatetimeIndex(start=sdate, freq='1B', periods=20)
        self.assertEquals(len(idx), 20)
        self.assertEquals(idx[0], sdate + 0 * dt.bday)
        self.assertEquals(idx.freq, 'B')

        idx = DatetimeIndex(end=edate, freq=('D', 5), periods=20)
        self.assertEquals(len(idx), 20)
        self.assertEquals(idx[-1], edate)
        self.assertEquals(idx.freq, '5D')

        idx1 = DatetimeIndex(start=sdate, end=edate, freq='W-SUN')
        idx2 = DatetimeIndex(start=sdate, end=edate,
                             freq=dt.Week(weekday=6))
        self.assertEquals(len(idx1), len(idx2))
        self.assertEquals(idx1.offset, idx2.offset)

        idx1 = DatetimeIndex(start=sdate, end=edate, freq='QS')
        idx2 = DatetimeIndex(start=sdate, end=edate,
                             freq=dt.QuarterBegin(startingMonth=1))
        self.assertEquals(len(idx1), len(idx2))
        self.assertEquals(idx1.offset, idx2.offset)

        idx1 = DatetimeIndex(start=sdate, end=edate, freq='BQ')
        idx2 = DatetimeIndex(start=sdate, end=edate,
                             freq=dt.BQuarterEnd(startingMonth=12))
        self.assertEquals(len(idx1), len(idx2))
        self.assertEquals(idx1.offset, idx2.offset)

    def test_dti_snap(self):
        dti = DatetimeIndex(['1/1/2002', '1/2/2002', '1/3/2002', '1/4/2002',
                             '1/5/2002', '1/6/2002', '1/7/2002'], freq='D')

        res = dti.snap(freq='W-MON')
        exp = date_range('12/31/2001', '1/7/2002', freq='w-mon')
        exp = exp.repeat([3, 4])
        self.assert_((res == exp).all())

        res = dti.snap(freq='B')

        exp = date_range('1/1/2002', '1/7/2002', freq='b')
        exp = exp.repeat([1, 1, 1, 2, 2])
        self.assert_((res == exp).all())

    def test_dti_reset_index_round_trip(self):
        dti = DatetimeIndex(start='1/1/2001', end='6/1/2001', freq='D')
        d1 = DataFrame({'v': np.random.rand(len(dti))}, index=dti)
        d2 = d1.reset_index()
        self.assert_(d2.dtypes[0] == np.dtype('M8[ns]'))
        d3 = d2.set_index('index')
        assert_frame_equal(d1, d3, check_names=False)

        # #2329
        stamp = datetime(2012, 11, 22)
        df = DataFrame([[stamp, 12.1]], columns=['Date', 'Value'])
        df = df.set_index('Date')

        self.assertEquals(df.index[0], stamp)
        self.assertEquals(df.reset_index()['Date'][0], stamp)

    def test_datetimeindex_union_join_empty(self):
        dti = DatetimeIndex(start='1/1/2001', end='2/1/2001', freq='D')
        empty = Index([])

        result = dti.union(empty)
        tm.assert_isinstance(result, DatetimeIndex)
        self.assert_(result is result)

        result = dti.join(empty)
        tm.assert_isinstance(result, DatetimeIndex)

    def test_series_set_value(self):
        # #1561

        dates = [datetime(2001, 1, 1), datetime(2001, 1, 2)]
        index = DatetimeIndex(dates)

        s = Series().set_value(dates[0], 1.)
        s2 = s.set_value(dates[1], np.nan)

        exp = Series([1., np.nan], index=index)

        assert_series_equal(s2, exp)

        # s = Series(index[:1], index[:1])
        # s2 = s.set_value(dates[1], index[1])
        # self.assert_(s2.values.dtype == 'M8[ns]')

    @slow
    def test_slice_locs_indexerror(self):
        times = [datetime(2000, 1, 1) + timedelta(minutes=i * 10)
                 for i in range(100000)]
        s = Series(lrange(100000), times)
        s.ix[datetime(1900, 1, 1):datetime(2100, 1, 1)]


class TestSeriesDatetime64(tm.TestCase):

    def setUp(self):
        self.series = Series(date_range('1/1/2000', periods=10))

    def test_auto_conversion(self):
        series = Series(list(date_range('1/1/2000', periods=10)))
        self.assert_(series.dtype == 'M8[ns]')

    def test_constructor_cant_cast_datetime64(self):
        self.assertRaises(TypeError, Series,
                          date_range('1/1/2000', periods=10), dtype=float)

    def test_series_comparison_scalars(self):
        val = datetime(2000, 1, 4)
        result = self.series > val
        expected = np.array([x > val for x in self.series])
        self.assert_(np.array_equal(result, expected))

        val = self.series[5]
        result = self.series > val
        expected = np.array([x > val for x in self.series])
        self.assert_(np.array_equal(result, expected))

    def test_between(self):
        left, right = self.series[[2, 7]]

        result = self.series.between(left, right)
        expected = (self.series >= left) & (self.series <= right)
        assert_series_equal(result, expected)

    #----------------------------------------------------------------------
    # NaT support

    def test_NaT_scalar(self):
        series = Series([0, 1000, 2000, iNaT], dtype='M8[ns]')

        val = series[3]
        self.assert_(com.isnull(val))

        series[2] = val
        self.assert_(com.isnull(series[2]))

    def test_set_none_nan(self):
        self.series[3] = None
        self.assert_(self.series[3] is NaT)

        self.series[3:5] = None
        self.assert_(self.series[4] is NaT)

        self.series[5] = np.nan
        self.assert_(self.series[5] is NaT)

        self.series[5:7] = np.nan
        self.assert_(self.series[6] is NaT)

    def test_intercept_astype_object(self):

        # this test no longer makes sense as series is by default already M8[ns]
        expected = self.series.astype('object')

        df = DataFrame({'a': self.series,
                        'b': np.random.randn(len(self.series))})

        result = df.values.squeeze()
        self.assert_((result[:, 0] == expected.values).all())

        df = DataFrame({'a': self.series,
                        'b': ['foo'] * len(self.series)})

        result = df.values.squeeze()
        self.assert_((result[:, 0] == expected.values).all())

    def test_union(self):
        rng1 = date_range('1/1/1999', '1/1/2012', freq='MS')
        s1 = Series(np.random.randn(len(rng1)), rng1)

        rng2 = date_range('1/1/1980', '12/1/2001', freq='MS')
        s2 = Series(np.random.randn(len(rng2)), rng2)
        df = DataFrame({'s1': s1, 's2': s2})
        self.assert_(df.index.values.dtype == np.dtype('M8[ns]'))

    def test_intersection(self):
        rng = date_range('6/1/2000', '6/15/2000', freq='D')
        rng = rng.delete(5)

        rng2 = date_range('5/15/2000', '6/20/2000', freq='D')
        rng2 = DatetimeIndex(rng2.values)

        result = rng.intersection(rng2)
        self.assert_(result.equals(rng))

        # empty same freq GH2129
        rng = date_range('6/1/2000', '6/15/2000', freq='T')
        result = rng[0:0].intersection(rng)
        self.assert_(len(result) == 0)

        result = rng.intersection(rng[0:0])
        self.assert_(len(result) == 0)

    def test_date_range_bms_bug(self):
        # #1645
        rng = date_range('1/1/2000', periods=10, freq='BMS')

        ex_first = Timestamp('2000-01-03')
        self.assertEquals(rng[0], ex_first)

    def test_string_index_series_name_converted(self):
        # #1644
        df = DataFrame(np.random.randn(10, 4),
                       index=date_range('1/1/2000', periods=10))

        result = df.ix['1/3/2000']
        self.assertEquals(result.name, df.index[2])

        result = df.T['1/3/2000']
        self.assertEquals(result.name, df.index[2])


class TestTimestamp(tm.TestCase):

    def test_class_ops(self):
        _skip_if_no_pytz()
        import pytz

        def compare(x,y):
            self.assertEqual(int(Timestamp(x).value/1e9), int(Timestamp(y).value/1e9))

        compare(Timestamp.now(),datetime.now())
        compare(Timestamp.now('UTC'),datetime.now(pytz.timezone('UTC')))
        compare(Timestamp.utcnow(),datetime.utcnow())
        compare(Timestamp.today(),datetime.today())

    def test_basics_nanos(self):
        val = np.int64(946684800000000000).view('M8[ns]')
        stamp = Timestamp(val.view('i8') + 500)
        self.assert_(stamp.year == 2000)
        self.assert_(stamp.month == 1)
        self.assert_(stamp.microsecond == 0)
        self.assert_(stamp.nanosecond == 500)

    def test_unit(self):
        def check(val,unit=None,h=1,s=1,us=0):
            stamp = Timestamp(val, unit=unit)
            self.assert_(stamp.year == 2000)
            self.assert_(stamp.month == 1)
            self.assert_(stamp.day == 1)
            self.assert_(stamp.hour == h)
            if unit != 'D':
                self.assert_(stamp.minute == 1)
                self.assert_(stamp.second == s)
                self.assert_(stamp.microsecond == us)
            else:
                self.assert_(stamp.minute == 0)
                self.assert_(stamp.second == 0)
                self.assert_(stamp.microsecond == 0)
            self.assert_(stamp.nanosecond == 0)

        ts = Timestamp('20000101 01:01:01')
        val = ts.value
        days = (ts - Timestamp('1970-01-01')).days

        check(val)
        check(val/long(1000),unit='us')
        check(val/long(1000000),unit='ms')
        check(val/long(1000000000),unit='s')
        check(days,unit='D',h=0)

        # using truediv, so these are like floats
        if compat.PY3:
            check((val+500000)/long(1000000000),unit='s',us=500)
            check((val+500000000)/long(1000000000),unit='s',us=500000)
            check((val+500000)/long(1000000),unit='ms',us=500)

        # get chopped in py2
        else:
            check((val+500000)/long(1000000000),unit='s')
            check((val+500000000)/long(1000000000),unit='s')
            check((val+500000)/long(1000000),unit='ms')

        # ok
        check((val+500000)/long(1000),unit='us',us=500)
        check((val+500000000)/long(1000000),unit='ms',us=500000)

        # floats
        check(val/1000.0 + 5,unit='us',us=5)
        check(val/1000.0 + 5000,unit='us',us=5000)
        check(val/1000000.0 + 0.5,unit='ms',us=500)
        check(val/1000000.0 + 0.005,unit='ms',us=5)
        check(val/1000000000.0 + 0.5,unit='s',us=500000)
        check(days + 0.5,unit='D',h=12)

        # nan
        result = Timestamp(np.nan)
        self.assert_(result is NaT)

        result = Timestamp(None)
        self.assert_(result is NaT)

        result = Timestamp(iNaT)
        self.assert_(result is NaT)

        result = Timestamp(NaT)
        self.assert_(result is NaT)

    def test_comparison(self):
        # 5-18-2012 00:00:00.000
        stamp = long(1337299200000000000)

        val = Timestamp(stamp)

        self.assert_(val == val)
        self.assert_(not val != val)
        self.assert_(not val < val)
        self.assert_(val <= val)
        self.assert_(not val > val)
        self.assert_(val >= val)

        other = datetime(2012, 5, 18)
        self.assert_(val == other)
        self.assert_(not val != other)
        self.assert_(not val < other)
        self.assert_(val <= other)
        self.assert_(not val > other)
        self.assert_(val >= other)

        other = Timestamp(stamp + 100)

        self.assert_(not val == other)
        self.assert_(val != other)
        self.assert_(val < other)
        self.assert_(val <= other)
        self.assert_(other > val)
        self.assert_(other >= val)

    def test_cant_compare_tz_naive_w_aware(self):
        _skip_if_no_pytz()
        # #1404
        a = Timestamp('3/12/2012')
        b = Timestamp('3/12/2012', tz='utc')

        self.assertRaises(Exception, a.__eq__, b)
        self.assertRaises(Exception, a.__ne__, b)
        self.assertRaises(Exception, a.__lt__, b)
        self.assertRaises(Exception, a.__gt__, b)
        self.assertRaises(Exception, b.__eq__, a)
        self.assertRaises(Exception, b.__ne__, a)
        self.assertRaises(Exception, b.__lt__, a)
        self.assertRaises(Exception, b.__gt__, a)

        if sys.version_info < (3, 3):
            self.assertRaises(Exception, a.__eq__, b.to_pydatetime())
            self.assertRaises(Exception, a.to_pydatetime().__eq__, b)
        else:
            self.assertFalse(a == b.to_pydatetime())
            self.assertFalse(a.to_pydatetime() == b)

    def test_delta_preserve_nanos(self):
        val = Timestamp(long(1337299200000000123))
        result = val + timedelta(1)
        self.assert_(result.nanosecond == val.nanosecond)

    def test_frequency_misc(self):
        self.assertEquals(fmod.get_freq_group('T'),
                          fmod.FreqGroup.FR_MIN)

        code, stride = fmod.get_freq_code(offsets.Hour())
        self.assertEquals(code, fmod.FreqGroup.FR_HR)

        code, stride = fmod.get_freq_code((5, 'T'))
        self.assertEquals(code, fmod.FreqGroup.FR_MIN)
        self.assertEquals(stride, 5)

        offset = offsets.Hour()
        result = fmod.to_offset(offset)
        self.assertEquals(result, offset)

        result = fmod.to_offset((5, 'T'))
        expected = offsets.Minute(5)
        self.assertEquals(result, expected)

        self.assertRaises(ValueError, fmod.get_freq_code, (5, 'baz'))

        self.assertRaises(ValueError, fmod.to_offset, '100foo')

        self.assertRaises(ValueError, fmod.to_offset, ('', ''))

        result = fmod.get_standard_freq(offsets.Hour())
        self.assertEquals(result, 'H')

    def test_hash_equivalent(self):
        d = {datetime(2011, 1, 1): 5}
        stamp = Timestamp(datetime(2011, 1, 1))
        self.assertEquals(d[stamp], 5)

    def test_timestamp_compare_scalars(self):
        # case where ndim == 0
        lhs = np.datetime64(datetime(2013, 12, 6))
        rhs = Timestamp('now')
        nat = Timestamp('nat')

        ops = {'gt': 'lt', 'lt': 'gt', 'ge': 'le', 'le': 'ge', 'eq': 'eq',
            'ne': 'ne'}

        for left, right in ops.items():
            left_f = getattr(operator, left)
            right_f = getattr(operator, right)

            if pd._np_version_under1p7:
                # you have to convert to timestamp for this to work with numpy
                # scalars
                expected = left_f(Timestamp(lhs), rhs)

                # otherwise a TypeError is thrown
                if left not in ('eq', 'ne'):
                    with tm.assertRaises(TypeError):
                        left_f(lhs, rhs)
            else:
                expected = left_f(lhs, rhs)

            result = right_f(rhs, lhs)
            self.assertEqual(result, expected)

            expected = left_f(rhs, nat)
            result = right_f(nat, rhs)
            self.assertEqual(result, expected)

    def test_timestamp_compare_series(self):
        # make sure we can compare Timestamps on the right AND left hand side
        # GH4982
        s = Series(date_range('20010101', periods=10), name='dates')
        s_nat = s.copy(deep=True)

        s[0] = pd.Timestamp('nat')
        s[3] = pd.Timestamp('nat')

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


class TestSlicing(tm.TestCase):

    def test_slice_year(self):
        dti = DatetimeIndex(freq='B', start=datetime(2005, 1, 1), periods=500)

        s = Series(np.arange(len(dti)), index=dti)
        result = s['2005']
        expected = s[s.index.year == 2005]
        assert_series_equal(result, expected)

        df = DataFrame(np.random.rand(len(dti), 5), index=dti)
        result = df.ix['2005']
        expected = df[df.index.year == 2005]
        assert_frame_equal(result, expected)

        rng = date_range('1/1/2000', '1/1/2010')

        result = rng.get_loc('2009')
        expected = slice(3288, 3653)
        self.assert_(result == expected)

    def test_slice_quarter(self):
        dti = DatetimeIndex(freq='D', start=datetime(2000, 6, 1), periods=500)

        s = Series(np.arange(len(dti)), index=dti)
        self.assertEquals(len(s['2001Q1']), 90)

        df = DataFrame(np.random.rand(len(dti), 5), index=dti)
        self.assertEquals(len(df.ix['1Q01']), 90)

    def test_slice_month(self):
        dti = DatetimeIndex(freq='D', start=datetime(2005, 1, 1), periods=500)
        s = Series(np.arange(len(dti)), index=dti)
        self.assertEquals(len(s['2005-11']), 30)

        df = DataFrame(np.random.rand(len(dti), 5), index=dti)
        self.assertEquals(len(df.ix['2005-11']), 30)

        assert_series_equal(s['2005-11'], s['11-2005'])

    def test_partial_slice(self):
        rng = DatetimeIndex(freq='D', start=datetime(2005, 1, 1), periods=500)
        s = Series(np.arange(len(rng)), index=rng)

        result = s['2005-05':'2006-02']
        expected = s['20050501':'20060228']
        assert_series_equal(result, expected)

        result = s['2005-05':]
        expected = s['20050501':]
        assert_series_equal(result, expected)

        result = s[:'2006-02']
        expected = s[:'20060228']
        assert_series_equal(result, expected)

        result = s['2005-1-1']
        self.assert_(result == s.irow(0))

        self.assertRaises(Exception, s.__getitem__, '2004-12-31')

    def test_partial_slice_daily(self):
        rng = DatetimeIndex(freq='H', start=datetime(2005, 1, 31), periods=500)
        s = Series(np.arange(len(rng)), index=rng)

        result = s['2005-1-31']
        assert_series_equal(result, s.ix[:24])

        self.assertRaises(Exception, s.__getitem__, '2004-12-31 00')

    def test_partial_slice_hourly(self):
        rng = DatetimeIndex(freq='T', start=datetime(2005, 1, 1, 20, 0, 0),
                            periods=500)
        s = Series(np.arange(len(rng)), index=rng)

        result = s['2005-1-1']
        assert_series_equal(result, s.ix[:60 * 4])

        result = s['2005-1-1 20']
        assert_series_equal(result, s.ix[:60])

        self.assert_(s['2005-1-1 20:00'] == s.ix[0])
        self.assertRaises(Exception, s.__getitem__, '2004-12-31 00:15')

    def test_partial_slice_minutely(self):
        rng = DatetimeIndex(freq='S', start=datetime(2005, 1, 1, 23, 59, 0),
                            periods=500)
        s = Series(np.arange(len(rng)), index=rng)

        result = s['2005-1-1 23:59']
        assert_series_equal(result, s.ix[:60])

        result = s['2005-1-1']
        assert_series_equal(result, s.ix[:60])

        self.assert_(s[Timestamp('2005-1-1 23:59:00')] == s.ix[0])
        self.assertRaises(Exception, s.__getitem__, '2004-12-31 00:00:00')

    def test_partial_slicing_with_multiindex(self):

        # GH 4758
        # partial string indexing with a multi-index buggy
        df = DataFrame({'ACCOUNT':["ACCT1", "ACCT1", "ACCT1", "ACCT2"],
                        'TICKER':["ABC", "MNP", "XYZ", "XYZ"],
                        'val':[1,2,3,4]},
                       index=date_range("2013-06-19 09:30:00", periods=4, freq='5T'))
        df_multi = df.set_index(['ACCOUNT', 'TICKER'], append=True)

        expected = DataFrame([[1]],index=Index(['ABC'],name='TICKER'),columns=['val'])
        result = df_multi.loc[('2013-06-19 09:30:00', 'ACCT1')]
        assert_frame_equal(result, expected)

        expected = df_multi.loc[(pd.Timestamp('2013-06-19 09:30:00', tz=None), 'ACCT1', 'ABC')]
        result = df_multi.loc[('2013-06-19 09:30:00', 'ACCT1', 'ABC')]
        assert_series_equal(result, expected)

        # this is a KeyError as we don't do partial string selection on multi-levels
        def f():
            df_multi.loc[('2013-06-19', 'ACCT1', 'ABC')]
        self.assertRaises(KeyError, f)

        # GH 4294
        # partial slice on a series mi
        s = pd.DataFrame(randn(1000, 1000), index=pd.date_range('2000-1-1', periods=1000)).stack()

        s2 = s[:-1].copy()
        expected = s2['2000-1-4']
        result = s2[pd.Timestamp('2000-1-4')]
        assert_series_equal(result, expected)

        result = s[pd.Timestamp('2000-1-4')]
        expected = s['2000-1-4']
        assert_series_equal(result, expected)

        df2 = pd.DataFrame(s)
        expected = df2.ix['2000-1-4']
        result = df2.ix[pd.Timestamp('2000-1-4')]
        assert_frame_equal(result, expected)

    def test_date_range_normalize(self):
        snap = datetime.today()
        n = 50

        rng = date_range(snap, periods=n, normalize=False, freq='2D')

        offset = timedelta(2)
        values = np.array([snap + i * offset for i in range(n)],
                          dtype='M8[ns]')

        self.assert_(np.array_equal(rng, values))

        rng = date_range(
            '1/1/2000 08:15', periods=n, normalize=False, freq='B')
        the_time = time(8, 15)
        for val in rng:
            self.assert_(val.time() == the_time)

    def test_timedelta(self):
        # this is valid too
        index = date_range('1/1/2000', periods=50, freq='B')
        shifted = index + timedelta(1)
        back = shifted + timedelta(-1)
        self.assert_(tm.equalContents(index, back))
        self.assertEqual(shifted.freq, index.freq)
        self.assertEqual(shifted.freq, back.freq)

        result = index - timedelta(1)
        expected = index + timedelta(-1)
        self.assert_(result.equals(expected))

        # GH4134, buggy with timedeltas
        rng = date_range('2013', '2014')
        s = Series(rng)
        result1 = rng - pd.offsets.Hour(1)
        result2 = DatetimeIndex(s - np.timedelta64(100000000))
        result3 = rng - np.timedelta64(100000000)
        result4 = DatetimeIndex(s - pd.offsets.Hour(1))
        self.assert_(result1.equals(result4))
        self.assert_(result2.equals(result3))

    def test_shift(self):
        ts = Series(np.random.randn(5),
                    index=date_range('1/1/2000', periods=5, freq='H'))

        result = ts.shift(1, freq='5T')
        exp_index = ts.index.shift(1, freq='5T')
        self.assert_(result.index.equals(exp_index))

        # GH #1063, multiple of same base
        result = ts.shift(1, freq='4H')
        exp_index = ts.index + datetools.Hour(4)
        self.assert_(result.index.equals(exp_index))

        idx = DatetimeIndex(['2000-01-01', '2000-01-02', '2000-01-04'])
        self.assertRaises(ValueError, idx.shift, 1)

    def test_setops_preserve_freq(self):
        rng = date_range('1/1/2000', '1/1/2002')

        result = rng[:50].union(rng[50:100])
        self.assert_(result.freq == rng.freq)

        result = rng[:50].union(rng[30:100])
        self.assert_(result.freq == rng.freq)

        result = rng[:50].union(rng[60:100])
        self.assert_(result.freq is None)

        result = rng[:50].intersection(rng[25:75])
        self.assert_(result.freqstr == 'D')

        nofreq = DatetimeIndex(list(rng[25:75]))
        result = rng[:50].union(nofreq)
        self.assert_(result.freq == rng.freq)

        result = rng[:50].intersection(nofreq)
        self.assert_(result.freq == rng.freq)

    def test_min_max(self):
        rng = date_range('1/1/2000', '12/31/2000')
        rng2 = rng.take(np.random.permutation(len(rng)))

        the_min = rng2.min()
        the_max = rng2.max()
        tm.assert_isinstance(the_min, Timestamp)
        tm.assert_isinstance(the_max, Timestamp)
        self.assertEqual(the_min, rng[0])
        self.assertEqual(the_max, rng[-1])

        self.assertEqual(rng.min(), rng[0])
        self.assertEqual(rng.max(), rng[-1])

    def test_min_max_series(self):
        rng = date_range('1/1/2000', periods=10, freq='4h')
        lvls = ['A', 'A', 'A', 'B', 'B', 'B', 'C', 'C', 'C', 'C']
        df = DataFrame({'TS': rng, 'V': np.random.randn(len(rng)),
                        'L': lvls})

        result = df.TS.max()
        exp = Timestamp(df.TS.iget(-1))
        self.assertTrue(isinstance(result, Timestamp))
        self.assertEqual(result, exp)

        result = df.TS.min()
        exp = Timestamp(df.TS.iget(0))
        self.assertTrue(isinstance(result, Timestamp))
        self.assertEqual(result, exp)

    def test_from_M8_structured(self):
        dates = [(datetime(2012, 9, 9, 0, 0),
                 datetime(2012, 9, 8, 15, 10))]
        arr = np.array(dates,
                       dtype=[('Date', 'M8[us]'), ('Forecasting', 'M8[us]')])
        df = DataFrame(arr)

        self.assertEqual(df['Date'][0], dates[0][0])
        self.assertEqual(df['Forecasting'][0], dates[0][1])

        s = Series(arr['Date'])
        self.assertTrue(s[0], Timestamp)
        self.assertEqual(s[0], dates[0][0])

        s = Series.from_array(arr['Date'], Index([0]))
        self.assertEqual(s[0], dates[0][0])

    def test_get_level_values_box(self):
        from pandas import MultiIndex

        dates = date_range('1/1/2000', periods=4)
        levels = [dates, [0, 1]]
        labels = [[0, 0, 1, 1, 2, 2, 3, 3],
                  [0, 1, 0, 1, 0, 1, 0, 1]]

        index = MultiIndex(levels=levels, labels=labels)

        self.assertTrue(isinstance(index.get_level_values(0)[0], Timestamp))

    def test_frame_apply_dont_convert_datetime64(self):
        from pandas.tseries.offsets import BDay
        df = DataFrame({'x1': [datetime(1996, 1, 1)]})

        df = df.applymap(lambda x: x + BDay())
        df = df.applymap(lambda x: x + BDay())

        self.assertTrue(df.x1.dtype == 'M8[ns]')

    def test_date_range_fy5252(self):
        dr = date_range(start="2013-01-01",
                           periods=2,
                           freq=offsets.FY5253(startingMonth=1,
                                               weekday=3,
                                               variation="nearest"))
        self.assertEqual(dr[0], Timestamp('2013-01-31'))
        self.assertEqual(dr[1], Timestamp('2014-01-30'))

class TimeConversionFormats(tm.TestCase):
    def test_to_datetime_format(self):
        values = ['1/1/2000', '1/2/2000', '1/3/2000']

        results1 = [ Timestamp('20000101'), Timestamp('20000201'),
                     Timestamp('20000301') ]
        results2 = [ Timestamp('20000101'), Timestamp('20000102'),
                     Timestamp('20000103') ]
        for vals, expecteds in [ (values, (Index(results1), Index(results2))),
                                 (Series(values),(Series(results1), Series(results2))),
                                 (values[0], (results1[0], results2[0])),
                                 (values[1], (results1[1], results2[1])),
                                 (values[2], (results1[2], results2[2])) ]:

            for i, fmt in enumerate(['%d/%m/%Y', '%m/%d/%Y']):
                result = to_datetime(vals, format=fmt)
                expected = expecteds[i]

                if isinstance(expected, Series):
                    assert_series_equal(result, Series(expected))
                elif isinstance(expected, Timestamp):
                    self.assert_(result == expected)
                else:
                    self.assert_(result.equals(expected))

    def test_to_datetime_format_YYYYMMDD(self):
        s = Series([19801222,19801222] + [19810105]*5)
        expected = Series([ Timestamp(x) for x in s.apply(str) ])

        result = to_datetime(s,format='%Y%m%d')
        assert_series_equal(result, expected)

        result = to_datetime(s.apply(str),format='%Y%m%d')
        assert_series_equal(result, expected)

        # with NaT
        expected = Series([Timestamp("19801222"),Timestamp("19801222")] + [Timestamp("19810105")]*5)
        expected[2] = np.nan
        s[2] = np.nan

        result = to_datetime(s,format='%Y%m%d')
        assert_series_equal(result, expected)

        # string with NaT
        s = s.apply(str)
        s[2] = 'nat'
        result = to_datetime(s,format='%Y%m%d')
        assert_series_equal(result, expected)


    def test_to_datetime_format_microsecond(self):
        val = '01-Apr-2011 00:00:01.978'
        format = '%d-%b-%Y %H:%M:%S.%f'
        result = to_datetime(val, format=format)
        exp = dt.datetime.strptime(val, format)
        self.assert_(result == exp)

    def test_to_datetime_format_time(self):
        data = [
                ['01/10/2010 15:20', '%m/%d/%Y %H:%M', Timestamp('2010-01-10 15:20')],
	            ['01/10/2010 05:43', '%m/%d/%Y %I:%M', Timestamp('2010-01-10 05:43')],
	            ['01/10/2010 13:56:01', '%m/%d/%Y %H:%M:%S', Timestamp('2010-01-10 13:56:01')]#,
	            #['01/10/2010 08:14 PM', '%m/%d/%Y %I:%M %p', Timestamp('2010-01-10 20:14')],
	            #['01/10/2010 07:40 AM', '%m/%d/%Y %I:%M %p', Timestamp('2010-01-10 07:40')],
	            #['01/10/2010 09:12:56 AM', '%m/%d/%Y %I:%M:%S %p', Timestamp('2010-01-10 09:12:56')]
            ]
        for s, format, dt in data:
            self.assertEqual(to_datetime(s, format=format), dt)

    def test_to_datetime_format_weeks(self):
        data = [
                ['2009324', '%Y%W%w', Timestamp('2009-08-13')],
                ['2013020', '%Y%U%w', Timestamp('2013-01-13')]
            ]
        for s, format, dt in data:
            self.assertEqual(to_datetime(s, format=format), dt)

class TestToDatetimeInferFormat(tm.TestCase):
    def test_to_datetime_infer_datetime_format_consistent_format(self):
        time_series = pd.Series(
            pd.date_range('20000101', periods=50, freq='H')
        )

        test_formats = [
            '%m-%d-%Y',
            '%m/%d/%Y %H:%M:%S.%f',
            '%Y-%m-%dT%H:%M:%S.%f',
        ]

        for test_format in test_formats:
            s_as_dt_strings = time_series.apply(
                lambda x: x.strftime(test_format)
            )

            with_format = pd.to_datetime(s_as_dt_strings, format=test_format)
            no_infer = pd.to_datetime(
                s_as_dt_strings, infer_datetime_format=False
            )
            yes_infer = pd.to_datetime(
                s_as_dt_strings, infer_datetime_format=True
            )

            # Whether the format is explicitly passed, it is inferred, or
            # it is not inferred, the results should all be the same
            self.assert_(np.array_equal(with_format, no_infer))
            self.assert_(np.array_equal(no_infer, yes_infer))

    def test_to_datetime_infer_datetime_format_inconsistent_format(self):
        test_series = pd.Series(
            np.array([
                '01/01/2011 00:00:00',
                '01-02-2011 00:00:00',
                '2011-01-03T00:00:00',
        ]))

        # When the format is inconsistent, infer_datetime_format should just
        # fallback to the default parsing
        self.assert_(np.array_equal(
            pd.to_datetime(test_series, infer_datetime_format=False),
            pd.to_datetime(test_series, infer_datetime_format=True)
        ))

        test_series = pd.Series(
            np.array([
                'Jan/01/2011',
                'Feb/01/2011',
                'Mar/01/2011',
        ]))

        self.assert_(np.array_equal(
            pd.to_datetime(test_series, infer_datetime_format=False),
            pd.to_datetime(test_series, infer_datetime_format=True)
        ))

    def test_to_datetime_infer_datetime_format_series_with_nans(self):
        test_series = pd.Series(
            np.array([
                '01/01/2011 00:00:00',
                np.nan,
                '01/03/2011 00:00:00',
                np.nan,
        ]))

        self.assert_(np.array_equal(
            pd.to_datetime(test_series, infer_datetime_format=False),
            pd.to_datetime(test_series, infer_datetime_format=True)
        ))

    def test_to_datetime_infer_datetime_format_series_starting_with_nans(self):
        test_series = pd.Series(
            np.array([
                np.nan,
                np.nan,
                '01/01/2011 00:00:00',
                '01/02/2011 00:00:00',
                '01/03/2011 00:00:00',
        ]))

        self.assert_(np.array_equal(
            pd.to_datetime(test_series, infer_datetime_format=False),
            pd.to_datetime(test_series, infer_datetime_format=True)
        ))


class TestGuessDatetimeFormat(tm.TestCase):
    def test_guess_datetime_format_with_parseable_formats(self):
        dt_string_to_format = (
            ('20111230', '%Y%m%d'),
            ('2011-12-30', '%Y-%m-%d'),
            ('30-12-2011', '%d-%m-%Y'),
            ('2011-12-30 00:00:00', '%Y-%m-%d %H:%M:%S'),
            ('2011-12-30T00:00:00', '%Y-%m-%dT%H:%M:%S'),
            ('2011-12-30 00:00:00.000000', '%Y-%m-%d %H:%M:%S.%f'),
        )

        for dt_string, dt_format in dt_string_to_format:
            self.assertEquals(
                tools._guess_datetime_format(dt_string),
                dt_format
            )

    def test_guess_datetime_format_with_dayfirst(self):
        ambiguous_string = '01/01/2011'
        self.assertEquals(
            tools._guess_datetime_format(ambiguous_string, dayfirst=True),
            '%d/%m/%Y'
        )
        self.assertEquals(
            tools._guess_datetime_format(ambiguous_string, dayfirst=False),
            '%m/%d/%Y'
        )

    def test_guess_datetime_format_with_locale_specific_formats(self):
        # The month names will vary depending on the locale, in which
        # case these wont be parsed properly (dateutil can't parse them)
        _skip_if_has_locale()

        dt_string_to_format = (
            ('30/Dec/2011', '%d/%b/%Y'),
            ('30/December/2011', '%d/%B/%Y'),
            ('30/Dec/2011 00:00:00', '%d/%b/%Y %H:%M:%S'),
        )

        for dt_string, dt_format in dt_string_to_format:
            self.assertEquals(
                tools._guess_datetime_format(dt_string),
                dt_format
            )

    def test_guess_datetime_format_invalid_inputs(self):
        # A datetime string must include a year, month and a day for it
        # to be guessable, in addition to being a string that looks like
        # a datetime
        invalid_dts = [
            '2013',
            '01/2013',
            '12:00:00',
            '1/1/1/1',
            'this_is_not_a_datetime',
            '51a',
            9,
            datetime(2011, 1, 1),
        ]

        for invalid_dt in invalid_dts:
            self.assertTrue(tools._guess_datetime_format(invalid_dt) is None)

    def test_guess_datetime_format_for_array(self):
        expected_format = '%Y-%m-%d %H:%M:%S.%f'
        dt_string = datetime(2011, 12, 30, 0, 0, 0).strftime(expected_format)

        test_arrays = [
            np.array([dt_string, dt_string, dt_string], dtype='O'),
            np.array([np.nan, np.nan, dt_string], dtype='O'),
            np.array([dt_string, 'random_string'], dtype='O'),
        ]

        for test_array in test_arrays:
            self.assertEqual(
                tools._guess_datetime_format_for_array(test_array),
                expected_format
            )

        format_for_string_of_nans = tools._guess_datetime_format_for_array(
            np.array([np.nan, np.nan, np.nan], dtype='O')
        )
        self.assertTrue(format_for_string_of_nans is None)

if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
