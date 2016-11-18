# pylint: disable-msg=E1101,W0612
import calendar
import operator
import sys
import warnings
from datetime import datetime, time, timedelta
from numpy.random import rand

import nose
import numpy as np
import pandas.index as _index
import pandas.lib as lib
import pandas.tslib as tslib

from pandas.types.common import is_datetime64_ns_dtype
import pandas as pd
import pandas.compat as compat
import pandas.core.common as com
import pandas.tseries.frequencies as frequencies
import pandas.tseries.offsets as offsets
import pandas.tseries.tools as tools
import pandas.util.testing as tm
from pandas import (
    Index, Series, DataFrame, isnull, date_range, Timestamp, Period,
    DatetimeIndex, Int64Index, to_datetime, bdate_range, Float64Index,
    NaT, timedelta_range, Timedelta, _np_version_under1p8, concat)
from pandas.compat import range, long, StringIO, lrange, lmap, zip, product
from pandas.compat.numpy import np_datetime64_compat
from pandas.core.common import PerformanceWarning
from pandas.tslib import iNaT
from pandas.util.testing import (
    assert_frame_equal, assert_series_equal, assert_almost_equal,
    _skip_if_has_locale, slow)

randn = np.random.randn


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
        tm.assertIsInstance(self.dups, Series)
        tm.assertIsInstance(self.dups.index, DatetimeIndex)

    def test_is_unique_monotonic(self):
        self.assertFalse(self.dups.index.is_unique)

    def test_index_unique(self):
        uniques = self.dups.index.unique()
        expected = DatetimeIndex([datetime(2000, 1, 2), datetime(2000, 1, 3),
                                  datetime(2000, 1, 4), datetime(2000, 1, 5)])
        self.assertEqual(uniques.dtype, 'M8[ns]')  # sanity
        tm.assert_index_equal(uniques, expected)
        self.assertEqual(self.dups.index.nunique(), 4)

        # #2563
        self.assertTrue(isinstance(uniques, DatetimeIndex))

        dups_local = self.dups.index.tz_localize('US/Eastern')
        dups_local.name = 'foo'
        result = dups_local.unique()
        expected = DatetimeIndex(expected, name='foo')
        expected = expected.tz_localize('US/Eastern')
        self.assertTrue(result.tz is not None)
        self.assertEqual(result.name, 'foo')
        tm.assert_index_equal(result, expected)

        # NaT, note this is excluded
        arr = [1370745748 + t for t in range(20)] + [iNaT]
        idx = DatetimeIndex(arr * 3)
        tm.assert_index_equal(idx.unique(), DatetimeIndex(arr))
        self.assertEqual(idx.nunique(), 20)
        self.assertEqual(idx.nunique(dropna=False), 21)

        arr = [Timestamp('2013-06-09 02:42:28') + timedelta(seconds=t)
               for t in range(20)] + [NaT]
        idx = DatetimeIndex(arr * 3)
        tm.assert_index_equal(idx.unique(), DatetimeIndex(arr))
        self.assertEqual(idx.nunique(), 20)
        self.assertEqual(idx.nunique(dropna=False), 21)

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
        ts[datetime(2000, 1, 6)] = 0
        self.assertEqual(ts[datetime(2000, 1, 6)], 0)

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
            self.assertIn(timestamp, df.index)

            # it works!
            df.ix[timestamp]
            self.assertTrue(len(df.ix[[timestamp]]) > 0)
        finally:
            _index._SIZE_CUTOFF = old_cutoff

    def test_indexing_unordered(self):
        # GH 2437
        rng = date_range(start='2011-01-01', end='2011-01-15')
        ts = Series(randn(len(rng)), index=rng)
        ts2 = concat([ts[0:4], ts[-4:], ts[4:-4]])

        for t in ts.index:
            # TODO: unused?
            s = str(t)  # noqa

            expected = ts[t]
            result = ts2[t]
            self.assertTrue(expected == result)

        # GH 3448 (ranges)
        def compare(slobj):
            result = ts2[slobj].copy()
            result = result.sort_index()
            expected = ts[slobj]
            assert_series_equal(result, expected)

        compare(slice('2011-01-01', '2011-01-15'))
        compare(slice('2010-12-30', '2011-01-15'))
        compare(slice('2011-01-01', '2011-01-16'))

        # partial ranges
        compare(slice('2011-01-01', '2011-01-6'))
        compare(slice('2011-01-06', '2011-01-8'))
        compare(slice('2011-01-06', '2011-01-12'))

        # single values
        result = ts2['2011'].sort_index()
        expected = ts['2011']
        assert_series_equal(result, expected)

        # diff freq
        rng = date_range(datetime(2005, 1, 1), periods=20, freq='M')
        ts = Series(np.arange(len(rng)), index=rng)
        ts = ts.take(np.random.permutation(20))

        result = ts['2005']
        for t in result.index:
            self.assertTrue(t.year == 2005)

    def test_indexing(self):

        idx = date_range("2001-1-1", periods=20, freq='M')
        ts = Series(np.random.rand(len(idx)), index=idx)

        # getting

        # GH 3070, make sure semantics work on Series/Frame
        expected = ts['2001']
        expected.name = 'A'

        df = DataFrame(dict(A=ts))
        result = df['2001']['A']
        assert_series_equal(expected, result)

        # setting
        ts['2001'] = 1
        expected = ts['2001']
        expected.name = 'A'

        df.loc['2001', 'A'] = 1

        result = df['2001']['A']
        assert_series_equal(expected, result)

        # GH3546 (not including times on the last day)
        idx = date_range(start='2013-05-31 00:00', end='2013-05-31 23:00',
                         freq='H')
        ts = Series(lrange(len(idx)), index=idx)
        expected = ts['2013-05']
        assert_series_equal(expected, ts)

        idx = date_range(start='2013-05-31 00:00', end='2013-05-31 23:59',
                         freq='S')
        ts = Series(lrange(len(idx)), index=idx)
        expected = ts['2013-05']
        assert_series_equal(expected, ts)

        idx = [Timestamp('2013-05-31 00:00'),
               Timestamp(datetime(2013, 5, 31, 23, 59, 59, 999999))]
        ts = Series(lrange(len(idx)), index=idx)
        expected = ts['2013']
        assert_series_equal(expected, ts)

        # GH 3925, indexing with a seconds resolution string / datetime object
        df = DataFrame(randn(5, 5),
                       columns=['open', 'high', 'low', 'close', 'volume'],
                       index=date_range('2012-01-02 18:01:00',
                                        periods=5, tz='US/Central', freq='s'))
        expected = df.loc[[df.index[2]]]
        result = df['2012-01-02 18:01:02']
        assert_frame_equal(result, expected)

        # this is a single date, so will raise
        self.assertRaises(KeyError, df.__getitem__, df.index[2], )

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


def assert_range_equal(left, right):
    assert (left.equals(right))
    assert (left.freq == right.freq)
    assert (left.tz == right.tz)


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

        self.assertEqual(v1, Timestamp('2/28/2005'))
        self.assertEqual(v2, Timestamp('4/30/2005'))
        self.assertEqual(v3, Timestamp('6/30/2005'))

        # don't carry freq through irregular slicing
        self.assertIsNone(dti2.freq)

    def test_pass_datetimeindex_to_index(self):
        # Bugs in #1396
        rng = date_range('1/1/2000', '3/1/2000')
        idx = Index(rng, dtype=object)

        expected = Index(rng.to_pydatetime(), dtype=object)

        self.assert_numpy_array_equal(idx.values, expected.values)

    def test_contiguous_boolean_preserve_freq(self):
        rng = date_range('1/1/2000', '3/1/2000', freq='B')

        mask = np.zeros(len(rng), dtype=bool)
        mask[10:20] = True

        masked = rng[mask]
        expected = rng[10:20]
        self.assertIsNotNone(expected.freq)
        assert_range_equal(masked, expected)

        mask[22] = True
        masked = rng[mask]
        self.assertIsNone(masked.freq)

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

        tm.assertIsInstance(s[5], Timestamp)

        rng = date_range('20090415', '20090519', freq='B')
        s = Series(rng, index=rng)
        tm.assertIsInstance(s[5], Timestamp)

        tm.assertIsInstance(s.iat[5], Timestamp)

    def test_series_box_timedelta(self):
        rng = timedelta_range('1 day 1 s', periods=5, freq='h')
        s = Series(rng)
        tm.assertIsInstance(s[1], Timedelta)
        tm.assertIsInstance(s.iat[2], Timedelta)

    def test_date_range_ambiguous_arguments(self):
        # #2538
        start = datetime(2011, 1, 1, 5, 3, 40)
        end = datetime(2011, 1, 1, 8, 9, 40)

        self.assertRaises(ValueError, date_range, start, end, freq='s',
                          periods=10)

    def test_timestamp_to_datetime(self):
        tm._skip_if_no_pytz()
        rng = date_range('20090415', '20090519', tz='US/Eastern')

        stamp = rng[0]
        dtval = stamp.to_pydatetime()
        self.assertEqual(stamp, dtval)
        self.assertEqual(stamp.tzinfo, dtval.tzinfo)

    def test_timestamp_to_datetime_dateutil(self):
        tm._skip_if_no_pytz()
        rng = date_range('20090415', '20090519', tz='dateutil/US/Eastern')

        stamp = rng[0]
        dtval = stamp.to_pydatetime()
        self.assertEqual(stamp, dtval)
        self.assertEqual(stamp.tzinfo, dtval.tzinfo)

    def test_timestamp_to_datetime_explicit_pytz(self):
        tm._skip_if_no_pytz()
        import pytz
        rng = date_range('20090415', '20090519',
                         tz=pytz.timezone('US/Eastern'))

        stamp = rng[0]
        dtval = stamp.to_pydatetime()
        self.assertEqual(stamp, dtval)
        self.assertEqual(stamp.tzinfo, dtval.tzinfo)

    def test_timestamp_to_datetime_explicit_dateutil(self):
        tm._skip_if_windows_python_3()
        tm._skip_if_no_dateutil()
        from pandas.tslib import _dateutil_gettz as gettz
        rng = date_range('20090415', '20090519', tz=gettz('US/Eastern'))

        stamp = rng[0]
        dtval = stamp.to_pydatetime()
        self.assertEqual(stamp, dtval)
        self.assertEqual(stamp.tzinfo, dtval.tzinfo)

    def test_index_convert_to_datetime_array(self):
        tm._skip_if_no_pytz()

        def _check_rng(rng):
            converted = rng.to_pydatetime()
            tm.assertIsInstance(converted, np.ndarray)
            for x, stamp in zip(converted, rng):
                tm.assertIsInstance(x, datetime)
                self.assertEqual(x, stamp.to_pydatetime())
                self.assertEqual(x.tzinfo, stamp.tzinfo)

        rng = date_range('20090415', '20090519')
        rng_eastern = date_range('20090415', '20090519', tz='US/Eastern')
        rng_utc = date_range('20090415', '20090519', tz='utc')

        _check_rng(rng)
        _check_rng(rng_eastern)
        _check_rng(rng_utc)

    def test_index_convert_to_datetime_array_explicit_pytz(self):
        tm._skip_if_no_pytz()
        import pytz

        def _check_rng(rng):
            converted = rng.to_pydatetime()
            tm.assertIsInstance(converted, np.ndarray)
            for x, stamp in zip(converted, rng):
                tm.assertIsInstance(x, datetime)
                self.assertEqual(x, stamp.to_pydatetime())
                self.assertEqual(x.tzinfo, stamp.tzinfo)

        rng = date_range('20090415', '20090519')
        rng_eastern = date_range('20090415', '20090519',
                                 tz=pytz.timezone('US/Eastern'))
        rng_utc = date_range('20090415', '20090519', tz=pytz.utc)

        _check_rng(rng)
        _check_rng(rng_eastern)
        _check_rng(rng_utc)

    def test_index_convert_to_datetime_array_dateutil(self):
        tm._skip_if_no_dateutil()
        import dateutil

        def _check_rng(rng):
            converted = rng.to_pydatetime()
            tm.assertIsInstance(converted, np.ndarray)
            for x, stamp in zip(converted, rng):
                tm.assertIsInstance(x, datetime)
                self.assertEqual(x, stamp.to_pydatetime())
                self.assertEqual(x.tzinfo, stamp.tzinfo)

        rng = date_range('20090415', '20090519')
        rng_eastern = date_range('20090415', '20090519',
                                 tz='dateutil/US/Eastern')
        rng_utc = date_range('20090415', '20090519', tz=dateutil.tz.tzutc())

        _check_rng(rng)
        _check_rng(rng_eastern)
        _check_rng(rng_utc)

    def test_ctor_str_intraday(self):
        rng = DatetimeIndex(['1-1-2000 00:00:01'])
        self.assertEqual(rng[0].second, 1)

    def test_series_ctor_plus_datetimeindex(self):
        rng = date_range('20090415', '20090519', freq='B')
        data = dict((k, 1) for k in rng)

        result = Series(data, index=rng)
        self.assertIs(result.index, rng)

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
                                freq=offsets.BDay())
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

        # neither monotonic increasing or decreasing
        rng2 = rng[[1, 0, 2]]

        self.assertRaises(ValueError, rng2.get_indexer, rng, method='pad')

    def test_frame_ctor_datetime64_column(self):
        rng = date_range('1/1/2000 00:00:00', '1/1/2000 1:59:50', freq='10s')
        dates = np.asarray(rng)

        df = DataFrame({'A': np.random.randn(len(rng)), 'B': dates})
        self.assertTrue(np.issubdtype(df['B'].dtype, np.dtype('M8[ns]')))

    def test_frame_add_datetime64_column(self):
        rng = date_range('1/1/2000 00:00:00', '1/1/2000 1:59:50', freq='10s')
        df = DataFrame(index=np.arange(len(rng)))

        df['A'] = rng
        self.assertTrue(np.issubdtype(df['A'].dtype, np.dtype('M8[ns]')))

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

            ex_vals = to_datetime(vals.astype('O')).values

            self.assertEqual(df[unit].dtype, ns_dtype)
            self.assertTrue((df[unit].values == ex_vals).all())

        # Test insertion into existing datetime64 column
        df = DataFrame({'ints': np.arange(n)}, index=np.arange(n))
        df['dates'] = np.arange(n, dtype=np.int64).view(ns_dtype)

        for unit in units:
            dtype = np.dtype('M8[%s]' % unit)
            vals = np.arange(n, dtype=np.int64).view(dtype)

            tmp = df.copy()

            tmp['dates'] = vals
            ex_vals = to_datetime(vals.astype('O')).values

            self.assertTrue((tmp['dates'].values == ex_vals).all())

    def test_to_datetime_unit(self):

        epoch = 1370745748
        s = Series([epoch + t for t in range(20)])
        result = to_datetime(s, unit='s')
        expected = Series([Timestamp('2013-06-09 02:42:28') + timedelta(
            seconds=t) for t in range(20)])
        assert_series_equal(result, expected)

        s = Series([epoch + t for t in range(20)]).astype(float)
        result = to_datetime(s, unit='s')
        expected = Series([Timestamp('2013-06-09 02:42:28') + timedelta(
            seconds=t) for t in range(20)])
        assert_series_equal(result, expected)

        s = Series([epoch + t for t in range(20)] + [iNaT])
        result = to_datetime(s, unit='s')
        expected = Series([Timestamp('2013-06-09 02:42:28') + timedelta(
            seconds=t) for t in range(20)] + [NaT])
        assert_series_equal(result, expected)

        s = Series([epoch + t for t in range(20)] + [iNaT]).astype(float)
        result = to_datetime(s, unit='s')
        expected = Series([Timestamp('2013-06-09 02:42:28') + timedelta(
            seconds=t) for t in range(20)] + [NaT])
        assert_series_equal(result, expected)

        # GH13834
        s = Series([epoch + t for t in np.arange(0, 2, .25)] +
                   [iNaT]).astype(float)
        result = to_datetime(s, unit='s')
        expected = Series([Timestamp('2013-06-09 02:42:28') + timedelta(
            seconds=t) for t in np.arange(0, 2, .25)] + [NaT])
        assert_series_equal(result, expected)

        s = concat([Series([epoch + t for t in range(20)]
                           ).astype(float), Series([np.nan])],
                   ignore_index=True)
        result = to_datetime(s, unit='s')
        expected = Series([Timestamp('2013-06-09 02:42:28') + timedelta(
            seconds=t) for t in range(20)] + [NaT])
        assert_series_equal(result, expected)

        result = to_datetime([1, 2, 'NaT', pd.NaT, np.nan], unit='D')
        expected = DatetimeIndex([Timestamp('1970-01-02'),
                                  Timestamp('1970-01-03')] + ['NaT'] * 3)
        tm.assert_index_equal(result, expected)

        with self.assertRaises(ValueError):
            to_datetime([1, 2, 'foo'], unit='D')
        with self.assertRaises(ValueError):
            to_datetime([1, 2, 111111111], unit='D')

        # coerce we can process
        expected = DatetimeIndex([Timestamp('1970-01-02'),
                                  Timestamp('1970-01-03')] + ['NaT'] * 1)
        result = to_datetime([1, 2, 'foo'], unit='D', errors='coerce')
        tm.assert_index_equal(result, expected)

        result = to_datetime([1, 2, 111111111], unit='D', errors='coerce')
        tm.assert_index_equal(result, expected)

    def test_series_ctor_datetime64(self):
        rng = date_range('1/1/2000 00:00:00', '1/1/2000 1:59:50', freq='10s')
        dates = np.asarray(rng)

        series = Series(dates)
        self.assertTrue(np.issubdtype(series.dtype, np.dtype('M8[ns]')))

    def test_index_cast_datetime64_other_units(self):
        arr = np.arange(0, 100, 10, dtype=np.int64).view('M8[D]')

        idx = Index(arr)

        self.assertTrue((idx.values == tslib.cast_to_nanoseconds(arr)).all())

    def test_reindex_series_add_nat(self):
        rng = date_range('1/1/2000 00:00:00', periods=10, freq='10s')
        series = Series(rng)

        result = series.reindex(lrange(15))
        self.assertTrue(np.issubdtype(result.dtype, np.dtype('M8[ns]')))

        mask = result.isnull()
        self.assertTrue(mask[-5:].all())
        self.assertFalse(mask[:-5].any())

    def test_reindex_frame_add_nat(self):
        rng = date_range('1/1/2000 00:00:00', periods=10, freq='10s')
        df = DataFrame({'A': np.random.randn(len(rng)), 'B': rng})

        result = df.reindex(lrange(15))
        self.assertTrue(np.issubdtype(result['B'].dtype, np.dtype('M8[ns]')))

        mask = com.isnull(result)['B']
        self.assertTrue(mask[-5:].all())
        self.assertFalse(mask[:-5].any())

    def test_series_repr_nat(self):
        series = Series([0, 1000, 2000, iNaT], dtype='M8[ns]')

        result = repr(series)
        expected = ('0   1970-01-01 00:00:00.000000\n'
                    '1   1970-01-01 00:00:00.000001\n'
                    '2   1970-01-01 00:00:00.000002\n'
                    '3                          NaT\n'
                    'dtype: datetime64[ns]')
        self.assertEqual(result, expected)

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
        tm.assertIsInstance(result2, DatetimeIndex)
        tm.assert_numpy_array_equal(result, result2.values)

        malformed = np.array(['1/100/2000', np.nan], dtype=object)

        # GH 10636, default is now 'raise'
        self.assertRaises(ValueError,
                          lambda: to_datetime(malformed, errors='raise'))

        result = to_datetime(malformed, errors='ignore')
        tm.assert_numpy_array_equal(result, malformed)

        self.assertRaises(ValueError, to_datetime, malformed, errors='raise')

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

        assert_series_equal(result, expected, check_names=False)
        self.assertEqual(result.name, 'foo')

        assert_series_equal(dresult, expected, check_names=False)
        self.assertEqual(dresult.name, 'foo')

    def test_to_datetime_iso8601(self):
        result = to_datetime(["2012-01-01 00:00:00"])
        exp = Timestamp("2012-01-01 00:00:00")
        self.assertEqual(result[0], exp)

        result = to_datetime(['20121001'])  # bad iso 8601
        exp = Timestamp('2012-10-01')
        self.assertEqual(result[0], exp)

    def test_to_datetime_default(self):
        rs = to_datetime('2001')
        xp = datetime(2001, 1, 1)
        self.assertTrue(rs, xp)

        # dayfirst is essentially broken

        # to_datetime('01-13-2012', dayfirst=True)
        # self.assertRaises(ValueError, to_datetime('01-13-2012',
        #                   dayfirst=True))

    def test_to_datetime_on_datetime64_series(self):
        # #2699
        s = Series(date_range('1/1/2000', periods=10))

        result = to_datetime(s)
        self.assertEqual(result[0], s[0])

    def test_to_datetime_with_apply(self):
        # this is only locale tested with US/None locales
        _skip_if_has_locale()

        # GH 5195
        # with a format and coerce a single item to_datetime fails
        td = Series(['May 04', 'Jun 02', 'Dec 11'], index=[1, 2, 3])
        expected = pd.to_datetime(td, format='%b %y')
        result = td.apply(pd.to_datetime, format='%b %y')
        assert_series_equal(result, expected)

        td = pd.Series(['May 04', 'Jun 02', ''], index=[1, 2, 3])
        self.assertRaises(ValueError,
                          lambda: pd.to_datetime(td, format='%b %y',
                                                 errors='raise'))
        self.assertRaises(ValueError,
                          lambda: td.apply(pd.to_datetime, format='%b %y',
                                           errors='raise'))
        expected = pd.to_datetime(td, format='%b %y', errors='coerce')

        result = td.apply(
            lambda x: pd.to_datetime(x, format='%b %y', errors='coerce'))
        assert_series_equal(result, expected)

    def test_nat_vector_field_access(self):
        idx = DatetimeIndex(['1/1/2000', None, None, '1/4/2000'])

        fields = ['year', 'quarter', 'month', 'day', 'hour', 'minute',
                  'second', 'microsecond', 'nanosecond', 'week', 'dayofyear',
                  'days_in_month', 'is_leap_year']

        for field in fields:
            result = getattr(idx, field)
            expected = [getattr(x, field) for x in idx]
            self.assert_numpy_array_equal(result, np.array(expected))

        s = pd.Series(idx)

        for field in fields:
            result = getattr(s.dt, field)
            expected = [getattr(x, field) for x in idx]
            self.assert_series_equal(result, pd.Series(expected))

    def test_nat_scalar_field_access(self):
        fields = ['year', 'quarter', 'month', 'day', 'hour', 'minute',
                  'second', 'microsecond', 'nanosecond', 'week', 'dayofyear',
                  'days_in_month', 'daysinmonth', 'dayofweek', 'weekday_name']
        for field in fields:
            result = getattr(NaT, field)
            self.assertTrue(np.isnan(result))

    def test_NaT_methods(self):
        # GH 9513
        raise_methods = ['astimezone', 'combine', 'ctime', 'dst',
                         'fromordinal', 'fromtimestamp', 'isocalendar',
                         'strftime', 'strptime', 'time', 'timestamp',
                         'timetuple', 'timetz', 'toordinal', 'tzname',
                         'utcfromtimestamp', 'utcnow', 'utcoffset',
                         'utctimetuple']
        nat_methods = ['date', 'now', 'replace', 'to_datetime', 'today']
        nan_methods = ['weekday', 'isoweekday']

        for method in raise_methods:
            if hasattr(NaT, method):
                self.assertRaises(ValueError, getattr(NaT, method))

        for method in nan_methods:
            if hasattr(NaT, method):
                self.assertTrue(np.isnan(getattr(NaT, method)()))

        for method in nat_methods:
            if hasattr(NaT, method):
                # see gh-8254
                exp_warning = None
                if method == 'to_datetime':
                    exp_warning = FutureWarning
                with tm.assert_produces_warning(
                        exp_warning, check_stacklevel=False):
                    self.assertIs(getattr(NaT, method)(), NaT)

        # GH 12300
        self.assertEqual(NaT.isoformat(), 'NaT')

    def test_to_datetime_types(self):

        # empty string
        result = to_datetime('')
        self.assertIs(result, NaT)

        result = to_datetime(['', ''])
        self.assertTrue(isnull(result).all())

        # ints
        result = Timestamp(0)
        expected = to_datetime(0)
        self.assertEqual(result, expected)

        # GH 3888 (strings)
        expected = to_datetime(['2012'])[0]
        result = to_datetime('2012')
        self.assertEqual(result, expected)

        # array = ['2012','20120101','20120101 12:01:01']
        array = ['20120101', '20120101 12:01:01']
        expected = list(to_datetime(array))
        result = lmap(Timestamp, array)
        tm.assert_almost_equal(result, expected)

        # currently fails ###
        # result = Timestamp('2012')
        # expected = to_datetime('2012')
        # self.assertEqual(result, expected)

    def test_to_datetime_unprocessable_input(self):
        # GH 4928
        self.assert_numpy_array_equal(
            to_datetime([1, '1'], errors='ignore'),
            np.array([1, '1'], dtype='O')
        )
        self.assertRaises(TypeError, to_datetime, [1, '1'], errors='raise')

    def test_to_datetime_other_datetime64_units(self):
        # 5/25/2012
        scalar = np.int64(1337904000000000).view('M8[us]')
        as_obj = scalar.astype('O')

        index = DatetimeIndex([scalar])
        self.assertEqual(index[0], scalar.astype('O'))

        value = Timestamp(scalar)
        self.assertEqual(value, as_obj)

    def test_to_datetime_list_of_integers(self):
        rng = date_range('1/1/2000', periods=20)
        rng = DatetimeIndex(rng.values)

        ints = list(rng.asi8)

        result = DatetimeIndex(ints)

        tm.assert_index_equal(rng, result)

    def test_to_datetime_freq(self):
        xp = bdate_range('2000-1-1', periods=10, tz='UTC')
        rs = xp.to_datetime()
        self.assertEqual(xp.freq, rs.freq)
        self.assertEqual(xp.tzinfo, rs.tzinfo)

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
            self.assertIn('2000', str(e))

    def test_reindex_with_datetimes(self):
        rng = date_range('1/1/2000', periods=20)
        ts = Series(np.random.randn(20), index=rng)

        result = ts.reindex(list(ts.index[5:10]))
        expected = ts[5:10]
        tm.assert_series_equal(result, expected)

        result = ts[list(ts.index[5:10])]
        tm.assert_series_equal(result, expected)

    def test_asfreq_keep_index_name(self):
        # GH #9854
        index_name = 'bar'
        index = pd.date_range('20130101', periods=20, name=index_name)
        df = pd.DataFrame([x for x in range(20)], columns=['foo'], index=index)

        self.assertEqual(index_name, df.index.name)
        self.assertEqual(index_name, df.asfreq('10D').index.name)

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
        self.assert_numpy_array_equal(result, expected)

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
        self.assertEqual(len(rng), 4)

    def test_date_range_negative_freq(self):
        # GH 11018
        rng = date_range('2011-12-31', freq='-2A', periods=3)
        exp = pd.DatetimeIndex(['2011-12-31', '2009-12-31',
                                '2007-12-31'], freq='-2A')
        tm.assert_index_equal(rng, exp)
        self.assertEqual(rng.freq, '-2A')

        rng = date_range('2011-01-31', freq='-2M', periods=3)
        exp = pd.DatetimeIndex(['2011-01-31', '2010-11-30',
                                '2010-09-30'], freq='-2M')
        tm.assert_index_equal(rng, exp)
        self.assertEqual(rng.freq, '-2M')

    def test_date_range_bms_bug(self):
        # #1645
        rng = date_range('1/1/2000', periods=10, freq='BMS')

        ex_first = Timestamp('2000-01-03')
        self.assertEqual(rng[0], ex_first)

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

    def test_first_subset(self):
        ts = _simple_ts('1/1/2000', '1/1/2010', freq='12h')
        result = ts.first('10d')
        self.assertEqual(len(result), 20)

        ts = _simple_ts('1/1/2000', '1/1/2010')
        result = ts.first('10d')
        self.assertEqual(len(result), 10)

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
        self.assertEqual(len(result), 20)

        ts = _simple_ts('1/1/2000', '1/1/2010')
        result = ts.last('10d')
        self.assertEqual(len(result), 10)

        result = ts.last('21D')
        expected = ts['12/12/2009':]
        assert_series_equal(result, expected)

        result = ts.last('21D')
        expected = ts[-21:]
        assert_series_equal(result, expected)

        result = ts[:0].last('3M')
        assert_series_equal(result, ts[:0])

    def test_format_pre_1900_dates(self):
        rng = date_range('1/1/1850', '1/1/1950', freq='A-DEC')
        rng.format()
        ts = Series(1, index=rng)
        repr(ts)

    def test_at_time(self):
        rng = date_range('1/1/2000', '1/5/2000', freq='5min')
        ts = Series(np.random.randn(len(rng)), index=rng)
        rs = ts.at_time(rng[1])
        self.assertTrue((rs.index.hour == rng[1].hour).all())
        self.assertTrue((rs.index.minute == rng[1].minute).all())
        self.assertTrue((rs.index.second == rng[1].second).all())

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
        self.assertEqual(len(rs), 0)

    def test_at_time_frame(self):
        rng = date_range('1/1/2000', '1/5/2000', freq='5min')
        ts = DataFrame(np.random.randn(len(rng), 2), index=rng)
        rs = ts.at_time(rng[1])
        self.assertTrue((rs.index.hour == rng[1].hour).all())
        self.assertTrue((rs.index.minute == rng[1].minute).all())
        self.assertTrue((rs.index.second == rng[1].second).all())

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
        self.assertEqual(len(rs), 0)

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

            self.assertEqual(len(filtered), exp_len)
            for rs in filtered.index:
                t = rs.time()
                if inc_start:
                    self.assertTrue(t >= stime)
                else:
                    self.assertTrue(t > stime)

                if inc_end:
                    self.assertTrue(t <= etime)
                else:
                    self.assertTrue(t < etime)

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

            self.assertEqual(len(filtered), exp_len)
            for rs in filtered.index:
                t = rs.time()
                if inc_start:
                    self.assertTrue((t >= stime) or (t <= etime))
                else:
                    self.assertTrue((t > stime) or (t <= etime))

                if inc_end:
                    self.assertTrue((t <= etime) or (t >= stime))
                else:
                    self.assertTrue((t < etime) or (t >= stime))

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

            self.assertEqual(len(filtered), exp_len)
            for rs in filtered.index:
                t = rs.time()
                if inc_start:
                    self.assertTrue(t >= stime)
                else:
                    self.assertTrue(t > stime)

                if inc_end:
                    self.assertTrue(t <= etime)
                else:
                    self.assertTrue(t < etime)

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

            self.assertEqual(len(filtered), exp_len)
            for rs in filtered.index:
                t = rs.time()
                if inc_start:
                    self.assertTrue((t >= stime) or (t <= etime))
                else:
                    self.assertTrue((t > stime) or (t <= etime))

                if inc_end:
                    self.assertTrue((t <= etime) or (t >= stime))
                else:
                    self.assertTrue((t < etime) or (t >= stime))

    def test_between_time_types(self):
        # GH11818
        rng = date_range('1/1/2000', '1/5/2000', freq='5min')
        self.assertRaises(ValueError, rng.indexer_between_time,
                          datetime(2010, 1, 2, 1), datetime(2010, 1, 2, 5))

        frame = DataFrame({'A': 0}, index=rng)
        self.assertRaises(ValueError, frame.between_time,
                          datetime(2010, 1, 2, 1), datetime(2010, 1, 2, 5))

        series = Series(0, index=rng)
        self.assertRaises(ValueError, series.between_time,
                          datetime(2010, 1, 2, 1), datetime(2010, 1, 2, 5))

    def test_between_time_formats(self):
        # GH11818
        _skip_if_has_locale()

        rng = date_range('1/1/2000', '1/5/2000', freq='5min')
        ts = DataFrame(np.random.randn(len(rng), 2), index=rng)

        strings = [("2:00", "2:30"), ("0200", "0230"), ("2:00am", "2:30am"),
                   ("0200am", "0230am"), ("2:00:00", "2:30:00"),
                   ("020000", "023000"), ("2:00:00am", "2:30:00am"),
                   ("020000am", "023000am")]
        expected_length = 28

        for time_string in strings:
            self.assertEqual(len(ts.between_time(*time_string)),
                             expected_length,
                             "%s - %s" % time_string)

    def test_dti_constructor_preserve_dti_freq(self):
        rng = date_range('1/1/2000', '1/2/2000', freq='5min')

        rng2 = DatetimeIndex(rng)
        self.assertEqual(rng.freq, rng2.freq)

    def test_dti_constructor_years_only(self):
        # GH 6961
        for tz in [None, 'UTC', 'Asia/Tokyo', 'dateutil/US/Pacific']:
            rng1 = date_range('2014', '2015', freq='M', tz=tz)
            expected1 = date_range('2014-01-31', '2014-12-31', freq='M', tz=tz)

            rng2 = date_range('2014', '2015', freq='MS', tz=tz)
            expected2 = date_range('2014-01-01', '2015-01-01', freq='MS',
                                   tz=tz)

            rng3 = date_range('2014', '2020', freq='A', tz=tz)
            expected3 = date_range('2014-12-31', '2019-12-31', freq='A', tz=tz)

            rng4 = date_range('2014', '2020', freq='AS', tz=tz)
            expected4 = date_range('2014-01-01', '2020-01-01', freq='AS',
                                   tz=tz)

            for rng, expected in [(rng1, expected1), (rng2, expected2),
                                  (rng3, expected3), (rng4, expected4)]:
                tm.assert_index_equal(rng, expected)

    def test_dti_constructor_small_int(self):
        # GH 13721
        exp = DatetimeIndex(['1970-01-01 00:00:00.00000000',
                             '1970-01-01 00:00:00.00000001',
                             '1970-01-01 00:00:00.00000002'])

        for dtype in [np.int64, np.int32, np.int16, np.int8]:
            arr = np.array([0, 10, 20], dtype=dtype)
            tm.assert_index_equal(DatetimeIndex(arr), exp)

    def test_dti_constructor_numpy_timeunits(self):
        # GH 9114
        base = pd.to_datetime(['2000-01-01T00:00', '2000-01-02T00:00', 'NaT'])

        for dtype in ['datetime64[h]', 'datetime64[m]', 'datetime64[s]',
                      'datetime64[ms]', 'datetime64[us]', 'datetime64[ns]']:
            values = base.values.astype(dtype)

            tm.assert_index_equal(DatetimeIndex(values), base)
            tm.assert_index_equal(to_datetime(values), base)

    def test_normalize(self):
        rng = date_range('1/1/2000 9:30', periods=10, freq='D')

        result = rng.normalize()
        expected = date_range('1/1/2000', periods=10, freq='D')
        tm.assert_index_equal(result, expected)

        rng_ns = pd.DatetimeIndex(np.array([1380585623454345752,
                                            1380585612343234312]).astype(
                                                "datetime64[ns]"))
        rng_ns_normalized = rng_ns.normalize()
        expected = pd.DatetimeIndex(np.array([1380585600000000000,
                                              1380585600000000000]).astype(
                                                  "datetime64[ns]"))
        tm.assert_index_equal(rng_ns_normalized, expected)

        self.assertTrue(result.is_normalized)
        self.assertFalse(rng.is_normalized)

    def test_to_period(self):
        from pandas.tseries.period import period_range

        ts = _simple_ts('1/1/2000', '1/1/2001')

        pts = ts.to_period()
        exp = ts.copy()
        exp.index = period_range('1/1/2000', '1/1/2001')
        assert_series_equal(pts, exp)

        pts = ts.to_period('M')
        exp.index = exp.index.asfreq('M')
        tm.assert_index_equal(pts.index, exp.index.asfreq('M'))
        assert_series_equal(pts, exp)

        # GH 7606 without freq
        idx = DatetimeIndex(['2011-01-01', '2011-01-02', '2011-01-03',
                             '2011-01-04'])
        exp_idx = pd.PeriodIndex(['2011-01-01', '2011-01-02', '2011-01-03',
                                  '2011-01-04'], freq='D')

        s = Series(np.random.randn(4), index=idx)
        expected = s.copy()
        expected.index = exp_idx
        assert_series_equal(s.to_period(), expected)

        df = DataFrame(np.random.randn(4, 4), index=idx, columns=idx)
        expected = df.copy()
        expected.index = exp_idx
        assert_frame_equal(df.to_period(), expected)

        expected = df.copy()
        expected.columns = exp_idx
        assert_frame_equal(df.to_period(axis=1), expected)

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

    def test_to_period_tz_pytz(self):
        tm._skip_if_no_pytz()
        from dateutil.tz import tzlocal
        from pytz import utc as UTC

        xp = date_range('1/1/2000', '4/1/2000').to_period()

        ts = date_range('1/1/2000', '4/1/2000', tz='US/Eastern')

        result = ts.to_period()[0]
        expected = ts[0].to_period()

        self.assertEqual(result, expected)
        tm.assert_index_equal(ts.to_period(), xp)

        ts = date_range('1/1/2000', '4/1/2000', tz=UTC)

        result = ts.to_period()[0]
        expected = ts[0].to_period()

        self.assertEqual(result, expected)
        tm.assert_index_equal(ts.to_period(), xp)

        ts = date_range('1/1/2000', '4/1/2000', tz=tzlocal())

        result = ts.to_period()[0]
        expected = ts[0].to_period()

        self.assertEqual(result, expected)
        tm.assert_index_equal(ts.to_period(), xp)

    def test_to_period_tz_explicit_pytz(self):
        tm._skip_if_no_pytz()
        import pytz
        from dateutil.tz import tzlocal

        xp = date_range('1/1/2000', '4/1/2000').to_period()

        ts = date_range('1/1/2000', '4/1/2000', tz=pytz.timezone('US/Eastern'))

        result = ts.to_period()[0]
        expected = ts[0].to_period()

        self.assertTrue(result == expected)
        tm.assert_index_equal(ts.to_period(), xp)

        ts = date_range('1/1/2000', '4/1/2000', tz=pytz.utc)

        result = ts.to_period()[0]
        expected = ts[0].to_period()

        self.assertTrue(result == expected)
        tm.assert_index_equal(ts.to_period(), xp)

        ts = date_range('1/1/2000', '4/1/2000', tz=tzlocal())

        result = ts.to_period()[0]
        expected = ts[0].to_period()

        self.assertTrue(result == expected)
        tm.assert_index_equal(ts.to_period(), xp)

    def test_to_period_tz_dateutil(self):
        tm._skip_if_no_dateutil()
        import dateutil
        from dateutil.tz import tzlocal

        xp = date_range('1/1/2000', '4/1/2000').to_period()

        ts = date_range('1/1/2000', '4/1/2000', tz='dateutil/US/Eastern')

        result = ts.to_period()[0]
        expected = ts[0].to_period()

        self.assertTrue(result == expected)
        tm.assert_index_equal(ts.to_period(), xp)

        ts = date_range('1/1/2000', '4/1/2000', tz=dateutil.tz.tzutc())

        result = ts.to_period()[0]
        expected = ts[0].to_period()

        self.assertTrue(result == expected)
        tm.assert_index_equal(ts.to_period(), xp)

        ts = date_range('1/1/2000', '4/1/2000', tz=tzlocal())

        result = ts.to_period()[0]
        expected = ts[0].to_period()

        self.assertTrue(result == expected)
        tm.assert_index_equal(ts.to_period(), xp)

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
        tm.assert_index_equal(pts.index, exp.index.asfreq('M'))

        df = df.T
        pts = df.to_period(axis=1)
        exp = df.copy()
        exp.columns = pr
        assert_frame_equal(pts, exp)

        pts = df.to_period('M', axis=1)
        tm.assert_index_equal(pts.columns, exp.columns.asfreq('M'))

        self.assertRaises(ValueError, df.to_period, axis=2)

    def test_timestamp_fields(self):
        # extra fields from DatetimeIndex like quarter and week
        idx = tm.makeDateIndex(100)

        fields = ['dayofweek', 'dayofyear', 'week', 'weekofyear', 'quarter',
                  'days_in_month', 'is_month_start', 'is_month_end',
                  'is_quarter_start', 'is_quarter_end', 'is_year_start',
                  'is_year_end', 'weekday_name']
        for f in fields:
            expected = getattr(idx, f)[-1]
            result = getattr(Timestamp(idx[-1]), f)
            self.assertEqual(result, expected)

        self.assertEqual(idx.freq, Timestamp(idx[-1], idx.freq).freq)
        self.assertEqual(idx.freqstr, Timestamp(idx[-1], idx.freq).freqstr)

    def test_woy_boundary(self):
        # make sure weeks at year boundaries are correct
        d = datetime(2013, 12, 31)
        result = Timestamp(d).week
        expected = 1  # ISO standard
        self.assertEqual(result, expected)

        d = datetime(2008, 12, 28)
        result = Timestamp(d).week
        expected = 52  # ISO standard
        self.assertEqual(result, expected)

        d = datetime(2009, 12, 31)
        result = Timestamp(d).week
        expected = 53  # ISO standard
        self.assertEqual(result, expected)

        d = datetime(2010, 1, 1)
        result = Timestamp(d).week
        expected = 53  # ISO standard
        self.assertEqual(result, expected)

        d = datetime(2010, 1, 3)
        result = Timestamp(d).week
        expected = 53  # ISO standard
        self.assertEqual(result, expected)

        result = np.array([Timestamp(datetime(*args)).week
                           for args in [(2000, 1, 1), (2000, 1, 2), (
                               2005, 1, 1), (2005, 1, 2)]])
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
        self.assertIn(iso8601, result)

    def test_timestamp_from_ordinal(self):

        # GH 3042
        dt = datetime(2011, 4, 16, 0, 0)
        ts = Timestamp.fromordinal(dt.toordinal())
        self.assertEqual(ts.to_pydatetime(), dt)

        # with a tzinfo
        stamp = Timestamp('2011-4-16', tz='US/Eastern')
        dt_tz = stamp.to_pydatetime()
        ts = Timestamp.fromordinal(dt_tz.toordinal(), tz='US/Eastern')
        self.assertEqual(ts.to_pydatetime(), dt_tz)

    def test_datetimeindex_integers_shift(self):
        rng = date_range('1/1/2000', periods=20)

        result = rng + 5
        expected = rng.shift(5)
        tm.assert_index_equal(result, expected)

        result = rng - 5
        expected = rng.shift(-5)
        tm.assert_index_equal(result, expected)

    def test_astype_object(self):
        # NumPy 1.6.1 weak ns support
        rng = date_range('1/1/2000', periods=20)

        casted = rng.astype('O')
        exp_values = list(rng)

        tm.assert_index_equal(casted, Index(exp_values, dtype=np.object_))
        self.assertEqual(casted.tolist(), exp_values)

    def test_catch_infinite_loop(self):
        offset = offsets.DateOffset(minute=5)
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
        tm.assert_index_equal(result.index, ex_index)
        tm.assert_index_equal(result_df.index, ex_index)

        appended = rng.append(rng)
        tm.assert_index_equal(appended, ex_index)

        appended = rng.append([rng, rng])
        ex_index = DatetimeIndex(np.tile(rng.values, 3))
        tm.assert_index_equal(appended, ex_index)

        # different index names
        rng1 = rng.copy()
        rng2 = rng.copy()
        rng1.name = 'foo'
        rng2.name = 'bar'
        self.assertEqual(rng1.append(rng1).name, 'foo')
        self.assertIsNone(rng1.append(rng2).name)

    def test_append_concat_tz(self):
        # GH 2938
        tm._skip_if_no_pytz()

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
        tm.assert_index_equal(result.index, rng3)
        tm.assert_index_equal(result_df.index, rng3)

        appended = rng.append(rng2)
        tm.assert_index_equal(appended, rng3)

    def test_append_concat_tz_explicit_pytz(self):
        # GH 2938
        tm._skip_if_no_pytz()
        from pytz import timezone as timezone

        rng = date_range('5/8/2012 1:45', periods=10, freq='5T',
                         tz=timezone('US/Eastern'))
        rng2 = date_range('5/8/2012 2:35', periods=10, freq='5T',
                          tz=timezone('US/Eastern'))
        rng3 = date_range('5/8/2012 1:45', periods=20, freq='5T',
                          tz=timezone('US/Eastern'))
        ts = Series(np.random.randn(len(rng)), rng)
        df = DataFrame(np.random.randn(len(rng), 4), index=rng)
        ts2 = Series(np.random.randn(len(rng2)), rng2)
        df2 = DataFrame(np.random.randn(len(rng2), 4), index=rng2)

        result = ts.append(ts2)
        result_df = df.append(df2)
        tm.assert_index_equal(result.index, rng3)
        tm.assert_index_equal(result_df.index, rng3)

        appended = rng.append(rng2)
        tm.assert_index_equal(appended, rng3)

    def test_append_concat_tz_dateutil(self):
        # GH 2938
        tm._skip_if_no_dateutil()
        rng = date_range('5/8/2012 1:45', periods=10, freq='5T',
                         tz='dateutil/US/Eastern')
        rng2 = date_range('5/8/2012 2:35', periods=10, freq='5T',
                          tz='dateutil/US/Eastern')
        rng3 = date_range('5/8/2012 1:45', periods=20, freq='5T',
                          tz='dateutil/US/Eastern')
        ts = Series(np.random.randn(len(rng)), rng)
        df = DataFrame(np.random.randn(len(rng), 4), index=rng)
        ts2 = Series(np.random.randn(len(rng2)), rng2)
        df2 = DataFrame(np.random.randn(len(rng2), 4), index=rng2)

        result = ts.append(ts2)
        result_df = df.append(df2)
        tm.assert_index_equal(result.index, rng3)
        tm.assert_index_equal(result_df.index, rng3)

        appended = rng.append(rng2)
        tm.assert_index_equal(appended, rng3)

    def test_set_dataframe_column_ns_dtype(self):
        x = DataFrame([datetime.now(), datetime.now()])
        self.assertEqual(x[0].dtype, np.dtype('M8[ns]'))

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
        self.assertTrue((index.asi8[50:100] == -1).all())

        arr = np.arange(1000, dtype=np.int64)
        index = DatetimeIndex(arr, copy=True)

        arr[50:100] = -1
        self.assertTrue((index.asi8[50:100] != -1).all())

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
        new_index = index.append(index + pd.DateOffset(days=1)).sort_values()

        exp = ts.reindex(new_index).interpolate(method='time')

        index = pd.date_range('1/1/2012', periods=4, freq='12H')
        ts = pd.Series([0, 12, 24, 36], index)
        new_index = index.append(index + pd.DateOffset(hours=1)).sort_values()
        result = ts.reindex(new_index).interpolate(method='time')

        self.assert_numpy_array_equal(result.values, exp.values)

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
        self.assertTrue(d['B'].isnull().all())

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
        self.assertTrue((-result).all())

        tst = DataFrame({'date': dates})
        result = tst.duplicated()
        self.assertTrue((-result).all())

    def test_timestamp_compare_with_early_datetime(self):
        # e.g. datetime.min
        stamp = Timestamp('2012-01-01')

        self.assertFalse(stamp == datetime.min)
        self.assertFalse(stamp == datetime(1600, 1, 1))
        self.assertFalse(stamp == datetime(2700, 1, 1))
        self.assertNotEqual(stamp, datetime.min)
        self.assertNotEqual(stamp, datetime(1600, 1, 1))
        self.assertNotEqual(stamp, datetime(2700, 1, 1))
        self.assertTrue(stamp > datetime(1600, 1, 1))
        self.assertTrue(stamp >= datetime(1600, 1, 1))
        self.assertTrue(stamp < datetime(2700, 1, 1))
        self.assertTrue(stamp <= datetime(2700, 1, 1))

    def test_to_html_timestamp(self):
        rng = date_range('2000-01-01', periods=10)
        df = DataFrame(np.random.randn(10, 4), index=rng)

        result = df.to_html()
        self.assertIn('2000-01-01', result)

    def test_to_csv_numpy_16_bug(self):
        frame = DataFrame({'a': date_range('1/1/2000', periods=10)})

        buf = StringIO()
        frame.to_csv(buf)

        result = buf.getvalue()
        self.assertIn('2000-01-01', result)

    def test_series_map_box_timestamps(self):
        # #2689, #2627
        s = Series(date_range('1/1/2000', periods=10))

        def f(x):
            return (x.hour, x.day, x.month)

        # it works!
        s.map(f)
        s.apply(f)
        DataFrame(s).applymap(f)

    def test_series_map_box_timedelta(self):
        # GH 11349
        s = Series(timedelta_range('1 day 1 s', periods=5, freq='h'))

        def f(x):
            return x.total_seconds()

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
        df1 = DataFrame({'date': ind, 'test': lrange(10)})

        # it works!
        pd.concat([df1, df2_obj])

    def test_asfreq_resample_set_correct_freq(self):
        # GH5613
        # we test if .asfreq() and .resample() set the correct value for .freq
        df = pd.DataFrame({'date': ["2012-01-01", "2012-01-02", "2012-01-03"],
                           'col': [1, 2, 3]})
        df = df.set_index(pd.to_datetime(df.date))

        # testing the settings before calling .asfreq() and .resample()
        self.assertEqual(df.index.freq, None)
        self.assertEqual(df.index.inferred_freq, 'D')

        # does .asfreq() set .freq correctly?
        self.assertEqual(df.asfreq('D').index.freq, 'D')

        # does .resample() set .freq correctly?
        self.assertEqual(df.resample('D').asfreq().index.freq, 'D')

    def test_pickle(self):

        # GH4606
        p = self.round_trip_pickle(NaT)
        self.assertTrue(p is NaT)

        idx = pd.to_datetime(['2013-01-01', NaT, '2014-01-06'])
        idx_p = self.round_trip_pickle(idx)
        self.assertTrue(idx_p[0] == idx[0])
        self.assertTrue(idx_p[1] is NaT)
        self.assertTrue(idx_p[2] == idx[2])

        # GH11002
        # don't infer freq
        idx = date_range('1750-1-1', '2050-1-1', freq='7D')
        idx_p = self.round_trip_pickle(idx)
        tm.assert_index_equal(idx, idx_p)

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


def _simple_ts(start, end, freq='D'):
    rng = date_range(start, end, freq=freq)
    return Series(np.random.randn(len(rng)), index=rng)


class TestToDatetime(tm.TestCase):
    _multiprocess_can_split_ = True

    def test_to_datetime_dt64s(self):
        in_bound_dts = [
            np.datetime64('2000-01-01'),
            np.datetime64('2000-01-02'),
        ]

        for dt in in_bound_dts:
            self.assertEqual(pd.to_datetime(dt), Timestamp(dt))

        oob_dts = [np.datetime64('1000-01-01'), np.datetime64('5000-01-02'), ]

        for dt in oob_dts:
            self.assertRaises(ValueError, pd.to_datetime, dt, errors='raise')
            self.assertRaises(ValueError, tslib.Timestamp, dt)
            self.assertIs(pd.to_datetime(dt, errors='coerce'), NaT)

    def test_to_datetime_array_of_dt64s(self):
        dts = [np.datetime64('2000-01-01'), np.datetime64('2000-01-02'), ]

        # Assuming all datetimes are in bounds, to_datetime() returns
        # an array that is equal to Timestamp() parsing
        self.assert_numpy_array_equal(
            pd.to_datetime(dts, box=False),
            np.array([Timestamp(x).asm8 for x in dts])
        )

        # A list of datetimes where the last one is out of bounds
        dts_with_oob = dts + [np.datetime64('9999-01-01')]

        self.assertRaises(ValueError, pd.to_datetime, dts_with_oob,
                          errors='raise')

        self.assert_numpy_array_equal(
            pd.to_datetime(dts_with_oob, box=False, errors='coerce'),
            np.array(
                [
                    Timestamp(dts_with_oob[0]).asm8,
                    Timestamp(dts_with_oob[1]).asm8,
                    iNaT,
                ],
                dtype='M8'
            )
        )

        # With errors='ignore', out of bounds datetime64s
        # are converted to their .item(), which depending on the version of
        # numpy is either a python datetime.datetime or datetime.date
        self.assert_numpy_array_equal(
            pd.to_datetime(dts_with_oob, box=False, errors='ignore'),
            np.array(
                [dt.item() for dt in dts_with_oob],
                dtype='O'
            )
        )

    def test_to_datetime_tz(self):

        # xref 8260
        # uniform returns a DatetimeIndex
        arr = [pd.Timestamp('2013-01-01 13:00:00-0800', tz='US/Pacific'),
               pd.Timestamp('2013-01-02 14:00:00-0800', tz='US/Pacific')]
        result = pd.to_datetime(arr)
        expected = DatetimeIndex(
            ['2013-01-01 13:00:00', '2013-01-02 14:00:00'], tz='US/Pacific')
        tm.assert_index_equal(result, expected)

        # mixed tzs will raise
        arr = [pd.Timestamp('2013-01-01 13:00:00', tz='US/Pacific'),
               pd.Timestamp('2013-01-02 14:00:00', tz='US/Eastern')]
        self.assertRaises(ValueError, lambda: pd.to_datetime(arr))

    def test_to_datetime_tz_pytz(self):

        # xref 8260
        tm._skip_if_no_pytz()
        import pytz

        us_eastern = pytz.timezone('US/Eastern')
        arr = np.array([us_eastern.localize(datetime(year=2000, month=1, day=1,
                                                     hour=3, minute=0)),
                        us_eastern.localize(datetime(year=2000, month=6, day=1,
                                                     hour=3, minute=0))],
                       dtype=object)
        result = pd.to_datetime(arr, utc=True)
        expected = DatetimeIndex(['2000-01-01 08:00:00+00:00',
                                  '2000-06-01 07:00:00+00:00'],
                                 dtype='datetime64[ns, UTC]', freq=None)
        tm.assert_index_equal(result, expected)

    def test_to_datetime_utc_is_true(self):
        # See gh-11934
        start = pd.Timestamp('2014-01-01', tz='utc')
        end = pd.Timestamp('2014-01-03', tz='utc')
        date_range = pd.bdate_range(start, end)

        result = pd.to_datetime(date_range, utc=True)
        expected = pd.DatetimeIndex(data=date_range)
        tm.assert_index_equal(result, expected)

    def test_to_datetime_tz_psycopg2(self):

        # xref 8260
        try:
            import psycopg2
        except ImportError:
            raise nose.SkipTest("no psycopg2 installed")

        # misc cases
        tz1 = psycopg2.tz.FixedOffsetTimezone(offset=-300, name=None)
        tz2 = psycopg2.tz.FixedOffsetTimezone(offset=-240, name=None)
        arr = np.array([datetime(2000, 1, 1, 3, 0, tzinfo=tz1),
                        datetime(2000, 6, 1, 3, 0, tzinfo=tz2)],
                       dtype=object)

        result = pd.to_datetime(arr, errors='coerce', utc=True)
        expected = DatetimeIndex(['2000-01-01 08:00:00+00:00',
                                  '2000-06-01 07:00:00+00:00'],
                                 dtype='datetime64[ns, UTC]', freq=None)
        tm.assert_index_equal(result, expected)

        # dtype coercion
        i = pd.DatetimeIndex([
            '2000-01-01 08:00:00+00:00'
        ], tz=psycopg2.tz.FixedOffsetTimezone(offset=-300, name=None))
        self.assertFalse(is_datetime64_ns_dtype(i))

        # tz coerceion
        result = pd.to_datetime(i, errors='coerce')
        tm.assert_index_equal(result, i)

        result = pd.to_datetime(i, errors='coerce', utc=True)
        expected = pd.DatetimeIndex(['2000-01-01 13:00:00'],
                                    dtype='datetime64[ns, UTC]')
        tm.assert_index_equal(result, expected)

    def test_datetime_bool(self):
        # GH13176
        with self.assertRaises(TypeError):
            to_datetime(False)
        self.assertTrue(to_datetime(False, errors="coerce") is tslib.NaT)
        self.assertEqual(to_datetime(False, errors="ignore"), False)
        with self.assertRaises(TypeError):
            to_datetime(True)
        self.assertTrue(to_datetime(True, errors="coerce") is tslib.NaT)
        self.assertEqual(to_datetime(True, errors="ignore"), True)
        with self.assertRaises(TypeError):
            to_datetime([False, datetime.today()])
        with self.assertRaises(TypeError):
            to_datetime(['20130101', True])
        tm.assert_index_equal(to_datetime([0, False, tslib.NaT, 0.0],
                                          errors="coerce"),
                              DatetimeIndex([to_datetime(0), tslib.NaT,
                                             tslib.NaT, to_datetime(0)]))

    def test_datetime_invalid_datatype(self):
        # GH13176

        with self.assertRaises(TypeError):
            pd.to_datetime(bool)
        with self.assertRaises(TypeError):
            pd.to_datetime(pd.to_datetime)

    def test_unit(self):
        # GH 11758
        # test proper behavior with erros

        with self.assertRaises(ValueError):
            to_datetime([1], unit='D', format='%Y%m%d')

        values = [11111111, 1, 1.0, tslib.iNaT, pd.NaT, np.nan,
                  'NaT', '']
        result = to_datetime(values, unit='D', errors='ignore')
        expected = Index([11111111, Timestamp('1970-01-02'),
                          Timestamp('1970-01-02'), pd.NaT,
                          pd.NaT, pd.NaT, pd.NaT, pd.NaT],
                         dtype=object)
        tm.assert_index_equal(result, expected)

        result = to_datetime(values, unit='D', errors='coerce')
        expected = DatetimeIndex(['NaT', '1970-01-02', '1970-01-02',
                                  'NaT', 'NaT', 'NaT', 'NaT', 'NaT'])
        tm.assert_index_equal(result, expected)

        with self.assertRaises(tslib.OutOfBoundsDatetime):
            to_datetime(values, unit='D', errors='raise')

        values = [1420043460000, tslib.iNaT, pd.NaT, np.nan, 'NaT']

        result = to_datetime(values, errors='ignore', unit='s')
        expected = Index([1420043460000, pd.NaT, pd.NaT,
                          pd.NaT, pd.NaT], dtype=object)
        tm.assert_index_equal(result, expected)

        result = to_datetime(values, errors='coerce', unit='s')
        expected = DatetimeIndex(['NaT', 'NaT', 'NaT', 'NaT', 'NaT'])
        tm.assert_index_equal(result, expected)

        with self.assertRaises(tslib.OutOfBoundsDatetime):
            to_datetime(values, errors='raise', unit='s')

        # if we have a string, then we raise a ValueError
        # and NOT an OutOfBoundsDatetime
        for val in ['foo', Timestamp('20130101')]:
            try:
                to_datetime(val, errors='raise', unit='s')
            except tslib.OutOfBoundsDatetime:
                raise AssertionError("incorrect exception raised")
            except ValueError:
                pass

    def test_unit_consistency(self):

        # consistency of conversions
        expected = Timestamp('1970-05-09 14:25:11')
        result = pd.to_datetime(11111111, unit='s', errors='raise')
        self.assertEqual(result, expected)
        self.assertIsInstance(result, Timestamp)

        result = pd.to_datetime(11111111, unit='s', errors='coerce')
        self.assertEqual(result, expected)
        self.assertIsInstance(result, Timestamp)

        result = pd.to_datetime(11111111, unit='s', errors='ignore')
        self.assertEqual(result, expected)
        self.assertIsInstance(result, Timestamp)

    def test_unit_with_numeric(self):

        # GH 13180
        # coercions from floats/ints are ok
        expected = DatetimeIndex(['2015-06-19 05:33:20',
                                  '2015-05-27 22:33:20'])
        arr1 = [1.434692e+18, 1.432766e+18]
        arr2 = np.array(arr1).astype('int64')
        for errors in ['ignore', 'raise', 'coerce']:
            result = pd.to_datetime(arr1, errors=errors)
            tm.assert_index_equal(result, expected)

            result = pd.to_datetime(arr2, errors=errors)
            tm.assert_index_equal(result, expected)

        # but we want to make sure that we are coercing
        # if we have ints/strings
        expected = DatetimeIndex(['NaT',
                                  '2015-06-19 05:33:20',
                                  '2015-05-27 22:33:20'])
        arr = ['foo', 1.434692e+18, 1.432766e+18]
        result = pd.to_datetime(arr, errors='coerce')
        tm.assert_index_equal(result, expected)

        expected = DatetimeIndex(['2015-06-19 05:33:20',
                                  '2015-05-27 22:33:20',
                                  'NaT',
                                  'NaT'])
        arr = [1.434692e+18, 1.432766e+18, 'foo', 'NaT']
        result = pd.to_datetime(arr, errors='coerce')
        tm.assert_index_equal(result, expected)

    def test_unit_mixed(self):

        # mixed integers/datetimes
        expected = DatetimeIndex(['2013-01-01', 'NaT', 'NaT'])
        arr = [pd.Timestamp('20130101'), 1.434692e+18, 1.432766e+18]
        result = pd.to_datetime(arr, errors='coerce')
        tm.assert_index_equal(result, expected)

        with self.assertRaises(ValueError):
            pd.to_datetime(arr, errors='raise')

        expected = DatetimeIndex(['NaT',
                                  'NaT',
                                  '2013-01-01'])
        arr = [1.434692e+18, 1.432766e+18, pd.Timestamp('20130101')]
        result = pd.to_datetime(arr, errors='coerce')
        tm.assert_index_equal(result, expected)

        with self.assertRaises(ValueError):
            pd.to_datetime(arr, errors='raise')

    def test_index_to_datetime(self):
        idx = Index(['1/1/2000', '1/2/2000', '1/3/2000'])

        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            result = idx.to_datetime()
            expected = DatetimeIndex(pd.to_datetime(idx.values))
            tm.assert_index_equal(result, expected)

        with tm.assert_produces_warning(FutureWarning,
                                        check_stacklevel=False):
            today = datetime.today()
            idx = Index([today], dtype=object)
            result = idx.to_datetime()
            expected = DatetimeIndex([today])
            tm.assert_index_equal(result, expected)

    def test_dataframe(self):

        df = DataFrame({'year': [2015, 2016],
                        'month': [2, 3],
                        'day': [4, 5],
                        'hour': [6, 7],
                        'minute': [58, 59],
                        'second': [10, 11],
                        'ms': [1, 1],
                        'us': [2, 2],
                        'ns': [3, 3]})

        result = to_datetime({'year': df['year'],
                              'month': df['month'],
                              'day': df['day']})
        expected = Series([Timestamp('20150204 00:00:00'),
                           Timestamp('20160305 00:0:00')])
        assert_series_equal(result, expected)

        # dict-like
        result = to_datetime(df[['year', 'month', 'day']].to_dict())
        assert_series_equal(result, expected)

        # dict but with constructable
        df2 = df[['year', 'month', 'day']].to_dict()
        df2['month'] = 2
        result = to_datetime(df2)
        expected2 = Series([Timestamp('20150204 00:00:00'),
                            Timestamp('20160205 00:0:00')])
        assert_series_equal(result, expected2)

        # unit mappings
        units = [{'year': 'years',
                  'month': 'months',
                  'day': 'days',
                  'hour': 'hours',
                  'minute': 'minutes',
                  'second': 'seconds'},
                 {'year': 'year',
                  'month': 'month',
                  'day': 'day',
                  'hour': 'hour',
                  'minute': 'minute',
                  'second': 'second'},
                 ]

        for d in units:
            result = to_datetime(df[list(d.keys())].rename(columns=d))
            expected = Series([Timestamp('20150204 06:58:10'),
                               Timestamp('20160305 07:59:11')])
            assert_series_equal(result, expected)

        d = {'year': 'year',
             'month': 'month',
             'day': 'day',
             'hour': 'hour',
             'minute': 'minute',
             'second': 'second',
             'ms': 'ms',
             'us': 'us',
             'ns': 'ns'}

        result = to_datetime(df.rename(columns=d))
        expected = Series([Timestamp('20150204 06:58:10.001002003'),
                           Timestamp('20160305 07:59:11.001002003')])
        assert_series_equal(result, expected)

        # coerce back to int
        result = to_datetime(df.astype(str))
        assert_series_equal(result, expected)

        # passing coerce
        df2 = DataFrame({'year': [2015, 2016],
                         'month': [2, 20],
                         'day': [4, 5]})
        with self.assertRaises(ValueError):
            to_datetime(df2)
        result = to_datetime(df2, errors='coerce')
        expected = Series([Timestamp('20150204 00:00:00'),
                           pd.NaT])
        assert_series_equal(result, expected)

        # extra columns
        with self.assertRaises(ValueError):
            df2 = df.copy()
            df2['foo'] = 1
            to_datetime(df2)

        # not enough
        for c in [['year'],
                  ['year', 'month'],
                  ['year', 'month', 'second'],
                  ['month', 'day'],
                  ['year', 'day', 'second']]:
            with self.assertRaises(ValueError):
                to_datetime(df[c])

        # duplicates
        df2 = DataFrame({'year': [2015, 2016],
                         'month': [2, 20],
                         'day': [4, 5]})
        df2.columns = ['year', 'year', 'day']
        with self.assertRaises(ValueError):
            to_datetime(df2)

        df2 = DataFrame({'year': [2015, 2016],
                         'month': [2, 20],
                         'day': [4, 5],
                         'hour': [4, 5]})
        df2.columns = ['year', 'month', 'day', 'day']
        with self.assertRaises(ValueError):
            to_datetime(df2)

    def test_dataframe_dtypes(self):
        # #13451
        df = DataFrame({'year': [2015, 2016],
                        'month': [2, 3],
                        'day': [4, 5]})

        # int16
        result = to_datetime(df.astype('int16'))
        expected = Series([Timestamp('20150204 00:00:00'),
                           Timestamp('20160305 00:00:00')])
        assert_series_equal(result, expected)

        # mixed dtypes
        df['month'] = df['month'].astype('int8')
        df['day'] = df['day'].astype('int8')
        result = to_datetime(df)
        expected = Series([Timestamp('20150204 00:00:00'),
                           Timestamp('20160305 00:00:00')])
        assert_series_equal(result, expected)

        # float
        df = DataFrame({'year': [2000, 2001],
                        'month': [1.5, 1],
                        'day': [1, 1]})
        with self.assertRaises(ValueError):
            to_datetime(df)


class TestDatetimeIndex(tm.TestCase):
    _multiprocess_can_split_ = True

    def test_hash_error(self):
        index = date_range('20010101', periods=10)
        with tm.assertRaisesRegexp(TypeError, "unhashable type: %r" %
                                   type(index).__name__):
            hash(index)

    def test_stringified_slice_with_tz(self):
        # GH2658
        import datetime
        start = datetime.datetime.now()
        idx = DatetimeIndex(start=start, freq="1d", periods=10)
        df = DataFrame(lrange(10), index=idx)
        df["2013-01-14 23:44:34.437768-05:00":]  # no exception here

    def test_append_join_nondatetimeindex(self):
        rng = date_range('1/1/2000', periods=10)
        idx = Index(['a', 'b', 'c', 'd'])

        result = rng.append(idx)
        tm.assertIsInstance(result[0], Timestamp)

        # it works
        rng.join(idx, how='outer')

    def test_to_period_nofreq(self):
        idx = DatetimeIndex(['2000-01-01', '2000-01-02', '2000-01-04'])
        self.assertRaises(ValueError, idx.to_period)

        idx = DatetimeIndex(['2000-01-01', '2000-01-02', '2000-01-03'],
                            freq='infer')
        self.assertEqual(idx.freqstr, 'D')
        expected = pd.PeriodIndex(['2000-01-01', '2000-01-02',
                                   '2000-01-03'], freq='D')
        tm.assert_index_equal(idx.to_period(), expected)

        # GH 7606
        idx = DatetimeIndex(['2000-01-01', '2000-01-02', '2000-01-03'])
        self.assertEqual(idx.freqstr, None)
        tm.assert_index_equal(idx.to_period(), expected)

    def test_000constructor_resolution(self):
        # 2252
        t1 = Timestamp((1352934390 * 1000000000) + 1000000 + 1000 + 1)
        idx = DatetimeIndex([t1])

        self.assertEqual(idx.nanosecond[0], t1.nanosecond)

    def test_constructor_coverage(self):
        rng = date_range('1/1/2000', periods=10.5)
        exp = date_range('1/1/2000', periods=10)
        tm.assert_index_equal(rng, exp)

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
        tm.assert_index_equal(result, expected)

        # NumPy string array
        strings = np.array(['2000-01-01', '2000-01-02', '2000-01-03'])
        result = DatetimeIndex(strings)
        expected = DatetimeIndex(strings.astype('O'))
        tm.assert_index_equal(result, expected)

        from_ints = DatetimeIndex(expected.asi8)
        tm.assert_index_equal(from_ints, expected)

        # string with NaT
        strings = np.array(['2000-01-01', '2000-01-02', 'NaT'])
        result = DatetimeIndex(strings)
        expected = DatetimeIndex(strings.astype('O'))
        tm.assert_index_equal(result, expected)

        from_ints = DatetimeIndex(expected.asi8)
        tm.assert_index_equal(from_ints, expected)

        # non-conforming
        self.assertRaises(ValueError, DatetimeIndex,
                          ['2000-01-01', '2000-01-02', '2000-01-04'], freq='D')

        self.assertRaises(ValueError, DatetimeIndex, start='2011-01-01',
                          freq='b')
        self.assertRaises(ValueError, DatetimeIndex, end='2011-01-01',
                          freq='B')
        self.assertRaises(ValueError, DatetimeIndex, periods=10, freq='D')

    def test_constructor_datetime64_tzformat(self):
        # GH 6572
        tm._skip_if_no_pytz()
        import pytz
        # ISO 8601 format results in pytz.FixedOffset
        for freq in ['AS', 'W-SUN']:
            idx = date_range('2013-01-01T00:00:00-05:00',
                             '2016-01-01T23:59:59-05:00', freq=freq)
            expected = date_range('2013-01-01T00:00:00', '2016-01-01T23:59:59',
                                  freq=freq, tz=pytz.FixedOffset(-300))
            tm.assert_index_equal(idx, expected)
            # Unable to use `US/Eastern` because of DST
            expected_i8 = date_range('2013-01-01T00:00:00',
                                     '2016-01-01T23:59:59', freq=freq,
                                     tz='America/Lima')
            self.assert_numpy_array_equal(idx.asi8, expected_i8.asi8)

            idx = date_range('2013-01-01T00:00:00+09:00',
                             '2016-01-01T23:59:59+09:00', freq=freq)
            expected = date_range('2013-01-01T00:00:00', '2016-01-01T23:59:59',
                                  freq=freq, tz=pytz.FixedOffset(540))
            tm.assert_index_equal(idx, expected)
            expected_i8 = date_range('2013-01-01T00:00:00',
                                     '2016-01-01T23:59:59', freq=freq,
                                     tz='Asia/Tokyo')
            self.assert_numpy_array_equal(idx.asi8, expected_i8.asi8)

        tm._skip_if_no_dateutil()

        # Non ISO 8601 format results in dateutil.tz.tzoffset
        for freq in ['AS', 'W-SUN']:
            idx = date_range('2013/1/1 0:00:00-5:00', '2016/1/1 23:59:59-5:00',
                             freq=freq)
            expected = date_range('2013-01-01T00:00:00', '2016-01-01T23:59:59',
                                  freq=freq, tz=pytz.FixedOffset(-300))
            tm.assert_index_equal(idx, expected)
            # Unable to use `US/Eastern` because of DST
            expected_i8 = date_range('2013-01-01T00:00:00',
                                     '2016-01-01T23:59:59', freq=freq,
                                     tz='America/Lima')
            self.assert_numpy_array_equal(idx.asi8, expected_i8.asi8)

            idx = date_range('2013/1/1 0:00:00+9:00',
                             '2016/1/1 23:59:59+09:00', freq=freq)
            expected = date_range('2013-01-01T00:00:00', '2016-01-01T23:59:59',
                                  freq=freq, tz=pytz.FixedOffset(540))
            tm.assert_index_equal(idx, expected)
            expected_i8 = date_range('2013-01-01T00:00:00',
                                     '2016-01-01T23:59:59', freq=freq,
                                     tz='Asia/Tokyo')
            self.assert_numpy_array_equal(idx.asi8, expected_i8.asi8)

    def test_constructor_dtype(self):

        # passing a dtype with a tz should localize
        idx = DatetimeIndex(['2013-01-01', '2013-01-02'],
                            dtype='datetime64[ns, US/Eastern]')
        expected = DatetimeIndex(['2013-01-01', '2013-01-02']
                                 ).tz_localize('US/Eastern')
        tm.assert_index_equal(idx, expected)

        idx = DatetimeIndex(['2013-01-01', '2013-01-02'],
                            tz='US/Eastern')
        tm.assert_index_equal(idx, expected)

        # if we already have a tz and its not the same, then raise
        idx = DatetimeIndex(['2013-01-01', '2013-01-02'],
                            dtype='datetime64[ns, US/Eastern]')

        self.assertRaises(ValueError,
                          lambda: DatetimeIndex(idx,
                                                dtype='datetime64[ns]'))

        # this is effectively trying to convert tz's
        self.assertRaises(TypeError,
                          lambda: DatetimeIndex(idx,
                                                dtype='datetime64[ns, CET]'))
        self.assertRaises(ValueError,
                          lambda: DatetimeIndex(
                              idx, tz='CET',
                              dtype='datetime64[ns, US/Eastern]'))
        result = DatetimeIndex(idx, dtype='datetime64[ns, US/Eastern]')
        tm.assert_index_equal(idx, result)

    def test_constructor_name(self):
        idx = DatetimeIndex(start='2000-01-01', periods=1, freq='A',
                            name='TEST')
        self.assertEqual(idx.name, 'TEST')

    def test_comparisons_coverage(self):
        rng = date_range('1/1/2000', periods=10)

        # raise TypeError for now
        self.assertRaises(TypeError, rng.__lt__, rng[3].value)

        result = rng == list(rng)
        exp = rng == rng
        self.assert_numpy_array_equal(result, exp)

    def test_comparisons_nat(self):

        fidx1 = pd.Index([1.0, np.nan, 3.0, np.nan, 5.0, 7.0])
        fidx2 = pd.Index([2.0, 3.0, np.nan, np.nan, 6.0, 7.0])

        didx1 = pd.DatetimeIndex(['2014-01-01', pd.NaT, '2014-03-01', pd.NaT,
                                  '2014-05-01', '2014-07-01'])
        didx2 = pd.DatetimeIndex(['2014-02-01', '2014-03-01', pd.NaT, pd.NaT,
                                  '2014-06-01', '2014-07-01'])
        darr = np.array([np_datetime64_compat('2014-02-01 00:00Z'),
                         np_datetime64_compat('2014-03-01 00:00Z'),
                         np_datetime64_compat('nat'), np.datetime64('nat'),
                         np_datetime64_compat('2014-06-01 00:00Z'),
                         np_datetime64_compat('2014-07-01 00:00Z')])

        if _np_version_under1p8:
            # cannot test array because np.datetime('nat') returns today's date
            cases = [(fidx1, fidx2), (didx1, didx2)]
        else:
            cases = [(fidx1, fidx2), (didx1, didx2), (didx1, darr)]

        # Check pd.NaT is handles as the same as np.nan
        with tm.assert_produces_warning(None):
            for idx1, idx2 in cases:

                result = idx1 < idx2
                expected = np.array([True, False, False, False, True, False])
                self.assert_numpy_array_equal(result, expected)

                result = idx2 > idx1
                expected = np.array([True, False, False, False, True, False])
                self.assert_numpy_array_equal(result, expected)

                result = idx1 <= idx2
                expected = np.array([True, False, False, False, True, True])
                self.assert_numpy_array_equal(result, expected)

                result = idx2 >= idx1
                expected = np.array([True, False, False, False, True, True])
                self.assert_numpy_array_equal(result, expected)

                result = idx1 == idx2
                expected = np.array([False, False, False, False, False, True])
                self.assert_numpy_array_equal(result, expected)

                result = idx1 != idx2
                expected = np.array([True, True, True, True, True, False])
                self.assert_numpy_array_equal(result, expected)

        with tm.assert_produces_warning(None):
            for idx1, val in [(fidx1, np.nan), (didx1, pd.NaT)]:
                result = idx1 < val
                expected = np.array([False, False, False, False, False, False])
                self.assert_numpy_array_equal(result, expected)
                result = idx1 > val
                self.assert_numpy_array_equal(result, expected)

                result = idx1 <= val
                self.assert_numpy_array_equal(result, expected)
                result = idx1 >= val
                self.assert_numpy_array_equal(result, expected)

                result = idx1 == val
                self.assert_numpy_array_equal(result, expected)

                result = idx1 != val
                expected = np.array([True, True, True, True, True, True])
                self.assert_numpy_array_equal(result, expected)

        # Check pd.NaT is handles as the same as np.nan
        with tm.assert_produces_warning(None):
            for idx1, val in [(fidx1, 3), (didx1, datetime(2014, 3, 1))]:
                result = idx1 < val
                expected = np.array([True, False, False, False, False, False])
                self.assert_numpy_array_equal(result, expected)
                result = idx1 > val
                expected = np.array([False, False, False, False, True, True])
                self.assert_numpy_array_equal(result, expected)

                result = idx1 <= val
                expected = np.array([True, False, True, False, False, False])
                self.assert_numpy_array_equal(result, expected)
                result = idx1 >= val
                expected = np.array([False, False, True, False, True, True])
                self.assert_numpy_array_equal(result, expected)

                result = idx1 == val
                expected = np.array([False, False, True, False, False, False])
                self.assert_numpy_array_equal(result, expected)

                result = idx1 != val
                expected = np.array([True, True, False, True, True, True])
                self.assert_numpy_array_equal(result, expected)

    def test_map(self):
        rng = date_range('1/1/2000', periods=10)

        f = lambda x: x.strftime('%Y%m%d')
        result = rng.map(f)
        exp = np.array([f(x) for x in rng], dtype='<U8')
        tm.assert_almost_equal(result, exp)

    def test_iteration_preserves_tz(self):

        tm._skip_if_no_dateutil()

        # GH 8890
        import dateutil
        index = date_range("2012-01-01", periods=3, freq='H', tz='US/Eastern')

        for i, ts in enumerate(index):
            result = ts
            expected = index[i]
            self.assertEqual(result, expected)

        index = date_range("2012-01-01", periods=3, freq='H',
                           tz=dateutil.tz.tzoffset(None, -28800))

        for i, ts in enumerate(index):
            result = ts
            expected = index[i]
            self.assertEqual(result._repr_base, expected._repr_base)
            self.assertEqual(result, expected)

        # 9100
        index = pd.DatetimeIndex(['2014-12-01 03:32:39.987000-08:00',
                                  '2014-12-01 04:12:34.987000-08:00'])
        for i, ts in enumerate(index):
            result = ts
            expected = index[i]
            self.assertEqual(result._repr_base, expected._repr_base)
            self.assertEqual(result, expected)

    def test_misc_coverage(self):
        rng = date_range('1/1/2000', periods=5)
        result = rng.groupby(rng.day)
        tm.assertIsInstance(list(result.values())[0][0], Timestamp)

        idx = DatetimeIndex(['2000-01-03', '2000-01-01', '2000-01-02'])
        self.assertFalse(idx.equals(list(idx)))

        non_datetime = Index(list('abc'))
        self.assertFalse(idx.equals(list(non_datetime)))

    def test_union_coverage(self):
        idx = DatetimeIndex(['2000-01-03', '2000-01-01', '2000-01-02'])
        ordered = DatetimeIndex(idx.sort_values(), freq='infer')
        result = ordered.union(idx)
        tm.assert_index_equal(result, ordered)

        result = ordered[:0].union(ordered)
        tm.assert_index_equal(result, ordered)
        self.assertEqual(result.freq, ordered.freq)

    def test_union_bug_1730(self):
        rng_a = date_range('1/1/2012', periods=4, freq='3H')
        rng_b = date_range('1/1/2012', periods=4, freq='4H')

        result = rng_a.union(rng_b)
        exp = DatetimeIndex(sorted(set(list(rng_a)) | set(list(rng_b))))
        tm.assert_index_equal(result, exp)

    def test_union_bug_1745(self):
        left = DatetimeIndex(['2012-05-11 15:19:49.695000'])
        right = DatetimeIndex(['2012-05-29 13:04:21.322000',
                               '2012-05-11 15:27:24.873000',
                               '2012-05-11 15:31:05.350000'])

        result = left.union(right)
        exp = DatetimeIndex(sorted(set(list(left)) | set(list(right))))
        tm.assert_index_equal(result, exp)

    def test_union_bug_4564(self):
        from pandas import DateOffset
        left = date_range("2013-01-01", "2013-02-01")
        right = left + DateOffset(minutes=15)

        result = left.union(right)
        exp = DatetimeIndex(sorted(set(list(left)) | set(list(right))))
        tm.assert_index_equal(result, exp)

    def test_union_freq_both_none(self):
        # GH11086
        expected = bdate_range('20150101', periods=10)
        expected.freq = None

        result = expected.union(expected)
        tm.assert_index_equal(result, expected)
        self.assertIsNone(result.freq)

    def test_union_dataframe_index(self):
        rng1 = date_range('1/1/1999', '1/1/2012', freq='MS')
        s1 = Series(np.random.randn(len(rng1)), rng1)

        rng2 = date_range('1/1/1980', '12/1/2001', freq='MS')
        s2 = Series(np.random.randn(len(rng2)), rng2)
        df = DataFrame({'s1': s1, 's2': s2})

        exp = pd.date_range('1/1/1980', '1/1/2012', freq='MS')
        tm.assert_index_equal(df.index, exp)

    def test_intersection_bug_1708(self):
        from pandas import DateOffset
        index_1 = date_range('1/1/2012', periods=4, freq='12H')
        index_2 = index_1 + DateOffset(hours=1)

        result = index_1 & index_2
        self.assertEqual(len(result), 0)

    def test_intersection(self):
        # GH 4690 (with tz)
        for tz in [None, 'Asia/Tokyo', 'US/Eastern', 'dateutil/US/Pacific']:
            base = date_range('6/1/2000', '6/30/2000', freq='D', name='idx')

            # if target has the same name, it is preserved
            rng2 = date_range('5/15/2000', '6/20/2000', freq='D', name='idx')
            expected2 = date_range('6/1/2000', '6/20/2000', freq='D',
                                   name='idx')

            # if target name is different, it will be reset
            rng3 = date_range('5/15/2000', '6/20/2000', freq='D', name='other')
            expected3 = date_range('6/1/2000', '6/20/2000', freq='D',
                                   name=None)

            rng4 = date_range('7/1/2000', '7/31/2000', freq='D', name='idx')
            expected4 = DatetimeIndex([], name='idx')

            for (rng, expected) in [(rng2, expected2), (rng3, expected3),
                                    (rng4, expected4)]:
                result = base.intersection(rng)
                tm.assert_index_equal(result, expected)
                self.assertEqual(result.name, expected.name)
                self.assertEqual(result.freq, expected.freq)
                self.assertEqual(result.tz, expected.tz)

            # non-monotonic
            base = DatetimeIndex(['2011-01-05', '2011-01-04',
                                  '2011-01-02', '2011-01-03'],
                                 tz=tz, name='idx')

            rng2 = DatetimeIndex(['2011-01-04', '2011-01-02',
                                  '2011-02-02', '2011-02-03'],
                                 tz=tz, name='idx')
            expected2 = DatetimeIndex(
                ['2011-01-04', '2011-01-02'], tz=tz, name='idx')

            rng3 = DatetimeIndex(['2011-01-04', '2011-01-02',
                                  '2011-02-02', '2011-02-03'],
                                 tz=tz, name='other')
            expected3 = DatetimeIndex(
                ['2011-01-04', '2011-01-02'], tz=tz, name=None)

            # GH 7880
            rng4 = date_range('7/1/2000', '7/31/2000', freq='D', tz=tz,
                              name='idx')
            expected4 = DatetimeIndex([], tz=tz, name='idx')

            for (rng, expected) in [(rng2, expected2), (rng3, expected3),
                                    (rng4, expected4)]:
                result = base.intersection(rng)
                tm.assert_index_equal(result, expected)
                self.assertEqual(result.name, expected.name)
                self.assertIsNone(result.freq)
                self.assertEqual(result.tz, expected.tz)

        # empty same freq GH2129
        rng = date_range('6/1/2000', '6/15/2000', freq='T')
        result = rng[0:0].intersection(rng)
        self.assertEqual(len(result), 0)

        result = rng.intersection(rng[0:0])
        self.assertEqual(len(result), 0)

    def test_string_index_series_name_converted(self):
        # #1644
        df = DataFrame(np.random.randn(10, 4),
                       index=date_range('1/1/2000', periods=10))

        result = df.ix['1/3/2000']
        self.assertEqual(result.name, df.index[2])

        result = df.T['1/3/2000']
        self.assertEqual(result.name, df.index[2])

    # GH 10699
    def test_datetime64_with_DateOffset(self):
        for klass, assert_func in zip([Series, DatetimeIndex],
                                      [self.assert_series_equal,
                                       tm.assert_index_equal]):
            s = klass(date_range('2000-01-01', '2000-01-31'), name='a')
            result = s + pd.DateOffset(years=1)
            result2 = pd.DateOffset(years=1) + s
            exp = klass(date_range('2001-01-01', '2001-01-31'), name='a')
            assert_func(result, exp)
            assert_func(result2, exp)

            result = s - pd.DateOffset(years=1)
            exp = klass(date_range('1999-01-01', '1999-01-31'), name='a')
            assert_func(result, exp)

            s = klass([Timestamp('2000-01-15 00:15:00', tz='US/Central'),
                       pd.Timestamp('2000-02-15', tz='US/Central')], name='a')
            result = s + pd.offsets.Day()
            result2 = pd.offsets.Day() + s
            exp = klass([Timestamp('2000-01-16 00:15:00', tz='US/Central'),
                         Timestamp('2000-02-16', tz='US/Central')], name='a')
            assert_func(result, exp)
            assert_func(result2, exp)

            s = klass([Timestamp('2000-01-15 00:15:00', tz='US/Central'),
                       pd.Timestamp('2000-02-15', tz='US/Central')], name='a')
            result = s + pd.offsets.MonthEnd()
            result2 = pd.offsets.MonthEnd() + s
            exp = klass([Timestamp('2000-01-31 00:15:00', tz='US/Central'),
                         Timestamp('2000-02-29', tz='US/Central')], name='a')
            assert_func(result, exp)
            assert_func(result2, exp)

            # array of offsets - valid for Series only
            if klass is Series:
                with tm.assert_produces_warning(PerformanceWarning):
                    s = klass([Timestamp('2000-1-1'), Timestamp('2000-2-1')])
                    result = s + Series([pd.offsets.DateOffset(years=1),
                                         pd.offsets.MonthEnd()])
                    exp = klass([Timestamp('2001-1-1'), Timestamp('2000-2-29')
                                 ])
                    assert_func(result, exp)

                    # same offset
                    result = s + Series([pd.offsets.DateOffset(years=1),
                                         pd.offsets.DateOffset(years=1)])
                    exp = klass([Timestamp('2001-1-1'), Timestamp('2001-2-1')])
                    assert_func(result, exp)

            s = klass([Timestamp('2000-01-05 00:15:00'),
                       Timestamp('2000-01-31 00:23:00'),
                       Timestamp('2000-01-01'),
                       Timestamp('2000-03-31'),
                       Timestamp('2000-02-29'),
                       Timestamp('2000-12-31'),
                       Timestamp('2000-05-15'),
                       Timestamp('2001-06-15')])

            # DateOffset relativedelta fastpath
            relative_kwargs = [('years', 2), ('months', 5), ('days', 3),
                               ('hours', 5), ('minutes', 10), ('seconds', 2),
                               ('microseconds', 5)]
            for i, kwd in enumerate(relative_kwargs):
                op = pd.DateOffset(**dict([kwd]))
                assert_func(klass([x + op for x in s]), s + op)
                assert_func(klass([x - op for x in s]), s - op)
                op = pd.DateOffset(**dict(relative_kwargs[:i + 1]))
                assert_func(klass([x + op for x in s]), s + op)
                assert_func(klass([x - op for x in s]), s - op)

            # assert these are equal on a piecewise basis
            offsets = ['YearBegin', ('YearBegin', {'month': 5}), 'YearEnd',
                       ('YearEnd', {'month': 5}), 'MonthBegin', 'MonthEnd',
                       'SemiMonthEnd', 'SemiMonthBegin',
                       'Week', ('Week', {
                           'weekday': 3
                       }), 'BusinessDay', 'BDay', 'QuarterEnd', 'QuarterBegin',
                       'CustomBusinessDay', 'CDay', 'CBMonthEnd',
                       'CBMonthBegin', 'BMonthBegin', 'BMonthEnd',
                       'BusinessHour', 'BYearBegin', 'BYearEnd',
                       'BQuarterBegin', ('LastWeekOfMonth', {
                           'weekday': 2
                       }), ('FY5253Quarter', {'qtr_with_extra_week': 1,
                                              'startingMonth': 1,
                                              'weekday': 2,
                                              'variation': 'nearest'}),
                       ('FY5253', {'weekday': 0,
                                   'startingMonth': 2,
                                   'variation':
                                   'nearest'}), ('WeekOfMonth', {'weekday': 2,
                                                                 'week': 2}),
                       'Easter', ('DateOffset', {'day': 4}),
                       ('DateOffset', {'month': 5})]

            with warnings.catch_warnings(record=True):
                for normalize in (True, False):
                    for do in offsets:
                        if isinstance(do, tuple):
                            do, kwargs = do
                        else:
                            do = do
                            kwargs = {}

                        for n in [0, 5]:
                            if (do in ['WeekOfMonth', 'LastWeekOfMonth',
                                       'FY5253Quarter', 'FY5253'] and n == 0):
                                continue
                            op = getattr(pd.offsets, do)(n,
                                                         normalize=normalize,
                                                         **kwargs)
                            assert_func(klass([x + op for x in s]), s + op)
                            assert_func(klass([x - op for x in s]), s - op)
                            assert_func(klass([op + x for x in s]), op + s)
    # def test_add_timedelta64(self):
    #     rng = date_range('1/1/2000', periods=5)
    #     delta = rng.values[3] - rng.values[1]

    #     result = rng + delta
    #     expected = rng + timedelta(2)
    #     self.assertTrue(result.equals(expected))

    def test_get_duplicates(self):
        idx = DatetimeIndex(['2000-01-01', '2000-01-02', '2000-01-02',
                             '2000-01-03', '2000-01-03', '2000-01-04'])

        result = idx.get_duplicates()
        ex = DatetimeIndex(['2000-01-02', '2000-01-03'])
        tm.assert_index_equal(result, ex)

    def test_argmin_argmax(self):
        idx = DatetimeIndex(['2000-01-04', '2000-01-01', '2000-01-02'])
        self.assertEqual(idx.argmin(), 1)
        self.assertEqual(idx.argmax(), 0)

    def test_sort_values(self):
        idx = DatetimeIndex(['2000-01-04', '2000-01-01', '2000-01-02'])

        ordered = idx.sort_values()
        self.assertTrue(ordered.is_monotonic)

        ordered = idx.sort_values(ascending=False)
        self.assertTrue(ordered[::-1].is_monotonic)

        ordered, dexer = idx.sort_values(return_indexer=True)
        self.assertTrue(ordered.is_monotonic)
        self.assert_numpy_array_equal(dexer,
                                      np.array([1, 2, 0], dtype=np.intp))

        ordered, dexer = idx.sort_values(return_indexer=True, ascending=False)
        self.assertTrue(ordered[::-1].is_monotonic)
        self.assert_numpy_array_equal(dexer,
                                      np.array([0, 2, 1], dtype=np.intp))

    def test_round(self):

        # round
        dt = Timestamp('20130101 09:10:11')
        result = dt.round('D')
        expected = Timestamp('20130101')
        self.assertEqual(result, expected)

        dt = Timestamp('20130101 19:10:11')
        result = dt.round('D')
        expected = Timestamp('20130102')
        self.assertEqual(result, expected)

        dt = Timestamp('20130201 12:00:00')
        result = dt.round('D')
        expected = Timestamp('20130202')
        self.assertEqual(result, expected)

        dt = Timestamp('20130104 12:00:00')
        result = dt.round('D')
        expected = Timestamp('20130105')
        self.assertEqual(result, expected)

        dt = Timestamp('20130104 12:32:00')
        result = dt.round('30Min')
        expected = Timestamp('20130104 12:30:00')
        self.assertEqual(result, expected)

        dti = date_range('20130101 09:10:11', periods=5)
        result = dti.round('D')
        expected = date_range('20130101', periods=5)
        tm.assert_index_equal(result, expected)

        # floor
        dt = Timestamp('20130101 09:10:11')
        result = dt.floor('D')
        expected = Timestamp('20130101')
        self.assertEqual(result, expected)

        # ceil
        dt = Timestamp('20130101 09:10:11')
        result = dt.ceil('D')
        expected = Timestamp('20130102')
        self.assertEqual(result, expected)

        # round with tz
        dt = Timestamp('20130101 09:10:11', tz='US/Eastern')
        result = dt.round('D')
        expected = Timestamp('20130101', tz='US/Eastern')
        self.assertEqual(result, expected)

        dt = Timestamp('20130101 09:10:11', tz='US/Eastern')
        result = dt.round('s')
        self.assertEqual(result, dt)

        dti = date_range('20130101 09:10:11',
                         periods=5).tz_localize('UTC').tz_convert('US/Eastern')
        result = dti.round('D')
        expected = date_range('20130101', periods=5).tz_localize('US/Eastern')
        tm.assert_index_equal(result, expected)

        result = dti.round('s')
        tm.assert_index_equal(result, dti)

        # invalid
        for freq in ['Y', 'M', 'foobar']:
            self.assertRaises(ValueError, lambda: dti.round(freq))

    def test_insert(self):
        idx = DatetimeIndex(
            ['2000-01-04', '2000-01-01', '2000-01-02'], name='idx')

        result = idx.insert(2, datetime(2000, 1, 5))
        exp = DatetimeIndex(['2000-01-04', '2000-01-01', '2000-01-05',
                             '2000-01-02'], name='idx')
        tm.assert_index_equal(result, exp)

        # insertion of non-datetime should coerce to object index
        result = idx.insert(1, 'inserted')
        expected = Index([datetime(2000, 1, 4), 'inserted',
                          datetime(2000, 1, 1),
                          datetime(2000, 1, 2)], name='idx')
        self.assertNotIsInstance(result, DatetimeIndex)
        tm.assert_index_equal(result, expected)
        self.assertEqual(result.name, expected.name)

        idx = date_range('1/1/2000', periods=3, freq='M', name='idx')

        # preserve freq
        expected_0 = DatetimeIndex(['1999-12-31', '2000-01-31', '2000-02-29',
                                    '2000-03-31'], name='idx', freq='M')
        expected_3 = DatetimeIndex(['2000-01-31', '2000-02-29', '2000-03-31',
                                    '2000-04-30'], name='idx', freq='M')

        # reset freq to None
        expected_1_nofreq = DatetimeIndex(['2000-01-31', '2000-01-31',
                                           '2000-02-29',
                                           '2000-03-31'], name='idx',
                                          freq=None)
        expected_3_nofreq = DatetimeIndex(['2000-01-31', '2000-02-29',
                                           '2000-03-31',
                                           '2000-01-02'], name='idx',
                                          freq=None)

        cases = [(0, datetime(1999, 12, 31), expected_0),
                 (-3, datetime(1999, 12, 31), expected_0),
                 (3, datetime(2000, 4, 30), expected_3),
                 (1, datetime(2000, 1, 31), expected_1_nofreq),
                 (3, datetime(2000, 1, 2), expected_3_nofreq)]

        for n, d, expected in cases:
            result = idx.insert(n, d)
            tm.assert_index_equal(result, expected)
            self.assertEqual(result.name, expected.name)
            self.assertEqual(result.freq, expected.freq)

        # reset freq to None
        result = idx.insert(3, datetime(2000, 1, 2))
        expected = DatetimeIndex(['2000-01-31', '2000-02-29', '2000-03-31',
                                  '2000-01-02'], name='idx', freq=None)
        tm.assert_index_equal(result, expected)
        self.assertEqual(result.name, expected.name)
        self.assertTrue(result.freq is None)

        # GH 7299
        tm._skip_if_no_pytz()
        import pytz

        idx = date_range('1/1/2000', periods=3, freq='D', tz='Asia/Tokyo',
                         name='idx')
        with tm.assertRaises(ValueError):
            result = idx.insert(3, pd.Timestamp('2000-01-04'))
        with tm.assertRaises(ValueError):
            result = idx.insert(3, datetime(2000, 1, 4))
        with tm.assertRaises(ValueError):
            result = idx.insert(3, pd.Timestamp('2000-01-04', tz='US/Eastern'))
        with tm.assertRaises(ValueError):
            result = idx.insert(3,
                                datetime(2000, 1, 4,
                                         tzinfo=pytz.timezone('US/Eastern')))

        for tz in ['US/Pacific', 'Asia/Singapore']:
            idx = date_range('1/1/2000 09:00', periods=6, freq='H', tz=tz,
                             name='idx')
            # preserve freq
            expected = date_range('1/1/2000 09:00', periods=7, freq='H', tz=tz,
                                  name='idx')
            for d in [pd.Timestamp('2000-01-01 15:00', tz=tz),
                      pytz.timezone(tz).localize(datetime(2000, 1, 1, 15))]:

                result = idx.insert(6, d)
                tm.assert_index_equal(result, expected)
                self.assertEqual(result.name, expected.name)
                self.assertEqual(result.freq, expected.freq)
                self.assertEqual(result.tz, expected.tz)

            expected = DatetimeIndex(['2000-01-01 09:00', '2000-01-01 10:00',
                                      '2000-01-01 11:00',
                                      '2000-01-01 12:00', '2000-01-01 13:00',
                                      '2000-01-01 14:00',
                                      '2000-01-01 10:00'], name='idx',
                                     tz=tz, freq=None)
            # reset freq to None
            for d in [pd.Timestamp('2000-01-01 10:00', tz=tz),
                      pytz.timezone(tz).localize(datetime(2000, 1, 1, 10))]:
                result = idx.insert(6, d)
                tm.assert_index_equal(result, expected)
                self.assertEqual(result.name, expected.name)
                self.assertTrue(result.freq is None)
                self.assertEqual(result.tz, expected.tz)

    def test_delete(self):
        idx = date_range(start='2000-01-01', periods=5, freq='M', name='idx')

        # prserve freq
        expected_0 = date_range(start='2000-02-01', periods=4, freq='M',
                                name='idx')
        expected_4 = date_range(start='2000-01-01', periods=4, freq='M',
                                name='idx')

        # reset freq to None
        expected_1 = DatetimeIndex(['2000-01-31', '2000-03-31', '2000-04-30',
                                    '2000-05-31'], freq=None, name='idx')

        cases = {0: expected_0,
                 -5: expected_0,
                 -1: expected_4,
                 4: expected_4,
                 1: expected_1}
        for n, expected in compat.iteritems(cases):
            result = idx.delete(n)
            tm.assert_index_equal(result, expected)
            self.assertEqual(result.name, expected.name)
            self.assertEqual(result.freq, expected.freq)

        with tm.assertRaises((IndexError, ValueError)):
            # either depeidnig on numpy version
            result = idx.delete(5)

        for tz in [None, 'Asia/Tokyo', 'US/Pacific']:
            idx = date_range(start='2000-01-01 09:00', periods=10, freq='H',
                             name='idx', tz=tz)

            expected = date_range(start='2000-01-01 10:00', periods=9,
                                  freq='H', name='idx', tz=tz)
            result = idx.delete(0)
            tm.assert_index_equal(result, expected)
            self.assertEqual(result.name, expected.name)
            self.assertEqual(result.freqstr, 'H')
            self.assertEqual(result.tz, expected.tz)

            expected = date_range(start='2000-01-01 09:00', periods=9,
                                  freq='H', name='idx', tz=tz)
            result = idx.delete(-1)
            tm.assert_index_equal(result, expected)
            self.assertEqual(result.name, expected.name)
            self.assertEqual(result.freqstr, 'H')
            self.assertEqual(result.tz, expected.tz)

    def test_delete_slice(self):
        idx = date_range(start='2000-01-01', periods=10, freq='D', name='idx')

        # prserve freq
        expected_0_2 = date_range(start='2000-01-04', periods=7, freq='D',
                                  name='idx')
        expected_7_9 = date_range(start='2000-01-01', periods=7, freq='D',
                                  name='idx')

        # reset freq to None
        expected_3_5 = DatetimeIndex(['2000-01-01', '2000-01-02', '2000-01-03',
                                      '2000-01-07', '2000-01-08', '2000-01-09',
                                      '2000-01-10'], freq=None, name='idx')

        cases = {(0, 1, 2): expected_0_2,
                 (7, 8, 9): expected_7_9,
                 (3, 4, 5): expected_3_5}
        for n, expected in compat.iteritems(cases):
            result = idx.delete(n)
            tm.assert_index_equal(result, expected)
            self.assertEqual(result.name, expected.name)
            self.assertEqual(result.freq, expected.freq)

            result = idx.delete(slice(n[0], n[-1] + 1))
            tm.assert_index_equal(result, expected)
            self.assertEqual(result.name, expected.name)
            self.assertEqual(result.freq, expected.freq)

        for tz in [None, 'Asia/Tokyo', 'US/Pacific']:
            ts = pd.Series(1, index=pd.date_range(
                '2000-01-01 09:00', periods=10, freq='H', name='idx', tz=tz))
            # preserve freq
            result = ts.drop(ts.index[:5]).index
            expected = pd.date_range('2000-01-01 14:00', periods=5, freq='H',
                                     name='idx', tz=tz)
            tm.assert_index_equal(result, expected)
            self.assertEqual(result.name, expected.name)
            self.assertEqual(result.freq, expected.freq)
            self.assertEqual(result.tz, expected.tz)

            # reset freq to None
            result = ts.drop(ts.index[[1, 3, 5, 7, 9]]).index
            expected = DatetimeIndex(['2000-01-01 09:00', '2000-01-01 11:00',
                                      '2000-01-01 13:00',
                                      '2000-01-01 15:00', '2000-01-01 17:00'],
                                     freq=None, name='idx', tz=tz)
            tm.assert_index_equal(result, expected)
            self.assertEqual(result.name, expected.name)
            self.assertEqual(result.freq, expected.freq)
            self.assertEqual(result.tz, expected.tz)

    def test_take(self):
        dates = [datetime(2010, 1, 1, 14), datetime(2010, 1, 1, 15),
                 datetime(2010, 1, 1, 17), datetime(2010, 1, 1, 21)]

        for tz in [None, 'US/Eastern', 'Asia/Tokyo']:
            idx = DatetimeIndex(start='2010-01-01 09:00',
                                end='2010-02-01 09:00', freq='H', tz=tz,
                                name='idx')
            expected = DatetimeIndex(dates, freq=None, name='idx', tz=tz)

            taken1 = idx.take([5, 6, 8, 12])
            taken2 = idx[[5, 6, 8, 12]]

            for taken in [taken1, taken2]:
                tm.assert_index_equal(taken, expected)
                tm.assertIsInstance(taken, DatetimeIndex)
                self.assertIsNone(taken.freq)
                self.assertEqual(taken.tz, expected.tz)
                self.assertEqual(taken.name, expected.name)

    def test_take_fill_value(self):
        # GH 12631
        idx = pd.DatetimeIndex(['2011-01-01', '2011-02-01', '2011-03-01'],
                               name='xxx')
        result = idx.take(np.array([1, 0, -1]))
        expected = pd.DatetimeIndex(['2011-02-01', '2011-01-01', '2011-03-01'],
                                    name='xxx')
        tm.assert_index_equal(result, expected)

        # fill_value
        result = idx.take(np.array([1, 0, -1]), fill_value=True)
        expected = pd.DatetimeIndex(['2011-02-01', '2011-01-01', 'NaT'],
                                    name='xxx')
        tm.assert_index_equal(result, expected)

        # allow_fill=False
        result = idx.take(np.array([1, 0, -1]), allow_fill=False,
                          fill_value=True)
        expected = pd.DatetimeIndex(['2011-02-01', '2011-01-01', '2011-03-01'],
                                    name='xxx')
        tm.assert_index_equal(result, expected)

        msg = ('When allow_fill=True and fill_value is not None, '
               'all indices must be >= -1')
        with tm.assertRaisesRegexp(ValueError, msg):
            idx.take(np.array([1, 0, -2]), fill_value=True)
        with tm.assertRaisesRegexp(ValueError, msg):
            idx.take(np.array([1, 0, -5]), fill_value=True)

        with tm.assertRaises(IndexError):
            idx.take(np.array([1, -5]))

    def test_take_fill_value_with_timezone(self):
        idx = pd.DatetimeIndex(['2011-01-01', '2011-02-01', '2011-03-01'],
                               name='xxx', tz='US/Eastern')
        result = idx.take(np.array([1, 0, -1]))
        expected = pd.DatetimeIndex(['2011-02-01', '2011-01-01', '2011-03-01'],
                                    name='xxx', tz='US/Eastern')
        tm.assert_index_equal(result, expected)

        # fill_value
        result = idx.take(np.array([1, 0, -1]), fill_value=True)
        expected = pd.DatetimeIndex(['2011-02-01', '2011-01-01', 'NaT'],
                                    name='xxx', tz='US/Eastern')
        tm.assert_index_equal(result, expected)

        # allow_fill=False
        result = idx.take(np.array([1, 0, -1]), allow_fill=False,
                          fill_value=True)
        expected = pd.DatetimeIndex(['2011-02-01', '2011-01-01', '2011-03-01'],
                                    name='xxx', tz='US/Eastern')
        tm.assert_index_equal(result, expected)

        msg = ('When allow_fill=True and fill_value is not None, '
               'all indices must be >= -1')
        with tm.assertRaisesRegexp(ValueError, msg):
            idx.take(np.array([1, 0, -2]), fill_value=True)
        with tm.assertRaisesRegexp(ValueError, msg):
            idx.take(np.array([1, 0, -5]), fill_value=True)

        with tm.assertRaises(IndexError):
            idx.take(np.array([1, -5]))

    def test_map_bug_1677(self):
        index = DatetimeIndex(['2012-04-25 09:30:00.393000'])
        f = index.asof

        result = index.map(f)
        expected = np.array([f(index[0])])
        self.assert_numpy_array_equal(result, expected)

    def test_groupby_function_tuple_1677(self):
        df = DataFrame(np.random.rand(100),
                       index=date_range("1/1/2000", periods=100))
        monthly_group = df.groupby(lambda x: (x.year, x.month))

        result = monthly_group.mean()
        tm.assertIsInstance(result.index[0], tuple)

    def test_append_numpy_bug_1681(self):
        # another datetime64 bug
        dr = date_range('2011/1/1', '2012/1/1', freq='W-FRI')
        a = DataFrame()
        c = DataFrame({'A': 'foo', 'B': dr}, index=dr)

        result = a.append(c)
        self.assertTrue((result['B'] == dr).all())

    def test_isin(self):
        index = tm.makeDateIndex(4)
        result = index.isin(index)
        self.assertTrue(result.all())

        result = index.isin(list(index))
        self.assertTrue(result.all())

        assert_almost_equal(index.isin([index[2], 5]),
                            np.array([False, False, True, False]))

    def test_union(self):
        i1 = Int64Index(np.arange(0, 20, 2))
        i2 = Int64Index(np.arange(10, 30, 2))
        result = i1.union(i2)
        expected = Int64Index(np.arange(0, 30, 2))
        tm.assert_index_equal(result, expected)

    def test_union_with_DatetimeIndex(self):
        i1 = Int64Index(np.arange(0, 20, 2))
        i2 = DatetimeIndex(start='2012-01-03 00:00:00', periods=10, freq='D')
        i1.union(i2)  # Works
        i2.union(i1)  # Fails with "AttributeError: can't set attribute"

    def test_time(self):
        rng = pd.date_range('1/1/2000', freq='12min', periods=10)
        result = pd.Index(rng).time
        expected = [t.time() for t in rng]
        self.assertTrue((result == expected).all())

    def test_date(self):
        rng = pd.date_range('1/1/2000', freq='12H', periods=10)
        result = pd.Index(rng).date
        expected = [t.date() for t in rng]
        self.assertTrue((result == expected).all())

    def test_does_not_convert_mixed_integer(self):
        df = tm.makeCustomDataframe(10, 10,
                                    data_gen_f=lambda *args, **kwargs: randn(),
                                    r_idx_type='i', c_idx_type='dt')
        cols = df.columns.join(df.index, how='outer')
        joined = cols.join(df.columns)
        self.assertEqual(cols.dtype, np.dtype('O'))
        self.assertEqual(cols.dtype, joined.dtype)
        tm.assert_numpy_array_equal(cols.values, joined.values)

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
            self.assertIs(index, joined)

    def assert_index_parameters(self, index):
        assert index.freq == '40960N'
        assert index.inferred_freq == '40960N'

    def test_ns_index(self):
        nsamples = 400
        ns = int(1e9 / 24414)
        dtstart = np.datetime64('2012-09-20T00:00:00')

        dt = dtstart + np.arange(nsamples) * np.timedelta64(ns, 'ns')
        freq = ns * offsets.Nano()
        index = pd.DatetimeIndex(dt, freq=freq, name='time')
        self.assert_index_parameters(index)

        new_index = pd.DatetimeIndex(start=index[0], end=index[-1],
                                     freq=index.freq)
        self.assert_index_parameters(new_index)

    def test_join_with_period_index(self):
        df = tm.makeCustomDataframe(
            10, 10, data_gen_f=lambda *args: np.random.randint(2),
            c_idx_type='p', r_idx_type='dt')
        s = df.iloc[:5, 0]
        joins = 'left', 'right', 'inner', 'outer'

        for join in joins:
            with tm.assertRaisesRegexp(ValueError, 'can only call with other '
                                       'PeriodIndex-ed objects'):
                df.columns.join(s.index, how=join)

    def test_factorize(self):
        idx1 = DatetimeIndex(['2014-01', '2014-01', '2014-02', '2014-02',
                              '2014-03', '2014-03'])

        exp_arr = np.array([0, 0, 1, 1, 2, 2], dtype=np.intp)
        exp_idx = DatetimeIndex(['2014-01', '2014-02', '2014-03'])

        arr, idx = idx1.factorize()
        self.assert_numpy_array_equal(arr, exp_arr)
        tm.assert_index_equal(idx, exp_idx)

        arr, idx = idx1.factorize(sort=True)
        self.assert_numpy_array_equal(arr, exp_arr)
        tm.assert_index_equal(idx, exp_idx)

        # tz must be preserved
        idx1 = idx1.tz_localize('Asia/Tokyo')
        exp_idx = exp_idx.tz_localize('Asia/Tokyo')

        arr, idx = idx1.factorize()
        self.assert_numpy_array_equal(arr, exp_arr)
        tm.assert_index_equal(idx, exp_idx)

        idx2 = pd.DatetimeIndex(['2014-03', '2014-03', '2014-02', '2014-01',
                                 '2014-03', '2014-01'])

        exp_arr = np.array([2, 2, 1, 0, 2, 0], dtype=np.intp)
        exp_idx = DatetimeIndex(['2014-01', '2014-02', '2014-03'])
        arr, idx = idx2.factorize(sort=True)
        self.assert_numpy_array_equal(arr, exp_arr)
        tm.assert_index_equal(idx, exp_idx)

        exp_arr = np.array([0, 0, 1, 2, 0, 2], dtype=np.intp)
        exp_idx = DatetimeIndex(['2014-03', '2014-02', '2014-01'])
        arr, idx = idx2.factorize()
        self.assert_numpy_array_equal(arr, exp_arr)
        tm.assert_index_equal(idx, exp_idx)

        # freq must be preserved
        idx3 = date_range('2000-01', periods=4, freq='M', tz='Asia/Tokyo')
        exp_arr = np.array([0, 1, 2, 3], dtype=np.intp)
        arr, idx = idx3.factorize()
        self.assert_numpy_array_equal(arr, exp_arr)
        tm.assert_index_equal(idx, idx3)

    def test_factorize_tz(self):
        # GH 13750
        for tz in [None, 'UTC', 'US/Eastern', 'Asia/Tokyo']:
            base = pd.date_range('2016-11-05', freq='H', periods=100, tz=tz)
            idx = base.repeat(5)

            exp_arr = np.arange(100, dtype=np.intp).repeat(5)

            for obj in [idx, pd.Series(idx)]:
                arr, res = obj.factorize()
                self.assert_numpy_array_equal(arr, exp_arr)
                tm.assert_index_equal(res, base)

    def test_factorize_dst(self):
        # GH 13750
        idx = pd.date_range('2016-11-06', freq='H', periods=12,
                            tz='US/Eastern')

        for obj in [idx, pd.Series(idx)]:
            arr, res = obj.factorize()
            self.assert_numpy_array_equal(arr, np.arange(12, dtype=np.intp))
            tm.assert_index_equal(res, idx)

        idx = pd.date_range('2016-06-13', freq='H', periods=12,
                            tz='US/Eastern')

        for obj in [idx, pd.Series(idx)]:
            arr, res = obj.factorize()
            self.assert_numpy_array_equal(arr, np.arange(12, dtype=np.intp))
            tm.assert_index_equal(res, idx)

    def test_slice_with_negative_step(self):
        ts = Series(np.arange(20),
                    date_range('2014-01-01', periods=20, freq='MS'))
        SLC = pd.IndexSlice

        def assert_slices_equivalent(l_slc, i_slc):
            assert_series_equal(ts[l_slc], ts.iloc[i_slc])
            assert_series_equal(ts.loc[l_slc], ts.iloc[i_slc])
            assert_series_equal(ts.ix[l_slc], ts.iloc[i_slc])

        assert_slices_equivalent(SLC[Timestamp('2014-10-01')::-1], SLC[9::-1])
        assert_slices_equivalent(SLC['2014-10-01'::-1], SLC[9::-1])

        assert_slices_equivalent(SLC[:Timestamp('2014-10-01'):-1], SLC[:8:-1])
        assert_slices_equivalent(SLC[:'2014-10-01':-1], SLC[:8:-1])

        assert_slices_equivalent(SLC['2015-02-01':'2014-10-01':-1],
                                 SLC[13:8:-1])
        assert_slices_equivalent(SLC[Timestamp('2015-02-01'):Timestamp(
            '2014-10-01'):-1], SLC[13:8:-1])
        assert_slices_equivalent(SLC['2015-02-01':Timestamp('2014-10-01'):-1],
                                 SLC[13:8:-1])
        assert_slices_equivalent(SLC[Timestamp('2015-02-01'):'2014-10-01':-1],
                                 SLC[13:8:-1])

        assert_slices_equivalent(SLC['2014-10-01':'2015-02-01':-1], SLC[:0])

    def test_slice_with_zero_step_raises(self):
        ts = Series(np.arange(20),
                    date_range('2014-01-01', periods=20, freq='MS'))
        self.assertRaisesRegexp(ValueError, 'slice step cannot be zero',
                                lambda: ts[::0])
        self.assertRaisesRegexp(ValueError, 'slice step cannot be zero',
                                lambda: ts.loc[::0])
        self.assertRaisesRegexp(ValueError, 'slice step cannot be zero',
                                lambda: ts.ix[::0])

    def test_slice_bounds_empty(self):
        # GH 14354
        empty_idx = DatetimeIndex(freq='1H', periods=0, end='2015')

        right = empty_idx._maybe_cast_slice_bound('2015-01-02', 'right', 'loc')
        exp = Timestamp('2015-01-02 23:59:59.999999999')
        self.assertEqual(right, exp)

        left = empty_idx._maybe_cast_slice_bound('2015-01-02', 'left', 'loc')
        exp = Timestamp('2015-01-02 00:00:00')
        self.assertEqual(left, exp)


class TestDatetime64(tm.TestCase):
    """
    Also test support for datetime64[ns] in Series / DataFrame
    """

    def setUp(self):
        dti = DatetimeIndex(start=datetime(2005, 1, 1),
                            end=datetime(2005, 1, 10), freq='Min')
        self.series = Series(rand(len(dti)), dti)

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

    def test_nanosecond_field(self):
        dti = DatetimeIndex(np.arange(10))

        self.assert_numpy_array_equal(dti.nanosecond,
                                      np.arange(10, dtype=np.int32))

    def test_datetimeindex_diff(self):
        dti1 = DatetimeIndex(freq='Q-JAN', start=datetime(1997, 12, 31),
                             periods=100)
        dti2 = DatetimeIndex(freq='Q-JAN', start=datetime(1997, 12, 31),
                             periods=98)
        self.assertEqual(len(dti1.difference(dti2)), 2)

    def test_fancy_getitem(self):
        dti = DatetimeIndex(freq='WOM-1FRI', start=datetime(2005, 1, 1),
                            end=datetime(2010, 1, 1))

        s = Series(np.arange(len(dti)), index=dti)

        self.assertEqual(s[48], 48)
        self.assertEqual(s['1/2/2009'], 48)
        self.assertEqual(s['2009-1-2'], 48)
        self.assertEqual(s[datetime(2009, 1, 2)], 48)
        self.assertEqual(s[lib.Timestamp(datetime(2009, 1, 2))], 48)
        self.assertRaises(KeyError, s.__getitem__, '2009-1-3')

        assert_series_equal(s['3/6/2009':'2009-06-05'],
                            s[datetime(2009, 3, 6):datetime(2009, 6, 5)])

    def test_fancy_setitem(self):
        dti = DatetimeIndex(freq='WOM-1FRI', start=datetime(2005, 1, 1),
                            end=datetime(2010, 1, 1))

        s = Series(np.arange(len(dti)), index=dti)
        s[48] = -1
        self.assertEqual(s[48], -1)
        s['1/2/2009'] = -2
        self.assertEqual(s[48], -2)
        s['1/2/2009':'2009-06-05'] = -3
        self.assertTrue((s[48:54] == -3).all())

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

    def test_dti_snap(self):
        dti = DatetimeIndex(['1/1/2002', '1/2/2002', '1/3/2002', '1/4/2002',
                             '1/5/2002', '1/6/2002', '1/7/2002'], freq='D')

        res = dti.snap(freq='W-MON')
        exp = date_range('12/31/2001', '1/7/2002', freq='w-mon')
        exp = exp.repeat([3, 4])
        self.assertTrue((res == exp).all())

        res = dti.snap(freq='B')

        exp = date_range('1/1/2002', '1/7/2002', freq='b')
        exp = exp.repeat([1, 1, 1, 2, 2])
        self.assertTrue((res == exp).all())

    def test_dti_reset_index_round_trip(self):
        dti = DatetimeIndex(start='1/1/2001', end='6/1/2001', freq='D')
        d1 = DataFrame({'v': np.random.rand(len(dti))}, index=dti)
        d2 = d1.reset_index()
        self.assertEqual(d2.dtypes[0], np.dtype('M8[ns]'))
        d3 = d2.set_index('index')
        assert_frame_equal(d1, d3, check_names=False)

        # #2329
        stamp = datetime(2012, 11, 22)
        df = DataFrame([[stamp, 12.1]], columns=['Date', 'Value'])
        df = df.set_index('Date')

        self.assertEqual(df.index[0], stamp)
        self.assertEqual(df.reset_index()['Date'][0], stamp)

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
        # self.assertEqual(s2.values.dtype, 'M8[ns]')

    @slow
    def test_slice_locs_indexerror(self):
        times = [datetime(2000, 1, 1) + timedelta(minutes=i * 10)
                 for i in range(100000)]
        s = Series(lrange(100000), times)
        s.ix[datetime(1900, 1, 1):datetime(2100, 1, 1)]

    def test_slicing_datetimes(self):

        # GH 7523

        # unique
        df = DataFrame(np.arange(4., dtype='float64'),
                       index=[datetime(2001, 1, i, 10, 00)
                              for i in [1, 2, 3, 4]])
        result = df.ix[datetime(2001, 1, 1, 10):]
        assert_frame_equal(result, df)
        result = df.ix[:datetime(2001, 1, 4, 10)]
        assert_frame_equal(result, df)
        result = df.ix[datetime(2001, 1, 1, 10):datetime(2001, 1, 4, 10)]
        assert_frame_equal(result, df)

        result = df.ix[datetime(2001, 1, 1, 11):]
        expected = df.iloc[1:]
        assert_frame_equal(result, expected)
        result = df.ix['20010101 11':]
        assert_frame_equal(result, expected)

        # duplicates
        df = pd.DataFrame(np.arange(5., dtype='float64'),
                          index=[datetime(2001, 1, i, 10, 00)
                                 for i in [1, 2, 2, 3, 4]])

        result = df.ix[datetime(2001, 1, 1, 10):]
        assert_frame_equal(result, df)
        result = df.ix[:datetime(2001, 1, 4, 10)]
        assert_frame_equal(result, df)
        result = df.ix[datetime(2001, 1, 1, 10):datetime(2001, 1, 4, 10)]
        assert_frame_equal(result, df)

        result = df.ix[datetime(2001, 1, 1, 11):]
        expected = df.iloc[1:]
        assert_frame_equal(result, expected)
        result = df.ix['20010101 11':]
        assert_frame_equal(result, expected)


class TestSeriesDatetime64(tm.TestCase):
    def setUp(self):
        self.series = Series(date_range('1/1/2000', periods=10))

    def test_auto_conversion(self):
        series = Series(list(date_range('1/1/2000', periods=10)))
        self.assertEqual(series.dtype, 'M8[ns]')

    def test_constructor_cant_cast_datetime64(self):
        msg = "Cannot cast datetime64 to "
        with tm.assertRaisesRegexp(TypeError, msg):
            Series(date_range('1/1/2000', periods=10), dtype=float)

        with tm.assertRaisesRegexp(TypeError, msg):
            Series(date_range('1/1/2000', periods=10), dtype=int)

    def test_constructor_cast_object(self):
        s = Series(date_range('1/1/2000', periods=10), dtype=object)
        exp = Series(date_range('1/1/2000', periods=10))
        tm.assert_series_equal(s, exp)

    def test_series_comparison_scalars(self):
        val = datetime(2000, 1, 4)
        result = self.series > val
        expected = Series([x > val for x in self.series])
        self.assert_series_equal(result, expected)

        val = self.series[5]
        result = self.series > val
        expected = Series([x > val for x in self.series])
        self.assert_series_equal(result, expected)

    def test_between(self):
        left, right = self.series[[2, 7]]

        result = self.series.between(left, right)
        expected = (self.series >= left) & (self.series <= right)
        assert_series_equal(result, expected)

    # ---------------------------------------------------------------------
    # NaT support

    def test_NaT_scalar(self):
        series = Series([0, 1000, 2000, iNaT], dtype='M8[ns]')

        val = series[3]
        self.assertTrue(com.isnull(val))

        series[2] = val
        self.assertTrue(com.isnull(series[2]))

    def test_NaT_cast(self):
        # GH10747
        result = Series([np.nan]).astype('M8[ns]')
        expected = Series([NaT])
        assert_series_equal(result, expected)

    def test_set_none_nan(self):
        self.series[3] = None
        self.assertIs(self.series[3], NaT)

        self.series[3:5] = None
        self.assertIs(self.series[4], NaT)

        self.series[5] = np.nan
        self.assertIs(self.series[5], NaT)

        self.series[5:7] = np.nan
        self.assertIs(self.series[6], NaT)

    def test_intercept_astype_object(self):

        # this test no longer makes sense as series is by default already
        # M8[ns]
        expected = self.series.astype('object')

        df = DataFrame({'a': self.series,
                        'b': np.random.randn(len(self.series))})
        exp_dtypes = pd.Series([np.dtype('datetime64[ns]'),
                                np.dtype('float64')], index=['a', 'b'])
        tm.assert_series_equal(df.dtypes, exp_dtypes)

        result = df.values.squeeze()
        self.assertTrue((result[:, 0] == expected.values).all())

        df = DataFrame({'a': self.series, 'b': ['foo'] * len(self.series)})

        result = df.values.squeeze()
        self.assertTrue((result[:, 0] == expected.values).all())

    def test_nat_operations(self):
        # GH 8617
        s = Series([0, pd.NaT], dtype='m8[ns]')
        exp = s[0]
        self.assertEqual(s.median(), exp)
        self.assertEqual(s.min(), exp)
        self.assertEqual(s.max(), exp)


class TestTimestamp(tm.TestCase):
    def test_class_ops_pytz(self):
        tm._skip_if_no_pytz()
        from pytz import timezone

        def compare(x, y):
            self.assertEqual(int(Timestamp(x).value / 1e9),
                             int(Timestamp(y).value / 1e9))

        compare(Timestamp.now(), datetime.now())
        compare(Timestamp.now('UTC'), datetime.now(timezone('UTC')))
        compare(Timestamp.utcnow(), datetime.utcnow())
        compare(Timestamp.today(), datetime.today())
        current_time = calendar.timegm(datetime.now().utctimetuple())
        compare(Timestamp.utcfromtimestamp(current_time),
                datetime.utcfromtimestamp(current_time))
        compare(Timestamp.fromtimestamp(current_time),
                datetime.fromtimestamp(current_time))

        date_component = datetime.utcnow()
        time_component = (date_component + timedelta(minutes=10)).time()
        compare(Timestamp.combine(date_component, time_component),
                datetime.combine(date_component, time_component))

    def test_class_ops_dateutil(self):
        tm._skip_if_no_dateutil()
        from dateutil.tz import tzutc

        def compare(x, y):
            self.assertEqual(int(np.round(Timestamp(x).value / 1e9)),
                             int(np.round(Timestamp(y).value / 1e9)))

        compare(Timestamp.now(), datetime.now())
        compare(Timestamp.now('UTC'), datetime.now(tzutc()))
        compare(Timestamp.utcnow(), datetime.utcnow())
        compare(Timestamp.today(), datetime.today())
        current_time = calendar.timegm(datetime.now().utctimetuple())
        compare(Timestamp.utcfromtimestamp(current_time),
                datetime.utcfromtimestamp(current_time))
        compare(Timestamp.fromtimestamp(current_time),
                datetime.fromtimestamp(current_time))

        date_component = datetime.utcnow()
        time_component = (date_component + timedelta(minutes=10)).time()
        compare(Timestamp.combine(date_component, time_component),
                datetime.combine(date_component, time_component))

    def test_basics_nanos(self):
        val = np.int64(946684800000000000).view('M8[ns]')
        stamp = Timestamp(val.view('i8') + 500)
        self.assertEqual(stamp.year, 2000)
        self.assertEqual(stamp.month, 1)
        self.assertEqual(stamp.microsecond, 0)
        self.assertEqual(stamp.nanosecond, 500)

        # GH 14415
        val = np.iinfo(np.int64).min + 80000000000000
        stamp = Timestamp(val)
        self.assertEqual(stamp.year, 1677)
        self.assertEqual(stamp.month, 9)
        self.assertEqual(stamp.day, 21)
        self.assertEqual(stamp.microsecond, 145224)
        self.assertEqual(stamp.nanosecond, 192)

    def test_unit(self):

        def check(val, unit=None, h=1, s=1, us=0):
            stamp = Timestamp(val, unit=unit)
            self.assertEqual(stamp.year, 2000)
            self.assertEqual(stamp.month, 1)
            self.assertEqual(stamp.day, 1)
            self.assertEqual(stamp.hour, h)
            if unit != 'D':
                self.assertEqual(stamp.minute, 1)
                self.assertEqual(stamp.second, s)
                self.assertEqual(stamp.microsecond, us)
            else:
                self.assertEqual(stamp.minute, 0)
                self.assertEqual(stamp.second, 0)
                self.assertEqual(stamp.microsecond, 0)
            self.assertEqual(stamp.nanosecond, 0)

        ts = Timestamp('20000101 01:01:01')
        val = ts.value
        days = (ts - Timestamp('1970-01-01')).days

        check(val)
        check(val / long(1000), unit='us')
        check(val / long(1000000), unit='ms')
        check(val / long(1000000000), unit='s')
        check(days, unit='D', h=0)

        # using truediv, so these are like floats
        if compat.PY3:
            check((val + 500000) / long(1000000000), unit='s', us=500)
            check((val + 500000000) / long(1000000000), unit='s', us=500000)
            check((val + 500000) / long(1000000), unit='ms', us=500)

        # get chopped in py2
        else:
            check((val + 500000) / long(1000000000), unit='s')
            check((val + 500000000) / long(1000000000), unit='s')
            check((val + 500000) / long(1000000), unit='ms')

        # ok
        check((val + 500000) / long(1000), unit='us', us=500)
        check((val + 500000000) / long(1000000), unit='ms', us=500000)

        # floats
        check(val / 1000.0 + 5, unit='us', us=5)
        check(val / 1000.0 + 5000, unit='us', us=5000)
        check(val / 1000000.0 + 0.5, unit='ms', us=500)
        check(val / 1000000.0 + 0.005, unit='ms', us=5)
        check(val / 1000000000.0 + 0.5, unit='s', us=500000)
        check(days + 0.5, unit='D', h=12)

        # nan
        result = Timestamp(np.nan)
        self.assertIs(result, NaT)

        result = Timestamp(None)
        self.assertIs(result, NaT)

        result = Timestamp(iNaT)
        self.assertIs(result, NaT)

        result = Timestamp(NaT)
        self.assertIs(result, NaT)

        result = Timestamp('NaT')
        self.assertIs(result, NaT)

        self.assertTrue(isnull(Timestamp('nat')))

    def test_roundtrip(self):

        # test value to string and back conversions
        # further test accessors
        base = Timestamp('20140101 00:00:00')

        result = Timestamp(base.value + pd.Timedelta('5ms').value)
        self.assertEqual(result, Timestamp(str(base) + ".005000"))
        self.assertEqual(result.microsecond, 5000)

        result = Timestamp(base.value + pd.Timedelta('5us').value)
        self.assertEqual(result, Timestamp(str(base) + ".000005"))
        self.assertEqual(result.microsecond, 5)

        result = Timestamp(base.value + pd.Timedelta('5ns').value)
        self.assertEqual(result, Timestamp(str(base) + ".000000005"))
        self.assertEqual(result.nanosecond, 5)
        self.assertEqual(result.microsecond, 0)

        result = Timestamp(base.value + pd.Timedelta('6ms 5us').value)
        self.assertEqual(result, Timestamp(str(base) + ".006005"))
        self.assertEqual(result.microsecond, 5 + 6 * 1000)

        result = Timestamp(base.value + pd.Timedelta('200ms 5us').value)
        self.assertEqual(result, Timestamp(str(base) + ".200005"))
        self.assertEqual(result.microsecond, 5 + 200 * 1000)

    def test_comparison(self):
        # 5-18-2012 00:00:00.000
        stamp = long(1337299200000000000)

        val = Timestamp(stamp)

        self.assertEqual(val, val)
        self.assertFalse(val != val)
        self.assertFalse(val < val)
        self.assertTrue(val <= val)
        self.assertFalse(val > val)
        self.assertTrue(val >= val)

        other = datetime(2012, 5, 18)
        self.assertEqual(val, other)
        self.assertFalse(val != other)
        self.assertFalse(val < other)
        self.assertTrue(val <= other)
        self.assertFalse(val > other)
        self.assertTrue(val >= other)

        other = Timestamp(stamp + 100)

        self.assertNotEqual(val, other)
        self.assertNotEqual(val, other)
        self.assertTrue(val < other)
        self.assertTrue(val <= other)
        self.assertTrue(other > val)
        self.assertTrue(other >= val)

    def test_compare_invalid(self):

        # GH 8058
        val = Timestamp('20130101 12:01:02')
        self.assertFalse(val == 'foo')
        self.assertFalse(val == 10.0)
        self.assertFalse(val == 1)
        self.assertFalse(val == long(1))
        self.assertFalse(val == [])
        self.assertFalse(val == {'foo': 1})
        self.assertFalse(val == np.float64(1))
        self.assertFalse(val == np.int64(1))

        self.assertTrue(val != 'foo')
        self.assertTrue(val != 10.0)
        self.assertTrue(val != 1)
        self.assertTrue(val != long(1))
        self.assertTrue(val != [])
        self.assertTrue(val != {'foo': 1})
        self.assertTrue(val != np.float64(1))
        self.assertTrue(val != np.int64(1))

        # ops testing
        df = DataFrame(randn(5, 2))
        a = df[0]
        b = Series(randn(5))
        b.name = Timestamp('2000-01-01')
        tm.assert_series_equal(a / b, 1 / (b / a))

    def test_cant_compare_tz_naive_w_aware(self):
        tm._skip_if_no_pytz()
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

    def test_cant_compare_tz_naive_w_aware_explicit_pytz(self):
        tm._skip_if_no_pytz()
        from pytz import utc
        # #1404
        a = Timestamp('3/12/2012')
        b = Timestamp('3/12/2012', tz=utc)

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

    def test_cant_compare_tz_naive_w_aware_dateutil(self):
        tm._skip_if_no_dateutil()
        from dateutil.tz import tzutc
        utc = tzutc()
        # #1404
        a = Timestamp('3/12/2012')
        b = Timestamp('3/12/2012', tz=utc)

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
        self.assertEqual(result.nanosecond, val.nanosecond)

    def test_frequency_misc(self):
        self.assertEqual(frequencies.get_freq_group('T'),
                         frequencies.FreqGroup.FR_MIN)

        code, stride = frequencies.get_freq_code(offsets.Hour())
        self.assertEqual(code, frequencies.FreqGroup.FR_HR)

        code, stride = frequencies.get_freq_code((5, 'T'))
        self.assertEqual(code, frequencies.FreqGroup.FR_MIN)
        self.assertEqual(stride, 5)

        offset = offsets.Hour()
        result = frequencies.to_offset(offset)
        self.assertEqual(result, offset)

        result = frequencies.to_offset((5, 'T'))
        expected = offsets.Minute(5)
        self.assertEqual(result, expected)

        self.assertRaises(ValueError, frequencies.get_freq_code, (5, 'baz'))

        self.assertRaises(ValueError, frequencies.to_offset, '100foo')

        self.assertRaises(ValueError, frequencies.to_offset, ('', ''))

        with tm.assert_produces_warning(FutureWarning, check_stacklevel=False):
            result = frequencies.get_standard_freq(offsets.Hour())
        self.assertEqual(result, 'H')

    def test_hash_equivalent(self):
        d = {datetime(2011, 1, 1): 5}
        stamp = Timestamp(datetime(2011, 1, 1))
        self.assertEqual(d[stamp], 5)

    def test_timestamp_compare_scalars(self):
        # case where ndim == 0
        lhs = np.datetime64(datetime(2013, 12, 6))
        rhs = Timestamp('now')
        nat = Timestamp('nat')

        ops = {'gt': 'lt',
               'lt': 'gt',
               'ge': 'le',
               'le': 'ge',
               'eq': 'eq',
               'ne': 'ne'}

        for left, right in ops.items():
            left_f = getattr(operator, left)
            right_f = getattr(operator, right)
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

    def test_is_leap_year(self):
        # GH 13727
        for tz in [None, 'UTC', 'US/Eastern', 'Asia/Tokyo']:
            dt = Timestamp('2000-01-01 00:00:00', tz=tz)
            self.assertTrue(dt.is_leap_year)
            self.assertIsInstance(dt.is_leap_year, bool)

            dt = Timestamp('1999-01-01 00:00:00', tz=tz)
            self.assertFalse(dt.is_leap_year)

            dt = Timestamp('2004-01-01 00:00:00', tz=tz)
            self.assertTrue(dt.is_leap_year)

            dt = Timestamp('2100-01-01 00:00:00', tz=tz)
            self.assertFalse(dt.is_leap_year)

        self.assertFalse(pd.NaT.is_leap_year)
        self.assertIsInstance(pd.NaT.is_leap_year, bool)


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
        self.assertEqual(result, expected)

    def test_slice_quarter(self):
        dti = DatetimeIndex(freq='D', start=datetime(2000, 6, 1), periods=500)

        s = Series(np.arange(len(dti)), index=dti)
        self.assertEqual(len(s['2001Q1']), 90)

        df = DataFrame(np.random.rand(len(dti), 5), index=dti)
        self.assertEqual(len(df.ix['1Q01']), 90)

    def test_slice_month(self):
        dti = DatetimeIndex(freq='D', start=datetime(2005, 1, 1), periods=500)
        s = Series(np.arange(len(dti)), index=dti)
        self.assertEqual(len(s['2005-11']), 30)

        df = DataFrame(np.random.rand(len(dti), 5), index=dti)
        self.assertEqual(len(df.ix['2005-11']), 30)

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
        self.assertEqual(result, s.iloc[0])

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

        self.assertEqual(s['2005-1-1 20:00'], s.ix[0])
        self.assertRaises(Exception, s.__getitem__, '2004-12-31 00:15')

    def test_partial_slice_minutely(self):
        rng = DatetimeIndex(freq='S', start=datetime(2005, 1, 1, 23, 59, 0),
                            periods=500)
        s = Series(np.arange(len(rng)), index=rng)

        result = s['2005-1-1 23:59']
        assert_series_equal(result, s.ix[:60])

        result = s['2005-1-1']
        assert_series_equal(result, s.ix[:60])

        self.assertEqual(s[Timestamp('2005-1-1 23:59:00')], s.ix[0])
        self.assertRaises(Exception, s.__getitem__, '2004-12-31 00:00:00')

    def test_partial_slice_second_precision(self):
        rng = DatetimeIndex(start=datetime(2005, 1, 1, 0, 0, 59,
                                           microsecond=999990),
                            periods=20, freq='US')
        s = Series(np.arange(20), rng)

        assert_series_equal(s['2005-1-1 00:00'], s.iloc[:10])
        assert_series_equal(s['2005-1-1 00:00:59'], s.iloc[:10])

        assert_series_equal(s['2005-1-1 00:01'], s.iloc[10:])
        assert_series_equal(s['2005-1-1 00:01:00'], s.iloc[10:])

        self.assertEqual(s[Timestamp('2005-1-1 00:00:59.999990')], s.iloc[0])
        self.assertRaisesRegexp(KeyError, '2005-1-1 00:00:00',
                                lambda: s['2005-1-1 00:00:00'])

    def test_partial_slicing_with_multiindex(self):

        # GH 4758
        # partial string indexing with a multi-index buggy
        df = DataFrame({'ACCOUNT': ["ACCT1", "ACCT1", "ACCT1", "ACCT2"],
                        'TICKER': ["ABC", "MNP", "XYZ", "XYZ"],
                        'val': [1, 2, 3, 4]},
                       index=date_range("2013-06-19 09:30:00",
                                        periods=4, freq='5T'))
        df_multi = df.set_index(['ACCOUNT', 'TICKER'], append=True)

        expected = DataFrame([
            [1]
        ], index=Index(['ABC'], name='TICKER'), columns=['val'])
        result = df_multi.loc[('2013-06-19 09:30:00', 'ACCT1')]
        assert_frame_equal(result, expected)

        expected = df_multi.loc[
            (pd.Timestamp('2013-06-19 09:30:00', tz=None), 'ACCT1', 'ABC')]
        result = df_multi.loc[('2013-06-19 09:30:00', 'ACCT1', 'ABC')]
        assert_series_equal(result, expected)

        # this is a KeyError as we don't do partial string selection on
        # multi-levels
        def f():
            df_multi.loc[('2013-06-19', 'ACCT1', 'ABC')]

        self.assertRaises(KeyError, f)

        # GH 4294
        # partial slice on a series mi
        s = pd.DataFrame(randn(1000, 1000), index=pd.date_range(
            '2000-1-1', periods=1000)).stack()

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
        values = DatetimeIndex([snap + i * offset for i in range(n)])

        tm.assert_index_equal(rng, values)

        rng = date_range('1/1/2000 08:15', periods=n, normalize=False,
                         freq='B')
        the_time = time(8, 15)
        for val in rng:
            self.assertEqual(val.time(), the_time)

    def test_timedelta(self):
        # this is valid too
        index = date_range('1/1/2000', periods=50, freq='B')
        shifted = index + timedelta(1)
        back = shifted + timedelta(-1)
        self.assertTrue(tm.equalContents(index, back))
        self.assertEqual(shifted.freq, index.freq)
        self.assertEqual(shifted.freq, back.freq)

        result = index - timedelta(1)
        expected = index + timedelta(-1)
        tm.assert_index_equal(result, expected)

        # GH4134, buggy with timedeltas
        rng = date_range('2013', '2014')
        s = Series(rng)
        result1 = rng - pd.offsets.Hour(1)
        result2 = DatetimeIndex(s - np.timedelta64(100000000))
        result3 = rng - np.timedelta64(100000000)
        result4 = DatetimeIndex(s - pd.offsets.Hour(1))
        tm.assert_index_equal(result1, result4)
        tm.assert_index_equal(result2, result3)

    def test_shift(self):
        ts = Series(np.random.randn(5),
                    index=date_range('1/1/2000', periods=5, freq='H'))

        result = ts.shift(1, freq='5T')
        exp_index = ts.index.shift(1, freq='5T')
        tm.assert_index_equal(result.index, exp_index)

        # GH #1063, multiple of same base
        result = ts.shift(1, freq='4H')
        exp_index = ts.index + offsets.Hour(4)
        tm.assert_index_equal(result.index, exp_index)

        idx = DatetimeIndex(['2000-01-01', '2000-01-02', '2000-01-04'])
        self.assertRaises(ValueError, idx.shift, 1)

    def test_setops_preserve_freq(self):
        for tz in [None, 'Asia/Tokyo', 'US/Eastern']:
            rng = date_range('1/1/2000', '1/1/2002', name='idx', tz=tz)

            result = rng[:50].union(rng[50:100])
            self.assertEqual(result.name, rng.name)
            self.assertEqual(result.freq, rng.freq)
            self.assertEqual(result.tz, rng.tz)

            result = rng[:50].union(rng[30:100])
            self.assertEqual(result.name, rng.name)
            self.assertEqual(result.freq, rng.freq)
            self.assertEqual(result.tz, rng.tz)

            result = rng[:50].union(rng[60:100])
            self.assertEqual(result.name, rng.name)
            self.assertIsNone(result.freq)
            self.assertEqual(result.tz, rng.tz)

            result = rng[:50].intersection(rng[25:75])
            self.assertEqual(result.name, rng.name)
            self.assertEqual(result.freqstr, 'D')
            self.assertEqual(result.tz, rng.tz)

            nofreq = DatetimeIndex(list(rng[25:75]), name='other')
            result = rng[:50].union(nofreq)
            self.assertIsNone(result.name)
            self.assertEqual(result.freq, rng.freq)
            self.assertEqual(result.tz, rng.tz)

            result = rng[:50].intersection(nofreq)
            self.assertIsNone(result.name)
            self.assertEqual(result.freq, rng.freq)
            self.assertEqual(result.tz, rng.tz)

    def test_min_max(self):
        rng = date_range('1/1/2000', '12/31/2000')
        rng2 = rng.take(np.random.permutation(len(rng)))

        the_min = rng2.min()
        the_max = rng2.max()
        tm.assertIsInstance(the_min, Timestamp)
        tm.assertIsInstance(the_max, Timestamp)
        self.assertEqual(the_min, rng[0])
        self.assertEqual(the_max, rng[-1])

        self.assertEqual(rng.min(), rng[0])
        self.assertEqual(rng.max(), rng[-1])

    def test_min_max_series(self):
        rng = date_range('1/1/2000', periods=10, freq='4h')
        lvls = ['A', 'A', 'A', 'B', 'B', 'B', 'C', 'C', 'C', 'C']
        df = DataFrame({'TS': rng, 'V': np.random.randn(len(rng)), 'L': lvls})

        result = df.TS.max()
        exp = Timestamp(df.TS.iat[-1])
        self.assertTrue(isinstance(result, Timestamp))
        self.assertEqual(result, exp)

        result = df.TS.min()
        exp = Timestamp(df.TS.iat[0])
        self.assertTrue(isinstance(result, Timestamp))
        self.assertEqual(result, exp)

    def test_from_M8_structured(self):
        dates = [(datetime(2012, 9, 9, 0, 0), datetime(2012, 9, 8, 15, 10))]
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
        labels = [[0, 0, 1, 1, 2, 2, 3, 3], [0, 1, 0, 1, 0, 1, 0, 1]]

        index = MultiIndex(levels=levels, labels=labels)

        self.assertTrue(isinstance(index.get_level_values(0)[0], Timestamp))

    def test_frame_apply_dont_convert_datetime64(self):
        from pandas.tseries.offsets import BDay
        df = DataFrame({'x1': [datetime(1996, 1, 1)]})

        df = df.applymap(lambda x: x + BDay())
        df = df.applymap(lambda x: x + BDay())

        self.assertTrue(df.x1.dtype == 'M8[ns]')

    def test_date_range_fy5252(self):
        dr = date_range(start="2013-01-01", periods=2, freq=offsets.FY5253(
            startingMonth=1, weekday=3, variation="nearest"))
        self.assertEqual(dr[0], Timestamp('2013-01-31'))
        self.assertEqual(dr[1], Timestamp('2014-01-30'))

    def test_partial_slice_doesnt_require_monotonicity(self):
        # For historical reasons.
        s = pd.Series(np.arange(10), pd.date_range('2014-01-01', periods=10))

        nonmonotonic = s[[3, 5, 4]]
        expected = nonmonotonic.iloc[:0]
        timestamp = pd.Timestamp('2014-01-10')

        assert_series_equal(nonmonotonic['2014-01-10':], expected)
        self.assertRaisesRegexp(KeyError,
                                r"Timestamp\('2014-01-10 00:00:00'\)",
                                lambda: nonmonotonic[timestamp:])

        assert_series_equal(nonmonotonic.ix['2014-01-10':], expected)
        self.assertRaisesRegexp(KeyError,
                                r"Timestamp\('2014-01-10 00:00:00'\)",
                                lambda: nonmonotonic.ix[timestamp:])


class TimeConversionFormats(tm.TestCase):
    def test_to_datetime_format(self):
        values = ['1/1/2000', '1/2/2000', '1/3/2000']

        results1 = [Timestamp('20000101'), Timestamp('20000201'),
                    Timestamp('20000301')]
        results2 = [Timestamp('20000101'), Timestamp('20000102'),
                    Timestamp('20000103')]
        for vals, expecteds in [(values, (Index(results1), Index(results2))),
                                (Series(values),
                                 (Series(results1), Series(results2))),
                                (values[0], (results1[0], results2[0])),
                                (values[1], (results1[1], results2[1])),
                                (values[2], (results1[2], results2[2]))]:

            for i, fmt in enumerate(['%d/%m/%Y', '%m/%d/%Y']):
                result = to_datetime(vals, format=fmt)
                expected = expecteds[i]

                if isinstance(expected, Series):
                    assert_series_equal(result, Series(expected))
                elif isinstance(expected, Timestamp):
                    self.assertEqual(result, expected)
                else:
                    tm.assert_index_equal(result, expected)

    def test_to_datetime_format_YYYYMMDD(self):
        s = Series([19801222, 19801222] + [19810105] * 5)
        expected = Series([Timestamp(x) for x in s.apply(str)])

        result = to_datetime(s, format='%Y%m%d')
        assert_series_equal(result, expected)

        result = to_datetime(s.apply(str), format='%Y%m%d')
        assert_series_equal(result, expected)

        # with NaT
        expected = Series([Timestamp("19801222"), Timestamp("19801222")] +
                          [Timestamp("19810105")] * 5)
        expected[2] = np.nan
        s[2] = np.nan

        result = to_datetime(s, format='%Y%m%d')
        assert_series_equal(result, expected)

        # string with NaT
        s = s.apply(str)
        s[2] = 'nat'
        result = to_datetime(s, format='%Y%m%d')
        assert_series_equal(result, expected)

        # coercion
        # GH 7930
        s = Series([20121231, 20141231, 99991231])
        result = pd.to_datetime(s, format='%Y%m%d', errors='ignore')
        expected = Series([datetime(2012, 12, 31),
                           datetime(2014, 12, 31), datetime(9999, 12, 31)],
                          dtype=object)
        self.assert_series_equal(result, expected)

        result = pd.to_datetime(s, format='%Y%m%d', errors='coerce')
        expected = Series(['20121231', '20141231', 'NaT'], dtype='M8[ns]')
        assert_series_equal(result, expected)

    # GH 10178
    def test_to_datetime_format_integer(self):
        s = Series([2000, 2001, 2002])
        expected = Series([Timestamp(x) for x in s.apply(str)])

        result = to_datetime(s, format='%Y')
        assert_series_equal(result, expected)

        s = Series([200001, 200105, 200206])
        expected = Series([Timestamp(x[:4] + '-' + x[4:]) for x in s.apply(str)
                           ])

        result = to_datetime(s, format='%Y%m')
        assert_series_equal(result, expected)

    def test_to_datetime_format_microsecond(self):
        val = '01-Apr-2011 00:00:01.978'
        format = '%d-%b-%Y %H:%M:%S.%f'
        result = to_datetime(val, format=format)
        exp = datetime.strptime(val, format)
        self.assertEqual(result, exp)

    def test_to_datetime_format_time(self):
        data = [
            ['01/10/2010 15:20', '%m/%d/%Y %H:%M',
             Timestamp('2010-01-10 15:20')],
            ['01/10/2010 05:43', '%m/%d/%Y %I:%M',
             Timestamp('2010-01-10 05:43')],
            ['01/10/2010 13:56:01', '%m/%d/%Y %H:%M:%S',
             Timestamp('2010-01-10 13:56:01')]  # ,
            # ['01/10/2010 08:14 PM', '%m/%d/%Y %I:%M %p',
            #  Timestamp('2010-01-10 20:14')],
            # ['01/10/2010 07:40 AM', '%m/%d/%Y %I:%M %p',
            #  Timestamp('2010-01-10 07:40')],
            # ['01/10/2010 09:12:56 AM', '%m/%d/%Y %I:%M:%S %p',
            #  Timestamp('2010-01-10 09:12:56')]
        ]
        for s, format, dt in data:
            self.assertEqual(to_datetime(s, format=format), dt)

    def test_to_datetime_with_non_exact(self):
        # GH 10834
        _skip_if_has_locale()

        # 8904
        # exact kw
        if sys.version_info < (2, 7):
            raise nose.SkipTest('on python version < 2.7')

        s = Series(['19MAY11', 'foobar19MAY11', '19MAY11:00:00:00',
                    '19MAY11 00:00:00Z'])
        result = to_datetime(s, format='%d%b%y', exact=False)
        expected = to_datetime(s.str.extract(r'(\d+\w+\d+)', expand=False),
                               format='%d%b%y')
        assert_series_equal(result, expected)

    def test_parse_nanoseconds_with_formula(self):

        # GH8989
        # trunctaing the nanoseconds when a format was provided
        for v in ["2012-01-01 09:00:00.000000001",
                  "2012-01-01 09:00:00.000001",
                  "2012-01-01 09:00:00.001",
                  "2012-01-01 09:00:00.001000",
                  "2012-01-01 09:00:00.001000000", ]:
            expected = pd.to_datetime(v)
            result = pd.to_datetime(v, format="%Y-%m-%d %H:%M:%S.%f")
            self.assertEqual(result, expected)

    def test_to_datetime_format_weeks(self):
        data = [
            ['2009324', '%Y%W%w', Timestamp('2009-08-13')],
            ['2013020', '%Y%U%w', Timestamp('2013-01-13')]
        ]
        for s, format, dt in data:
            self.assertEqual(to_datetime(s, format=format), dt)


class TestToDatetimeInferFormat(tm.TestCase):

    def test_to_datetime_infer_datetime_format_consistent_format(self):
        s = pd.Series(pd.date_range('20000101', periods=50, freq='H'))

        test_formats = ['%m-%d-%Y', '%m/%d/%Y %H:%M:%S.%f',
                        '%Y-%m-%dT%H:%M:%S.%f']

        for test_format in test_formats:
            s_as_dt_strings = s.apply(lambda x: x.strftime(test_format))

            with_format = pd.to_datetime(s_as_dt_strings, format=test_format)
            no_infer = pd.to_datetime(s_as_dt_strings,
                                      infer_datetime_format=False)
            yes_infer = pd.to_datetime(s_as_dt_strings,
                                       infer_datetime_format=True)

            # Whether the format is explicitly passed, it is inferred, or
            # it is not inferred, the results should all be the same
            self.assert_series_equal(with_format, no_infer)
            self.assert_series_equal(no_infer, yes_infer)

    def test_to_datetime_infer_datetime_format_inconsistent_format(self):
        s = pd.Series(np.array(['01/01/2011 00:00:00',
                                '01-02-2011 00:00:00',
                                '2011-01-03T00:00:00']))

        # When the format is inconsistent, infer_datetime_format should just
        # fallback to the default parsing
        tm.assert_series_equal(pd.to_datetime(s, infer_datetime_format=False),
                               pd.to_datetime(s, infer_datetime_format=True))

        s = pd.Series(np.array(['Jan/01/2011', 'Feb/01/2011', 'Mar/01/2011']))

        tm.assert_series_equal(pd.to_datetime(s, infer_datetime_format=False),
                               pd.to_datetime(s, infer_datetime_format=True))

    def test_to_datetime_infer_datetime_format_series_with_nans(self):
        s = pd.Series(np.array(['01/01/2011 00:00:00', np.nan,
                                '01/03/2011 00:00:00', np.nan]))
        tm.assert_series_equal(pd.to_datetime(s, infer_datetime_format=False),
                               pd.to_datetime(s, infer_datetime_format=True))

    def test_to_datetime_infer_datetime_format_series_starting_with_nans(self):
        s = pd.Series(np.array([np.nan, np.nan, '01/01/2011 00:00:00',
                                '01/02/2011 00:00:00', '01/03/2011 00:00:00']))

        tm.assert_series_equal(pd.to_datetime(s, infer_datetime_format=False),
                               pd.to_datetime(s, infer_datetime_format=True))

    def test_to_datetime_iso8601_noleading_0s(self):
        # GH 11871
        s = pd.Series(['2014-1-1', '2014-2-2', '2015-3-3'])
        expected = pd.Series([pd.Timestamp('2014-01-01'),
                              pd.Timestamp('2014-02-02'),
                              pd.Timestamp('2015-03-03')])
        tm.assert_series_equal(pd.to_datetime(s), expected)
        tm.assert_series_equal(pd.to_datetime(s, format='%Y-%m-%d'), expected)


class TestGuessDatetimeFormat(tm.TestCase):
    def test_guess_datetime_format_with_parseable_formats(self):
        dt_string_to_format = (('20111230', '%Y%m%d'),
                               ('2011-12-30', '%Y-%m-%d'),
                               ('30-12-2011', '%d-%m-%Y'),
                               ('2011-12-30 00:00:00', '%Y-%m-%d %H:%M:%S'),
                               ('2011-12-30T00:00:00', '%Y-%m-%dT%H:%M:%S'),
                               ('2011-12-30 00:00:00.000000',
                                '%Y-%m-%d %H:%M:%S.%f'), )

        for dt_string, dt_format in dt_string_to_format:
            self.assertEqual(
                tools._guess_datetime_format(dt_string),
                dt_format
            )

    def test_guess_datetime_format_with_dayfirst(self):
        ambiguous_string = '01/01/2011'
        self.assertEqual(
            tools._guess_datetime_format(ambiguous_string, dayfirst=True),
            '%d/%m/%Y'
        )
        self.assertEqual(
            tools._guess_datetime_format(ambiguous_string, dayfirst=False),
            '%m/%d/%Y'
        )

    def test_guess_datetime_format_with_locale_specific_formats(self):
        # The month names will vary depending on the locale, in which
        # case these wont be parsed properly (dateutil can't parse them)
        _skip_if_has_locale()

        dt_string_to_format = (('30/Dec/2011', '%d/%b/%Y'),
                               ('30/December/2011', '%d/%B/%Y'),
                               ('30/Dec/2011 00:00:00', '%d/%b/%Y %H:%M:%S'), )

        for dt_string, dt_format in dt_string_to_format:
            self.assertEqual(
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

    def test_guess_datetime_format_nopadding(self):
        # GH 11142
        dt_string_to_format = (('2011-1-1', '%Y-%m-%d'),
                               ('30-1-2011', '%d-%m-%Y'),
                               ('1/1/2011', '%m/%d/%Y'),
                               ('2011-1-1 00:00:00', '%Y-%m-%d %H:%M:%S'),
                               ('2011-1-1 0:0:0', '%Y-%m-%d %H:%M:%S'),
                               ('2011-1-3T00:00:0', '%Y-%m-%dT%H:%M:%S'))

        for dt_string, dt_format in dt_string_to_format:
            self.assertEqual(
                tools._guess_datetime_format(dt_string),
                dt_format
            )

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
            np.array(
                [np.nan, np.nan, np.nan], dtype='O'))
        self.assertTrue(format_for_string_of_nans is None)


class TestTimestampToJulianDate(tm.TestCase):
    def test_compare_1700(self):
        r = Timestamp('1700-06-23').to_julian_date()
        self.assertEqual(r, 2342145.5)

    def test_compare_2000(self):
        r = Timestamp('2000-04-12').to_julian_date()
        self.assertEqual(r, 2451646.5)

    def test_compare_2100(self):
        r = Timestamp('2100-08-12').to_julian_date()
        self.assertEqual(r, 2488292.5)

    def test_compare_hour01(self):
        r = Timestamp('2000-08-12T01:00:00').to_julian_date()
        self.assertEqual(r, 2451768.5416666666666666)

    def test_compare_hour13(self):
        r = Timestamp('2000-08-12T13:00:00').to_julian_date()
        self.assertEqual(r, 2451769.0416666666666666)


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


class TestDaysInMonth(tm.TestCase):
    def test_coerce_deprecation(self):

        # deprecation of coerce
        with tm.assert_produces_warning(FutureWarning):
            to_datetime('2015-02-29', coerce=True)
        with tm.assert_produces_warning(FutureWarning):
            self.assertRaises(ValueError,
                              lambda: to_datetime('2015-02-29', coerce=False))

        # multiple arguments
        for e, c in zip(['raise', 'ignore', 'coerce'], [True, False]):
            with tm.assert_produces_warning(FutureWarning):
                self.assertRaises(TypeError,
                                  lambda: to_datetime('2015-02-29', errors=e,
                                                      coerce=c))

    # tests for issue #10154
    def test_day_not_in_month_coerce(self):
        self.assertTrue(isnull(to_datetime('2015-02-29', errors='coerce')))
        self.assertTrue(isnull(to_datetime('2015-02-29', format="%Y-%m-%d",
                                           errors='coerce')))
        self.assertTrue(isnull(to_datetime('2015-02-32', format="%Y-%m-%d",
                                           errors='coerce')))
        self.assertTrue(isnull(to_datetime('2015-04-31', format="%Y-%m-%d",
                                           errors='coerce')))

    def test_day_not_in_month_raise(self):
        self.assertRaises(ValueError, to_datetime, '2015-02-29',
                          errors='raise')
        self.assertRaises(ValueError, to_datetime, '2015-02-29',
                          errors='raise', format="%Y-%m-%d")
        self.assertRaises(ValueError, to_datetime, '2015-02-32',
                          errors='raise', format="%Y-%m-%d")
        self.assertRaises(ValueError, to_datetime, '2015-04-31',
                          errors='raise', format="%Y-%m-%d")

    def test_day_not_in_month_ignore(self):
        self.assertEqual(to_datetime(
            '2015-02-29', errors='ignore'), '2015-02-29')
        self.assertEqual(to_datetime(
            '2015-02-29', errors='ignore', format="%Y-%m-%d"), '2015-02-29')
        self.assertEqual(to_datetime(
            '2015-02-32', errors='ignore', format="%Y-%m-%d"), '2015-02-32')
        self.assertEqual(to_datetime(
            '2015-04-31', errors='ignore', format="%Y-%m-%d"), '2015-04-31')


if __name__ == '__main__':
    nose.runmodule(argv=[__file__, '-vvs', '-x', '--pdb', '--pdb-failure'],
                   exit=False)
