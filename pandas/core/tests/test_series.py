# pylint: disable-msg=E1101,W0612

from datetime import datetime, timedelta
import os
import operator
import unittest

import numpy as np

from pandas.core.api import (Index, Series, TimeSeries, DataFrame, isnull)
import pandas.core.datetools as datetools
from pandas.util.testing import assert_series_equal
import pandas.util.testing as common

#-------------------------------------------------------------------------------
# Series test cases

class TestSeries(unittest.TestCase):
    def setUp(self):
        self.ts = common.makeTimeSeries()
        self.series = common.makeStringSeries()
        self.objSeries = common.makeObjectSeries()

        self.empty = Series([], index=[])

    def test_constructor(self):
        # Recognize TimeSeries
        self.assert_(isinstance(self.ts, TimeSeries))

        # Pass in Series
        derived = Series(self.ts)
        self.assert_(isinstance(derived, TimeSeries))

        self.assert_(common.equalContents(derived.index, self.ts.index))
        # Ensure new index is not created
        self.assertEquals(id(self.ts.index), id(derived.index))

        # Pass in scalar
        scalar = Series(0.5)
        self.assert_(isinstance(scalar, float))

        # Mixed type Series
        mixed = Series(['hello', np.NaN], index=[0, 1])
        self.assert_(mixed.dtype == np.object_)
        self.assert_(mixed[1] is np.NaN)

        self.assertRaises(Exception, Series, [0, 1, 2], index=None)

        self.assert_(not isinstance(self.empty, TimeSeries))
        self.assert_(not isinstance(Series({}), TimeSeries))

        self.assertRaises(Exception, Series, np.random.randn(3, 3),
                          index=np.arange(3))

    def test_constructor_corner(self):
        df = common.makeTimeDataFrame()
        objs = [df, df]
        s = Series(objs, index=[0, 1])
        self.assert_(isinstance(s, Series))

    def test_fromDict(self):
        data = {'a' : 0, 'b' : 1, 'c' : 2, 'd' : 3}

        series = Series(data)
        self.assert_(common.is_sorted(series.index))

        data = {'a' : 0, 'b' : '1', 'c' : '2', 'd' : datetime.now()}
        series = Series(data)
        self.assert_(series.dtype == np.object_)

        data = {'a' : 0, 'b' : '1', 'c' : '2', 'd' : '3'}
        series = Series(data)
        self.assert_(series.dtype == np.object_)

        data = {'a' : '0', 'b' : '1'}
        series = Series(data, dtype=float)
        self.assert_(series.dtype == np.float64)

    def test_setindex(self):
        # wrong type
        series = self.series.copy()
        self.assertRaises(TypeError, series._set_index, None)

        # wrong length
        series = self.series.copy()
        self.assertRaises(AssertionError, series._set_index,
                          np.arange(len(series) - 1))

        # works
        series = self.series.copy()
        series.index = np.arange(len(series))
        self.assert_(isinstance(series.index, Index))

    def test_array_finalize(self):
        pass

    def test_fromValue(self):
        nans = Series.fromValue(np.NaN, index=self.ts.index)
        self.assert_(nans.dtype == np.float_)

        strings = Series.fromValue('foo', index=self.ts.index)
        self.assert_(strings.dtype == np.object_)

        d = datetime.now()
        dates = Series.fromValue(d, index=self.ts.index)
        self.assert_(dates.dtype == np.object_)

    def test_contains(self):
        common.assert_contains_all(self.ts.index, self.ts)

    def test_save_load(self):
        self.series.save('tmp1')
        self.ts.save('tmp3')
        unp_series = Series.load('tmp1')
        unp_ts = Series.load('tmp3')
        os.remove('tmp1')
        os.remove('tmp3')
        assert_series_equal(unp_series, self.series)
        assert_series_equal(unp_ts, self.ts)

    def test_getitem_get(self):
        idx1 = self.series.index[5]
        idx2 = self.objSeries.index[5]

        self.assertEqual(self.series[idx1], self.series.get(idx1))
        self.assertEqual(self.objSeries[idx2], self.objSeries.get(idx2))

        self.assertEqual(self.series[idx1], self.series[5])
        self.assertEqual(self.objSeries[idx2], self.objSeries[5])

        self.assert_(self.series.get(-1) is None)
        self.assertEqual(self.series[5], self.series.get(self.series.index[5]))

        # missing
        d = self.ts.index[0] - datetools.bday
        self.assertRaises(Exception, self.ts.__getitem__, d),

    def test_fancy(self):
        slice1 = self.series[[1,2,3]]
        slice2 = self.objSeries[[1,2,3]]
        self.assertEqual(self.series.index[2], slice1.index[1])
        self.assertEqual(self.objSeries.index[2], slice2.index[1])
        self.assertEqual(self.series[2], slice1[1])
        self.assertEqual(self.objSeries[2], slice2[1])

    def test_slice(self):
        numSlice = self.series[10:20]
        numSliceEnd = self.series[-10:]
        objSlice = self.objSeries[10:20]

        self.assert_(self.series.index[9] not in numSlice.index)
        self.assert_(self.objSeries.index[9] not in objSlice.index)

        self.assertEqual(len(numSlice), len(numSlice.index))
        self.assertEqual(self.series[numSlice.index[0]],
                         numSlice[numSlice.index[0]])

        self.assertEqual(numSlice.index[1], self.series.index[11])

        self.assert_(common.equalContents(numSliceEnd,
                                          np.array(self.series)[-10:]))

    def test_setitem(self):
        self.ts[self.ts.index[5]] = np.NaN
        self.ts[[1,2,17]] = np.NaN
        self.ts[6] = np.NaN
        self.assert_(np.isnan(self.ts[6]))
        self.assert_(np.isnan(self.ts[2]))
        self.ts[np.isnan(self.ts)] = 5
        self.assert_(not np.isnan(self.ts[2]))

        # caught this bug when writing tests
        series = Series(common.makeIntIndex(20).astype(float),
                        index=common.makeIntIndex(20))

        series[::2] = 0
        self.assert_((series[::2] == 0).all())

    def test_setslice(self):
        sl = self.ts[5:20]
        self.assertEqual(len(sl), len(sl.index))
        self.assertEqual(len(sl.index.indexMap), len(sl.index))

    def test_repr(self):
        str(self.ts)
        str(self.series)
        str(self.objSeries)

        str(Series(common.randn(1000), index=np.arange(1000)))

        # empty
        str(self.empty)

        # with NaNs
        self.series[5:7] = np.NaN
        str(self.series)

    def test_toString(self):
        from cStringIO import StringIO
        self.ts.toString(buffer=StringIO())

    def test_iter(self):
        for i, val in enumerate(self.series):
            self.assertEqual(val, self.series[i])

        for i, val in enumerate(self.ts):
            self.assertEqual(val, self.ts[i])

    def test_keys(self):
        self.assert_(self.ts.keys() is self.ts.index)

    def test_values(self):
        self.assert_(np.array_equal(self.ts, self.ts.values()))

    def test_iteritems(self):
        for idx, val in self.series.iteritems():
            self.assertEqual(val, self.series[idx])

        for idx, val in self.ts.iteritems():
            self.assertEqual(val, self.ts[idx])

    def test_stats(self):
        self.series[5:15] = np.NaN

        s1 = np.array(self.series)
        s1 = s1[-np.isnan(s1)]

        self.assertEquals(np.min(s1), self.series.min())
        self.assertEquals(np.max(s1), self.series.max())
        self.assertEquals(np.sum(s1), self.series.sum())
        self.assertEquals(np.mean(s1), self.series.mean())
        self.assertEquals(np.std(s1, ddof=1), self.series.std())
        self.assertEquals(np.var(s1, ddof=1), self.series.var())

        try:
            from scipy.stats import skew
            common.assert_almost_equal(skew(s1, bias=False),
                                       self.series.skew())
        except ImportError:
            pass

        self.assert_(not np.isnan(np.sum(self.series)))
        self.assert_(not np.isnan(np.mean(self.series)))
        self.assert_(not np.isnan(np.std(self.series)))
        self.assert_(not np.isnan(np.var(self.series)))
        self.assert_(not np.isnan(np.min(self.series)))
        self.assert_(not np.isnan(np.max(self.series)))

        self.assert_(np.isnan(Series([1.], index=[1]).std()))
        self.assert_(np.isnan(Series([1.], index=[1]).var()))
        self.assert_(np.isnan(Series([1.], index=[1]).skew()))

    def test_append(self):
        appendedSeries = self.series.append(self.ts)
        for idx, value in appendedSeries.iteritems():
            if idx in self.series.index:
                self.assertEqual(value, self.series[idx])
            elif idx in self.ts.index:
                self.assertEqual(value, self.ts[idx])
            else:
                self.fail("orphaned index!")

        self.assertRaises(Exception, self.ts.append, self.ts)

    def test_operators(self):
        series = self.ts
        other = self.ts[::2]

        def _check_op(other, op):
            cython_or_numpy = op(series, other)
            python = series.combineFunc(other, op)

            common.assert_almost_equal(cython_or_numpy, python)

        def check(other):
            _check_op(other, operator.add)
            _check_op(other, operator.sub)
            _check_op(other, operator.div)
            _check_op(other, operator.mul)
            _check_op(other, operator.pow)

            _check_op(other, lambda x, y: operator.add(y, x))
            _check_op(other, lambda x, y: operator.sub(y, x))
            _check_op(other, lambda x, y: operator.div(y, x))
            _check_op(other, lambda x, y: operator.mul(y, x))
            _check_op(other, lambda x, y: operator.pow(y, x))

        check(self.ts * 2)
        check(self.ts[::2])
        check(5)

        def check_comparators(other):
            _check_op(other, operator.gt)
            _check_op(other, operator.ge)
            _check_op(other, operator.eq)
            _check_op(other, operator.lt)
            _check_op(other, operator.le)

        check_comparators(5)
        check_comparators(self.ts + 1)

    def test_operators_date(self):
        result = self.objSeries + timedelta(1)
        result = self.objSeries - timedelta(1)

    def test_operators_corner(self):
        series = self.ts

        empty = Series([], index=Index([]))

        result = series + empty
        self.assert_(np.isnan(result).all())

        result = empty + Series([], index=Index([]))
        self.assert_(len(result) == 0)

        deltas = Series([timedelta(1)] * 5, index=np.arange(5))
        sub_deltas = deltas[::2]

        deltas5 = deltas * 5
        deltas = deltas + sub_deltas

    def test_operators_frame(self):
        # rpow does not work with DataFrame
        df = DataFrame({'A' : self.ts})

        common.assert_almost_equal(self.ts + self.ts, (self.ts + df)['A'])
        self.assertRaises(Exception, self.ts.__pow__, df)

    def test_combineFirst(self):
        series = Series(common.makeIntIndex(20).astype(float),
                        index=common.makeIntIndex(20))

        series_copy = series * 2
        series_copy[::2] = np.NaN

        # nothing used from the input
        combined = series.combineFirst(series_copy)

        self.assert_(np.array_equal(combined, series))

        # Holes filled from input
        combined = series_copy.combineFirst(series)
        self.assert_(np.isfinite(combined).all())

        self.assert_(np.array_equal(combined[::2], series[::2]))
        self.assert_(np.array_equal(combined[1::2], series_copy[1::2]))

        # mixed types
        index = common.makeStringIndex(20)
        floats = Series(common.randn(20), index=index)
        strings = Series(common.makeStringIndex(10), index=index[::2])

        combined = strings.combineFirst(floats)

        common.assert_dict_equal(strings, combined, compare_keys=False)
        common.assert_dict_equal(floats[1::2], combined, compare_keys=False)

        # corner case
        s = Series([1., 2, 3], index=[0, 1, 2])
        result = s.combineFirst(Series([], index=[]))
        assert_series_equal(s, result)

    def test_overloads(self):
        methods = ['argsort', 'cumsum', 'cumprod']

        for method in methods:
            func = getattr(np, method)

            self.assert_(np.array_equal(func(self.ts), func(np.array(self.ts))))

            # with missing values
            ts = self.ts.copy()
            ts[::2] = np.NaN

            result = func(ts)[1::2]
            expected = func(np.array(ts.valid()))

            self.assert_(np.array_equal(result, expected))

        argsorted = self.ts.argsort()
        self.assert_(argsorted.dtype == np.int_)

    def test_median(self):
        self.assertAlmostEqual(np.median(self.ts), self.ts.median())

        ts = self.ts.copy()
        ts[::2] = np.NaN

        self.assertAlmostEqual(np.median(ts.valid()), ts.median())

    def test_corr(self):
        # full overlap
        self.assertAlmostEqual(self.ts.corr(self.ts), 1)

        # partial overlap
        self.assertAlmostEqual(self.ts[:15].corr(self.ts[5:]), 1)

        # No overlap
        self.assert_(np.isnan(self.ts[::2].corr(self.ts[1::2])))

        # additional checks?

    def test_copy(self):
        ts = self.ts.copy()

        ts[::2] = np.NaN

        # Did not modify original Series
        self.assertFalse(np.isnan(self.ts[0]))

    def test_count(self):
        self.assertEqual(self.ts.count(), len(self.ts))

        self.ts[::2] = np.NaN

        self.assertEqual(self.ts.count(), np.isfinite(self.ts).sum())

    def test_sort(self):
        ts = self.ts.copy()
        ts.sort()

        self.assert_(np.array_equal(ts, self.ts.order()))
        self.assert_(np.array_equal(ts.index, self.ts.order().index))

    def test_order(self):

        ts = self.ts.copy()
        ts[:5] = np.NaN
        vals = ts.values()

        result = ts.order()
        self.assert_(np.isnan(result[-5:]).all())
        self.assert_(np.array_equal(result[:-5], np.sort(vals[5:])))

        result = ts.order(missingAtEnd=False)
        self.assert_(np.isnan(result[:5]).all())
        self.assert_(np.array_equal(result[5:], np.sort(vals[5:])))

    def test_map(self):
        result = self.ts.map(lambda x: x * 2)

        self.assert_(np.array_equal(result, self.ts * 2))

    def test_toCSV(self):
        self.ts.toCSV('_foo')
        os.remove('_foo')

    def test_toDict(self):
        self.assert_(np.array_equal(Series(self.ts.toDict()), self.ts))

    def test_cap(self):
        val = self.ts.median()

        self.assertEqual(self.ts.cap(val).max(), val)

    def test_floor(self):
        val = self.ts.median()

        self.assertEqual(self.ts.floor(val).min(), val)

    def test_valid(self):
        ts = self.ts.copy()
        ts[::2] = np.NaN

        result = ts.valid()
        self.assertEqual(len(result), ts.count())

        common.assert_dict_equal(result, ts, compare_keys=False)

    def test_shift(self):
        shifted = self.ts.shift(1)
        unshifted = shifted.shift(-1)

        common.assert_dict_equal(unshifted.valid(), self.ts, compare_keys=False)

        offset = datetools.bday
        shifted = self.ts.shift(1, offset=offset)
        unshifted = shifted.shift(-1, offset=offset)

        assert_series_equal(unshifted, self.ts)

        unshifted = self.ts.shift(0, offset=offset)
        assert_series_equal(unshifted, self.ts)

        shifted = self.ts.shift(1, timeRule='WEEKDAY')
        unshifted = shifted.shift(-1, timeRule='WEEKDAY')

        assert_series_equal(unshifted, self.ts)

        # corner case
        unshifted = self.ts.shift(0)
        assert_series_equal(unshifted, self.ts)

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

        truncated = ts.truncate(before=self.ts.index[-1] + offset,
                                after=self.ts.index[0] - offset)
        assert(len(truncated) == 0)

    def test_asOf(self):
        self.ts[5:10] = np.NaN
        self.ts[15:20] = np.NaN

        val1 = self.ts.asOf(self.ts.index[7])
        val2 = self.ts.asOf(self.ts.index[19])

        self.assertEqual(val1, self.ts[4])
        self.assertEqual(val2, self.ts[14])

        # accepts strings
        val1 = self.ts.asOf(str(self.ts.index[7]))
        self.assertEqual(val1, self.ts[4])

        # in there
        self.assertEqual(self.ts.asOf(self.ts.index[3]), self.ts[3])

        # no as of value
        d = self.ts.index[0] - datetools.bday
        self.assert_(np.isnan(self.ts.asOf(d)))

    def test_merge(self):
        index, data = common.getMixedTypeDict()

        source = Series(data['B'], index=data['C'])
        target = Series(data['C'][:4], index=data['D'][:4])

        merged = target.merge(source)

        for k, v in merged.iteritems():
            self.assertEqual(v, source[target[k]])

        # input could be a dict
        merged = target.merge(source.toDict())

        for k, v in merged.iteritems():
            self.assertEqual(v, source[target[k]])

    def test_reindex(self):
        identity = self.series.reindex(self.series.index)
        self.assertEqual(id(self.series.index), id(identity.index))

        subIndex = self.series.index[10:20]
        subSeries = self.series.reindex(subIndex)

        for idx, val in subSeries.iteritems():
            self.assertEqual(val, self.series[idx])

        subIndex2 = self.ts.index[10:20]
        subTS = self.ts.reindex(subIndex2)

        for idx, val in subTS.iteritems():
            self.assertEqual(val, self.ts[idx])
        crapSeries = self.ts.reindex(subIndex)

        self.assert_(np.isnan(crapSeries).all())

        # This is extremely important for the Cython code to not screw up
        nonContigIndex = self.ts.index[::2]
        subNonContig = self.ts.reindex(nonContigIndex)
        for idx, val in subNonContig.iteritems():
            self.assertEqual(val, self.ts[idx])

        # bad fill method
        ts = self.ts[::2]
        self.assertRaises(Exception, ts.reindex, self.ts.index, fillMethod='foo')

        # corner case: pad empty series
        reindexed = self.empty.reindex(self.ts.index, fillMethod='pad')

        # pass non-Index
        reindexed = self.ts.reindex(list(self.ts.index))
        assert_series_equal(self.ts, reindexed)

    def test_reindex_int(self):
        ts = self.ts[::2]
        int_ts = Series(np.zeros(len(ts), dtype=int), index=ts.index)

        # this should work fine
        reindexed_int = int_ts.reindex(self.ts.index)

        # if NaNs introduced
        self.assert_(reindexed_int.dtype == np.float_)

        # NO NaNs introduced
        reindexed_int = int_ts.reindex(int_ts.index[::2])
        self.assert_(reindexed_int.dtype == np.int_)

    def test_reindex_bool(self):

        # A series other than float, int, string, or object
        ts = self.ts[::2]
        bool_ts = Series(np.zeros(len(ts), dtype=bool), index=ts.index)

        # this should work fine
        reindexed_bool = bool_ts.reindex(self.ts.index)

        # if NaNs introduced
        self.assert_(reindexed_bool.dtype == np.object_)

        # NO NaNs introduced
        reindexed_bool = bool_ts.reindex(bool_ts.index[::2])
        self.assert_(reindexed_bool.dtype == np.bool_)

    def test_reindex_bool_pad(self):
        # fail
        ts = self.ts[5:]
        bool_ts = Series(np.zeros(len(ts), dtype=bool), index=ts.index)
        filled_bool = bool_ts.reindex(self.ts.index, fillMethod='pad')
        self.assert_(isnull(filled_bool[:5]).all())

    def test_preserveRefs(self):
        sl = self.ts[5:10]
        seq = self.ts[[5,10,15]]
        sl[4] = np.NaN
        seq[1] = np.NaN
        self.assertFalse(np.isnan(self.ts[9]))
        self.assertFalse(np.isnan(self.ts[10]))

#-------------------------------------------------------------------------------
# TimeSeries-specific

    def test_fill(self):
        ts = Series([0., 1., 2., 3., 4.], index=common.makeDateIndex(5))

        self.assert_(np.array_equal(ts, ts.fill()))

        ts[2] = np.NaN

        self.assert_(np.array_equal(ts.fill(), [0., 1., 1., 3., 4.]))
        self.assert_(np.array_equal(ts.fill(method='backfill'), [0., 1., 3., 3., 4.]))

        self.assert_(np.array_equal(ts.fill(value=5), [0., 1., 5., 3., 4.]))

    def test_asfreq(self):
        ts = Series([0., 1., 2.], index=[datetime(2009, 10, 30),
                                         datetime(2009, 11, 30),
                                         datetime(2009, 12, 31)])

        daily_ts = ts.asfreq('WEEKDAY')
        monthly_ts = daily_ts.asfreq('EOM')
        self.assert_(np.array_equal(monthly_ts, ts))

        daily_ts = ts.asfreq('WEEKDAY', fillMethod='pad')
        monthly_ts = daily_ts.asfreq('EOM')
        self.assert_(np.array_equal(monthly_ts, ts))

        daily_ts = ts.asfreq(datetools.bday)
        monthly_ts = daily_ts.asfreq(datetools.bmonthEnd)
        self.assert_(np.array_equal(monthly_ts, ts))

    def test_interpolate(self):
        ts = Series(np.arange(len(self.ts), dtype=float), self.ts.index)

        ts_copy = ts.copy()
        ts_copy[5:10] = np.NaN

        linear_interp = ts_copy.interpolate(method='linear')
        self.assert_(np.array_equal(linear_interp, ts))

        ord_ts = Series([d.toordinal() for d in self.ts.index],
                        index=self.ts.index).astype(float)

        ord_ts_copy = ord_ts.copy()
        ord_ts_copy[5:10] = np.NaN

        time_interp = ord_ts_copy.interpolate(method='time')
        self.assert_(np.array_equal(time_interp, ord_ts))

        # try time interpolation on a non-TimeSeries
        self.assertRaises(Exception, self.series.interpolate, method='time')

    def test_weekday(self):
        # Just run the function
        weekdays = self.ts.weekday

    def test_diff(self):
        # Just run the function
        self.ts.diff()

    def test_autocorr(self):
        # Just run the function
        self.ts.autocorr()

    def test_firstValid(self):
        ts = self.ts.copy()
        ts[:5] = np.NaN

        index = ts._firstTimeWithValue()
        self.assertEqual(index, ts.index[5])

        ts[-5:] = np.NaN
        index = ts._lastTimeWithValue()
        self.assertEqual(index, ts.index[-6])

        ser = Series([], index=[])
        self.assert_(ser._lastTimeWithValue() is None)
        self.assert_(ser._firstTimeWithValue() is None)

    def test_lastValid(self):
        pass

#-------------------------------------------------------------------------------
# GroupBy

    def test_groupby(self):
        data = Series(np.arange(9) / 3, index=np.arange(9))

        index = np.arange(9)
        np.random.shuffle(index)
        data = data.reindex(index)

        grouped = data.groupby(lambda x: x // 3)

        repr(grouped.groups) # nothing else here

        for k, v in grouped:
            self.assertEqual(len(v), 3)

        agged = grouped.aggregate(np.mean)
        self.assertEqual(agged[1], 1)

        assert_series_equal(agged, grouped.agg(np.mean)) # shorthand

        transformed = grouped.transform(lambda x: x * x.sum())
        self.assertEqual(transformed[7], 12)

        value_grouped = data.groupby(data)
        assert_series_equal(value_grouped.aggregate(np.mean), agged)

        # complex agg
        agged = grouped.aggregate([np.mean, np.std])
        agged = grouped.aggregate({'one' : np.mean,
                                   'two' : np.std})

        group_constants = {
            0 : 10,
            1 : 20,
            2 : 30
        }
        agged = grouped.agg(lambda x: group_constants[x.groupName] + x.mean())
        self.assertEqual(agged[1], 21)

        # corner cases
        self.assertRaises(Exception, grouped._aggregate_named,
                          lambda x: x * 2)

    def test_groupby_transform(self):
        data = Series(np.arange(9) / 3, index=np.arange(9))

        index = np.arange(9)
        np.random.shuffle(index)
        data = data.reindex(index)

        grouped = data.groupby(lambda x: x // 3)

        transformed = grouped.transform(lambda x: x * x.sum())
        self.assertEqual(transformed[7], 12)

        # corner cases
        self.assertRaises(Exception, grouped.transform,
                          lambda x: x.mean())

if __name__ == '__main__':
    unittest.main()
