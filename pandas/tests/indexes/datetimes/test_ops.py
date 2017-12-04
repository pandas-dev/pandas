import pytz
import pytest
import dateutil
import warnings
import numpy as np
from datetime import datetime

from itertools import product
import pandas as pd
import pandas._libs.tslib as tslib
from pandas._libs.tslibs.offsets import shift_months
import pandas.util.testing as tm
from pandas import (DatetimeIndex, PeriodIndex, Series, Timestamp,
                    date_range, _np_version_under1p10, Index,
                    bdate_range)
from pandas.tseries.offsets import BMonthEnd, CDay, BDay
from pandas.tests.test_base import Ops


START, END = datetime(2009, 1, 1), datetime(2010, 1, 1)


class TestDatetimeIndexOps(Ops):
    tz = [None, 'UTC', 'Asia/Tokyo', 'US/Eastern', 'dateutil/Asia/Singapore',
          'dateutil/US/Pacific']

    def setup_method(self, method):
        super(TestDatetimeIndexOps, self).setup_method(method)
        mask = lambda x: (isinstance(x, DatetimeIndex) or
                          isinstance(x, PeriodIndex))
        self.is_valid_objs = [o for o in self.objs if mask(o)]
        self.not_valid_objs = [o for o in self.objs if not mask(o)]

    def test_ops_properties(self):
        f = lambda x: isinstance(x, DatetimeIndex)
        self.check_ops_properties(DatetimeIndex._field_ops, f)
        self.check_ops_properties(DatetimeIndex._object_ops, f)
        self.check_ops_properties(DatetimeIndex._bool_ops, f)

    def test_ops_properties_basic(self):

        # sanity check that the behavior didn't change
        # GH7206
        for op in ['year', 'day', 'second', 'weekday']:
            pytest.raises(TypeError, lambda x: getattr(self.dt_series, op))

        # attribute access should still work!
        s = Series(dict(year=2000, month=1, day=10))
        assert s.year == 2000
        assert s.month == 1
        assert s.day == 10
        pytest.raises(AttributeError, lambda: s.weekday)

    def test_astype_object(self):
        idx = pd.date_range(start='2013-01-01', periods=4, freq='M',
                            name='idx')
        expected_list = [Timestamp('2013-01-31'),
                         Timestamp('2013-02-28'),
                         Timestamp('2013-03-31'),
                         Timestamp('2013-04-30')]
        expected = pd.Index(expected_list, dtype=object, name='idx')
        result = idx.astype(object)
        assert isinstance(result, Index)

        assert result.dtype == object
        tm.assert_index_equal(result, expected)
        assert result.name == expected.name
        assert idx.tolist() == expected_list

        idx = pd.date_range(start='2013-01-01', periods=4, freq='M',
                            name='idx', tz='Asia/Tokyo')
        expected_list = [Timestamp('2013-01-31', tz='Asia/Tokyo'),
                         Timestamp('2013-02-28', tz='Asia/Tokyo'),
                         Timestamp('2013-03-31', tz='Asia/Tokyo'),
                         Timestamp('2013-04-30', tz='Asia/Tokyo')]
        expected = pd.Index(expected_list, dtype=object, name='idx')
        result = idx.astype(object)
        assert isinstance(result, Index)
        assert result.dtype == object
        tm.assert_index_equal(result, expected)
        assert result.name == expected.name
        assert idx.tolist() == expected_list

        idx = DatetimeIndex([datetime(2013, 1, 1), datetime(2013, 1, 2),
                             pd.NaT, datetime(2013, 1, 4)], name='idx')
        expected_list = [Timestamp('2013-01-01'),
                         Timestamp('2013-01-02'), pd.NaT,
                         Timestamp('2013-01-04')]
        expected = pd.Index(expected_list, dtype=object, name='idx')
        result = idx.astype(object)
        assert isinstance(result, Index)
        assert result.dtype == object
        tm.assert_index_equal(result, expected)
        assert result.name == expected.name
        assert idx.tolist() == expected_list

    def test_minmax(self):
        for tz in self.tz:
            # monotonic
            idx1 = pd.DatetimeIndex(['2011-01-01', '2011-01-02',
                                     '2011-01-03'], tz=tz)
            assert idx1.is_monotonic

            # non-monotonic
            idx2 = pd.DatetimeIndex(['2011-01-01', pd.NaT, '2011-01-03',
                                     '2011-01-02', pd.NaT], tz=tz)
            assert not idx2.is_monotonic

            for idx in [idx1, idx2]:
                assert idx.min() == Timestamp('2011-01-01', tz=tz)
                assert idx.max() == Timestamp('2011-01-03', tz=tz)
                assert idx.argmin() == 0
                assert idx.argmax() == 2

        for op in ['min', 'max']:
            # Return NaT
            obj = DatetimeIndex([])
            assert pd.isna(getattr(obj, op)())

            obj = DatetimeIndex([pd.NaT])
            assert pd.isna(getattr(obj, op)())

            obj = DatetimeIndex([pd.NaT, pd.NaT, pd.NaT])
            assert pd.isna(getattr(obj, op)())

    def test_numpy_minmax(self):
        dr = pd.date_range(start='2016-01-15', end='2016-01-20')

        assert np.min(dr) == Timestamp('2016-01-15 00:00:00', freq='D')
        assert np.max(dr) == Timestamp('2016-01-20 00:00:00', freq='D')

        errmsg = "the 'out' parameter is not supported"
        tm.assert_raises_regex(ValueError, errmsg, np.min, dr, out=0)
        tm.assert_raises_regex(ValueError, errmsg, np.max, dr, out=0)

        assert np.argmin(dr) == 0
        assert np.argmax(dr) == 5

        if not _np_version_under1p10:
            errmsg = "the 'out' parameter is not supported"
            tm.assert_raises_regex(
                ValueError, errmsg, np.argmin, dr, out=0)
            tm.assert_raises_regex(
                ValueError, errmsg, np.argmax, dr, out=0)

    def test_round(self):
        for tz in self.tz:
            rng = pd.date_range(start='2016-01-01', periods=5,
                                freq='30Min', tz=tz)
            elt = rng[1]

            expected_rng = DatetimeIndex([
                Timestamp('2016-01-01 00:00:00', tz=tz, freq='30T'),
                Timestamp('2016-01-01 00:00:00', tz=tz, freq='30T'),
                Timestamp('2016-01-01 01:00:00', tz=tz, freq='30T'),
                Timestamp('2016-01-01 02:00:00', tz=tz, freq='30T'),
                Timestamp('2016-01-01 02:00:00', tz=tz, freq='30T'),
            ])
            expected_elt = expected_rng[1]

            tm.assert_index_equal(rng.round(freq='H'), expected_rng)
            assert elt.round(freq='H') == expected_elt

            msg = pd.tseries.frequencies._INVALID_FREQ_ERROR
            with tm.assert_raises_regex(ValueError, msg):
                rng.round(freq='foo')
            with tm.assert_raises_regex(ValueError, msg):
                elt.round(freq='foo')

            msg = "<MonthEnd> is a non-fixed frequency"
            tm.assert_raises_regex(ValueError, msg, rng.round, freq='M')
            tm.assert_raises_regex(ValueError, msg, elt.round, freq='M')

            # GH 14440 & 15578
            index = pd.DatetimeIndex(['2016-10-17 12:00:00.0015'], tz=tz)
            result = index.round('ms')
            expected = pd.DatetimeIndex(['2016-10-17 12:00:00.002000'], tz=tz)
            tm.assert_index_equal(result, expected)

            for freq in ['us', 'ns']:
                tm.assert_index_equal(index, index.round(freq))

            index = pd.DatetimeIndex(['2016-10-17 12:00:00.00149'], tz=tz)
            result = index.round('ms')
            expected = pd.DatetimeIndex(['2016-10-17 12:00:00.001000'], tz=tz)
            tm.assert_index_equal(result, expected)

            index = pd.DatetimeIndex(['2016-10-17 12:00:00.001501031'])
            result = index.round('10ns')
            expected = pd.DatetimeIndex(['2016-10-17 12:00:00.001501030'])
            tm.assert_index_equal(result, expected)

            with tm.assert_produces_warning():
                ts = '2016-10-17 12:00:00.001501031'
                pd.DatetimeIndex([ts]).round('1010ns')

    def test_repeat_range(self):
        rng = date_range('1/1/2000', '1/1/2001')

        result = rng.repeat(5)
        assert result.freq is None
        assert len(result) == 5 * len(rng)

        for tz in self.tz:
            index = pd.date_range('2001-01-01', periods=2, freq='D', tz=tz)
            exp = pd.DatetimeIndex(['2001-01-01', '2001-01-01',
                                    '2001-01-02', '2001-01-02'], tz=tz)
            for res in [index.repeat(2), np.repeat(index, 2)]:
                tm.assert_index_equal(res, exp)
                assert res.freq is None

            index = pd.date_range('2001-01-01', periods=2, freq='2D', tz=tz)
            exp = pd.DatetimeIndex(['2001-01-01', '2001-01-01',
                                    '2001-01-03', '2001-01-03'], tz=tz)
            for res in [index.repeat(2), np.repeat(index, 2)]:
                tm.assert_index_equal(res, exp)
                assert res.freq is None

            index = pd.DatetimeIndex(['2001-01-01', 'NaT', '2003-01-01'],
                                     tz=tz)
            exp = pd.DatetimeIndex(['2001-01-01', '2001-01-01', '2001-01-01',
                                    'NaT', 'NaT', 'NaT',
                                    '2003-01-01', '2003-01-01', '2003-01-01'],
                                   tz=tz)
            for res in [index.repeat(3), np.repeat(index, 3)]:
                tm.assert_index_equal(res, exp)
                assert res.freq is None

    def test_repeat(self):
        reps = 2
        msg = "the 'axis' parameter is not supported"

        for tz in self.tz:
            rng = pd.date_range(start='2016-01-01', periods=2,
                                freq='30Min', tz=tz)

            expected_rng = DatetimeIndex([
                Timestamp('2016-01-01 00:00:00', tz=tz, freq='30T'),
                Timestamp('2016-01-01 00:00:00', tz=tz, freq='30T'),
                Timestamp('2016-01-01 00:30:00', tz=tz, freq='30T'),
                Timestamp('2016-01-01 00:30:00', tz=tz, freq='30T'),
            ])

            res = rng.repeat(reps)
            tm.assert_index_equal(res, expected_rng)
            assert res.freq is None

            tm.assert_index_equal(np.repeat(rng, reps), expected_rng)
            tm.assert_raises_regex(ValueError, msg, np.repeat,
                                   rng, reps, axis=1)

    def test_representation(self):

        idx = []
        idx.append(DatetimeIndex([], freq='D'))
        idx.append(DatetimeIndex(['2011-01-01'], freq='D'))
        idx.append(DatetimeIndex(['2011-01-01', '2011-01-02'], freq='D'))
        idx.append(DatetimeIndex(
            ['2011-01-01', '2011-01-02', '2011-01-03'], freq='D'))
        idx.append(DatetimeIndex(
            ['2011-01-01 09:00', '2011-01-01 10:00', '2011-01-01 11:00'
             ], freq='H', tz='Asia/Tokyo'))
        idx.append(DatetimeIndex(
            ['2011-01-01 09:00', '2011-01-01 10:00', pd.NaT], tz='US/Eastern'))
        idx.append(DatetimeIndex(
            ['2011-01-01 09:00', '2011-01-01 10:00', pd.NaT], tz='UTC'))

        exp = []
        exp.append("""DatetimeIndex([], dtype='datetime64[ns]', freq='D')""")
        exp.append("DatetimeIndex(['2011-01-01'], dtype='datetime64[ns]', "
                   "freq='D')")
        exp.append("DatetimeIndex(['2011-01-01', '2011-01-02'], "
                   "dtype='datetime64[ns]', freq='D')")
        exp.append("DatetimeIndex(['2011-01-01', '2011-01-02', '2011-01-03'], "
                   "dtype='datetime64[ns]', freq='D')")
        exp.append("DatetimeIndex(['2011-01-01 09:00:00+09:00', "
                   "'2011-01-01 10:00:00+09:00', '2011-01-01 11:00:00+09:00']"
                   ", dtype='datetime64[ns, Asia/Tokyo]', freq='H')")
        exp.append("DatetimeIndex(['2011-01-01 09:00:00-05:00', "
                   "'2011-01-01 10:00:00-05:00', 'NaT'], "
                   "dtype='datetime64[ns, US/Eastern]', freq=None)")
        exp.append("DatetimeIndex(['2011-01-01 09:00:00+00:00', "
                   "'2011-01-01 10:00:00+00:00', 'NaT'], "
                   "dtype='datetime64[ns, UTC]', freq=None)""")

        with pd.option_context('display.width', 300):
            for indx, expected in zip(idx, exp):
                for func in ['__repr__', '__unicode__', '__str__']:
                    result = getattr(indx, func)()
                    assert result == expected

    def test_representation_to_series(self):
        idx1 = DatetimeIndex([], freq='D')
        idx2 = DatetimeIndex(['2011-01-01'], freq='D')
        idx3 = DatetimeIndex(['2011-01-01', '2011-01-02'], freq='D')
        idx4 = DatetimeIndex(
            ['2011-01-01', '2011-01-02', '2011-01-03'], freq='D')
        idx5 = DatetimeIndex(['2011-01-01 09:00', '2011-01-01 10:00',
                              '2011-01-01 11:00'], freq='H', tz='Asia/Tokyo')
        idx6 = DatetimeIndex(['2011-01-01 09:00', '2011-01-01 10:00', pd.NaT],
                             tz='US/Eastern')
        idx7 = DatetimeIndex(['2011-01-01 09:00', '2011-01-02 10:15'])

        exp1 = """Series([], dtype: datetime64[ns])"""

        exp2 = ("0   2011-01-01\n"
                "dtype: datetime64[ns]")

        exp3 = ("0   2011-01-01\n"
                "1   2011-01-02\n"
                "dtype: datetime64[ns]")

        exp4 = ("0   2011-01-01\n"
                "1   2011-01-02\n"
                "2   2011-01-03\n"
                "dtype: datetime64[ns]")

        exp5 = ("0   2011-01-01 09:00:00+09:00\n"
                "1   2011-01-01 10:00:00+09:00\n"
                "2   2011-01-01 11:00:00+09:00\n"
                "dtype: datetime64[ns, Asia/Tokyo]")

        exp6 = ("0   2011-01-01 09:00:00-05:00\n"
                "1   2011-01-01 10:00:00-05:00\n"
                "2                         NaT\n"
                "dtype: datetime64[ns, US/Eastern]")

        exp7 = ("0   2011-01-01 09:00:00\n"
                "1   2011-01-02 10:15:00\n"
                "dtype: datetime64[ns]")

        with pd.option_context('display.width', 300):
            for idx, expected in zip([idx1, idx2, idx3, idx4,
                                      idx5, idx6, idx7],
                                     [exp1, exp2, exp3, exp4,
                                      exp5, exp6, exp7]):
                result = repr(Series(idx))
                assert result == expected

    def test_summary(self):
        # GH9116
        idx1 = DatetimeIndex([], freq='D')
        idx2 = DatetimeIndex(['2011-01-01'], freq='D')
        idx3 = DatetimeIndex(['2011-01-01', '2011-01-02'], freq='D')
        idx4 = DatetimeIndex(
            ['2011-01-01', '2011-01-02', '2011-01-03'], freq='D')
        idx5 = DatetimeIndex(['2011-01-01 09:00', '2011-01-01 10:00',
                              '2011-01-01 11:00'],
                             freq='H', tz='Asia/Tokyo')
        idx6 = DatetimeIndex(['2011-01-01 09:00', '2011-01-01 10:00', pd.NaT],
                             tz='US/Eastern')

        exp1 = ("DatetimeIndex: 0 entries\n"
                "Freq: D")

        exp2 = ("DatetimeIndex: 1 entries, 2011-01-01 to 2011-01-01\n"
                "Freq: D")

        exp3 = ("DatetimeIndex: 2 entries, 2011-01-01 to 2011-01-02\n"
                "Freq: D")

        exp4 = ("DatetimeIndex: 3 entries, 2011-01-01 to 2011-01-03\n"
                "Freq: D")

        exp5 = ("DatetimeIndex: 3 entries, 2011-01-01 09:00:00+09:00 "
                "to 2011-01-01 11:00:00+09:00\n"
                "Freq: H")

        exp6 = """DatetimeIndex: 3 entries, 2011-01-01 09:00:00-05:00 to NaT"""

        for idx, expected in zip([idx1, idx2, idx3, idx4, idx5, idx6],
                                 [exp1, exp2, exp3, exp4, exp5, exp6]):
            result = idx.summary()
            assert result == expected

    def test_resolution(self):
        for freq, expected in zip(['A', 'Q', 'M', 'D', 'H', 'T',
                                   'S', 'L', 'U'],
                                  ['day', 'day', 'day', 'day', 'hour',
                                   'minute', 'second', 'millisecond',
                                   'microsecond']):
            for tz in self.tz:
                idx = pd.date_range(start='2013-04-01', periods=30, freq=freq,
                                    tz=tz)
                assert idx.resolution == expected

    def test_comp_nat(self):
        left = pd.DatetimeIndex([pd.Timestamp('2011-01-01'), pd.NaT,
                                 pd.Timestamp('2011-01-03')])
        right = pd.DatetimeIndex([pd.NaT, pd.NaT, pd.Timestamp('2011-01-03')])

        for lhs, rhs in [(left, right),
                         (left.astype(object), right.astype(object))]:
            result = rhs == lhs
            expected = np.array([False, False, True])
            tm.assert_numpy_array_equal(result, expected)

            result = lhs != rhs
            expected = np.array([True, True, False])
            tm.assert_numpy_array_equal(result, expected)

            expected = np.array([False, False, False])
            tm.assert_numpy_array_equal(lhs == pd.NaT, expected)
            tm.assert_numpy_array_equal(pd.NaT == rhs, expected)

            expected = np.array([True, True, True])
            tm.assert_numpy_array_equal(lhs != pd.NaT, expected)
            tm.assert_numpy_array_equal(pd.NaT != lhs, expected)

            expected = np.array([False, False, False])
            tm.assert_numpy_array_equal(lhs < pd.NaT, expected)
            tm.assert_numpy_array_equal(pd.NaT > lhs, expected)

    def test_value_counts_unique(self):
        # GH 7735
        for tz in self.tz:
            idx = pd.date_range('2011-01-01 09:00', freq='H', periods=10)
            # create repeated values, 'n'th element is repeated by n+1 times
            idx = DatetimeIndex(np.repeat(idx.values, range(1, len(idx) + 1)),
                                tz=tz)

            exp_idx = pd.date_range('2011-01-01 18:00', freq='-1H', periods=10,
                                    tz=tz)
            expected = Series(range(10, 0, -1), index=exp_idx, dtype='int64')

            for obj in [idx, Series(idx)]:
                tm.assert_series_equal(obj.value_counts(), expected)

            expected = pd.date_range('2011-01-01 09:00', freq='H', periods=10,
                                     tz=tz)
            tm.assert_index_equal(idx.unique(), expected)

            idx = DatetimeIndex(['2013-01-01 09:00', '2013-01-01 09:00',
                                 '2013-01-01 09:00', '2013-01-01 08:00',
                                 '2013-01-01 08:00', pd.NaT], tz=tz)

            exp_idx = DatetimeIndex(['2013-01-01 09:00', '2013-01-01 08:00'],
                                    tz=tz)
            expected = Series([3, 2], index=exp_idx)

            for obj in [idx, Series(idx)]:
                tm.assert_series_equal(obj.value_counts(), expected)

            exp_idx = DatetimeIndex(['2013-01-01 09:00', '2013-01-01 08:00',
                                     pd.NaT], tz=tz)
            expected = Series([3, 2, 1], index=exp_idx)

            for obj in [idx, Series(idx)]:
                tm.assert_series_equal(obj.value_counts(dropna=False),
                                       expected)

            tm.assert_index_equal(idx.unique(), exp_idx)

    def test_nonunique_contains(self):
        # GH 9512
        for idx in map(DatetimeIndex,
                       ([0, 1, 0], [0, 0, -1], [0, -1, -1],
                        ['2015', '2015', '2016'], ['2015', '2015', '2014'])):
            assert idx[0] in idx

    def test_order(self):
        # with freq
        idx1 = DatetimeIndex(['2011-01-01', '2011-01-02',
                              '2011-01-03'], freq='D', name='idx')
        idx2 = DatetimeIndex(['2011-01-01 09:00', '2011-01-01 10:00',
                              '2011-01-01 11:00'], freq='H',
                             tz='Asia/Tokyo', name='tzidx')

        for idx in [idx1, idx2]:
            ordered = idx.sort_values()
            tm.assert_index_equal(ordered, idx)
            assert ordered.freq == idx.freq

            ordered = idx.sort_values(ascending=False)
            expected = idx[::-1]
            tm.assert_index_equal(ordered, expected)
            assert ordered.freq == expected.freq
            assert ordered.freq.n == -1

            ordered, indexer = idx.sort_values(return_indexer=True)
            tm.assert_index_equal(ordered, idx)
            tm.assert_numpy_array_equal(indexer, np.array([0, 1, 2]),
                                        check_dtype=False)
            assert ordered.freq == idx.freq

            ordered, indexer = idx.sort_values(return_indexer=True,
                                               ascending=False)
            expected = idx[::-1]
            tm.assert_index_equal(ordered, expected)
            tm.assert_numpy_array_equal(indexer,
                                        np.array([2, 1, 0]),
                                        check_dtype=False)
            assert ordered.freq == expected.freq
            assert ordered.freq.n == -1

        # without freq
        for tz in self.tz:
            idx1 = DatetimeIndex(['2011-01-01', '2011-01-03', '2011-01-05',
                                  '2011-01-02', '2011-01-01'],
                                 tz=tz, name='idx1')
            exp1 = DatetimeIndex(['2011-01-01', '2011-01-01', '2011-01-02',
                                  '2011-01-03', '2011-01-05'],
                                 tz=tz, name='idx1')

            idx2 = DatetimeIndex(['2011-01-01', '2011-01-03', '2011-01-05',
                                  '2011-01-02', '2011-01-01'],
                                 tz=tz, name='idx2')

            exp2 = DatetimeIndex(['2011-01-01', '2011-01-01', '2011-01-02',
                                  '2011-01-03', '2011-01-05'],
                                 tz=tz, name='idx2')

            idx3 = DatetimeIndex([pd.NaT, '2011-01-03', '2011-01-05',
                                  '2011-01-02', pd.NaT], tz=tz, name='idx3')
            exp3 = DatetimeIndex([pd.NaT, pd.NaT, '2011-01-02', '2011-01-03',
                                  '2011-01-05'], tz=tz, name='idx3')

            for idx, expected in [(idx1, exp1), (idx2, exp2), (idx3, exp3)]:
                ordered = idx.sort_values()
                tm.assert_index_equal(ordered, expected)
                assert ordered.freq is None

                ordered = idx.sort_values(ascending=False)
                tm.assert_index_equal(ordered, expected[::-1])
                assert ordered.freq is None

                ordered, indexer = idx.sort_values(return_indexer=True)
                tm.assert_index_equal(ordered, expected)

                exp = np.array([0, 4, 3, 1, 2])
                tm.assert_numpy_array_equal(indexer, exp, check_dtype=False)
                assert ordered.freq is None

                ordered, indexer = idx.sort_values(return_indexer=True,
                                                   ascending=False)
                tm.assert_index_equal(ordered, expected[::-1])

                exp = np.array([2, 1, 3, 4, 0])
                tm.assert_numpy_array_equal(indexer, exp, check_dtype=False)
                assert ordered.freq is None

    def test_drop_duplicates_metadata(self):
        # GH 10115
        idx = pd.date_range('2011-01-01', '2011-01-31', freq='D', name='idx')
        result = idx.drop_duplicates()
        tm.assert_index_equal(idx, result)
        assert idx.freq == result.freq

        idx_dup = idx.append(idx)
        assert idx_dup.freq is None  # freq is reset
        result = idx_dup.drop_duplicates()
        tm.assert_index_equal(idx, result)
        assert result.freq is None

    def test_drop_duplicates(self):
        # to check Index/Series compat
        base = pd.date_range('2011-01-01', '2011-01-31', freq='D', name='idx')
        idx = base.append(base[:5])

        res = idx.drop_duplicates()
        tm.assert_index_equal(res, base)
        res = Series(idx).drop_duplicates()
        tm.assert_series_equal(res, Series(base))

        res = idx.drop_duplicates(keep='last')
        exp = base[5:].append(base[:5])
        tm.assert_index_equal(res, exp)
        res = Series(idx).drop_duplicates(keep='last')
        tm.assert_series_equal(res, Series(exp, index=np.arange(5, 36)))

        res = idx.drop_duplicates(keep=False)
        tm.assert_index_equal(res, base[5:])
        res = Series(idx).drop_duplicates(keep=False)
        tm.assert_series_equal(res, Series(base[5:], index=np.arange(5, 31)))

    def test_infer_freq(self):
        # GH 11018
        for freq in ['A', '2A', '-2A', 'Q', '-1Q', 'M', '-1M', 'D', '3D',
                     '-3D', 'W', '-1W', 'H', '2H', '-2H', 'T', '2T', 'S',
                     '-3S']:
            idx = pd.date_range('2011-01-01 09:00:00', freq=freq, periods=10)
            result = pd.DatetimeIndex(idx.asi8, freq='infer')
            tm.assert_index_equal(idx, result)
            assert result.freq == freq

    def test_nat_new(self):
        idx = pd.date_range('2011-01-01', freq='D', periods=5, name='x')
        result = idx._nat_new()
        exp = pd.DatetimeIndex([pd.NaT] * 5, name='x')
        tm.assert_index_equal(result, exp)

        result = idx._nat_new(box=False)
        exp = np.array([tslib.iNaT] * 5, dtype=np.int64)
        tm.assert_numpy_array_equal(result, exp)

    def test_shift(self):
        # GH 9903
        for tz in self.tz:
            idx = pd.DatetimeIndex([], name='xxx', tz=tz)
            tm.assert_index_equal(idx.shift(0, freq='H'), idx)
            tm.assert_index_equal(idx.shift(3, freq='H'), idx)

            idx = pd.DatetimeIndex(['2011-01-01 10:00', '2011-01-01 11:00'
                                    '2011-01-01 12:00'], name='xxx', tz=tz)
            tm.assert_index_equal(idx.shift(0, freq='H'), idx)
            exp = pd.DatetimeIndex(['2011-01-01 13:00', '2011-01-01 14:00'
                                    '2011-01-01 15:00'], name='xxx', tz=tz)
            tm.assert_index_equal(idx.shift(3, freq='H'), exp)
            exp = pd.DatetimeIndex(['2011-01-01 07:00', '2011-01-01 08:00'
                                    '2011-01-01 09:00'], name='xxx', tz=tz)
            tm.assert_index_equal(idx.shift(-3, freq='H'), exp)

    def test_nat(self):
        assert pd.DatetimeIndex._na_value is pd.NaT
        assert pd.DatetimeIndex([])._na_value is pd.NaT

        for tz in [None, 'US/Eastern', 'UTC']:
            idx = pd.DatetimeIndex(['2011-01-01', '2011-01-02'], tz=tz)
            assert idx._can_hold_na

            tm.assert_numpy_array_equal(idx._isnan, np.array([False, False]))
            assert not idx.hasnans
            tm.assert_numpy_array_equal(idx._nan_idxs,
                                        np.array([], dtype=np.intp))

            idx = pd.DatetimeIndex(['2011-01-01', 'NaT'], tz=tz)
            assert idx._can_hold_na

            tm.assert_numpy_array_equal(idx._isnan, np.array([False, True]))
            assert idx.hasnans
            tm.assert_numpy_array_equal(idx._nan_idxs,
                                        np.array([1], dtype=np.intp))

    def test_equals(self):
        # GH 13107
        for tz in [None, 'UTC', 'US/Eastern', 'Asia/Tokyo']:
            idx = pd.DatetimeIndex(['2011-01-01', '2011-01-02', 'NaT'])
            assert idx.equals(idx)
            assert idx.equals(idx.copy())
            assert idx.equals(idx.astype(object))
            assert idx.astype(object).equals(idx)
            assert idx.astype(object).equals(idx.astype(object))
            assert not idx.equals(list(idx))
            assert not idx.equals(pd.Series(idx))

            idx2 = pd.DatetimeIndex(['2011-01-01', '2011-01-02', 'NaT'],
                                    tz='US/Pacific')
            assert not idx.equals(idx2)
            assert not idx.equals(idx2.copy())
            assert not idx.equals(idx2.astype(object))
            assert not idx.astype(object).equals(idx2)
            assert not idx.equals(list(idx2))
            assert not idx.equals(pd.Series(idx2))

            # same internal, different tz
            idx3 = pd.DatetimeIndex._simple_new(idx.asi8, tz='US/Pacific')
            tm.assert_numpy_array_equal(idx.asi8, idx3.asi8)
            assert not idx.equals(idx3)
            assert not idx.equals(idx3.copy())
            assert not idx.equals(idx3.astype(object))
            assert not idx.astype(object).equals(idx3)
            assert not idx.equals(list(idx3))
            assert not idx.equals(pd.Series(idx3))


@pytest.mark.parametrize('years,months', product([-1, 0, 1], [-2, 0, 2]))
def test_shift_months(years, months):
    s = DatetimeIndex([Timestamp('2000-01-05 00:15:00'),
                       Timestamp('2000-01-31 00:23:00'),
                       Timestamp('2000-01-01'),
                       Timestamp('2000-02-29'),
                       Timestamp('2000-12-31')])
    actual = DatetimeIndex(shift_months(s.asi8, years * 12 + months))
    expected = DatetimeIndex([x + pd.offsets.DateOffset(
        years=years, months=months) for x in s])
    tm.assert_index_equal(actual, expected)


class TestBusinessDatetimeIndex(object):

    def setup_method(self, method):
        self.rng = bdate_range(START, END)

    def test_comparison(self):
        d = self.rng[10]

        comp = self.rng > d
        assert comp[11]
        assert not comp[9]

    def test_pickle_unpickle(self):
        unpickled = tm.round_trip_pickle(self.rng)
        assert unpickled.offset is not None

    def test_copy(self):
        cp = self.rng.copy()
        repr(cp)
        tm.assert_index_equal(cp, self.rng)

    def test_repr(self):
        # only really care that it works
        repr(self.rng)

    def test_shift(self):
        shifted = self.rng.shift(5)
        assert shifted[0] == self.rng[5]
        assert shifted.offset == self.rng.offset

        shifted = self.rng.shift(-5)
        assert shifted[5] == self.rng[0]
        assert shifted.offset == self.rng.offset

        shifted = self.rng.shift(0)
        assert shifted[0] == self.rng[0]
        assert shifted.offset == self.rng.offset

        rng = date_range(START, END, freq=BMonthEnd())
        shifted = rng.shift(1, freq=BDay())
        assert shifted[0] == rng[0] + BDay()

    def test_summary(self):
        self.rng.summary()
        self.rng[2:2].summary()

    def test_summary_pytz(self):
        bdate_range('1/1/2005', '1/1/2009', tz=pytz.utc).summary()

    def test_summary_dateutil(self):
        bdate_range('1/1/2005', '1/1/2009', tz=dateutil.tz.tzutc()).summary()

    def test_equals(self):
        assert not self.rng.equals(list(self.rng))

    def test_identical(self):
        t1 = self.rng.copy()
        t2 = self.rng.copy()
        assert t1.identical(t2)

        # name
        t1 = t1.rename('foo')
        assert t1.equals(t2)
        assert not t1.identical(t2)
        t2 = t2.rename('foo')
        assert t1.identical(t2)

        # freq
        t2v = Index(t2.values)
        assert t1.equals(t2v)
        assert not t1.identical(t2v)


class TestCustomDatetimeIndex(object):
    def setup_method(self, method):
        self.rng = bdate_range(START, END, freq='C')

    def test_comparison(self):
        d = self.rng[10]

        comp = self.rng > d
        assert comp[11]
        assert not comp[9]

    def test_copy(self):
        cp = self.rng.copy()
        repr(cp)
        tm.assert_index_equal(cp, self.rng)

    def test_repr(self):
        # only really care that it works
        repr(self.rng)

    def test_shift(self):

        shifted = self.rng.shift(5)
        assert shifted[0] == self.rng[5]
        assert shifted.offset == self.rng.offset

        shifted = self.rng.shift(-5)
        assert shifted[5] == self.rng[0]
        assert shifted.offset == self.rng.offset

        shifted = self.rng.shift(0)
        assert shifted[0] == self.rng[0]
        assert shifted.offset == self.rng.offset

        # PerformanceWarning
        with warnings.catch_warnings(record=True):
            rng = date_range(START, END, freq=BMonthEnd())
            shifted = rng.shift(1, freq=CDay())
            assert shifted[0] == rng[0] + CDay()

    def test_pickle_unpickle(self):
        unpickled = tm.round_trip_pickle(self.rng)
        assert unpickled.offset is not None

    def test_summary(self):
        self.rng.summary()
        self.rng[2:2].summary()

    def test_summary_pytz(self):
        bdate_range('1/1/2005', '1/1/2009', freq='C', tz=pytz.utc).summary()

    def test_summary_dateutil(self):
        bdate_range('1/1/2005', '1/1/2009', freq='C',
                    tz=dateutil.tz.tzutc()).summary()

    def test_equals(self):
        assert not self.rng.equals(list(self.rng))
