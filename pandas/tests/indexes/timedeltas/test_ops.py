import pytest

import numpy as np
from datetime import timedelta

import pandas as pd
import pandas.util.testing as tm
from pandas import to_timedelta
from pandas import (Series, Timedelta, Timestamp, TimedeltaIndex,
                    timedelta_range,
                    _np_version_under1p10, Index)
from pandas._libs.tslib import iNaT
from pandas.tests.test_base import Ops


class TestTimedeltaIndexOps(Ops):
    def setup_method(self, method):
        super(TestTimedeltaIndexOps, self).setup_method(method)
        mask = lambda x: isinstance(x, TimedeltaIndex)
        self.is_valid_objs = [o for o in self.objs if mask(o)]
        self.not_valid_objs = []

    def test_ops_properties(self):
        f = lambda x: isinstance(x, TimedeltaIndex)
        self.check_ops_properties(TimedeltaIndex._field_ops, f)
        self.check_ops_properties(TimedeltaIndex._object_ops, f)

    def test_astype_object(self):
        idx = timedelta_range(start='1 days', periods=4, freq='D', name='idx')
        expected_list = [Timedelta('1 days'), Timedelta('2 days'),
                         Timedelta('3 days'), Timedelta('4 days')]
        expected = pd.Index(expected_list, dtype=object, name='idx')
        result = idx.astype(object)
        assert isinstance(result, Index)

        assert result.dtype == object
        tm.assert_index_equal(result, expected)
        assert result.name == expected.name
        assert idx.tolist() == expected_list

        idx = TimedeltaIndex([timedelta(days=1), timedelta(days=2), pd.NaT,
                              timedelta(days=4)], name='idx')
        expected_list = [Timedelta('1 days'), Timedelta('2 days'), pd.NaT,
                         Timedelta('4 days')]
        expected = pd.Index(expected_list, dtype=object, name='idx')
        result = idx.astype(object)
        assert isinstance(result, Index)
        assert result.dtype == object
        tm.assert_index_equal(result, expected)
        assert result.name == expected.name
        assert idx.tolist() == expected_list

    def test_minmax(self):

        # monotonic
        idx1 = TimedeltaIndex(['1 days', '2 days', '3 days'])
        assert idx1.is_monotonic

        # non-monotonic
        idx2 = TimedeltaIndex(['1 days', np.nan, '3 days', 'NaT'])
        assert not idx2.is_monotonic

        for idx in [idx1, idx2]:
            assert idx.min() == Timedelta('1 days')
            assert idx.max() == Timedelta('3 days')
            assert idx.argmin() == 0
            assert idx.argmax() == 2

        for op in ['min', 'max']:
            # Return NaT
            obj = TimedeltaIndex([])
            assert pd.isna(getattr(obj, op)())

            obj = TimedeltaIndex([pd.NaT])
            assert pd.isna(getattr(obj, op)())

            obj = TimedeltaIndex([pd.NaT, pd.NaT, pd.NaT])
            assert pd.isna(getattr(obj, op)())

    def test_numpy_minmax(self):
        dr = pd.date_range(start='2016-01-15', end='2016-01-20')
        td = TimedeltaIndex(np.asarray(dr))

        assert np.min(td) == Timedelta('16815 days')
        assert np.max(td) == Timedelta('16820 days')

        errmsg = "the 'out' parameter is not supported"
        tm.assert_raises_regex(ValueError, errmsg, np.min, td, out=0)
        tm.assert_raises_regex(ValueError, errmsg, np.max, td, out=0)

        assert np.argmin(td) == 0
        assert np.argmax(td) == 5

        if not _np_version_under1p10:
            errmsg = "the 'out' parameter is not supported"
            tm.assert_raises_regex(
                ValueError, errmsg, np.argmin, td, out=0)
            tm.assert_raises_regex(
                ValueError, errmsg, np.argmax, td, out=0)

    def test_round(self):
        td = pd.timedelta_range(start='16801 days', periods=5, freq='30Min')
        elt = td[1]

        expected_rng = TimedeltaIndex([
            Timedelta('16801 days 00:00:00'),
            Timedelta('16801 days 00:00:00'),
            Timedelta('16801 days 01:00:00'),
            Timedelta('16801 days 02:00:00'),
            Timedelta('16801 days 02:00:00'),
        ])
        expected_elt = expected_rng[1]

        tm.assert_index_equal(td.round(freq='H'), expected_rng)
        assert elt.round(freq='H') == expected_elt

        msg = pd.tseries.frequencies._INVALID_FREQ_ERROR
        with tm.assert_raises_regex(ValueError, msg):
            td.round(freq='foo')
        with tm.assert_raises_regex(ValueError, msg):
            elt.round(freq='foo')

        msg = "<MonthEnd> is a non-fixed frequency"
        tm.assert_raises_regex(ValueError, msg, td.round, freq='M')
        tm.assert_raises_regex(ValueError, msg, elt.round, freq='M')

    def test_representation(self):
        idx1 = TimedeltaIndex([], freq='D')
        idx2 = TimedeltaIndex(['1 days'], freq='D')
        idx3 = TimedeltaIndex(['1 days', '2 days'], freq='D')
        idx4 = TimedeltaIndex(['1 days', '2 days', '3 days'], freq='D')
        idx5 = TimedeltaIndex(['1 days 00:00:01', '2 days', '3 days'])

        exp1 = """TimedeltaIndex([], dtype='timedelta64[ns]', freq='D')"""

        exp2 = ("TimedeltaIndex(['1 days'], dtype='timedelta64[ns]', "
                "freq='D')")

        exp3 = ("TimedeltaIndex(['1 days', '2 days'], "
                "dtype='timedelta64[ns]', freq='D')")

        exp4 = ("TimedeltaIndex(['1 days', '2 days', '3 days'], "
                "dtype='timedelta64[ns]', freq='D')")

        exp5 = ("TimedeltaIndex(['1 days 00:00:01', '2 days 00:00:00', "
                "'3 days 00:00:00'], dtype='timedelta64[ns]', freq=None)")

        with pd.option_context('display.width', 300):
            for idx, expected in zip([idx1, idx2, idx3, idx4, idx5],
                                     [exp1, exp2, exp3, exp4, exp5]):
                for func in ['__repr__', '__unicode__', '__str__']:
                    result = getattr(idx, func)()
                    assert result == expected

    def test_representation_to_series(self):
        idx1 = TimedeltaIndex([], freq='D')
        idx2 = TimedeltaIndex(['1 days'], freq='D')
        idx3 = TimedeltaIndex(['1 days', '2 days'], freq='D')
        idx4 = TimedeltaIndex(['1 days', '2 days', '3 days'], freq='D')
        idx5 = TimedeltaIndex(['1 days 00:00:01', '2 days', '3 days'])

        exp1 = """Series([], dtype: timedelta64[ns])"""

        exp2 = """0   1 days
dtype: timedelta64[ns]"""

        exp3 = """0   1 days
1   2 days
dtype: timedelta64[ns]"""

        exp4 = """0   1 days
1   2 days
2   3 days
dtype: timedelta64[ns]"""

        exp5 = """0   1 days 00:00:01
1   2 days 00:00:00
2   3 days 00:00:00
dtype: timedelta64[ns]"""

        with pd.option_context('display.width', 300):
            for idx, expected in zip([idx1, idx2, idx3, idx4, idx5],
                                     [exp1, exp2, exp3, exp4, exp5]):
                result = repr(pd.Series(idx))
                assert result == expected

    def test_summary(self):
        # GH9116
        idx1 = TimedeltaIndex([], freq='D')
        idx2 = TimedeltaIndex(['1 days'], freq='D')
        idx3 = TimedeltaIndex(['1 days', '2 days'], freq='D')
        idx4 = TimedeltaIndex(['1 days', '2 days', '3 days'], freq='D')
        idx5 = TimedeltaIndex(['1 days 00:00:01', '2 days', '3 days'])

        exp1 = ("TimedeltaIndex: 0 entries\n"
                "Freq: D")

        exp2 = ("TimedeltaIndex: 1 entries, 1 days to 1 days\n"
                "Freq: D")

        exp3 = ("TimedeltaIndex: 2 entries, 1 days to 2 days\n"
                "Freq: D")

        exp4 = ("TimedeltaIndex: 3 entries, 1 days to 3 days\n"
                "Freq: D")

        exp5 = ("TimedeltaIndex: 3 entries, 1 days 00:00:01 to 3 days "
                "00:00:00")

        for idx, expected in zip([idx1, idx2, idx3, idx4, idx5],
                                 [exp1, exp2, exp3, exp4, exp5]):
            result = idx.summary()
            assert result == expected

    def test_comp_nat(self):
        left = pd.TimedeltaIndex([pd.Timedelta('1 days'), pd.NaT,
                                  pd.Timedelta('3 days')])
        right = pd.TimedeltaIndex([pd.NaT, pd.NaT, pd.Timedelta('3 days')])

        for lhs, rhs in [(left, right),
                         (left.astype(object), right.astype(object))]:
            result = rhs == lhs
            expected = np.array([False, False, True])
            tm.assert_numpy_array_equal(result, expected)

            result = rhs != lhs
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

        idx = timedelta_range('1 days 09:00:00', freq='H', periods=10)
        # create repeated values, 'n'th element is repeated by n+1 times
        idx = TimedeltaIndex(np.repeat(idx.values, range(1, len(idx) + 1)))

        exp_idx = timedelta_range('1 days 18:00:00', freq='-1H', periods=10)
        expected = Series(range(10, 0, -1), index=exp_idx, dtype='int64')

        for obj in [idx, Series(idx)]:
            tm.assert_series_equal(obj.value_counts(), expected)

        expected = timedelta_range('1 days 09:00:00', freq='H', periods=10)
        tm.assert_index_equal(idx.unique(), expected)

        idx = TimedeltaIndex(['1 days 09:00:00', '1 days 09:00:00',
                              '1 days 09:00:00', '1 days 08:00:00',
                              '1 days 08:00:00', pd.NaT])

        exp_idx = TimedeltaIndex(['1 days 09:00:00', '1 days 08:00:00'])
        expected = Series([3, 2], index=exp_idx)

        for obj in [idx, Series(idx)]:
            tm.assert_series_equal(obj.value_counts(), expected)

        exp_idx = TimedeltaIndex(['1 days 09:00:00', '1 days 08:00:00',
                                  pd.NaT])
        expected = Series([3, 2, 1], index=exp_idx)

        for obj in [idx, Series(idx)]:
            tm.assert_series_equal(obj.value_counts(dropna=False), expected)

        tm.assert_index_equal(idx.unique(), exp_idx)

    def test_nonunique_contains(self):
        # GH 9512
        for idx in map(TimedeltaIndex, ([0, 1, 0], [0, 0, -1], [0, -1, -1],
                                        ['00:01:00', '00:01:00', '00:02:00'],
                                        ['00:01:00', '00:01:00', '00:00:01'])):
            assert idx[0] in idx

    def test_unknown_attribute(self):
        # see gh-9680
        tdi = pd.timedelta_range(start=0, periods=10, freq='1s')
        ts = pd.Series(np.random.normal(size=10), index=tdi)
        assert 'foo' not in ts.__dict__.keys()
        pytest.raises(AttributeError, lambda: ts.foo)

    def test_order(self):
        # GH 10295
        idx1 = TimedeltaIndex(['1 day', '2 day', '3 day'], freq='D',
                              name='idx')
        idx2 = TimedeltaIndex(
            ['1 hour', '2 hour', '3 hour'], freq='H', name='idx')

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
            tm.assert_index_equal(ordered, idx[::-1])
            assert ordered.freq == expected.freq
            assert ordered.freq.n == -1

        idx1 = TimedeltaIndex(['1 hour', '3 hour', '5 hour',
                               '2 hour ', '1 hour'], name='idx1')
        exp1 = TimedeltaIndex(['1 hour', '1 hour', '2 hour',
                               '3 hour', '5 hour'], name='idx1')

        idx2 = TimedeltaIndex(['1 day', '3 day', '5 day',
                               '2 day', '1 day'], name='idx2')

        # TODO(wesm): unused?
        # exp2 = TimedeltaIndex(['1 day', '1 day', '2 day',
        #                        '3 day', '5 day'], name='idx2')

        # idx3 = TimedeltaIndex([pd.NaT, '3 minute', '5 minute',
        #                        '2 minute', pd.NaT], name='idx3')
        # exp3 = TimedeltaIndex([pd.NaT, pd.NaT, '2 minute', '3 minute',
        #                        '5 minute'], name='idx3')

        for idx, expected in [(idx1, exp1), (idx1, exp1), (idx1, exp1)]:
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
        idx = pd.timedelta_range('1 day', '31 day', freq='D', name='idx')
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
        base = pd.timedelta_range('1 day', '31 day', freq='D', name='idx')
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
        for freq in ['D', '3D', '-3D', 'H', '2H', '-2H', 'T', '2T', 'S', '-3S'
                     ]:
            idx = pd.timedelta_range('1', freq=freq, periods=10)
            result = pd.TimedeltaIndex(idx.asi8, freq='infer')
            tm.assert_index_equal(idx, result)
            assert result.freq == freq

    def test_nat_new(self):

        idx = pd.timedelta_range('1', freq='D', periods=5, name='x')
        result = idx._nat_new()
        exp = pd.TimedeltaIndex([pd.NaT] * 5, name='x')
        tm.assert_index_equal(result, exp)

        result = idx._nat_new(box=False)
        exp = np.array([iNaT] * 5, dtype=np.int64)
        tm.assert_numpy_array_equal(result, exp)

    def test_shift(self):
        # GH 9903
        idx = pd.TimedeltaIndex([], name='xxx')
        tm.assert_index_equal(idx.shift(0, freq='H'), idx)
        tm.assert_index_equal(idx.shift(3, freq='H'), idx)

        idx = pd.TimedeltaIndex(['5 hours', '6 hours', '9 hours'], name='xxx')
        tm.assert_index_equal(idx.shift(0, freq='H'), idx)
        exp = pd.TimedeltaIndex(['8 hours', '9 hours', '12 hours'], name='xxx')
        tm.assert_index_equal(idx.shift(3, freq='H'), exp)
        exp = pd.TimedeltaIndex(['2 hours', '3 hours', '6 hours'], name='xxx')
        tm.assert_index_equal(idx.shift(-3, freq='H'), exp)

        tm.assert_index_equal(idx.shift(0, freq='T'), idx)
        exp = pd.TimedeltaIndex(['05:03:00', '06:03:00', '9:03:00'],
                                name='xxx')
        tm.assert_index_equal(idx.shift(3, freq='T'), exp)
        exp = pd.TimedeltaIndex(['04:57:00', '05:57:00', '8:57:00'],
                                name='xxx')
        tm.assert_index_equal(idx.shift(-3, freq='T'), exp)

    def test_repeat(self):
        index = pd.timedelta_range('1 days', periods=2, freq='D')
        exp = pd.TimedeltaIndex(['1 days', '1 days', '2 days', '2 days'])
        for res in [index.repeat(2), np.repeat(index, 2)]:
            tm.assert_index_equal(res, exp)
            assert res.freq is None

        index = TimedeltaIndex(['1 days', 'NaT', '3 days'])
        exp = TimedeltaIndex(['1 days', '1 days', '1 days',
                              'NaT', 'NaT', 'NaT',
                              '3 days', '3 days', '3 days'])
        for res in [index.repeat(3), np.repeat(index, 3)]:
            tm.assert_index_equal(res, exp)
            assert res.freq is None

    def test_nat(self):
        assert pd.TimedeltaIndex._na_value is pd.NaT
        assert pd.TimedeltaIndex([])._na_value is pd.NaT

        idx = pd.TimedeltaIndex(['1 days', '2 days'])
        assert idx._can_hold_na

        tm.assert_numpy_array_equal(idx._isnan, np.array([False, False]))
        assert not idx.hasnans
        tm.assert_numpy_array_equal(idx._nan_idxs,
                                    np.array([], dtype=np.intp))

        idx = pd.TimedeltaIndex(['1 days', 'NaT'])
        assert idx._can_hold_na

        tm.assert_numpy_array_equal(idx._isnan, np.array([False, True]))
        assert idx.hasnans
        tm.assert_numpy_array_equal(idx._nan_idxs,
                                    np.array([1], dtype=np.intp))

    def test_equals(self):
        # GH 13107
        idx = pd.TimedeltaIndex(['1 days', '2 days', 'NaT'])
        assert idx.equals(idx)
        assert idx.equals(idx.copy())
        assert idx.equals(idx.astype(object))
        assert idx.astype(object).equals(idx)
        assert idx.astype(object).equals(idx.astype(object))
        assert not idx.equals(list(idx))
        assert not idx.equals(pd.Series(idx))

        idx2 = pd.TimedeltaIndex(['2 days', '1 days', 'NaT'])
        assert not idx.equals(idx2)
        assert not idx.equals(idx2.copy())
        assert not idx.equals(idx2.astype(object))
        assert not idx.astype(object).equals(idx2)
        assert not idx.astype(object).equals(idx2.astype(object))
        assert not idx.equals(list(idx2))
        assert not idx.equals(pd.Series(idx2))


class TestTimedeltas(object):
    _multiprocess_can_split_ = True

    def test_ops_error_str(self):
        # GH 13624
        tdi = TimedeltaIndex(['1 day', '2 days'])

        for l, r in [(tdi, 'a'), ('a', tdi)]:
            with pytest.raises(TypeError):
                l + r

            with pytest.raises(TypeError):
                l > r

            with pytest.raises(TypeError):
                l == r

            with pytest.raises(TypeError):
                l != r

    def test_timedelta_ops(self):
        # GH4984
        # make sure ops return Timedelta
        s = Series([Timestamp('20130101') + timedelta(seconds=i * i)
                    for i in range(10)])
        td = s.diff()

        result = td.mean()
        expected = to_timedelta(timedelta(seconds=9))
        assert result == expected

        result = td.to_frame().mean()
        assert result[0] == expected

        result = td.quantile(.1)
        expected = Timedelta(np.timedelta64(2600, 'ms'))
        assert result == expected

        result = td.median()
        expected = to_timedelta('00:00:09')
        assert result == expected

        result = td.to_frame().median()
        assert result[0] == expected

        # GH 6462
        # consistency in returned values for sum
        result = td.sum()
        expected = to_timedelta('00:01:21')
        assert result == expected

        result = td.to_frame().sum()
        assert result[0] == expected

        # std
        result = td.std()
        expected = to_timedelta(Series(td.dropna().values).std())
        assert result == expected

        result = td.to_frame().std()
        assert result[0] == expected

        # invalid ops
        for op in ['skew', 'kurt', 'sem', 'prod']:
            pytest.raises(TypeError, getattr(td, op))

        # GH 10040
        # make sure NaT is properly handled by median()
        s = Series([Timestamp('2015-02-03'), Timestamp('2015-02-07')])
        assert s.diff().median() == timedelta(days=4)

        s = Series([Timestamp('2015-02-03'), Timestamp('2015-02-07'),
                    Timestamp('2015-02-15')])
        assert s.diff().median() == timedelta(days=6)

    def test_compare_timedelta_series(self):
        # regresssion test for GH5963
        s = pd.Series([timedelta(days=1), timedelta(days=2)])
        actual = s > timedelta(days=1)
        expected = pd.Series([False, True])
        tm.assert_series_equal(actual, expected)

    def test_compare_timedelta_ndarray(self):
        # GH11835
        periods = [Timedelta('0 days 01:00:00'), Timedelta('0 days 01:00:00')]
        arr = np.array(periods)
        result = arr[0] > arr
        expected = np.array([False, False])
        tm.assert_numpy_array_equal(result, expected)
