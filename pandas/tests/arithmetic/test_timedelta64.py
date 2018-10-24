# -*- coding: utf-8 -*-
# Arithmetc tests for DataFrame/Series/Index/Array classes that should
# behave identically.
from datetime import datetime, timedelta
import operator

import pytest
import numpy as np

import pandas as pd
import pandas.util.testing as tm

from pandas.core import ops
from pandas.errors import NullFrequencyError, PerformanceWarning
from pandas import (
    timedelta_range,
    Timedelta, Timestamp, NaT, Series, TimedeltaIndex, DatetimeIndex,
    DataFrame)


# ------------------------------------------------------------------
# Timedelta64[ns] dtype Comparisons

class TestTimedelta64ArrayComparisons(object):
    # TODO: All of these need to be parametrized over box

    def test_compare_timedelta_series(self):
        # regresssion test for GH5963
        s = pd.Series([timedelta(days=1), timedelta(days=2)])
        actual = s > timedelta(days=1)
        expected = pd.Series([False, True])
        tm.assert_series_equal(actual, expected)

    def test_tdi_cmp_str_invalid(self):
        # GH#13624
        tdi = TimedeltaIndex(['1 day', '2 days'])

        for left, right in [(tdi, 'a'), ('a', tdi)]:
            with pytest.raises(TypeError):
                left > right

            with pytest.raises(TypeError):
                left == right

            with pytest.raises(TypeError):
                left != right

    @pytest.mark.parametrize('dtype', [None, object])
    def test_comp_nat(self, dtype):
        left = pd.TimedeltaIndex([pd.Timedelta('1 days'), pd.NaT,
                                  pd.Timedelta('3 days')])
        right = pd.TimedeltaIndex([pd.NaT, pd.NaT, pd.Timedelta('3 days')])

        lhs, rhs = left, right
        if dtype is object:
            lhs, rhs = left.astype(object), right.astype(object)

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

    def test_comparisons_nat(self):
        tdidx1 = pd.TimedeltaIndex(['1 day', pd.NaT, '1 day 00:00:01', pd.NaT,
                                    '1 day 00:00:01', '5 day 00:00:03'])
        tdidx2 = pd.TimedeltaIndex(['2 day', '2 day', pd.NaT, pd.NaT,
                                    '1 day 00:00:02', '5 days 00:00:03'])
        tdarr = np.array([np.timedelta64(2, 'D'),
                          np.timedelta64(2, 'D'), np.timedelta64('nat'),
                          np.timedelta64('nat'),
                          np.timedelta64(1, 'D') + np.timedelta64(2, 's'),
                          np.timedelta64(5, 'D') + np.timedelta64(3, 's')])

        cases = [(tdidx1, tdidx2), (tdidx1, tdarr)]

        # Check pd.NaT is handles as the same as np.nan
        for idx1, idx2 in cases:

            result = idx1 < idx2
            expected = np.array([True, False, False, False, True, False])
            tm.assert_numpy_array_equal(result, expected)

            result = idx2 > idx1
            expected = np.array([True, False, False, False, True, False])
            tm.assert_numpy_array_equal(result, expected)

            result = idx1 <= idx2
            expected = np.array([True, False, False, False, True, True])
            tm.assert_numpy_array_equal(result, expected)

            result = idx2 >= idx1
            expected = np.array([True, False, False, False, True, True])
            tm.assert_numpy_array_equal(result, expected)

            result = idx1 == idx2
            expected = np.array([False, False, False, False, False, True])
            tm.assert_numpy_array_equal(result, expected)

            result = idx1 != idx2
            expected = np.array([True, True, True, True, True, False])
            tm.assert_numpy_array_equal(result, expected)

    # TODO: better name
    def test_comparisons_coverage(self):
        rng = timedelta_range('1 days', periods=10)

        result = rng < rng[3]
        expected = np.array([True, True, True] + [False] * 7)
        tm.assert_numpy_array_equal(result, expected)

        # raise TypeError for now
        with pytest.raises(TypeError):
            rng < rng[3].value

        result = rng == list(rng)
        exp = rng == rng
        tm.assert_numpy_array_equal(result, exp)


# ------------------------------------------------------------------
# Timedelta64[ns] dtype Arithmetic Operations

class TestAddSubNaTMasking(object):
    # TODO: parametrize over boxes

    def test_tdi_add_timestamp_nat_masking(self):
        # GH#17991 checking for overflow-masking with NaT
        tdinat = pd.to_timedelta(['24658 days 11:15:00', 'NaT'])

        tsneg = Timestamp('1950-01-01')
        ts_neg_variants = [tsneg,
                           tsneg.to_pydatetime(),
                           tsneg.to_datetime64().astype('datetime64[ns]'),
                           tsneg.to_datetime64().astype('datetime64[D]')]

        tspos = Timestamp('1980-01-01')
        ts_pos_variants = [tspos,
                           tspos.to_pydatetime(),
                           tspos.to_datetime64().astype('datetime64[ns]'),
                           tspos.to_datetime64().astype('datetime64[D]')]

        for variant in ts_neg_variants + ts_pos_variants:
            res = tdinat + variant
            assert res[1] is pd.NaT

    def test_tdi_add_overflow(self):
        # See GH#14068
        msg = "too (big|large) to convert"
        with tm.assert_raises_regex(OverflowError, msg):
            pd.to_timedelta(106580, 'D') + Timestamp('2000')
        with tm.assert_raises_regex(OverflowError, msg):
            Timestamp('2000') + pd.to_timedelta(106580, 'D')

        _NaT = int(pd.NaT) + 1
        msg = "Overflow in int64 addition"
        with tm.assert_raises_regex(OverflowError, msg):
            pd.to_timedelta([106580], 'D') + Timestamp('2000')
        with tm.assert_raises_regex(OverflowError, msg):
            Timestamp('2000') + pd.to_timedelta([106580], 'D')
        with tm.assert_raises_regex(OverflowError, msg):
            pd.to_timedelta([_NaT]) - Timedelta('1 days')
        with tm.assert_raises_regex(OverflowError, msg):
            pd.to_timedelta(['5 days', _NaT]) - Timedelta('1 days')
        with tm.assert_raises_regex(OverflowError, msg):
            (pd.to_timedelta([_NaT, '5 days', '1 hours']) -
             pd.to_timedelta(['7 seconds', _NaT, '4 hours']))

        # These should not overflow!
        exp = TimedeltaIndex([pd.NaT])
        result = pd.to_timedelta([pd.NaT]) - Timedelta('1 days')
        tm.assert_index_equal(result, exp)

        exp = TimedeltaIndex(['4 days', pd.NaT])
        result = pd.to_timedelta(['5 days', pd.NaT]) - Timedelta('1 days')
        tm.assert_index_equal(result, exp)

        exp = TimedeltaIndex([pd.NaT, pd.NaT, '5 hours'])
        result = (pd.to_timedelta([pd.NaT, '5 days', '1 hours']) +
                  pd.to_timedelta(['7 seconds', pd.NaT, '4 hours']))
        tm.assert_index_equal(result, exp)


class TestTimedeltaArraylikeAddSubOps(object):
    # Tests for timedelta64[ns] __add__, __sub__, __radd__, __rsub__

    # TODO: moved from tests.series.test_operators, needs splitting, cleanup,
    # de-duplication, box-parametrization...
    def test_operators_timedelta64(self):
        # series ops
        v1 = pd.date_range('2012-1-1', periods=3, freq='D')
        v2 = pd.date_range('2012-1-2', periods=3, freq='D')
        rs = Series(v2) - Series(v1)
        xp = Series(1e9 * 3600 * 24,
                    rs.index).astype('int64').astype('timedelta64[ns]')
        tm.assert_series_equal(rs, xp)
        assert rs.dtype == 'timedelta64[ns]'

        df = DataFrame(dict(A=v1))
        td = Series([timedelta(days=i) for i in range(3)])
        assert td.dtype == 'timedelta64[ns]'

        # series on the rhs
        result = df['A'] - df['A'].shift()
        assert result.dtype == 'timedelta64[ns]'

        result = df['A'] + td
        assert result.dtype == 'M8[ns]'

        # scalar Timestamp on rhs
        maxa = df['A'].max()
        assert isinstance(maxa, Timestamp)

        resultb = df['A'] - df['A'].max()
        assert resultb.dtype == 'timedelta64[ns]'

        # timestamp on lhs
        result = resultb + df['A']
        values = [Timestamp('20111230'), Timestamp('20120101'),
                  Timestamp('20120103')]
        expected = Series(values, name='A')
        tm.assert_series_equal(result, expected)

        # datetimes on rhs
        result = df['A'] - datetime(2001, 1, 1)
        expected = Series(
            [timedelta(days=4017 + i) for i in range(3)], name='A')
        tm.assert_series_equal(result, expected)
        assert result.dtype == 'm8[ns]'

        d = datetime(2001, 1, 1, 3, 4)
        resulta = df['A'] - d
        assert resulta.dtype == 'm8[ns]'

        # roundtrip
        resultb = resulta + d
        tm.assert_series_equal(df['A'], resultb)

        # timedeltas on rhs
        td = timedelta(days=1)
        resulta = df['A'] + td
        resultb = resulta - td
        tm.assert_series_equal(resultb, df['A'])
        assert resultb.dtype == 'M8[ns]'

        # roundtrip
        td = timedelta(minutes=5, seconds=3)
        resulta = df['A'] + td
        resultb = resulta - td
        tm.assert_series_equal(df['A'], resultb)
        assert resultb.dtype == 'M8[ns]'

        # inplace
        value = rs[2] + np.timedelta64(timedelta(minutes=5, seconds=1))
        rs[2] += np.timedelta64(timedelta(minutes=5, seconds=1))
        assert rs[2] == value

    def test_timedelta64_ops_nat(self):
        # GH 11349
        timedelta_series = Series([NaT, Timedelta('1s')])
        nat_series_dtype_timedelta = Series([NaT, NaT],
                                            dtype='timedelta64[ns]')
        single_nat_dtype_timedelta = Series([NaT], dtype='timedelta64[ns]')

        # subtraction
        tm.assert_series_equal(timedelta_series - NaT,
                               nat_series_dtype_timedelta)
        tm.assert_series_equal(-NaT + timedelta_series,
                               nat_series_dtype_timedelta)

        tm.assert_series_equal(timedelta_series - single_nat_dtype_timedelta,
                               nat_series_dtype_timedelta)
        tm.assert_series_equal(-single_nat_dtype_timedelta + timedelta_series,
                               nat_series_dtype_timedelta)

        # addition
        tm.assert_series_equal(nat_series_dtype_timedelta + NaT,
                               nat_series_dtype_timedelta)
        tm.assert_series_equal(NaT + nat_series_dtype_timedelta,
                               nat_series_dtype_timedelta)

        tm.assert_series_equal(nat_series_dtype_timedelta +
                               single_nat_dtype_timedelta,
                               nat_series_dtype_timedelta)
        tm.assert_series_equal(single_nat_dtype_timedelta +
                               nat_series_dtype_timedelta,
                               nat_series_dtype_timedelta)

        tm.assert_series_equal(timedelta_series + NaT,
                               nat_series_dtype_timedelta)
        tm.assert_series_equal(NaT + timedelta_series,
                               nat_series_dtype_timedelta)

        tm.assert_series_equal(timedelta_series + single_nat_dtype_timedelta,
                               nat_series_dtype_timedelta)
        tm.assert_series_equal(single_nat_dtype_timedelta + timedelta_series,
                               nat_series_dtype_timedelta)

        tm.assert_series_equal(nat_series_dtype_timedelta + NaT,
                               nat_series_dtype_timedelta)
        tm.assert_series_equal(NaT + nat_series_dtype_timedelta,
                               nat_series_dtype_timedelta)

        tm.assert_series_equal(nat_series_dtype_timedelta +
                               single_nat_dtype_timedelta,
                               nat_series_dtype_timedelta)
        tm.assert_series_equal(single_nat_dtype_timedelta +
                               nat_series_dtype_timedelta,
                               nat_series_dtype_timedelta)

        # multiplication
        tm.assert_series_equal(nat_series_dtype_timedelta * 1.0,
                               nat_series_dtype_timedelta)
        tm.assert_series_equal(1.0 * nat_series_dtype_timedelta,
                               nat_series_dtype_timedelta)

        tm.assert_series_equal(timedelta_series * 1, timedelta_series)
        tm.assert_series_equal(1 * timedelta_series, timedelta_series)

        tm.assert_series_equal(timedelta_series * 1.5,
                               Series([NaT, Timedelta('1.5s')]))
        tm.assert_series_equal(1.5 * timedelta_series,
                               Series([NaT, Timedelta('1.5s')]))

        tm.assert_series_equal(timedelta_series * np.nan,
                               nat_series_dtype_timedelta)
        tm.assert_series_equal(np.nan * timedelta_series,
                               nat_series_dtype_timedelta)

        # division
        tm.assert_series_equal(timedelta_series / 2,
                               Series([NaT, Timedelta('0.5s')]))
        tm.assert_series_equal(timedelta_series / 2.0,
                               Series([NaT, Timedelta('0.5s')]))
        tm.assert_series_equal(timedelta_series / np.nan,
                               nat_series_dtype_timedelta)

    # -------------------------------------------------------------
    # Invalid Operations

    def test_td64arr_add_str_invalid(self, box):
        # GH#13624
        tdi = TimedeltaIndex(['1 day', '2 days'])
        tdi = tm.box_expected(tdi, box)

        with pytest.raises(TypeError):
            tdi + 'a'
        with pytest.raises(TypeError):
            'a' + tdi

    @pytest.mark.parametrize('other', [3.14, np.array([2.0, 3.0])])
    @pytest.mark.parametrize('op', [operator.add, ops.radd,
                                    operator.sub, ops.rsub],
                             ids=lambda x: x.__name__)
    def test_td64arr_add_sub_float(self, box, op, other):
        tdi = TimedeltaIndex(['-1 days', '-1 days'])
        tdi = tm.box_expected(tdi, box, transpose=True)

        with pytest.raises(TypeError):
            op(tdi, other)

    @pytest.mark.parametrize('freq', [None, 'H'])
    def test_td64arr_sub_period(self, box, freq):
        # GH#13078
        # not supported, check TypeError
        p = pd.Period('2011-01-01', freq='D')
        idx = TimedeltaIndex(['1 hours', '2 hours'], freq=freq)
        idx = tm.box_expected(idx, box)

        with pytest.raises(TypeError):
            idx - p

        with pytest.raises(TypeError):
            p - idx

    @pytest.mark.parametrize('pi_freq', ['D', 'W', 'Q', 'H'])
    @pytest.mark.parametrize('tdi_freq', [None, 'H'])
    def test_td64arr_sub_pi(self, box, tdi_freq, pi_freq):
        # GH#20049 subtracting PeriodIndex should raise TypeError
        tdi = TimedeltaIndex(['1 hours', '2 hours'], freq=tdi_freq)
        dti = Timestamp('2018-03-07 17:16:40') + tdi
        pi = dti.to_period(pi_freq)

        # TODO: parametrize over box for pi?
        tdi = tm.box_expected(tdi, box, transpose=True)
        with pytest.raises(TypeError):
            tdi - pi

    # -------------------------------------------------------------
    # Binary operations td64 arraylike and datetime-like

    def test_td64arr_sub_timestamp_raises(self, box):
        idx = TimedeltaIndex(['1 day', '2 day'])
        idx = tm.box_expected(idx, box)

        msg = ("cannot subtract a datelike from|"
               "Could not operate|"
               "cannot perform operation")
        with tm.assert_raises_regex(TypeError, msg):
            idx - Timestamp('2011-01-01')

    def test_td64arr_add_timestamp(self, box):
        idx = TimedeltaIndex(['1 day', '2 day'])
        expected = DatetimeIndex(['2011-01-02', '2011-01-03'])

        idx = tm.box_expected(idx, box)
        expected = tm.box_expected(expected, box)

        result = idx + Timestamp('2011-01-01')
        tm.assert_equal(result, expected)

        result = Timestamp('2011-01-01') + idx
        tm.assert_equal(result, expected)

    def test_td64arr_add_sub_timestamp(self, box):
        # GH#11925
        ts = Timestamp('2012-01-01')
        # TODO: parametrize over types of datetime scalar?

        tdser = Series(timedelta_range('1 day', periods=3))
        expected = Series(pd.date_range('2012-01-02', periods=3))

        tdser = tm.box_expected(tdser, box)
        expected = tm.box_expected(expected, box)

        tm.assert_equal(ts + tdser, expected)
        tm.assert_equal(tdser + ts, expected)

        expected2 = Series(pd.date_range('2011-12-31',
                                         periods=3, freq='-1D'))
        expected2 = tm.box_expected(expected2, box)

        tm.assert_equal(ts - tdser, expected2)
        tm.assert_equal(ts + (-tdser), expected2)

        with pytest.raises(TypeError):
            tdser - ts

    def test_tdi_sub_dt64_array(self, box):
        dti = pd.date_range('2016-01-01', periods=3)
        tdi = dti - dti.shift(1)
        dtarr = dti.values
        expected = pd.DatetimeIndex(dtarr) - tdi

        tdi = tm.box_expected(tdi, box, transpose=True)
        expected = tm.box_expected(expected, box, transpose=True)

        with pytest.raises(TypeError):
            tdi - dtarr

        # TimedeltaIndex.__rsub__
        result = dtarr - tdi
        tm.assert_equal(result, expected)

    def test_tdi_add_dt64_array(self, box):
        dti = pd.date_range('2016-01-01', periods=3)
        tdi = dti - dti.shift(1)
        dtarr = dti.values
        expected = pd.DatetimeIndex(dtarr) + tdi

        tdi = tm.box_expected(tdi, box, transpose=True)
        expected = tm.box_expected(expected, box, transpose=True)

        result = tdi + dtarr
        tm.assert_equal(result, expected)
        result = dtarr + tdi
        tm.assert_equal(result, expected)

    # ------------------------------------------------------------------
    # Operations with int-like others

    def test_td64arr_add_int_series_invalid(self, box, tdser):
        tdser = tm.box_expected(tdser, box)
        err = TypeError if box is not pd.Index else NullFrequencyError
        int_ser = Series([2, 3, 4])

        with pytest.raises(err):
            tdser + int_ser
        with pytest.raises(err):
            int_ser + tdser
        with pytest.raises(err):
            tdser - int_ser
        with pytest.raises(err):
            int_ser - tdser

    def test_td64arr_add_intlike(self, box):
        # GH#19123
        tdi = TimedeltaIndex(['59 days', '59 days', 'NaT'])
        ser = tm.box_expected(tdi, box, transpose=True)
        err = TypeError if box is not pd.Index else NullFrequencyError

        other = Series([20, 30, 40], dtype='uint8')

        # TODO: separate/parametrize
        with pytest.raises(err):
            ser + 1
        with pytest.raises(err):
            ser - 1

        with pytest.raises(err):
            ser + other
        with pytest.raises(err):
            ser - other

        with pytest.raises(err):
            ser + np.array(other)
        with pytest.raises(err):
            ser - np.array(other)

        with pytest.raises(err):
            ser + pd.Index(other)
        with pytest.raises(err):
            ser - pd.Index(other)

    @pytest.mark.parametrize('scalar', [1, 1.5, np.array(2)])
    def test_td64arr_add_sub_numeric_scalar_invalid(self, box, scalar, tdser):
        tdser = tm.box_expected(tdser, box)
        err = TypeError
        if box is pd.Index and not isinstance(scalar, float):
            err = NullFrequencyError

        with pytest.raises(err):
            tdser + scalar
        with pytest.raises(err):
            scalar + tdser
        with pytest.raises(err):
            tdser - scalar
        with pytest.raises(err):
            scalar - tdser

    @pytest.mark.parametrize('dtype', ['int64', 'int32', 'int16',
                                       'uint64', 'uint32', 'uint16', 'uint8',
                                       'float64', 'float32', 'float16'])
    @pytest.mark.parametrize('vec', [
        np.array([1, 2, 3]),
        pd.Index([1, 2, 3]),
        Series([1, 2, 3])
        # TODO: Add DataFrame in here?
    ], ids=lambda x: type(x).__name__)
    def test_td64arr_add_sub_numeric_arr_invalid(self, box, vec, dtype, tdser):
        tdser = tm.box_expected(tdser, box, transpose=True)
        err = TypeError
        if box is pd.Index and not dtype.startswith('float'):
            err = NullFrequencyError

        vector = vec.astype(dtype)
        # TODO: parametrize over these four ops?
        with pytest.raises(err):
            tdser + vector
        with pytest.raises(err):
            vector + tdser
        with pytest.raises(err):
            tdser - vector
        with pytest.raises(err):
            vector - tdser

    # ------------------------------------------------------------------
    # Operations with timedelta-like others

    # TODO: this was taken from tests.series.test_ops; de-duplicate
    @pytest.mark.parametrize('scalar_td', [timedelta(minutes=5, seconds=4),
                                           Timedelta(minutes=5, seconds=4),
                                           Timedelta('5m4s').to_timedelta64()])
    def test_operators_timedelta64_with_timedelta(self, scalar_td):
        # smoke tests
        td1 = Series([timedelta(minutes=5, seconds=3)] * 3)
        td1.iloc[2] = np.nan

        td1 + scalar_td
        scalar_td + td1
        td1 - scalar_td
        scalar_td - td1
        td1 / scalar_td
        scalar_td / td1

    # TODO: this was taken from tests.series.test_ops; de-duplicate
    def test_timedelta64_operations_with_timedeltas(self):
        # td operate with td
        td1 = Series([timedelta(minutes=5, seconds=3)] * 3)
        td2 = timedelta(minutes=5, seconds=4)
        result = td1 - td2
        expected = (Series([timedelta(seconds=0)] * 3) -
                    Series([timedelta(seconds=1)] * 3))
        assert result.dtype == 'm8[ns]'
        tm.assert_series_equal(result, expected)

        result2 = td2 - td1
        expected = (Series([timedelta(seconds=1)] * 3) -
                    Series([timedelta(seconds=0)] * 3))
        tm.assert_series_equal(result2, expected)

        # roundtrip
        tm.assert_series_equal(result + td2, td1)

        # Now again, using pd.to_timedelta, which should build
        # a Series or a scalar, depending on input.
        td1 = Series(pd.to_timedelta(['00:05:03'] * 3))
        td2 = pd.to_timedelta('00:05:04')
        result = td1 - td2
        expected = (Series([timedelta(seconds=0)] * 3) -
                    Series([timedelta(seconds=1)] * 3))
        assert result.dtype == 'm8[ns]'
        tm.assert_series_equal(result, expected)

        result2 = td2 - td1
        expected = (Series([timedelta(seconds=1)] * 3) -
                    Series([timedelta(seconds=0)] * 3))
        tm.assert_series_equal(result2, expected)

        # roundtrip
        tm.assert_series_equal(result + td2, td1)

    def test_td64arr_add_td64_array(self, box):
        dti = pd.date_range('2016-01-01', periods=3)
        tdi = dti - dti.shift(1)
        tdarr = tdi.values

        expected = 2 * tdi
        tdi = tm.box_expected(tdi, box, transpose=True)
        expected = tm.box_expected(expected, box, transpose=True)

        result = tdi + tdarr
        tm.assert_equal(result, expected)
        result = tdarr + tdi
        tm.assert_equal(result, expected)

    def test_td64arr_sub_td64_array(self, box):
        dti = pd.date_range('2016-01-01', periods=3)
        tdi = dti - dti.shift(1)
        tdarr = tdi.values

        expected = 0 * tdi
        tdi = tm.box_expected(tdi, box, transpose=True)
        expected = tm.box_expected(expected, box, transpose=True)

        result = tdi - tdarr
        tm.assert_equal(result, expected)
        result = tdarr - tdi
        tm.assert_equal(result, expected)

    # TODO: parametrize over [add, sub, radd, rsub]?
    @pytest.mark.parametrize('names', [(None, None, None),
                                       ('Egon', 'Venkman', None),
                                       ('NCC1701D', 'NCC1701D', 'NCC1701D')])
    def test_td64arr_add_sub_tdi(self, box, names):
        # GH#17250 make sure result dtype is correct
        # GH#19043 make sure names are propagated correctly
        if box is pd.DataFrame and names[0] != names[1]:
            return

        tdi = TimedeltaIndex(['0 days', '1 day'], name=names[0])
        ser = Series([Timedelta(hours=3), Timedelta(hours=4)], name=names[1])
        expected = Series([Timedelta(hours=3), Timedelta(days=1, hours=4)],
                          name=names[2])

        ser = tm.box_expected(ser, box, transpose=True)
        expected = tm.box_expected(expected, box, transpose=True)

        result = tdi + ser
        tm.assert_equal(result, expected)
        if box is not pd.DataFrame:
            assert result.dtype == 'timedelta64[ns]'
        else:
            assert result.dtypes[0] == 'timedelta64[ns]'

        result = ser + tdi
        tm.assert_equal(result, expected)
        if box is not pd.DataFrame:
            assert result.dtype == 'timedelta64[ns]'
        else:
            assert result.dtypes[0] == 'timedelta64[ns]'

        expected = Series([Timedelta(hours=-3), Timedelta(days=1, hours=-4)],
                          name=names[2])
        expected = tm.box_expected(expected, box, transpose=True)

        result = tdi - ser
        tm.assert_equal(result, expected)
        if box is not pd.DataFrame:
            assert result.dtype == 'timedelta64[ns]'
        else:
            assert result.dtypes[0] == 'timedelta64[ns]'

        result = ser - tdi
        tm.assert_equal(result, -expected)
        if box is not pd.DataFrame:
            assert result.dtype == 'timedelta64[ns]'
        else:
            assert result.dtypes[0] == 'timedelta64[ns]'

    def test_td64arr_sub_NaT(self, box):
        # GH#18808
        ser = Series([NaT, Timedelta('1s')])
        expected = Series([NaT, NaT], dtype='timedelta64[ns]')

        ser = tm.box_expected(ser, box)
        expected = tm.box_expected(expected, box)

        res = ser - pd.NaT
        tm.assert_equal(res, expected)

    def test_td64arr_add_timedeltalike(self, two_hours, box):
        # only test adding/sub offsets as + is now numeric
        rng = timedelta_range('1 days', '10 days')
        expected = timedelta_range('1 days 02:00:00', '10 days 02:00:00',
                                   freq='D')
        rng = tm.box_expected(rng, box)
        expected = tm.box_expected(expected, box)

        result = rng + two_hours
        tm.assert_equal(result, expected)

    def test_td64arr_sub_timedeltalike(self, two_hours, box):
        # only test adding/sub offsets as - is now numeric
        rng = timedelta_range('1 days', '10 days')
        expected = timedelta_range('0 days 22:00:00', '9 days 22:00:00')

        rng = tm.box_expected(rng, box)
        expected = tm.box_expected(expected, box)

        result = rng - two_hours
        tm.assert_equal(result, expected)

    # ------------------------------------------------------------------
    # __add__/__sub__ with DateOffsets and arrays of DateOffsets

    # TODO: this was taken from tests.series.test_operators; de-duplicate
    def test_timedelta64_operations_with_DateOffset(self):
        # GH#10699
        td = Series([timedelta(minutes=5, seconds=3)] * 3)
        result = td + pd.offsets.Minute(1)
        expected = Series([timedelta(minutes=6, seconds=3)] * 3)
        tm.assert_series_equal(result, expected)

        result = td - pd.offsets.Minute(1)
        expected = Series([timedelta(minutes=4, seconds=3)] * 3)
        tm.assert_series_equal(result, expected)

        with tm.assert_produces_warning(PerformanceWarning):
            result = td + Series([pd.offsets.Minute(1), pd.offsets.Second(3),
                                  pd.offsets.Hour(2)])
        expected = Series([timedelta(minutes=6, seconds=3),
                           timedelta(minutes=5, seconds=6),
                           timedelta(hours=2, minutes=5, seconds=3)])
        tm.assert_series_equal(result, expected)

        result = td + pd.offsets.Minute(1) + pd.offsets.Second(12)
        expected = Series([timedelta(minutes=6, seconds=15)] * 3)
        tm.assert_series_equal(result, expected)

        # valid DateOffsets
        for do in ['Hour', 'Minute', 'Second', 'Day', 'Micro', 'Milli',
                   'Nano']:
            op = getattr(pd.offsets, do)
            td + op(5)
            op(5) + td
            td - op(5)
            op(5) - td

    @pytest.mark.parametrize('names', [(None, None, None),
                                       ('foo', 'bar', None),
                                       ('foo', 'foo', 'foo')])
    def test_td64arr_add_offset_index(self, names, box):
        # GH#18849, GH#19744
        if box is pd.DataFrame and names[0] != names[1]:
            return

        tdi = TimedeltaIndex(['1 days 00:00:00', '3 days 04:00:00'],
                             name=names[0])
        other = pd.Index([pd.offsets.Hour(n=1), pd.offsets.Minute(n=-2)],
                         name=names[1])

        expected = TimedeltaIndex([tdi[n] + other[n] for n in range(len(tdi))],
                                  freq='infer', name=names[2])
        tdi = tm.box_expected(tdi, box, transpose=True)
        expected = tm.box_expected(expected, box, transpose=True)

        exwarning = PerformanceWarning if box is not pd.DataFrame else None
        with tm.assert_produces_warning(exwarning):
            res = tdi + other
        tm.assert_equal(res, expected)

        with tm.assert_produces_warning(exwarning):
            res2 = other + tdi
        tm.assert_equal(res2, expected)

    # TODO: combine with test_td64arr_add_offset_index by parametrizing
    # over second box?
    def test_td64arr_add_offset_array(self, box):
        # GH#18849
        tdi = TimedeltaIndex(['1 days 00:00:00', '3 days 04:00:00'])
        other = np.array([pd.offsets.Hour(n=1), pd.offsets.Minute(n=-2)])

        expected = TimedeltaIndex([tdi[n] + other[n] for n in range(len(tdi))],
                                  freq='infer')

        tdi = tm.box_expected(tdi, box, transpose=True)
        expected = tm.box_expected(expected, box, transpose=True)

        exwarning = PerformanceWarning if box is not pd.DataFrame else None
        with tm.assert_produces_warning(exwarning):
            res = tdi + other
        tm.assert_equal(res, expected)

        with tm.assert_produces_warning(exwarning):
            res2 = other + tdi
        tm.assert_equal(res2, expected)

    @pytest.mark.parametrize('names', [(None, None, None),
                                       ('foo', 'bar', None),
                                       ('foo', 'foo', 'foo')])
    def test_td64arr_sub_offset_index(self, names, box):
        # GH#18824, GH#19744
        if box is pd.DataFrame and names[0] != names[1]:
            return

        tdi = TimedeltaIndex(['1 days 00:00:00', '3 days 04:00:00'],
                             name=names[0])
        other = pd.Index([pd.offsets.Hour(n=1), pd.offsets.Minute(n=-2)],
                         name=names[1])

        expected = TimedeltaIndex([tdi[n] - other[n] for n in range(len(tdi))],
                                  freq='infer', name=names[2])

        tdi = tm.box_expected(tdi, box, transpose=True)
        expected = tm.box_expected(expected, box, transpose=True)

        exwarning = PerformanceWarning if box is not pd.DataFrame else None
        with tm.assert_produces_warning(exwarning):
            res = tdi - other
        tm.assert_equal(res, expected)

    def test_td64arr_sub_offset_array(self, box):
        # GH#18824
        tdi = TimedeltaIndex(['1 days 00:00:00', '3 days 04:00:00'])
        other = np.array([pd.offsets.Hour(n=1), pd.offsets.Minute(n=-2)])

        expected = TimedeltaIndex([tdi[n] - other[n] for n in range(len(tdi))],
                                  freq='infer')

        tdi = tm.box_expected(tdi, box, transpose=True)
        expected = tm.box_expected(expected, box, transpose=True)

        exwarning = PerformanceWarning if box is not pd.DataFrame else None
        with tm.assert_produces_warning(exwarning):
            res = tdi - other
        tm.assert_equal(res, expected)

    @pytest.mark.parametrize('names', [(None, None, None),
                                       ('foo', 'bar', None),
                                       ('foo', 'foo', 'foo')])
    def test_td64arr_with_offset_series(self, names, box):
        # GH#18849
        if box is pd.DataFrame and names[0] != names[1]:
            return

        box2 = Series if box is pd.Index else box
        exwarning = PerformanceWarning if box is not pd.DataFrame else None

        tdi = TimedeltaIndex(['1 days 00:00:00', '3 days 04:00:00'],
                             name=names[0])
        other = Series([pd.offsets.Hour(n=1), pd.offsets.Minute(n=-2)],
                       name=names[1])

        expected_add = Series([tdi[n] + other[n] for n in range(len(tdi))],
                              name=names[2])
        expected_sub = Series([tdi[n] - other[n] for n in range(len(tdi))],
                              name=names[2])

        tdi = tm.box_expected(tdi, box, transpose=True)
        expected_add = tm.box_expected(expected_add, box2, transpose=True)

        with tm.assert_produces_warning(exwarning):
            res = tdi + other
        tm.assert_equal(res, expected_add)

        with tm.assert_produces_warning(exwarning):
            res2 = other + tdi
        tm.assert_equal(res2, expected_add)

        expected_sub = tm.box_expected(expected_sub, box2, transpose=True)

        with tm.assert_produces_warning(exwarning):
            res3 = tdi - other
        tm.assert_equal(res3, expected_sub)

    @pytest.mark.parametrize('obox', [np.array, pd.Index, pd.Series])
    def test_td64arr_addsub_anchored_offset_arraylike(self, obox, box):
        # GH#18824
        tdi = TimedeltaIndex(['1 days 00:00:00', '3 days 04:00:00'])
        tdi = tm.box_expected(tdi, box, transpose=True)

        anchored = obox([pd.offsets.MonthEnd(), pd.offsets.Day(n=2)])

        # addition/subtraction ops with anchored offsets should issue
        # a PerformanceWarning and _then_ raise a TypeError.
        with pytest.raises(TypeError):
            with tm.assert_produces_warning(PerformanceWarning):
                tdi + anchored
        with pytest.raises(TypeError):
            with tm.assert_produces_warning(PerformanceWarning):
                anchored + tdi
        with pytest.raises(TypeError):
            with tm.assert_produces_warning(PerformanceWarning):
                tdi - anchored
        with pytest.raises(TypeError):
            with tm.assert_produces_warning(PerformanceWarning):
                anchored - tdi


class TestTimedeltaArraylikeMulDivOps(object):
    # Tests for timedelta64[ns]
    # __mul__, __rmul__, __div__, __rdiv__, __floordiv__, __rfloordiv__

    # TODO: Moved from tests.series.test_operators; needs cleanup
    @pytest.mark.parametrize("m", [1, 3, 10])
    @pytest.mark.parametrize("unit", ['D', 'h', 'm', 's', 'ms', 'us', 'ns'])
    def test_timedelta64_conversions(self, m, unit):
        startdate = Series(pd.date_range('2013-01-01', '2013-01-03'))
        enddate = Series(pd.date_range('2013-03-01', '2013-03-03'))

        ser = enddate - startdate
        ser[2] = np.nan

        # op
        expected = Series([x / np.timedelta64(m, unit) for x in ser])
        result = ser / np.timedelta64(m, unit)
        tm.assert_series_equal(result, expected)

        # reverse op
        expected = Series([Timedelta(np.timedelta64(m, unit)) / x
                           for x in ser])
        result = np.timedelta64(m, unit) / ser
        tm.assert_series_equal(result, expected)

    # ------------------------------------------------------------------
    # Multiplication
    # organized with scalar others first, then array-like

    def test_td64arr_mul_int(self, box):
        idx = TimedeltaIndex(np.arange(5, dtype='int64'))
        idx = tm.box_expected(idx, box)

        result = idx * 1
        tm.assert_equal(result, idx)

        result = 1 * idx
        tm.assert_equal(result, idx)

    def test_td64arr_mul_tdlike_scalar_raises(self, two_hours, box):
        rng = timedelta_range('1 days', '10 days', name='foo')
        rng = tm.box_expected(rng, box)
        with pytest.raises(TypeError):
            rng * two_hours

    def test_tdi_mul_int_array_zerodim(self, box):
        rng5 = np.arange(5, dtype='int64')
        idx = TimedeltaIndex(rng5)
        expected = TimedeltaIndex(rng5 * 5)

        idx = tm.box_expected(idx, box)
        expected = tm.box_expected(expected, box)

        result = idx * np.array(5, dtype='int64')
        tm.assert_equal(result, expected)

    def test_tdi_mul_int_array(self, box):
        rng5 = np.arange(5, dtype='int64')
        idx = TimedeltaIndex(rng5)
        expected = TimedeltaIndex(rng5 ** 2)

        idx = tm.box_expected(idx, box, transpose=True)
        expected = tm.box_expected(expected, box, transpose=True)

        result = idx * rng5
        tm.assert_equal(result, expected)

    def test_tdi_mul_int_series(self, box):
        idx = TimedeltaIndex(np.arange(5, dtype='int64'))
        expected = TimedeltaIndex(np.arange(5, dtype='int64') ** 2)

        idx = tm.box_expected(idx, box, transpose=True)

        box2 = pd.Series if box is pd.Index else box
        expected = tm.box_expected(expected, box2, transpose=True)

        result = idx * pd.Series(np.arange(5, dtype='int64'))
        tm.assert_equal(result, expected)

    def test_tdi_mul_float_series(self, box):
        idx = TimedeltaIndex(np.arange(5, dtype='int64'))
        idx = tm.box_expected(idx, box, transpose=True)

        rng5f = np.arange(5, dtype='float64')
        expected = TimedeltaIndex(rng5f * (rng5f + 0.1))
        box2 = pd.Series if box is pd.Index else box
        expected = tm.box_expected(expected, box2, transpose=True)

        result = idx * Series(rng5f + 0.1)
        tm.assert_equal(result, expected)

    # TODO: Put Series/DataFrame in others?
    @pytest.mark.parametrize('other', [
        np.arange(1, 11),
        pd.Int64Index(range(1, 11)),
        pd.UInt64Index(range(1, 11)),
        pd.Float64Index(range(1, 11)),
        pd.RangeIndex(1, 11)
    ], ids=lambda x: type(x).__name__)
    def test_tdi_rmul_arraylike(self, other, box):
        tdi = TimedeltaIndex(['1 Day'] * 10)
        expected = timedelta_range('1 days', '10 days')

        tdi = tm.box_expected(tdi, box, transpose=True)
        expected = tm.box_expected(expected, box, transpose=True)

        result = other * tdi
        tm.assert_equal(result, expected)
        commute = tdi * other
        tm.assert_equal(commute, expected)

    # ------------------------------------------------------------------
    # __div__

    def test_td64arr_div_nat_invalid(self, box):
        # don't allow division by NaT (maybe could in the future)
        rng = timedelta_range('1 days', '10 days', name='foo')
        rng = tm.box_expected(rng, box)
        with pytest.raises(TypeError):
            rng / pd.NaT

    def test_td64arr_div_int(self, box):
        idx = TimedeltaIndex(np.arange(5, dtype='int64'))
        idx = tm.box_expected(idx, box)

        result = idx / 1
        tm.assert_equal(result, idx)

    def test_tdi_div_tdlike_scalar(self, two_hours, box):
        # GH#20088, GH#22163 ensure DataFrame returns correct dtype
        rng = timedelta_range('1 days', '10 days', name='foo')
        expected = pd.Float64Index((np.arange(10) + 1) * 12, name='foo')

        rng = tm.box_expected(rng, box)
        expected = tm.box_expected(expected, box)

        result = rng / two_hours
        tm.assert_equal(result, expected)

    def test_tdi_div_tdlike_scalar_with_nat(self, two_hours, box):
        rng = TimedeltaIndex(['1 days', pd.NaT, '2 days'], name='foo')
        expected = pd.Float64Index([12, np.nan, 24], name='foo')

        rng = tm.box_expected(rng, box)
        expected = tm.box_expected(expected, box)

        result = rng / two_hours
        tm.assert_equal(result, expected)

    # ------------------------------------------------------------------
    # __floordiv__, __rfloordiv__

    def test_td64arr_floordiv_tdscalar(self, box, scalar_td):
        # GH#18831
        td1 = Series([timedelta(minutes=5, seconds=3)] * 3)
        td1.iloc[2] = np.nan

        expected = Series([0, 0, np.nan])

        td1 = tm.box_expected(td1, box)
        expected = tm.box_expected(expected, box)

        result = td1 // scalar_td
        tm.assert_equal(result, expected)

    def test_td64arr_rfloordiv_tdscalar(self, box, scalar_td):
        # GH#18831
        td1 = Series([timedelta(minutes=5, seconds=3)] * 3)
        td1.iloc[2] = np.nan

        expected = Series([1, 1, np.nan])

        td1 = tm.box_expected(td1, box)
        expected = tm.box_expected(expected, box)

        result = scalar_td // td1
        tm.assert_equal(result, expected)

    def test_td64arr_rfloordiv_tdscalar_explicit(self, box, scalar_td):
        # GH#18831
        td1 = Series([timedelta(minutes=5, seconds=3)] * 3)
        td1.iloc[2] = np.nan

        expected = Series([1, 1, np.nan])

        td1 = tm.box_expected(td1, box)
        expected = tm.box_expected(expected, box)

        # We can test __rfloordiv__ using this syntax,
        # see `test_timedelta_rfloordiv`
        result = td1.__rfloordiv__(scalar_td)
        tm.assert_equal(result, expected)

    def test_td64arr_floordiv_int(self, box):
        idx = TimedeltaIndex(np.arange(5, dtype='int64'))
        idx = tm.box_expected(idx, box)
        result = idx // 1
        tm.assert_equal(result, idx)

    def test_td64arr_floordiv_tdlike_scalar(self, two_hours, box):
        tdi = timedelta_range('1 days', '10 days', name='foo')
        expected = pd.Int64Index((np.arange(10) + 1) * 12, name='foo')

        tdi = tm.box_expected(tdi, box)
        expected = tm.box_expected(expected, box)

        result = tdi // two_hours
        tm.assert_equal(result, expected)

    # TODO: Is this redundant with test_td64arr_floordiv_tdlike_scalar?
    @pytest.mark.parametrize('scalar_td', [
        timedelta(minutes=10, seconds=7),
        Timedelta('10m7s'),
        Timedelta('10m7s').to_timedelta64()
    ], ids=lambda x: type(x).__name__)
    def test_td64arr_rfloordiv_tdlike_scalar(self, scalar_td, box):
        # GH#19125
        tdi = TimedeltaIndex(['00:05:03', '00:05:03', pd.NaT], freq=None)
        expected = pd.Index([2.0, 2.0, np.nan])

        tdi = tm.box_expected(tdi, box)
        expected = tm.box_expected(expected, box)

        res = tdi.__rfloordiv__(scalar_td)
        tm.assert_equal(res, expected)

        expected = pd.Index([0.0, 0.0, np.nan])
        expected = tm.box_expected(expected, box)

        res = tdi // (scalar_td)
        tm.assert_equal(res, expected)

    # ------------------------------------------------------------------
    # Operations with invalid others

    def test_td64arr_mul_tdscalar_invalid(self, box, scalar_td):
        td1 = Series([timedelta(minutes=5, seconds=3)] * 3)
        td1.iloc[2] = np.nan

        td1 = tm.box_expected(td1, box)

        # check that we are getting a TypeError
        # with 'operate' (from core/ops.py) for the ops that are not
        # defined
        pattern = 'operate|unsupported|cannot|not supported'
        with tm.assert_raises_regex(TypeError, pattern):
            td1 * scalar_td
        with tm.assert_raises_regex(TypeError, pattern):
            scalar_td * td1

    def test_td64arr_mul_too_short_raises(self, box):
        idx = TimedeltaIndex(np.arange(5, dtype='int64'))
        idx = tm.box_expected(idx, box)
        with pytest.raises(TypeError):
            idx * idx[:3]
        with pytest.raises(ValueError):
            idx * np.array([1, 2])

    def test_td64arr_mul_td64arr_raises(self, box):
        idx = TimedeltaIndex(np.arange(5, dtype='int64'))
        idx = tm.box_expected(idx, box)
        with pytest.raises(TypeError):
            idx * idx

    # ------------------------------------------------------------------
    # Operations with numeric others

    @pytest.mark.parametrize('one', [1, np.array(1), 1.0, np.array(1.0)])
    def test_td64arr_mul_numeric_scalar(self, box, one, tdser):
        # GH#4521
        # divide/multiply by integers
        expected = Series(['-59 Days', '-59 Days', 'NaT'],
                          dtype='timedelta64[ns]')

        tdser = tm.box_expected(tdser, box)
        expected = tm.box_expected(expected, box)

        result = tdser * (-one)
        tm.assert_equal(result, expected)
        result = (-one) * tdser
        tm.assert_equal(result, expected)

        expected = Series(['118 Days', '118 Days', 'NaT'],
                          dtype='timedelta64[ns]')
        expected = tm.box_expected(expected, box)

        result = tdser * (2 * one)
        tm.assert_equal(result, expected)
        result = (2 * one) * tdser
        tm.assert_equal(result, expected)

    @pytest.mark.parametrize('two', [2, 2.0, np.array(2), np.array(2.0)])
    def test_td64arr_div_numeric_scalar(self, box, two, tdser):
        # GH#4521
        # divide/multiply by integers
        expected = Series(['29.5D', '29.5D', 'NaT'], dtype='timedelta64[ns]')

        tdser = tm.box_expected(tdser, box)
        expected = tm.box_expected(expected, box)

        result = tdser / two
        tm.assert_equal(result, expected)

    @pytest.mark.parametrize('dtype', ['int64', 'int32', 'int16',
                                       'uint64', 'uint32', 'uint16', 'uint8',
                                       'float64', 'float32', 'float16'])
    @pytest.mark.parametrize('vector', [np.array([20, 30, 40]),
                                        pd.Index([20, 30, 40]),
                                        Series([20, 30, 40])],
                             ids=lambda x: type(x).__name__)
    @pytest.mark.parametrize('op', [operator.mul, ops.rmul])
    def test_td64arr_rmul_numeric_array(self, op, box, vector, dtype, tdser):
        # GH#4521
        # divide/multiply by integers
        vector = vector.astype(dtype)

        expected = Series(['1180 Days', '1770 Days', 'NaT'],
                          dtype='timedelta64[ns]')

        tdser = tm.box_expected(tdser, box, transpose=True)
        # TODO: Make this up-casting more systematic?
        box2 = Series if (box is pd.Index and type(vector) is Series) else box
        expected = tm.box_expected(expected, box2, transpose=True)

        result = op(vector, tdser)
        tm.assert_equal(result, expected)

    @pytest.mark.parametrize('dtype', ['int64', 'int32', 'int16',
                                       'uint64', 'uint32', 'uint16', 'uint8',
                                       'float64', 'float32', 'float16'])
    @pytest.mark.parametrize('vector', [np.array([20, 30, 40]),
                                        pd.Index([20, 30, 40]),
                                        Series([20, 30, 40])],
                             ids=lambda x: type(x).__name__)
    def test_td64arr_div_numeric_array(self, box, vector, dtype, tdser):
        # GH#4521
        # divide/multiply by integers
        vector = vector.astype(dtype)
        expected = Series(['2.95D', '1D 23H 12m', 'NaT'],
                          dtype='timedelta64[ns]')

        tdser = tm.box_expected(tdser, box, transpose=True)
        box = Series if (box is pd.Index and type(vector) is Series) else box
        expected = tm.box_expected(expected, box, transpose=True)

        result = tdser / vector
        tm.assert_equal(result, expected)

        with pytest.raises(TypeError):
            vector / tdser

    @pytest.mark.parametrize('names', [(None, None, None),
                                       ('Egon', 'Venkman', None),
                                       ('NCC1701D', 'NCC1701D', 'NCC1701D')])
    def test_td64arr_mul_int_series(self, box, names):
        # GH#19042 test for correct name attachment
        if box is pd.DataFrame and names[0] != names[1]:
            return

        tdi = TimedeltaIndex(['0days', '1day', '2days', '3days', '4days'],
                             name=names[0])
        # TODO: Should we be parametrizing over types for `ser` too?
        ser = Series([0, 1, 2, 3, 4], dtype=np.int64, name=names[1])

        expected = Series(['0days', '1day', '4days', '9days', '16days'],
                          dtype='timedelta64[ns]',
                          name=names[2])

        tdi = tm.box_expected(tdi, box, transpose=True)
        box = Series if (box is pd.Index and type(ser) is Series) else box
        expected = tm.box_expected(expected, box, transpose=True)

        result = ser * tdi
        tm.assert_equal(result, expected)

        result = tdi * ser
        tm.assert_equal(result, expected)

    # TODO: Should we be parametrizing over types for `ser` too?
    @pytest.mark.parametrize('names', [(None, None, None),
                                       ('Egon', 'Venkman', None),
                                       ('NCC1701D', 'NCC1701D', 'NCC1701D')])
    def test_float_series_rdiv_td64arr(self, box, names):
        # GH#19042 test for correct name attachment
        # TODO: the direct operation TimedeltaIndex / Series still
        # needs to be fixed.
        tdi = TimedeltaIndex(['0days', '1day', '2days', '3days', '4days'],
                             name=names[0])
        ser = Series([1.5, 3, 4.5, 6, 7.5], dtype=np.float64, name=names[1])

        expected = Series([tdi[n] / ser[n] for n in range(len(ser))],
                          dtype='timedelta64[ns]',
                          name=names[2])

        tdi = tm.box_expected(tdi, box)
        box = Series if (box is pd.Index and type(ser) is Series) else box
        expected = tm.box_expected(expected, box)

        result = ser.__rdiv__(tdi)
        if box is pd.DataFrame:
            # TODO: Should we skip this case sooner or test something else?
            assert result is NotImplemented
        else:
            tm.assert_equal(result, expected)


class TestTimedeltaArraylikeInvalidArithmeticOps(object):

    def test_td64arr_pow_invalid(self, scalar_td, box):
        td1 = Series([timedelta(minutes=5, seconds=3)] * 3)
        td1.iloc[2] = np.nan

        td1 = tm.box_expected(td1, box)

        # check that we are getting a TypeError
        # with 'operate' (from core/ops.py) for the ops that are not
        # defined
        pattern = 'operate|unsupported|cannot|not supported'
        with tm.assert_raises_regex(TypeError, pattern):
            scalar_td ** td1

        with tm.assert_raises_regex(TypeError, pattern):
            td1 ** scalar_td
